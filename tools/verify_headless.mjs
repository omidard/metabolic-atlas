// Headless verification: real flux search (LP + FVA + sampling) at 1440 and
// 390, expanded pathway card, interactive length filter, keyboard walk.
// Run: node tools/verify_headless.mjs   (serves docs/ itself)
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { extname, join } from 'node:path';
import pw from '/data/node_modules/playwright-core/index.js';
const { chromium } = pw;

const DOCS = '/data/metabolic_atlas_platform/docs';
const OUT = '/data/metabolic_atlas_platform/tools/shots';
const MIME = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.json': 'application/json', '.md': 'text/plain', '.tsv': 'text/plain' };

const server = createServer(async (req, res) => {
  try {
    let p = decodeURIComponent(new URL(req.url, 'http://x').pathname);
    if (p === '/') p = '/index.html';
    const data = await readFile(join(DOCS, p));
    res.writeHead(200, { 'content-type': MIME[extname(p)] || 'application/octet-stream' });
    res.end(data);
  } catch { res.writeHead(404); res.end('nf'); }
});
await new Promise(r => server.listen(8931, r));

const browser = await chromium.launch({
  executablePath: '/home/omidard/.cache/ms-playwright/chromium-1134/chrome-linux/chrome',
  args: ['--no-sandbox', '--use-angle=swiftshader'],
});
const errors = [];
const findings = [];
const ok = (name, cond, detail = '') => {
  findings.push(`${cond ? 'PASS' : 'FAIL'} ${name}${detail ? ' :: ' + detail : ''}`);
  console.log(findings[findings.length - 1]);
};

async function newPage(w, h) {
  const page = await browser.newPage({ viewport: { width: w, height: h }, deviceScaleFactor: 2 });
  page.on('console', m => { if (m.type() === 'error') errors.push(`[console ${w}px] ${m.text()}`); });
  page.on('pageerror', e => errors.push(`[pageerror ${w}px] ${e.message}`));
  return page;
}

// ---------- 1440px: full flux search ----------
const page = await newPage(1440, 900);
await page.goto('http://localhost:8931/index.html');
await page.waitForFunction(() => document.querySelector('#dataset-line')?.textContent.includes('union map'), null, { timeout: 60000 });
await page.waitForFunction(() => document.querySelector('#map-status')?.hidden === true, null, { timeout: 60000 }).catch(() => {});
const mapUp = await page.evaluate(() => document.querySelector('#map-status')?.hidden === true);
ok('3D map initialised (WebGL)', mapUp);

// run the regression search: pyr_e -> lac__D_e was falsely infeasible before
// the per-GEM step-realization fix (union id collision, D_LACtex); it must
// now test feasible on GCF_000178395.2 + Parageobacillus - glucose.
async function pick(kind, mid) {
  await page.fill(`#${kind}-input`, mid);
  await page.waitForSelector(`#${kind}-listbox li[data-mid="${mid}"]`, { timeout: 10000 });
  await page.click(`#${kind}-listbox li[data-mid="${mid}"]`);
}
await pick('sub', 'pyr_e');
await pick('prod', 'lac__D_e');
await page.click('#run-search');
await page.waitForFunction(() => document.querySelector('#results-summary')?.textContent.includes('pathway'), null, { timeout: 30000 });
ok('search ran', true, await page.$eval('#results-summary', el => el.textContent.trim().slice(0, 140)));
ok('length histogram rendered', await page.locator('.lenbar').count() > 0, `${await page.locator('.lenbar').count()} length bars`);
await page.screenshot({ path: `${OUT}/1440_discover.png`, fullPage: false });

// length filter interaction
const firstLen = await page.locator('.lenbar').first().getAttribute('data-len');
await page.locator('.lenbar').first().click();
const filteredTxt = await page.$eval('#pathway-list .morebar span', el => el.textContent);
ok('length filter filters the list', /matching pathways \(of/.test(filteredTxt), filteredTxt.trim());
await page.click('#len-clear');
ok('length filter clears', !(await page.$('#len-clear')));

// Simulate: choose GEM + medium, small sample count, run
await page.click('#nav-simulate');
await page.waitForSelector('#m2-gem', { timeout: 60000 });
await page.selectOption('#m2-gem', 'GCF_000178395.2');
await page.waitForFunction(() => document.querySelector('#m2-status')?.textContent.includes('GCF_000178395.2'), null, { timeout: 60000 });
await page.selectOption('#m2-medium', 'Parageobacillus - glucose');
await page.fill('#sim-samples', '10');
await page.waitForFunction(() => !document.querySelector('#sim-run').disabled, null, { timeout: 30000 });
ok('cancel hidden before run', await page.$eval('#sim-cancel', el => el.hidden));
const t0 = Date.now();
await page.click('#sim-run');
await page.waitForFunction(() => !document.querySelector('#sim-cancel') || document.querySelector('#sim-cancel').hidden === false, null, { timeout: 10000 }).catch(() => {});
await page.waitForFunction(() => /FVA \+ sampling on \d+ of \d+ feasible/.test(document.querySelector('#results-summary .mode2-status')?.textContent || ''), null, { timeout: 420000 });
const runLine = await page.$eval('#results-summary .mode2-status', el => el.textContent.replace(/\s+/g, ' ').trim());
ok('flux run completed', true, `${Math.round((Date.now() - t0) / 1000)}s :: ${runLine.slice(0, 240)}`);
ok('cancel re-hidden after run', await page.$eval('#sim-cancel', el => el.hidden));

// feasibility chips + tiers + legend
const chips = await page.$$eval('#results-body .feas-chip', els => els.map(e => e.textContent.trim()));
ok('feasible chip with flux number', chips.some(c => /feasible · [\d.]+ mmol/.test(c)), chips.slice(0, 6).join(' | '));
const legendTxt = await page.$eval('#map-legend', el => el.textContent);
ok('map legend has flux-search tiers', /Flux search/.test(legendTxt) && /shortest feasible/.test(legendTxt));
ok('map legend has biomass core', /representative biomass/.test(legendTxt));
await page.screenshot({ path: `${OUT}/1440_simulate_run.png` });

// expand the best feasible card: FVA + sampling chart, genes, subsystems
const feasCard = page.locator('#results-body .pcard', { has: page.locator('.feas-chip.ok') }).first();
await feasCard.locator('summary').click();
await page.waitForTimeout(1500);   // enrichment fetch (cached gem) + chart render
const cardHtml = await feasCard.evaluate(el => el.innerHTML);
ok('FVA vs sampling chart in card', /FVA range vs sampled flux/.test(cardHtml));
ok('sampling denominator stated', /of 10 requested samples solved|sampling not computed/.test(cardHtml));
ok('yield with basis in card', /yield .* mol\/mol|yield not computed/.test(cardHtml));
ok('gene symbols in card', /<strong>[a-zA-Z0-9]+<\/strong> <span class="mono locus">/.test(cardHtml) || /no gene rule/.test(cardHtml));
ok('subsystem chips in card', /subchip/.test(cardHtml));
ok('per-reaction pFBA + FVA in step rows', /v [\d.\-e]+ · FVA \[/.test(cardHtml));
await feasCard.screenshot({ path: `${OUT}/1440_card_expanded.png` });
await page.locator('#map-pane').screenshot({ path: `${OUT}/1440_map_tiers.png` });
ok('cancel visually hidden after run (CSS)', await page.$eval('#sim-cancel', el => getComputedStyle(el).display === 'none'));

// histogram now shows feasible split + feasible-only + flux sort work
const legendCtl = await page.$eval('#list-controls', el => el.textContent);
ok('histogram legend has run denominators', /tested feasible/.test(legendCtl), legendCtl.match(/feasibility on[^)]+\)/)?.[0] || '');
await page.check('#feas-only');
const feasOnlyTxt = await page.$eval('#pathway-list .morebar span', el => el.textContent);
ok('feasible-only filter honest', /of .* found/.test(feasOnlyTxt) || /pathways\./.test(feasOnlyTxt), feasOnlyTxt.trim());
await page.selectOption('#sort-select', 'flux');
const firstChip = await page.$eval('#pathway-list .pcard .feas-chip', el => el.textContent);
ok('flux sort puts a feasible pathway first', /feasible/.test(firstChip), firstChip.trim());
await page.uncheck('#feas-only');
await page.screenshot({ path: `${OUT}/1440_sorted.png` });

// per-card deep run button exists for an untested pathway
const naBtn = await page.$$eval('#pathway-list .pcard .mode2-slot button', els => els.map(e => e.textContent.trim()));
ok('per-pathway test/deep buttons exist', naBtn.some(t => /Test this pathway|Run FVA \+ flux sampling|Re-run FVA/.test(t)), naBtn.slice(0, 4).join(' | '));

// Model stage: gene-symbol card + GPR symbols
await page.click('#nav-model');
await page.waitForFunction(() => document.querySelector('#gem-body') && !document.querySelector('#gem-body').hidden, null, { timeout: 30000 });
const dashTxt = await page.$eval('#gem-body', el => el.textContent);
ok('gene-symbol coverage card', /Gene symbols/.test(dashTxt));
const gprHtml = await page.$eval('#rxn-table tbody', el => el.innerHTML);
ok('GPR cells carry symbol + locus', /<strong>[a-zA-Z0-9]+<\/strong>/.test(gprHtml));
await page.screenshot({ path: `${OUT}/1440_model_dashboard.png`, fullPage: false });

// keyboard walk: tab reaches search input, a length bar, and a card summary
await page.click('#nav-discover');
const kb = await page.evaluate(() => {
  const els = ['#sub-input', '.lenbar', '#sort-select', '#feas-only', '.pcard summary', '#run-search'];
  return els.map(sel => {
    const el = document.querySelector(sel);
    if (!el) return `${sel}: MISSING`;
    el.focus();
    return `${sel}: ${document.activeElement === el ? 'focusable' : 'NOT focusable'}`;
  });
});
ok('keyboard: all list controls focusable', kb.every(s => /focusable/.test(s) && !/NOT/.test(s)), kb.join('; '));
await page.close();

// ---------- 390px ----------
const p390 = await newPage(390, 844);
await p390.goto('http://localhost:8931/index.html');
await p390.waitForFunction(() => document.querySelector('#dataset-line')?.textContent.includes('union map'), null, { timeout: 60000 });
for (const sel of ['#sub-input', '#prod-input', '#run-search', '#nav-simulate', '#nav-model']) {
  const box = await p390.locator(sel).boundingBox();
  ok(`390px: ${sel} operable in-viewport`, !!box && box.x >= 0 && box.x + box.width <= 391, box ? `x=${Math.round(box.x)} w=${Math.round(box.width)}` : 'no box');
}
await p390.click('#example-row button');
await p390.waitForFunction(() => document.querySelector('#results-summary')?.textContent.includes('pathway'), null, { timeout: 30000 });
const lb = await p390.locator('.lenbar').first().boundingBox();
ok('390px: length bars reachable', !!lb && lb.x + lb.width <= 391);
await p390.screenshot({ path: `${OUT}/390_discover.png` });
await p390.click('#nav-simulate');
await p390.waitForSelector('#sim-run', { timeout: 60000 });
const sb = await p390.locator('#sim-run').boundingBox();
ok('390px: #sim-run operable in-viewport', !!sb && sb.x >= 0 && sb.x + sb.width <= 391, sb ? `x=${Math.round(sb.x)} w=${Math.round(sb.width)}` : 'no box');
const ss = await p390.locator('#sim-samples').boundingBox();
ok('390px: #sim-samples operable in-viewport', !!ss && ss.x >= 0 && ss.x + ss.width <= 391);
await p390.screenshot({ path: `${OUT}/390_simulate.png` });
await p390.close();

await browser.close();
server.close();

console.log('\n---- console/page errors:', errors.length);
errors.slice(0, 12).forEach(e => console.log(e));
const fails = findings.filter(f => f.startsWith('FAIL')).length;
console.log(`\n==== ${findings.length - fails} PASS / ${fails} FAIL`);
process.exit(fails || errors.length ? 1 : 0);
