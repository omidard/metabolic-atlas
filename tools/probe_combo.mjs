// Probe: which substrate->product x GEM x medium combos yield feasible
// pathways (real solves through the UI). Prints the run summary per combo.
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { extname, join } from 'node:path';
import pw from '/data/node_modules/playwright-core/index.js';
const { chromium } = pw;

const DOCS = '/data/metabolic_atlas_platform/docs';
const MIME = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json' };
const server = createServer(async (req, res) => {
  try {
    let p = decodeURIComponent(new URL(req.url, 'http://x').pathname);
    if (p === '/') p = '/index.html';
    const data = await readFile(join(DOCS, p));
    res.writeHead(200, { 'content-type': MIME[extname(p)] || 'application/octet-stream' });
    res.end(data);
  } catch { res.writeHead(404); res.end('nf'); }
});
await new Promise(r => server.listen(8932, r));

const combos = [
  ['pyr_e', 'lac__D_e', 'GCF_000178395.2', 'Parageobacillus - glucose'],
  ['meoh_e', 'ac_e', 'GCF_000807675.2', 'Eubacterium - methanol'],
  ['pyr_e', 'ac_e', 'GCF_000007565.2', 'Pseudomonas - glucose'],
  ['glc__D_e', 'lac__D_e', 'GCF_000178395.2', 'Parageobacillus - glucose'],
];

const browser = await chromium.launch({
  executablePath: '/home/omidard/.cache/ms-playwright/chromium-1134/chrome-linux/chrome',
  args: ['--no-sandbox', '--use-angle=swiftshader'],
});
const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
await page.goto('http://localhost:8932/index.html');
await page.waitForFunction(() => document.querySelector('#dataset-line')?.textContent.includes('union map'), null, { timeout: 60000 });
// load simulate module once
await page.click('#nav-simulate');
await page.waitForSelector('#m2-gem', { timeout: 60000 });
await page.fill('#sim-samples', '10');

async function pick(kind, mid) {
  await page.click('#nav-discover');
  await page.fill(`#${kind}-input`, mid);
  await page.waitForSelector(`#${kind}-listbox li[data-mid="${mid}"]`, { timeout: 10000 });
  await page.click(`#${kind}-listbox li[data-mid="${mid}"]`);
}

for (const [sub, prod, gem, medium] of combos) {
  try {
    await pick('sub', sub);
    await pick('prod', prod);
    await page.click('#run-search');
    await page.waitForFunction(() => /pathway/.test(document.querySelector('#results-summary')?.textContent || ''), null, { timeout: 30000 });
    const found = await page.$eval('#results-summary', el => el.textContent.replace(/\s+/g, ' ').trim().slice(0, 90));
    if (/No pathway found/.test(found)) { console.log(`-- ${sub}->${prod} ${gem}: NO PATHWAYS (${found})`); continue; }
    await page.click('#nav-simulate');
    await page.selectOption('#m2-gem', gem);
    await page.waitForFunction((g) => document.querySelector('#m2-status')?.textContent.includes(g), gem, { timeout: 60000 });
    await page.selectOption('#m2-medium', medium);
    await page.waitForFunction(() => !document.querySelector('#sim-run').disabled, null, { timeout: 30000 });
    await page.click('#sim-run');
    await page.waitForFunction(() => /feasible|does not contain/.test(document.querySelector('#results-summary .mode2-status')?.textContent || '') && /FVA \+ sampling on|does not contain/.test(document.querySelector('#results-summary .mode2-status')?.textContent || ''), null, { timeout: 420000 });
    const line = await page.$eval('#results-summary .mode2-status', el => el.textContent.replace(/\s+/g, ' ').trim());
    console.log(`== ${sub}->${prod} ${gem} on ${medium}\n   ${found}\n   ${line.slice(0, 260)}`);
  } catch (e) {
    console.log(`!! ${sub}->${prod} ${gem}: ${e.message.split('\n')[0]}`);
  }
}
await browser.close();
server.close();
