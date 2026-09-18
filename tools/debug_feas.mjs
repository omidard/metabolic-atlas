// Debug: why is a trivial pathway infeasible? Runs the UI search, then digs
// into the LP pieces in the page context with the real modules.
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
await new Promise(r => server.listen(8934, r));

const browser = await chromium.launch({
  executablePath: '/home/omidard/.cache/ms-playwright/chromium-1134/chrome-linux/chrome',
  args: ['--no-sandbox', '--use-angle=swiftshader'],
});
const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
page.on('console', m => console.log('[pg]', m.text().slice(0, 400)));
await page.goto('http://localhost:8934/index.html');
await page.waitForFunction(() => document.querySelector('#dataset-line')?.textContent.includes('union map'), null, { timeout: 60000 });

const out = await page.evaluate(async () => {
  const { loadGraph, loadGraphMeta, loadGem, loadMedia } = await import('./assets/data.js');
  const { searchPathways } = await import('./assets/search.js');
  const fba = await import('./assets/fba.js');
  const graph = await loadGraph();
  const meta = await loadGraphMeta();
  const gem = await loadGem('GCF_000178395.2');
  const media = await loadMedia();
  const res = searchPathways(graph, meta, 'pyr_e', 'lac__D_e',
    ['Parageobacillus thermoglucosidasius', 'Pseudomonas putida', 'Cupriavidus necator', 'Eubacterium limosum']);
  const lines = [];
  lines.push(`pathways: ${res.pathways.length}`);
  const gemRxns = new Set(gem.reactions.map(r => r.id));
  // medium bounds like mode2.effectiveMedium with substrate swap
  const def = media['Parageobacillus - glucose'];
  const bounds = {};
  for (const [ex, lb] of Object.entries(def.components || {})) if (!ex.startsWith('_')) bounds[ex] = lb;
  for (const [ex, lb] of Object.entries(def.supplements || {})) if (!ex.startsWith('_') && !(ex in bounds)) bounds[ex] = lb;
  bounds['EX_glc__D_e'] = 0;
  bounds['EX_pyr_e'] = def.carbon_cap ?? -10;
  const target = fba.productTarget(gem, 'lac__D_e');
  lines.push(`target: ${JSON.stringify(target && { kind: target.kind, id: target.id })}`);
  for (let i = 0; i < Math.min(10, res.pathways.length); i++) {
    const pw = res.pathways[i];
    const { stepCons, missingStep } = fba.stepConstraints(gem, pw);
    const steps = pw.steps.map((st, k) => `${st.from}->${st.to}[${st.rxns.map(r => r.id + (gemRxns.has(r.id) ? '' : '(absent)') + ':' + r.dir).join(',')}]`).join(' ');
    if (!stepCons) { lines.push(`#${i + 1} len ${pw.len} NOT CARRIED (step ${missingStep + 1}) :: ${steps}`); continue; }
    const r = await fba.pathwayFeasibility(gem, 'GCF_000178395.2', bounds, pw, target);
    lines.push(`#${i + 1} len ${pw.len} feasible=${r.feasible} z=${r.productFlux} status=${r.status} :: ${steps}`);
    if (!r.feasible && i < 3) {
      const diag = await fba.diagnoseSteps(gem, 'GCF_000178395.2', bounds, pw, target);
      lines.push(`   diag: ${diag.map(d => `step${d.step + 1} ok=${d.ok} flux=${d.flux}`).join('; ')}`);
    }
  }
  // sanity: growth and plain product max without step constraints
  const g = await fba.maxGrowth(gem, 'GCF_000178395.2', bounds);
  lines.push(`growth: optimal=${g.optimal} mu=${g.mu}`);
  const glpk = await fba.getGLPK();
  const lp = fba.buildLP(glpk, gem, bounds, { acc: 'GCF_000178395.2', extraCols: target.extraCols, objective: { direction: 'max', vars: [{ name: target.id, coef: 1 }] } });
  const s = await fba.solveLP(glpk, lp);
  lines.push(`product max (no step constraints): optimal=${s.optimal} z=${s.z} status=${s.status}`);
  return lines.join('\n');
});
console.log(out);
await browser.close();
server.close();
