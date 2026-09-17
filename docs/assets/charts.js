// Shared chart helpers: small designed SVG/HTML charts rendered from data the
// page has already loaded or solved. One fixed palette; every count carries its
// denominator in the caller's title line; absent values are skipped, never
// drawn as zero. All text at or above 12px at the rendered size (SVG widths are
// fixed at build time, never scaled down by the viewport).

export const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

// ---- one palette for the 14 pathway groups (union map, legends, all charts)
export const GROUP_COLORS = {
  'Glycolysis': '#A93F32',
  'TCA': '#A65D1E',
  'PPP': '#8F7A1F',
  'ETC': '#5F8A26',
  'Fatty acid': '#2E8A50',
  'Amino acid': '#2A8A80',
  'Nucleotide': '#2E7A9E',
  'Cofactor': '#3A62B0',
  'Lipid': '#6551B5',
  'Cell envelope': '#8A4BA8',
  'Carbon': '#A8449A',
  'Ion transport': '#AE3F6B',
  'Transport': '#7A6A52',
  'Other': '#98948C',
};
export const NEUTRAL_BAR = '#98948C';

const nf = new Intl.NumberFormat('en-US');
const fmtN = (x) => nf.format(x);
const pct = (n, d) => (d > 0 ? Math.round(n / d * 100) : null);

// A chart block: topic-phrase title, the chart, an optional caption line.
export function chartBlock(title, inner, caption = '') {
  return `<figure class="chart">
    <figcaption class="chart-title">${title}</figcaption>
    ${inner}
    ${caption ? `<div class="chart-cap">${caption}</div>` : ''}
  </figure>`;
}

// ---- horizontal bars (HTML): label · bar · value (+ share of total)
// rows: [{label, value, color, note}], total = denominator for the share text.
export function hBars(rows, { total = null, valueFmt = fmtN, showPct = true, barMax = null } = {}) {
  const max = barMax ?? Math.max(...rows.map(r => Math.abs(r.value)), 1e-12);
  return `<div class="ch-bars">${rows.map(r => {
    const w = Math.max(Math.abs(r.value) / max * 100, 0.5);
    const share = (showPct && total != null) ? pct(r.value, total) : null;
    return `<div class="ch-row">
      <span class="ch-label" title="${esc(r.label)}">${esc(r.label)}</span>
      <span class="ch-track"><span class="ch-bar" style="width:${w.toFixed(1)}%;background:${r.color || NEUTRAL_BAR}"></span></span>
      <span class="ch-val">${valueFmt(r.value)}${share != null ? ` <span class="ch-share">(${share}%)</span>` : ''}${r.note ? ` <span class="ch-share">${esc(r.note)}</span>` : ''}</span>
    </div>`;
  }).join('')}</div>`;
}

// ---- single 100% stacked bar with legend. segs: [{label, value, color}]
export function stackBar(segs, { total = null } = {}) {
  const sum = segs.reduce((a, s) => a + s.value, 0) || 1;
  const den = total ?? sum;
  return `<div class="ch-stack" role="img" aria-label="${esc(segs.map(s => `${s.label} ${s.value} of ${den}`).join(', '))}">
    <div class="ch-stack-track">${segs.filter(s => s.value > 0).map(s =>
      `<span class="ch-stack-seg" style="width:${(s.value / sum * 100).toFixed(2)}%;background:${s.color}" title="${esc(s.label)}: ${fmtN(s.value)} of ${fmtN(den)}"></span>`).join('')}
    </div>
    <div class="ch-legend">${segs.map(s =>
      `<span class="ch-lg"><span class="ch-sw" style="background:${s.color}"></span>${esc(s.label)} ${fmtN(s.value)}<span class="ch-share"> (${pct(s.value, den)}%)</span></span>`).join('')}</div>
  </div>`;
}

// ---- donut (SVG, fixed size) with legend. segs: [{label, value, color}]
export function donut(segs, { size = 132, centerTop = '', centerBottom = '' } = {}) {
  const sum = segs.reduce((a, s) => a + s.value, 0);
  const R = size / 2, r0 = R - 13, C = 2 * Math.PI * r0;
  let off = 0;
  const rings = segs.filter(s => s.value > 0).map(s => {
    const frac = s.value / (sum || 1);
    const dash = `${(frac * C).toFixed(2)} ${(C - frac * C).toFixed(2)}`;
    const el = `<circle r="${r0}" cx="${R}" cy="${R}" fill="none" stroke="${s.color}" stroke-width="13"
      stroke-dasharray="${dash}" stroke-dashoffset="${(-off * C).toFixed(2)}" transform="rotate(-90 ${R} ${R})"/>`;
    off += frac;
    return el;
  }).join('');
  return `<div class="ch-donut">
    <svg width="${size}" height="${size}" role="img" aria-label="${esc(segs.map(s => `${s.label} ${s.value}`).join(', '))}">
      ${rings}
      <text x="${R}" y="${R - 2}" text-anchor="middle" class="ch-svg-strong">${esc(centerTop)}</text>
      <text x="${R}" y="${R + 14}" text-anchor="middle" class="ch-svg-dim">${esc(centerBottom)}</text>
    </svg>
    <div class="ch-legend ch-legend-col">${segs.map(s =>
      `<span class="ch-lg"><span class="ch-sw" style="background:${s.color}"></span>${esc(s.label)} ${fmtN(s.value)}<span class="ch-share"> (${sum ? pct(s.value, sum) : 0}%)</span></span>`).join('')}</div>
  </div>`;
}

// ---- fraction bar: num of den, with a % readout beside it.
export function fractionBar(num, den, { color = 'var(--accent)', label = '' } = {}) {
  if (num == null || den == null || den <= 0) {
    return `<div class="ch-fraction"><span class="ch-val">not computed</span></div>`;
  }
  const p = num / den * 100;
  return `<div class="ch-fraction" role="img" aria-label="${esc(label)} ${fmtN(num)} of ${fmtN(den)} (${p.toFixed(1)}%)">
    <span class="ch-track ch-track-wide"><span class="ch-bar" style="width:${Math.min(p, 100).toFixed(1)}%;background:${color}"></span></span>
    <span class="ch-val">${fmtN(num)} of ${fmtN(den)} <span class="ch-share">(${p.toFixed(1)}%)</span></span>
  </div>`;
}

// ---- dumbbell / interval rows (SVG): rows [{label, lo, hi}] in data units.
export function intervals(rows, { width = 360, unit = '', min = 0, max = null, color = 'var(--accent-ink)' } = {}) {
  const usable = rows.filter(r => r.lo != null && r.hi != null);
  if (!usable.length) return '<p class="status">not computed</p>';
  const hiMax = max ?? ((Math.max(...usable.map(r => r.hi)) * 1.15) || 1);
  const mL = 110, mR = 66, rowH = 26, mT = 6, mB = 22;
  const H = mT + rows.length * rowH + mB;
  const iw = width - mL - mR;
  const X = (v) => mL + (v - min) / (hiMax - min) * iw;
  const ticks = niceTicks(min, hiMax, 4);
  return `<svg width="${width}" height="${H}" role="img" aria-label="${esc(rows.map(r => `${r.label}: ${r.lo} to ${r.hi} ${unit}`).join('; '))}">
    ${ticks.map(t => `<line x1="${X(t)}" x2="${X(t)}" y1="${mT}" y2="${H - mB}" stroke="var(--line)" stroke-width="1"/>
      <text x="${X(t)}" y="${H - 7}" text-anchor="middle" class="ch-svg-dim">${tickLabel(t)}</text>`).join('')}
    ${rows.map((r, i) => {
      const y = mT + i * rowH + rowH / 2;
      if (r.lo == null || r.hi == null) {
        return `<text x="0" y="${y + 4}" class="ch-svg-label">${esc(r.label)}</text>
          <text x="${mL}" y="${y + 4}" class="ch-svg-dim">not computed</text>`;
      }
      return `<text x="0" y="${y + 4}" class="ch-svg-label">${esc(r.label)}</text>
        <line x1="${X(r.lo)}" x2="${X(r.hi)}" y1="${y}" y2="${y}" stroke="${color}" stroke-width="3" stroke-linecap="round"/>
        <circle cx="${X(r.lo)}" cy="${y}" r="4" fill="${color}"/>
        <circle cx="${X(r.hi)}" cy="${y}" r="4" fill="${color}"/>
        <text x="${X(r.hi) + 8}" y="${y + 4}" class="ch-svg-val">${r.lo === r.hi ? r.hi : `${r.lo} to ${r.hi}`}</text>`;
    }).join('')}
  </svg>`;
}

// ---- histogram over integer or binned x. bins: [{x, count, color}]
export function histogram(bins, { width = 360, height = 130, xLabel = '', yLabel = 'count', barLabel = null } = {}) {
  const mL = 40, mR = 8, mT = 8, mB = 34;
  const iw = width - mL - mR, ih = height - mT - mB;
  const cMax = Math.max(...bins.map(b => b.count), 1);
  const bw = iw / bins.length;
  return `<svg width="${width}" height="${height}" role="img" aria-label="Histogram of ${esc(xLabel)}: ${esc(bins.map(b => `${b.x}: ${b.count}`).join(', '))}">
    <line x1="${mL}" x2="${width - mR}" y1="${mT + ih}" y2="${mT + ih}" stroke="var(--line-strong)"/>
    ${bins.map((b, i) => {
      const h = b.count / cMax * ih;
      const x = mL + i * bw;
      return `${b.count ? `<rect x="${(x + 1).toFixed(1)}" y="${(mT + ih - h).toFixed(1)}" width="${Math.max(bw - 2, 1).toFixed(1)}" height="${h.toFixed(1)}" fill="${b.color || 'var(--accent)'}"/>` : ''}
        ${b.count && (barLabel ? barLabel(b) : true) ? `<text x="${(x + bw / 2).toFixed(1)}" y="${(mT + ih - h - 4).toFixed(1)}" text-anchor="middle" class="ch-svg-dim">${fmtN(b.count)}</text>` : ''}
        ${b.tick !== false ? `<text x="${(x + bw / 2).toFixed(1)}" y="${height - 18}" text-anchor="middle" class="ch-svg-dim">${esc(String(b.x))}</text>` : ''}`;
    }).join('')}
    <text x="${mL + iw / 2}" y="${height - 3}" text-anchor="middle" class="ch-svg-dim">${esc(xLabel)}</text>
  </svg>`;
}

// ---- heatmap: rows x cols with per-cell count and shading toward `color`.
// cells[i][j] = number or null (renders as a dot for "none").
export function heatmap(rowLabels, colLabels, cells, {
  cellW = 46, cellH = 26, labelW = 150, color = [42, 98, 176], showValues = true,
  colColors = null, topLabelH = 84, valueMax = null,
} = {}) {
  const flat = cells.flat().filter(v => v != null);
  const vMax = valueMax ?? Math.max(...flat, 1);
  const W = labelW + colLabels.length * cellW + 8;
  const H = topLabelH + rowLabels.length * cellH + 6;
  const shade = (v) => {
    const t = Math.pow(v / vMax, 0.6);
    return `rgba(${color[0]},${color[1]},${color[2]},${(0.08 + 0.84 * t).toFixed(3)})`;
  };
  return `<svg width="${W}" height="${H}" role="img" aria-label="Heatmap, ${rowLabels.length} rows by ${colLabels.length} columns.">
    ${colLabels.map((c, j) => `<text transform="translate(${labelW + j * cellW + cellW / 2 + 4} ${topLabelH - 6}) rotate(-45)"
      class="ch-svg-dim" text-anchor="start" ${colColors && colColors[j] ? `fill="${colColors[j]}" style="fill:${colColors[j]}"` : ''}>${esc(c)}</text>`).join('')}
    ${rowLabels.map((r, i) => `<text x="0" y="${topLabelH + i * cellH + cellH / 2 + 4}" class="ch-svg-label">${esc(r)}</text>`).join('')}
    ${cells.map((row, i) => row.map((v, j) => {
      const x = labelW + j * cellW, y = topLabelH + i * cellH;
      if (v == null || v === 0) {
        return `<rect x="${x + 1}" y="${y + 1}" width="${cellW - 2}" height="${cellH - 2}" fill="var(--surface-2)"/>
          <circle cx="${x + cellW / 2}" cy="${y + cellH / 2}" r="1.5" fill="var(--line-strong)"/>`;
      }
      const dark = v / vMax > 0.5;
      return `<rect x="${x + 1}" y="${y + 1}" width="${cellW - 2}" height="${cellH - 2}" fill="${shade(v)}"/>
        ${showValues ? `<text x="${x + cellW / 2}" y="${y + cellH / 2 + 4}" text-anchor="middle" class="${dark ? 'ch-svg-oncolor' : 'ch-svg-dim'}">${fmtN(v)}</text>` : ''}`;
    }).join('')).join('')}
  </svg>`;
}

// ---- boolean presence matrix: rows x cols, on/off cells coloured per column.
export function presenceMatrix(rowLabels, colLabels, cells, colColors, { cellW = 15, cellH = 18, labelW = 92, topLabelH = 8 } = {}) {
  const W = labelW + colLabels.length * cellW + 6;
  const H = topLabelH + rowLabels.length * cellH + 4;
  return `<svg width="${W}" height="${H}" role="img" aria-label="Presence matrix, ${rowLabels.length} pathways by ${colLabels.length} GEMs.">
    ${rowLabels.map((r, i) => `<text x="0" y="${topLabelH + i * cellH + cellH / 2 + 4}" class="ch-svg-label">${esc(r)}</text>`).join('')}
    ${cells.map((row, i) => row.map((on, j) => {
      const x = labelW + j * cellW, y = topLabelH + i * cellH;
      return on
        ? `<rect x="${x + 1}" y="${y + 2}" width="${cellW - 3}" height="${cellH - 5}" rx="2" fill="${colColors[j]}"><title>${esc(colLabels[j])}: carries every step</title></rect>`
        : `<rect x="${x + 1}" y="${y + 2}" width="${cellW - 3}" height="${cellH - 5}" rx="2" fill="var(--surface-2)" stroke="var(--line)"><title>${esc(colLabels[j])}: missing at least one step</title></rect>`;
    }).join('')).join('')}
  </svg>`;
}

// ---- diverging horizontal bars: positive right (accent), negative left (grey).
// rows: [{label, value, note}] sorted by caller.
export function divergingBars(rows, { width = 520, posColor = 'var(--accent)', negColor = '#98948C', unit = '', valueFmt = (v) => v.toFixed(3), posLabel = '', negLabel = '' } = {}) {
  const rowH = 24, mT = 20, mB = 8, labelW = 128, valW = 74;
  const H = mT + rows.length * rowH + mB;
  const iw = width - labelW - valW - 16;
  const vMax = Math.max(...rows.map(r => Math.abs(r.value)), 1e-12);
  const x0 = labelW + iw / 2;
  return `<svg width="${width}" height="${H}" role="img" aria-label="Diverging bars${unit ? ', ' + esc(unit) : ''}.">
    <line x1="${x0}" x2="${x0}" y1="${mT - 6}" y2="${H - mB}" stroke="var(--line-strong)"/>
    ${negLabel ? `<text x="${x0 - 6}" y="${mT - 8}" text-anchor="end" class="ch-svg-dim">${esc(negLabel)}</text>` : ''}
    ${posLabel ? `<text x="${x0 + 6}" y="${mT - 8}" class="ch-svg-dim">${esc(posLabel)}</text>` : ''}
    ${rows.map((r, i) => {
      const y = mT + i * rowH;
      const w = Math.abs(r.value) / vMax * (iw / 2 - 4);
      const pos = r.value >= 0;
      return `<text x="0" y="${y + rowH / 2 + 4}" class="ch-svg-label">${esc(r.label)}</text>
        <rect x="${pos ? x0 : x0 - w}" y="${y + 4}" width="${Math.max(w, 1).toFixed(1)}" height="${rowH - 9}" fill="${pos ? posColor : negColor}"/>
        <text x="${pos ? x0 + w + 5 : x0 - w - 5}" y="${y + rowH / 2 + 4}" ${pos ? '' : 'text-anchor="end"'} class="ch-svg-val">${valueFmt(r.value)}</text>`;
    }).join('')}
  </svg>`;
}

// ---- quantile rows: min…max whisker, p5-p95 band, median tick, per row.
// rows: [{label, min, p5, median, p95, max}]
export function quantileRows(rows, { width = 520, unit = '', color = 'var(--accent)' } = {}) {
  const rowH = 26, mT = 24, mB = 24, labelW = 128, mR = 12;
  const H = mT + rows.length * rowH + mB;
  const iw = width - labelW - mR;
  const lo = Math.min(...rows.map(r => r.min), 0);
  const hi = Math.max(...rows.map(r => r.max), 1e-9);
  const X = (v) => labelW + (v - lo) / (hi - lo || 1) * iw;
  const ticks = niceTicks(lo, hi, 5);
  return `<svg width="${width}" height="${H}" role="img" aria-label="Per-reaction flux distributions: 5th to 95th percentile band, median tick, minimum to maximum whisker.">
    ${ticks.map(t => `<line x1="${X(t)}" x2="${X(t)}" y1="${mT - 4}" y2="${H - mB}" stroke="var(--line)"/>
      <text x="${X(t)}" y="${H - 9}" text-anchor="middle" class="ch-svg-dim">${tickLabel(t)}</text>`).join('')}
    ${lo < 0 ? `<line x1="${X(0)}" x2="${X(0)}" y1="${mT - 4}" y2="${H - mB}" stroke="var(--line-strong)"/>` : ''}
    <text x="${labelW}" y="${mT - 10}" class="ch-svg-dim">flux (${esc(unit)}) · band = 5-95%, tick = median, whisker = min to max</text>
    ${rows.map((r, i) => {
      const y = mT + i * rowH + rowH / 2;
      return `<text x="0" y="${y + 4}" class="ch-svg-label">${esc(r.label)}</text>
        <line x1="${X(r.min)}" x2="${X(r.max)}" y1="${y}" y2="${y}" stroke="var(--line-strong)" stroke-width="1"/>
        <rect x="${X(r.p5)}" y="${y - 6}" width="${Math.max(X(r.p95) - X(r.p5), 1.5).toFixed(1)}" height="12" rx="2" fill="${color}" fill-opacity="0.75"/>
        <line x1="${X(r.median)}" x2="${X(r.median)}" y1="${y - 8}" y2="${y + 8}" stroke="var(--ink)" stroke-width="2"/>`;
    }).join('')}
  </svg>`;
}

// ---- grouped/paired vertical bars: groups: [{label, bars: [{label, value, color}]}]
export function pairedBars(groups, { width = 420, height = 170, unit = '', valueFmt = (v) => v == null ? 'not computed' : v.toFixed(3) } = {}) {
  const mL = 8, mR = 8, mT = 24, mB = 40;
  const iw = width - mL - mR, ih = height - mT - mB;
  const all = groups.flatMap(g => g.bars.map(b => b.value)).filter(v => v != null);
  const vMax = Math.max(...all, 1e-9);
  const gw = iw / groups.length;
  return `<svg width="${width}" height="${height}" role="img" aria-label="${esc(groups.map(g => `${g.label}: ${g.bars.map(b => `${b.label} ${valueFmt(b.value)}`).join(', ')}`).join('; '))}">
    <line x1="${mL}" x2="${width - mR}" y1="${mT + ih}" y2="${mT + ih}" stroke="var(--line-strong)"/>
    ${groups.map((g, i) => {
      const bx = mL + i * gw;
      const bw = Math.min((gw - 18) / g.bars.length, 44);
      const start = bx + (gw - bw * g.bars.length) / 2;
      return `${g.bars.map((b, j) => {
        if (b.value == null) {
          return `<text x="${start + j * bw + bw / 2}" y="${mT + ih - 6}" text-anchor="middle" class="ch-svg-dim" transform="rotate(-90 ${start + j * bw + bw / 2} ${mT + ih - 6})">not computed</text>`;
        }
        const h = Math.abs(b.value) / vMax * ih;
        return `<rect x="${(start + j * bw + 2).toFixed(1)}" y="${(mT + ih - h).toFixed(1)}" width="${(bw - 4).toFixed(1)}" height="${Math.max(h, 1).toFixed(1)}" fill="${b.color}"/>
          <text x="${start + j * bw + bw / 2}" y="${(mT + ih - h - 5).toFixed(1)}" text-anchor="middle" class="ch-svg-val">${valueFmt(b.value)}</text>`;
      }).join('')}
      <text x="${bx + gw / 2}" y="${height - 22}" text-anchor="middle" class="ch-svg-label" style="text-anchor:middle">${esc(g.label)}</text>
      <text x="${bx + gw / 2}" y="${height - 6}" text-anchor="middle" class="ch-svg-dim">${esc(g.sub || '')}</text>`;
    }).join('')}
    <text x="${mL}" y="14" class="ch-svg-dim">${esc(unit)}</text>
  </svg>`;
}

// ---- headline stat card with an optional comparison bullet vs a mean.
export function statCard(k, v, d, { mean = null, meanLabel = '', color = 'var(--accent)' } = {}) {
  let bullet = '';
  if (mean != null && typeof v === 'number') {
    const max = Math.max(v, mean) * 1.1;
    const dpct = mean > 0 ? ((v - mean) / mean * 100) : null;
    bullet = `<div class="stat-bullet" role="img" aria-label="${fmtN(v)} versus ${meanLabel} ${fmtN(Math.round(mean))}">
      <span class="ch-track"><span class="ch-bar" style="width:${(v / max * 100).toFixed(1)}%;background:${color}"></span>
        <span class="stat-mean" style="left:${(mean / max * 100).toFixed(1)}%"></span></span>
      <span class="stat-delta">${dpct == null ? '' : `${dpct >= 0 ? '+' : ''}${dpct.toFixed(0)}% vs ${esc(meanLabel)} ${fmtN(Math.round(mean))}`}</span>
    </div>`;
  }
  return `<div class="stat"><div class="k">${k}</div><div class="v">${typeof v === 'number' ? fmtN(v) : v}</div><div class="d">${d}</div>${bullet}</div>`;
}

// ---- shared tick helpers
export function tickLabel(t) {
  const a = Math.abs(t);
  if (a >= 100 || a === 0) return String(Math.round(t));
  if (a >= 1) return String(+t.toFixed(1));
  return t.toPrecision(2);
}
export function niceTicks(lo, hi, n) {
  const span = hi - lo || 1;
  const step0 = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 5, 10].map(m => m * mag).find(s => span / s <= n) || 10 * mag;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(v);
  return out;
}

// Width helper: fit a chart to its host, floor 300, cap `max`.
export function fitWidth(el, max = 640) {
  const w = el ? el.clientWidth : 0;
  return Math.max(300, Math.min(w || max, max));
}
