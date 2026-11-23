/* app.js
 * LFQ Explorer: client-side differential proteomics for Proteome Discoverer exports.
 * Everything runs in memory in this tab. No network requests are made with user data,
 * and nothing is written to cookies, localStorage, IndexedDB or any server.
 */
'use strict';

(() => {
  /* ===================================================================
   * Helpers
   * =================================================================== */
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));
  const tick = () => new Promise(res => setTimeout(res, 20));
  const isNum = v => typeof v === 'number' && Number.isFinite(v);
  const norm = s => String(s ?? '').toLowerCase().replace(/\s+/g, ' ').trim();
  const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const safeName = s => String(s).replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, '') || 'file';
  const unique = arr => [...new Set(arr)];

  function toNum(v) {
    if (typeof v === 'number') return v;
    if (v == null) return NaN;
    const s = String(v).trim();
    return s === '' ? NaN : Number(s);
  }
  const fmt = (v, d = 2) => (isNum(v) ? v.toFixed(d) : '–');
  function fmtP(v) {
    if (!isNum(v)) return '–';
    if (v === 0) return '< 1e-300';
    return v < 1e-3 ? v.toExponential(2) : v.toPrecision(3);
  }

  const FONT = '"IBM Plex Sans", system-ui, sans-serif';
  const PALETTE = ['#0B8A99', '#C0392B', '#7A5BB5', '#C98A12', '#2F72C0', '#3E8E4E', '#B04A86', '#6F8187', '#8C5A2B', '#1F9A8A'];
  const SYMBOLS = ['circle', 'square', 'diamond', 'triangle-up', 'cross', 'star', 'hexagon', 'triangle-down'];

  /* ===================================================================
   * State (in memory only)
   * =================================================================== */
  const state = {
    wb: null,
    fileName: '',
    parsed: null,     // header, data rows, detected columns
    samples: [],      // editable sample design
    res: null,        // analysis results
    sel: null,        // selected protein index into res.proteins
    activeTab: null,
    de: { comp: 0, swap: false, page: 0, sortKey: 'p', sortDir: 1, filter: '', sigOnly: false }
  };

  /* ===================================================================
   * UI utilities
   * =================================================================== */
  function notice(msg, kind = 'error') {
    const el = $('#notice');
    if (!msg) { el.hidden = true; return; }
    el.className = 'notice' + (kind === 'error' ? '' : ' notice-' + kind);
    el.textContent = msg;
    el.hidden = false;
  }
  function setStatus(msg, busy = false) {
    const el = $('#run-status');
    el.textContent = msg || '';
    el.classList.toggle('is-busy', !!busy);
  }
  function markStale() { if (state.res) $('#stale').hidden = false; }

  function theme() {
    const cs = getComputedStyle(document.documentElement);
    const v = n => cs.getPropertyValue(n).trim();
    return { ink: v('--ink'), muted: v('--ink-muted'), rule: v('--rule'), up: v('--up'), down: v('--down'), ns: v('--ns'), accent: v('--accent'), mid: v('--heat-mid'), surface: v('--surface') };
  }

  function merge(a, b) {
    for (const k of Object.keys(b)) {
      const bv = b[k], av = a[k];
      if (bv && typeof bv === 'object' && !Array.isArray(bv) && av && typeof av === 'object' && !Array.isArray(av)) merge(av, bv);
      else a[k] = bv;
    }
    return a;
  }

  function baseLayout(extra = {}) {
    const t = theme();
    const axis = () => ({ gridcolor: t.rule, linecolor: t.rule, zerolinecolor: t.rule, tickfont: { color: t.muted }, title: { font: { color: t.ink, size: 13 } }, automargin: true });
    return merge({
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: FONT, size: 12, color: t.ink },
      margin: { l: 64, r: 24, t: 44, b: 56 },
      hoverlabel: { font: { family: FONT, size: 12 } },
      legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom', font: { size: 12 } },
      xaxis: axis(),
      yaxis: axis()
    }, extra);
  }

  function cfg(name) {
    return {
      responsive: true,
      displaylogo: false,
      modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines'],
      toImageButtonOptions: { format: 'png', filename: safeName(name), scale: 3 },
      modeBarButtonsToAdd: [{
        name: 'Download plot as SVG',
        icon: Plotly.Icons.disk,
        click: gd => Plotly.downloadImage(gd, { format: 'svg', filename: safeName(name) })
      }]
    };
  }

  function onClick(gd, handler) {
    if (gd.removeAllListeners) gd.removeAllListeners('plotly_click');
    gd.on('plotly_click', handler);
  }

  function groupColor(g) {
    const groups = state.res ? state.res.groups : unique(state.samples.filter(s => s.include).map(s => s.group));
    const i = groups.indexOf(g);
    return PALETTE[(i < 0 ? 0 : i) % PALETTE.length];
  }

  /* ===================================================================
   * Hero illustration (decorative)
   * =================================================================== */
  function drawHero() {
    const svg = $('#hero-volcano');
    if (!svg) return;
    let seed = 11;
    const rnd = () => {
      seed = (seed + 0x6D2B79F5) | 0;
      let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    const gauss = () => { let u = 0, v = 0; while (!u) u = rnd(); while (!v) v = rnd(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
    const cx0 = 200, base = 240, sx = 56, sy = 40;
    let dots = '';
    for (let i = 0; i < 560; i++) {
      const x = gauss() * 1.05;
      const y = Math.abs(x) * (0.5 + rnd() * 1.7) + rnd() * rnd() * 1.0;
      const sig = Math.abs(x) > 1 && y > 1.3;
      const cls = sig ? (x > 0 ? 'h-up' : 'h-down') : 'h-ns';
      const cx = cx0 + x * sx, cy = base - y * sy;
      if (cx < 8 || cx > 392 || cy < 8) continue;
      dots += `<circle class="${cls}" cx="${cx.toFixed(1)}" cy="${cy.toFixed(1)}" r="${sig ? 3.3 : 2.5}"/>`;
    }
    const thrY = base - 1.3 * sy;
    svg.innerHTML =
      `<line class="thr" x1="${cx0 - sx}" y1="8" x2="${cx0 - sx}" y2="${base}"/>` +
      `<line class="thr" x1="${cx0 + sx}" y1="8" x2="${cx0 + sx}" y2="${base}"/>` +
      `<line class="thr" x1="8" y1="${thrY}" x2="392" y2="${thrY}"/>` +
      `<line class="ax" x1="8" y1="${base + 6}" x2="392" y2="${base + 6}"/>` + dots;
  }

  /* ===================================================================
   * File loading and parsing
   * =================================================================== */
  const QUANT_TYPES = [
    { label: 'Normalized abundances', re: /Abundances?\s*\(Normalized\)\s*:/i },
    { label: 'Scaled abundances', re: /Abundances?\s*\(Scaled\)\s*:/i },
    { label: 'Raw abundances', re: /^Abundance\s*:/i }
  ];

  function findCol(header, names) {
    const H = header.map(norm);
    for (const n of names) {
      const i = H.indexOf(norm(n));
      if (i >= 0) return i;
    }
    return -1;
  }

  function detectColumns(header) {
    let fdr = header.findIndex(h => /^protein fdr confidence:\s*combined$/i.test(h.trim()));
    if (fdr < 0) fdr = header.findIndex(h => /protein fdr confidence/i.test(h));
    return {
      accession: findCol(header, ['Accession', 'Protein Accession', 'Master Protein Accessions']),
      gene: findCol(header, ['Gene Symbol', 'Gene Symbols', 'Gene Name', 'Gene Names', 'Gene']),
      description: findCol(header, ['Description', 'Protein Description']),
      master: findCol(header, ['Master']),
      contaminant: findCol(header, ['Contaminant']),
      fdr,
      uniquePep: findCol(header, ['# Unique Peptides', 'Unique Peptides', '# Unique peptides'])
    };
  }

  function detectQuant(header) {
    return QUANT_TYPES
      .map(t => ({ label: t.label, idx: header.map((h, i) => (t.re.test(h) ? i : -1)).filter(i => i >= 0) }))
      .filter(t => t.idx.length);
  }

  function headerRowOf(ws) {
    if (!ws || !ws['!ref']) return -1;
    const r = XLSX.utils.decode_range(ws['!ref']);
    for (let R = r.s.r; R <= Math.min(r.e.r, r.s.r + 9); R++) {
      for (let C = r.s.c; C <= r.e.c; C++) {
        const cell = ws[XLSX.utils.encode_cell({ r: R, c: C })];
        if (cell && norm(cell.v) === 'accession') return R - r.s.r;
      }
    }
    return -1;
  }

  async function loadFile(file) {
    if (!file) return;
    notice(null);
    try {
      await needLibs('XLSX', 'Plotly');
      setStatus(`Reading ${file.name}…`, true);
      await tick();
      const lower = file.name.toLowerCase();
      let wb;
      if (/\.(csv|tsv|txt)$/.test(lower)) wb = XLSX.read(await file.text(), { type: 'string' });
      else wb = XLSX.read(new Uint8Array(await file.arrayBuffer()), { type: 'array' });

      resetResults();
      state.wb = wb;
      state.fileName = file.name;

      const sel = $('#sheet-select');
      sel.innerHTML = wb.SheetNames.map(n => `<option>${esc(n)}</option>`).join('');
      $('#sheet-row').hidden = wb.SheetNames.length < 2;
      const pick = wb.SheetNames.find(n => headerRowOf(wb.Sheets[n]) >= 0) || wb.SheetNames[0];
      sel.value = pick;
      parseSheet(pick);
      setStatus('File loaded. Check the samples, then run the analysis.');
    } catch (e) {
      console.error(e);
      setStatus('');
      notice(e.message || 'This file could not be read.');
    }
  }

  function parseSheet(name) {
    const ws = state.wb.Sheets[name];
    const h = Math.max(0, headerRowOf(ws));
    const rows = XLSX.utils.sheet_to_json(ws, { header: 1, raw: true, defval: null });
    const header = (rows[h] || []).map(c => (c == null ? '' : String(c).trim()));
    const data = rows.slice(h + 1).filter(r => r && r.some(v => v !== null && v !== ''));
    const cols = detectColumns(header);
    if (cols.accession < 0) throw new Error('No “Accession” column was found. Export the Proteins table from Proteome Discoverer with the Accession column included.');
    const quant = detectQuant(header);
    if (!quant.length) throw new Error('No per-sample abundance columns were found. Expected headers such as “Abundances (Normalized): F3: Sample, 11E8_KO”.');

    state.parsed = { sheet: name, header, data, cols, quant, quantIdx: 0 };

    const qs = $('#quant-select');
    qs.innerHTML = quant.map((q, i) => `<option value="${i}">${esc(q.label)} (${q.idx.length} columns)</option>`).join('');
    $('#quant-row').hidden = quant.length < 2;

    buildSamples();
    renderFileInfo();
    renderFilterAvailability();
    renderSampleTable();

    $('#intro').hidden = true;
    $('#tabs').hidden = false;
    $('#open-samples').disabled = false;
    $('#run-btn').disabled = false;
    setTabsEnabled(false);
    showTab('samples');
  }

  function buildSamples() {
    const { header, quant, quantIdx } = state.parsed;
    state.samples = quant[quantIdx].idx.map((ci, k) => {
      const h = header[ci];
      const m = h.match(/:\s*([^:]+):\s*Sample,\s*(.+)$/i);
      let sample, group;
      if (m) { sample = m[1].trim(); group = m[2].trim(); }
      else {
        const parts = h.split(':');
        sample = (parts[1] || `S${k + 1}`).trim();
        group = '';
      }
      return { col: ci, header: h, sample, group, batch: '', include: true };
    });
  }

  function renderFileInfo() {
    const P = state.parsed;
    const groups = unique(state.samples.map(s => s.group).filter(Boolean));
    const el = $('#file-info');
    el.innerHTML =
      `<b>${esc(state.fileName)}</b>` +
      `<span>${P.data.length.toLocaleString()} protein rows, ${state.samples.length} samples, ${groups.length} groups</span>`;
    el.hidden = false;
  }

  function renderFilterAvailability() {
    const c = state.parsed.cols;
    const missing = [];
    const setAvail = (id, ok, label) => {
      const input = $(id);
      input.disabled = !ok;
      input.closest('label').classList.toggle('is-disabled', !ok);
      if (!ok) missing.push(label);
    };
    setAvail('#f-master', c.master >= 0, 'Master');
    setAvail('#f-contam', c.contaminant >= 0, 'Contaminant');
    setAvail('#f-highconf', c.fdr >= 0, 'Protein FDR Confidence');
    setAvail('#f-minpep', c.uniquePep >= 0, '# Unique Peptides');
    const note = $('#filter-notes');
    note.hidden = !missing.length;
    note.textContent = missing.length ? `Not in this file, so skipped: ${missing.join(', ')}.` : '';
  }

  /* ===================================================================
   * Sample table
   * =================================================================== */
  function renderSampleTable() {
    const batchOn = $('#b-enable').checked;
    $('#sample-tbody').innerHTML = state.samples.map((s, k) => `
      <tr data-k="${k}" class="${s.include ? '' : 'is-off'}">
        <td><input type="checkbox" data-f="include" ${s.include ? 'checked' : ''} aria-label="Include ${esc(s.sample)}"></td>
        <td><input type="text" data-f="sample" value="${esc(s.sample)}" aria-label="Sample name"></td>
        <td><input type="text" data-f="group" value="${esc(s.group)}" list="group-options" aria-label="Group for ${esc(s.sample)}"></td>
        <td><input type="text" data-f="batch" value="${esc(s.batch)}" ${batchOn ? '' : 'disabled'} placeholder="${batchOn ? 'e.g. 1' : 'Off'}" aria-label="Batch for ${esc(s.sample)}"></td>
        <td class="muted small">${esc(s.header)}</td>
      </tr>`).join('');
    renderGroupSummary();
  }

  function renderGroupSummary() {
    const inc = state.samples.filter(s => s.include);
    const groups = unique(inc.map(s => s.group.trim()).filter(Boolean));
    const minPer = Math.max(2, +$('#s-minper').value || 2);
    let html = groups.map(g => {
      const n = inc.filter(s => s.group.trim() === g).length;
      const warn = n < minPer;
      return `<span class="chip${warn ? ' is-warn' : ''}" title="${warn ? 'Too few samples to test' : ''}"><span class="dot" style="background:${PALETTE[groups.indexOf(g) % PALETTE.length]}"></span>${esc(g)}: ${n}</span>`;
    }).join('');
    const ungrouped = inc.filter(s => !s.group.trim()).length;
    if (ungrouped) html += `<span class="chip is-warn">${ungrouped} included sample${ungrouped > 1 ? 's' : ''} without a group</span>`;
    if ($('#b-enable').checked) {
      const batches = unique(inc.map(s => s.batch.trim()).filter(Boolean));
      const noBatch = inc.filter(s => !s.batch.trim()).length;
      html += batches.map(b => `<span class="chip">Batch ${esc(b)}: ${inc.filter(s => s.batch.trim() === b).length}</span>`).join('');
      if (noBatch) html += `<span class="chip is-warn">${noBatch} sample${noBatch > 1 ? 's' : ''} without a batch</span>`;
    }
    $('#group-summary').innerHTML = html;
    $('#group-options').innerHTML = groups.map(g => `<option value="${esc(g)}">`).join('');
  }

  function onSampleEdit(e) {
    const input = e.target;
    const tr = input.closest('tr');
    if (!tr || !input.dataset.f) return;
    const s = state.samples[+tr.dataset.k];
    if (input.dataset.f === 'include') {
      s.include = input.checked;
      tr.classList.toggle('is-off', !input.checked);
    } else {
      s[input.dataset.f] = input.value;
    }
    renderGroupSummary();
    markStale();
  }

  async function loadBatchFile(file) {
    const msg = $('#b-file-msg');
    try {
      await needLibs('XLSX');
      setStatus('');
      const lower = file.name.toLowerCase();
      const wb = /\.(csv|tsv|txt)$/.test(lower)
        ? XLSX.read(await file.text(), { type: 'string' })
        : XLSX.read(new Uint8Array(await file.arrayBuffer()), { type: 'array' });
      const rows = XLSX.utils.sheet_to_json(wb.Sheets[wb.SheetNames[0]], { header: 1, raw: false, defval: '' });
      const head = (rows[0] || []).map(norm);
      const si = head.indexOf('sample'), bi = head.indexOf('batch');
      if (si < 0 || bi < 0) throw new Error('The batch file needs columns named Sample and Batch.');
      const map = new Map(rows.slice(1).map(r => [String(r[si]).trim(), String(r[bi]).trim()]));
      let matched = 0;
      state.samples.forEach(s => { if (map.has(s.sample.trim())) { s.batch = map.get(s.sample.trim()); matched++; } });
      renderSampleTable();
      markStale();
      msg.textContent = `Batch labels set for ${matched} of ${state.samples.length} samples.`;
      msg.classList.toggle('warn-text', matched < state.samples.length);
      msg.hidden = false;
    } catch (e) {
      msg.textContent = e.message;
      msg.classList.add('warn-text');
      msg.hidden = false;
    }
  }

  /* ===================================================================
   * Analysis pipeline
   * =================================================================== */
  function readSettings() {
    const on = id => $(id).checked && !$(id).disabled;
    return {
      masterOnly: on('#f-master'),
      removeContam: on('#f-contam'),
      highConf: on('#f-highconf'),
      minPep: $('#f-minpep').disabled ? 0 : Math.max(0, +$('#f-minpep').value || 0),
      minRep: Math.max(1, Math.round(+$('#f-minrep').value || 1)),
      batch: $('#b-enable').checked,
      protect: $('#b-protect').checked,
      remask: $('#b-remask').checked,
      minPer: Math.max(2, Math.round(+$('#s-minper').value || 2)),
      quant: state.parsed.quant[state.parsed.quantIdx].label
    };
  }

  function geneOf(r, c) {
    let g = '';
    if (c.gene >= 0 && r[c.gene] != null) g = String(r[c.gene]).trim();
    if (!g && c.description >= 0) {
      const m = String(r[c.description] ?? '').match(/\bGN=([^\s;]+)/);
      if (m) g = m[1];
    }
    return g;
  }

  function buildMatrix(S, inc) {
    const { data, cols } = state.parsed;
    const counts = [];
    let rows = data.filter(r => String(r[cols.accession] ?? '').trim() !== '');
    counts.push(['Protein rows in file', rows.length]);

    if (S.masterOnly) {
      rows = rows.filter(r => norm(r[cols.master]) === 'master protein');
      counts.push(['Master proteins', rows.length]);
    }
    if (S.removeContam) {
      rows = rows.filter(r => !['true', '1', 'yes'].includes(norm(r[cols.contaminant])));
      counts.push(['Contaminants removed', rows.length]);
    }
    if (S.highConf) {
      rows = rows.filter(r => norm(r[cols.fdr]) === 'high');
      counts.push(['High FDR confidence', rows.length]);
    }
    if (S.minPep > 0) {
      rows = rows.filter(r => { const v = toNum(r[cols.uniquePep]); return isNum(v) && v >= S.minPep; });
      counts.push([`At least ${S.minPep} unique peptides`, rows.length]);
    }

    // Collapse duplicated accessions by median, mask non-positive values, log2
    const byAcc = new Map();
    for (const r of rows) {
      const acc = String(r[cols.accession]).trim();
      let e = byAcc.get(acc);
      if (!e) {
        e = { acc, gene: geneOf(r, cols), desc: cols.description >= 0 ? String(r[cols.description] ?? '').trim() : '', vals: inc.map(() => []) };
        byAcc.set(acc, e);
      }
      inc.forEach((s, j) => { const v = toNum(r[s.col]); if (isNum(v)) e.vals[j].push(v); });
    }
    const proteins = [];
    for (const e of byAcc.values()) {
      const raw = new Float64Array(inc.length);
      e.vals.forEach((arr, j) => {
        const m = arr.length ? Stats.median(arr) : NaN;
        raw[j] = m > 0 ? Math.log2(m) : NaN;
      });
      proteins.push({ acc: e.acc, gene: e.gene, desc: e.desc, raw });
    }
    if (byAcc.size !== rows.length) counts.push(['Unique accessions', proteins.length]);
    return { proteins, counts };
  }

  function computePCA(M, n) {
    if (n < 3) return null;
    const keep = [];
    for (let i = 0; i < M.length; i++) {
      const seen = new Set();
      for (const v of M[i]) if (isNum(v)) { seen.add(v); if (seen.size > 1) break; }
      if (seen.size > 1) keep.push(i);
    }
    if (keep.length < 2) return null;
    const X = Array.from({ length: n }, () => new Float64Array(keep.length));
    keep.forEach((pi, c) => {
      const row = M[pi];
      const med = Stats.median(row);
      for (let j = 0; j < n; j++) X[j][c] = isNum(row[j]) ? row[j] : med;
    });
    Stats.standardizeColumns(X);
    const pc = Stats.pca(X, Math.min(n - 1, 8));
    return { scores: pc.scores, explained: pc.explained, corr: Stats.corrRows(X), nProteins: keep.length };
  }

  function differential(M, groups, groupIdx, minPer) {
    const comps = [];
    for (let a = 0; a < groups.length; a++) {
      for (let b = a + 1; b < groups.length; b++) {
        const i1 = groupIdx[a], i2 = groupIdx[b];
        if (i1.length < minPer || i2.length < minPer) continue;
        const results = [];
        for (let p = 0; p < M.length; p++) {
          const row = M[p];
          const v1 = i1.map(j => row[j]).filter(isNum);
          const v2 = i2.map(j => row[j]).filter(isNum);
          if (v1.length < minPer || v2.length < minPer) continue;
          const w = Stats.welch(v1, v2);
          const m1 = Stats.mean(v1), m2 = Stats.mean(v2);
          results.push({ i: p, n1: v1.length, n2: v2.length, m1, m2, lfc: m1 - m2, t: w.t, p: w.p, adj: NaN });
        }
        const adj = Stats.bh(results.map(r => r.p));
        results.forEach((r, k) => { r.adj = adj[k]; });
        comps.push({ g1: groups[a], g2: groups[b], results, byProtein: new Map(results.map((r, k) => [r.i, k])) });
      }
    }
    return comps;
  }

  async function runAnalysis() {
    if (!state.parsed) return;
    notice(null);
    const S = readSettings();
    const inc = state.samples.filter(s => s.include).map(s => ({ ...s, sample: s.sample.trim(), group: s.group.trim(), batch: s.batch.trim() }));

    if (inc.length < 2) return notice('Include at least two samples.');
    const noGroup = inc.filter(s => !s.group);
    if (noGroup.length) return notice(`Assign a group to every included sample. Missing for: ${noGroup.map(s => s.sample).join(', ')}.`);
    const dupNames = inc.map(s => s.sample).filter((v, i, a) => a.indexOf(v) !== i);
    if (dupNames.length) return notice(`Sample names must be unique. Repeated: ${unique(dupNames).join(', ')}.`);
    const groups = unique(inc.map(s => s.group));
    if (groups.length < 2) return notice('At least two groups are needed for differential analysis.');
    const groupIdx = groups.map(g => inc.map((s, j) => (s.group === g ? j : -1)).filter(j => j >= 0));
    const testable = groupIdx.filter(ix => ix.length >= S.minPer).length;
    if (testable < 2) return notice(`At least two groups need ${S.minPer} or more samples to run the t-tests.`);
    if (S.batch) {
      const noBatch = inc.filter(s => !s.batch);
      if (noBatch.length) return notice(`Batch correction is on, but these samples have no batch label: ${noBatch.map(s => s.sample).join(', ')}. Enter them in the sample table or load a batch CSV.`);
    }

    $('#run-btn').disabled = true;
    try {
      setStatus('Filtering proteins…', true); await tick();
      const { proteins, counts } = buildMatrix(S, inc);

      setStatus('Applying the detection filter…', true); await tick();
      const kept = proteins.filter(p => groupIdx.some(ix => {
        const need = Math.min(S.minRep, ix.length);
        let c = 0;
        for (const j of ix) if (isNum(p.raw[j])) c++;
        return c > 0 && c >= need;
      }));
      counts.push([`Detected in ≥ ${S.minRep} replicates of a group`, kept.length]);
      if (kept.length < 3) throw new Error('Fewer than three proteins passed the filters. Loosen the filters or check the sample selection.');

      const raw = kept.map(p => p.raw);
      let corrected = null;
      if (S.batch) {
        setStatus('Running ComBat batch correction…', true); await tick();
        const Y = raw.map(row => { const med = Stats.median(row); return Float64Array.from(row, v => (isNum(v) ? v : med)); });
        corrected = Stats.combat(Y, inc.map(s => s.batch), S.protect ? inc.map(s => s.group) : null);
        if (S.remask) corrected.forEach((row, i) => { for (let j = 0; j < row.length; j++) if (!isNum(raw[i][j])) row[j] = NaN; });
      }
      const M = corrected || raw;

      setStatus('Computing PCA and correlations…', true); await tick();
      const pca = computePCA(M, inc.length);
      const pcaRaw = corrected ? computePCA(raw, inc.length) : null;

      setStatus('Testing differential abundance…', true); await tick();
      const comps = differential(M, groups, groupIdx, S.minPer);
      if (!comps.length) throw new Error('No comparison had enough samples to test.');

      state.res = {
        S, samples: inc, groups, groupIdx, counts,
        proteins: kept.map(p => ({ acc: p.acc, gene: p.gene, desc: p.desc })),
        raw, corrected, M, pca, pcaRaw, comps, ranAt: new Date()
      };
      state.sel = null;
      Object.assign(state.de, { comp: 0, swap: false, page: 0 });

      populateResultControls();
      setTabsEnabled(true);
      $('#stale').hidden = true;
      setStatus(`Done. ${kept.length.toLocaleString()} proteins, ${comps.length} comparison${comps.length > 1 ? 's' : ''}.`);
      showTab('overview');
    } catch (e) {
      console.error(e);
      setStatus('');
      notice(e.message || 'The analysis stopped because of an unexpected error.');
    } finally {
      $('#run-btn').disabled = false;
    }
  }

  function resetResults() {
    state.res = null;
    state.sel = null;
    $('#stale').hidden = true;
    ['plot-missing', 'plot-dist', 'plot-pca', 'plot-scree', 'plot-corr', 'plot-heatmap', 'plot-volcano'].forEach(id => {
      const el = document.getElementById(id);
      if (el && window.Plotly) Plotly.purge(el);
    });
  }

  /* ===================================================================
   * Tabs
   * =================================================================== */
  function setTabsEnabled(on) {
    $$('#tabs button').forEach(b => { if (b.dataset.tab !== 'samples') b.disabled = !on; });
  }

  function showTab(name) {
    state.activeTab = name;
    $$('#tabs button').forEach(b => b.setAttribute('aria-selected', String(b.dataset.tab === name)));
    $$('.panel').forEach(p => { p.hidden = p.id !== `panel-${name}`; });
    renderActive();
  }

  function renderActive() {
    const name = state.activeTab;
    if (name === 'samples') return;
    if (!state.res) return;
    try {
      ({ overview: renderOverview, pca: renderPCA, corr: renderCorrelation, heatmap: renderHeatmap, de: renderDE, downloads: renderDownloads })[name]?.();
    } catch (e) {
      console.error(e);
      notice('A plot could not be drawn: ' + e.message);
    }
  }

  function populateResultControls() {
    const R = state.res;
    const k = R.pca ? R.pca.explained.length : 0;
    const opts = [...Array(k).keys()].map(i => `<option value="${i}">PC${i + 1}</option>`).join('');
    $('#pca-x').innerHTML = opts;
    $('#pca-y').innerHTML = opts;
    $('#pca-x').value = '0';
    $('#pca-y').value = k > 1 ? '1' : '0';
    $('#pca-before-wrap').hidden = !R.pcaRaw;
    $('#pca-before').checked = false;
    $('#ov-corr-wrap').hidden = !R.corrected;

    $('#de-comp').innerHTML = R.comps.map((c, i) => `<option value="${i}">${esc(c.g1)} vs ${esc(c.g2)}</option>`).join('');
    $('#de-comp').value = '0';
    $('#protein-list').innerHTML = R.proteins.map(p => `<option value="${esc(p.gene || p.acc)}">${esc(p.acc)}</option>`).join('');
  }

  /* ===================================================================
   * Overview
   * =================================================================== */
  function renderOverview() {
    const R = state.res;
    const t = theme();
    $('#ov-counts').innerHTML = R.counts.map(([k, v]) => `<li><span>${esc(k)}</span><b>${v.toLocaleString()}</b></li>`).join('');
    const design = R.groups.map((g, i) => `${g} (${R.groupIdx[i].length})`).join(', ');
    $('#ov-design').textContent = `${R.samples.length} samples in ${R.groups.length} groups: ${design}.` +
      (R.corrected ? ' ComBat batch correction applied.' : '');

    const names = R.samples.map(s => s.sample);
    const colors = R.samples.map(s => groupColor(s.group));
    const miss = names.map((_, j) => R.raw.reduce((c, row) => c + (isNum(row[j]) ? 0 : 1), 0));
    Plotly.react('plot-missing', [{
      type: 'bar', x: names, y: miss, marker: { color: colors },
      customdata: R.samples.map(s => s.group),
      hovertemplate: '%{x} (%{customdata})<br>%{y} missing<extra></extra>'
    }], baseLayout({
      margin: { t: 16, b: 70 },
      xaxis: { title: { text: '' }, tickangle: -45, type: 'category' },
      yaxis: { title: { text: 'Proteins without a value' } },
      bargap: 0.25
    }), cfg('missing_values_per_sample'));

    const useCorr = R.corrected && $('#ov-corr').checked;
    const mat = useCorr ? R.corrected : R.raw;
    const traces = R.samples.map((s, j) => ({
      type: 'box', name: s.sample, y: mat.map(row => row[j]).filter(isNum),
      marker: { color: groupColor(s.group), size: 2 }, line: { width: 1.2 },
      boxpoints: false, legendgroup: s.group, showlegend: false,
      hovertemplate: `${esc(s.sample)} (${esc(s.group)})<extra></extra>`
    }));
    Plotly.react('plot-dist', traces, baseLayout({
      margin: { t: 16, b: 70 },
      xaxis: { tickangle: -45, type: 'category' },
      yaxis: { title: { text: useCorr ? 'log2 abundance (batch corrected)' : 'log2 abundance' } }
    }), cfg('abundance_distribution'));
  }

  /* ===================================================================
   * PCA
   * =================================================================== */
  function renderPCA() {
    const R = state.res;
    const before = R.pcaRaw && $('#pca-before').checked;
    const P = before ? R.pcaRaw : R.pca;
    const note = $('#pca-note');
    if (!P) {
      note.textContent = 'PCA needs at least three samples and two variable proteins.';
      Plotly.purge('plot-pca'); Plotly.purge('plot-scree');
      return;
    }
    const cx = +$('#pca-x').value, cy = +$('#pca-y').value;
    const batches = R.S.batch ? unique(R.samples.map(s => s.batch)) : null;
    note.textContent = `${P.nProteins.toLocaleString()} proteins, median-imputed and standardized.` +
      (batches ? ' Marker shape shows the batch.' : '') +
      (before ? ' Showing values before batch correction.' : '');

    const traces = R.groups.map((g, gi) => {
      const idx = R.groupIdx[gi];
      return {
        type: 'scatter', mode: 'markers+text', name: g,
        x: idx.map(j => P.scores[j][cx]), y: idx.map(j => P.scores[j][cy]),
        text: idx.map(j => R.samples[j].sample), textposition: 'top center',
        textfont: { size: 11, color: theme().muted },
        customdata: idx.map(j => R.samples[j].batch || '–'),
        marker: {
          size: 14, color: groupColor(g), line: { width: 1, color: theme().surface },
          symbol: batches ? idx.map(j => SYMBOLS[batches.indexOf(R.samples[j].batch) % SYMBOLS.length]) : 'circle'
        },
        hovertemplate: `<b>%{text}</b><br>${esc(g)}` + (batches ? '<br>Batch %{customdata}' : '') + '<extra></extra>'
      };
    });
    Plotly.react('plot-pca', traces, baseLayout({
      xaxis: { title: { text: `PC${cx + 1} (${(P.explained[cx] * 100).toFixed(1)}%)` } },
      yaxis: { title: { text: `PC${cy + 1} (${(P.explained[cy] * 100).toFixed(1)}%)` } }
    }), cfg('PCA'));

    Plotly.react('plot-scree', [{
      type: 'bar', x: P.explained.map((_, i) => `PC${i + 1}`), y: P.explained.map(v => v * 100),
      marker: { color: P.explained.map((_, i) => (i === cx || i === cy ? theme().accent : theme().ns)) },
      hovertemplate: '%{x}: %{y:.1f}%<extra></extra>'
    }], baseLayout({
      title: { text: 'Variance explained', font: { size: 13 }, x: 0, xanchor: 'left' },
      margin: { t: 40 },
      yaxis: { title: { text: '% of variance' } }
    }), cfg('PCA_scree'));
  }

  /* ===================================================================
   * Clustered heatmap with dendrograms
   * =================================================================== */
  function clusteredHeatmap(elId, o) {
    const t = theme();
    const nr = o.z.length, nc = o.z[0].length;
    const rc = o.rowCluster && nr > 2 ? Stats.hclust(o.z) : null;
    const cc = o.colCluster && nc > 2 ? Stats.hclust(o.z[0].map((_, j) => o.z.map(r => r[j]))) : null;
    const ro = rc ? rc.order : [...Array(nr).keys()];
    const co = cc ? cc.order : [...Array(nc).keys()];
    const z = ro.map(i => co.map(j => o.z[i][j]));
    const xs = [...Array(nc).keys()], ys = [...Array(nr).keys()];
    const hover = ro.map(i => co.map(j => o.hover(i, j)));

    const colorbarCorner = !!rc;
    const traces = [{
      type: 'heatmap', z, x: xs, y: ys, xaxis: 'x', yaxis: 'y',
      colorscale: o.colorscale, zmid: o.zmid, zmin: o.zmin, zmax: o.zmax,
      text: hover, hovertemplate: '%{text}<extra></extra>',
      xgap: o.gap ? 1 : 0, ygap: o.gap ? 1 : 0,
      texttemplate: o.annotate ? '%{z:.2f}' : '', textfont: { size: 10 },
      colorbar: colorbarCorner
        ? { title: { text: o.cbTitle, side: 'top', font: { size: 11 } }, x: 0, xanchor: 'left', y: 1, yanchor: 'top', len: 0.15, thickness: 10, tickfont: { size: 10, color: t.muted }, outlinewidth: 0 }
        : { title: { text: o.cbTitle, font: { size: 11 } }, thickness: 12, len: 0.5, tickfont: { size: 10, color: t.muted }, outlinewidth: 0 }
    }];

    const longestRow = Math.max(...o.rowLabels.map(l => l.length));
    const longestCol = Math.max(...o.colLabels.map(l => l.length));
    const rowLabelSpace = Math.min(240, 24 + longestRow * 7);
    const layout = baseLayout({
      showlegend: false,
      height: o.height,
      margin: { t: 10, b: Math.min(200, 40 + longestCol * 6.5), l: rc ? 10 : rowLabelSpace, r: rc ? rowLabelSpace : 90 },
      xaxis: {
        domain: [rc ? 0.14 : 0, 1], tickvals: xs, ticktext: co.map(j => o.colTickText(j)),
        tickangle: -50, showgrid: false, zeroline: false, showline: false, range: [-0.5, nc - 0.5], anchor: 'y', automargin: false
      },
      yaxis: {
        domain: [0, cc ? 0.84 : 1], tickvals: ys, ticktext: ro.map(i => o.rowTickText(i)),
        side: rc ? 'right' : 'left', showgrid: false, zeroline: false, showline: false,
        range: [nr - 0.5, -0.5], anchor: 'x', automargin: false, tickfont: { size: nr > 120 ? 8 : nr > 60 ? 9 : 11 }
      }
    });

    if (cc) {
      traces.push({ type: 'scatter', mode: 'lines', x: cc.seg.pos, y: cc.seg.h, xaxis: 'x', yaxis: 'y2', line: { color: t.muted, width: 1 }, hoverinfo: 'skip' });
      layout.yaxis2 = { domain: [0.855, 1], anchor: 'x', visible: false, range: [0, cc.maxH * 1.05 || 1] };
    }
    if (rc) {
      traces.push({ type: 'scatter', mode: 'lines', x: rc.seg.h, y: rc.seg.pos, xaxis: 'x2', yaxis: 'y', line: { color: t.muted, width: 1 }, hoverinfo: 'skip' });
      layout.xaxis2 = { domain: [0, 0.13], anchor: 'y', visible: false, range: [rc.maxH * 1.05 || 1, 0] };
    }

    const gd = document.getElementById(elId);
    return Plotly.react(gd, traces, layout, cfg(o.fileName)).then(() => {
      onClick(gd, ev => {
        const p = ev.points && ev.points[0];
        if (!p || p.data.type !== 'heatmap' || !o.onRowClick) return;
        o.onRowClick(ro[p.y]);
      });
    });
  }

  const dot = c => `<span style="color:${c}">●</span>`;

  function divergingScale() {
    const t = theme();
    return [[0, t.down], [0.5, t.mid], [1, t.up]];
  }

  /* ===================================================================
   * Sample correlation
   * =================================================================== */
  function renderCorrelation() {
    const R = state.res;
    if (!R.pca) { Plotly.purge('plot-corr'); return; }
    const C = R.pca.corr;
    const names = R.samples.map(s => s.sample);
    const n = names.length;
    clusteredHeatmap('plot-corr', {
      z: C, rowCluster: true, colCluster: true, zmid: 0, gap: true,
      colorscale: divergingScale(), cbTitle: 'r',
      annotate: $('#corr-annot').checked && n <= 30,
      rowLabels: names, colLabels: names,
      rowTickText: i => `${esc(names[i])} ${dot(groupColor(R.samples[i].group))}`,
      colTickText: j => `${dot(groupColor(R.samples[j].group))} ${esc(names[j])}`,
      hover: (i, j) => `${esc(names[i])} (${esc(R.samples[i].group)})<br>${esc(names[j])} (${esc(R.samples[j].group)})<br>r = ${C[i][j].toFixed(3)}`,
      height: Math.max(560, Math.min(1000, 260 + n * 36)),
      fileName: 'sample_correlation_clustering'
    });
  }

  /* ===================================================================
   * Protein heatmap
   * =================================================================== */
  function renderHeatmap() {
    const R = state.res;
    const src = $('#hm-source').value;
    const N = +$('#hm-n').value;
    const note = $('#hm-note');
    let idx;
    if (src === 'variable') {
      const v = R.M.map((row, i) => {
        const vals = Array.from(row).filter(isNum);
        return [i, vals.length >= 2 ? Stats.variance(vals, 1) : -1];
      }).filter(d => d[1] >= 0).sort((a, b) => b[1] - a[1]);
      idx = v.slice(0, N).map(d => d[0]);
      note.textContent = `The ${idx.length} proteins with the highest variance across samples, shown as row z-scores (missing values set to 0). Click a row to open the protein in the differential abundance view.`;
    } else {
      const { c, g1, g2, sign } = compView();
      const T = thresholds();
      const sig = c.results.filter(r => statusOf(r, T, sign) !== 0).sort((a, b) => (T.adj ? a.adj - b.adj : a.p - b.p));
      idx = sig.slice(0, N).map(r => r.i);
      note.textContent = idx.length
        ? `${idx.length} of ${sig.length} significant proteins for ${g1} vs ${g2}, ranked by p-value, shown as row z-scores. Click a row to inspect it.`
        : `No proteins pass the current thresholds for ${g1} vs ${g2}. Loosen them in the differential abundance tab.`;
    }
    if (idx.length < 2) { Plotly.purge('plot-heatmap'); return; }

    const Z = idx.map(i => {
      const row = Array.from(R.M[i]);
      const vals = row.filter(isNum);
      const m = Stats.mean(vals), sd = Math.sqrt(Stats.variance(vals, 1));
      return row.map(v => (isNum(v) && sd > 0 ? (v - m) / sd : 0));
    });
    const lim = Math.max(2, Math.min(4, Math.max(...Z.flat().map(Math.abs))));
    const labels = idx.map(i => R.proteins[i].gene || R.proteins[i].acc);
    const samples = R.samples;

    clusteredHeatmap('plot-heatmap', {
      z: Z, rowCluster: true, colCluster: $('#hm-colclust').checked,
      zmid: 0, zmin: -lim, zmax: lim, colorscale: divergingScale(), cbTitle: 'z',
      rowLabels: labels, colLabels: samples.map(s => s.sample),
      rowTickText: i => esc(labels[i]),
      colTickText: j => `${dot(groupColor(samples[j].group))} ${esc(samples[j].sample)}`,
      hover: (i, j) => {
        const p = R.proteins[idx[i]];
        const v = R.M[idx[i]][j];
        return `<b>${esc(p.gene || p.acc)}</b> ${esc(p.acc)}<br>${esc(samples[j].sample)} (${esc(samples[j].group)})<br>z = ${Z[i][j].toFixed(2)}<br>log2 = ${isNum(v) ? v.toFixed(2) : 'missing'}`;
      },
      height: Math.max(520, 180 + idx.length * (idx.length > 120 ? 9 : 16)),
      fileName: 'protein_heatmap',
      onRowClick: i => { selectProtein(idx[i]); showTab('de'); }
    });
  }

  /* ===================================================================
   * Differential abundance: volcano, inspector, table
   * =================================================================== */
  function compView() {
    const c = state.res.comps[state.de.comp];
    const sw = state.de.swap;
    return { c, g1: sw ? c.g2 : c.g1, g2: sw ? c.g1 : c.g2, sign: sw ? -1 : 1 };
  }
  function thresholds() {
    const fc = +$('#de-fc').value;
    const p = +$('#de-p').value;
    return { fc: isNum(fc) && fc >= 0 ? fc : 1, p: isNum(p) && p > 0 ? p : 0.05, adj: $('#de-yaxis').value === 'adj' };
  }
  function statusOf(r, T, sign) {
    const pv = T.adj ? r.adj : r.p;
    if (!(pv < T.p) || !(Math.abs(r.lfc) >= T.fc)) return 0;
    return sign * r.lfc > 0 ? 1 : -1;
  }
  const label = i => state.res.proteins[i].gene || state.res.proteins[i].acc;

  function renderDE() {
    renderVolcano();
    renderInspector();
    renderTable();
  }

  function renderVolcano() {
    const R = state.res;
    const { c, g1, g2, sign } = compView();
    const T = thresholds();
    const t = theme();
    const yv = r => -Math.log10(Math.max(T.adj ? r.adj : r.p, 1e-300));

    const buckets = { up: [], down: [], ns: [] };
    c.results.forEach((r, k) => {
      if (!isNum(T.adj ? r.adj : r.p)) return;
      const s = statusOf(r, T, sign);
      (s > 0 ? buckets.up : s < 0 ? buckets.down : buckets.ns).push(k);
    });

    const trace = (ks, name, color, size, opacity) => ({
      type: 'scattergl', mode: 'markers', name: `${name} (${ks.length})`,
      x: ks.map(k => sign * c.results[k].lfc), y: ks.map(k => yv(c.results[k])),
      customdata: ks, text: ks.map(k => label(c.results[k].i)),
      hovertext: ks.map(k => R.proteins[c.results[k].i].acc),
      hovertemplate: '<b>%{text}</b> %{hovertext}<br>log2 FC %{x:.2f}<br>−log10 p %{y:.2f}<extra></extra>',
      marker: { color, size, opacity }
    });
    const traces = [
      trace(buckets.ns, 'Not significant', t.ns, 6, 0.65),
      trace(buckets.down, `Higher in ${g2}`, t.down, 8, 0.9),
      trace(buckets.up, `Higher in ${g1}`, t.up, 8, 0.9)
    ];

    const annotations = [];
    const nLab = Math.max(0, Math.round(+$('#de-labels').value || 0));
    const sigSorted = [...buckets.up, ...buckets.down].sort((a, b) => yv(c.results[b]) - yv(c.results[a])).slice(0, nLab);
    for (const k of sigSorted) {
      const r = c.results[k];
      annotations.push({ x: sign * r.lfc, y: yv(r), text: esc(label(r.i)), showarrow: true, arrowhead: 0, arrowwidth: 1, arrowcolor: t.muted, ax: sign * r.lfc > 0 ? 18 : -18, ay: -16, font: { size: 10, color: t.ink } });
    }

    if (state.sel != null && c.byProtein.has(state.sel)) {
      const k = c.byProtein.get(state.sel);
      const r = c.results[k];
      if (isNum(T.adj ? r.adj : r.p)) {
        traces.push({
          type: 'scattergl', mode: 'markers', name: 'Selected', showlegend: false,
          x: [sign * r.lfc], y: [yv(r)], customdata: [k], text: [label(r.i)],
          hovertemplate: '<b>%{text}</b><extra></extra>',
          marker: { size: 18, color: 'rgba(0,0,0,0)', line: { color: t.ink, width: 2.5 }, symbol: 'circle' }
        });
        annotations.push({ x: sign * r.lfc, y: yv(r), text: `<b>${esc(label(r.i))}</b>`, showarrow: true, arrowhead: 0, ax: 0, ay: -30, font: { size: 12, color: t.ink }, bgcolor: t.surface, bordercolor: t.ink, borderwidth: 1, borderpad: 3 });
      }
    }

    const yLine = -Math.log10(T.p);
    const shapes = [{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: yLine, y1: yLine, line: { color: t.muted, width: 1, dash: 'dash' } }];
    if (T.fc > 0) {
      for (const x of [-T.fc, T.fc]) shapes.push({ type: 'line', yref: 'paper', y0: 0, y1: 1, x0: x, x1: x, line: { color: t.muted, width: 1, dash: 'dash' } });
    }

    const gd = document.getElementById('plot-volcano');
    Plotly.react(gd, traces, baseLayout({
      shapes, annotations,
      margin: { t: 56 },
      xaxis: { title: { text: `log2 fold change (${g1} / ${g2})` } },
      yaxis: { title: { text: T.adj ? '−log10 adjusted p' : '−log10 p' }, rangemode: 'tozero' },
      uirevision: `${state.de.comp}-${state.de.swap}-${T.adj}`
    }), cfg(`volcano_${g1}_vs_${g2}`)).then(() => {
      onClick(gd, ev => {
        const p = ev.points && ev.points[0];
        if (!p || p.customdata == null) return;
        selectProtein(c.results[p.customdata].i);
      });
    });

    const tested = c.results.length;
    $('#de-summary').innerHTML =
      `<b>${esc(g1)}</b> vs <b>${esc(g2)}</b>: ${tested.toLocaleString()} proteins tested, ` +
      `<span class="n-up">${buckets.up.length} higher in ${esc(g1)}</span> and ` +
      `<span class="n-down">${buckets.down.length} higher in ${esc(g2)}</span> ` +
      `at ${T.adj ? 'adjusted ' : ''}p &lt; ${T.p} and |log2 FC| ≥ ${T.fc}. Click any point to inspect it.`;
  }

  function selectProtein(i) {
    state.sel = i;
    if (state.activeTab === 'de') renderDE();
  }

  function renderInspector() {
    const el = $('#inspector');
    const R = state.res;
    if (state.sel == null) {
      el.innerHTML = '<p class="insp-empty">Click a point on the volcano plot, a row in the table below, or a row in the protein heatmap to see that protein here.</p>';
      return;
    }
    const p = R.proteins[state.sel];
    const { c, g1, g2, sign } = compView();
    const T = thresholds();
    const k = c.byProtein.get(state.sel);
    const r = k != null ? c.results[k] : null;
    const shortDesc = p.desc.split(/\s+OS=/)[0];

    let stats;
    if (r) {
      const lfc = sign * r.lfc;
      const st = statusOf(r, T, sign);
      stats = `
        <dl class="insp-stats">
          <div><dt>log2 FC</dt><dd class="${st > 0 ? 'up' : st < 0 ? 'down' : ''}">${fmt(lfc)}</dd></div>
          <div><dt>Fold change</dt><dd>${fmt(Math.pow(2, lfc))}</dd></div>
          <div><dt>p</dt><dd>${fmtP(r.p)}</dd></div>
          <div><dt>Adjusted p</dt><dd>${fmtP(r.adj)}</dd></div>
          <div><dt>Values in ${esc(g1)}</dt><dd>${sign > 0 ? r.n1 : r.n2}</dd></div>
          <div><dt>Values in ${esc(g2)}</dt><dd>${sign > 0 ? r.n2 : r.n1}</dd></div>
        </dl>`;
    } else {
      stats = `<p class="hint">Not tested in ${esc(g1)} vs ${esc(g2)}: fewer than ${R.S.minPer} values in one of the groups.</p>`;
    }

    const rows = R.comps.map((cc, ci) => {
      const kk = cc.byProtein.get(state.sel);
      const rr = kk != null ? cc.results[kk] : null;
      const cur = ci === state.de.comp;
      const s = cur ? sign : 1;
      const a = cur ? g1 : cc.g1, b = cur ? g2 : cc.g2;
      return `<tr class="${cur ? 'is-current' : ''}"><td>${esc(a)} vs ${esc(b)}</td><td class="num">${rr ? fmt(s * rr.lfc) : '–'}</td><td class="num">${rr ? fmtP(rr.adj) : '–'}</td></tr>`;
    }).join('');

    el.innerHTML = `
      <h3 class="insp-gene">${esc(p.gene || p.acc)}</h3>
      <p class="insp-acc">UniProt <a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(p.acc)}/entry" target="_blank" rel="noopener noreferrer">${esc(p.acc)}</a></p>
      ${shortDesc ? `<p class="insp-desc">${esc(shortDesc)}</p>` : ''}
      ${stats}
      <div id="insp-plot"></div>
      <table class="insp-comps">
        <thead><tr><th>Comparison</th><th class="num">log2 FC</th><th class="num">Adj. p</th></tr></thead>
        <tbody>${rows}</tbody>
      </table>`;

    const t = theme();
    const vals = R.M[state.sel];
    const traces = R.groups.map((g, gi) => {
      const ix = R.groupIdx[gi].filter(j => isNum(vals[j]));
      return {
        type: 'box', name: g, y: ix.map(j => vals[j]), text: ix.map(j => R.samples[j].sample),
        boxpoints: 'all', jitter: 0.35, pointpos: 0,
        marker: { color: groupColor(g), size: 8, line: { width: 1, color: t.surface } },
        line: { color: groupColor(g), width: 1.2 }, fillcolor: 'rgba(0,0,0,0)',
        hovertemplate: '%{text}: %{y:.2f}<extra></extra>'
      };
    });
    Plotly.react('insp-plot', traces, baseLayout({
      height: 260, showlegend: false,
      margin: { t: 8, b: 60, l: 48, r: 8 },
      xaxis: { tickangle: -30 },
      yaxis: { title: { text: R.corrected ? 'log2 (corrected)' : 'log2 abundance' } }
    }), cfg(`protein_${p.gene || p.acc}`));
  }

  const PAGE = 50;
  function tableRows() {
    const R = state.res;
    const { c, sign } = compView();
    const T = thresholds();
    const q = state.de.filter.trim().toLowerCase();
    let rows = c.results.map(r => {
      const p = R.proteins[r.i];
      return { r, p, lfc: sign * r.lfc, m1: sign > 0 ? r.m1 : r.m2, m2: sign > 0 ? r.m2 : r.m1, st: statusOf(r, T, sign) };
    });
    if (q) rows = rows.filter(o => o.p.gene.toLowerCase().includes(q) || o.p.acc.toLowerCase().includes(q) || o.p.desc.toLowerCase().includes(q));
    if (state.de.sigOnly) rows = rows.filter(o => o.st !== 0);
    const key = state.de.sortKey, dir = state.de.sortDir;
    const val = o => ({ gene: o.p.gene || o.p.acc, acc: o.p.acc, lfc: o.lfc, m1: o.m1, m2: o.m2, p: o.r.p, adj: o.r.adj }[key]);
    rows.sort((a, b) => {
      const va = val(a), vb = val(b);
      if (typeof va === 'string') return dir * va.localeCompare(vb);
      const na = !isNum(va), nb = !isNum(vb);
      if (na || nb) return na === nb ? 0 : na ? 1 : -1;
      return dir * (va - vb);
    });
    return rows;
  }

  function renderTable() {
    const { g1, g2 } = compView();
    $('#th-m1').textContent = `Mean ${g1}`;
    $('#th-m2').textContent = `Mean ${g2}`;
    $$('#de-table th').forEach(th => {
      const b = th.querySelector('button');
      if (!b) return;
      th.setAttribute('aria-sort', b.dataset.sort === state.de.sortKey ? (state.de.sortDir > 0 ? 'ascending' : 'descending') : 'none');
    });

    const rows = tableRows();
    const pages = Math.max(1, Math.ceil(rows.length / PAGE));
    state.de.page = Math.min(state.de.page, pages - 1);
    const slice = rows.slice(state.de.page * PAGE, (state.de.page + 1) * PAGE);
    $('#de-table tbody').innerHTML = slice.length ? slice.map(o => `
      <tr data-i="${o.r.i}" tabindex="0" class="${o.r.i === state.sel ? 'is-selected' : ''}">
        <td class="gene">${esc(o.p.gene || '–')}</td>
        <td>${esc(o.p.acc)}</td>
        <td class="num">${fmt(o.lfc)}</td>
        <td class="num">${fmt(o.m1)}</td>
        <td class="num">${fmt(o.m2)}</td>
        <td class="num">${fmtP(o.r.p)}</td>
        <td class="num">${fmtP(o.r.adj)}</td>
        <td><span class="call ${o.st > 0 ? 'call-up' : o.st < 0 ? 'call-down' : ''}">${o.st > 0 ? `Higher in ${esc(g1)}` : o.st < 0 ? `Higher in ${esc(g2)}` : 'Not significant'}</span></td>
      </tr>`).join('') : '<tr><td colspan="8" class="muted">No proteins match this filter.</td></tr>';

    $('#de-pager').innerHTML = `
      <button class="btn btn-quiet" type="button" data-page="-1" ${state.de.page === 0 ? 'disabled' : ''}>Previous</button>
      <span>Page ${state.de.page + 1} of ${pages}, ${rows.length.toLocaleString()} proteins</span>
      <button class="btn btn-quiet" type="button" data-page="1" ${state.de.page >= pages - 1 ? 'disabled' : ''}>Next</button>`;
  }

  function findProtein(q) {
    const R = state.res;
    q = q.trim().toLowerCase();
    if (!q || !R) return -1;
    let i = R.proteins.findIndex(p => p.gene.toLowerCase() === q || p.acc.toLowerCase() === q);
    if (i < 0) i = R.proteins.findIndex(p => p.gene.toLowerCase().startsWith(q) || p.acc.toLowerCase().startsWith(q));
    return i;
  }

  /* ===================================================================
   * Downloads (built in memory, handed to the browser, then released)
   * =================================================================== */
  function csvCell(v) {
    if (v == null) return '';
    if (typeof v === 'number') return Number.isFinite(v) ? String(v) : '';
    const s = String(v);
    return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  }
  const toCSV = rows => rows.map(r => r.map(csvCell).join(',')).join('\r\n') + '\r\n';

  function saveBlob(name, blob) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1500);
  }
  const saveText = (name, text, type = 'text/csv') => saveBlob(name, new Blob([text], { type: `${type};charset=utf-8` }));

  function matrixCSV(mat) {
    const R = state.res;
    const head = ['Accession', 'Gene_Symbol', 'Description', ...R.samples.map(s => s.sample)];
    return toCSV([head, ...mat.map((row, i) => [R.proteins[i].acc, R.proteins[i].gene, R.proteins[i].desc, ...Array.from(row)])]);
  }

  function deRows(ci, swap) {
    const R = state.res;
    const c = R.comps[ci];
    const sign = swap ? -1 : 1;
    const g1 = swap ? c.g2 : c.g1, g2 = swap ? c.g1 : c.g2;
    const T = thresholds();
    return c.results.map(r => {
      const p = R.proteins[r.i];
      const lfc = sign * r.lfc;
      const st = statusOf(r, T, sign);
      return [p.acc, p.gene, p.desc, g1, g2, swap ? r.n2 : r.n1, swap ? r.n1 : r.n2, swap ? r.m2 : r.m1, swap ? r.m1 : r.m2,
        lfc, Math.pow(2, lfc), sign * r.t, r.p, r.adj, st !== 0, st > 0 ? 'Up' : st < 0 ? 'Down' : ''];
    });
  }
  const DE_HEAD = ['Protein_ID', 'Gene_Symbol', 'Description', 'Group1', 'Group2', 'N_Group1', 'N_Group2', 'Mean_Group1', 'Mean_Group2', 'Log2FC', 'FoldChange', 'T_Statistic', 'P_Value', 'Adjusted_P_Value', 'Significant', 'Regulation'];
  const deCSV = (ci, swap) => toCSV([DE_HEAD, ...deRows(ci, swap)]);
  const deName = (ci, swap) => { const c = state.res.comps[ci]; return `DE_${safeName(swap ? c.g2 : c.g1)}_vs_${safeName(swap ? c.g1 : c.g2)}.csv`; };

  function pcaCSV(P) {
    const R = state.res;
    const k = P.explained.length;
    return toCSV([
      ['Sample', 'Group', 'Batch', ...P.explained.map((_, i) => `PC${i + 1}`)],
      ...R.samples.map((s, j) => [s.sample, s.group, s.batch, ...P.scores[j].slice(0, k)]),
      [],
      ['Explained_variance_ratio', '', '', ...P.explained]
    ]);
  }

  function corrCSV() {
    const R = state.res;
    const names = R.samples.map(s => s.sample);
    return toCSV([['Sample', ...names], ...R.pca.corr.map((row, i) => [names[i], ...row])]);
  }

  function metadataCSV() {
    const R = state.res;
    return toCSV([['Sample', 'Group', 'Batch', 'Source_Column'], ...R.samples.map(s => [s.sample, s.group, s.batch, s.header])]);
  }

  function parametersText() {
    const R = state.res, S = R.S, T = thresholds();
    return [
      'LFQ Explorer analysis parameters',
      `Run: ${R.ranAt.toISOString()}`,
      `Input file: ${state.fileName} (sheet: ${state.parsed.sheet})`,
      `Quantity: ${S.quant}`,
      '',
      'Filters',
      `  Master proteins only: ${S.masterOnly}`,
      `  Remove contaminants: ${S.removeContam}`,
      `  High FDR confidence only: ${S.highConf}`,
      `  Minimum unique peptides: ${S.minPep}`,
      `  Detected in at least ${S.minRep} replicates of one group`,
      ...R.counts.map(([k, v]) => `  ${k}: ${v}`),
      '',
      'Processing',
      '  Duplicated accessions collapsed by median; values <= 0 treated as missing; log2 transform.',
      `  Batch correction: ${S.batch ? `parametric ComBat, group protected: ${S.protect}, original missing values kept missing: ${S.remask}` : 'none'}`,
      '  PCA and clustering: per-protein median imputation, standardization (z-score), average linkage on Euclidean distance.',
      '',
      'Statistics',
      `  Welch's t-test on log2 values, minimum ${S.minPer} values per group; Benjamini-Hochberg FDR per comparison.`,
      `  Significance thresholds used in exported files: ${T.adj ? 'adjusted' : 'raw'} p < ${T.p}, |log2FC| >= ${T.fc}`,
      '',
      'Groups',
      ...R.groups.map((g, i) => `  ${g}: ${R.groupIdx[i].map(j => R.samples[j].sample).join(', ')}`)
    ].join('\n') + '\n';
  }

  function downloadItems() {
    const R = state.res;
    const items = [
      { name: 'Log2_Abundance_Filtered.csv', title: 'Filtered log2 abundance matrix', desc: `${R.proteins.length.toLocaleString()} proteins × ${R.samples.length} samples, before any batch correction`, make: () => matrixCSV(R.raw) }
    ];
    if (R.corrected) items.push({ name: 'ComBat_Corrected_Abundance.csv', title: 'Batch-corrected log2 abundance matrix', desc: 'ComBat output used for PCA, heatmaps and statistics', make: () => matrixCSV(R.corrected) });
    items.push({ name: 'All_Differential_Expression.csv', title: 'All comparisons combined', desc: 'Every pairwise test in one table; Significant uses the current thresholds', make: () => toCSV([DE_HEAD, ...R.comps.flatMap((_, ci) => deRows(ci, false))]) });
    R.comps.forEach((c, ci) => items.push({ name: deName(ci, false), title: `${c.g1} vs ${c.g2}`, desc: 'Differential abundance results', make: () => deCSV(ci, false) }));
    if (R.pca) {
      items.push({ name: 'PCA_Scores.csv', title: 'PCA scores', desc: 'Sample coordinates and explained variance', make: () => pcaCSV(R.pca) });
      items.push({ name: 'Sample_Correlation.csv', title: 'Sample correlation matrix', desc: 'Pearson r between samples', make: corrCSV });
    }
    items.push({ name: 'Sample_Metadata.csv', title: 'Sample design', desc: 'Sample, group and batch as used in this run', make: metadataCSV });
    items.push({ name: 'Analysis_Parameters.txt', title: 'Analysis parameters', desc: 'Filters, methods and thresholds for your records', make: parametersText, type: 'text/plain' });
    return items;
  }

  function renderDownloads() {
    const items = downloadItems();
    $('#dl-list').innerHTML = items.map((it, i) => `
      <li>
        <div class="dl-text"><b>${esc(it.title)}</b><span>${esc(it.desc)}. ${esc(it.name)}</span></div>
        <button class="btn btn-quiet" type="button" data-dl="${i}">Download</button>
      </li>`).join('');
  }

  async function downloadZip() {
    const btn = $('#dl-zip');
    btn.disabled = true;
    try {
      await needLibs('JSZip');
      setStatus('');
      const zip = new JSZip();
      for (const it of downloadItems()) zip.file(it.name, it.make());
      const blob = await zip.generateAsync({ type: 'blob', compression: 'DEFLATE' });
      saveBlob(`LFQ_Explorer_results_${safeName(state.fileName.replace(/\.[^.]+$/, ''))}.zip`, blob);
    } catch (e) {
      setStatus('');
      notice(e.message);
    } finally {
      btn.disabled = false;
    }
  }

  /* ===================================================================
   * Event wiring
   * =================================================================== */
  function wire() {
    const fi = $('#file-input');
    fi.addEventListener('change', () => {
      const file = fi.files && fi.files[0];
      fi.value = '';           // allows picking the same file again
      if (file) loadFile(file);
    });

    $('#sheet-select').addEventListener('change', e => {
      try { resetResults(); parseSheet(e.target.value); } catch (err) { notice(err.message); }
    });
    $('#quant-select').addEventListener('change', e => {
      state.parsed.quantIdx = +e.target.value;
      buildSamples(); renderSampleTable(); markStale();
    });

    $('#open-samples').addEventListener('click', () => showTab('samples'));
    $('#sample-tbody').addEventListener('input', onSampleEdit);
    $('#sample-tbody').addEventListener('change', onSampleEdit);

    ['#f-master', '#f-contam', '#f-highconf', '#f-minpep', '#f-minrep', '#b-protect', '#b-remask'].forEach(id => $(id).addEventListener('change', markStale));
    $('#s-minper').addEventListener('change', () => { renderGroupSummary(); markStale(); });
    $('#b-enable').addEventListener('change', e => {
      $('#batch-fields').hidden = !e.target.checked;
      if (state.parsed) { renderSampleTable(); if (e.target.checked) showTab('samples'); }
      markStale();
    });
    $('#b-file').addEventListener('change', e => { if (e.target.files[0]) loadBatchFile(e.target.files[0]); e.target.value = ''; });

    $('#run-btn').addEventListener('click', runAnalysis);

    $('#tabs').addEventListener('click', e => {
      const b = e.target.closest('button[data-tab]');
      if (b && !b.disabled) showTab(b.dataset.tab);
    });
    $('#tabs').addEventListener('keydown', e => {
      if (e.key !== 'ArrowRight' && e.key !== 'ArrowLeft') return;
      const tabs = $$('#tabs button:not(:disabled)');
      const i = tabs.indexOf(document.activeElement);
      if (i < 0) return;
      const next = tabs[(i + (e.key === 'ArrowRight' ? 1 : -1) + tabs.length) % tabs.length];
      next.focus(); showTab(next.dataset.tab);
    });

    // Overview / PCA / correlation / heatmap controls
    $('#ov-corr').addEventListener('change', renderOverview);
    ['#pca-x', '#pca-y', '#pca-before'].forEach(id => $(id).addEventListener('change', renderPCA));
    $('#corr-annot').addEventListener('change', renderCorrelation);
    ['#hm-source', '#hm-n', '#hm-colclust'].forEach(id => $(id).addEventListener('change', renderHeatmap));

    // Differential abundance controls
    $('#de-comp').addEventListener('change', e => { state.de.comp = +e.target.value; state.de.swap = false; state.de.page = 0; renderDE(); });
    $('#de-swap').addEventListener('click', () => { state.de.swap = !state.de.swap; renderDE(); });
    ['#de-fc', '#de-p', '#de-yaxis', '#de-labels'].forEach(id => $(id).addEventListener('change', () => { state.de.page = 0; renderDE(); }));
    $('#de-search').addEventListener('change', e => {
      const i = findProtein(e.target.value);
      if (i >= 0) selectProtein(i);
      else if (e.target.value.trim()) notice(`No protein matching “${e.target.value.trim()}” passed the filters.`, 'warn');
    });
    let filterTimer;
    $('#de-filter').addEventListener('input', e => {
      clearTimeout(filterTimer);
      filterTimer = setTimeout(() => { state.de.filter = e.target.value; state.de.page = 0; renderTable(); }, 150);
    });
    $('#de-sigonly').addEventListener('change', e => { state.de.sigOnly = e.target.checked; state.de.page = 0; renderTable(); });
    $('#de-table thead').addEventListener('click', e => {
      const b = e.target.closest('button[data-sort]');
      if (!b) return;
      const key = b.dataset.sort;
      if (state.de.sortKey === key) state.de.sortDir *= -1;
      else { state.de.sortKey = key; state.de.sortDir = key === 'gene' || key === 'acc' || key === 'p' || key === 'adj' ? 1 : -1; }
      renderTable();
    });
    const pickRow = e => {
      const tr = e.target.closest('tr[data-i]');
      if (tr) selectProtein(+tr.dataset.i);
    };
    $('#de-table tbody').addEventListener('click', pickRow);
    $('#de-table tbody').addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); pickRow(e); } });
    $('#de-pager').addEventListener('click', e => {
      const b = e.target.closest('button[data-page]');
      if (!b) return;
      state.de.page += +b.dataset.page;
      renderTable();
    });
    $('#de-dl').addEventListener('click', () => saveText(deName(state.de.comp, state.de.swap), deCSV(state.de.comp, state.de.swap)));

    // Downloads
    $('#dl-list').addEventListener('click', e => {
      const b = e.target.closest('button[data-dl]');
      if (!b) return;
      const it = downloadItems()[+b.dataset.dl];
      saveText(it.name, it.make(), it.type || 'text/csv');
    });
    $('#dl-zip').addEventListener('click', downloadZip);

    // Re-theme plots when the system colour scheme changes
    window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', renderActive);
  }

  /* ===================================================================
   * Library loading with fallbacks: local lib/ folder first, then CDNs.
   * =================================================================== */
  const LIBS = {
    XLSX: { label: 'SheetJS (Excel reader)', urls: [
      'lib/xlsx.full.min.js',
      'https://cdn.jsdelivr.net/npm/xlsx@0.18.5/dist/xlsx.full.min.js',
      'https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js',
      'https://unpkg.com/xlsx@0.18.5/dist/xlsx.full.min.js'
    ] },
    Plotly: { label: 'Plotly (plots)', urls: [
      'lib/plotly.min.js',
      'https://cdn.plot.ly/plotly-2.35.2.min.js',
      'https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js',
      'https://unpkg.com/plotly.js-dist-min@2.35.2/plotly.min.js'
    ] },
    JSZip: { label: 'JSZip (ZIP downloads)', urls: [
      'lib/jszip.min.js',
      'https://cdn.jsdelivr.net/npm/jszip@3.10.1/dist/jszip.min.js',
      'https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js',
      'https://unpkg.com/jszip@3.10.1/dist/jszip.min.js'
    ] }
  };

  function loadScript(url) {
    return new Promise((resolve, reject) => {
      const el = document.createElement('script');
      el.src = url;
      el.async = true;
      el.onload = () => resolve();
      el.onerror = () => { el.remove(); reject(new Error(url)); };
      document.head.appendChild(el);
    });
  }

  const libPromises = {};
  function ensureLib(name) {
    if (window[name]) return Promise.resolve();
    if (!libPromises[name]) {
      libPromises[name] = (async () => {
        for (const url of LIBS[name].urls) {
          try { await loadScript(url); if (window[name]) return; } catch (_) { /* try the next source */ }
        }
        throw new Error(`${LIBS[name].label} could not be loaded. Check the internet connection or firewall, or place the library file in a lib/ folder next to index.html (see README).`);
      })();
      libPromises[name].catch(() => { delete libPromises[name]; });   // allow a retry later
    }
    return libPromises[name];
  }

  async function needLibs(...names) {
    if (names.every(n => window[n])) return;
    setStatus('Loading components…', true);
    await Promise.all(names.map(ensureLib));
  }

  document.addEventListener('DOMContentLoaded', () => {
    drawHero();
    wire();
    // Start fetching libraries straight away so they are ready when a file is chosen.
    needLibs('XLSX', 'Plotly')
      .then(() => setStatus(''))
      .catch(e => { setStatus(''); notice(e.message); });
  });
})();
