/* stats.js
 * Numerical routines for LFQ Explorer. No dependencies, runs entirely in the browser.
 * Mirrors the Python pipeline: Welch's t-test, Benjamini–Hochberg FDR, median imputation,
 * standardisation, PCA, average-linkage hierarchical clustering and parametric ComBat.
 */
'use strict';

const Stats = (() => {
  const isNum = v => typeof v === 'number' && Number.isFinite(v);

  function mean(a) {
    let s = 0, n = 0;
    for (const v of a) { s += v; n++; }
    return n ? s / n : NaN;
  }

  function variance(a, ddof = 1) {
    const n = a.length;
    if (n - ddof <= 0) return NaN;
    const m = mean(a);
    let s = 0;
    for (const v of a) { const d = v - m; s += d * d; }
    return s / (n - ddof);
  }

  function median(a) {
    const b = Array.from(a).filter(isNum).sort((x, y) => x - y);
    const n = b.length;
    if (!n) return NaN;
    const h = n >> 1;
    return n % 2 ? b[h] : (b[h - 1] + b[h]) / 2;
  }

  /* ---------- Student t distribution via the regularised incomplete beta ---------- */
  const LANCZOS = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];

  function lgamma(x) {
    if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
    x -= 1;
    let a = LANCZOS[0];
    const t = x + 7.5;
    for (let i = 1; i < 9; i++) a += LANCZOS[i] / (x + i);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }

  function betacf(a, b, x) {
    const MAXIT = 300, EPS = 3e-14, FPMIN = 1e-300;
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }

  function ibeta(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
    if (x < (a + 1) / (a + b + 2)) return bt * betacf(a, b, x) / a;
    return 1 - bt * betacf(b, a, 1 - x) / b;
  }

  function tTwoSided(t, df) {
    if (!isNum(t) || !isNum(df) || df <= 0) return NaN;
    const p = ibeta(df / (df + t * t), df / 2, 0.5);
    return Math.min(1, Math.max(0, p));
  }

  /* Welch's unequal-variance t-test (scipy.stats.ttest_ind, equal_var=False) */
  function welch(a, b) {
    const n1 = a.length, n2 = b.length;
    const m1 = mean(a), m2 = mean(b);
    const s1 = variance(a, 1) / n1, s2 = variance(b, 1) / n2;
    const se2 = s1 + s2;
    if (!(se2 > 0)) return { t: NaN, df: NaN, p: NaN };
    const t = (m1 - m2) / Math.sqrt(se2);
    const df = se2 * se2 / (s1 * s1 / (n1 - 1) + s2 * s2 / (n2 - 1));
    return { t, df, p: tTwoSided(t, df) };
  }

  /* Benjamini–Hochberg; NaN p-values are ignored and stay NaN */
  function bh(p) {
    const out = new Array(p.length).fill(NaN);
    const idx = [];
    for (let i = 0; i < p.length; i++) if (isNum(p[i])) idx.push(i);
    const m = idx.length;
    idx.sort((a, b) => p[a] - p[b]);
    let prev = 1;
    for (let r = m - 1; r >= 0; r--) {
      const i = idx[r];
      const v = Math.min(prev, p[i] * m / (r + 1));
      out[i] = v;
      prev = v;
    }
    return out;
  }

  /* ---------- linear algebra ---------- */
  function invert(A) {
    const n = A.length;
    const M = A.map((r, i) => [...r, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
    for (let c = 0; c < n; c++) {
      let piv = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
      if (Math.abs(M[piv][c]) < 1e-10) return null;
      [M[c], M[piv]] = [M[piv], M[c]];
      const d = M[c][c];
      for (let j = 0; j < 2 * n; j++) M[c][j] /= d;
      for (let r = 0; r < n; r++) {
        if (r === c) continue;
        const f = M[r][c];
        if (f === 0) continue;
        for (let j = 0; j < 2 * n; j++) M[r][j] -= f * M[c][j];
      }
    }
    return M.map(r => r.slice(n));
  }

  /* Jacobi eigen-decomposition of a small symmetric matrix. Values sorted descending. */
  function eigSym(S) {
    const n = S.length;
    const a = S.map(r => Array.from(r));
    const V = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] * a[p][q];
      if (off < 1e-22) break;
      for (let p = 0; p < n; p++) {
        for (let q = p + 1; q < n; q++) {
          if (Math.abs(a[p][q]) < 1e-300) continue;
          const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
          const t = theta === 0 ? 1 : Math.sign(theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          const c = 1 / Math.sqrt(t * t + 1), s = t * c;
          for (let k = 0; k < n; k++) {
            const akp = a[k][p], akq = a[k][q];
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
          }
          for (let k = 0; k < n; k++) {
            const apk = a[p][k], aqk = a[q][k];
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
          }
          for (let k = 0; k < n; k++) {
            const vkp = V[k][p], vkq = V[k][q];
            V[k][p] = c * vkp - s * vkq;
            V[k][q] = s * vkp + c * vkq;
          }
        }
      }
    }
    const order = [...Array(n).keys()].sort((i, j) => a[j][j] - a[i][i]);
    return {
      values: order.map(i => a[i][i]),
      vectors: V.map(row => order.map(i => row[i]))
    };
  }

  /* Column-wise z-score (sklearn StandardScaler, ddof = 0). X: rows = samples. Mutates X. */
  function standardizeColumns(X) {
    const n = X.length, p = X[0].length;
    for (let j = 0; j < p; j++) {
      let m = 0;
      for (let i = 0; i < n; i++) m += X[i][j];
      m /= n;
      let ss = 0;
      for (let i = 0; i < n; i++) { const d = X[i][j] - m; ss += d * d; }
      const sd = Math.sqrt(ss / n);
      for (let i = 0; i < n; i++) X[i][j] = sd > 0 ? (X[i][j] - m) / sd : 0;
    }
    return X;
  }

  /* PCA on a centred matrix (rows = samples) via the n×n Gram matrix. */
  function pca(X, k) {
    const n = X.length;
    const G = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        let s = 0;
        const xi = X[i], xj = X[j];
        for (let c = 0; c < xi.length; c++) s += xi[c] * xj[c];
        G[i][j] = G[j][i] = s;
      }
    }
    const { values, vectors } = eigSym(G);
    const total = values.reduce((s, v) => s + Math.max(v, 0), 0) || 1;
    k = Math.min(k, n);
    const scores = Array.from({ length: n }, () => new Array(k).fill(0));
    const explained = [];
    for (let c = 0; c < k; c++) {
      const lam = Math.max(values[c], 0);
      // deterministic sign: largest-magnitude score is positive
      let maxAbs = 0, sign = 1;
      for (let i = 0; i < n; i++) if (Math.abs(vectors[i][c]) > maxAbs) { maxAbs = Math.abs(vectors[i][c]); sign = Math.sign(vectors[i][c]) || 1; }
      for (let i = 0; i < n; i++) scores[i][c] = sign * vectors[i][c] * Math.sqrt(lam);
      explained.push(lam / total);
    }
    return { scores, explained };
  }

  /* Pearson correlation between rows (np.corrcoef). */
  function corrRows(X) {
    const n = X.length;
    const Z = X.map(r => {
      const m = mean(r);
      const z = new Float64Array(r.length);
      let ss = 0;
      for (let j = 0; j < r.length; j++) { z[j] = r[j] - m; ss += z[j] * z[j]; }
      const s = Math.sqrt(ss) || 1;
      for (let j = 0; j < r.length; j++) z[j] /= s;
      return z;
    });
    const C = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        let s = 0;
        for (let c = 0; c < Z[i].length; c++) s += Z[i][c] * Z[j][c];
        C[i][j] = C[j][i] = s;
      }
    }
    return C;
  }

  /* Average-linkage agglomerative clustering on Euclidean distance (seaborn clustermap default).
   * Returns the leaf order and dendrogram line segments aligned to positions 0..n-1. */
  function hclust(rows) {
    const n = rows.length;
    if (n < 2) return { order: n ? [0] : [], seg: { pos: [], h: [] }, maxH: 0 };
    const N = 2 * n - 1;
    const D = Array.from({ length: N }, () => new Float64Array(N));
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        let s = 0;
        const a = rows[i], b = rows[j];
        for (let c = 0; c < a.length; c++) { const d = a[c] - b[c]; s += d * d; }
        D[i][j] = D[j][i] = Math.sqrt(s);
      }
    }
    const size = new Array(N).fill(1);
    const height = new Float64Array(N);
    const kids = new Array(N);
    let active = [...Array(n).keys()];
    for (let k = 0; k < n - 1; k++) {
      let best = Infinity, bi = -1, bj = -1;
      for (let x = 0; x < active.length; x++) {
        const ax = active[x];
        for (let y = x + 1; y < active.length; y++) {
          const d = D[ax][active[y]];
          if (d < best) { best = d; bi = ax; bj = active[y]; }
        }
      }
      const id = n + k;
      size[id] = size[bi] + size[bj];
      height[id] = best;
      kids[id] = [bi, bj];
      active = active.filter(a => a !== bi && a !== bj);
      for (const c of active) {
        const d = (D[bi][c] * size[bi] + D[bj][c] * size[bj]) / size[id];
        D[id][c] = D[c][id] = d;
      }
      active.push(id);
    }
    const root = N - 1;
    const order = [];
    const stack = [root];
    while (stack.length) {
      const id = stack.pop();
      if (id < n) order.push(id);
      else { stack.push(kids[id][1]); stack.push(kids[id][0]); }
    }
    const pos = new Float64Array(N);
    order.forEach((leaf, i) => { pos[leaf] = i; });
    const segPos = [], segH = [];
    for (let id = n; id < N; id++) {
      const [a, b] = kids[id];
      pos[id] = (pos[a] + pos[b]) / 2;
      segPos.push(pos[a], pos[a], pos[b], pos[b], null);
      segH.push(height[a], height[id], height[id], height[b], null);
    }
    return { order, seg: { pos: segPos, h: segH }, maxH: height[root] };
  }

  /* Parametric empirical-Bayes ComBat (Johnson, Li & Rabinovic 2007), as in sva / neuroCombat.
   * Y: array of Float64Array (features × samples), no missing values.
   * batch: array of batch labels. covariate: array of group labels to protect, or null. */
  function combat(Y, batch, covariate) {
    const G = Y.length, n = batch.length;
    const bLevels = [...new Set(batch)];
    const B = bLevels.length;
    if (B < 2) throw new Error('Batch correction needs at least two different batch labels.');
    const bIdx = batch.map(b => bLevels.indexOf(b));
    const members = bLevels.map((_, k) => bIdx.map((v, j) => (v === k ? j : -1)).filter(j => j >= 0));
    const tooSmall = bLevels.filter((_, k) => members[k].length < 2);
    if (tooSmall.length) throw new Error(`Each batch needs at least two samples. Batch “${tooSmall.join('”, “')}” has only one.`);

    const cLevels = covariate ? [...new Set(covariate)] : [];
    const P = B + Math.max(cLevels.length - 1, 0);
    if (n <= P) throw new Error('Too few samples for this batch and group design. Add samples or turn off group protection.');
    const X = batch.map((_, j) => {
      const row = new Array(P).fill(0);
      row[bIdx[j]] = 1;
      if (covariate) {
        const ci = cLevels.indexOf(covariate[j]);
        if (ci > 0) row[B + ci - 1] = 1;
      }
      return row;
    });
    const XtX = Array.from({ length: P }, (_, a) => Array.from({ length: P }, (_, b) => {
      let s = 0;
      for (let j = 0; j < n; j++) s += X[j][a] * X[j][b];
      return s;
    }));
    const inv = invert(XtX);
    if (!inv) throw new Error('Batches and groups are completely confounded, so batch effects cannot be separated from biology. Check the batch labels or turn off group protection.');
    const H = Array.from({ length: P }, (_, a) => {
      const r = new Float64Array(n);
      for (let j = 0; j < n; j++) {
        let s = 0;
        for (let b = 0; b < P; b++) s += inv[a][b] * X[j][b];
        r[j] = s;
      }
      return r;
    });
    const nb = members.map(m => m.length);

    const sdat = new Array(G), stand = new Array(G), sdv = new Float64Array(G);
    for (let g = 0; g < G; g++) {
      const y = Y[g];
      const beta = new Float64Array(P);
      for (let a = 0; a < P; a++) {
        let s = 0;
        const h = H[a];
        for (let j = 0; j < n; j++) s += h[j] * y[j];
        beta[a] = s;
      }
      let grand = 0;
      for (let k = 0; k < B; k++) grand += (nb[k] / n) * beta[k];
      let ss = 0;
      const st = new Float64Array(n);
      for (let j = 0; j < n; j++) {
        let fit = 0, cov = 0;
        for (let a = 0; a < P; a++) {
          fit += X[j][a] * beta[a];
          if (a >= B) cov += X[j][a] * beta[a];
        }
        const r = y[j] - fit;
        ss += r * r;
        st[j] = grand + cov;
      }
      const sd = Math.sqrt(ss / n);
      sdv[g] = sd;
      stand[g] = st;
      const s = new Float64Array(n);
      for (let j = 0; j < n; j++) s[j] = sd > 1e-12 ? (y[j] - st[j]) / sd : 0;
      sdat[g] = s;
    }

    const gHat = members.map(() => new Float64Array(G));
    const dHat = members.map(() => new Float64Array(G));
    for (let g = 0; g < G; g++) {
      for (let k = 0; k < B; k++) {
        const vals = members[k].map(j => sdat[g][j]);
        gHat[k][g] = mean(vals);
        dHat[k][g] = variance(vals, 1);
      }
    }

    const gStar = [], dStar = [];
    for (let k = 0; k < B; k++) {
      const gbar = mean(gHat[k]);
      const t2 = variance(gHat[k], 1);
      const m = mean(dHat[k]);
      const s2 = Math.max(variance(dHat[k], 1), 1e-12);
      const aPrior = (2 * s2 + m * m) / s2;
      const bPrior = (m * s2 + m * m * m) / s2;
      const nk = nb[k];
      const gs = new Float64Array(G), ds = new Float64Array(G);
      for (let g = 0; g < G; g++) {
        let gOld = gHat[k][g], dOld = dHat[k][g];
        let gNew = gOld, dNew = dOld;
        for (let it = 0; it < 1000; it++) {
          const denom = nk * t2 + dOld;
          gNew = denom > 0 ? (nk * t2 * gHat[k][g] + dOld * gbar) / denom : gHat[k][g];
          let sum2 = 0;
          for (const j of members[k]) { const d = sdat[g][j] - gNew; sum2 += d * d; }
          dNew = (0.5 * sum2 + bPrior) / (nk / 2 + aPrior - 1);
          const change = Math.max(
            Math.abs(gNew - gOld) / (Math.abs(gOld) + 1e-12),
            Math.abs(dNew - dOld) / (Math.abs(dOld) + 1e-12)
          );
          gOld = gNew; dOld = dNew;
          if (change < 1e-4) break;
        }
        gs[g] = gNew;
        ds[g] = dNew > 0 ? dNew : 1;
      }
      gStar.push(gs);
      dStar.push(ds);
    }

    return Y.map((y, g) => {
      const out = new Float64Array(n);
      for (let j = 0; j < n; j++) {
        const k = bIdx[j];
        out[j] = sdv[g] > 1e-12
          ? ((sdat[g][j] - gStar[k][g]) / Math.sqrt(dStar[k][g])) * sdv[g] + stand[g][j]
          : y[j];
      }
      return out;
    });
  }

  return { isNum, mean, variance, median, welch, bh, tTwoSided, eigSym, standardizeColumns, pca, corrRows, hclust, combat, invert };
})();
