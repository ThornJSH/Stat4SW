/**
 * 관계 분석: 상관분석(피어슨, 스피어만, 통제변수가 있으면 편상관), 교차분석(카이제곱 등). 순수 함수: (dataset, request) → 결과 객체.
 * 결과표 구성: 상관 행렬(변수마다 상관계수 / 유의확률(양측) / N).
 *
 * 켄달 타우-b: 동순위를 보정한 순위상관. 타우-b = (일치쌍 − 불일치쌍)/√((n₀ − n₁)(n₀ − n₂)), n₀ = n(n−1)/2, n₁·n₂ = 각 변수의 동순위 쌍 수.
 * 불일치쌍은 x, y 순으로 정렬한 뒤 y 의 역전 수를 병합 정렬로 센다 (O(n log n), Knight 1966). 유의확률은 scipy 와 같은 규칙: 동순위가 없고
 * n ≤ 33 (또는 불일치쌍이 극단적으로 적음)이면 정확 분포, 아니면 동순위 보정 분산의 정규 근사. 통제변수(편상관)와는 함께 쓰지 않는다.
 *
 * 편상관(통제변수 칸을 채운 경우): 통제변수 Z 의 선형 영향을 뺀 상관. 상관행렬에서 Σ = R_vv − R_vz R_zz⁻¹ R_zv 를 구해
 * r_ij·Z = Σ_ij/√(Σ_ii Σ_jj). 유의확률은 t = r√(df/(1 − r²)), df = N − 2 − (통제변수 수). Fisher z 신뢰구간은 SE = 1/√(N − 통제변수 수 − 3).
 * 스피어만은 모든 변수를 순위(동순위 평균)로 바꾼 뒤 같은 계산. 결측은 항상 목록별 제거.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function relDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var desc = require('./Stats_Descriptive.js');
    var ch = require('./Charts.js');
    return {
      getVars: data.getVars, selectComplete: data.selectComplete, compareValues_: data.compareValues_,
      pt: dist.pt, pValueT: dist.pValueT, qnorm: dist.qnorm, pnorm: dist.pnorm, pchisq: dist.pchisq, lgamma: dist.lgamma,
      centralSums: desc.centralSums, scatterBlock: ch.scatterBlock, CHART_MAX_PAIRS_: ch.CHART_MAX_PAIRS_,
      olsFit: require('./Stats_Regression.js').olsFit,
      advCorrelation_: require('./Stats_Assumptions.js').advCorrelation_, advAppend_: require('./Stats_Assumptions.js').advAppend_
    };
  }
  return {
    getVars: getVars, selectComplete: selectComplete, compareValues_: compareValues_,
    pt: pt, pValueT: pValueT, qnorm: qnorm, pnorm: pnorm, pchisq: pchisq, lgamma: lgamma, centralSums: centralSums,
    scatterBlock: scatterBlock, CHART_MAX_PAIRS_: CHART_MAX_PAIRS_, olsFit: olsFit,
    advCorrelation_: advCorrelation_, advAppend_: advAppend_
  };
}

function relName_(v) { return v.label ? v.label : v.name; }

/* ---------- 상관계수 ---------- */

/**
 * 피어슨 상관계수와 양측 p값 (n ≥ 3). 표준편차가 0인 변수가 있으면 r 은 null.
 * 편차의 곱합은 평균 반올림 오차를 보정한다.
 */
function pearsonR(x, y) {
  var d = relDeps_();
  var n = x.length;
  var res = { n: n, r: null, df: n - 2, p: null };
  if (n < 3) return res;
  var cx = d.centralSums(x), cy = d.centralSums(y);
  var sxy = 0, ex = 0, ey = 0, i;
  for (i = 0; i < n; i++) {
    var dx = x[i] - cx.mean, dy = y[i] - cy.mean;
    sxy += dx * dy; ex += dx; ey += dy;
  }
  sxy -= ex * ey / n;
  if (!(cx.c2 > 0) || !(cy.c2 > 0)) return res;
  res.r = Math.max(-1, Math.min(1, sxy / Math.sqrt(cx.c2 * cy.c2)));
  res.p = Math.abs(res.r) === 1 ? 0 : d.pValueT(res.r * Math.sqrt((n - 2) / (1 - res.r * res.r)), n - 2);
  return res;
}

/** 순위 (동순위는 평균순위) */
function rankAverage_(values) {
  var idx = values.map(function (v, i) { return i; }).sort(function (a, b) { return values[a] - values[b]; });
  var ranks = new Array(values.length);
  var i = 0;
  while (i < idx.length) {
    var j = i;
    while (j + 1 < idx.length && values[idx[j + 1]] === values[idx[i]]) j++;
    var avg = (i + j) / 2 + 1;
    for (var k = i; k <= j; k++) ranks[idx[k]] = avg;
    i = j + 1;
  }
  return ranks;
}

/** 스피어만 순위상관: 순위의 피어슨 상관, p값은 t 근사 (자유도 n−2) */
function spearmanR(x, y) {
  return pearsonR(rankAverage_(x), rankAverage_(y));
}

/* ---------- 켄달 타우-b ---------- */

/** 정렬된 값에서 동순위 묶음 통계: t = Σ t(t−1)/2, t0 = Σ t(t−1)(t−2), t1 = Σ t(t−1)(2t+5) */
function relTieStats_(sorted) {
  var out = { t: 0, t0: 0, t1: 0 }, i = 0, n = sorted.length;
  while (i < n) {
    var j = i;
    while (j + 1 < n && sorted[j + 1] === sorted[i]) j++;
    var g = j - i + 1;
    if (g > 1) { out.t += g * (g - 1) / 2; out.t0 += g * (g - 1) * (g - 2); out.t1 += g * (g - 1) * (2 * g + 5); }
    i = j + 1;
  }
  return out;
}

/** 배열의 역전 수(i < j 이고 a[i] > a[j] 인 쌍의 수): 아래에서 위로 병합 정렬 */
function relInversions_(arr) {
  var n = arr.length, a = new Float64Array(arr), b = new Float64Array(n), inv = 0, w, i, lo, mid, hi, p, q, k;
  for (w = 1; w < n; w *= 2) {
    for (lo = 0; lo < n; lo += 2 * w) {
      mid = Math.min(lo + w, n); hi = Math.min(lo + 2 * w, n);
      p = lo; q = mid; k = lo;
      while (p < mid && q < hi) {
        if (a[q] < a[p]) { b[k++] = a[q++]; inv += mid - p; } else b[k++] = a[p++];
      }
      while (p < mid) b[k++] = a[p++];
      while (q < hi) b[k++] = a[q++];
    }
    var tmp = a; a = b; b = tmp;
  }
  return inv;
}

/** 동순위가 없을 때 불일치쌍 수 D (순열의 역전 수)의 분포에서 양측 정확 유의확률: 2 P(D ≤ min(dis, tot − dis)) (최대 1) */
function relKendallExactP_(n, dis, tot) {
  var c = Math.min(dis, tot - dis);
  if (2 * c === tot) return 1;
  if (c <= 1) {   // 닫힌 식: P(D=0) = 1/n!, P(D=1) = (n−1)/n!
    return Math.min(1, 2 * (c === 0 ? 1 : n) * Math.exp(-relDeps_().lgamma(n + 1)));
  }
  // P_k(j) = (1/k) Σ_{i=0}^{k−1} P_{k−1}(j−i), 누적합으로 O(c) 씩
  var prob = new Float64Array(c + 1), cum = new Float64Array(c + 1), k, j;
  prob[0] = 1;
  for (k = 2; k <= n; k++) {
    var s = 0;
    for (j = 0; j <= c; j++) { s += prob[j]; cum[j] = s; }
    for (j = 0; j <= c; j++) prob[j] = (cum[j] - (j - k >= 0 ? cum[j - k] : 0)) / k;
  }
  var cdf = 0;
  for (j = 0; j <= c; j++) cdf += prob[j];
  return Math.min(1, 2 * cdf);
}

/**
 * 켄달 타우-b 와 양측 유의확률 (n ≥ 3). 한 변수의 값이 모두 같으면 r 은 null.
 * method: 'auto'(기본, scipy 와 같은 규칙) | 'asymptotic' | 'exact'(동순위가 없을 때만, 시험용)
 * 반환: { n, r, p, method('exact' | 'asymptotic'), dis }
 */
function kendallTauB(x, y, method) {
  var d = relDeps_();
  var n = x.length;
  var res = { n: n, r: null, p: null, method: null, dis: null };
  if (n < 3) return res;
  var idx = new Array(n), i;
  for (i = 0; i < n; i++) idx[i] = i;
  idx.sort(function (a, b) { return x[a] - x[b] || y[a] - y[b]; });
  var xs = new Array(n), ys = new Array(n), joint = 0, run = 1;
  for (i = 0; i < n; i++) { xs[i] = x[idx[i]]; ys[i] = y[idx[i]]; }
  for (i = 1; i <= n; i++) {
    if (i < n && xs[i] === xs[i - 1] && ys[i] === ys[i - 1]) run++;
    else { joint += run * (run - 1) / 2; run = 1; }
  }
  var dis = relInversions_(ys);
  var tx = relTieStats_(xs), ty = relTieStats_(ys.slice().sort(function (a, b) { return a - b; }));
  var tot = n * (n - 1) / 2;
  var denom = (tot - tx.t) * (tot - ty.t);
  if (!(denom > 0)) return res;
  var cmd = tot - tx.t - ty.t + joint - 2 * dis;
  res.r = Math.max(-1, Math.min(1, cmd / Math.sqrt(denom)));
  res.dis = dis;
  var noTies = tx.t === 0 && ty.t === 0;
  var useExact = method === 'exact' || ((method === undefined || method === 'auto') && noTies && (n <= 33 || Math.min(dis, tot - dis) <= 1));
  if (useExact && !noTies) useExact = false;
  if (useExact) {
    res.p = relKendallExactP_(n, dis, tot);
    res.method = 'exact';
  } else {
    var m = n * (n - 1);
    var v = (m * (2 * n + 5) - tx.t1 - ty.t1) / 18 + 2 * tx.t * ty.t / m + (n > 2 ? tx.t0 * ty.t0 / (9 * m * (n - 2)) : 0);
    res.p = v > 0 ? 2 * d.pnorm(-Math.abs(cmd) / Math.sqrt(v)) : null;
    res.method = 'asymptotic';
  }
  return res;
}

/** 피어슨 r 의 Fisher z 신뢰구간 (n > 3, |r| < 1). conf: 신뢰수준(%) */
function pearsonCI(r, n, conf) {
  var out = { lo: null, hi: null };
  if (r === null || n <= 3 || Math.abs(r) >= 1) return out;
  var z = Math.atanh(r);
  var se = 1 / Math.sqrt(n - 3);
  var zc = relDeps_().qnorm((1 - conf / 100) / 2, false);
  out.lo = Math.tanh(z - zc * se);
  out.hi = Math.tanh(z + zc * se);
  return out;
}

/** 두 변수가 모두 유효한 사례의 값 쌍 */
function validPairs_(a, b) {
  var x = [], y = [];
  for (var i = 0; i < a.length; i++) {
    if (a[i] !== null && b[i] !== null) { x.push(a[i]); y.push(b[i]); }
  }
  return { x: x, y: y };
}

/* ---------- 편상관 ---------- */

var PART_MAX_CONTROLS_ = 10;

/** 대칭 양정치 행렬의 역행렬(Cholesky). 양정치가 아니면 null */
function relInvert_(A) {
  var n = A.length, i, j, k, L = [];
  for (i = 0; i < n; i++) L.push(new Array(n).fill(0));
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      var s = A[i][j];
      for (k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) { if (!(s > 1e-12)) return null; L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j];
    }
  }
  var inv = [];
  for (i = 0; i < n; i++) inv.push(new Array(n).fill(0));
  for (var col = 0; col < n; col++) {
    var y = new Array(n).fill(0), x = new Array(n).fill(0);
    for (i = 0; i < n; i++) { var s2 = i === col ? 1 : 0; for (k = 0; k < i; k++) s2 -= L[i][k] * y[k]; y[i] = s2 / L[i][i]; }
    for (i = n - 1; i >= 0; i--) { var s3 = y[i]; for (k = i + 1; k < n; k++) s3 -= L[k][i] * x[k]; x[i] = s3 / L[i][i]; }
    for (j = 0; j < n; j++) inv[j][col] = x[j];
  }
  return inv;
}

/**
 * 편상관. varCols: 분석 변수별 값 배열, ctrlCols: 통제변수별 값 배열 (모두 같은 길이, 결측 없음).
 * 반환: { n, c, df, r0 (통제 전 상관 행렬), r (편상관 행렬), p, lo, hi } 또는 { constant: 'var'|'ctrl', index } / { collinear: true }
 * 분산이 통제변수로 완전히 설명되는 변수가 있으면 그 변수의 편상관은 null.
 */
function partialCorrelation(varCols, ctrlCols, conf) {
  var d = relDeps_();
  var k = varCols.length, c = ctrlCols.length, n = varCols[0].length, m = k + c, i, j, t;
  var cols = varCols.concat(ctrlCols);
  var cen = [], sd = [];
  for (j = 0; j < m; j++) {
    var s = 0;
    for (t = 0; t < n; t++) s += cols[j][t];
    var mean = s / n, ss = 0, cc = new Array(n);
    for (t = 0; t < n; t++) { cc[t] = cols[j][t] - mean; ss += cc[t] * cc[t]; }
    if (!(ss > 0)) return { constant: j < k ? 'var' : 'ctrl', index: j < k ? j : j - k };
    var sq = Math.sqrt(ss);
    for (t = 0; t < n; t++) cc[t] /= sq;
    cen.push(cc);
  }
  var R = [];
  for (i = 0; i < m; i++) R.push(new Array(m).fill(0));
  for (i = 0; i < m; i++) {
    R[i][i] = 1;
    for (j = i + 1; j < m; j++) {
      var dot = 0;
      for (t = 0; t < n; t++) dot += cen[i][t] * cen[j][t];
      R[i][j] = dot; R[j][i] = dot;
    }
  }
  var Rcc = [];
  for (i = 0; i < c; i++) Rcc.push(R[k + i].slice(k));
  var inv = relInvert_(Rcc);
  if (!inv) return { collinear: true };
  var df = n - 2 - c;
  var Sig = [];
  for (i = 0; i < k; i++) {
    Sig.push(new Array(k));
    for (j = 0; j < k; j++) {
      var acc = R[i][j];
      for (var a = 0; a < c; a++) for (var b = 0; b < c; b++) acc -= R[i][k + a] * inv[a][b] * R[k + b][j];
      Sig[i][j] = acc;
    }
  }
  var out = { n: n, c: c, df: df, r0: [], r: [], p: [], lo: [], hi: [] };
  var zc = df >= 1 && n - c - 3 >= 1 ? d.qnorm((1 - conf / 100) / 2, false) : null;
  for (i = 0; i < k; i++) {
    out.r0.push(R[i].slice(0, k));
    out.r.push(new Array(k).fill(null)); out.p.push(new Array(k).fill(null)); out.lo.push(new Array(k).fill(null)); out.hi.push(new Array(k).fill(null));
  }
  for (i = 0; i < k; i++) {
    for (j = 0; j < k; j++) {
      if (i === j) { out.r[i][j] = 1; continue; }
      if (!(Sig[i][i] > 1e-12) || !(Sig[j][j] > 1e-12)) continue;
      var r = Math.max(-1, Math.min(1, Sig[i][j] / Math.sqrt(Sig[i][i] * Sig[j][j])));
      out.r[i][j] = r;
      if (df >= 1) out.p[i][j] = Math.abs(r) === 1 ? 0 : d.pValueT(r * Math.sqrt(df / (1 - r * r)), df);
      if (zc !== null && Math.abs(r) < 1) {
        var z = Math.atanh(r), se = 1 / Math.sqrt(n - c - 3);
        out.lo[i][j] = Math.tanh(z - zc * se); out.hi[i][j] = Math.tanh(z + zc * se);
      }
    }
  }
  return out;
}

function runPartialCorrelation_(dataset, request) {
  var d = relDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < 2) throw new Error('상관분석에는 변수를 2개 이상 선택하세요.');
  var ctrls = d.getVars(dataset, request.controls, { numeric: true });
  if (ctrls.length > PART_MAX_CONTROLS_) throw new Error('통제변수는 ' + PART_MAX_CONTROLS_ + '개까지 선택할 수 있습니다.');
  var seen = {};
  vars.forEach(function (v) { seen[v.name] = true; });
  ctrls.forEach(function (v) { if (seen[v.name]) throw new Error('변수 "' + v.name + '"를 분석 변수와 통제변수에 동시에 쓸 수 없습니다.'); });
  var opt = request.options || {};
  if (opt.kendall) throw new Error('켄달 타우-b는 통제변수(편상관)와 함께 쓸 수 없습니다. 켄달 옵션을 끄거나 통제변수를 비우세요.');
  var methods = CORR_METHODS_.filter(function (m) { return opt[m.key]; });
  if (methods.length === 0) throw new Error('상관계수 종류(피어슨, 스피어만)를 하나 이상 선택하세요.');
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');

  var names = vars.map(function (v) { return v.name; }).concat(ctrls.map(function (v) { return v.name; }));
  var sel = d.selectComplete(dataset, names);
  var n = sel.columns[names[0]].length, k = vars.length, c = ctrls.length;
  if (n - 2 - c < 1) throw new Error('모든 변수가 유효한 사례가 ' + n + '건이라 편상관을 계산할 수 없습니다 (통제변수 ' + c + '개일 때 ' + (c + 3) + '건 이상 필요). 결측이 많은 변수를 빼 보세요.');
  var warnings = [];
  if (sel.dropped > 0) warnings.push('목록별 결측 제거로 ' + sel.dropped + '건을 제외했습니다 (분석 변수와 통제변수 중 하나라도 비어 있는 사례). 사용한 사례는 ' + n + '건입니다.');
  var ctrlNames = ctrls.map(relName_).join(', ');
  var blocks = [];
  var scat = 0, scatSkipped = false;

  methods.forEach(function (m) {
    var rank = m.key === 'spearman';
    var vc = vars.map(function (v) { return rank ? rankAverage_(sel.columns[v.name]) : sel.columns[v.name]; });
    var cc = ctrls.map(function (v) { return rank ? rankAverage_(sel.columns[v.name]) : sel.columns[v.name]; });
    var res = partialCorrelation(vc, cc, conf);
    if (res.constant) throw new Error('변수 "' + relName_(res.constant === 'var' ? vars[res.index] : ctrls[res.index]) + '"의 값이 모두 같아 편상관을 계산할 수 없습니다.');
    if (res.collinear) throw new Error('통제변수들이 서로 선형결합이라(거의 같은 정보) 편상관을 계산할 수 없습니다. 같은 정보를 담은 통제변수를 하나 빼 보세요.');
    var label = rank ? '스피어만 순위 편상관' : '피어슨 편상관';
    var columns = [{ key: 'name', label: '' }, { key: 'stat', label: '' }];
    vars.forEach(function (v, idx) { columns.push({ key: 'c' + idx, label: relName_(v), format: { type: 'fixed', digits: 3 } }); });
    var rows = [], i, j;
    for (i = 0; i < k; i++) {
      var rRow = { _kind: 'group', _fmt: { type: 'fixed', digits: 3 }, name: relName_(vars[i]), stat: label + ' 계수' };
      var pRow = { _fmt: { type: 'pvalue' }, name: '', stat: '유의확률(양측)' };
      var dRow = { _fmt: { type: 'int' }, name: '', stat: '자유도' };
      for (j = 0; j < k; j++) {
        rRow['c' + j] = res.r[i][j] === null ? NaN : res.r[i][j];
        pRow['c' + j] = i === j ? null : (res.p[i][j] === null ? NaN : res.p[i][j]);
        dRow['c' + j] = i === j ? 0 : res.df;
      }
      rows.push(rRow, pRow, dRow);
    }
    blocks.push({
      type: 'table', caption: label + ' 행렬 (통제변수: ' + ctrlNames + ')', columns: columns, rows: rows,
      notes: [rank
        ? '통제변수의 영향을 뺀 순위 상관이다: 모든 변수를 순위(동순위는 평균순위)로 바꾼 뒤 편상관을 구했다. 유의확률은 t 근사(자유도 N − 2 − 통제변수 수 = ' + res.df + ')이다.'
        : '통제변수(' + ctrlNames + ')의 선형 영향을 두 변수에서 각각 뺀 뒤의 상관계수이다. 유의확률은 상관이 0이라는 귀무가설의 양측 t검정(자유도 N − 2 − 통제변수 수 = ' + res.df + ')이다.',
        '목록별 결측 제거: 분석 변수와 통제변수가 모두 유효한 ' + n + '건을 사용했다.']
    });
    var cRows = [];
    for (i = 0; i < k; i++) {
      for (j = i + 1; j < k; j++) {
        cRows.push({ a: relName_(vars[i]), b: relName_(vars[j]), r0: res.r0[i][j], r: res.r[i][j], diff: res.r[i][j] === null ? null : res.r[i][j] - res.r0[i][j],
          df: res.df, p: res.p[i][j], lo: res.lo[i][j], hi: res.hi[i][j] });
      }
    }
    var ccols = [
      { key: 'a', label: '변수1' }, { key: 'b', label: '변수2' },
      { key: 'r0', label: '통제 전 상관', format: { type: 'fixed', digits: 3 } }, { key: 'r', label: '편상관', format: { type: 'fixed', digits: 3 } },
      { key: 'diff', label: '변화(편상관 − 통제 전)', format: { type: 'fixed', digits: 3 } }, { key: 'df', label: '자유도', format: { type: 'int' } },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
    ];
    if (opt.ci) {
      ccols.push({ key: 'lo', label: conf + '% 신뢰구간 하한', format: { type: 'fixed', digits: 3 } }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: { type: 'fixed', digits: 3 } });
    }
    blocks.push({
      type: 'table', caption: label + ': 통제 전 → 통제 후 비교',
      columns: ccols, rows: cRows,
      notes: ['"통제 전 상관"은 같은 ' + n + '건으로 구한 단순 상관이다. 편상관이 통제 전보다 0에 가까워지면 통제변수가 두 변수의 관계를 상당 부분 설명(또는 매개)한다는 뜻이고, 부호가 바뀌거나 커지면 억제 효과이다. 인과의 증거는 아니다.']
        .concat(opt.ci ? ['신뢰구간은 Fisher z 변환(표준오차 1/√(N − 통제변수 수 − 3))으로 구했다.'] : [])
    });
    // 통제 후 잔차 산점도
    if (opt.scatter) {
      for (i = 0; i < k && !scatSkipped; i++) {
        for (j = i + 1; j < k; j++) {
          if (scat >= d.CHART_MAX_PAIRS_) { scatSkipped = true; break; }
          var preds = cc.map(function (col, q) { return { name: 'z' + q, values: col }; });
          var fx = d.olsFit(preds, vc[i], conf), fy = d.olsFit(preds, vc[j], conf);
          if (!fx.resid || !fy.resid) continue;
          blocks.push(d.scatterBlock((rank ? '순위의 ' : '') + '통제 후 잔차 산점도: ' + relName_(vars[i]) + ' × ' + relName_(vars[j]), fx.resid, fy.resid,
            relName_(vars[i]) + ' (통제 후 잔차)', relName_(vars[j]) + ' (통제 후 잔차)', res.r[i][j]));
          scat++;
        }
      }
    }
  });
  if (scatSkipped) warnings.push('변수 쌍이 많아 산점도는 처음 ' + d.CHART_MAX_PAIRS_ + '쌍만 그렸습니다.');
  return { analysis: 'correlation', title: '편상관분석 (통제변수: ' + ctrlNames + ')', blocks: blocks, warnings: warnings };
}

/* ---------- 상관분석 ---------- */

var CORR_METHODS_ = [
  { key: 'pearson', label: '피어슨 상관', fn: pearsonR },
  { key: 'spearman', label: '스피어만 순위상관', fn: spearmanR },
  { key: 'kendall', label: '켄달 타우-b 상관', fn: kendallTauB }
];

function runCorrelation(dataset, request) {
  if (request.controls && request.controls.length) return runPartialCorrelation_(dataset, request);
  var d = relDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < 2) throw new Error('상관분석에는 변수를 2개 이상 선택하세요.');
  var opt = request.options || {};
  var methods = CORR_METHODS_.filter(function (m) { return opt[m.key]; });
  if (methods.length === 0) throw new Error('상관계수 종류(피어슨, 스피어만, 켄달 타우-b)를 하나 이상 선택하세요.');
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');

  var names = vars.map(function (v) { return v.name; });
  var cols = {};
  var warnings = [];
  if (opt.listwise) {
    var sel = d.selectComplete(dataset, names);
    cols = sel.columns;
    if (sel.dropped > 0) warnings.push('목록별 결측 제거로 ' + sel.dropped + '건을 제외했습니다.');
  } else {
    names.forEach(function (nm) { cols[nm] = dataset.data[nm]; });
  }

  var k = vars.length;
  var blocks = [];
  var ciRows = [];

  methods.forEach(function (m) {
    // 쌍별 결과 행렬
    var res = [];
    for (var i = 0; i < k; i++) res.push(new Array(k));
    for (i = 0; i < k; i++) {
      for (var j = i; j < k; j++) {
        var pr = validPairs_(cols[names[i]], cols[names[j]]);
        var r;
        if (i === j) {
          r = { n: pr.x.length, r: 1, p: null };
        } else {
          r = m.fn(pr.x, pr.y);
          if (r.r === null && pr.x.length >= 3) {
            warnings.push('"' + relName_(vars[i]) + '"와 "' + relName_(vars[j]) + '": 한 변수의 값이 모두 같아 상관계수를 계산할 수 없습니다.');
          }
        }
        res[i][j] = res[j][i] = r;
        if (m.key === 'pearson' && opt.ci && i < j) {
          var ci = pearsonCI(r.r, r.n, conf);
          ciRows.push({ a: relName_(vars[i]), b: relName_(vars[j]), n: r.n, r: r.r, lo: ci.lo, hi: ci.hi });
        }
      }
    }

    var columns = [{ key: 'name', label: '' }, { key: 'stat', label: '' }];
    vars.forEach(function (v, idx) {
      columns.push({ key: 'c' + idx, label: relName_(v), format: { type: 'fixed', digits: 3 } });
    });
    var rows = [];
    for (i = 0; i < k; i++) {
      var rRow = { _kind: 'group', _fmt: { type: 'fixed', digits: 3 }, name: relName_(vars[i]), stat: m.key === 'pearson' ? 'Pearson 상관계수' : m.key === 'spearman' ? 'Spearman 상관계수' : 'Kendall 타우-b' };
      var pRow = { _fmt: { type: 'pvalue' }, name: '', stat: '유의확률(양측)' };
      var nRow = { _fmt: { type: 'int' }, name: '', stat: 'N' };
      for (j = 0; j < k; j++) {
        var c = res[i][j];
        rRow['c' + j] = c.r === null ? NaN : c.r;
        pRow['c' + j] = i === j ? null : (c.r === null ? NaN : c.p);
        nRow['c' + j] = c.n;
      }
      rows.push(rRow, pRow, nRow);
    }
    var exactUsed = false;
    for (i = 0; i < k; i++) for (j = i + 1; j < k; j++) if (res[i][j].method === 'exact') exactUsed = true;
    var notes = [m.key === 'pearson'
      ? '피어슨 상관계수: 두 변수의 선형 관계의 세기(-1 ~ 1). 유의확률은 상관이 0이라는 귀무가설의 양측 t검정(자유도 N−2)이다.'
      : m.key === 'spearman'
        ? '스피어만 상관계수: 순위(동순위는 평균순위)의 피어슨 상관이다. 유의확률은 t 근사(자유도 N−2)이다.'
        : '켄달 타우-b: 동순위를 보정한 순위상관으로, 같은 방향(일치)인 쌍과 반대 방향(불일치)인 쌍의 차이를 전체 쌍으로 나눈 값이다 (-1 ~ 1). 보통 스피어만보다 절댓값이 작다. 표본이 작거나 동순위가 많은 서열 자료에서 스피어만의 대안으로 쓴다. 유의확률은 동순위 보정 분산의 정규 근사(양측)이고' +
          (exactUsed ? ', 동순위가 없고 N이 33 이하인 쌍은 정확 분포로 구했다.' : ' (동순위가 없고 N이 33 이하이면 정확 분포를 쓴다).')];
    notes.push(opt.listwise
      ? '목록별 결측 제거: 선택한 모든 변수가 유효한 사례만 사용했다.'
      : '쌍별 결측 제거: 각 변수 쌍에서 두 변수가 모두 유효한 사례를 사용했다 (N은 쌍마다 다를 수 있다).');
    blocks.push({ type: 'table', caption: m.label + ' 행렬', columns: columns, rows: rows, notes: notes });
  });

  if (ciRows.length) {
    blocks.push({
      type: 'table',
      caption: '피어슨 상관계수의 ' + conf + '% 신뢰구간 (Fisher z 변환)',
      columns: [
        { key: 'a', label: '변수1' }, { key: 'b', label: '변수2' },
        { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'r', label: '상관계수', format: { type: 'fixed', digits: 3 } },
        { key: 'lo', label: '하한', format: { type: 'fixed', digits: 3 } },
        { key: 'hi', label: '상한', format: { type: 'fixed', digits: 3 } }
      ],
      rows: ciRows,
      notes: ['N이 4 미만이거나 상관계수가 ±1이면 구간을 계산하지 않는다.']
    });
  }

  if (opt.assumptions && opt.pearson) {
    var advVars = vars.map(function (v) {
      return { name: relName_(v), values: cols[v.name].filter(function (x) { return x !== null; }), ordinal: v.scale === 'ordinal' };
    });
    var advPart = d.advCorrelation_(advVars, !!opt.spearman);
    blocks = d.advAppend_(blocks, { rows: advPart.rows, recs: advPart.recs }, false);
  }

  if (opt.scatter) {
    var made = 0, skipped = false;
    for (var a1 = 0; a1 < k && !skipped; a1++) {
      for (var b1 = a1 + 1; b1 < k; b1++) {
        if (made >= d.CHART_MAX_PAIRS_) { skipped = true; break; }
        var pr = validPairs_(cols[names[a1]], cols[names[b1]]);
        if (pr.x.length < 2) continue;
        blocks.push(d.scatterBlock('산점도: ' + relName_(vars[a1]) + ' × ' + relName_(vars[b1]), pr.x, pr.y,
          relName_(vars[a1]), relName_(vars[b1]), pearsonR(pr.x, pr.y).r));
        made++;
      }
    }
    if (skipped) warnings.push('변수 쌍이 많아 산점도는 처음 ' + d.CHART_MAX_PAIRS_ + '쌍만 그렸습니다.');
  }

  return { analysis: 'correlation', title: '상관분석', blocks: blocks, warnings: warnings };
}

/* ---------- 교차분석 ---------- */

var CROSSTAB_MAX_CATEGORIES_ = 30;

/**
 * 2×2 표 [[a, b], [c, d]] 의 Fisher 정확검정.
 * two: 양측 (관측 표와 같거나 더 작은 확률의 표를 모두 더함, scipy와 같은 허용오차 1e-7)
 * one: 단측 (관측값 쪽의 더 작은 꼬리 확률)
 */
function fisherExact2x2(a, b, c, d) {
  var lg = relDeps_().lgamma;
  var r1 = a + b, r2 = c + d, c1 = a + c, N = r1 + r2;
  var lo = Math.max(0, c1 - r2), hi = Math.min(r1, c1);
  var constant = lg(r1 + 1) + lg(r2 + 1) + lg(c1 + 1) + lg(N - c1 + 1) - lg(N + 1);
  function logP(x) {
    return constant - lg(x + 1) - lg(r1 - x + 1) - lg(c1 - x + 1) - lg(r2 - c1 + x + 1);
  }
  var obs = logP(a), tol = Math.log1p(1e-7);
  var two = 0, le = 0, ge = 0;
  for (var x = lo; x <= hi; x++) {
    var lp = logP(x), p = Math.exp(lp);
    if (x <= a) le += p;
    if (x >= a) ge += p;
    if (lp <= obs + tol) two += p;
  }
  return { two: Math.min(1, two), one: Math.min(1, le, ge) };
}

/** 범주(값 순)와 레이블: [{ value, label }] */
function categories_(values, variable) {
  var d = relDeps_();
  var seen = {}, cats = [];
  values.forEach(function (v) {
    var key = typeof v + ':' + v;
    if (!seen[key]) { seen[key] = true; cats.push({ value: v }); }
  });
  cats.sort(d.compareValues_);
  cats.forEach(function (c) {
    var label = variable.valueLabels[String(c.value)];
    c.label = label !== undefined ? label : String(c.value);
  });
  return cats;
}

function runCrosstab(dataset, request) {
  var d = relDeps_();
  var rv = d.getVars(dataset, request.rowvar);
  var cv = d.getVars(dataset, request.colvar);
  if (rv.length !== 1 || cv.length !== 1) throw new Error('행 변수와 열 변수에 변수를 하나씩 선택하세요.');
  var rowVar = rv[0], colVar = cv[0];
  if (rowVar.name === colVar.name) throw new Error('행 변수와 열 변수에 같은 변수 "' + rowVar.name + '"를 쓸 수 없습니다.');
  var opt = request.options || {};
  var rowName = relName_(rowVar), colName = relName_(colVar);

  var rowVals = [], colVals = [], i, j;
  for (i = 0; i < dataset.n; i++) {
    var a = dataset.data[rowVar.name][i], b = dataset.data[colVar.name][i];
    if (a !== null && b !== null) { rowVals.push(a); colVals.push(b); }
  }
  var N = rowVals.length;
  if (N === 0) throw new Error('두 변수가 모두 유효한 사례가 없습니다.');
  var rc = categories_(rowVals, rowVar), cc = categories_(colVals, colVar);
  var nr = rc.length, nc = cc.length;
  if (nr < 2 || nc < 2) {
    throw new Error('교차분석에는 각 변수에 범주가 2개 이상 필요합니다 (행 변수 "' + rowName + '" ' + nr + '개, 열 변수 "' + colName + '" ' + nc + '개).');
  }
  if (nr > CROSSTAB_MAX_CATEGORIES_ || nc > CROSSTAB_MAX_CATEGORIES_) {
    throw new Error('범주가 너무 많습니다 (행 ' + nr + '개, 열 ' + nc + '개, 최대 ' + CROSSTAB_MAX_CATEGORIES_ + '개). 범주가 적은 변수를 선택하세요.');
  }

  // 빈도표
  var rowIndex = {}, colIndex = {};
  rc.forEach(function (c, k) { rowIndex[typeof c.value + ':' + c.value] = k; });
  cc.forEach(function (c, k) { colIndex[typeof c.value + ':' + c.value] = k; });
  var obs = [], R = new Array(nr), C = new Array(nc);
  for (i = 0; i < nr; i++) {
    obs.push(new Array(nc));
    R[i] = 0;
    for (j = 0; j < nc; j++) obs[i][j] = 0;
  }
  for (j = 0; j < nc; j++) C[j] = 0;
  for (var t = 0; t < N; t++) {
    var ri = rowIndex[typeof rowVals[t] + ':' + rowVals[t]], ci = colIndex[typeof colVals[t] + ':' + colVals[t]];
    obs[ri][ci]++; R[ri]++; C[ci]++;
  }

  // 기대빈도와 통계량
  var exp = [], adj = [], chi2 = 0, lr = 0, yates = 0, low = 0, minExp = Infinity;
  for (i = 0; i < nr; i++) {
    exp.push([]); adj.push([]);
    for (j = 0; j < nc; j++) {
      var e = R[i] * C[j] / N, o = obs[i][j];
      exp[i].push(e);
      chi2 += (o - e) * (o - e) / e;
      if (o > 0) lr += o * Math.log(o / e);
      var yd = Math.max(0, Math.abs(o - e) - 0.5);
      yates += yd * yd / e;
      adj[i].push((o - e) / Math.sqrt(e * (1 - R[i] / N) * (1 - C[j] / N)));
      if (e < 5) low++;
      if (e < minExp) minExp = e;
    }
  }
  lr *= 2;
  var df = (nr - 1) * (nc - 1);
  var is2x2 = nr === 2 && nc === 2;
  var warnings = [];
  var lowPct = low / (nr * nc) * 100;
  if (lowPct > 20) {
    warnings.push('기대빈도가 5 미만인 셀이 ' + lowPct.toFixed(1) + '%입니다. 카이제곱 근사가 부정확할 수 있습니다' +
      (is2x2 ? ' (Fisher 정확검정을 참고하세요).' : ' (범주를 합치는 것을 고려하세요).'));
  }

  /* --- 교차표 --- */
  var statDefs = [{ key: 'count', label: '빈도', fmt: { type: 'int' } }];
  if (opt.expected) statDefs.push({ key: 'exp', label: '기대빈도', fmt: { type: 'fixed', digits: 1 } });
  if (opt.rowPct) statDefs.push({ key: 'row', label: '행 %', fmt: { type: 'fixed', digits: 1 } });
  if (opt.colPct) statDefs.push({ key: 'col', label: '열 %', fmt: { type: 'fixed', digits: 1 } });
  if (opt.totalPct) statDefs.push({ key: 'tot', label: '전체 %', fmt: { type: 'fixed', digits: 1 } });
  if (opt.residual) statDefs.push({ key: 'res', label: '조정된 표준화 잔차', fmt: { type: 'fixed', digits: 1 } });

  var columns = [{ key: 'name', label: rowName }, { key: 'stat', label: '' }];
  cc.forEach(function (c, k) { columns.push({ key: 'c' + k, label: c.label, format: { type: 'int' } }); });
  columns.push({ key: 'total', label: '합계', format: { type: 'int' } });

  /** 셀 (ri, cj) 의 값. ri 또는 cj 가 -1 이면 그쪽의 합계. */
  function cellValue(key, ri, cj) {
    var o = ri < 0 ? (cj < 0 ? N : C[cj]) : (cj < 0 ? R[ri] : obs[ri][cj]);
    var rowBase = ri < 0 ? N : R[ri];
    var colBase = cj < 0 ? N : C[cj];
    switch (key) {
      case 'count': return o;
      case 'exp': return ri >= 0 && cj >= 0 ? exp[ri][cj] : o; // 합계의 기대빈도는 관측 합계와 같다
      case 'row': return o / rowBase * 100;
      case 'col': return o / colBase * 100;
      case 'tot': return o / N * 100;
      case 'res': return ri >= 0 && cj >= 0 ? adj[ri][cj] : null;
      default: return null;
    }
  }

  var rows = [];
  var groups = rc.map(function (c, k) { return { label: c.label, ri: k }; });
  groups.push({ label: '합계', ri: -1 });
  groups.forEach(function (g) {
    statDefs.forEach(function (sd, si) {
      var row = { name: si === 0 ? g.label : '', stat: sd.label, _fmt: sd.fmt };
      if (si === 0) row._kind = g.ri < 0 ? 'total' : 'group';
      for (var cj = -1; cj < nc; cj++) row[cj < 0 ? 'total' : 'c' + cj] = cellValue(sd.key, g.ri, cj);
      rows.push(row);
    });
  });

  var blocks = [{
    type: 'table',
    caption: '교차표: ' + rowName + ' × ' + colName + ' (열 변수: ' + colName + ')',
    columns: columns,
    rows: rows,
    notes: ['두 변수가 모두 유효한 ' + N + '건을 사용했다' + (dataset.n > N ? ' (' + (dataset.n - N) + '건 제외).' : '.')]
  }];

  /* --- 카이제곱 검정 --- */
  var pAsym = function (v, dof) { return d.pchisq(v, dof, false); };
  var testRows = [];
  testRows.push({ test: 'Pearson 카이제곱', value: chi2, df: df, p: pAsym(chi2, df), p2: null, p1: null });
  if (is2x2) testRows.push({ test: '연속성 수정(Yates)', value: yates, df: 1, p: pAsym(yates, 1), p2: null, p1: null });
  testRows.push({ test: '우도비', value: lr, df: df, p: pAsym(lr, df), p2: null, p1: null });
  if (is2x2) {
    var fe = fisherExact2x2(obs[0][0], obs[0][1], obs[1][0], obs[1][1]);
    testRows.push({ test: 'Fisher의 정확검정', value: null, df: null, p: null, p2: fe.two, p1: fe.one });
  }
  var numeric = rc.every(function (c) { return typeof c.value === 'number'; }) &&
    cc.every(function (c) { return typeof c.value === 'number'; });
  if (numeric) {
    var lin = pearsonR(rowVals, colVals);
    if (lin.r !== null) {
      var m2 = (N - 1) * lin.r * lin.r;
      testRows.push({ test: '선형 대 선형 연관', value: m2, df: 1, p: pAsym(m2, 1), p2: null, p1: null });
    }
  }
  testRows.push({ _kind: 'total', _fmt: { type: 'int' }, test: '유효 케이스 수', value: N, df: null, p: null, p2: null, p1: null });

  blocks.push({
    type: 'table',
    caption: '카이제곱 검정',
    columns: [
      { key: 'test', label: '' },
      { key: 'value', label: '값', format: { type: 'fixed', digits: 3 } },
      { key: 'df', label: '자유도', format: { type: 'int' } },
      { key: 'p', label: '점근 유의확률(양측)', format: { type: 'pvalue' } },
      { key: 'p2', label: '정확 유의확률(양측)', format: { type: 'pvalue' } },
      { key: 'p1', label: '정확 유의확률(단측)', format: { type: 'pvalue' } }
    ],
    rows: testRows,
    notes: [
      '기대빈도가 5 미만인 셀: ' + low + '개(' + lowPct.toFixed(1) + '%), 최소 기대빈도 ' + minExp.toFixed(2) + '.',
      '귀무가설: 두 변수는 서로 독립이다. 유의확률이 작으면(예: .05 미만) 독립성을 기각한다.' +
        (numeric ? ' 선형 대 선형 연관은 범주의 숫자 코드를 점수로 쓴다.' : '')
    ]
  });

  /* --- 효과크기 --- */
  if (opt.effect) {
    var effRows = [];
    if (is2x2) {
      var sign = obs[0][0] * obs[1][1] - obs[0][1] * obs[1][0] >= 0 ? 1 : -1;
      effRows.push({ name: '파이', value: sign * Math.sqrt(chi2 / N), p: pAsym(chi2, df) });
    }
    effRows.push({ name: 'Cramér의 V', value: Math.sqrt(chi2 / (N * Math.min(nr - 1, nc - 1))), p: pAsym(chi2, df) });
    effRows.push({ name: '분할계수', value: Math.sqrt(chi2 / (chi2 + N)), p: pAsym(chi2, df) });
    blocks.push({
      type: 'table',
      caption: '대칭적 측도',
      columns: [
        { key: 'name', label: '' },
        { key: 'value', label: '값', format: { type: 'fixed', digits: 3 } },
        { key: 'p', label: '근사 유의확률', format: { type: 'pvalue' } }
      ],
      rows: effRows,
      notes: ['파이는 2×2 표에서만 계산하며 부호는 연관의 방향이다. 근사 유의확률은 Pearson 카이제곱의 값과 같다.']
    });
  }

  return { analysis: 'crosstab', title: '교차분석', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = {
    pearsonR: pearsonR,
    spearmanR: spearmanR,
    kendallTauB: kendallTauB,
    rankAverage_: rankAverage_,
    pearsonCI: pearsonCI,
    partialCorrelation: partialCorrelation,
    fisherExact2x2: fisherExact2x2,
    runCorrelation: runCorrelation,
    runCrosstab: runCrosstab
  };
}
