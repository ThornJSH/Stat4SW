/**
 * 평가자 일치도: Cohen kappa(평가자 2명, 가중 선택), Fleiss kappa(3명 이상), 급내상관계수 ICC(6유형).
 * 순수 함수: (dataset, request) → 결과 객체. 평가자(또는 측정)마다 변수 한 열, 사례(대상)가 행이다. 하나라도 비어 있는 사례는 통째로 제외한다.
 *
 * - Cohen kappa: κ = (p_o − p_e)/(1 − p_e) (가중: κ = 1 − Σ d_ij p_ij / Σ d_ij p_i· p_·j, d = 선형 |i−j|/(r−1) 또는 이차 ((i−j)/(r−1))²).
 *   표준오차는 Fleiss·Cohen·Everitt(1969)의 비영가설 식(신뢰구간용), 영가설 표준오차로 z = κ/SE₀ 검정.
 * - Fleiss kappa: P̄ = 사례별 일치도의 평균, P̄e = Σ p_j², κ = (P̄ − P̄e)/(1 − P̄e). 범주별 κ_j = 1 − Σ x_ij(n − x_ij)/(N n (n − 1) p_j q_j).
 *   영가설의 표준오차는 Fleiss(1971)식 (비교할 다른 도구가 없어 식만 구현).
 * - ICC (McGraw & Wong 1996, Shrout & Fleiss 1979): 이원 분산분석의 평균제곱(대상 MSB, 평가자 MSJ, 잔차 MSE, 대상 안 MSW)으로
 *   ICC(1)=(MSB−MSW)/(MSB+(k−1)MSW), ICC(2)=(MSB−MSE)/(MSB+(k−1)MSE+k(MSJ−MSE)/n), ICC(3)=(MSB−MSE)/(MSB+(k−1)MSE), 평균 측정은 ICC(1,k)=(MSB−MSW)/MSB,
 *   ICC(2,k)=(MSB−MSE)/(MSB+(MSJ−MSE)/n), ICC(3,k)=(MSB−MSE)/MSB. 유의확률은 F(n−1, 자유도), 신뢰구간은 F 분포 식 (pingouin과 같은 식, 신뢰수준만 일반화).
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function agrDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return { getVars: data.getVars, compareValues_: data.compareValues_, pf: dist.pf, qf: dist.qf, qnorm: dist.qnorm, pnorm: dist.pnorm };
  }
  return { getVars: getVars, compareValues_: compareValues_, pf: pf, qf: qf, qnorm: qnorm, pnorm: pnorm };
}

var AGR_MAX_CATS_ = 20;
var AGR_MAX_RATERS_ = 20;
var AGR_F3_ = { type: 'fixed', digits: 3 };
function agrName_(v) { return v.label ? v.label : v.name; }

/* ---------- Cohen kappa ---------- */

/** 가중 거리 행렬 d (0 이면 완전 일치). kind: 'none' | 'linear' | 'quadratic' */
function agrDistance_(r, kind) {
  var d = [], i, j;
  for (i = 0; i < r; i++) {
    d.push([]);
    for (j = 0; j < r; j++) {
      var x = r > 1 ? Math.abs(i - j) / (r - 1) : 0;
      d[i].push(kind === 'linear' ? x : kind === 'quadratic' ? x * x : (i === j ? 0 : 1));
    }
  }
  return d;
}

/**
 * 2명 평가자의 교차표 T (r × r 개수)에서 Cohen kappa. weights: 'none' | 'linear' | 'quadratic', conf: 신뢰수준(%)
 * 반환: { n, po, pe, kappa, se, se0, z, p, lo, hi, kmax }  (pe = 1 이면 kappa 등은 null)
 */
function cohenKappa(T, weights, conf) {
  var dd = agrDeps_();
  var r = T.length, i, j, N = 0;
  T.forEach(function (row) { row.forEach(function (v) { N += v; }); });
  var P = T.map(function (row) { return row.map(function (v) { return v / N; }); });
  var row = new Array(r).fill(0), col = new Array(r).fill(0);
  for (i = 0; i < r; i++) for (j = 0; j < r; j++) { row[i] += P[i][j]; col[j] += P[i][j]; }
  var d = agrDistance_(r, weights);
  var obsD = 0, expD = 0;
  for (i = 0; i < r; i++) for (j = 0; j < r; j++) { obsD += d[i][j] * P[i][j]; expD += d[i][j] * row[i] * col[j]; }
  var res = { n: N, po: 1 - obsD, pe: 1 - expD, kappa: null, se: null, se0: null, z: null, p: null, lo: null, hi: null, kmax: null };
  if (!(expD > 0)) return res;
  var kappa = 1 - obsD / expD;
  res.kappa = kappa;
  // 표준오차 (가중치 w = 1 − d)
  var w = d.map(function (rw) { return rw.map(function (v) { return 1 - v; }); });
  var wRow = new Array(r).fill(0), wCol = new Array(r).fill(0), wExp = 0;
  for (i = 0; i < r; i++) for (j = 0; j < r; j++) { wRow[i] += col[j] * w[i][j]; wCol[j] += row[i] * w[i][j]; wExp += w[i][j] * row[i] * col[j]; }
  var fac = 1 / ((1 - wExp) * (1 - wExp) * N), v1, v0, a, b;
  if (weights === 'none') {
    var ta = 0, tb = 0, tc0 = 0;
    for (i = 0; i < r; i++) {
      ta += P[i][i] * Math.pow(1 - (row[i] + col[i]) * (1 - kappa), 2);
      tc0 += row[i] * col[i] * (row[i] + col[i]);
      for (j = 0; j < r; j++) if (i !== j) tb += P[i][j] * Math.pow(col[i] + row[j], 2);   // P[i][j]: i행 j열, (i열 비율 + j행 비율)
    }
    tb *= (1 - kappa) * (1 - kappa);
    v1 = (ta + tb - Math.pow(kappa - wExp * (1 - kappa), 2)) * fac;
    v0 = (wExp + wExp * wExp - tc0) * fac;
  } else {
    var s1 = 0, s0 = 0;
    for (a = 0; a < r; a++) {
      for (b = 0; b < r; b++) {
        s1 += P[a][b] * Math.pow(w[a][b] - (wCol[b] + wRow[a]) * (1 - kappa), 2);
        s0 += row[a] * col[b] * Math.pow(w[a][b] - (wCol[b] + wRow[a]), 2);
      }
    }
    v1 = (s1 - Math.pow(kappa - wExp * (1 - kappa), 2)) * fac;
    v0 = (s0 - wExp * wExp) * fac;
  }
  if (v1 > 0) res.se = Math.sqrt(v1);
  if (v0 > 0) res.se0 = Math.sqrt(v0);
  var zc = dd.qnorm((1 - conf / 100) / 2, false);
  if (res.se !== null) { res.lo = kappa - zc * res.se; res.hi = kappa + zc * res.se; }
  if (res.se0 !== null) { res.z = kappa / res.se0; res.p = 2 * dd.pnorm(-Math.abs(res.z), true); }
  var pe0 = 0, mn = 0;
  for (i = 0; i < r; i++) { pe0 += row[i] * col[i]; mn += Math.min(row[i], col[i]); }
  if (pe0 < 1) res.kmax = (mn - pe0) / (1 - pe0);
  return res;
}

/**
 * Fleiss kappa. X: 사례 × 범주 개수 (각 행의 합 = 평가자 수 n). 반환:
 * { N, n, p: [범주 비율], counts, pbar, pe, kappa, se, z, pval, kappaJ, seJ, zJ, pJ }
 */
function fleissKappa(X) {
  var dd = agrDeps_();
  var N = X.length, r = X[0].length, n = 0, i, j;
  X[0].forEach(function (v) { n += v; });
  var counts = new Array(r).fill(0), pbar = 0;
  for (i = 0; i < N; i++) {
    var s2 = 0;
    for (j = 0; j < r; j++) { counts[j] += X[i][j]; s2 += X[i][j] * X[i][j]; }
    pbar += (s2 - n) / (n * (n - 1));
  }
  pbar /= N;
  var p = counts.map(function (c) { return c / (N * n); });
  var pe = 0, sp = 0, spq = 0;
  p.forEach(function (v) { pe += v * v; sp += v * (1 - v); });
  p.forEach(function (v) { spq += v * (1 - v) * ((1 - v) - v); });
  var res = { N: N, n: n, p: p, counts: counts, pbar: pbar, pe: pe, kappa: null, se: null, z: null, pval: null, kappaJ: [], seJ: null, zJ: [], pJ: [] };
  if (!(1 - pe > 0)) return res;
  res.kappa = (pbar - pe) / (1 - pe);
  var base = 2 / (N * n * (n - 1));
  var inner = sp * sp - spq;
  if (inner > 0 && sp > 0) {
    res.se = Math.sqrt(base) * Math.sqrt(inner) / sp;
    res.z = res.kappa / res.se;
    res.pval = 2 * dd.pnorm(-Math.abs(res.z), true);
  }
  res.seJ = Math.sqrt(base);
  for (j = 0; j < r; j++) {
    var pq = p[j] * (1 - p[j]);
    if (!(pq > 0)) { res.kappaJ.push(null); res.zJ.push(null); res.pJ.push(null); continue; }
    var sx = 0;
    for (i = 0; i < N; i++) sx += X[i][j] * (n - X[i][j]);
    var kj = 1 - sx / (N * n * (n - 1) * pq);
    res.kappaJ.push(kj);
    res.zJ.push(kj / res.seJ);
    res.pJ.push(2 * dd.pnorm(-Math.abs(kj / res.seJ), true));
  }
  return res;
}

/* ---------- ICC ---------- */

var AGR_ICC_TYPES_ = [
  { key: 'icc1', label: 'ICC(1,1)', desc: '일원 확률 · 단일 측정', mean: false },
  { key: 'icc2', label: 'ICC(2,1)', desc: '이원 확률 · 절대 일치도 · 단일 측정', mean: false },
  { key: 'icc3', label: 'ICC(3,1)', desc: '이원 혼합 · 일관성 · 단일 측정', mean: false },
  { key: 'icc1k', label: 'ICC(1,k)', desc: '일원 확률 · 평균 측정', mean: true },
  { key: 'icc2k', label: 'ICC(2,k)', desc: '이원 확률 · 절대 일치도 · 평균 측정', mean: true },
  { key: 'icc3k', label: 'ICC(3,k)', desc: '이원 혼합 · 일관성 · 평균 측정', mean: true }
];

/**
 * columns: 평가자(측정)별 값 배열 (같은 길이 n, 결측 없음), conf: 신뢰수준(%)
 * 반환: { n, k, ssb, ssj, sse, sst, dfb, dfj, dfe, msb, msj, mse, msw, means, sds, rows: [{ key, icc, F, df1, df2, p, lo, hi }] }
 * 변동이 없으면(msb = mse = 0 등) 값이 null.
 */
function iccCompute(columns, conf) {
  var d = agrDeps_();
  var k = columns.length, n = columns[0].length, i, j;
  var grand = 0, subj = new Array(n).fill(0), means = [], sds = [];
  for (j = 0; j < k; j++) {
    var s = 0;
    for (i = 0; i < n; i++) { s += columns[j][i]; subj[i] += columns[j][i]; }
    means.push(s / n); grand += s;
  }
  grand /= n * k;
  var ssb = 0, ssj = 0, sst = 0;
  for (i = 0; i < n; i++) ssb += k * (subj[i] / k - grand) * (subj[i] / k - grand);
  for (j = 0; j < k; j++) {
    ssj += n * (means[j] - grand) * (means[j] - grand);
    var ss = 0;
    for (i = 0; i < n; i++) { sst += (columns[j][i] - grand) * (columns[j][i] - grand); ss += (columns[j][i] - means[j]) * (columns[j][i] - means[j]); }
    sds.push(n > 1 ? Math.sqrt(ss / (n - 1)) : null);
  }
  var sse = Math.max(0, sst - ssb - ssj);
  var dfb = n - 1, dfj = k - 1, dfe = (n - 1) * (k - 1);
  var msb = ssb / dfb, msj = ssj / dfj, mse = sse / dfe, msw = (ssj + sse) / (dfj + dfe);
  var r = { n: n, k: k, ssb: ssb, ssj: ssj, sse: sse, sst: sst, dfb: dfb, dfj: dfj, dfe: dfe, msb: msb, msj: msj, mse: mse, msw: msw, means: means, sds: sds, rows: [] };
  var alpha = 1 - conf / 100;
  var ppf = function (q, a, b) { return d.qf(1 - q, a, b, false); };      // F 분포의 q 분위수
  var empty = function (key) { return { key: key, icc: null, F: null, df1: dfb, df2: null, p: null, lo: null, hi: null }; };
  var ok = msb > 0 && msw > 0 && mse > 0 && isFinite(msb + msw + mse);
  if (!ok) { AGR_ICC_TYPES_.forEach(function (t) { r.rows.push(empty(t.key)); }); return r; }
  var icc1 = (msb - msw) / (msb + (k - 1) * msw);
  var icc2 = (msb - mse) / (msb + (k - 1) * mse + k * (msj - mse) / n);
  var icc3 = (msb - mse) / (msb + (k - 1) * mse);
  var icc1k = (msb - msw) / msb, icc2k = (msb - mse) / (msb + (msj - mse) / n), icc3k = (msb - mse) / msb;
  var F1 = msb / msw, F2 = msb / mse, df1kd = n * (k - 1), df2kd = dfe;
  var q = 1 - alpha / 2;
  var f1l = F1 / ppf(q, dfb, df1kd), f1u = F1 * ppf(q, df1kd, dfb);
  var l1 = (f1l - 1) / (f1l + (k - 1)), u1 = (f1u - 1) / (f1u + (k - 1));
  var f3l = F2 / ppf(q, dfb, df2kd), f3u = F2 * ppf(q, df2kd, dfb);
  var l3 = (f3l - 1) / (f3l + (k - 1)), u3 = (f3u - 1) / (f3u + (k - 1));
  var fj = msj / mse;
  var vn = df2kd * Math.pow(k * icc2 * fj + n * (1 + (k - 1) * icc2) - k * icc2, 2);
  var vd = dfb * k * k * icc2 * icc2 * fj * fj + Math.pow(n * (1 + (k - 1) * icc2) - k * icc2, 2);
  var v = vn / vd;
  var f2u = ppf(q, n - 1, v), f2l = ppf(q, v, n - 1);
  var l2 = n * (msb - f2u * mse) / (f2u * (k * msj + (k * n - k - n) * mse) + n * msb);
  var u2 = n * (f2l * msb - mse) / (k * msj + (k * n - k - n) * mse + n * f2l * msb);
  var mk = function (key, icc, F, df2, lo, hi) { return { key: key, icc: icc, F: F, df1: dfb, df2: df2, p: d.pf(F, dfb, df2, false), lo: lo, hi: hi }; };
  r.rows = [
    mk('icc1', icc1, F1, df1kd, l1, u1), mk('icc2', icc2, F2, df2kd, l2, u2), mk('icc3', icc3, F2, df2kd, l3, u3),
    mk('icc1k', icc1k, F1, df1kd, 1 - 1 / f1l, 1 - 1 / f1u),
    mk('icc2k', icc2k, F2, df2kd, l2 * k / (1 + l2 * (k - 1)), u2 * k / (1 + u2 * (k - 1))),
    mk('icc3k', icc3k, F2, df2kd, 1 - 1 / f3l, 1 - 1 / f3u)
  ];
  return r;
}

/* ---------- 실행 ---------- */

function agrConf_(opt) {
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  return conf;
}

/** 목록별 결측 제거 (평가자 모두에 값이 있는 사례) */
function agrComplete_(dataset, vars) {
  var rows = [];
  for (var i = 0; i < dataset.n; i++) {
    var ok = true;
    for (var a = 0; a < vars.length; a++) if (dataset.data[vars[a].name][i] === null) { ok = false; break; }
    if (ok) rows.push(i);
  }
  return rows;
}

function agrPstr_(p) { return p === null ? '.' : (p < 0.001 ? '<.001' : p.toFixed(3).replace(/^0/, '')); }

/** Landis & Koch 기준을 쉬운 말로 */
function agrKappaWord_(k) {
  if (k === null) return '계산 불가';
  if (k < 0) return '우연보다 낮음';
  if (k <= 0.2) return '매우 낮음';
  if (k <= 0.4) return '낮음';
  if (k <= 0.6) return '보통';
  if (k <= 0.8) return '높음';
  return '거의 완벽';
}

function runKappa(dataset, request) {
  var d = agrDeps_();
  var vars = d.getVars(dataset, request.vars);
  if (vars.length < 2) throw new Error('kappa 는 평가자 변수를 2개 이상 선택하세요 (평가자마다 열 하나).');
  if (vars.length > AGR_MAX_RATERS_) throw new Error('평가자 변수는 ' + AGR_MAX_RATERS_ + '개까지 선택할 수 있습니다.');
  var opt = request.options || {};
  var conf = agrConf_(opt);
  var weights = opt.weights === 'linear' || opt.weights === 'quadratic' ? opt.weights : 'none';
  var rowsOk = agrComplete_(dataset, vars);
  var N = rowsOk.length, k = vars.length;
  if (N < 2) throw new Error('모든 평가자가 평가한 사례가 ' + N + '건입니다 (2건 이상 필요).');
  var names = vars.map(agrName_);
  var cols = vars.map(function (v) { return rowsOk.map(function (r) { return dataset.data[v.name][r]; }); });
  // 범주: 모든 평가자의 값의 합집합 (값 순)
  var cats = [], seen = {};
  cols.forEach(function (c) {
    c.forEach(function (x) {
      var key = typeof x + ':' + x;
      if (!seen[key]) { seen[key] = true; cats.push({ value: x }); }
    });
  });
  cats.sort(d.compareValues_);
  var pos = {};
  cats.forEach(function (c, i) {
    pos[typeof c.value + ':' + c.value] = i;
    c.label = String(c.value);
    for (var q = 0; q < vars.length; q++) {
      var lab = vars[q].valueLabels[String(c.value)];
      if (lab !== undefined) { c.label = lab; break; }
    }
  });
  var r = cats.length;
  if (r < 2) throw new Error('평가에 쓰인 범주가 ' + r + '개입니다. 2개 이상이어야 kappa 를 계산할 수 있습니다.');
  if (r > AGR_MAX_CATS_) throw new Error('평가 범주가 ' + r + '개입니다 (최대 ' + AGR_MAX_CATS_ + '개). 점수가 연속이면 급내상관(ICC)을 쓰세요.');
  var idx = cols.map(function (c) { return c.map(function (x) { return pos[typeof x + ':' + x]; }); });

  var warnings = [];
  if (dataset.n > N) warnings.push('목록별 결측 제거로 ' + (dataset.n - N) + '건을 제외했습니다 (평가자 중 하나라도 비어 있는 사례). 사용한 사례는 ' + N + '건입니다.');
  if (N < 30) warnings.push('사례가 ' + N + '건으로 적어 kappa 의 신뢰구간이 매우 넓을 수 있습니다.');
  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [{ label: '분석에 포함 (사례)', n: N }, { label: '결측으로 제외', n: dataset.n - N }, { label: '평가자 수 (변수)', n: k }, { label: '평가 범주 수', n: r }, { _kind: 'total', label: '합계', n: dataset.n }],
    notes: ['평가자: ' + names.join(', ') + '. 범주 ' + r + '개: ' + cats.map(function (c) { return c.label; }).join(', ') + '.']
  });

  var F = AGR_F3_;
  if (k === 2) {
    var T = [];
    for (var i = 0; i < r; i++) T.push(new Array(r).fill(0));
    for (var t = 0; t < N; t++) T[idx[0][t]][idx[1][t]]++;
    var cols2 = [{ key: 'name', label: names[0] + ' ＼ ' + names[1] }];
    cats.forEach(function (c, j) { cols2.push({ key: 'c' + j, label: c.label, format: { type: 'int' } }); });
    cols2.push({ key: 'sum', label: '합계', format: { type: 'int' } });
    var xrows = T.map(function (rw, a) {
      var row = { name: cats[a].label, sum: rw.reduce(function (s, v) { return s + v; }, 0) };
      rw.forEach(function (v, b) { row['c' + b] = v; });
      return row;
    });
    var tot = { _kind: 'total', name: '합계', sum: N };
    cats.forEach(function (c, b) { var s = 0; T.forEach(function (rw) { s += rw[b]; }); tot['c' + b] = s; });
    xrows.push(tot);
    blocks.push({ type: 'table', caption: '평가 일치 교차표 (행: ' + names[0] + ', 열: ' + names[1] + ')', columns: cols2, rows: xrows, notes: ['대각선의 값이 두 평가자가 같은 범주로 평가한 사례 수이다.'] });
    var ck = cohenKappa(T, weights, conf);
    var wname = weights === 'linear' ? '선형 가중 kappa' : weights === 'quadratic' ? '이차 가중 kappa' : 'Cohen kappa';
    blocks.push({
      type: 'table', caption: wname,
      columns: [
        { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'po', label: weights === 'none' ? '관찰 일치율' : '관찰 일치도(가중)', format: F },
        { key: 'pe', label: weights === 'none' ? '우연 일치율' : '우연 일치도(가중)', format: F },
        { key: 'kappa', label: 'kappa', format: F },
        { key: 'se', label: '표준오차', format: F },
        { key: 'lo', label: conf + '% 신뢰구간 하한', format: F },
        { key: 'hi', label: conf + '% 신뢰구간 상한', format: F },
        { key: 'z', label: 'z (kappa = 0 검정)', format: F },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
      ].concat(weights === 'none' ? [{ key: 'kmax', label: 'kappa 최대값', format: F }] : []),
      rows: [ck],
      notes: [
        'kappa = (관찰 일치 − 우연 일치) / (1 − 우연 일치). 일치도 수준(Landis & Koch 기준): 0.2 이하 매우 낮음, 0.2~0.4 낮음, 0.4~0.6 보통, 0.6~0.8 높음, 0.8 초과 거의 완벽 → 이 자료는 "' + agrKappaWord_(ck.kappa) + '".',
        '표준오차와 신뢰구간은 Fleiss·Cohen·Everitt(1969)의 식이고, z와 유의확률은 kappa = 0 이라는 귀무가설의 표준오차로 구했다. 유의해도 일치가 "높다"는 뜻은 아니니 kappa의 크기와 신뢰구간을 보세요.' +
          (weights === 'none' ? ' kappa 최대값은 각 평가자의 범주 비율이 같을 때 가능한 kappa의 최댓값이다.' : ' 가중 kappa는 범주의 순서를 고려해 먼 범주끼리의 불일치를 더 크게 벌점을 준다 (범주 순서는 값의 크기 순).')
      ]
    });
    if (ck.kappa === null) warnings.push('우연 일치가 100%라(한 범주로만 평가) kappa 를 계산할 수 없습니다.');
    else {
      if (ck.kappa < 0) warnings.push('kappa 가 음수입니다. 두 평가자의 일치가 우연보다도 낮습니다 (평가 기준이나 변수 선택을 확인하세요).');
      var maxMarg = 0, rs = [], cs = [];
      T.forEach(function (rw, a) { rs.push(rw.reduce(function (s, v) { return s + v; }, 0) / N); });
      for (var b2 = 0; b2 < r; b2++) { var s2 = 0; T.forEach(function (rw) { s2 += rw[b2]; }); cs.push(s2 / N); }
      maxMarg = Math.max(Math.max.apply(null, rs), Math.max.apply(null, cs));
      if (weights === 'none' && ck.po >= 0.8 && ck.kappa < 0.4) warnings.push('관찰 일치율은 ' + (ck.po * 100).toFixed(1) + '%로 높은데 kappa 는 ' + ck.kappa.toFixed(3) + '로 낮습니다. 한쪽 범주에 평가가 몰려(최대 ' + (maxMarg * 100).toFixed(0) + '%) 우연 일치가 높아서 생기는 "kappa 역설"일 수 있습니다. 일치율도 함께 보고하세요.');
    }
  } else {
    var X = [];
    for (var s3 = 0; s3 < N; s3++) {
      var xr = new Array(r).fill(0);
      for (var q2 = 0; q2 < k; q2++) xr[idx[q2][s3]]++;
      X.push(xr);
    }
    var fk = fleissKappa(X);
    var frows = cats.map(function (c, j) {
      return { name: c.label, cnt: fk.counts[j], p: fk.p[j], kappa: fk.kappaJ[j], z: fk.zJ[j], pv: fk.pJ[j] };
    });
    frows.push({ _kind: 'total', name: '전체 (Fleiss kappa)', cnt: N * k, p: 1, kappa: fk.kappa, se: fk.se, z: fk.z, pv: fk.pval });
    blocks.push({
      type: 'table', caption: 'Fleiss kappa (평가자 ' + k + '명)',
      columns: [
        { key: 'name', label: '범주' }, { key: 'cnt', label: '평가 횟수', format: { type: 'int' } }, { key: 'p', label: '비율', format: F },
        { key: 'kappa', label: 'kappa', format: F }, { key: 'se', label: '표준오차 (영가설)', format: F }, { key: 'z', label: 'z', format: F }, { key: 'pv', label: '유의확률', format: { type: 'pvalue' } }
      ],
      rows: frows.map(function (rw, q3) { if (q3 < frows.length - 1) rw.se = fk.seJ; return rw; }),
      notes: [
        '평균 사례 일치도 P̄ = ' + fk.pbar.toFixed(3) + ', 우연 일치 P̄e = ' + fk.pe.toFixed(3) + ', kappa = (P̄ − P̄e)/(1 − P̄e). 일치도 수준(Landis & Koch 기준): 0.2 이하 매우 낮음, 0.2~0.4 낮음, 0.4~0.6 보통, 0.6~0.8 높음, 0.8 초과 거의 완벽 → 이 자료는 "' + agrKappaWord_(fk.kappa) + '".',
        '범주별 kappa 는 그 범주를 다른 범주와 구분하는 데서의 일치도이다. 표준오차와 z는 kappa = 0 귀무가설에서 Fleiss(1971)의 식으로 구했다 (신뢰구간은 제공하지 않는다). 가중 kappa 옵션은 평가자가 2명일 때만 쓴다.'
      ]
    });
    if (fk.kappa === null) warnings.push('우연 일치가 100%라(한 범주로만 평가) kappa 를 계산할 수 없습니다.');
    // 평가자 쌍별 Cohen kappa 행렬 (옵션의 가중치 적용)
    var pairCols = [{ key: 'name', label: '' }];
    names.forEach(function (nm, j) { pairCols.push({ key: 'c' + j, label: nm, format: F }); });
    var pairRows = [], sumK = 0, cntK = 0;
    for (var a2 = 0; a2 < k; a2++) {
      var prow = { name: names[a2] };
      for (var b3 = 0; b3 < k; b3++) {
        if (a2 === b3) { prow['c' + b3] = 1; continue; }
        var Tp = [];
        for (var u = 0; u < r; u++) Tp.push(new Array(r).fill(0));
        for (var t2 = 0; t2 < N; t2++) Tp[idx[a2][t2]][idx[b3][t2]]++;
        var kp = cohenKappa(Tp, weights, conf).kappa;
        prow['c' + b3] = kp;
        if (a2 < b3 && kp !== null) { sumK += kp; cntK++; }
      }
      pairRows.push(prow);
    }
    blocks.push({
      type: 'table', caption: '평가자 쌍별 Cohen kappa' + (weights === 'none' ? '' : ' (' + (weights === 'linear' ? '선형' : '이차') + ' 가중)'),
      columns: pairCols, rows: pairRows,
      notes: ['쌍별 kappa의 평균(Light kappa) = ' + (cntK ? (sumK / cntK).toFixed(3) : '.') + '. 어느 평가자 쌍이 일치가 낮은지 보는 용도이다. 대각선은 자기 자신이라 1이다.']
    });
  }
  return { analysis: 'kappa', title: k === 2 ? '평가자 일치도 (Cohen kappa)' : '평가자 일치도 (Fleiss kappa)', blocks: blocks, warnings: warnings };
}

function runIcc(dataset, request) {
  var d = agrDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < 2) throw new Error('급내상관(ICC)은 평가자(또는 측정) 변수를 2개 이상 선택하세요 (평가자마다 열 하나).');
  if (vars.length > AGR_MAX_RATERS_) throw new Error('변수는 ' + AGR_MAX_RATERS_ + '개까지 선택할 수 있습니다.');
  var opt = request.options || {};
  var conf = agrConf_(opt);
  var rowsOk = agrComplete_(dataset, vars);
  var N = rowsOk.length, k = vars.length;
  if (N < 3) throw new Error('모든 변수에 값이 있는 사례가 ' + N + '건입니다 (3건 이상 필요).');
  var names = vars.map(agrName_);
  var cols = vars.map(function (v) { return rowsOk.map(function (r) { return dataset.data[v.name][r]; }); });
  var res = iccCompute(cols, conf);
  var warnings = [];
  if (dataset.n > N) warnings.push('목록별 결측 제거로 ' + (dataset.n - N) + '건을 제외했습니다 (변수 중 하나라도 비어 있는 대상). 사용한 대상은 ' + N + '건입니다.');
  if (N < 20) warnings.push('대상이 ' + N + '건으로 적어 ICC의 신뢰구간이 매우 넓을 수 있습니다 (대상 30건 이상을 권장하는 기준이 흔합니다).');
  if (res.rows[0].icc === null) warnings.push('값에 변동이 없거나 오차가 0이라 ICC 를 계산할 수 없습니다.');
  var F = AGR_F3_;
  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [{ label: '분석에 포함 (대상)', n: N }, { label: '결측으로 제외', n: dataset.n - N }, { label: '평가자·측정 수 (변수)', n: k }, { _kind: 'total', label: '합계', n: dataset.n }],
    notes: ['평가자(또는 측정): ' + names.join(', ') + '. 한 대상이 모든 변수에 값이 있어야 하며 비어 있는 대상은 통째로 제외했다.']
  });
  blocks.push({
    type: 'table', caption: '기술통계',
    columns: [{ key: 'name', label: '평가자(변수)' }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'mean', label: '평균', format: F }, { key: 'sd', label: '표준편차', format: F }],
    rows: names.map(function (nm, j) { return { name: nm, n: N, mean: res.means[j], sd: res.sds[j] }; })
  });
  if (opt.anova !== false) blocks.push({
    type: 'table', caption: '이원 분산분석표 (대상 × 평가자)',
    columns: [{ key: 'source', label: '' }, { key: 'ss', label: '제곱합', format: F }, { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'ms', label: '평균제곱', format: F }],
    rows: [
      { source: '대상(사례) 사이', ss: res.ssb, df: res.dfb, ms: res.msb },
      { source: '평가자 사이', ss: res.ssj, df: res.dfj, ms: res.msj },
      { source: '잔차', ss: res.sse, df: res.dfe, ms: res.mse },
      { _kind: 'total', source: '합계', ss: res.sst, df: N * k - 1 }
    ],
    notes: ['ICC 는 이 평균제곱들로 계산한다. 대상 사이 변동이 평가자·잔차 변동보다 클수록 ICC 가 1에 가깝다. 대상 안 평균제곱(일원 모형) = ' + res.msw.toFixed(3) + '.']
  });
  blocks.push({
    type: 'table', caption: '급내상관계수 (ICC)',
    columns: [
      { key: 'label', label: '유형' }, { key: 'desc', label: '모형 · 일치 기준 · 측정' }, { key: 'icc', label: 'ICC', format: F },
      { key: 'lo', label: conf + '% 신뢰구간 하한', format: F }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: F },
      { key: 'F', label: 'F', format: F }, { key: 'df1', label: '자유도1', format: { type: 'int' } }, { key: 'df2', label: '자유도2', format: { type: 'int' } },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'word', label: '해석(추정값)' }
    ],
    rows: AGR_ICC_TYPES_.map(function (tp, q) {
      var rw = res.rows[q];
      var word = rw.icc === null ? '계산 불가' : rw.icc < 0.5 ? '낮음' : rw.icc < 0.75 ? '보통' : rw.icc < 0.9 ? '좋음' : '매우 좋음';
      return { _bold: tp.key === 'icc2' || tp.key === 'icc3', label: tp.label, desc: tp.desc, icc: rw.icc, lo: rw.lo, hi: rw.hi, F: rw.F, df1: rw.df1, df2: rw.df2, p: rw.p, word: word };
    }),
    notes: [
      '어느 값을 볼지: 평가자가 같은 사람들에 고정돼 있고 점수의 순서(일관성)만 중요하면 ICC(3,1), 평가자가 더 큰 집단에서 뽑힌 표본이고 점수 자체의 일치(절대 일치도)가 중요하면 ICC(2,1)을 봅니다 (굵게 표시한 단일 측정 값). 실제 측정에서 여러 평가자의 평균을 쓸 때만 평균 측정(k) 값을 봅니다. 대상마다 서로 다른 평가자가 평가했다면 ICC(1,1)입니다.',
      '해석 기준(Koo & Li 2016): 0.5 미만 낮음, 0.5~0.75 보통, 0.75~0.9 좋음, 0.9 초과 매우 좋음. 이 기준은 추정값이 아니라 신뢰구간으로 판단하는 것이 권장됩니다 (하한이 어느 구간인지 보세요).',
      '유의확률은 ICC = 0 검정이라 ICC 가 낮아도 유의할 수 있으니 크기와 신뢰구간을 보세요. 모형은 일원 확률, 이원 확률, 이원 혼합이 있고, 유형은 절대 일치도와 일관성입니다 (ICC(1), ICC(2), ICC(3)의 단일·평균 측정).'
    ]
  });
  return { analysis: 'icc', title: '급내상관계수 (ICC)', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = { cohenKappa: cohenKappa, fleissKappa: fleissKappa, iccCompute: iccCompute, runKappa: runKappa, runIcc: runIcc };
}
