/**
 * 이원 분산분석(집단 요인 2개, 개체-간 설계). 순수 함수: (dataset, request) → 결과 객체.
 * 결과표 구성: 개체-간 효과 검정, 추정 주변평균, 사후검정, 단순 주효과 (추정 주변평균의 쌍별 비교와 같은 검정).
 *
 * 제곱합: 유형 III(기본) 또는 유형 II. 효과(합 0) 코딩한 설계의 축소 모형과 전체 모형의 잔차제곱합 차이로 구한다.
 * 설계가 범주형이라 모든 모형의 잔차제곱합 = 셀 안 제곱합 + Σ n_ij (셀 평균 − 적합값)² 이다. 그래서 사례 단위가 아니라
 * 셀 평균(a×b개)에 대한 가중(n_ij) 최소제곱으로 구한다 (사례 수가 아무리 커도 계산량이 셀 수에만 달려 있다).
 *   유형 III: 항 하나를 뺀 모형(나머지 항은 모두 포함) − 전체 모형
 *   유형 II : A 는 B 만 있는 모형 − (A + B), B 는 A 만 있는 모형 − (A + B), A×B 는 (A + B) − 전체
 * 사례 수가 셀마다 같으면 두 유형은 일치한다. 모든 셀(A×B 조합)에 사례가 있어야 한다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function an2Deps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return {
      getVars: data.getVars, compareValues_: data.compareValues_,
      pf: dist.pf, pt: dist.pt, qt: dist.qt, ptukey: dist.ptukey, qtukey: dist.qtukey,
      leveneTest: require('./Stats_Assumptions.js').leveneTest, qrLeastSquares_: require('./Stats_Regression.js').qrLeastSquares_,
      lineChartBlock: require('./Charts.js').lineChartBlock
    };
  }
  return {
    getVars: getVars, compareValues_: compareValues_, pf: pf, pt: pt, qt: qt, ptukey: ptukey, qtukey: qtukey,
    leveneTest: leveneTest, qrLeastSquares_: qrLeastSquares_, lineChartBlock: lineChartBlock
  };
}

var AN2_MAX_LEVELS_ = 10;
var AN2_SMALL_CELL_ = 5;
function an2Name_(v) { return v.label ? v.label : v.name; }

/* ---------- 계산 ---------- */

/** 범주 변수의 수준 목록(값 순, 레이블)과 사례별 수준 번호 */
function an2Levels_(v, vals) {
  var d = an2Deps_();
  var cats = [], seen = {};
  vals.forEach(function (x) {
    var key = typeof x + ':' + x;
    if (!seen[key]) { seen[key] = true; cats.push({ value: x }); }
  });
  cats.sort(d.compareValues_);
  var pos = {};
  cats.forEach(function (c, i) {
    pos[typeof c.value + ':' + c.value] = i;
    var label = v.valueLabels[String(c.value)];
    c.label = label !== undefined ? label : String(c.value);
  });
  return { levels: cats, index: vals.map(function (x) { return pos[typeof x + ':' + x]; }) };
}

/** 효과(합 0) 코딩: 수준 i 의 a−1 개 값 (마지막 수준은 모두 −1) */
function an2Effect_(i, k) {
  var r = [];
  for (var l = 0; l < k - 1; l++) r.push(i === l ? 1 : (i === k - 1 ? -1 : 0));
  return r;
}

/**
 * 셀 단위 가중 최소제곱의 잔차제곱합. useA/useB/useAB: 모형에 넣을 항 (절편은 항상). cells: 셀 평균 배열 [i][j] {n, mean}, ssWithin: 셀 안 제곱합
 */
function an2CellSse_(cells, a, b, ssWithin, useA, useB, useAB) {
  var d = an2Deps_();
  var list = [], i, j, k;
  for (i = 0; i < a; i++) for (j = 0; j < b; j++) list.push({ i: i, j: j, n: cells[i][j].n, mean: cells[i][j].mean });
  var X = list.map(function (c) {
    var x = [1], ea = an2Effect_(c.i, a), eb = an2Effect_(c.j, b);
    if (useA) x = x.concat(ea);
    if (useB) x = x.concat(eb);
    if (useAB) ea.forEach(function (va) { eb.forEach(function (vb) { x.push(va * vb); }); });
    return x;
  });
  var p = X[0].length;
  var cols = [];
  for (k = 0; k < p; k++) cols.push(list.map(function (c, r) { return X[r][k] * Math.sqrt(c.n); }));
  var fit = d.qrLeastSquares_(cols, list.map(function (c) { return c.mean * Math.sqrt(c.n); }));
  if (!fit.ok) throw new Error('설계 행렬을 계산할 수 없습니다 (셀 구성을 확인하세요).');
  var sse = ssWithin;
  list.forEach(function (c, r) {
    var f = 0;
    for (k = 0; k < p; k++) f += fit.beta[k] * X[r][k];
    sse += c.n * (c.mean - f) * (c.mean - f);
  });
  return sse;
}

/**
 * y: 숫자 배열, ia/ib: 사례별 요인 수준 번호(0부터), a/b: 수준 수, ssType: 'III' | 'II'.
 * 반환: { n, a, b, cells[i][j] {n, mean, sd}, rowN, rowMean, colN, colMean, terms: [{key, ss, df, ms, F, p, peta}], intercept, model, sse, dfe, mse, sst, r2, adjR2, emptyCell }
 */
function anova2Analyze(y, ia, ib, a, b, ssType) {
  var d = an2Deps_();
  var N = y.length, i, j, t;
  var cn = [], cs = [], cq = [];
  for (i = 0; i < a; i++) { cn.push(new Array(b).fill(0)); cs.push(new Array(b).fill(0)); cq.push(new Array(b).fill(0)); }
  for (t = 0; t < N; t++) { cn[ia[t]][ib[t]]++; cs[ia[t]][ib[t]] += y[t]; }
  for (i = 0; i < a; i++) for (j = 0; j < b; j++) if (cn[i][j] === 0) return { emptyCell: [i, j] };
  var cells = [];
  for (i = 0; i < a; i++) {
    cells.push([]);
    for (j = 0; j < b; j++) cells[i].push({ n: cn[i][j], mean: cs[i][j] / cn[i][j], sd: null });
  }
  for (t = 0; t < N; t++) {
    var cell = cells[ia[t]][ib[t]];
    cq[ia[t]][ib[t]] += (y[t] - cell.mean) * (y[t] - cell.mean);
  }
  for (i = 0; i < a; i++) for (j = 0; j < b; j++) if (cells[i][j].n > 1) cells[i][j].sd = Math.sqrt(cq[i][j] / (cells[i][j].n - 1));

  var res = { n: N, a: a, b: b, cells: cells };
  var dfe = N - a * b;
  if (dfe < 1) return Object.assign(res, { tooFew: true });
  res.dfe = dfe;

  var rowN = new Array(a).fill(0), rowS = new Array(a).fill(0), colN = new Array(b).fill(0), colS = new Array(b).fill(0);
  for (i = 0; i < a; i++) for (j = 0; j < b; j++) { rowN[i] += cn[i][j]; rowS[i] += cs[i][j]; colN[j] += cn[i][j]; colS[j] += cs[i][j]; }
  res.rowN = rowN; res.colN = colN;
  res.rowMean = rowS.map(function (s, k) { return s / rowN[k]; });
  res.colMean = colS.map(function (s, k) { return s / colN[k]; });

  var mean = 0, sumSq = 0;
  y.forEach(function (v) { mean += v; sumSq += v * v; });
  mean /= N;
  var sst = 0;
  y.forEach(function (v) { sst += (v - mean) * (v - mean); });
  res.sst = sst; res.sumSq = sumSq;

  var ssWithin = 0;
  for (i = 0; i < a; i++) for (j = 0; j < b; j++) ssWithin += cq[i][j];
  var sseFull = ssWithin;                                    // 셀 평균을 모두 추정하는 전체(포화) 모형의 잔차제곱합
  res.sse = sseFull;
  res.mse = sseFull / dfe;
  var ss = {};
  if (ssType === 'II') {
    var sseAB = an2CellSse_(cells, a, b, ssWithin, true, true, false);
    ss.A = an2CellSse_(cells, a, b, ssWithin, false, true, false) - sseAB;
    ss.B = an2CellSse_(cells, a, b, ssWithin, true, false, false) - sseAB;
    ss.AB = sseAB - sseFull;
  } else {
    ss.A = an2CellSse_(cells, a, b, ssWithin, false, true, true) - sseFull;
    ss.B = an2CellSse_(cells, a, b, ssWithin, true, false, true) - sseFull;
    ss.AB = an2CellSse_(cells, a, b, ssWithin, true, true, false) - sseFull;
  }
  var dfs = { A: a - 1, B: b - 1, AB: (a - 1) * (b - 1) };
  res.terms = ['A', 'B', 'AB'].map(function (key) {
    var s = Math.max(0, ss[key]);
    var ms = s / dfs[key];
    var F = res.mse > 0 ? ms / res.mse : null;
    return { key: key, ss: s, df: dfs[key], ms: ms, F: F, p: F === null ? null : d.pf(F, dfs[key], dfe, false), peta: s + sseFull > 0 ? s / (s + sseFull) : null };
  });
  var msModel = (sst - sseFull) / (a * b - 1), Fm = res.mse > 0 ? msModel / res.mse : null;
  res.model = { ss: sst - sseFull, df: a * b - 1, ms: msModel, F: Fm, p: Fm === null ? null : d.pf(Fm, a * b - 1, dfe, false), peta: sst > 0 ? (sst - sseFull) / sst : null };
  if (ssType !== 'II') {
    var mu = 0, inv = 0;
    for (i = 0; i < a; i++) for (j = 0; j < b; j++) { mu += cells[i][j].mean; inv += 1 / cells[i][j].n; }
    mu /= a * b;
    var ssI = mu * mu * (a * b) * (a * b) / inv, Fi = res.mse > 0 ? ssI / res.mse : null;
    res.intercept = { ss: ssI, df: 1, ms: ssI, F: Fi, p: Fi === null ? null : d.pf(Fi, 1, dfe, false) };
  }
  res.r2 = sst > 0 ? 1 - sseFull / sst : null;
  res.adjR2 = sst > 0 ? 1 - (1 - res.r2) * (N - 1) / dfe : null;
  return res;
}

/** 추정 주변평균: 요인 A·B 수준별(셀 평균의 비가중 평균)과 셀별. SE = √(MSE · Σ 1/n / 개수²), 자유도 = 오차 자유도 */
function anova2Emm(res, conf) {
  var d = an2Deps_();
  var a = res.a, b = res.b, i, j;
  var tc = d.qt((1 - conf / 100) / 2, res.dfe, false);
  var mk = function (est, invSum, count) {
    var se = res.mse > 0 ? Math.sqrt(res.mse * invSum) / count : null;
    return { mean: est, se: se, lo: se === null ? null : est - tc * se, hi: se === null ? null : est + tc * se };
  };
  var rowE = [], colE = [], cellE = [];
  for (i = 0; i < a; i++) {
    var s = 0, inv = 0;
    for (j = 0; j < b; j++) { s += res.cells[i][j].mean; inv += 1 / res.cells[i][j].n; }
    rowE.push(mk(s / b, inv, b));
  }
  for (j = 0; j < b; j++) {
    var s2 = 0, inv2 = 0;
    for (i = 0; i < a; i++) { s2 += res.cells[i][j].mean; inv2 += 1 / res.cells[i][j].n; }
    colE.push(mk(s2 / a, inv2, a));
  }
  for (i = 0; i < a; i++) {
    cellE.push([]);
    for (j = 0; j < b; j++) cellE[i].push(mk(res.cells[i][j].mean, 1 / res.cells[i][j].n, 1));
  }
  return { rows: rowE, cols: colE, cells: cellE };
}

/**
 * 주효과의 사후검정(관측 주변평균 비교, 오차평균제곱·오차자유도는 이원 모형의 것). 평균차이 = 수준 i − 수준 j.
 * method: 'tukey' (표본 크기가 다르면 Tukey-Kramer) | 'bonferroni'. 반환: [{ i, j, diff, se, p, lo, hi }]
 */
function anova2PostHoc(method, means, ns, mse, dfe, conf) {
  var d = an2Deps_();
  var k = means.length, alpha = 1 - conf / 100, pairs = k * (k - 1) / 2, rows = [];
  var qc = method === 'tukey' ? d.qtukey(alpha, k, dfe, false) : null;
  var tc = method === 'bonferroni' ? d.qt(alpha / (2 * pairs), dfe, false) : null;
  for (var i = 0; i < k - 1; i++) {
    for (var j = i + 1; j < k; j++) {
      var row = { i: i, j: j, diff: means[i] - means[j], se: null, p: null, lo: null, hi: null };
      rows.push(row);
      if (!(mse > 0)) continue;
      row.se = Math.sqrt(mse * (1 / ns[i] + 1 / ns[j]));
      var t = row.diff / row.se, hw;
      if (method === 'tukey') {
        row.p = d.ptukey(Math.abs(t) * Math.SQRT2, k, dfe, false);
        hw = qc * row.se / Math.SQRT2;
      } else {
        row.p = Math.min(1, pairs * d.pt(-Math.abs(t), dfe) * 2);
        hw = tc * row.se;
      }
      row.lo = row.diff - hw;
      row.hi = row.diff + hw;
    }
  }
  return rows;
}

/**
 * 단순 주효과: 한 요인의 각 수준에서 다른 요인의 효과. SS = Σ n_ij (셀 평균 − 그 수준의 가중평균)², 자유도 = 수준 수 − 1,
 * F = (SS/df) / MSE(이원 모형), 자유도2 = 오차 자유도 (추정 주변평균 비교의 일변량 검정).
 * 반환: { aAtB: [B 수준 j 마다 A 효과], bAtA: [A 수준 i 마다 B 효과] }  각 항목 { ss, df, ms, F, p }
 */
function anova2Simple(res) {
  var d = an2Deps_();
  var a = res.a, b = res.b, i, j;
  var one = function (means, ns) {
    var tot = 0, sw = 0;
    means.forEach(function (m, k) { tot += ns[k]; sw += ns[k] * m; });
    var g = sw / tot, ss = 0;
    means.forEach(function (m, k) { ss += ns[k] * (m - g) * (m - g); });
    var df = means.length - 1, ms = ss / df, F = res.mse > 0 ? ms / res.mse : null;
    return { ss: ss, df: df, ms: ms, F: F, p: F === null ? null : d.pf(F, df, res.dfe, false) };
  };
  var aAtB = [], bAtA = [];
  for (j = 0; j < b; j++) {
    var m1 = [], n1 = [];
    for (i = 0; i < a; i++) { m1.push(res.cells[i][j].mean); n1.push(res.cells[i][j].n); }
    aAtB.push(one(m1, n1));
  }
  for (i = 0; i < a; i++) {
    var m2 = [], n2 = [];
    for (j = 0; j < b; j++) { m2.push(res.cells[i][j].mean); n2.push(res.cells[i][j].n); }
    bAtA.push(one(m2, n2));
  }
  return { aAtB: aAtB, bAtA: bAtA };
}

/* ---------- 실행 ---------- */

function runAnova2(dataset, request) {
  var d = an2Deps_();
  var dep = d.getVars(dataset, request.dependent, { numeric: true });
  if (dep.length !== 1) throw new Error('종속변수는 하나만 선택하세요.');
  var f1 = d.getVars(dataset, request.factor1), f2 = d.getVars(dataset, request.factor2);
  if (f1.length !== 1 || f2.length !== 1) throw new Error('요인 1과 요인 2는 각각 변수 하나씩 선택하세요.');
  var seen = {};
  seen[dep[0].name] = '종속변수';
  [f1[0], f2[0]].forEach(function (v, k) {
    if (seen[v.name]) throw new Error('변수 "' + v.name + '"를 ' + seen[v.name] + '와(과) 요인 ' + (k + 1) + '에 동시에 쓸 수 없습니다.');
    seen[v.name] = '요인 ' + (k + 1);
  });
  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  var ssType = opt.ss === 'II' ? 'II' : 'III';
  var posthoc = opt.posthoc === 'bonferroni' ? 'bonferroni' : (opt.posthoc === 'none' || opt.posthoc === false ? null : 'tukey');
  var nameA = an2Name_(f1[0]), nameB = an2Name_(f2[0]);

  // 목록별 결측 제거
  var rowsOk = [];
  for (var r = 0; r < dataset.n; r++) {
    if (dataset.data[dep[0].name][r] !== null && dataset.data[f1[0].name][r] !== null && dataset.data[f2[0].name][r] !== null) rowsOk.push(r);
  }
  var N = rowsOk.length;
  var pick = function (name) { return rowsOk.map(function (x) { return dataset.data[name][x]; }); };
  var y = pick(dep[0].name);
  var LA = an2Levels_(f1[0], pick(f1[0].name)), LB = an2Levels_(f2[0], pick(f2[0].name));
  var a = LA.levels.length, b = LB.levels.length;
  [[a, nameA], [b, nameB]].forEach(function (x) {
    if (x[0] < 2) throw new Error('요인 "' + x[1] + '"의 수준이 ' + x[0] + '개입니다. 2개 이상이어야 합니다.');
    if (x[0] > AN2_MAX_LEVELS_) throw new Error('요인 "' + x[1] + '"의 수준이 ' + x[0] + '개입니다 (최대 ' + AN2_MAX_LEVELS_ + '개). 범주가 적은 변수를 선택하세요.');
  });

  var res = anova2Analyze(y, LA.index, LB.index, a, b, ssType);
  if (res.emptyCell) {
    throw new Error('요인 조합 "' + nameA + ' = ' + LA.levels[res.emptyCell[0]].label + ', ' + nameB + ' = ' + LB.levels[res.emptyCell[1]].label +
      '"에 사례가 없습니다. 이원 분산분석은 모든 조합에 사례가 있어야 합니다. 범주를 합치거나 다른 변수를 선택하세요.');
  }
  if (res.tooFew) throw new Error('사례 수(' + N + ')가 셀 수(' + (a * b) + ')보다 많아야 합니다. 셀마다 평균 1건보다 많아야 오차를 추정할 수 있습니다.');

  var warnings = [];
  if (dataset.n > N) warnings.push('목록별 결측 제거로 ' + (dataset.n - N) + '건을 제외했습니다 (종속변수와 두 요인 중 하나라도 빠진 사례). 사용한 사례는 ' + N + '건입니다.');
  var balanced = true, small = 0, firstN = res.cells[0][0].n;
  res.cells.forEach(function (row) { row.forEach(function (c) { if (c.n !== firstN) balanced = false; if (c.n < AN2_SMALL_CELL_) small++; }); });
  if (!balanced) warnings.push('셀마다 사례 수가 다릅니다 (불균형 설계). 제곱합 유형에 따라 주효과 결과가 달라질 수 있습니다 (지금: 유형 ' + ssType + '). 보통 유형 III을 씁니다.');
  if (small > 0) warnings.push('사례가 ' + AN2_SMALL_CELL_ + '건 미만인 셀이 ' + small + '개 있습니다. 결과가 불안정할 수 있습니다.');
  if (!(res.mse > 0)) warnings.push('모든 셀에서 값이 같아 오차가 0입니다. F와 유의확률을 계산할 수 없습니다.');

  var F3 = { type: 'fixed', digits: 3 };
  var depName = an2Name_(dep[0]);
  var blocks = [];

  // 기술통계
  var descRows = [];
  var labA = function (i) { return LA.levels[i].label; }, labB = function (j) { return LB.levels[j].label; };
  var cellsOf = function (i, j) { return res.cells[i][j]; };
  var allByCell = [];
  var sdOf = function (vals) {
    if (vals.length < 2) return null;
    var m = 0, s = 0;
    vals.forEach(function (v) { m += v; });
    m /= vals.length;
    vals.forEach(function (v) { s += (v - m) * (v - m); });
    return Math.sqrt(s / (vals.length - 1));
  };
  var valsOf = function (fa, fb) {
    var out = [];
    for (var t = 0; t < N; t++) if ((fa === null || LA.index[t] === fa) && (fb === null || LB.index[t] === fb)) out.push(y[t]);
    return out;
  };
  var meanOf = function (vals) { var s = 0; vals.forEach(function (v) { s += v; }); return s / vals.length; };
  var pushDesc = function (la, lb, fa, fb, bold) {
    var vals = valsOf(fa, fb);
    var row = { a: la, b: lb, mean: meanOf(vals), sd: sdOf(vals), n: vals.length };
    if (bold) row._kind = 'total';
    descRows.push(row);
  };
  var i, j;
  for (i = 0; i < a; i++) {
    for (j = 0; j < b; j++) pushDesc(labA(i), labB(j), i, j, false);
    pushDesc(labA(i), '전체', i, null, false);
  }
  for (j = 0; j < b; j++) pushDesc('전체', labB(j), null, j, false);
  pushDesc('전체', '전체', null, null, true);
  blocks.push({
    type: 'table',
    caption: '기술통계 (종속변수: ' + depName + ')',
    columns: [
      { key: 'a', label: nameA },
      { key: 'b', label: nameB },
      { key: 'mean', label: '평균', format: F3 },
      { key: 'sd', label: '표준편차', format: F3 },
      { key: 'n', label: 'N', format: { type: 'int' } }
    ],
    rows: descRows,
    notes: ['관측한 평균이다 (가중). 추정 주변평균은 옵션으로 본다. 요인: ' + nameA + ' ' + a + '수준, ' + nameB + ' ' + b + '수준, 사례 ' + N + '건.']
  });

  // Levene (평균 기준)
  if (opt.levene !== false) {
    for (i = 0; i < a; i++) for (j = 0; j < b; j++) allByCell.push(valsOf(i, j));
    var lv = d.leveneTest(allByCell, false);
    if (lv.p !== null && lv.p < 0.05) warnings.push('Levene 등분산 검정이 유의합니다 (p = ' + (lv.p < 0.001 ? '<.001' : lv.p.toFixed(3).replace(/^0/, '')) + '). 셀 사이 분산이 같다는 가정이 의심됩니다. 셀 크기가 비슷하면 영향이 작지만, 다르면 결과를 조심해서 해석하세요.');
    blocks.push({
      type: 'table',
      caption: 'Levene 등분산 검정 (종속변수: ' + depName + ')',
      columns: [
        { key: 'W', label: 'Levene 통계량', format: F3 },
        { key: 'df1', label: '자유도1', format: { type: 'int' } },
        { key: 'df2', label: '자유도2', format: { type: 'int' } },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
      ],
      rows: [{ W: lv.W, df1: lv.df1, df2: lv.df2, p: lv.p }],
      notes: ['각 셀(' + nameA + ' × ' + nameB + ' 조합, ' + (a * b) + '개)의 오차분산이 같다는 귀무가설을 평균 기준으로 검정한다. 유의하면(p < .05) 등분산 가정이 의심된다.']
    });
  }

  // 개체-간 효과 검정
  var tRows = [{ source: '수정된 모형', ss: res.model.ss, df: res.model.df, ms: res.model.ms, F: res.model.F, p: res.model.p, peta: res.model.peta }];
  if (res.intercept) tRows.push({ source: '절편', ss: res.intercept.ss, df: res.intercept.df, ms: res.intercept.ms, F: res.intercept.F, p: res.intercept.p });
  var names3 = { A: nameA, B: nameB, AB: nameA + ' × ' + nameB };
  res.terms.forEach(function (tm) { tRows.push({ source: names3[tm.key], ss: tm.ss, df: tm.df, ms: tm.ms, F: tm.F, p: tm.p, peta: tm.peta }); });
  tRows.push({ source: '오차', ss: res.sse, df: res.dfe, ms: res.mse });
  tRows.push({ _kind: 'total', source: '합계', ss: res.sumSq, df: N });
  tRows.push({ _kind: 'total', source: '수정 합계', ss: res.sst, df: N - 1 });
  var tNotes = ['제곱합 유형 ' + ssType + (ssType === 'III' ? ' (기본).' : ': 주효과는 다른 주효과를 보정하고 상호작용은 포함하지 않는다.') + ' 부분 η² = 제곱합 /(제곱합 + 오차 제곱합).'];
  if (res.r2 !== null) tNotes.push('R 제곱 = ' + res.r2.toFixed(3).replace(/^0/, '') + ' (수정된 R 제곱 = ' + res.adjR2.toFixed(3).replace(/^(-?)0/, '$1') + ').');
  var inter = res.terms[2];
  if (inter.p !== null && inter.p < 0.05) tNotes.push('상호작용이 유의합니다 (p = ' + (inter.p < 0.001 ? '<.001' : inter.p.toFixed(3).replace(/^0/, '')) + '). 이때는 주효과를 그대로 해석하지 말고 상호작용 그림과 단순 주효과를 보세요.');
  else if (inter.p !== null) tNotes.push('상호작용이 유의하지 않습니다. 주효과를 해석하세요.');
  blocks.push({
    type: 'table',
    caption: '개체-간 효과 검정 (종속변수: ' + depName + ')',
    columns: [
      { key: 'source', label: '' },
      { key: 'ss', label: '제곱합', format: F3 },
      { key: 'df', label: '자유도', format: { type: 'int' } },
      { key: 'ms', label: '평균제곱', format: F3 },
      { key: 'F', label: 'F', format: F3 },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
      { key: 'peta', label: '부분 η²', format: F3 }
    ],
    rows: tRows,
    notes: tNotes
  });

  // 추정 주변평균
  if (opt.emmeans) {
    var em = anova2Emm(res, conf);
    var eRows = [];
    em.rows.forEach(function (e, k) { eRows.push({ f: nameA, lev: labA(k), mean: e.mean, se: e.se, lo: e.lo, hi: e.hi }); });
    em.cols.forEach(function (e, k) { eRows.push({ f: nameB, lev: labB(k), mean: e.mean, se: e.se, lo: e.lo, hi: e.hi }); });
    for (i = 0; i < a; i++) for (j = 0; j < b; j++) eRows.push({ f: nameA + ' × ' + nameB, lev: labA(i) + ', ' + labB(j), mean: em.cells[i][j].mean, se: em.cells[i][j].se, lo: em.cells[i][j].lo, hi: em.cells[i][j].hi });
    blocks.push({
      type: 'table',
      caption: '추정 주변평균 (종속변수: ' + depName + ')',
      columns: [
        { key: 'f', label: '요인' },
        { key: 'lev', label: '수준' },
        { key: 'mean', label: '평균', format: F3 },
        { key: 'se', label: '표준오차', format: F3 },
        { key: 'lo', label: conf + '% 신뢰구간 하한', format: F3 },
        { key: 'hi', label: conf + '% 신뢰구간 상한', format: F3 }
      ],
      rows: eRows,
      notes: ['요인별 값은 다른 요인 수준의 셀 평균을 같은 비중으로 평균한 값이다 (사례 수가 다르면 기술통계의 관측 평균과 다를 수 있다). 표준오차와 신뢰구간은 이원 모형의 오차평균제곱과 오차 자유도(' + res.dfe + ')를 쓴다.']
    });
  }

  // 주효과 사후검정
  var phName = posthoc === 'tukey' ? 'Tukey HSD' : 'Bonferroni';
  var addPost = function (name, labelOf, means, ns, k) {
    var rows = anova2PostHoc(posthoc, means, ns, res.mse, res.dfe, conf);
    blocks.push({
      type: 'table',
      caption: '사후검정: ' + name + ' (' + phName + ')',
      columns: [
        { key: 'i', label: '(I) ' + name },
        { key: 'j', label: '(J) ' + name },
        { key: 'diff', label: '평균차이 (I−J)', format: F3 },
        { key: 'se', label: '표준오차', format: F3 },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
        { key: 'lo', label: conf + '% 신뢰구간 하한', format: F3 },
        { key: 'hi', label: conf + '% 신뢰구간 상한', format: F3 }
      ],
      rows: rows.map(function (rw) { return { i: labelOf(rw.i), j: labelOf(rw.j), diff: rw.diff, se: rw.se, p: rw.p, lo: rw.lo, hi: rw.hi }; }),
      notes: ['관측 평균(다른 요인을 무시한 ' + name + ' 수준별 평균)을 비교하고 오차평균제곱·자유도는 이원 모형의 것을 쓴다.' +
        (posthoc === 'tukey' ? ' 사례 수가 다르면 Tukey-Kramer.' : ' 비교 횟수 ' + (k * (k - 1) / 2) + '번만큼 유의확률을 보정한다.')]
    });
  };
  if (posthoc) {
    if (a >= 3) addPost(nameA, labA, res.rowMean, res.rowN, a);
    if (b >= 3) addPost(nameB, labB, res.colMean, res.colN, b);
  }

  // 단순 주효과
  if (opt.simple !== false) {
    var se = anova2Simple(res);
    var sRows = [];
    se.aAtB.forEach(function (e, k) { sRows.push({ eff: nameA + '의 효과 (' + nameB + ' = ' + labB(k) + ')', ss: e.ss, df: e.df, ms: e.ms, F: e.F, p: e.p }); });
    se.bAtA.forEach(function (e, k) { sRows.push({ eff: nameB + '의 효과 (' + nameA + ' = ' + labA(k) + ')', ss: e.ss, df: e.df, ms: e.ms, F: e.F, p: e.p }); });
    blocks.push({
      type: 'table',
      caption: '단순 주효과 (종속변수: ' + depName + ')',
      columns: [
        { key: 'eff', label: '' },
        { key: 'ss', label: '제곱합', format: F3 },
        { key: 'df', label: '자유도', format: { type: 'int' } },
        { key: 'ms', label: '평균제곱', format: F3 },
        { key: 'F', label: 'F', format: F3 },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
      ],
      rows: sRows,
      notes: ['한 요인의 각 수준에서 다른 요인의 효과를 검정한다 (F = 그 수준의 평균제곱 / 이원 모형의 오차평균제곱, 자유도2 = ' + res.dfe + '). 상호작용이 유의할 때 해석하는 표이고, 검정을 여러 개 하므로 유의수준을 나눠서(예: .05/수준 수) 보는 것이 안전하다.' +
        (inter.p !== null && inter.p >= 0.05 ? ' 지금은 상호작용이 유의하지 않으니 참고만 하세요.' : '')]
    });
  }

  // 상호작용 그림: x = 요인 1 수준, 선 = 요인 2 수준 (관측 셀 평균)
  if (opt.plot !== false) {
    var series = [];
    for (j = 0; j < b; j++) {
      var pts = [];
      for (i = 0; i < a; i++) pts.push([i + 1, res.cells[i][j].mean]);
      series.push({ label: nameB + ' = ' + labB(j), points: pts });
    }
    blocks.push(d.lineChartBlock('상호작용 그림: ' + nameA + '별 ' + nameB + ' 평균', nameA, '평균 (' + depName + ')', series, LA.levels.map(function (l) { return l.label; })));
  }

  return { analysis: 'anova2', title: '이원 분산분석', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = {
    anova2Analyze: anova2Analyze,
    anova2Emm: anova2Emm,
    anova2PostHoc: anova2PostHoc,
    anova2Simple: anova2Simple,
    runAnova2: runAnova2
  };
}
