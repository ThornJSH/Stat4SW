/**
 * 이분형 로지스틱 회귀. 순수 함수: (dataset, request) → 결과 객체. 결과표 구성은 일반적인 보고 형식(처리 요약, 종속변수 코딩,
 * 모형 계수의 omnibus 검정, 모형 요약, Hosmer-Lemeshow, 분류표, 방정식의 변수)를 따르고 ROC/AUC를 더했다.
 *
 * 추정: 최대우도, IRLS(Newton-Raphson)를 가중 최소제곱 Householder QR로 푼다 (qrLeastSquares_ 재사용).
 * 독립변수를 평균으로 중심화해 계산한 뒤 원래 척도로 되돌린다. 단계마다 가능도가 줄면 보폭을 절반으로 줄인다.
 * 완전·준완전 분리(추정이 무한대로 발산)와 미수렴은 경고로 알린다.
 * 종속변수는 값이 정확히 2가지여야 하고, 기본은 큰 값이 사건(1)이다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function logDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var reg = require('./Stats_Regression.js');
    return {
      getVars: data.getVars, compareValues_: data.compareValues_, pchisq: dist.pchisq, pnorm: dist.pnorm, qnorm: dist.qnorm,
      qrLeastSquares_: reg.qrLeastSquares_, lineChartBlock: require('./Charts.js').lineChartBlock
    };
  }
  return {
    getVars: getVars, compareValues_: compareValues_, pchisq: pchisq, pnorm: pnorm, qnorm: qnorm,
    qrLeastSquares_: qrLeastSquares_, lineChartBlock: lineChartBlock
  };
}

var LOG_MAX_ITER_ = 100;
var LOG_MAX_CATEGORIES_ = 10;
var LOG_ROC_POINTS_ = 400;
function logName_(v) { return v.label ? v.label : v.name; }

function logSigmoid_(eta) {
  if (eta >= 0) { var e = Math.exp(-eta); return 1 / (1 + e); }
  var e2 = Math.exp(eta);
  return e2 / (1 + e2);
}

/** 로그가능도 Σ[y·η − ln(1 + e^η)] (수치 안정) */
function logLikelihood_(eta, y) {
  var ll = 0;
  for (var i = 0; i < y.length; i++) {
    var e = eta[i];
    ll += y[i] * e - (e > 0 ? e + Math.log1p(Math.exp(-e)) : Math.log1p(Math.exp(e)));
  }
  return ll;
}

/**
 * 로지스틱 회귀 적합. predictors: [{ name, values }], y: 0/1 배열. 절편 포함.
 * 반환: { n, k, B(절편 먼저), cov, se, ll, ll0, p(예측확률), iterations, converged, separation, ... } 또는 { constant|collinear|tooFew|singular }
 */
function logisticFit(predictors, y) {
  var d = logDeps_();
  var n = y.length, k = predictors.length, j, i;
  var n1 = 0;
  for (i = 0; i < n; i++) n1 += y[i];
  var n0 = n - n1;
  if (n1 === 0 || n0 === 0) return { oneClass: true };
  if (n - (k + 1) < 1) return { tooFew: true, n: n, p: k + 1 };

  var means = [], sds = [], cols = [];
  for (j = 0; j < k; j++) {
    var s = 0, v = predictors[j].values;
    for (i = 0; i < n; i++) s += v[i];
    var m = s / n, ss = 0;
    var c = new Array(n);
    for (i = 0; i < n; i++) { c[i] = v[i] - m; ss += c[i] * c[i]; }
    if (!(ss > 0)) return { constant: predictors[j].name };
    means.push(m);
    sds.push(Math.sqrt(ss / Math.max(1, n - 1)));
    cols.push(c);
  }
  var X = [new Array(n).fill(1)].concat(cols);           // 열 단위 설계 행렬 (상수항 포함)
  var p1 = k + 1;

  var beta = new Array(p1).fill(0);
  beta[0] = Math.log(n1 / n0);
  var eta = new Array(n).fill(beta[0]);
  var ll = logLikelihood_(eta, y);
  var ll0 = n1 * Math.log(n1 / n) + n0 * Math.log(n0 / n);
  var converged = false, it, minv = null, singular = false;
  var p = new Array(n);

  for (it = 1; it <= LOG_MAX_ITER_; it++) {
    var w = new Array(n), z = new Array(n);
    for (i = 0; i < n; i++) {
      p[i] = logSigmoid_(eta[i]);
      var wi = Math.max(p[i] * (1 - p[i]), 1e-12);
      w[i] = wi;
      z[i] = eta[i] + (y[i] - p[i]) / wi;
    }
    var sw = w.map(Math.sqrt);
    var wcols = X.map(function (col) { return col.map(function (v, t) { return v * sw[t]; }); });
    var wz = z.map(function (v, t) { return v * sw[t]; });
    var fit = d.qrLeastSquares_(wcols, wz);
    if (!fit.ok) { singular = true; break; }
    minv = fit.minv;
    var step = 1, newBeta, newEta, newLL, tries = 0;
    do {
      newBeta = fit.beta.map(function (b, t) { return beta[t] + step * (b - beta[t]); });
      newEta = new Array(n);
      for (i = 0; i < n; i++) {
        var e = newBeta[0];
        for (j = 0; j < k; j++) e += newBeta[j + 1] * cols[j][i];
        newEta[i] = e;
      }
      newLL = logLikelihood_(newEta, y);
      step /= 2; tries++;
    } while (newLL < ll - 1e-10 && tries < 12);
    var maxDelta = 0;
    for (j = 0; j < p1; j++) maxDelta = Math.max(maxDelta, Math.abs(newBeta[j] - beta[j]));
    var dLL = Math.abs(newLL - ll);
    beta = newBeta; eta = newEta;
    var prevLL = ll;
    ll = newLL;
    if (dLL < 1e-11 * (1 + Math.abs(ll)) && maxDelta < 1e-8) { converged = true; break; }
    if (it > 1 && !(ll > prevLL - 1e-12) && maxDelta < 1e-8) { converged = true; break; }
  }
  for (i = 0; i < n; i++) p[i] = logSigmoid_(eta[i]);
  if (singular) return { singular: true };

  // 수렴한 점의 공분산 (X'WX)^-1 를 다시 계산
  var w2 = p.map(function (pi) { return Math.max(pi * (1 - pi), 1e-300); });
  var sw2 = w2.map(Math.sqrt);
  var fin = d.qrLeastSquares_(X.map(function (col) { return col.map(function (v, t) { return v * sw2[t]; }); }), eta.map(function (e, t) { return e * sw2[t]; }));
  var cov = null;
  if (fin.ok) cov = fin.minv;

  // 원래 척도로 되돌리기: B_raw = A·B_c, A = [[1, −m'], [0, I]]
  var B = beta.slice();
  for (j = 0; j < k; j++) B[0] -= beta[j + 1] * means[j];
  var covRaw = null, se = null;
  if (cov) {
    var A = [];
    for (i = 0; i < p1; i++) { A.push(new Array(p1).fill(0)); A[i][i] = 1; }
    for (j = 0; j < k; j++) A[0][j + 1] = -means[j];
    var AC = A.map(function (row) {
      return cov[0].map(function (_, c2) { var s2 = 0; for (var r = 0; r < p1; r++) s2 += row[r] * cov[r][c2]; return s2; });
    });
    covRaw = AC.map(function (row) {
      return A.map(function (arow) { var s3 = 0; for (var r = 0; r < p1; r++) s3 += row[r] * arow[r]; return s3; });
    });
    se = covRaw.map(function (row, t) { return Math.sqrt(Math.max(0, row[t])); });
  }

  // 분리 의심: 계수가 매우 크거나 예측확률이 0 또는 1에 붙음, 표준오차가 비정상적으로 큼
  // (상수는 중심화한 값, 기울기는 독립변수 1 표준편차당 효과로 비교한다 — 단위가 큰 변수나 큰 평균 때문에 잘못 경고하지 않도록)
  var maxAbs = Math.abs(beta[0]);
  for (j = 0; j < k; j++) maxAbs = Math.max(maxAbs, Math.abs(beta[j + 1]) * sds[j]);
  var extreme = 0;
  for (i = 0; i < n; i++) if ((y[i] === 1 && p[i] > 1 - 1e-9) || (y[i] === 0 && p[i] < 1e-9)) extreme++;
  var maxSe = 0;
  if (se) for (j = 0; j < p1; j++) if (isFinite(se[j])) maxSe = Math.max(maxSe, se[j]);
  var separation = !converged || maxAbs > 25 || extreme > 0 || maxSe > 1e4;

  return {
    n: n, k: k, n1: n1, n0: n0, B: B, cov: covRaw, se: se, ll: ll, ll0: ll0, p: p, iterations: it > LOG_MAX_ITER_ ? LOG_MAX_ITER_ : it,
    converged: converged, separation: separation
  };
}

/** 분위수 (선형 보간, numpy 기본과 같음). sorted 는 오름차순, q 는 0~1 */
function logQuantile_(sorted, q) {
  var pos = (sorted.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
}

/**
 * Hosmer-Lemeshow 적합도: 예측확률의 분위수(기본 10분위)로 집단을 나눠 관측과 기대 빈도를 비교한다 (R의 hoslem.test와 같은 방식).
 * 경계는 구간 (b0, b1], (b1, b2], … 이고 첫 구간은 최솟값을 포함한다. 같은 경계는 합친다.
 * 반환: { groups: [{ n, o1, e1, o0, e0 }], chi2, df, p } (집단이 3개 미만이면 chi2, df, p 는 null)
 */
function hosmerLemeshow(p, y, g, chisqUpper) {
  var n = p.length, i;
  var sorted = p.slice().sort(function (a, b) { return a - b; });
  var breaks = [];
  for (i = 0; i <= g; i++) {
    var b = logQuantile_(sorted, i / g);
    if (!breaks.length || b > breaks[breaks.length - 1]) breaks.push(b);
  }
  var ng = breaks.length - 1;
  var groups = [];
  for (i = 0; i < ng; i++) groups.push({ n: 0, o1: 0, e1: 0 });
  for (i = 0; i < n; i++) {
    var gi = 0;
    while (gi < ng - 1 && p[i] > breaks[gi + 1]) gi++;
    groups[gi].n++; groups[gi].o1 += y[i]; groups[gi].e1 += p[i];
  }
  groups.forEach(function (gr) { gr.o0 = gr.n - gr.o1; gr.e0 = gr.n - gr.e1; });
  var res = { groups: groups, chi2: null, df: null, p: null };
  if (ng >= 3) {
    var chi = 0;
    groups.forEach(function (gr) {
      if (gr.e1 > 1e-12) chi += (gr.o1 - gr.e1) * (gr.o1 - gr.e1) / gr.e1;
      if (gr.e0 > 1e-12) chi += (gr.o0 - gr.e0) * (gr.o0 - gr.e0) / gr.e0;
    });
    res.chi2 = chi; res.df = ng - 2;
    res.p = chisqUpper(chi, res.df);
  }
  return res;
}

/** 분류표: 예측확률이 절단값을 넘으면(>) 사건으로 분류. 반환 { tp, fn, fp, tn, ... } */
function classificationTable(p, y, cutoff) {
  var tp = 0, fn = 0, fp = 0, tn = 0;
  for (var i = 0; i < p.length; i++) {
    var pred = p[i] > cutoff ? 1 : 0;
    if (y[i] === 1) { if (pred === 1) tp++; else fn++; } else { if (pred === 1) fp++; else tn++; }
  }
  var n = p.length;
  return {
    tp: tp, fn: fn, fp: fp, tn: tn,
    sens: tp + fn ? tp / (tp + fn) : null, spec: tn + fp ? tn / (tn + fp) : null,
    pctEvent: tp + fn ? tp / (tp + fn) * 100 : null, pctNon: tn + fp ? tn / (tn + fp) * 100 : null,
    overall: n ? (tp + tn) / n * 100 : null
  };
}

/** 평균 순위 (동순위 평균) */
function logRanks_(values) {
  var n = values.length, idx = new Array(n), i;
  for (i = 0; i < n; i++) idx[i] = i;
  idx.sort(function (a, b) { return values[a] - values[b]; });
  var ranks = new Array(n);
  i = 0;
  while (i < n) {
    var j = i;
    while (j + 1 < n && values[idx[j + 1]] === values[idx[i]]) j++;
    var avg = (i + j) / 2 + 1;
    for (var k = i; k <= j; k++) ranks[idx[k]] = avg;
    i = j + 1;
  }
  return ranks;
}

/**
 * ROC 분석: AUC (동점은 0.5), DeLong(1988) 표준오차, 정규 근사 신뢰구간, 귀무가설 AUC = 0.5 검정(Mann-Whitney 귀무분산).
 * 반환: { auc, se, lo, hi, p, points: [[1−특이도, 민감도], …] } (사건·비사건이 각각 1건 이상이어야 한다)
 */
function rocAnalysis(p, y, conf, d) {
  var n = p.length, i;
  var pos = [], neg = [];
  for (i = 0; i < n; i++) (y[i] === 1 ? pos : neg).push(p[i]);
  var m = pos.length, q = neg.length;
  if (!m || !q) return null;
  var all = logRanks_(p);
  var rpos = logRanks_(pos), rneg = logRanks_(neg);
  var v10 = [], v01 = [], pi = 0, ni = 0;
  for (i = 0; i < n; i++) {
    if (y[i] === 1) { v10.push((all[i] - rpos[pi]) / q); pi++; }
    else { v01.push(1 - (all[i] - rneg[ni]) / m); ni++; }
  }
  var auc = 0;
  v10.forEach(function (v) { auc += v; });
  auc /= m;
  function svar(arr) {
    var mean = 0, s = 0, t;
    for (t = 0; t < arr.length; t++) mean += arr[t];
    mean /= arr.length;
    for (t = 0; t < arr.length; t++) s += (arr[t] - mean) * (arr[t] - mean);
    return arr.length > 1 ? s / (arr.length - 1) : 0;
  }
  var se = Math.sqrt(svar(v10) / m + svar(v01) / q);
  var z = d.qnorm(1 - (1 - conf / 100) / 2);
  var se0 = Math.sqrt((m + q + 1) / (12 * m * q));
  var zt = (auc - 0.5) / se0;
  // 곡선: 임계값을 큰 것부터 내려 가며 (1−특이도, 민감도)
  var order = [];
  for (i = 0; i < n; i++) order.push(i);
  order.sort(function (a, b) { return p[b] - p[a]; });
  var points = [[0, 0]], tp = 0, fp = 0;
  for (i = 0; i < n; i++) {
    if (y[order[i]] === 1) tp++; else fp++;
    if (i === n - 1 || p[order[i + 1]] !== p[order[i]]) points.push([fp / q, tp / m]);
  }
  if (points.length > LOG_ROC_POINTS_) {
    var stride = Math.ceil(points.length / LOG_ROC_POINTS_), thin = [];
    for (i = 0; i < points.length; i += stride) thin.push(points[i]);
    if (thin[thin.length - 1] !== points[points.length - 1]) thin.push(points[points.length - 1]);
    points = thin;
  }
  return {
    auc: auc, se: se, lo: Math.max(0, auc - z * se), hi: Math.min(1, auc + z * se),
    p: 2 * d.pnorm(Math.abs(zt), false), points: points, nPos: m, nNeg: q
  };
}

function runLogistic(dataset, request) {
  var d = logDeps_();
  var dep = d.getVars(dataset, request.dependent);
  if (dep.length !== 1) throw new Error('종속변수는 하나만 선택하세요.');
  var covs = request.covariates && request.covariates.length ? d.getVars(dataset, request.covariates, { numeric: true }) : [];
  var facs = request.factors && request.factors.length ? d.getVars(dataset, request.factors) : [];
  if (covs.length + facs.length === 0) throw new Error('독립변수(공변량 또는 범주형 독립변수)를 하나 이상 선택하세요.');
  var seen = {};
  seen[dep[0].name] = '종속변수';
  covs.concat(facs).forEach(function (v) {
    if (seen[v.name]) throw new Error('변수 "' + v.name + '"를 ' + seen[v.name] + '와(과) 독립변수에 동시에 쓸 수 없습니다.');
    seen[v.name] = '독립변수';
  });
  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  var cutoff = opt.cutoff === undefined || opt.cutoff === null || opt.cutoff === '' ? 0.5 : Number(opt.cutoff);
  if (!isFinite(cutoff) || cutoff <= 0 || cutoff >= 1) throw new Error('"분류 절단값"은 0보다 크고 1보다 작아야 합니다.');

  var all = dep.concat(covs, facs);
  var rowsOk = [];
  for (var i = 0; i < dataset.n; i++) {
    var ok = true;
    for (var a = 0; a < all.length; a++) if (dataset.data[all[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length;
  var pick = function (name) { return rowsOk.map(function (r) { return dataset.data[name][r]; }); };
  var warnings = [];
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건).');

  // 종속변수 코딩 (사용자가 시트에서 0/1로 미리 바꿀 필요가 없다):
  //  ① 기본: 값이 정확히 2가지이면 큰 값이 사건(1) (옵션으로 작은 값)
  //  ② 사건으로 볼 값(eventValues)을 고르면 그 값들이 1, 나머지가 0 (3가지 이상인 변수도 가능)
  //  ③ eventAtLeast 를 주면 그 값 이상이 1 (숫자 변수)
  var evVals = Array.isArray(opt.eventValues) ? opt.eventValues.map(String) : [];
  var atLeast = null;
  if (opt.eventAtLeast !== undefined && opt.eventAtLeast !== null && String(opt.eventAtLeast).trim() !== '') {
    atLeast = Number(opt.eventAtLeast);
    if (!isFinite(atLeast)) throw new Error('"이 값 이상이면 사건"은 숫자여야 합니다.');
  }
  if (evVals.length && atLeast !== null) throw new Error('"사건으로 볼 값"과 "이 값 이상이면 사건"은 둘 중 하나만 쓰세요.');
  var yRaw = pick(dep[0].name), levels = [], seenL = {}, tally = {};
  yRaw.forEach(function (x) {
    var key = typeof x + ':' + x;
    if (!seenL[key]) { seenL[key] = true; levels.push({ value: x, key: key }); tally[key] = 0; }
    tally[key]++;
  });
  levels.sort(d.compareValues_);
  levels.forEach(function (l) {
    var lab = dep[0].valueLabels[String(l.value)];
    l.label = lab !== undefined ? lab : String(l.value);
    l.count = tally[l.key];
  });
  var joinLabels = function (arr) {
    var s = arr.map(function (l) { return l.label; }).join(' / ');
    return s.length > 60 ? s.slice(0, 59) + '…' : s;
  };
  var cleanLab = function (s) { return String(s).replace(/[,;=\n]+/g, ' '); };
  var y, eventLabel, nonEventLabel, codingRows, codingNote, saved, eventIdx = 1, nonIdx = 0;
  if (atLeast !== null) {
    yRaw.forEach(function (x) {
      if (typeof x !== 'number' && !isFinite(Number(x))) throw new Error('"이 값 이상이면 사건"은 숫자 변수에만 쓸 수 있습니다. 종속변수 "' + logName_(dep[0]) + '"에 숫자가 아닌 값("' + x + '")이 있습니다.');
    });
    y = yRaw.map(function (x) { return Number(x) >= atLeast ? 1 : 0; });
    var aTxt = String(Number(atLeast.toPrecision(12)));
    eventLabel = aTxt + ' 이상'; nonEventLabel = aTxt + ' 미만';
    codingRows = [{ orig: nonEventLabel, code: 0 }, { orig: eventLabel, code: 1 }];
    codingNote = '"' + logName_(dep[0]) + '"의 값이 ' + aTxt + ' 이상이면 사건(1), 미만이면 0으로 나눴다.';
    saved = { event: 1, non: 0, labels: '0=' + cleanLab(nonEventLabel) + ', 1=' + cleanLab(eventLabel) };
  } else if (evVals.length) {
    var evSet = {};
    evVals.forEach(function (k) { evSet[k] = true; });
    var evLevels = levels.filter(function (l) { return evSet[String(l.value)]; });
    var nonLevels = levels.filter(function (l) { return !evSet[String(l.value)]; });
    if (!evLevels.length) throw new Error('고른 "사건으로 볼 값"이 사용한 사례의 종속변수 값에 없습니다. 분석 창에서 값 목록을 새로 불러와 다시 고르세요.');
    if (!nonLevels.length) throw new Error('모든 값을 사건으로 골랐습니다. 사건이 아닌 값이 하나 이상 있어야 합니다.');
    y = yRaw.map(function (x) { return evSet[String(x)] ? 1 : 0; });
    eventLabel = joinLabels(evLevels); nonEventLabel = joinLabels(nonLevels);
    codingRows = nonLevels.map(function (l) { return { orig: l.label, code: 0, n: l.count }; }).concat(evLevels.map(function (l) { return { orig: l.label, code: 1, n: l.count }; }));
    codingNote = '"사건으로 볼 값"으로 고른 값(' + eventLabel + ')이 사건(1)이고 나머지 값은 0이다.';
    saved = { event: 1, non: 0, labels: '0=' + cleanLab(nonEventLabel) + ', 1=' + cleanLab(eventLabel) };
  } else {
    if (levels.length !== 2) {
      throw new Error('종속변수 "' + logName_(dep[0]) + '"의 값이 ' + levels.length + '가지입니다. 이분형 로지스틱 회귀는 값이 정확히 2가지여야 합니다.' +
        (levels.length > 2 ? ' 분석 창에서 종속변수를 보내면 나오는 "사건(1)으로 볼 값"에서 사건으로 볼 값을 고르거나, 숫자 변수는 "이 값 이상이면 사건"을 입력하세요.' : ''));
    }
    eventIdx = opt.flip ? 0 : 1; nonIdx = 1 - eventIdx;
    y = yRaw.map(function (x) { return x === levels[eventIdx].value ? 1 : 0; });
    eventLabel = levels[eventIdx].label; nonEventLabel = levels[nonIdx].label;
    codingRows = [{ orig: nonEventLabel, code: 0, n: levels[nonIdx].count }, { orig: eventLabel, code: 1, n: levels[eventIdx].count }];
    codingNote = '기본은 큰 값이 사건이다. 사건으로 볼 값을 직접 고르거나 "작은 값을 사건으로" 옵션을 쓸 수 있다.';
    if (opt.flip) codingNote = '옵션 "작은 값을 사건으로"를 썼다.';
    var vlAll = levels.map(function (l) { return cleanLab(l.value) + '=' + cleanLab(l.label); }).join(', ');
    saved = { event: levels[eventIdx].value, non: levels[nonIdx].value, labels: vlAll };
  }
  if (opt.flip && (atLeast !== null || evVals.length)) warnings.push('옵션 "작은 값을 사건으로"는 사건 값을 직접 고르지 않았을 때만 적용됩니다. 지금은 고른 기준을 따랐습니다.');

  var predictors = covs.map(function (v) { return { name: logName_(v), values: pick(v.name) }; });
  var refNotes = [];
  facs.forEach(function (v) {
    var vals = pick(v.name), cats = [], seenC = {};
    vals.forEach(function (x) { var key = typeof x + ':' + x; if (!seenC[key]) { seenC[key] = true; cats.push({ value: x }); } });
    cats.sort(d.compareValues_);
    cats.forEach(function (c) { var lab = v.valueLabels[String(c.value)]; c.label = lab !== undefined ? lab : String(c.value); });
    if (cats.length < 2) throw new Error('범주형 독립변수 "' + logName_(v) + '"의 범주가 ' + cats.length + '개입니다. 2개 이상이어야 합니다.');
    if (cats.length > LOG_MAX_CATEGORIES_) throw new Error('범주형 독립변수 "' + logName_(v) + '"의 범주가 ' + cats.length + '개입니다 (최대 ' + LOG_MAX_CATEGORIES_ + '개).');
    var refKey = opt.refCategories && opt.refCategories[v.name] !== undefined && opt.refCategories[v.name] !== '' && opt.refCategories[v.name] !== null ? String(opt.refCategories[v.name]) : null;
    var refIdx = opt.lastRef ? cats.length - 1 : 0;
    if (refKey !== null) {
      refIdx = -1;
      cats.forEach(function (c, ci) { if (String(c.value) === refKey) refIdx = ci; });
      if (refIdx < 0) throw new Error('범주형 독립변수 "' + logName_(v) + '"의 기준 범주 "' + refKey + '"이(가) 사용한 사례에 없습니다. 분석 창에서 기준 범주를 다시 고르세요.');
    }
    refNotes.push('"' + logName_(v) + '"의 기준 범주: ' + cats[refIdx].label);
    cats.forEach(function (c, idx) {
      if (idx === refIdx) return;
      predictors.push({ name: logName_(v) + ' = ' + c.label, values: vals.map(function (x) { return x === c.value ? 1 : 0; }) });
    });
  });

  var fit = logisticFit(predictors, y);
  if (fit.oneClass) throw new Error('종속변수의 사건과 비사건이 모두 있어야 합니다 (사용한 사례가 한 범주뿐입니다).');
  if (fit.tooFew) throw new Error('사례 수(' + fit.n + ')가 계수의 수(' + fit.p + ', 상수 포함)보다 많아야 합니다.');
  if (fit.constant) throw new Error('변수 "' + fit.constant + '"의 값이 모두 같아 분석에 쓸 수 없습니다.');
  if (fit.singular) throw new Error('독립변수들이 서로 선형결합이거나 완전히 분리돼 계수를 계산할 수 없습니다. 같은 정보를 담은 변수를 함께 쓰지 않았는지, 범주 하나가 사건과 완전히 겹치지 않는지 확인하세요.');

  if (fit.separation) {
    warnings.push(fit.converged
      ? '계수 또는 표준오차가 매우 크거나 예측확률이 0 또는 1에 붙어 있습니다. 완전·준완전 분리(어떤 변수가 결과를 거의 완벽하게 가름)가 의심됩니다. 이 경우 계수와 표준오차는 믿을 수 없으니 해당 변수를 빼거나 범주를 합쳐 보세요.'
      : '계산이 ' + fit.iterations + '번 반복 안에 수렴하지 않았습니다. 완전·준완전 분리 또는 변수 사이의 공선성이 의심되며 결과를 믿을 수 없습니다.');
  }
  var epv = Math.min(fit.n1, fit.n0) / fit.k;
  if (epv < 10) warnings.push('사건 수와 비사건 수 중 작은 쪽(' + Math.min(fit.n1, fit.n0) + '건)을 독립변수 항 ' + fit.k + '개로 나누면 ' + epv.toFixed(1) + '건으로 관행(10건 이상)에 못 미칩니다. 계수 추정이 불안정할 수 있습니다.');

  var F3 = { type: 'fixed', digits: 3 };
  var blocks = [];
  var lrChi = 2 * (fit.ll - fit.ll0);

  blocks.push({
    type: 'table', caption: '케이스 처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'pct', label: '%', format: { type: 'fixed', digits: 1 } }],
    rows: [
      { label: '분석에 포함', n: n, pct: dataset.n ? n / dataset.n * 100 : null },
      { label: '결측으로 제외', n: dataset.n - n, pct: dataset.n ? (dataset.n - n) / dataset.n * 100 : null },
      { _kind: 'total', label: '합계', n: dataset.n, pct: dataset.n ? 100 : null }
    ],
    notes: ['목록별 결측 제거: 종속변수와 모든 독립변수가 유효한 사례만 사용했다.']
  });

  codingRows.forEach(function (r) { if (r.n === undefined) r.n = r.code === 1 ? fit.n1 : fit.n0; });
  blocks.push({
    type: 'table', caption: '종속변수 코딩 (' + logName_(dep[0]) + ')',
    columns: [{ key: 'orig', label: '원래 값' }, { key: 'code', label: '내부 값', format: { type: 'int' } }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: codingRows,
    notes: ['내부 값 1이 "사건"이다. 계수와 오즈비는 사건(1)이 일어날 로그 오즈에 대한 것이다. ' + codingNote].concat(refNotes)
  });

  blocks.push({
    type: 'table', caption: '모형 검정 (모형 계수의 omnibus 검정)',
    columns: [{ key: 'chi2', label: '카이제곱', format: F3 }, { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }],
    rows: [{ chi2: lrChi, df: fit.k, p: d.pchisq(Math.max(0, lrChi), fit.k, false) }],
    notes: ['귀무가설: 모든 독립변수의 계수가 0이다 (상수만 있는 모형 대비 가능도비 검정). 유의확률이 작으면 모형이 상수만 있는 모형보다 낫다.']
  });

  var cs = 1 - Math.exp(2 * (fit.ll0 - fit.ll) / n);
  var maxCs = 1 - Math.exp(2 * fit.ll0 / n);
  blocks.push({
    type: 'table', caption: '모형 요약',
    columns: [
      { key: 'll0', label: '−2 로그가능도 (상수만)', format: F3 }, { key: 'll1', label: '−2 로그가능도 (모형)', format: F3 },
      { key: 'cs', label: 'Cox와 Snell R 제곱', format: F3 }, { key: 'nk', label: 'Nagelkerke R 제곱', format: F3 },
      { key: 'mf', label: 'McFadden R 제곱', format: F3 }, { key: 'it', label: '반복 횟수', format: { type: 'int' } }
    ],
    rows: [{ ll0: -2 * fit.ll0, ll1: -2 * fit.ll, cs: cs, nk: maxCs > 0 ? cs / maxCs : null, mf: 1 - fit.ll / fit.ll0, it: fit.iterations }],
    notes: ['R 제곱 계열은 선형 회귀의 설명력과 같은 뜻이 아니라 "유사 R 제곱"이다. Nagelkerke는 Cox와 Snell을 최댓값 1로 맞춘 것, McFadden = 1 − ln L(모형) / ln L(상수만).']
  });

  if (opt.hl !== false) {
    var hl = hosmerLemeshow(fit.p, y, 10, function (c, df) { return d.pchisq(c, df, false); });
    blocks.push({
      type: 'table', caption: 'Hosmer-Lemeshow 검정',
      columns: [{ key: 'chi2', label: '카이제곱', format: F3 }, { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }],
      rows: [{ chi2: hl.chi2, df: hl.df, p: hl.p }],
      notes: ['귀무가설: 모형이 자료에 잘 맞는다. 유의확률이 **크면**(예: .05 이상) 적합이 나쁘다는 증거가 없다. 표본이 크면 사소한 차이에도 기각되고 작으면 검정력이 낮다.',
        '예측확률의 10분위(분위수 경계, R의 hoslem.test와 같은 방식)로 집단을 나눴다. 같은 확률의 처리 방식에 따라 다른 통계 프로그램과 값이 조금 다를 수 있다.'].concat(hl.chi2 === null ? ['집단이 3개 미만이라 검정을 계산할 수 없다.'] : [])
    });
    blocks.push({
      type: 'table', caption: 'Hosmer-Lemeshow 검정의 분할표',
      columns: [{ key: 'g', label: '집단', format: { type: 'int' } }, { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'o1', label: '관측 (사건)', format: { type: 'fixed', digits: 0 } }, { key: 'e1', label: '기대 (사건)', format: { type: 'fixed', digits: 2 } },
        { key: 'o0', label: '관측 (비사건)', format: { type: 'fixed', digits: 0 } }, { key: 'e0', label: '기대 (비사건)', format: { type: 'fixed', digits: 2 } }],
      rows: hl.groups.map(function (g, idx) { return { g: idx + 1, n: g.n, o1: g.o1, e1: g.e1, o0: g.o0, e0: g.e0 }; }),
      notes: []
    });
  }

  var cls = classificationTable(fit.p, y, cutoff);
  blocks.push({
    type: 'table', caption: '분류표 (절단값 ' + cutoff + ')',
    columns: [
      { key: 'obs', label: '관측' }, { key: 'p0', label: '예측: ' + nonEventLabel, format: { type: 'int' } },
      { key: 'p1', label: '예측: ' + eventLabel, format: { type: 'int' } }, { key: 'pct', label: '정확도(%)', format: { type: 'fixed', digits: 1 } }
    ],
    rows: [
      { obs: nonEventLabel, p0: cls.tn, p1: cls.fp, pct: cls.pctNon },
      { obs: eventLabel, p0: cls.fn, p1: cls.tp, pct: cls.pctEvent },
      { _kind: 'total', obs: '전체', p0: null, p1: null, pct: cls.overall }
    ],
    notes: ['예측확률이 절단값을 넘으면 사건으로 분류했다. 민감도(사건을 맞힌 비율) = ' + (cls.sens === null ? '-' : (cls.sens * 100).toFixed(1) + '%') + ', 특이도(비사건을 맞힌 비율) = ' + (cls.spec === null ? '-' : (cls.spec * 100).toFixed(1) + '%') + '. 사건이 드물면 절단값 0.5로는 사건을 거의 못 맞힐 수 있어 절단값을 바꿔 본다.']
  });

  var z = d.qnorm(1 - (1 - conf / 100) / 2);
  var coefRows = [];
  var coefNames = ['상수'].concat(predictors.map(function (pr) { return pr.name; }));
  coefNames.forEach(function (nm, j) {
    var b = fit.B[j], se = fit.se ? fit.se[j] : null;
    var wald = se !== null && se > 0 ? (b / se) * (b / se) : null;
    coefRows.push({
      name: nm, B: b, se: se, wald: wald, df: 1, p: wald !== null ? d.pchisq(wald, 1, false) : null,
      or: Math.exp(b), lo: se !== null ? Math.exp(b - z * se) : null, hi: se !== null ? Math.exp(b + z * se) : null
    });
  });
  blocks.push({
    type: 'table', caption: '방정식의 변수',
    columns: [
      { key: 'name', label: '' }, { key: 'B', label: 'B', format: F3 }, { key: 'se', label: '표준오차', format: F3 },
      { key: 'wald', label: 'Wald', format: F3 }, { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
      { key: 'or', label: 'Exp(B) (오즈비)', format: F3 }, { key: 'lo', label: 'Exp(B)의 ' + conf + '% 신뢰구간 하한', format: F3 }, { key: 'hi', label: 'Exp(B)의 ' + conf + '% 신뢰구간 상한', format: F3 }
    ],
    rows: coefRows,
    notes: ['Wald = (B / 표준오차)², 유의확률은 χ²(1). Exp(B)는 독립변수가 1 늘 때(더미 변수는 기준 범주 대비) 사건의 오즈가 몇 배가 되는지이고 1이면 효과가 없다. 신뢰구간은 Wald 방식이다.']
  });

  if (opt.roc !== false) {
    var roc = rocAnalysis(fit.p, y, conf, d);
    if (roc) {
      blocks.push({
        type: 'table', caption: 'ROC 곡선 아래 면적 (AUC)',
        columns: [
          { key: 'auc', label: 'AUC', format: F3 }, { key: 'se', label: '표준오차 (DeLong)', format: F3 },
          { key: 'lo', label: conf + '% 신뢰구간 하한', format: F3 }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: F3 },
          { key: 'p', label: '유의확률 (AUC = .5)', format: { type: 'pvalue' } }
        ],
        rows: [{ auc: roc.auc, se: roc.se, lo: roc.lo, hi: roc.hi, p: roc.p }],
        notes: ['AUC는 사건 사례가 비사건 사례보다 예측확률이 높을 확률(동점은 0.5)이다. .5는 무작위, 1은 완벽한 구분이다 (대략 .7 이상 쓸 만, .8 이상 좋음, .9 이상 매우 좋음). 표준오차는 DeLong(1988), 유의확률은 귀무가설 AUC = .5에 대한 검정이다.',
          '이 AUC는 모형을 만든 같은 자료로 계산해 낙관적일 수 있다 (새 자료에서는 더 낮을 수 있다).']
      });
      blocks.push(d.lineChartBlock('ROC 곡선', '1 − 특이도 (위양성률)', '민감도', [
        { label: 'ROC 곡선 (AUC ' + Number(roc.auc.toFixed(3)) + ')', points: roc.points },
        { label: '기준선', points: [[0, 0], [1, 1]] }
      ]));
    } else {
      warnings.push('사건 또는 비사건이 없어 ROC를 계산하지 않았습니다.');
    }
  }

  var result = { analysis: 'logistic', title: '로지스틱 회귀', blocks: blocks, warnings: warnings };
  if (opt.save) {
    // 사례 순서(데이터셋 전체)에 맞춘 새 열: 제외된 사례는 결측
    var prob = new Array(dataset.n).fill(null), grp = new Array(dataset.n).fill(null);
    rowsOk.forEach(function (r, idx) {
      prob[r] = fit.p[idx];
      grp[r] = fit.p[idx] > cutoff ? saved.event : saved.non;
    });
    result.saveColumns = [
      { name: '예측확률', label: '예측확률 (' + logName_(dep[0]) + ' = ' + eventLabel + ')', scale: '연속', values: prob, fixed: 4 },
      { name: '예측소속', label: '예측 소속 (' + logName_(dep[0]) + ', 절단값 ' + cutoff + ')', scale: '명목', values: grp, valueLabels: saved.labels }
    ];
  }
  // 성향점수 매칭이 같은 적합 결과를 다시 쓰도록 내부 값을 돌려준다 (요청에 returnInternals 가 있을 때만, 호출한 쪽이 지운다)
  if (request.returnInternals) result.internals = { rowsOk: rowsOk, predictors: predictors, y: y, fit: fit };
  return result;
}

if (typeof module !== 'undefined') {
  module.exports = {
    logisticFit: logisticFit,
    hosmerLemeshow: hosmerLemeshow,
    classificationTable: classificationTable,
    rocAnalysis: rocAnalysis,
    runLogistic: runLogistic
  };
}
