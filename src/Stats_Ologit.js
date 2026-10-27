/**
 * 서열(비례오즈) 로지스틱 회귀 — 누적 로짓 모형 logit P(Y ≤ j) = θ_j − xβ. 순수 함수: (dataset, request) → 결과 객체.
 * 부호 규약은 R polr·statsmodels OrderedModel 과 같다: β 가 양수이면 변수가 클수록 더 높은 범주가 될 가능성이 크다.
 * 결과표 구성은 일반적인 보고 형식(처리 요약, 모형 적합 정보, 유사 R², 모수 추정값)를 따르고 비례오즈 가정 점검(Brant 검정)과 분류표를 더했다.
 *
 * 추정: 최대우도. 독립변수를 표준화해 뉴턴-랩슨(해석적 기울기·헤세 행렬)으로 풀고, 문턱값이 오름차순이 아니거나 로그가능도가 줄면 보폭을 줄인다.
 *   표준오차는 관측 정보 행렬의 역행렬, 원래 척도로 되돌릴 때 선형 변환으로 공분산을 옮긴다. 완전·준완전 분리와 미수렴은 경고한다.
 * Brant 검정(1990): 문턱값마다 이분형 로지스틱(Y > j)을 따로 적합하고 계수 차이의 Wald 검정을 한다. 계수 사이의 공분산은
 *   V_ml = (X'W_m X)⁻¹ X'[π_l(1 − π_m)]X (X'W_l X)⁻¹ (m ≤ l).
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function olDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var reg = require('./Stats_Regression.js');
    return {
      getVars: data.getVars, compareValues_: data.compareValues_, pchisq: dist.pchisq, pnorm: dist.pnorm, qnorm: dist.qnorm,
      regDummies_: reg.regDummies_, logisticFit: require('./Stats_Logistic.js').logisticFit
    };
  }
  return {
    getVars: getVars, compareValues_: compareValues_, pchisq: pchisq, pnorm: pnorm, qnorm: qnorm,
    regDummies_: regDummies_, logisticFit: logisticFit
  };
}

var OL_MIN_CATS_ = 3;
var OL_MAX_CATS_ = 10;
var OL_MAX_ITER_ = 100;
var OL_SMALL_CAT_ = 5;
function olName_(v) { return v.label ? v.label : v.name; }

function olF_(x) { return x >= 0 ? 1 / (1 + Math.exp(-x)) : Math.exp(x) / (1 + Math.exp(x)); }

/** 대칭 양정치 행렬의 Cholesky 분해로 A⁻¹ 을 구한다. 양정치가 아니면 null */
function olInvert_(A) {
  var n = A.length, i, j, k, L = [];
  for (i = 0; i < n; i++) L.push(new Array(n).fill(0));
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      var s = A[i][j];
      for (k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) { if (!(s > 1e-300)) return null; L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j];
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

/** 로그가능도, 기울기, 헤세 행렬(요청 시). phi = [β(k), θ(J−1)], X: 사례 × k (표준화), y: 범주 번호 0..J−1 */
function olEval_(phi, X, y, J, k, needDeriv) {
  var n = y.length, p = k + J - 1, i, a, b;
  var theta = phi.slice(k);
  var ll = 0, g = null, H = null;
  if (needDeriv) {
    g = new Array(p).fill(0);
    H = [];
    for (a = 0; a < p; a++) H.push(new Array(p).fill(0));
  }
  for (i = 0; i < n; i++) {
    var eta = 0;
    for (a = 0; a < k; a++) eta += phi[a] * X[i][a];
    var c = y[i];
    var up = c < J - 1, lo = c > 0;
    var Fa = up ? olF_(theta[c] - eta) : 1, Fb = lo ? olF_(theta[c - 1] - eta) : 0;
    var P = Fa - Fb;
    if (!(P > 1e-300)) return { ll: -Infinity };
    ll += Math.log(P);
    if (!needDeriv) continue;
    var ga = up ? Fa * (1 - Fa) : 0, gb = lo ? Fb * (1 - Fb) : 0;
    var ha = up ? ga * (1 - 2 * Fa) : 0, hb = lo ? gb * (1 - 2 * Fb) : 0;
    // 편미분 (P 의 η, θ_c, θ_{c-1} 에 대한 1·2계 도함수)
    var Pe = -(ga - gb), Pc = ga, Pd = -gb;
    var Pee = ha - hb, Pec = -ha, Ped = hb, Pcc = ha, Pdd = -hb;
    var gE = Pe / P, gC = Pc / P, gD = Pd / P;
    var hEE = Pee / P - gE * gE, hEC = Pec / P - gE * gC, hED = Ped / P - gE * gD, hCC = Pcc / P - gC * gC, hDD = Pdd / P - gD * gD, hCD = -gC * gD;
    for (a = 0; a < k; a++) {
      g[a] += gE * X[i][a];
      for (b = 0; b <= a; b++) H[a][b] += hEE * X[i][a] * X[i][b];
    }
    var ic = k + c, id = k + c - 1;
    if (up) g[ic] += gC;
    if (lo) g[id] += gD;
    for (a = 0; a < k; a++) {
      if (up) H[ic][a] += hEC * X[i][a];
      if (lo) H[id][a] += hED * X[i][a];
    }
    if (up) H[ic][ic] += hCC;
    if (lo) H[id][id] += hDD;
    if (up && lo) { H[ic][id] += hCD; H[id][ic] += hCD; }
  }
  if (needDeriv) {
    for (a = 0; a < p; a++) for (b = a + 1; b < p; b++) if (b < k && a < k) H[a][b] = H[b][a];
    // 아래삼각만 채운 β-β 블록을 대칭으로, θ-β 블록은 이미 [θ][β] 에 있으므로 [β][θ] 로 복사
    for (a = 0; a < k; a++) for (b = k; b < p; b++) H[a][b] = H[b][a];
  }
  return { ll: ll, g: g, H: H };
}

/**
 * 서열 로지스틱 적합. predictors: [{ name, values }], y: 범주 번호 배열(0..J−1), J: 범주 수.
 * 반환: { n, k, J, beta, theta, cov (모수 [β, θ] 의 공분산, 원래 척도), se (같은 순서), ll, ll0, probs (n × J), iterations, converged, separation }
 *      또는 { constant|tooFew|singular|emptyCategory }
 */
function ordinalFit(predictors, y, J) {
  var n = y.length, k = predictors.length, i, j, a, b;
  var counts = new Array(J).fill(0);
  y.forEach(function (c) { counts[c]++; });
  for (j = 0; j < J; j++) if (counts[j] === 0) return { emptyCategory: j };
  var p = k + J - 1;
  if (n - p < 1) return { tooFew: true, n: n, p: p };
  var means = [], sds = [], X = [];
  for (i = 0; i < n; i++) X.push(new Array(k));
  for (a = 0; a < k; a++) {
    var s = 0, v = predictors[a].values;
    for (i = 0; i < n; i++) s += v[i];
    var m = s / n, ss = 0;
    for (i = 0; i < n; i++) ss += (v[i] - m) * (v[i] - m);
    if (!(ss > 0)) return { constant: predictors[a].name };
    var sd = Math.sqrt(ss / Math.max(1, n - 1));
    means.push(m); sds.push(sd);
    for (i = 0; i < n; i++) X[i][a] = (v[i] - m) / sd;
  }
  // 시작값: β = 0, θ = 누적 비율의 로짓
  var phi = new Array(p).fill(0), cum = 0, ll0 = 0;
  for (j = 0; j < J - 1; j++) { cum += counts[j]; phi[k + j] = Math.log(cum / (n - cum)); }
  counts.forEach(function (c) { ll0 += c * Math.log(c / n); });
  var cur = olEval_(phi, X, y, J, k, true);
  var converged = false, it, singular = false;
  for (it = 1; it <= OL_MAX_ITER_; it++) {
    var negH = cur.H.map(function (row) { return row.map(function (v) { return -v; }); });
    var inv = olInvert_(negH), ridge = 0;
    while (!inv && ridge < 1e6) {
      ridge = ridge === 0 ? 1e-8 : ridge * 100;
      inv = olInvert_(negH.map(function (row, r) { return row.map(function (v, c) { return r === c ? v + ridge : v; }); }));
    }
    if (!inv) { singular = true; break; }
    var dir = new Array(p).fill(0);
    for (a = 0; a < p; a++) for (b = 0; b < p; b++) dir[a] += inv[a][b] * cur.g[b];
    var step = 1, tries = 0, next = null, nphi = null;
    do {
      nphi = phi.map(function (v, t) { return v + step * dir[t]; });
      var ordered = true;
      for (j = 1; j < J - 1; j++) if (!(nphi[k + j] - nphi[k + j - 1] > 1e-9)) { ordered = false; break; }
      next = ordered ? olEval_(nphi, X, y, J, k, false) : { ll: -Infinity };
      if (next.ll >= cur.ll - 1e-12) break;
      step /= 2; tries++;
    } while (tries < 30);
    if (!(next.ll >= cur.ll - 1e-12)) { converged = it > 1; break; }
    var maxDelta = 0;
    for (a = 0; a < p; a++) maxDelta = Math.max(maxDelta, Math.abs(nphi[a] - phi[a]));
    var dLL = next.ll - cur.ll;
    phi = nphi;
    cur = olEval_(phi, X, y, J, k, true);
    if (maxDelta < 1e-9 && Math.abs(dLL) < 1e-11 * (1 + Math.abs(cur.ll))) { converged = true; break; }
  }
  if (singular) return { singular: true };
  var negH2 = cur.H.map(function (row) { return row.map(function (v) { return -v; }); });
  var covS = olInvert_(negH2);
  // 원래 척도: β_raw = β_s / sd, θ_raw = θ_s + Σ β_raw m
  var beta = [], theta = [];
  for (a = 0; a < k; a++) beta.push(phi[a] / sds[a]);
  for (j = 0; j < J - 1; j++) {
    var t = phi[k + j];
    for (a = 0; a < k; a++) t += beta[a] * means[a];
    theta.push(t);
  }
  var cov = null, se = null;
  if (covS) {
    var A = [];
    for (a = 0; a < p; a++) { A.push(new Array(p).fill(0)); A[a][a] = 1; }
    for (a = 0; a < k; a++) A[a][a] = 1 / sds[a];
    for (j = 0; j < J - 1; j++) for (a = 0; a < k; a++) A[k + j][a] = means[a] / sds[a];
    var AC = A.map(function (row) { return covS[0].map(function (_, c2) { var s2 = 0; for (var r = 0; r < p; r++) s2 += row[r] * covS[r][c2]; return s2; }); });
    cov = AC.map(function (row) { return A.map(function (arow) { var s3 = 0; for (var r = 0; r < p; r++) s3 += row[r] * arow[r]; return s3; }); });
    se = cov.map(function (row, t2) { return Math.sqrt(Math.max(0, row[t2])); });
  }
  // 범주별 확률
  var probs = [];
  for (i = 0; i < n; i++) {
    var eta = 0;
    for (a = 0; a < k; a++) eta += phi[a] * X[i][a];
    var row = [], prev = 0;
    for (j = 0; j < J; j++) {
      var Fj = j < J - 1 ? olF_(phi[k + j] - eta) : 1;
      row.push(Fj - prev); prev = Fj;
    }
    probs.push(row);
  }
  var maxStd = 0;
  for (a = 0; a < k; a++) maxStd = Math.max(maxStd, Math.abs(phi[a]));
  var maxSe = 0;
  if (se) se.forEach(function (v) { if (isFinite(v)) maxSe = Math.max(maxSe, v); });
  var separation = !converged || maxStd > 25 || maxSe > 1e4 || !covS;
  return { n: n, k: k, J: J, beta: beta, theta: theta, cov: cov, se: se, ll: cur.ll, ll0: ll0, probs: probs, counts: counts,
    iterations: it > OL_MAX_ITER_ ? OL_MAX_ITER_ : it, converged: converged, separation: separation };
}

/**
 * Brant 검정. predictors, y(0..J−1), J. 문턱값 j = 0..J−2 마다 이분형 로지스틱(Y > j)을 적합한다.
 * 반환: { ok, splitB: [문턱값][변수] 계수, omnibus: { chi2, df, p }, perVar: [{ chi2, df, p }] } 또는 { ok: false, reason }
 */
function brantTest(predictors, y, J) {
  var d = olDeps_();
  var n = y.length, k = predictors.length, S = J - 1, m, l, i, a, b;
  var fits = [];
  for (m = 0; m < S; m++) {
    var yb = y.map(function (c) { return c > m ? 1 : 0; });
    var f = d.logisticFit(predictors, yb);
    if (f.oneClass || f.tooFew || f.constant || f.singular || !f.cov) return { ok: false, reason: '문턱값 ' + (m + 1) + '의 이분형 로지스틱을 계산할 수 없습니다.' };
    if (!f.converged || f.separation) return { ok: false, reason: '문턱값 ' + (m + 1) + '의 이분형 로지스틱이 수렴하지 않았거나 분리가 의심됩니다.' };
    fits.push(f);
  }
  var q = k + 1;
  var Xrow = function (t) { var r = [1]; for (var c = 0; c < k; c++) r.push(predictors[c].values[t]); return r; };
  // 사례마다 x x' 의 위삼각 성분을 한 번만 만들어 두고 (q(q+1)/2 개), 문턱값 쌍마다 가중합만 계산한다
  var npk = q * (q + 1) / 2, pack = new Float64Array(n * npk);
  for (i = 0; i < n; i++) {
    var xr = Xrow(i), o = i * npk, pp = 0;
    for (a = 0; a < q; a++) for (b = a; b < q; b++) pack[o + pp++] = xr[a] * xr[b];
  }
  // 블록 공분산 V[m][l] (q × q)
  var V = [];
  for (m = 0; m < S; m++) V.push(new Array(S));
  for (m = 0; m < S; m++) {
    for (l = m; l < S; l++) {
      var acc = new Float64Array(npk), pmv = fits[m].p, plv = fits[l].p;
      for (i = 0; i < n; i++) {
        var w = plv[i] - pmv[i] * plv[i];          // m ≤ l 이면 π_l (1 − π_m)
        var o2 = i * npk;
        for (var pp2 = 0; pp2 < npk; pp2++) acc[pp2] += w * pack[o2 + pp2];
      }
      var Sml = [], ppi = 0;
      for (a = 0; a < q; a++) Sml.push(new Array(q).fill(0));
      for (a = 0; a < q; a++) for (b = a; b < q; b++) { Sml[a][b] = acc[ppi]; Sml[b][a] = acc[ppi]; ppi++; }
      var left = fits[m].cov, right = fits[l].cov;
      var T1 = left.map(function (row) { return Sml[0].map(function (_, c2) { var s = 0; for (var r = 0; r < q; r++) s += row[r] * Sml[r][c2]; return s; }); });
      V[m][l] = T1.map(function (row) { return right[0].map(function (_, c2) { var s = 0; for (var r = 0; r < q; r++) s += row[r] * right[r][c2]; return s; }); });
      V[l][m] = V[m][l][0].map(function (_, c2) { return V[m][l].map(function (row) { return row[c2]; }); });
    }
  }
  var splitB = fits.map(function (f) { return f.B.slice(1); });
  // 기울기만의 공분산과 차이 대비
  var beta = [], idxOf = function (mm, v) { return mm * k + v; };
  for (m = 0; m < S; m++) for (a = 0; a < k; a++) beta.push(splitB[m][a]);
  var cov = [];
  for (a = 0; a < S * k; a++) cov.push(new Array(S * k).fill(0));
  for (m = 0; m < S; m++) for (l = 0; l < S; l++) for (a = 0; a < k; a++) for (b = 0; b < k; b++) cov[idxOf(m, a)][idxOf(l, b)] = V[m][l][a + 1][b + 1];
  var wald = function (vars) {
    var rowsD = [];
    vars.forEach(function (v) {
      for (var mm = 0; mm < S - 1; mm++) {
        var dr = new Array(S * k).fill(0);
        dr[idxOf(mm, v)] = 1; dr[idxOf(mm + 1, v)] = -1;
        rowsD.push(dr);
      }
    });
    var diff = rowsD.map(function (dr) { var s = 0; for (var t = 0; t < dr.length; t++) s += dr[t] * beta[t]; return s; });
    var DV = rowsD.map(function (dr) { return cov[0].map(function (_, c2) { var s = 0; for (var t = 0; t < dr.length; t++) s += dr[t] * cov[t][c2]; return s; }); });
    var M = rowsD.map(function (dr, r1) { return rowsD.map(function (dr2, r2) { var s = 0; for (var t = 0; t < dr2.length; t++) s += DV[r1][t] * dr2[t]; return s; }); });
    var Mi = olInvert_(M);
    if (!Mi) return { chi2: null, df: rowsD.length, p: null };
    var chi2 = 0;
    for (var r1 = 0; r1 < diff.length; r1++) for (var r2 = 0; r2 < diff.length; r2++) chi2 += diff[r1] * Mi[r1][r2] * diff[r2];
    return { chi2: chi2, df: rowsD.length, p: d.pchisq(chi2, rowsD.length, false) };
  };
  var all = [];
  for (a = 0; a < k; a++) all.push(a);
  var perVar = all.map(function (v) { return wald([v]); });
  return { ok: true, splitB: splitB, omnibus: wald(all), perVar: perVar };
}

/* ---------- 실행 ---------- */

function runOrdinalLogit(dataset, request) {
  var d = olDeps_();
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
  var depName = olName_(dep[0]);

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

  // 종속변수의 범주: 값 순서 = 범주 순서
  var yRaw = pick(dep[0].name), levels = [], seenL = {}, tally = {};
  yRaw.forEach(function (x) {
    var key = typeof x + ':' + x;
    if (!seenL[key]) { seenL[key] = true; levels.push({ value: x, key: key }); tally[key] = 0; }
    tally[key]++;
  });
  levels.sort(d.compareValues_);
  var J = levels.length;
  if (J < OL_MIN_CATS_ || J > OL_MAX_CATS_) {
    throw new Error('종속변수 "' + depName + '"의 값이 ' + J + '가지입니다. 서열 로지스틱 회귀는 값이 ' + OL_MIN_CATS_ + '~' + OL_MAX_CATS_ + '가지여야 합니다.' +
      (J === 2 ? ' 2가지이면 이분형 로지스틱 회귀를 쓰세요.' : J > OL_MAX_CATS_ ? ' 점수가 연속에 가까우면 선형 회귀를 쓰거나 범주를 합쳐 보세요.' : ''));
  }
  var pos = {};
  levels.forEach(function (l, idx) {
    pos[l.key] = idx;
    var lab = dep[0].valueLabels[String(l.value)];
    l.label = lab !== undefined ? lab : String(l.value);
    l.count = tally[l.key];
  });
  var y = yRaw.map(function (x) { return pos[typeof x + ':' + x]; });
  var small = levels.filter(function (l) { return l.count < OL_SMALL_CAT_; });
  if (small.length) warnings.push('사례가 ' + OL_SMALL_CAT_ + '건 미만인 범주가 있습니다 (' + small.map(function (l) { return l.label + ' ' + l.count + '건'; }).join(', ') + '). 추정이 불안정할 수 있으니 이웃한 범주를 합치는 것을 고려하세요.');

  var predictors = covs.map(function (v) { return { name: olName_(v), values: pick(v.name) }; });
  var refNotes = [];
  facs.forEach(function (v) {
    var dm = d.regDummies_(v, pick(v.name), opt.lastRef);
    refNotes.push(dm.note);
    dm.predictors.forEach(function (pr) { predictors.push(pr); });
  });

  var fit = ordinalFit(predictors, y, J);
  if (fit.tooFew) throw new Error('사례 수(' + fit.n + ')가 모수의 수(' + fit.p + ', 문턱값 포함)보다 많아야 합니다.');
  if (fit.constant) throw new Error('변수 "' + fit.constant + '"의 값이 모두 같아 분석에 쓸 수 없습니다.');
  if (fit.singular) throw new Error('독립변수들이 서로 선형결합이거나 완전히 분리돼 계수를 계산할 수 없습니다. 같은 정보를 담은 변수를 함께 쓰지 않았는지 확인하세요.');
  if (fit.emptyCategory !== undefined) throw new Error('범주에 사례가 없습니다.');
  if (fit.separation) {
    warnings.push(fit.converged
      ? '계수 또는 표준오차가 매우 큽니다. 완전·준완전 분리(어떤 변수가 결과를 거의 완벽하게 가름)가 의심됩니다. 이 경우 계수와 표준오차는 믿을 수 없으니 해당 변수를 빼거나 범주를 합쳐 보세요.'
      : '계산이 ' + fit.iterations + '번 반복 안에 수렴하지 않았습니다. 완전·준완전 분리 또는 변수 사이의 공선성이 의심되며 결과를 믿을 수 없습니다.');
  }
  var nParams = fit.k + J - 1;
  if (n < 10 * nParams) warnings.push('사례 수(' + n + '건)가 모수 수(' + nParams + '개, 문턱값 포함)의 10배에 못 미칩니다. 추정이 불안정할 수 있습니다.');

  var F3 = { type: 'fixed', digits: 3 };
  var blocks = [];
  var lrChi = 2 * (fit.ll - fit.ll0), lrDf = fit.k;

  var sumRows = levels.map(function (l) { return { label: depName + ' = ' + l.label, n: l.count, pct: l.count / n * 100 }; });
  sumRows.push({ label: '분석에 포함', n: n, pct: dataset.n ? n / dataset.n * 100 : null });
  sumRows.push({ label: '결측으로 제외', n: dataset.n - n, pct: dataset.n ? (dataset.n - n) / dataset.n * 100 : null });
  sumRows.push({ _kind: 'total', label: '합계', n: dataset.n, pct: dataset.n ? 100 : null });
  blocks.push({
    type: 'table', caption: '케이스 처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'pct', label: '%', format: { type: 'fixed', digits: 1 } }],
    rows: sumRows,
    notes: ['종속변수의 범주는 값이 작은 것부터 큰 순서(' + levels.map(function (l) { return l.label; }).join(' < ') + ')이다. 앞 세 줄의 % 는 분석에 포함한 사례 중 비율이다.']
  });

  blocks.push({
    type: 'table', caption: '모형 적합 정보',
    columns: [{ key: 'model', label: '모형' }, { key: 'm2ll', label: '−2 로그 가능도', format: F3 }, { key: 'chi2', label: '카이제곱', format: F3 }, { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }],
    rows: [
      { model: '절편만', m2ll: -2 * fit.ll0 },
      { model: '최종', m2ll: -2 * fit.ll, chi2: lrChi, df: lrDf, p: d.pchisq(Math.max(0, lrChi), lrDf, false) }
    ],
    notes: ['카이제곱은 두 모형의 −2 로그 가능도 차이(우도비 검정)이다. 유의하면 독립변수들이 함께 결과 범주를 설명한다. 반복 ' + fit.iterations + '회' + (fit.converged ? ' 만에 수렴' : ', 수렴 안 됨') + '.']
  });

  var cs = 1 - Math.exp(2 * (fit.ll0 - fit.ll) / n);
  var csMax = 1 - Math.exp(2 * fit.ll0 / n);
  blocks.push({
    type: 'table', caption: '유사 R 제곱',
    columns: [{ key: 'cs', label: 'Cox와 Snell', format: F3 }, { key: 'nk', label: 'Nagelkerke', format: F3 }, { key: 'mf', label: 'McFadden', format: F3 }],
    rows: [{ cs: cs, nk: csMax > 0 ? cs / csMax : null, mf: 1 - fit.ll / fit.ll0 }],
    notes: ['선형 회귀의 R²와 같은 뜻이 아니라 모형의 설명력을 가늠하는 참고 값이다.']
  });

  // 모수 추정값
  var z = d.qnorm((1 - conf / 100) / 2, false);
  var paramRows = [], idx = 0;
  var addRow = function (name, b, se, isLoc) {
    var w = se > 0 ? b / se : null;
    paramRows.push({
      name: name, B: b, se: se, wald: w === null ? null : w * w, p: w === null ? null : 2 * d.pnorm(-Math.abs(w), true),
      lo: se > 0 ? b - z * se : null, hi: se > 0 ? b + z * se : null,
      or: isLoc ? Math.exp(b) : null, orLo: isLoc && se > 0 ? Math.exp(b - z * se) : null, orHi: isLoc && se > 0 ? Math.exp(b + z * se) : null
    });
  };
  for (var j = 0; j < J - 1; j++) addRow('문턱값 [' + depName + ' ≤ ' + levels[j].label + ']', fit.theta[j], fit.se ? fit.se[fit.k + j] : null, false);
  predictors.forEach(function (pr, t) { addRow(pr.name, fit.beta[t], fit.se ? fit.se[t] : null, true); });
  var pNotes = ['모형: logit P(' + depName + ' ≤ 범주 j) = 문턱값 j − (독립변수 × B). B가 양수이면 그 변수가 클수록 더 높은 범주가 될 가능성이 큽니다. Exp(B)는 변수가 1 늘 때 "더 높은 범주(누적)에 속할 오즈"가 몇 배가 되는지이고 1이면 효과가 없다.',
    'Wald = (B / 표준오차)², 유의확률은 정규 근사의 양측 검정이다.'].concat(refNotes);
  blocks.push({
    type: 'table', caption: '모수 추정값',
    columns: [
      { key: 'name', label: '' }, { key: 'B', label: 'B', format: F3 }, { key: 'se', label: '표준오차', format: F3 }, { key: 'wald', label: 'Wald', format: F3 }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
      { key: 'lo', label: 'B의 ' + conf + '% 신뢰구간 하한', format: F3 }, { key: 'hi', label: 'B의 ' + conf + '% 신뢰구간 상한', format: F3 },
      { key: 'or', label: 'Exp(B)', format: F3 }, { key: 'orLo', label: 'Exp(B) 하한', format: F3 }, { key: 'orHi', label: 'Exp(B) 상한', format: F3 }
    ],
    rows: paramRows, notes: pNotes
  });

  // 비례오즈 가정 점검 (Brant)
  var brant = null;
  if (opt.parallel !== false) {
    brant = brantTest(predictors, y, J);
    if (brant.ok) {
      var bcols = [{ key: 'name', label: '' }];
      for (var s = 0; s < J - 1; s++) bcols.push({ key: 'b' + s, label: 'B [' + depName + ' > ' + levels[s].label + ']', format: F3 });
      bcols.push({ key: 'chi2', label: 'Wald 카이제곱', format: F3 }, { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } });
      var brows = predictors.map(function (pr, t) {
        var row = { name: pr.name, chi2: brant.perVar[t].chi2, df: brant.perVar[t].df, p: brant.perVar[t].p };
        brant.splitB.forEach(function (bs, s2) { row['b' + s2] = bs[t]; });
        return row;
      });
      brows.push({ _kind: 'total', name: '전체 (Brant)', chi2: brant.omnibus.chi2, df: brant.omnibus.df, p: brant.omnibus.p });
      blocks.push({
        type: 'table', caption: '비례오즈 가정 점검 (Brant 검정)',
        columns: bcols, rows: brows,
        notes: ['문턱값마다 "그 범주보다 높은가"를 이분형 로지스틱 회귀로 따로 추정한 계수(B)를 나란히 놓았다. 가정이 맞으면 계수가 문턱값마다 비슷하고, Wald 검정의 유의확률이 .05 이상이다. 유의하면 그 변수의 효과가 문턱값마다 달라 비례오즈 가정이 의심된다.',
          'Brant(1990)의 Wald 검정이다 (우도비 방식의 평행선 검정과는 값이 다를 수 있다). 표본이 크면 사소한 차이에도 유의해지므로 계수의 차이가 실제로 큰지도 보세요.']
      });
      if (brant.omnibus.p !== null && brant.omnibus.p < 0.05) {
        var bad = [];
        predictors.forEach(function (pr, t) { if (brant.perVar[t].p !== null && brant.perVar[t].p < 0.05) bad.push(pr.name); });
        warnings.push('비례오즈 가정이 의심됩니다 (Brant 전체 p = ' + (brant.omnibus.p < 0.001 ? '<.001' : brant.omnibus.p.toFixed(3)) + (bad.length ? ', 문제가 되는 변수: ' + bad.join(', ') : '') + '). 계수가 문턱값마다 크게 다르면 이 모형의 단일 B 대신 범주를 합친 이분형 로지스틱 회귀를 함께 보세요.');
      }
    } else {
      warnings.push('비례오즈 가정 점검(Brant)을 계산하지 못했습니다: ' + brant.reason);
    }
  }

  // 분류표와 예측
  var predCat = fit.probs.map(function (row) {
    var best = 0;
    row.forEach(function (v, t) { if (v > row[best] + 1e-15) best = t; });
    return best;
  });
  if (opt.classify) {
    var T = [];
    for (var u = 0; u < J; u++) T.push(new Array(J).fill(0));
    for (var t2 = 0; t2 < n; t2++) T[y[t2]][predCat[t2]]++;
    var hit = 0;
    T.forEach(function (row, c) { hit += row[c]; });
    var ccols = [{ key: 'name', label: '실제 ＼ 예측' }];
    levels.forEach(function (l, c) { ccols.push({ key: 'c' + c, label: l.label, format: { type: 'int' } }); });
    ccols.push({ key: 'sum', label: '합계', format: { type: 'int' } }, { key: 'pct', label: '정분류 %', format: { type: 'fixed', digits: 1 } });
    var crows = T.map(function (row, c) {
      var sum = row.reduce(function (s2, v) { return s2 + v; }, 0);
      var r = { name: levels[c].label, sum: sum, pct: sum ? row[c] / sum * 100 : null };
      row.forEach(function (v, c2) { r['c' + c2] = v; });
      return r;
    });
    var crow = { _kind: 'total', name: '전체', sum: n, pct: hit / n * 100 };
    levels.forEach(function (l, c) { var s3 = 0; T.forEach(function (row) { s3 += row[c]; }); crow['c' + c] = s3; });
    crows.push(crow);
    blocks.push({
      type: 'table', caption: '분류표 (가장 가능성이 높은 범주로 예측)',
      columns: ccols, rows: crows,
      notes: ['각 사례를 확률이 가장 큰 범주로 예측한 결과이다. 전체 정분류 % 는 ' + (hit / n * 100).toFixed(1) + '이다. 서열 범주에서는 정확히 맞히지 못해도 이웃 범주로 틀리는 경우가 많으니 대각선 근처에 몰리는지도 보세요.']
    });
  }

  var result = { analysis: 'ologit', title: '서열 로지스틱 회귀 (비례오즈)', blocks: blocks, warnings: warnings };
  if (opt.save) {
    var vlAll = levels.map(function (l) { return String(l.value).replace(/[,;=\n]+/g, ' ') + '=' + String(l.label).replace(/[,;=\n]+/g, ' '); }).join(', ');
    var predCol = new Array(dataset.n).fill(null);
    var probCols = levels.map(function () { return new Array(dataset.n).fill(null); });
    rowsOk.forEach(function (r, t3) {
      predCol[r] = levels[predCat[t3]].value;
      fit.probs[t3].forEach(function (pv, c) { probCols[c][r] = pv; });
    });
    result.saveColumns = [{ name: '예측범주', label: '예측 범주 (' + depName + ', 확률이 가장 큰 범주)', scale: '서열', values: predCol, valueLabels: vlAll }]
      .concat(levels.map(function (l, c) { return { name: '확률_' + String(l.value).replace(/[^0-9A-Za-z가-힣_]+/g, '_'), label: '예측확률 (' + depName + ' = ' + l.label + ')', scale: '연속', values: probCols[c], fixed: 4 }; }));
  }
  return result;
}

if (typeof module !== 'undefined') {
  module.exports = { ordinalFit: ordinalFit, brantTest: brantTest, runOrdinalLogit: runOrdinalLogit };
}
