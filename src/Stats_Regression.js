/**
 * OLS 다중회귀. 순수 함수: (dataset, request) → 결과 객체.
 * 결과표 구성: 모형 요약, 분산분석, 계수, 공선성 진단.
 *
 * 수치 방식: 독립변수를 평균으로 중심화한 뒤 Householder QR로 푼다 (정규방정식을 직접 풀지 않는다).
 * 상수항과 그 표준오차는 중심화한 결과에서 복원한다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function regDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var desc = require('./Stats_Descriptive.js');
    return {
      getVars: data.getVars, compareValues_: data.compareValues_,
      pf: dist.pf, pt: dist.pt, qt: dist.qt, centralSums: desc.centralSums
    };
  }
  return { getVars: getVars, compareValues_: compareValues_, pf: pf, pt: pt, qt: qt, centralSums: centralSums };
}

function regName_(v) { return v.label ? v.label : v.name; }

var REG_RANK_TOL_ = 1e-7;       // 중심화한 열의 잔여 크기 / 원래 크기 가 이보다 작으면 완전 공선성
var REG_MAX_CATEGORIES_ = 10;   // 범주형 독립변수의 최대 범주 수
var REG_VIF_WARN_ = 10;

/* ---------- 선형대수 ---------- */

/**
 * 열 단위 행렬 cols (n × p) 와 y 의 Householder QR 최소제곱 (절편 없음).
 * 반환: { ok, badCol } 또는 { ok, beta, minv }  (minv = (X'X)^-1, p × p)
 * 열 j 의 잔여 크기가 원래 크기의 REG_RANK_TOL_ 이하이면 { ok: false, badCol: j }.
 */
function qrLeastSquares_(cols, y) {
  var p = cols.length, n = y.length, i, j, k;
  var A = cols.map(function (c) { return c.slice(); });
  var b = y.slice();
  var norms0 = cols.map(function (c) {
    var s = 0;
    for (var t = 0; t < n; t++) s += c[t] * c[t];
    return Math.sqrt(s);
  });
  for (j = 0; j < p; j++) {
    var col = A[j], norm = 0;
    for (i = j; i < n; i++) norm += col[i] * col[i];
    norm = Math.sqrt(norm);
    if (!(norm > REG_RANK_TOL_ * norms0[j])) return { ok: false, badCol: j };
    var alpha = col[j] > 0 ? -norm : norm;
    var v = new Array(n - j);
    for (i = j; i < n; i++) v[i - j] = col[i];
    v[0] -= alpha;
    var vv = 0;
    for (i = 0; i < v.length; i++) vv += v[i] * v[i];
    for (k = j + 1; k < p; k++) {
      var dot = 0;
      for (i = j; i < n; i++) dot += v[i - j] * A[k][i];
      var f = 2 * dot / vv;
      for (i = j; i < n; i++) A[k][i] -= f * v[i - j];
    }
    var db = 0;
    for (i = j; i < n; i++) db += v[i - j] * b[i];
    var fb = 2 * db / vv;
    for (i = j; i < n; i++) b[i] -= fb * v[i - j];
    A[j][j] = alpha;
  }
  // R[i][j] = A[j][i] (i ≤ j). R beta = b[0..p-1] 를 뒤에서부터 푼다.
  var beta = new Array(p);
  for (i = p - 1; i >= 0; i--) {
    var s = b[i];
    for (k = i + 1; k < p; k++) s -= A[k][i] * beta[k];
    beta[i] = s / A[i][i];
  }
  // R^-1 (상삼각) 을 구하고 (X'X)^-1 = R^-1 R^-T
  var rinv = [];
  for (i = 0; i < p; i++) rinv.push(new Array(p).fill(0));
  for (j = 0; j < p; j++) {
    rinv[j][j] = 1 / A[j][j];
    for (i = j - 1; i >= 0; i--) {
      var t2 = 0;
      for (k = i + 1; k <= j; k++) t2 += A[k][i] * rinv[k][j];
      rinv[i][j] = -t2 / A[i][i];
    }
  }
  var minv = [];
  for (i = 0; i < p; i++) {
    minv.push(new Array(p));
    for (j = 0; j < p; j++) {
      var m = 0;
      for (k = Math.max(i, j); k < p; k++) m += rinv[i][k] * rinv[j][k];
      minv[i][j] = m;
    }
  }
  return { ok: true, beta: beta, minv: minv };
}

/** 대칭행렬의 고유값(내림차순)과 고유벡터 (순환 Jacobi). vectors[k] 는 k번째 고유벡터. */
function jacobiEigen_(S) {
  var n = S.length, i, j, k;
  var A = S.map(function (r) { return r.slice(); });
  var V = [];
  for (i = 0; i < n; i++) { V.push(new Array(n).fill(0)); V[i][i] = 1; }
  for (var sweep = 0; sweep < 100; sweep++) {
    var off = 0, diag = 0;
    for (i = 0; i < n; i++) {
      diag += A[i][i] * A[i][i];
      for (j = i + 1; j < n; j++) off += A[i][j] * A[i][j];
    }
    if (off <= 1e-30 * diag) break;
    for (var p = 0; p < n - 1; p++) {
      for (var q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-300) continue;
        var theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        var t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        var c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (k = 0; k < n; k++) {
          var akp = A[k][p], akq = A[k][q];
          A[k][p] = c * akp - s * akq;
          A[k][q] = s * akp + c * akq;
        }
        for (k = 0; k < n; k++) {
          var apk = A[p][k], aqk = A[q][k];
          A[p][k] = c * apk - s * aqk;
          A[q][k] = s * apk + c * aqk;
        }
        for (k = 0; k < n; k++) {
          var vkp = V[k][p], vkq = V[k][q];
          V[k][p] = c * vkp - s * vkq;
          V[k][q] = s * vkp + c * vkq;
        }
      }
    }
  }
  var order = [];
  for (i = 0; i < n; i++) order.push(i);
  order.sort(function (a, b) { return A[b][b] - A[a][a]; });
  return {
    values: order.map(function (o) { return A[o][o]; }),
    vectors: order.map(function (o) { return V.map(function (row) { return row[o]; }); })
  };
}

/* ---------- OLS ---------- */

/**
 * OLS 적합. predictors: [{ name, values }], y: 숫자 배열. 절편은 항상 포함한다.
 * 완전 공선성이면 { collinear: 이름 }, 사례가 부족하면 { tooFew: true }.
 */
function olsFit(predictors, y, conf) {
  var d = regDeps_();
  var n = y.length, k = predictors.length, p = k + 1, i, j;
  var res = { n: n, k: k, p: p };
  if (n - p < 1) return { tooFew: true, n: n, p: p };

  var cy = d.centralSums(y);
  var means = [], c2x = [], centered = [];
  for (j = 0; j < k; j++) {
    var cs = d.centralSums(predictors[j].values);
    if (!(cs.c2 > 0)) return { constant: predictors[j].name };
    means.push(cs.mean);
    c2x.push(cs.c2);
    centered.push(predictors[j].values.map(function (x) { return x - cs.mean; }));
  }
  var yc = y.map(function (v) { return v - cy.mean; });
  var fit = qrLeastSquares_(centered, yc);
  if (!fit.ok) return { collinear: predictors[fit.badCol].name };

  var b0 = cy.mean;
  for (j = 0; j < k; j++) b0 -= fit.beta[j] * means[j];
  var B = [b0].concat(fit.beta);

  // 잔차 (원래 척도에서 계산)
  var resid = new Array(n), sse = 0;
  for (i = 0; i < n; i++) {
    var yhat = b0;
    for (j = 0; j < k; j++) yhat += fit.beta[j] * predictors[j].values[i];
    resid[i] = y[i] - yhat;
    sse += resid[i] * resid[i];
  }
  var sst = cy.c2;
  // 완전한 적합: 잔차 제곱합이 실수 오차 수준(상대 1e-24 이하)이면 0으로 본다
  if (sst > 0 && sse <= 1e-24 * sst) sse = 0;
  var ssr = Math.max(0, sst - sse);
  var dfr = n - p, dfm = k;
  res.B = B; res.resid = resid;
  res.sse = sse; res.ssr = ssr; res.sst = sst;
  res.dfm = dfm; res.dfr = dfr; res.dft = n - 1;
  res.msr = ssr / dfm; res.mse = sse / dfr;
  res.s = Math.sqrt(res.mse);
  res.r2 = sst > 0 ? ssr / sst : null;
  res.r = res.r2 === null ? null : Math.sqrt(res.r2);
  res.adjR2 = sst > 0 ? 1 - (1 - res.r2) * (n - 1) / dfr : null;
  res.F = null; res.pF = null;
  if (res.mse > 0) {
    res.F = res.msr / res.mse;
    res.pF = d.pf(res.F, dfm, dfr, false);
  }

  // 계수의 표준오차, t, p, 표준화 계수, 공차, VIF
  var minv = fit.minv;
  // 기울기 계수들의 공분산 행렬 (조절효과의 단순기울기 표준오차에 쓴다). 중심화한 설계에서 기울기의 공분산은 mse · (X'X)^-1
  res.slopeCov = res.mse > 0 ? minv.map(function (row) { return row.map(function (v) { return res.mse * v; }); }) : null;
  res.means = means;
  var tc = d.qt((1 - conf / 100) / 2, dfr, false);
  res.coef = [];
  for (j = -1; j < k; j++) {
    var idx = j + 1, se = null;
    if (res.mse > 0) {
      if (j < 0) {
        var q = 0, a, b;
        for (a = 0; a < k; a++) for (b = 0; b < k; b++) q += means[a] * minv[a][b] * means[b];
        se = res.s * Math.sqrt(1 / n + q);
      } else {
        se = res.s * Math.sqrt(minv[j][j]);
      }
    }
    var t = se !== null ? B[idx] / se : null;
    var c = {
      name: j < 0 ? '(상수)' : predictors[j].name,
      B: B[idx], se: se, t: t,
      p: t !== null ? d.pt(-Math.abs(t), dfr) * 2 : null,
      lo: se !== null ? B[idx] - tc * se : null,
      hi: se !== null ? B[idx] + tc * se : null,
      beta: null, tol: null, vif: null
    };
    if (j >= 0) {
      if (sst > 0) c.beta = B[idx] * Math.sqrt(c2x[j] / sst);
      c.vif = c2x[j] * minv[j][j];
      c.tol = 1 / c.vif;
    }
    res.coef.push(c);
  }

  // 더빈-왓슨 (사례 순서대로)
  var num = 0;
  for (i = 1; i < n; i++) num += (resid[i] - resid[i - 1]) * (resid[i] - resid[i - 1]);
  res.dw = sse > 0 ? num / sse : null;
  return res;
}

/**
 * 공선성 진단 (조건지수·분산 비율 방식): 상수항을 포함한 설계 행렬의 각 열을 길이 1로 맞추고 X'X 의 고유분해를 한다.
 * 조건지수 = sqrt(최대 고유값 / 고유값), 분산 비율 π_ij = (v_ij²/λ_j) / Σ_j (v_ij²/λ_j).
 * 반환: { eig: [고유값], cond: [조건지수], prop: [차원][계수] }
 */
function collinearityDiagnostics(predictors, n) {
  var k = predictors.length, p = k + 1, i, j;
  var cols = [new Array(n).fill(1)];
  predictors.forEach(function (pr) { cols.push(pr.values); });
  var scaled = cols.map(function (c) {
    var s = 0;
    for (var t = 0; t < n; t++) s += c[t] * c[t];
    s = Math.sqrt(s);
    return c.map(function (x) { return x / s; });
  });
  var S = [];
  for (i = 0; i < p; i++) {
    S.push(new Array(p));
    for (j = 0; j < p; j++) {
      var dot = 0;
      for (var t2 = 0; t2 < n; t2++) dot += scaled[i][t2] * scaled[j][t2];
      S[i][j] = dot;
    }
  }
  var e = jacobiEigen_(S);
  var out = { eig: e.values, cond: e.values.map(function (v) { return Math.sqrt(e.values[0] / v); }), prop: [] };
  // phi[i][j] = v_ij² / λ_j (i: 계수, j: 차원)
  var denom = new Array(p).fill(0), phi = [];
  for (i = 0; i < p; i++) {
    phi.push(new Array(p));
    for (j = 0; j < p; j++) {
      phi[i][j] = e.vectors[j][i] * e.vectors[j][i] / e.values[j];
      denom[i] += phi[i][j];
    }
  }
  for (j = 0; j < p; j++) {
    out.prop.push([]);
    for (i = 0; i < p; i++) out.prop[j].push(phi[i][j] / denom[i]);
  }
  return out;
}

/* ---------- 실행 ---------- */

function runRegression(dataset, request) {
  var d = regDeps_();
  var dep = d.getVars(dataset, request.dependent, { numeric: true });
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

  // 목록별 결측 제거
  var all = dep.concat(covs, facs);
  var rowsOk = [];
  for (var i = 0; i < dataset.n; i++) {
    var ok = true;
    for (var a = 0; a < all.length; a++) if (dataset.data[all[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length;
  var pick = function (name) { return rowsOk.map(function (r) { return dataset.data[name][r]; }); };
  var y = pick(dep[0].name);

  var predictors = covs.map(function (v) { return { name: regName_(v), values: pick(v.name) }; });
  var refNotes = [];
  facs.forEach(function (v) {
    var dm = regDummies_(v, pick(v.name), opt.lastRef);
    refNotes.push(dm.note);
    dm.predictors.forEach(function (pr) { predictors.push(pr); });
  });

  var fit = olsFit(predictors, y, conf);
  if (fit.tooFew) throw new Error('사례 수(' + fit.n + ')가 계수의 수(' + fit.p + ', 상수 포함)보다 많아야 합니다.');
  if (fit.constant) throw new Error('변수 "' + fit.constant + '"의 값이 모두 같아 회귀에 쓸 수 없습니다.');
  if (fit.collinear) throw new Error('변수 "' + fit.collinear + '"이(가) 앞선 독립변수들의 선형결합이라 회귀를 계산할 수 없습니다 (완전 공선성). 이 변수를 빼 보세요.');

  var warnings = [];
  if (fit.mse === 0) warnings.push('잔차가 모두 0입니다 (완전한 적합). 표준오차, t, F를 계산할 수 없습니다.');
  fit.coef.forEach(function (c) {
    if (c.vif !== null && c.vif > REG_VIF_WARN_) {
      warnings.push('"' + c.name + '"의 VIF가 ' + c.vif.toFixed(1) + '로 ' + REG_VIF_WARN_ + '을 넘습니다. 독립변수 사이의 공선성이 큽니다.');
    }
  });
  var outliers = fit.s > 0 ? fit.resid.filter(function (e) { return Math.abs(e / fit.s) > 3; }).length : 0;

  var F3 = { type: 'fixed', digits: 3 };
  var blocks = [];

  var summaryCols = [
    { key: 'r', label: 'R', format: F3 },
    { key: 'r2', label: 'R 제곱', format: F3 },
    { key: 'adj', label: '수정된 R 제곱', format: F3 },
    { key: 's', label: '추정값의 표준오차', format: F3 }
  ];
  var summaryRow = { r: fit.r, r2: fit.r2, adj: fit.adjR2, s: fit.s };
  if (opt.dw) { summaryCols.push({ key: 'dw', label: '더빈-왓슨', format: F3 }); summaryRow.dw = fit.dw; }
  blocks.push({
    type: 'table',
    caption: '모형 요약 (종속변수: ' + regName_(dep[0]) + ')',
    columns: summaryCols,
    rows: [summaryRow],
    notes: ['독립변수: 상수 포함 ' + fit.p + '개 계수. 유효 사례 ' + n + '건' + (dataset.n > n ? ' (결측으로 ' + (dataset.n - n) + '건 제외).' : '.')]
  });

  blocks.push({
    type: 'table',
    caption: '분산분석',
    columns: [
      { key: 'source', label: '' },
      { key: 'ss', label: '제곱합', format: F3 },
      { key: 'df', label: '자유도', format: { type: 'int' } },
      { key: 'ms', label: '평균제곱', format: F3 },
      { key: 'F', label: 'F', format: F3 },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
    ],
    rows: [
      { source: '회귀', ss: fit.ssr, df: fit.dfm, ms: fit.msr, F: fit.F, p: fit.pF },
      { source: '잔차', ss: fit.sse, df: fit.dfr, ms: fit.mse },
      { _kind: 'total', source: '합계', ss: fit.sst, df: fit.dft }
    ],
    notes: []
  });

  var coefCols = [
    { key: 'name', label: '' },
    { key: 'B', label: 'B', format: F3 },
    { key: 'se', label: '표준오차', format: F3 },
    { key: 'beta', label: '베타', format: F3 },
    { key: 't', label: 't', format: F3 },
    { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
  ];
  if (opt.ci) {
    coefCols.push({ key: 'lo', label: 'B의 ' + conf + '% 신뢰구간 하한', format: F3 });
    coefCols.push({ key: 'hi', label: 'B의 ' + conf + '% 신뢰구간 상한', format: F3 });
  }
  coefCols.push({ key: 'tol', label: '공차', format: F3 });
  coefCols.push({ key: 'vif', label: 'VIF', format: F3 });
  var coefNotes = ['베타는 표준화 계수이고, 공차 = 1/VIF 이다. 상수에는 베타·공차·VIF가 없다.'].concat(refNotes);
  if (outliers > 0) coefNotes.push('표준화 잔차(잔차/추정값의 표준오차)의 절댓값이 3을 넘는 사례가 ' + outliers + '건 있다.');
  blocks.push({ type: 'table', caption: '계수', columns: coefCols, rows: fit.coef, notes: coefNotes });

  if (opt.collinearity) blocks.push(regCollinearityBlock_(predictors, n, fit.coef, '공선성 진단'));

  return { analysis: 'regression', title: '선형 회귀', blocks: blocks, warnings: warnings };
}

/** 범주형 변수 하나를 더미 열로 바꾼다 (기준 범주 제외). 반환: { predictors: [{ name, values }], note } */
function regDummies_(v, vals, lastRef) {
  var d = regDeps_();
  var cats = [], seenC = {};
  vals.forEach(function (x) {
    var key = typeof x + ':' + x;
    if (!seenC[key]) { seenC[key] = true; cats.push({ value: x }); }
  });
  cats.sort(d.compareValues_);
  cats.forEach(function (c) {
    var label = v.valueLabels[String(c.value)];
    c.label = label !== undefined ? label : String(c.value);
  });
  if (cats.length < 2) throw new Error('범주형 독립변수 "' + regName_(v) + '"의 범주가 ' + cats.length + '개입니다. 2개 이상이어야 합니다.');
  if (cats.length > REG_MAX_CATEGORIES_) {
    throw new Error('범주형 독립변수 "' + regName_(v) + '"의 범주가 ' + cats.length + '개입니다 (최대 ' + REG_MAX_CATEGORIES_ + '개).');
  }
  var refIdx = lastRef ? cats.length - 1 : 0;
  var predictors = [];
  cats.forEach(function (c, idx) {
    if (idx === refIdx) return;
    predictors.push({
      name: regName_(v) + ' = ' + c.label,
      values: vals.map(function (x) { return x === c.value ? 1 : 0; })
    });
  });
  return { predictors: predictors, note: '"' + regName_(v) + '"의 기준 범주: ' + cats[refIdx].label };
}

/** 공선성 진단 표 블록. coefs: 적합 결과의 coef (상수 포함, 열 이름용) */
function regCollinearityBlock_(predictors, n, coefs, caption) {
  var F3 = { type: 'fixed', digits: 3 };
  var diag = collinearityDiagnostics(predictors, n);
  var dcols = [
    { key: 'dim', label: '차원', format: { type: 'int' } },
    { key: 'eig', label: '고유값', format: F3 },
    { key: 'cond', label: '조건지수', format: F3 }
  ];
  coefs.forEach(function (c, idx) { dcols.push({ key: 'v' + idx, label: '분산 비율: ' + c.name, format: F3 }); });
  var drows = diag.eig.map(function (ev, j) {
    var row = { dim: j + 1, eig: ev, cond: diag.cond[j] };
    diag.prop[j].forEach(function (pr, idx) { row['v' + idx] = pr; });
    return row;
  });
  return {
    type: 'table',
    caption: caption,
    columns: dcols,
    rows: drows,
    notes: ['상수항을 포함한 설계 행렬의 각 열을 길이 1로 맞춘 뒤 계산한다. 조건지수가 30을 넘고 한 차원에서 둘 이상의 분산 비율이 0.5를 넘으면 공선성이 의심된다.']
  };
}

if (typeof module !== 'undefined') {
  module.exports = {
    qrLeastSquares_: qrLeastSquares_,
    jacobiEigen_: jacobiEigen_,
    olsFit: olsFit,
    collinearityDiagnostics: collinearityDiagnostics,
    regDummies_: regDummies_,
    regCollinearityBlock_: regCollinearityBlock_,
    regName_: regName_,
    runRegression: runRegression
  };
}
