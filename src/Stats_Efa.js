/**
 * 탐색적 요인분석(EFA). 순수 함수: (dataset, request) → 결과 객체.
 * 입력: 숫자·서열 변수 3~50개 (피어슨 상관행렬로 분석). 사례는 목록별 결측 제거.
 * 추출: 주성분(PCA) 또는 주축분해(PAF, 초기 공통성 = 제곱 다중상관 SMC, 공통성이 수렴할 때까지 반복).
 * 회전: 없음, 베리맥스(직각, Kaiser 정규화, SVD 방식 반복; R 의 stats::varimax 와 같은 알고리즘이지만 끝까지 수렴시킨다),
 *       프로맥스(사각, Kaiser 정규화: 행을 길이 1로 맞춘 뒤 베리맥스, |적재값|^(κ−1) 목표에 최소제곱 회전(κ = 4), 행 길이 복원;
 *       factor_analyzer 처럼 정규화한다 — R 의 stats::promax 는 정규화하지 않는다. 패턴·구조 행렬과 요인 상관).
 * 보조 통계: KMO(전체와 변수별 MSA), Bartlett 구형성 검정, 고유값(스크리 그림), 설명된 총분산, 공통성,
 *           병렬분석(Horn: 같은 크기의 난수 자료의 상관행렬 고유값 평균·95백분위와 비교, 시드 고정),
 *           요인 점수(회귀 방법, Thurstone: Z·R⁻¹·구조행렬)를 새 열로 저장.
 *
 * 표기 규약: 요인(성분)은 회전 뒤 제곱합 큰 순으로 정렬하고 각 열의 적재값 합이 양수가 되도록 부호를 맞춘다.
 * PAF 반복은 1e-6 에 수렴할 때까지 돈다 (반복을 일찍 멈추는 프로그램과는 소수 셋째 자리가 다를 수 있다).
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function efaDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var reg = require('./Stats_Regression.js');
    return {
      getVars: data.getVars, pchisq: dist.pchisq, jacobiEigen_: reg.jacobiEigen_,
      lineChartBlock: require('./Charts.js').lineChartBlock, mulberry32_: require('./Stats_Mediation.js').mulberry32_
    };
  }
  return { getVars: getVars, pchisq: pchisq, jacobiEigen_: jacobiEigen_, lineChartBlock: lineChartBlock, mulberry32_: mulberry32_ };
}

var EFA_MIN_VARS_ = 3;
var EFA_MAX_VARS_ = 50;
var EFA_PAF_TOL_ = 1e-6;
var EFA_PAF_MAX_ITER_ = 300;
var EFA_VARIMAX_MAX_ITER_ = 2000;
var EFA_VARIMAX_EPS_ = 1e-12;
var EFA_CROSS_ = 0.32;          // 교차적재로 보는 적재값 크기
var EFA_SALIENT_ = 0.40;        // 요인에 속한다고 보는 적재값 크기
var EFA_LOW_COMM_ = 0.3;
var EFA_PROMAX_KAPPA_ = 4;
var EFA_PAR_SIMS_MAX_ = 100;     // 병렬분석 모의 횟수 (사례 수 × 변수 수가 크면 줄인다)
var EFA_PAR_SIMS_MIN_ = 20;
var EFA_PAR_BUDGET_ = 3e7;       // 모의 횟수 × 사례 수 × 변수 수² / 2 의 상한
var EFA_PAR_DEFAULT_SEED_ = 2026;
var EFA_PAR_PERCENTILE_ = 0.95;
function efaName_(v) { return v.label ? v.label : v.name; }

/* ---------- 상관행렬과 선형대수 ---------- */

/** 열 단위 자료 cols (p개, 각 길이 n) 의 피어슨 상관행렬. 분산이 0인 열이 있으면 { constant: 열 번호 } */
function efaCorrelation(cols) {
  var p = cols.length, n = cols[0].length, i, j, t;
  var mean = [], c = [];
  for (j = 0; j < p; j++) {
    var s = 0;
    for (t = 0; t < n; t++) s += cols[j][t];
    mean.push(s / n);
  }
  for (j = 0; j < p; j++) {
    var cj = new Array(n), ss = 0;
    for (t = 0; t < n; t++) { cj[t] = cols[j][t] - mean[j]; ss += cj[t] * cj[t]; }
    if (!(ss > 0)) return { constant: j };
    var sd = Math.sqrt(ss);
    for (t = 0; t < n; t++) cj[t] /= sd;
    c.push(cj);
  }
  var R = [];
  for (i = 0; i < p; i++) R.push(new Array(p).fill(0));
  for (i = 0; i < p; i++) {
    R[i][i] = 1;
    for (j = i + 1; j < p; j++) {
      var dot = 0;
      for (t = 0; t < n; t++) dot += c[i][t] * c[j][t];
      R[i][j] = dot; R[j][i] = dot;
    }
  }
  return { R: R };
}

/** 대칭 양정치 행렬의 Cholesky 분해 (하삼각 L). 양정치가 아니면 null */
function efaCholesky_(A) {
  var n = A.length, i, j, k, L = [];
  for (i = 0; i < n; i++) L.push(new Array(n).fill(0));
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      var s = A[i][j];
      for (k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) {
        if (!(s > 1e-12)) return null;
        L[i][i] = Math.sqrt(s);
      } else {
        L[i][j] = s / L[j][j];
      }
    }
  }
  return L;
}

/** Cholesky 인수 L 로 (LL')^-1 과 log det 를 구한다 */
function efaInvertChol_(L) {
  var n = L.length, i, j, k, logdet = 0;
  for (i = 0; i < n; i++) logdet += 2 * Math.log(L[i][i]);
  var inv = [];
  for (i = 0; i < n; i++) inv.push(new Array(n).fill(0));
  for (var col = 0; col < n; col++) {
    var y = new Array(n).fill(0);
    for (i = 0; i < n; i++) {
      var s = i === col ? 1 : 0;
      for (k = 0; k < i; k++) s -= L[i][k] * y[k];
      y[i] = s / L[i][i];
    }
    var x = new Array(n).fill(0);
    for (i = n - 1; i >= 0; i--) {
      var s2 = y[i];
      for (k = i + 1; k < n; k++) s2 -= L[k][i] * x[k];
      x[i] = s2 / L[i][i];
    }
    for (j = 0; j < n; j++) inv[j][col] = x[j];
  }
  return { inv: inv, logdet: logdet };
}

function efaMul_(A, B) {
  var n = A.length, k = B.length, m = B[0].length, out = [], i, j, q;
  for (i = 0; i < n; i++) {
    out.push(new Array(m).fill(0));
    for (q = 0; q < k; q++) {
      var a = A[i][q];
      if (a === 0) continue;
      for (j = 0; j < m; j++) out[i][j] += a * B[q][j];
    }
  }
  return out;
}
function efaT_(A) {
  return A[0].map(function (_, j) { return A.map(function (row) { return row[j]; }); });
}

/* ---------- KMO, Bartlett ---------- */

/** KMO 전체값과 변수별 MSA. 편상관은 역행렬에서 q_ij = -Rinv_ij / sqrt(Rinv_ii Rinv_jj) */
function efaKmo(R, Rinv) {
  var p = R.length, i, j;
  var rs = new Array(p).fill(0), qs = new Array(p).fill(0), rt = 0, qt = 0;
  for (i = 0; i < p; i++) {
    for (j = 0; j < p; j++) {
      if (i === j) continue;
      var r2 = R[i][j] * R[i][j];
      var q = -Rinv[i][j] / Math.sqrt(Rinv[i][i] * Rinv[j][j]);
      rs[i] += r2; qs[i] += q * q; rt += r2; qt += q * q;
    }
  }
  return { kmo: rt / (rt + qt), msa: rs.map(function (r, k) { return r / (r + qs[k]); }) };
}

/** Bartlett 구형성: χ² = −(n − 1 − (2p + 5)/6) ln|R|, 자유도 p(p − 1)/2 */
function efaBartlett(logdet, n, p) {
  var d = efaDeps_();
  var chi2 = -(n - 1 - (2 * p + 5) / 6) * logdet;
  var df = p * (p - 1) / 2;
  return { chi2: chi2, df: df, p: d.pchisq(chi2, df, false) };
}

/* ---------- 추출 ---------- */

/** 고유값 분해 결과에서 열 m 개의 적재값 행렬 (p × m): 고유벡터 × √고유값 (음수 고유값은 0) */
function efaLoadingsFromEigen_(e, m) {
  var p = e.vectors[0].length, L = [], i, k;
  for (i = 0; i < p; i++) {
    L.push([]);
    for (k = 0; k < m; k++) L[i].push(e.vectors[k][i] * Math.sqrt(Math.max(e.values[k], 0)));
  }
  return L;
}

function efaRowSS_(L) {
  return L.map(function (row) { var s = 0; row.forEach(function (v) { s += v * v; }); return s; });
}
function efaColSS_(L) {
  var m = L[0].length, out = new Array(m).fill(0);
  L.forEach(function (row) { row.forEach(function (v, k) { out[k] += v * v; }); });
  return out;
}

/** 주축분해: 대각을 공통성으로 바꾼 상관행렬의 고유분해를 공통성이 수렴할 때까지 반복 */
function efaPaf(R, Rinv, m) {
  var d = efaDeps_();
  var p = R.length, i, j;
  var h2 = [];
  for (i = 0; i < p; i++) h2.push(1 - 1 / Rinv[i][i]);
  var init = h2.slice();
  var L = null, iter = 0, converged = false;
  while (iter < EFA_PAF_MAX_ITER_) {
    iter++;
    var Rr = R.map(function (row) { return row.slice(); });
    for (i = 0; i < p; i++) Rr[i][i] = h2[i];
    var e = d.jacobiEigen_(Rr);
    L = efaLoadingsFromEigen_(e, m);
    var nh = efaRowSS_(L), change = 0;
    for (i = 0; i < p; i++) change = Math.max(change, Math.abs(nh[i] - h2[i]));
    h2 = nh;
    if (change < EFA_PAF_TOL_) { converged = true; break; }
  }
  return { L: L, h2: h2, initial: init, iterations: iter, converged: converged };
}

/* ---------- 회전 ---------- */

/**
 * 베리맥스 (Kaiser 정규화). 행렬 B = x'(z³ − z·diag(열 제곱합)/p) 의 극분해 인수 B(B'B)^(-1/2) 를 회전행렬로 쓴다
 * (R 의 stats::varimax 가 SVD 로 U V' 를 구하는 것과 같다).
 */
function efaVarimax(L) {
  var d = efaDeps_();
  var p = L.length, m = L[0].length, i, j, k, a, b;
  var sc = L.map(function (row) { var s = 0; row.forEach(function (v) { s += v * v; }); return Math.sqrt(s) || 1; });
  var x = L.map(function (row, ri) { return row.map(function (v) { return v / sc[ri]; }); });
  var T = [];
  for (i = 0; i < m; i++) { T.push(new Array(m).fill(0)); T[i][i] = 1; }
  var dsum = 0;
  for (var it = 0; it < EFA_VARIMAX_MAX_ITER_; it++) {
    var z = x.map(function (row) {
      return T[0].map(function (_, c) { var s = 0; for (var q = 0; q < m; q++) s += row[q] * T[q][c]; return s; });
    });
    var cs = new Array(m).fill(0);
    z.forEach(function (row) { row.forEach(function (v, c) { cs[c] += v * v; }); });
    var G = z.map(function (row) { return row.map(function (v, c) { return v * v * v - v * cs[c] / p; }); });
    var B = [];
    for (a = 0; a < m; a++) {
      B.push(new Array(m).fill(0));
      for (b = 0; b < m; b++) { var s = 0; for (k = 0; k < p; k++) s += x[k][a] * G[k][b]; B[a][b] = s; }
    }
    var BtB = [];
    for (a = 0; a < m; a++) {
      BtB.push(new Array(m).fill(0));
      for (b = 0; b < m; b++) { var s2 = 0; for (k = 0; k < m; k++) s2 += B[k][a] * B[k][b]; BtB[a][b] = s2; }
    }
    var e = d.jacobiEigen_(BtB);
    var sv = e.values.map(function (v) { return Math.sqrt(Math.max(v, 0)); });
    var W = [];
    for (a = 0; a < m; a++) {
      W.push(new Array(m).fill(0));
      for (b = 0; b < m; b++) {
        var s3 = 0;
        for (k = 0; k < m; k++) if (sv[k] > 1e-14) s3 += e.vectors[k][a] * e.vectors[k][b] / sv[k];
        W[a][b] = s3;
      }
    }
    for (a = 0; a < m; a++) {
      for (b = 0; b < m; b++) { var s4 = 0; for (k = 0; k < m; k++) s4 += B[a][k] * W[k][b]; T[a][b] = s4; }
    }
    var dpast = dsum;
    dsum = 0;
    sv.forEach(function (v) { dsum += v; });
    if (dsum < dpast * (1 + EFA_VARIMAX_EPS_)) break;
  }
  var rotated = x.map(function (row, ri) {
    return T[0].map(function (_, c) { var s = 0; for (var q = 0; q < m; q++) s += row[q] * T[q][c]; return s * sc[ri]; });
  });
  return { L: rotated, T: T };
}

/** 열 정렬 순서(제곱합 큰 순)와 부호(각 열의 합이 양수가 되도록) */
function efaCanonInfo_(L) {
  var m = L[0].length, ss = efaColSS_(L);
  var order = [];
  for (var k = 0; k < m; k++) order.push(k);
  order.sort(function (a, b) { return ss[b] - ss[a]; });
  var sign = order.map(function (k) { var s = 0; L.forEach(function (row) { s += row[k]; }); return s < 0 ? -1 : 1; });
  return { order: order, sign: sign };
}

/** 열을 제곱합 큰 순으로 정렬하고 각 열의 합이 양수가 되도록 부호를 맞춘다 */
function efaCanonical(L) {
  var info = efaCanonInfo_(L);
  return L.map(function (row) { return info.order.map(function (k, q) { return row[k] * info.sign[q]; }); });
}

/**
 * 프로맥스 (Kaiser 정규화, κ = 4): 행을 길이 1로 맞춘 적재값의 베리맥스 결과 x 로 목표 Q = x·|x|^(κ−1) 를 만들고
 * U = (x'x)⁻¹ x'Q 를 열마다 (U'U)⁻¹ 의 대각이 1이 되도록 정규화한다. 패턴 행렬 P = (x U)·행 길이(공통성의 제곱근), 요인 상관 Φ = (U'U)⁻¹.
 * 정준화(열 순서, 부호)는 패턴 행렬의 열 제곱합 기준이고 Φ 에도 같이 적용한다. 계산 불가면 null
 */
function efaPromax(L, kappa) {
  var m = L[0].length, j, k;
  var rl = L.map(function (row) { var s = 0; row.forEach(function (v) { s += v * v; }); return Math.sqrt(s) || 1; });
  var vm = efaVarimax(L.map(function (row, i) { return row.map(function (v) { return v / rl[i]; }); })).L;
  var Q = vm.map(function (row) { return row.map(function (v) { return v * Math.pow(Math.abs(v), kappa - 1); }); });
  var xt = efaT_(vm);
  var ch = efaCholesky_(efaMul_(xt, vm));
  if (!ch) return null;
  var U = efaMul_(efaInvertChol_(ch).inv, efaMul_(xt, Q));
  var ch2 = efaCholesky_(efaMul_(efaT_(U), U));
  if (!ch2) return null;
  var pinv = efaInvertChol_(ch2).inv;
  for (j = 0; j < m; j++) {
    var sd = Math.sqrt(pinv[j][j]);
    for (k = 0; k < m; k++) U[k][j] *= sd;
  }
  var ch3 = efaCholesky_(efaMul_(efaT_(U), U));
  if (!ch3) return null;
  var phi = efaInvertChol_(ch3).inv;
  var P = efaMul_(vm, U).map(function (row, i) { return row.map(function (v) { return v * rl[i]; }); });
  var info = efaCanonInfo_(P);
  var Pc = P.map(function (row) { return info.order.map(function (c, q) { return row[c] * info.sign[q]; }); });
  var phiC = info.order.map(function (a, qa) {
    return info.order.map(function (b, qb) { return phi[a][b] * info.sign[qa] * info.sign[qb]; });
  });
  return { P: Pc, phi: phiC };
}

/** 구조 행렬 S = P Φ (패턴 행렬과 요인 상관) */
function efaStructure(P, phi) { return efaMul_(P, phi); }

/** 요인 점수 가중치 (회귀 방법): W = R⁻¹ S. 점수 = 표준화한 변수 × W */
function efaScoreWeights(Rinv, S) { return efaMul_(Rinv, S); }

/** 표준정규 난수 한 쌍씩 만드는 Box-Muller */
function efaGauss_(rnd) {
  var spare = null;
  return function () {
    if (spare !== null) { var s = spare; spare = null; return s; }
    var u = 0, v = 0;
    while (u === 0) u = rnd();
    v = rnd();
    var r = Math.sqrt(-2 * Math.log(u));
    spare = r * Math.sin(2 * Math.PI * v);
    return r * Math.cos(2 * Math.PI * v);
  };
}

/** 정렬된 배열의 분위수 (선형 보간, numpy 기본과 같음) */
function efaQuantile_(sorted, q) {
  var pos = (sorted.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
}

/**
 * 병렬분석(Horn): n × p 표준정규 난수 자료의 상관행렬 고유값을 sims 번 구해 순번별 평균과 95백분위를 낸다.
 * 반환: { sims, mean: [p], p95: [p] }. 시드가 같으면 항상 같은 결과.
 */
function efaParallel(n, p, seed, sims) {
  var d = efaDeps_();
  var gauss = efaGauss_(d.mulberry32_(seed));
  var all = [], s, j, t;
  for (s = 0; s < sims; s++) {
    var cols = [];
    for (j = 0; j < p; j++) {
      var c = new Array(n);
      for (t = 0; t < n; t++) c[t] = gauss();
      cols.push(c);
    }
    var cr = efaCorrelation(cols);
    if (cr.constant !== undefined) { s--; continue; }
    all.push(d.jacobiEigen_(cr.R).values);
  }
  var mean = [], p95 = [];
  for (j = 0; j < p; j++) {
    var col = all.map(function (e) { return e[j]; });
    var sum = 0;
    col.forEach(function (v) { sum += v; });
    mean.push(sum / sims);
    p95.push(efaQuantile_(col.slice().sort(function (a, b) { return a - b; }), EFA_PAR_PERCENTILE_));
  }
  return { sims: sims, mean: mean, p95: p95 };
}

/** 병렬분석 모의 횟수: 크기가 크면 줄인다 (최소 20회) */
function efaParallelSims(n, p) {
  var per = n * p * p / 2;
  return Math.max(EFA_PAR_SIMS_MIN_, Math.min(EFA_PAR_SIMS_MAX_, Math.floor(EFA_PAR_BUDGET_ / per)));
}

/** 크기순 정렬 순서: 절댓값이 가장 큰 요인별로 묶고(요인 번호 순), 묶음 안에서는 적재값 큰 순 */
function efaSortOrder(L) {
  var rows = L.map(function (row, i) {
    var best = 0;
    row.forEach(function (v, k) { if (Math.abs(v) > Math.abs(row[best])) best = k; });
    return { i: i, f: best, v: Math.abs(row[best]) };
  });
  rows.sort(function (a, b) { return a.f - b.f || b.v - a.v || a.i - b.i; });
  return rows.map(function (r) { return r.i; });
}

/* ---------- 분석 본체 (표를 만들지 않는 계산) ---------- */

/**
 * columns: 변수별 숫자 배열 (결측 없음). opts: { method: 'paf'|'pca', rotation: 'varimax'|'promax'|'none', nfactors: 정수 또는 null(고유값 1 이상) }
 * 반환: { n, p, R, Rinv, eig, kmo, msa, bartlett, m, mAuto, L(회전 전), h2, h2Initial, rotated(회전 뒤 또는 null: 프로맥스는 패턴 행렬),
 *        oblique, phi, structure(프로맥스), ssInit, ssExt, ssRot(직각 회전만), paf }
 * 특이 상관행렬이면 { singular: true }
 */
function efaAnalyze(columns, opts) {
  var d = efaDeps_();
  var p = columns.length, n = columns[0].length;
  var cr = efaCorrelation(columns);
  if (cr.constant !== undefined) return { constant: cr.constant };
  var R = cr.R;
  var chol = efaCholesky_(R);
  if (!chol) return { singular: true };
  var ic = efaInvertChol_(chol);
  var e = d.jacobiEigen_(R);
  var kk = efaKmo(R, ic.inv);
  var out = { n: n, p: p, R: R, Rinv: ic.inv, eig: e.values, kmo: kk.kmo, msa: kk.msa, bartlett: efaBartlett(ic.logdet, n, p) };

  var kaiser = 0;
  e.values.forEach(function (v) { if (v > 1 + 1e-9) kaiser++; });
  out.kaiser = kaiser;
  var m = opts.nfactors !== null && opts.nfactors !== undefined ? opts.nfactors : kaiser;
  out.mAuto = opts.nfactors === null || opts.nfactors === undefined;
  out.m = m;
  if (m < 1) return Object.assign(out, { noFactor: true });

  var L, h2, initial;
  if (opts.method === 'pca') {
    L = efaLoadingsFromEigen_(e, m);
    h2 = efaRowSS_(L);
    initial = new Array(p).fill(1);
  } else {
    var paf = efaPaf(R, ic.inv, m);
    L = paf.L; h2 = paf.h2; initial = paf.initial;
    out.paf = { iterations: paf.iterations, converged: paf.converged };
  }
  L = efaCanonical(L);
  out.L = L;
  out.h2 = h2;
  out.h2Initial = initial;
  out.ssInit = e.values.slice();
  out.ssExt = efaColSS_(L);
  out.rotated = null;
  out.oblique = false;
  if (opts.rotation === 'varimax' && m > 1) {
    out.rotated = efaCanonical(efaVarimax(L).L);
    out.ssRot = efaColSS_(out.rotated);
  } else if (opts.rotation === 'promax' && m > 1) {
    var pm = efaPromax(L, EFA_PROMAX_KAPPA_);
    if (!pm) return Object.assign(out, { promaxFailed: true });
    out.oblique = true;
    out.rotated = pm.P;
    out.phi = pm.phi;
    out.structure = efaStructure(pm.P, pm.phi);
  }
  return out;
}

/* ---------- 실행 ---------- */

function runEfa(dataset, request) {
  var d = efaDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < EFA_MIN_VARS_) throw new Error('요인분석에는 변수를 ' + EFA_MIN_VARS_ + '개 이상 선택하세요.');
  if (vars.length > EFA_MAX_VARS_) throw new Error('변수는 ' + EFA_MAX_VARS_ + '개까지 선택할 수 있습니다.');
  var opt = request.options || {};
  var method = opt.method === 'pca' ? 'pca' : 'paf';
  var rotation = opt.rotation === 'none' ? 'none' : opt.rotation === 'promax' ? 'promax' : 'varimax';
  var seed = opt.seed === undefined || opt.seed === null || String(opt.seed).trim() === '' ? EFA_PAR_DEFAULT_SEED_ : Number(opt.seed);
  if (!isFinite(seed) || seed !== Math.floor(seed) || Math.abs(seed) > 2147483647) throw new Error('"병렬분석 시드"는 정수여야 합니다.');
  var p = vars.length;
  var nfRaw = opt.nfactors;
  var nfactors = null;
  if (nfRaw !== undefined && nfRaw !== null && String(nfRaw).trim() !== '') {
    nfactors = Number(nfRaw);
    var maxF = method === 'pca' ? p : p - 1;
    if (!isFinite(nfactors) || nfactors !== Math.floor(nfactors) || nfactors < 1 || nfactors > maxF) {
      throw new Error('"요인 수"는 1 이상 ' + maxF + ' 이하의 정수여야 합니다 (비우면 고유값 1 이상으로 정합니다).');
    }
  }
  var suppress = opt.suppress === undefined || opt.suppress === null || String(opt.suppress).trim() === '' ? 0.3 : Number(opt.suppress);
  if (!isFinite(suppress) || suppress < 0 || suppress >= 1) throw new Error('"작은 적재값 숨기기 기준"은(는) 0 이상 1 미만이어야 합니다.');

  // 목록별 결측 제거
  var rowsOk = [];
  for (var i = 0; i < dataset.n; i++) {
    var ok = true;
    for (var a = 0; a < p; a++) if (dataset.data[vars[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length;
  if (n <= p) throw new Error('모든 변수에 값이 있는 사례가 ' + n + '건이라 계산할 수 없습니다 (변수 수 ' + p + '개보다 많아야 합니다). 결측이 많은 변수를 빼 보세요.');
  var columns = vars.map(function (v) { return rowsOk.map(function (r) { return dataset.data[v.name][r]; }); });
  var names = vars.map(efaName_);

  var res = efaAnalyze(columns, { method: method, rotation: rotation, nfactors: nfactors });
  if (res.constant !== undefined) throw new Error('변수 "' + names[res.constant] + '"의 값이 모두 같아 요인분석에 쓸 수 없습니다.');
  if (res.singular) throw new Error('상관행렬이 특이해 계산할 수 없습니다. 거의 같은 값을 갖는 변수(예: 합계 변수와 그 구성 문항)가 함께 들어 있는지 확인하고 하나를 빼 보세요.');
  if (res.promaxFailed) throw new Error('프로맥스 회전을 계산하지 못했습니다. 요인 수를 줄이거나 회전을 베리맥스로 바꿔 보세요.');
  if (res.noFactor) throw new Error('고유값이 1 이상인 요인이 없습니다. 변수 사이의 상관이 매우 낮습니다. "요인 수"를 직접 지정하거나 변수를 바꿔 보세요.');
  var m = res.m;
  var word = method === 'pca' ? '성분' : '요인';

  var warnings = [];
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (선택한 변수 중 하나라도 빠진 사례). 사용한 사례는 ' + n + '건입니다.');
  if (n < 100 || n < 5 * p) warnings.push('사례 수(' + n + '건)가 적습니다. 변수 수의 5배 이상이면서 100건 이상이 흔한 권장입니다. 요인 구조가 불안정할 수 있습니다.');
  if (res.kmo < 0.5) warnings.push('KMO가 ' + res.kmo.toFixed(3) + '로 0.5 미만이라 요인분석에 적합하지 않을 수 있습니다. 공통성 표의 변수별 MSA가 낮은 변수를 빼 보세요.');
  if (res.bartlett.p >= 0.05) warnings.push('Bartlett 구형성 검정이 유의하지 않습니다 (유의확률 ≥ .05). 변수 사이에 상관이 거의 없어 요인으로 묶기 어려울 수 있습니다.');
  if (res.paf && !res.paf.converged) warnings.push('주축분해의 반복이 ' + res.paf.iterations + '회까지 수렴하지 않았습니다. 요인 수를 줄이거나 변수를 바꿔 보세요.');
  var heywood = [];
  res.h2.forEach(function (h, k) { if (h >= 1 - 1e-9 && method === 'paf') heywood.push(names[k]); });
  if (heywood.length) warnings.push('공통성이 1 이상인 변수가 있습니다 (Heywood 사례): ' + heywood.join(', ') + '. 요인 수가 너무 많거나 사례가 적을 때 생기니 요인 수를 줄여 보세요.');
  var lowComm = [];
  res.h2.forEach(function (h, k) { if (h < EFA_LOW_COMM_) lowComm.push(names[k]); });
  if (lowComm.length) warnings.push('추출 뒤 공통성이 ' + EFA_LOW_COMM_ + ' 미만인 변수: ' + lowComm.join(', ') + '. 공통요인으로 잘 설명되지 않으니 제외를 고려해 보세요.');

  var finalL = res.rotated || res.L;                   // 해석에 쓰는 적재 행렬
  var judged = res.rotated || (m === 1 ? res.L : null);
  if (judged) {
    var cross = [], cnt = new Array(m).fill(0);
    judged.forEach(function (row, k) {
      var big = 0;
      row.forEach(function (v, f) { if (Math.abs(v) >= EFA_CROSS_) big++; });
      if (big >= 2) cross.push(names[k]);
      var best = 0;
      row.forEach(function (v, f) { if (Math.abs(v) > Math.abs(row[best])) best = f; });
      if (Math.abs(row[best]) >= EFA_SALIENT_) cnt[best]++;
    });
    if (cross.length) warnings.push('둘 이상의 ' + word + '에 ' + EFA_CROSS_ + ' 이상으로 걸치는 변수(교차적재): ' + cross.join(', ') + '.');
    var thin = [];
    cnt.forEach(function (c, f) { if (c < 2) thin.push((f + 1)); });
    if (thin.length) warnings.push(word + ' ' + thin.join(', ') + '에 속하는 변수(적재값 ' + EFA_SALIENT_.toFixed(2) + ' 이상)가 2개 미만입니다. 요인 수를 줄여 보세요.');
  }
  if (m === 1 && rotation !== 'none') warnings.push(word + '이 1개라 회전하지 않았습니다.');

  var F3 = { type: 'fixed', digits: 3 };
  var blocks = [];

  var obl = res.oblique;
  var rotText = !res.rotated ? '회전 없음' : obl ? '프로맥스 (사각, Kaiser 정규화, κ=' + EFA_PROMAX_KAPPA_ + ')' : '베리맥스 (Kaiser 정규화)';
  var kmoWord = res.kmo >= 0.9 ? '매우 좋음' : res.kmo >= 0.8 ? '좋음' : res.kmo >= 0.7 ? '보통' : res.kmo >= 0.6 ? '평범' : res.kmo >= 0.5 ? '나쁨' : '부적합';
  blocks.push({
    type: 'table',
    caption: 'KMO와 Bartlett 검정',
    columns: [
      { key: 'kmo', label: 'KMO 표본적절성', format: F3 },
      { key: 'chi2', label: 'Bartlett 근사 카이제곱', format: F3 },
      { key: 'df', label: '자유도', format: { type: 'int' } },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
    ],
    rows: [{ kmo: res.kmo, chi2: res.bartlett.chi2, df: res.bartlett.df, p: res.bartlett.p }],
    notes: [
      '분석 설정: 추출 ' + (method === 'pca' ? '주성분' : '주축분해') + ', ' + rotText + ', ' + word + ' 수 ' + m + (res.mAuto ? ' (고유값 1 이상)' : ' (직접 지정)') + ', 변수 ' + p + '개, 사례 ' + n + '건.',
      'KMO 해석(Kaiser): 0.9 이상 매우 좋음, 0.8 좋음, 0.7 보통, 0.6 평범, 0.5 나쁨, 0.5 미만 부적합 → 이 자료는 "' + kmoWord + '". Bartlett 검정이 유의하면 변수 사이에 상관이 있어 요인분석을 할 수 있다.'
    ]
  });

  blocks.push({
    type: 'table',
    caption: '공통성',
    columns: [
      { key: 'name', label: '' },
      { key: 'init', label: '초기', format: F3 },
      { key: 'ext', label: '추출', format: F3 },
      { key: 'msa', label: '표본적절성(MSA)', format: F3 }
    ],
    rows: names.map(function (nm, k) { return { name: nm, init: res.h2Initial[k], ext: res.h2[k], msa: res.msa[k] }; }),
    notes: [(method === 'pca'
      ? '주성분은 초기 공통성이 1이고, 추출 공통성은 뽑은 ' + word + '이 변수의 분산을 설명하는 비율이다.'
      : '주축분해의 초기 공통성은 다른 변수들로 그 변수를 예측한 R²(제곱 다중상관)이고, 추출 공통성은 뽑은 ' + word + '이 설명하는 비율이다.') +
      ' MSA는 변수별 KMO로 0.5 미만이면 그 변수가 요인분석에 맞지 않을 수 있다.']
  });

  var tvCols = [
    { key: 'comp', label: word, format: { type: 'int' } },
    { key: 'iTot', label: '초기 고유값: 합계', format: F3 },
    { key: 'iPct', label: '초기 고유값: 분산 %', format: F3 },
    { key: 'iCum', label: '초기 고유값: 누적 %', format: F3 },
    { key: 'eTot', label: '추출 제곱합 적재값: 합계', format: F3 },
    { key: 'ePct', label: '추출 제곱합 적재값: 분산 %', format: F3 },
    { key: 'eCum', label: '추출 제곱합 적재값: 누적 %', format: F3 }
  ];
  if (res.rotated && !obl) {
    tvCols.push({ key: 'rTot', label: '회전 제곱합 적재값: 합계', format: F3 });
    tvCols.push({ key: 'rPct', label: '회전 제곱합 적재값: 분산 %', format: F3 });
    tvCols.push({ key: 'rCum', label: '회전 제곱합 적재값: 누적 %', format: F3 });
  }
  var iCum = 0, eCum = 0, rCum = 0;
  var tvRows = res.ssInit.map(function (ev, k) {
    iCum += 100 * ev / p;
    var row = { comp: k + 1, iTot: ev, iPct: 100 * ev / p, iCum: iCum };
    if (k < m) {
      eCum += 100 * res.ssExt[k] / p;
      row.eTot = res.ssExt[k]; row.ePct = 100 * res.ssExt[k] / p; row.eCum = eCum;
      if (res.rotated && !obl) { rCum += 100 * res.ssRot[k] / p; row.rTot = res.ssRot[k]; row.rPct = 100 * res.ssRot[k] / p; row.rCum = rCum; }
    }
    return row;
  });
  var tvNotes = ['분산 %는 합계를 변수 수(' + p + ')로 나눈 값이다 (변수를 표준화했으므로 전체 분산 = 변수 수).'];
  if (res.mAuto) tvNotes.push('고유값 1 이상 기준은 ' + word + ' 수를 많게 잡는 경향이 있다. 스크리 그림에서 꺾이는 지점과 해석 가능성도 함께 보고 정하세요.');
  if (res.rotated && !obl) tvNotes.push('직각 회전(베리맥스)에서는 회전하면 ' + word + '별 합계는 달라지지만 합은 같다.');
  if (obl) tvNotes.push('사각 회전(프로맥스)에서는 ' + word + '이 서로 상관되어 회전 뒤 제곱합을 더해 분산 비율을 낼 수 없다. 회전 뒤 합계는 표시하지 않는다.');
  blocks.push({ type: 'table', caption: '설명된 총분산', columns: tvCols, rows: tvRows, notes: tvNotes });

  var par = null;
  if (opt.parallel) {
    var sims = efaParallelSims(n, p);
    par = efaParallel(n, p, seed, sims);
    var rec = 0;
    while (rec < p && res.ssInit[rec] > par.p95[rec]) rec++;
    par.rec = rec;
    var shown = Math.min(p, Math.max(m, rec, res.kaiser) + 3);
    var parRows = [];
    for (var q = 0; q < shown; q++) {
      parRows.push({ comp: q + 1, obs: res.ssInit[q], mean: par.mean[q], p95: par.p95[q], judge: res.ssInit[q] > par.p95[q] && q < rec ? '유지' : '제외' });
    }
    blocks.push({
      type: 'table', caption: '병렬분석',
      columns: [
        { key: 'comp', label: word, format: { type: 'int' } },
        { key: 'obs', label: '실제 자료의 고유값', format: F3 },
        { key: 'mean', label: '난수 자료 고유값의 평균', format: F3 },
        { key: 'p95', label: '난수 자료 고유값의 95백분위', format: F3 },
        { key: 'judge', label: '판정' }
      ],
      rows: parRows,
      notes: [
        '병렬분석(Horn)은 같은 크기(사례 ' + n + '건, 변수 ' + p + '개)의 난수 자료에서 상관행렬의 고유값을 ' + sims + '번 구해 평균과 95백분위를 낸 뒤, 실제 고유값이 난수의 95백분위보다 큰 ' + word + '만 유지하는 방법이다 (시드 ' + seed + ', 시드가 같으면 같은 결과). 고유값 1 이상 기준보다 ' + word + ' 수를 덜 부풀린다고 알려져 있다.',
        '병렬분석이 권하는 ' + word + ' 수: ' + rec + '개 (고유값 1 이상 기준: ' + res.kaiser + '개, 이 분석에 쓴 수: ' + m + '개). 표에는 앞의 ' + shown + '개만 보인다.'
      ]
    });
    if (res.mAuto && rec !== res.kaiser && rec >= 1) {
      warnings.push('병렬분석은 ' + word + ' ' + rec + '개를 권하지만 고유값 1 이상 기준으로는 ' + res.kaiser + '개가 나와 ' + m + '개로 분석했습니다. "요인 수" 칸에 ' + rec + '을 넣어 다시 분석해 비교해 보세요.');
    }
  }

  if (opt.scree !== false) {
    var pts = res.ssInit.map(function (ev, k) { return [k + 1, ev]; });
    var series = [{ label: '고유값', points: pts }, { label: '기준 1', points: [[1, 1], [p, 1]] }];
    if (par) series.push({ label: '난수 자료의 95백분위', points: par.p95.map(function (v, k) { return [k + 1, v]; }) });
    blocks.push(d.lineChartBlock('스크리 그림', word + ' 번호', '고유값', series));
  }

  var sortRows = opt.sort === false ? names.map(function (_, k) { return k; }) : efaSortOrder(finalL);
  var loadingCols = function (L, caption, note) {
    var cols = [{ key: 'name', label: '' }];
    for (var f = 0; f < m; f++) cols.push({ key: 'f' + f, label: word + ' ' + (f + 1), format: F3 });
    var rows = sortRows.map(function (k) {
      var row = { name: names[k] };
      L[k].forEach(function (v, f) { row['f' + f] = Math.abs(v) < suppress ? null : v; });
      return row;
    });
    var notes = [note];
    if (suppress > 0) notes.push('절댓값이 ' + suppress + ' 미만인 적재값은 비워 보여 준다 (세부 옵션 "작은 적재값 숨기기 기준").');
    if (opt.sort !== false) notes.push('변수는 가장 크게 적재된 ' + word + '별로 묶고 적재값 큰 순으로 정렬했다.');
    return { type: 'table', caption: caption, columns: cols, rows: rows, notes: notes };
  };
  blocks.push(loadingCols(res.L, (res.rotated ? word + ' 행렬 (회전 전)' : word + ' 행렬'),
    '추출 방법: ' + (method === 'pca' ? '주성분' : '주축분해') + (res.paf ? ' (' + res.paf.iterations + '회 반복 후 수렴)' : '') + '. 각 ' + word + '의 적재값 합이 양수가 되도록 부호를 맞췄다.'));
  if (res.rotated && !obl) {
    blocks.push(loadingCols(res.rotated, '회전된 ' + word + ' 행렬', '회전 방법: 베리맥스 (Kaiser 정규화). ' + word + '은 회전 뒤 제곱합이 큰 순으로 번호를 붙였다. 해석에는 이 행렬을 쓴다.'));
  }
  if (obl) {
    blocks.push(loadingCols(res.rotated, '패턴 행렬', '회전 방법: 프로맥스 (사각, Kaiser 정규화, κ=' + EFA_PROMAX_KAPPA_ + '). 패턴 행렬은 다른 ' + word + '의 영향을 뺀 각 ' + word + '의 고유한 기여라 ' + word + ' 해석에 쓴다. ' + word + '은 회전 전 적재값의 제곱합이 큰 순으로 번호를 붙였다.'));
    blocks.push(loadingCols(res.structure, '구조 행렬', '구조 행렬은 변수와 ' + word + ' 사이의 상관 (= 패턴 행렬 × ' + word + ' 상관 행렬)이다. 사각 회전에서는 패턴 행렬과 다르다.'));
    var phiCols = [{ key: 'name', label: '' }];
    for (var pf = 0; pf < m; pf++) phiCols.push({ key: 'f' + pf, label: word + ' ' + (pf + 1), format: F3 });
    blocks.push({
      type: 'table', caption: word + ' 상관 행렬', columns: phiCols,
      rows: res.phi.map(function (row, a) {
        var r = { name: word + ' ' + (a + 1) };
        row.forEach(function (v, b) { r['f' + b] = v; });
        return r;
      }),
      notes: [word + ' 사이의 상관이 모두 0.3보다 작으면 직각 회전(베리맥스)과 결과가 비슷하다. 상관이 크면 사각 회전이 알맞다.']
    });
  }

  var result = { analysis: 'efa', title: '탐색적 요인분석', blocks: blocks, warnings: warnings };
  if (opt.save) {
    var S = obl ? res.structure : finalL;
    var W = efaScoreWeights(res.Rinv, S);
    var stats = columns.map(function (c) {
      var mean = 0, ss = 0;
      c.forEach(function (v) { mean += v; });
      mean /= n;
      c.forEach(function (v) { ss += (v - mean) * (v - mean); });
      return { mean: mean, sd: Math.sqrt(ss / (n - 1)) };
    });
    var full = [];
    for (var f2 = 0; f2 < m; f2++) full.push(new Array(dataset.n).fill(null));
    rowsOk.forEach(function (r, idx) {
      for (var f3 = 0; f3 < m; f3++) {
        var s = 0;
        for (var j = 0; j < p; j++) s += (columns[j][idx] - stats[j].mean) / stats[j].sd * W[j][f3];
        full[f3][r] = s;
      }
    });
    var tag = (method === 'pca' ? '주성분' : '주축분해') + '·' + (obl ? '프로맥스' : res.rotated ? '베리맥스' : '회전 없음');
    result.saveColumns = full.map(function (vals, f4) {
      return { name: word + '점수' + (f4 + 1), label: word + ' 점수 ' + (f4 + 1) + ' (회귀 방법, ' + tag + ')', scale: '연속', values: vals, fixed: 4 };
    });
    result.warnings.push('요인 점수는 회귀 방법(Thurstone)으로 구했습니다. 표준화한 변수에 R⁻¹·' + (obl ? '구조 행렬' : '적재 행렬') + '을 곱한 값이라 평균 0이지만 분산은 1이 아닐 수 있습니다. 보고서에는 요인별 문항의 평균·합계 점수를 쓰는 것이 더 흔합니다.');
  }
  return result;
}

if (typeof module !== 'undefined') {
  module.exports = {
    efaCorrelation: efaCorrelation,
    efaKmo: efaKmo,
    efaBartlett: efaBartlett,
    efaPaf: efaPaf,
    efaVarimax: efaVarimax,
    efaCanonical: efaCanonical,
    efaPromax: efaPromax,
    efaStructure: efaStructure,
    efaScoreWeights: efaScoreWeights,
    efaParallel: efaParallel,
    efaParallelSims: efaParallelSims,
    efaSortOrder: efaSortOrder,
    efaAnalyze: efaAnalyze,
    runEfa: runEfa
  };
}
