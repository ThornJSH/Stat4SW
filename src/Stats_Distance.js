/**
 * 거리·유사도(근접도 행렬)와 마할라노비스 이상치. 순수 함수: (dataset, request) → 결과 객체.
 *
 * - 사례 간 거리(기본): 사례(행)마다 선택한 변수 값의 벡터. 변수 간 거리: 변수(열)마다 사례 값의 벡터.
 * - 측도: 거리(유클리드, 제곱 유클리드, 맨해튼, 체비쇼프, 민코프스키, 마할라노비스, 코사인 거리, 상관 거리)와
 *   유사도(코사인, 피어슨 상관), 이분(0/1) 자료의 단순일치·Jaccard·Dice(유사도)와 Hamming(거리)
 * - 표준화는 변수(열)별로 한다: z점수, 0~1, −1~1, 최댓값 1, 평균 1, 표준편차 1 (표준편차는 n−1)
 * - 마할라노비스 이상치: D² = (x−x̄)ᵀ S⁻¹ (x−x̄), S 는 표본 공분산(n−1). 자유도 k(변수 수)의 χ² 근사로 p < .001 을 이상치 후보로 본다.
 * - 결측은 목록별 제거. 사례 간 거리는 n²개라 한도(DIST_MAX_PAIRS_CASES_)를 둔다. 화면에는 앞 50×50까지, 전체는 새 탭에 저장(500사례까지).
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function distDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return { getVars: data.getVars, pchisq: dist.pchisq, qchisq: dist.qchisq, choleskySolve_: require('./Stats_Mediation.js').choleskySolve_ };
  }
  return { getVars: getVars, pchisq: pchisq, qchisq: qchisq, choleskySolve_: choleskySolve_ };
}

var DIST_SHOW_MAX_ = 50;             // 화면에 보이는 행렬의 최대 크기
var DIST_SAVE_MAX_ = 500;            // 새 탭에 저장할 수 있는 최대 사례(또는 변수) 수
var DIST_MAX_PAIRS_CASES_ = 1000;    // 사례 간 거리를 계산할 수 있는 최대 사례 수 (n² 쌍)
var DIST_OUTLIER_P_ = 0.001;         // 이상치 후보 기준 (관행)
var DIST_OUTLIER_LIST_ = 20;         // 이상치 표에 보이는 최대 사례 수

/** [키, 이름, 종류('distance'|'similarity'), 이분 자료 전용 여부, 단위에 민감한지] */
var DIST_MEASURES_ = {
  euclid: ['유클리드 거리', 'distance', false, true],
  sqeuclid: ['제곱 유클리드 거리', 'distance', false, true],
  cityblock: ['맨해튼(블록) 거리', 'distance', false, true],
  chebyshev: ['체비쇼프 거리', 'distance', false, true],
  minkowski: ['민코프스키 거리', 'distance', false, true],
  mahalanobis: ['마할라노비스 거리', 'distance', false, false],
  cosineDist: ['코사인 거리 (1 − 코사인 유사도)', 'distance', false, true],
  cosine: ['코사인 유사도', 'similarity', false, true],
  corrDist: ['상관 거리 (1 − r)', 'distance', false, false],
  correlation: ['피어슨 상관 (유사도)', 'similarity', false, false],
  matching: ['단순일치 계수 (유사도)', 'similarity', true, false],
  jaccard: ['Jaccard 계수 (유사도)', 'similarity', true, false],
  dice: ['Dice 계수 (유사도)', 'similarity', true, false],
  hamming: ['Hamming 거리 (다른 변수의 비율)', 'distance', true, false]
};
var DIST_STD_ = {
  none: '표준화 없음', z: 'z점수', range01: '범위 0~1', range11: '범위 −1~1', max1: '최댓값 1', mean1: '평균 1', sd1: '표준편차 1'
};

function distName_(v) { return v.label ? v.label : v.name; }

/**
 * 변수(열)별 표준화. cols: [[값...]] (결측 없음), names: 오류 안내용 변수 이름. 새 배열을 돌려준다.
 */
function standardizeColumns(cols, kind, names) {
  if (kind === 'none' || kind === undefined || kind === null) return cols.map(function (c) { return c.slice(); });
  return cols.map(function (c, j) {
    var n = c.length, i, sum = 0, mn = Infinity, mx = -Infinity;
    for (i = 0; i < n; i++) { sum += c[i]; if (c[i] < mn) mn = c[i]; if (c[i] > mx) mx = c[i]; }
    var mean = sum / n, ss = 0;
    for (i = 0; i < n; i++) ss += (c[i] - mean) * (c[i] - mean);
    var sd = n > 1 ? Math.sqrt(ss / (n - 1)) : 0;
    var nm = names && names[j] ? '"' + names[j] + '"' : '변수';
    var out = new Array(n);
    if (kind === 'z' || kind === 'sd1') {
      if (!(sd > 0)) throw new Error(nm + '의 값이 모두 같아 표준화할 수 없습니다 (표준편차 0).');
      for (i = 0; i < n; i++) out[i] = kind === 'z' ? (c[i] - mean) / sd : c[i] / sd;
    } else if (kind === 'range01' || kind === 'range11') {
      if (!(mx > mn)) throw new Error(nm + '의 값이 모두 같아 표준화할 수 없습니다 (범위 0).');
      for (i = 0; i < n; i++) out[i] = kind === 'range01' ? (c[i] - mn) / (mx - mn) : 2 * (c[i] - mn) / (mx - mn) - 1;
    } else if (kind === 'max1') {
      if (!(mx !== 0)) throw new Error(nm + '의 최댓값이 0이라 최댓값 1로 표준화할 수 없습니다.');
      for (i = 0; i < n; i++) out[i] = c[i] / mx;
    } else if (kind === 'mean1') {
      if (!(mean !== 0)) throw new Error(nm + '의 평균이 0이라 평균 1로 표준화할 수 없습니다.');
      for (i = 0; i < n; i++) out[i] = c[i] / mean;
    } else {
      throw new Error('알 수 없는 표준화 방식입니다: ' + kind);
    }
    return out;
  });
}

/** 열들로 표본 공분산 행렬 (n−1) */
function distCovariance_(cols) {
  var k = cols.length, n = cols[0].length, means = [], i, j, t;
  for (j = 0; j < k; j++) { var s = 0; for (t = 0; t < n; t++) s += cols[j][t]; means.push(s / n); }
  var C = [];
  for (i = 0; i < k; i++) {
    C.push(new Array(k));
    for (j = 0; j <= i; j++) {
      var acc = 0;
      for (t = 0; t < n; t++) acc += (cols[i][t] - means[i]) * (cols[j][t] - means[j]);
      C[i][j] = acc / (n - 1);
      if (j < i) C[j][i] = C[i][j];
    }
  }
  return { C: C, means: means };
}

/** 공분산 행렬의 역행렬 (촐레스키). 특이하면 null */
function distInverse_(C) {
  var d = distDeps_(), k = C.length, id = [], i;
  for (i = 0; i < k; i++) { var e = new Array(k).fill(0); e[i] = 1; id.push(e); }
  var cols = d.choleskySolve_(C, id);
  if (!cols) return null;
  var inv = [];
  for (i = 0; i < k; i++) inv.push(cols.map(function (c) { return c[i]; }));
  return inv;
}

/**
 * 마할라노비스 D² (각 사례에서 평균 벡터까지). cols: 변수별 값(결측 없음).
 * 반환: { d2:[...], k, n } 또는 { singular: true } / { tooFew: true }
 */
function mahalanobisD2(cols) {
  var k = cols.length, n = k ? cols[0].length : 0;
  if (n < k + 1) return { tooFew: true };
  var cv = distCovariance_(cols), inv = distInverse_(cv.C);
  if (!inv) return { singular: true };
  var d2 = new Array(n), t, i, j;
  for (t = 0; t < n; t++) {
    var q = 0;
    for (i = 0; i < k; i++) {
      var di = cols[i][t] - cv.means[i];
      for (j = 0; j < k; j++) q += di * inv[i][j] * (cols[j][t] - cv.means[j]);
    }
    d2[t] = q < 0 ? 0 : q;
  }
  return { d2: d2, k: k, n: n };
}

/**
 * 벡터 쌍의 측도 함수를 만든다. 계산할 수 없으면(길이 0 벡터의 코사인, 분산 0의 상관 등) null 을 돌려주는 함수.
 * opts: { p (민코프스키), inv (마할라노비스의 공분산 역행렬) }
 */
function distMeasureFn_(kind, opts) {
  var p = opts && opts.p, inv = opts && opts.inv;
  var sim = function (a, b) {
    var i, dot = 0, na = 0, nb = 0;
    for (i = 0; i < a.length; i++) { dot += a[i] * b[i]; na += a[i] * a[i]; nb += b[i] * b[i]; }
    return na > 0 && nb > 0 ? dot / Math.sqrt(na * nb) : null;
  };
  var corr = function (a, b) {
    var i, n = a.length, ma = 0, mb = 0;
    for (i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
    ma /= n; mb /= n;
    var sab = 0, saa = 0, sbb = 0;
    for (i = 0; i < n; i++) { sab += (a[i] - ma) * (b[i] - mb); saa += (a[i] - ma) * (a[i] - ma); sbb += (b[i] - mb) * (b[i] - mb); }
    return saa > 0 && sbb > 0 ? sab / Math.sqrt(saa * sbb) : null;
  };
  var counts = function (a, b) {
    var i, n11 = 0, n00 = 0, mis = 0;
    for (i = 0; i < a.length; i++) { if (a[i] === 1 && b[i] === 1) n11++; else if (a[i] === 0 && b[i] === 0) n00++; else mis++; }
    return { a: n11, d: n00, bc: mis, n: a.length };
  };
  switch (kind) {
    case 'euclid': case 'sqeuclid':
      return function (a, b) { var s = 0, i; for (i = 0; i < a.length; i++) s += (a[i] - b[i]) * (a[i] - b[i]); return kind === 'euclid' ? Math.sqrt(s) : s; };
    case 'cityblock':
      return function (a, b) { var s = 0, i; for (i = 0; i < a.length; i++) s += Math.abs(a[i] - b[i]); return s; };
    case 'chebyshev':
      return function (a, b) { var m = 0, i; for (i = 0; i < a.length; i++) { var v = Math.abs(a[i] - b[i]); if (v > m) m = v; } return m; };
    case 'minkowski':
      return function (a, b) { var s = 0, i; for (i = 0; i < a.length; i++) s += Math.pow(Math.abs(a[i] - b[i]), p); return Math.pow(s, 1 / p); };
    case 'mahalanobis':
      return function (a, b) {
        var k = a.length, i, j, q = 0, df = new Array(k);
        for (i = 0; i < k; i++) df[i] = a[i] - b[i];
        for (i = 0; i < k; i++) for (j = 0; j < k; j++) q += df[i] * inv[i][j] * df[j];
        return Math.sqrt(q < 0 ? 0 : q);
      };
    case 'cosine': return sim;
    case 'cosineDist': return function (a, b) { var s = sim(a, b); return s === null ? null : 1 - s; };
    case 'correlation': return corr;
    case 'corrDist': return function (a, b) { var r = corr(a, b); return r === null ? null : 1 - r; };
    case 'matching': return function (a, b) { var c = counts(a, b); return (c.a + c.d) / c.n; };
    case 'jaccard': return function (a, b) { var c = counts(a, b); return c.a + c.bc > 0 ? c.a / (c.a + c.bc) : null; };
    case 'dice': return function (a, b) { var c = counts(a, b); return 2 * c.a + c.bc > 0 ? 2 * c.a / (2 * c.a + c.bc) : null; };
    case 'hamming': return function (a, b) { var c = counts(a, b); return c.bc / c.n; };
    default: throw new Error('알 수 없는 측도입니다: ' + kind);
  }
}

/**
 * 근접도 행렬. vectors: 벡터들의 배열(벡터 하나는 숫자 배열). 대칭이라 위 삼각만 계산한다.
 * 대각선: 거리는 0, 유사도는 자기 자신과의 값(코사인·상관·이분 계수는 1; 계산할 수 없으면 null).
 */
function proximityMatrix(vectors, kind, opts) {
  var f = distMeasureFn_(kind, opts), m = vectors.length, M = [], i, j;
  var isSim = DIST_MEASURES_[kind][1] === 'similarity';
  for (i = 0; i < m; i++) M.push(new Array(m));
  for (i = 0; i < m; i++) {
    M[i][i] = isSim ? f(vectors[i], vectors[i]) : 0;
    for (j = i + 1; j < m; j++) { var v = f(vectors[i], vectors[j]); M[i][j] = v; M[j][i] = v; }
  }
  return M;
}

function distAllBinary_(cols) {
  return cols.every(function (c) { return c.every(function (x) { return x === 0 || x === 1; }); });
}

function runDistance(dataset, request) {
  var d = distDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < 2) throw new Error('변수를 2개 이상 선택하세요.');
  var opt = request.options || {};
  var target = opt.target === 'var' ? 'var' : 'case';
  var measure = opt.measure || 'euclid';
  if (!DIST_MEASURES_[measure]) throw new Error('알 수 없는 측도입니다: ' + measure);
  var std = opt.standardize || 'none';
  if (!DIST_STD_[std]) throw new Error('알 수 없는 표준화 방식입니다: ' + std);
  var p = opt.p === undefined || opt.p === null || opt.p === '' ? 3 : Number(opt.p);
  if (measure === 'minkowski' && !(isFinite(p) && p >= 1)) throw new Error('민코프스키 거리의 차수 p는 1 이상의 숫자여야 합니다.');
  if (measure === 'mahalanobis' && target === 'var') throw new Error('마할라노비스 거리는 사례 간 거리에만 쓸 수 있습니다 (변수 간 거리에는 "사례 간"을 선택하세요).');
  var wantOutlier = !!opt.outlier || !!opt.saveD2;

  var rowsOk = [], i, a;
  for (i = 0; i < dataset.n; i++) {
    var ok = true;
    for (a = 0; a < vars.length; a++) if (dataset.data[vars[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length, k = vars.length;
  var warnings = [];
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건).');
  if (n < 2) throw new Error('사용할 수 있는 사례가 2건 이상 있어야 합니다 (현재 ' + n + '건).');
  var names = vars.map(distName_);
  var cols = vars.map(function (v) { return rowsOk.map(function (r) { return dataset.data[v.name][r]; }); });
  if (DIST_MEASURES_[measure][2] && !distAllBinary_(cols)) throw new Error('"' + DIST_MEASURES_[measure][0] + '"은(는) 값이 0 또는 1인 변수에만 쓸 수 있습니다. 다른 값이 있는 변수가 포함돼 있습니다.');
  // 사례 간 거리는 n²개라 사례가 많으면 행렬을 계산하지 않는다 (마할라노비스 이상치 점검만 가능)
  var skipMatrix = target === 'case' && n > DIST_MAX_PAIRS_CASES_;
  if (skipMatrix && !wantOutlier) {
    throw new Error('사례 간 거리는 사례가 ' + DIST_MAX_PAIRS_CASES_ + '건 이하일 때 계산합니다 (현재 ' + n + '건, 쌍이 ' + Math.round(n * (n - 1) / 2) + '개). 변수 간 거리를 쓰거나 사례를 줄이세요. 마할라노비스 이상치 점검은 사례가 많아도 할 수 있습니다.');
  }
  if (skipMatrix && opt.saveMatrix) throw new Error('사례가 ' + n + '건이라 행렬을 저장할 수 없습니다 (최대 ' + DIST_SAVE_MAX_ + '개).');
  if (skipMatrix) warnings.push('사례가 ' + n + '건으로 많아(한도 ' + DIST_MAX_PAIRS_CASES_ + '건) 거리 행렬은 계산하지 않고 마할라노비스 이상치 점검만 했습니다.');
  var rowLabel = function (idx) { return '행 ' + dataset.rowNumbers[rowsOk[idx]]; };

  // 단위 경고: 표준화하지 않았는데 단위에 민감한 측도이고 변수의 표준편차가 크게 다르면
  if (std === 'none' && DIST_MEASURES_[measure][3]) {
    var sds = cols.map(function (c) { var m = 0, s = 0, t; for (t = 0; t < c.length; t++) m += c[t]; m /= c.length; for (t = 0; t < c.length; t++) s += (c[t] - m) * (c[t] - m); return Math.sqrt(s / (c.length - 1)); });
    var pos = sds.filter(function (s) { return s > 0; });
    if (pos.length >= 2 && Math.max.apply(null, pos) / Math.min.apply(null, pos) > 10 && !DIST_MEASURES_[measure][2]) {
      warnings.push('변수의 표준편차가 크게 다릅니다 (최대/최소 ' + (Math.max.apply(null, pos) / Math.min.apply(null, pos)).toFixed(1) + '배). 표준화하지 않으면 단위가 큰 변수가 거리를 지배합니다. 옵션의 표준화(z점수 등)를 고려하세요.');
    }
  }

  var stdCols = DIST_MEASURES_[measure][2] ? cols : standardizeColumns(cols, std, names);
  var vectors, itemLabels;
  if (target === 'case') {
    vectors = rowsOk.map(function (r, idx) { return stdCols.map(function (c) { return c[idx]; }); });
    itemLabels = rowsOk.map(function (r, idx) { return rowLabel(idx); });
  } else {
    vectors = stdCols;
    itemLabels = names;
  }
  var mOpts = { p: p };
  if (measure === 'mahalanobis' && !skipMatrix) {
    var cv = distCovariance_(stdCols), inv = distInverse_(cv.C);
    if (!inv || n < k + 1) throw new Error('변수들의 공분산 행렬을 계산할 수 없습니다 (변수가 서로 완전히 겹치거나 사례가 변수보다 충분히 많지 않습니다). 마할라노비스 거리에는 사례 수가 변수 수보다 많아야 하고 변수 사이 완전한 공선성이 없어야 합니다.');
    mOpts.inv = inv;
  }
  var M = skipMatrix ? null : proximityMatrix(vectors, measure, mOpts);
  var m = vectors.length;
  var kindWord = DIST_MEASURES_[measure][1] === 'similarity' ? '유사도' : '거리';
  var mName = DIST_MEASURES_[measure][0] + (measure === 'minkowski' ? ' (p = ' + p + ')' : '');
  var targetWord = target === 'case' ? '사례 간' : '변수 간';

  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [
      { label: '분석에 포함한 사례', n: n }, { label: '결측으로 제외', n: dataset.n - n }, { label: '변수', n: k },
      { _kind: 'total', label: targetWord + ' 쌍 수', n: m * (m - 1) / 2 }
    ],
    notes: ['측도: ' + mName + ', ' + targetWord + ', 표준화: ' + (DIST_MEASURES_[measure][2] ? '쓰지 않음(이분 자료)' : DIST_STD_[std]) + '. 목록별 결측 제거: 선택한 모든 변수가 유효한 사례만 썼다.']
  });

  var result = { analysis: 'distance', title: '거리·유사도', blocks: blocks, warnings: warnings };
  if (!skipMatrix) distMatrixBlocks_(blocks, M, m, itemLabels, mName, kindWord, targetWord);
  return distFinish_(result, { dataset: dataset, rowsOk: rowsOk, cols: cols, names: names, n: n, k: k, opt: opt, rowLabel: rowLabel, wantOutlier: wantOutlier, M: M, m: m, itemLabels: itemLabels, mName: mName, kindWord: kindWord, targetWord: targetWord, target: target });
}

/** 행렬 표, 요약 표를 blocks 에 더한다 */
function distMatrixBlocks_(blocks, M, m, itemLabels, mName, kindWord, targetWord) {
  var i, j;
  // 행렬 표시 (앞 DIST_SHOW_MAX_ 개까지)
  var shown = Math.min(m, DIST_SHOW_MAX_);
  var mcols = [{ key: 'label', label: '' }];
  for (var cj = 0; cj < shown; cj++) mcols.push({ key: 'c' + cj, label: itemLabels[cj], format: { type: 'fixed', digits: 3 } });
  var mrows = [];
  for (i = 0; i < shown; i++) {
    var row = { label: itemLabels[i] };
    for (var cj2 = 0; cj2 < shown; cj2++) row['c' + cj2] = M[i][cj2];
    mrows.push(row);
  }
  var mnotes = [kindWord === '거리' ? '값이 클수록 서로 다르다 (대각선은 0).' : '값이 클수록 서로 비슷하다 (대각선은 자기 자신과의 값).'];
  if (m > shown) mnotes.push('화면에는 앞 ' + shown + '×' + shown + '만 보인다 (전체 ' + m + '×' + m + '). 전체는 옵션 "근접도 행렬을 새 탭에 저장"으로 저장할 수 있다 (' + DIST_SAVE_MAX_ + '개까지). 아래 요약 통계는 전체를 계산한 값이다.');
  blocks.push({ type: 'table', caption: mName + ' 행렬 (' + targetWord + ')', columns: mcols, rows: mrows, notes: mnotes });

  // 요약: 쌍의 최소·평균·최대, 가장 가까운(비슷한)·먼 쌍
  var cnt = 0, sum = 0, lo = null, hi = null, loPair = null, hiPair = null;
  for (i = 0; i < m; i++) for (j = i + 1; j < m; j++) {
    var v = M[i][j];
    if (v === null) continue;
    cnt++; sum += v;
    if (lo === null || v < lo) { lo = v; loPair = [i, j]; }
    if (hi === null || v > hi) { hi = v; hiPair = [i, j]; }
  }
  var closeIdx = kindWord === '거리' ? loPair : hiPair, farIdx = kindWord === '거리' ? hiPair : loPair;
  var pairText = function (pr) { return pr ? itemLabels[pr[0]] + ' — ' + itemLabels[pr[1]] : ''; };
  blocks.push({
    type: 'table', caption: mName + ' 요약 (쌍 ' + cnt + '개)',
    columns: [{ key: 'label', label: '' }, { key: 'value', label: '값', format: { type: 'fixed', digits: 3 } }, { key: 'pair', label: '쌍' }],
    rows: [
      { label: '가장 ' + (kindWord === '거리' ? '가까운' : '비슷한') + ' 쌍', value: kindWord === '거리' ? lo : hi, pair: pairText(closeIdx) },
      { label: '가장 ' + (kindWord === '거리' ? '먼' : '다른') + ' 쌍', value: kindWord === '거리' ? hi : lo, pair: pairText(farIdx) },
      { label: '평균', value: cnt ? sum / cnt : null, pair: '' }
    ],
    notes: cnt < m * (m - 1) / 2 ? ['계산할 수 없는 쌍(길이 0인 벡터의 코사인, 변하지 않는 벡터의 상관 등) ' + (m * (m - 1) / 2 - cnt) + '개는 제외했다.'] : []
  });
}

/** 마할라노비스 이상치 점검, 저장 열, 새 탭에 저장할 행렬을 결과에 더한다 */
function distFinish_(result, c) {
  var d = distDeps_(), blocks = result.blocks, opt = c.opt, n = c.n, k = c.k, dataset = c.dataset, rowsOk = c.rowsOk, rowLabel = c.rowLabel;
  var M = c.M, m = c.m, itemLabels = c.itemLabels, names = c.names, i;

  // 마할라노비스 이상치 점검
  if (c.wantOutlier) {
    var md = mahalanobisD2(c.cols);       // 이상치는 표준화와 무관하다 (마할라노비스는 단위에 불변)
    if (md.tooFew) throw new Error('마할라노비스 이상치 점검에는 사례 수가 변수 수보다 많아야 합니다 (사례 ' + n + '건, 변수 ' + k + '개).');
    if (md.singular) throw new Error('변수들의 공분산 행렬을 계산할 수 없습니다 (변수가 서로 완전히 겹칩니다). 같은 정보를 담은 변수를 함께 쓰지 않았는지 확인하세요.');
    var crit = d.qchisq(DIST_OUTLIER_P_, k, false);
    var pv = md.d2.map(function (q) { return d.pchisq(q, k, false); });
    var order = md.d2.map(function (q, idx) { return idx; }).sort(function (x, y) { return md.d2[y] - md.d2[x]; });
    var flagged = order.filter(function (idx) { return pv[idx] < DIST_OUTLIER_P_; });
    var listIdx = flagged.length ? flagged.slice(0, DIST_OUTLIER_LIST_) : order.slice(0, Math.min(5, order.length));
    blocks.push({
      type: 'table', caption: '마할라노비스 이상치 점검 (D²)',
      columns: [{ key: 'label', label: '' }, { key: 'value', label: '값', format: { type: 'fixed', digits: 3 } }],
      rows: [
        { label: '변수 수 (자유도)', value: k }, { label: '사례 수', value: n },
        { label: '임계값 (χ², p = ' + DIST_OUTLIER_P_ + ')', value: crit }, { label: '최대 D²', value: md.d2[order[0]] },
        { _kind: 'total', label: '이상치 후보 (p < ' + DIST_OUTLIER_P_ + ')', value: flagged.length }
      ].map(function (r) { if (r.label === '변수 수 (자유도)' || r.label === '사례 수' || r.label.indexOf('이상치 후보') === 0) r.value = Math.round(r.value); return r; }),
      notes: ['D²는 각 사례가 변수 평균 벡터에서 얼마나 먼지를 변수 사이 상관까지 반영해 잰 값이다. 정규분포에 가까우면 D²는 자유도 k인 χ² 분포를 따른다. p < .001 이면 이상치 후보이며, 바로 지우지 말고 입력 오류나 다른 집단인지 먼저 확인한다.']
    });
    blocks.push({
      type: 'table', caption: flagged.length ? '이상치 후보 (D²가 큰 순' + (flagged.length > listIdx.length ? ', 앞 ' + listIdx.length + '건' : '') + ')' : '이상치 후보 없음 — D²가 가장 큰 사례 ' + listIdx.length + '건',
      columns: [{ key: 'label', label: '사례' }, { key: 'd2', label: 'D²', format: { type: 'fixed', digits: 3 } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }],
      rows: listIdx.map(function (idx) { return { label: rowLabel(idx), d2: md.d2[idx], p: pv[idx] }; }),
      notes: ['사례 이름은 데이터 탭의 행 번호이다. 유의확률은 χ² 근사의 오른쪽 꼬리 확률이다.']
    });
    if (opt.saveD2) {
      var full = new Array(dataset.n).fill(null);
      rowsOk.forEach(function (r, idx) { full[r] = md.d2[idx]; });
      result.saveColumns = [{ name: '마할라노비스_D2', label: '마할라노비스 D² (' + names.join(', ') + ')', scale: '연속', values: full, fixed: 4 }];
    }
  }

  // 근접도 행렬을 새 탭에 저장
  if (opt.saveMatrix) {
    if (m > DIST_SAVE_MAX_) throw new Error((c.target === 'case' ? '사례' : '변수') + '가 ' + m + '개라 행렬을 새 탭에 저장할 수 없습니다 (최대 ' + DIST_SAVE_MAX_ + '개). 사례를 줄이거나 변수 간 거리를 쓰세요.');
    var sheetRows = [[c.mName + ' (' + c.targetWord + ')'].concat(itemLabels)];
    for (i = 0; i < m; i++) sheetRows.push([itemLabels[i]].concat(M[i].map(function (x) { return x === null ? '' : x; })));
    result.saveSheet = { name: c.kindWord === '거리' ? '거리행렬' : '유사도행렬', rows: sheetRows };
  }
  return result;
}

if (typeof module !== 'undefined') {
  module.exports = {
    distCovariance_: distCovariance_, distInverse_: distInverse_,
    standardizeColumns: standardizeColumns, mahalanobisD2: mahalanobisD2, proximityMatrix: proximityMatrix, runDistance: runDistance,
    DIST_MEASURES_: DIST_MEASURES_, DIST_STD_: DIST_STD_
  };
}
