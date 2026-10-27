/**
 * k-평균 군집분석. 순수 함수: (dataset, request) → 결과 객체.
 *
 * - 변수는 선택에 따라 표준화(기본 z점수)한 공간에서 군집한다. 단위가 다르면 표준화하지 않은 변수가 거리를 지배한다.
 * - 초기 중심: k-means++ (첫 중심은 floor(u·n), 다음은 D²에 비례해 누적합 > u·합계인 첫 사례, 중심마다 난수 1개).
 *   난수는 mulberry32(시드)라 같은 시드이면 항상 같은 결과이다.
 * - 반복: Lloyd. 가장 가까운 중심(동점은 번호가 작은 쪽)에 배정하고 평균으로 중심을 갱신, 배정이 바뀌지 않으면 수렴.
 *   빈 군집이 생기면 자기 중심에서 가장 먼 사례(크기 2 이상인 군집에서)로 채운다.
 * - 여러 번 시작해 군집 내 제곱합이 가장 작은 결과를 쓴다. 군집 번호는 크기 내림차순(동점은 원래 번호 순).
 * - 군집 수 판단 보조(엘보): k = 2~10 의 군집 내 제곱합, 평균 실루엣(사례 2,000건 이하). 사례 5,000건 초과이면 생략.
 * - 변수별 F(일원 분산분석)는 군집이 데이터에서 만들어졌으므로 유의확률은 기술적 참고일 뿐이다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function kmDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var dd = require('./Stats_Distance.js');
    return {
      getVars: data.getVars, pf: dist.pf, mulberry32_: require('./Stats_Mediation.js').mulberry32_,
      standardizeColumns: dd.standardizeColumns, lineChartBlock: require('./Charts.js').lineChartBlock
    };
  }
  return { getVars: getVars, pf: pf, mulberry32_: mulberry32_, standardizeColumns: standardizeColumns, lineChartBlock: lineChartBlock };
}

var KM_MAX_K_ = 10;
var KM_MAX_STARTS_ = 50;
var KM_MAX_ITER_CAP_ = 500;
var KM_ELBOW_MAX_N_ = 5000;
var KM_SIL_MAX_N_ = 2000;
var KM_SMALL_CLUSTER_ = 0.05;
function kmName_(v) { return v.label ? v.label : v.name; }

function kmSq_(a, b) { var s = 0, j; for (j = 0; j < a.length; j++) s += (a[j] - b[j]) * (a[j] - b[j]); return s; }

/** k-means++ 초기 중심 (사례 인덱스 배열). X: 사례별 벡터, rnd: 난수 함수 */
function kmeansPlusPlus(X, K, rnd) {
  var n = X.length, i, c;
  var chosen = [Math.floor(rnd() * n)];
  var d2 = new Array(n);
  for (i = 0; i < n; i++) d2[i] = kmSq_(X[i], X[chosen[0]]);
  for (c = 1; c < K; c++) {
    var total = 0, pick = n - 1;
    for (i = 0; i < n; i++) total += d2[i];
    if (total === 0) {
      pick = Math.floor(rnd() * n);
    } else {
      var u = rnd() * total, cum = 0;
      for (i = 0; i < n; i++) { cum += d2[i]; if (cum > u) { pick = i; break; } }
    }
    chosen.push(pick);
    for (i = 0; i < n; i++) { var v = kmSq_(X[i], X[pick]); if (v < d2[i]) d2[i] = v; }
  }
  return chosen;
}

/** Lloyd 반복. centers: 초기 중심 벡터들(복사해서 쓴다). 반환: { labels, centers, wss, iterations, converged } */
function kmeansLloyd(X, centers, maxIter) {
  var n = X.length, K = centers.length, dim = X[0].length;
  var cen = centers.map(function (c) { return c.slice(); });
  var labels = new Array(n).fill(-1), dist = new Array(n), counts, i, j, d, it, converged = false;
  for (it = 1; it <= maxIter; it++) {
    var changed = false;
    counts = new Array(K).fill(0);
    for (i = 0; i < n; i++) {
      var best = 0, bd = kmSq_(X[i], cen[0]);
      for (j = 1; j < K; j++) { d = kmSq_(X[i], cen[j]); if (d < bd) { bd = d; best = j; } }
      if (labels[i] !== best) changed = true;
      labels[i] = best; dist[i] = bd; counts[best]++;
    }
    for (j = 0; j < K; j++) {
      if (counts[j] > 0) continue;
      var far = -1;
      for (i = 0; i < n; i++) if (counts[labels[i]] > 1 && (far < 0 || dist[i] > dist[far])) far = i;
      counts[labels[far]]--; labels[far] = j; counts[j] = 1; dist[far] = 0; changed = true;
    }
    for (j = 0; j < K; j++) {
      var s = new Array(dim).fill(0), t;
      for (i = 0; i < n; i++) if (labels[i] === j) for (t = 0; t < dim; t++) s[t] += X[i][t];
      for (t = 0; t < dim; t++) cen[j][t] = s[t] / counts[j];
    }
    if (!changed) { converged = true; break; }
  }
  var wss = 0;
  for (i = 0; i < n; i++) wss += kmSq_(X[i], cen[labels[i]]);
  return { labels: labels, centers: cen, wss: wss, iterations: it > maxIter ? maxIter : it, converged: converged };
}

/**
 * 여러 번 시작해 가장 좋은 결과. 반환: { labels(크기 내림차순 번호 0부터), centers, wss, iterations, converged, init(사례 인덱스), bestStart(1부터), sizes }
 */
function kmeansRun(X, K, seed, starts, maxIter) {
  var d = kmDeps_();
  var rnd = d.mulberry32_(seed), best = null, s;
  for (s = 0; s < starts; s++) {
    var init = kmeansPlusPlus(X, K, rnd);
    var r = kmeansLloyd(X, init.map(function (i) { return X[i]; }), maxIter);
    if (best === null || r.wss < best.wss) { r.init = init; r.bestStart = s + 1; best = r; }
  }
  var sizes = new Array(K).fill(0), j;
  best.labels.forEach(function (l) { sizes[l]++; });
  var order = [];
  for (j = 0; j < K; j++) order.push(j);
  order.sort(function (a, b) { return sizes[b] - sizes[a] || a - b; });
  var remap = new Array(K);
  order.forEach(function (old, nw) { remap[old] = nw; });
  best.labels = best.labels.map(function (l) { return remap[l]; });
  best.centers = order.map(function (old) { return best.centers[old]; });
  best.sizes = order.map(function (old) { return sizes[old]; });
  return best;
}

/** 평균 실루엣 (유클리드). 군집이 1개 사례뿐이면 그 사례의 실루엣은 0 */
function kmeansSilhouette(X, labels, K) {
  var n = X.length, i, j, sum = 0;
  var sums = new Array(K), cnt = new Array(K).fill(0);
  for (i = 0; i < n; i++) cnt[labels[i]]++;
  for (i = 0; i < n; i++) {
    for (j = 0; j < K; j++) sums[j] = 0;
    for (j = 0; j < n; j++) if (j !== i) sums[labels[j]] += Math.sqrt(kmSq_(X[i], X[j]));
    var own = labels[i];
    if (cnt[own] <= 1) continue;
    var a = sums[own] / (cnt[own] - 1), b = Infinity;
    for (j = 0; j < K; j++) if (j !== own && cnt[j] > 0) { var m = sums[j] / cnt[j]; if (m < b) b = m; }
    var mx = Math.max(a, b);
    sum += mx > 0 ? (b - a) / mx : 0;
  }
  return sum / n;
}

function runKmeans(dataset, request) {
  var d = kmDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < 2) throw new Error('변수를 2개 이상 선택하세요.');
  var opt = request.options || {};
  var num = function (v, def) { return v === undefined || v === null || v === '' ? def : Number(v); };
  var K = num(opt.k, 3), seed = num(opt.seed, 12345), starts = num(opt.starts, 10), maxIter = num(opt.maxIter, 100);
  if (!isFinite(K) || K !== Math.floor(K) || K < 2 || K > KM_MAX_K_) throw new Error('군집 수는 2~' + KM_MAX_K_ + '의 정수여야 합니다.');
  if (!isFinite(seed) || seed !== Math.floor(seed)) throw new Error('시드는 정수여야 합니다.');
  if (!isFinite(starts) || starts !== Math.floor(starts) || starts < 1 || starts > KM_MAX_STARTS_) throw new Error('시작 횟수는 1~' + KM_MAX_STARTS_ + '의 정수여야 합니다.');
  if (!isFinite(maxIter) || maxIter !== Math.floor(maxIter) || maxIter < 1 || maxIter > KM_MAX_ITER_CAP_) throw new Error('최대 반복 횟수는 1~' + KM_MAX_ITER_CAP_ + '의 정수여야 합니다.');
  var std = opt.standardize === 'none' || opt.standardize === 'range01' ? opt.standardize : 'z';

  var rowsOk = [], i, a;
  for (i = 0; i < dataset.n; i++) {
    var ok = true;
    for (a = 0; a < vars.length; a++) if (dataset.data[vars[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length, k = vars.length;
  var warnings = [];
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건).');
  if (n < K + 1) throw new Error('사례 수(' + n + ')가 군집 수(' + K + ')보다 충분히 많아야 합니다.');
  var names = vars.map(kmName_);
  var cols = vars.map(function (v) { return rowsOk.map(function (r) { return dataset.data[v.name][r]; }); });
  var work = d.standardizeColumns(cols, std, names);
  var X = rowsOk.map(function (r, idx) { return work.map(function (c) { return c[idx]; }); });

  if (std === 'none') {
    var sds = cols.map(function (c) { var m = 0, s = 0, t; for (t = 0; t < n; t++) m += c[t]; m /= n; for (t = 0; t < n; t++) s += (c[t] - m) * (c[t] - m); return Math.sqrt(s / (n - 1)); });
    var pos = sds.filter(function (s) { return s > 0; });
    if (pos.length >= 2 && Math.max.apply(null, pos) / Math.min.apply(null, pos) > 10) {
      warnings.push('변수의 표준편차가 크게 다릅니다 (최대/최소 ' + (Math.max.apply(null, pos) / Math.min.apply(null, pos)).toFixed(1) + '배). 표준화하지 않으면 단위가 큰 변수가 군집을 좌우합니다. 표준화(z점수)를 고려하세요.');
    }
  }

  var res = kmeansRun(X, K, seed, starts, maxIter);
  var mean = new Array(k).fill(0), t;
  X.forEach(function (x) { for (t = 0; t < k; t++) mean[t] += x[t] / n; });
  var tss = 0;
  X.forEach(function (x) { tss += kmSq_(x, mean); });
  if (!res.converged) warnings.push('최대 반복(' + maxIter + '회) 안에 수렴하지 않았습니다. 최대 반복 횟수를 늘려 보세요.');
  var minShare = Math.min.apply(null, res.sizes) / n;
  if (minShare < KM_SMALL_CLUSTER_) warnings.push('가장 작은 군집이 전체의 ' + (minShare * 100).toFixed(1) + '%(' + Math.min.apply(null, res.sizes) + '건)뿐입니다. 이상치 몇 개가 군집이 됐거나 군집 수가 많을 수 있습니다.');

  var F3 = { type: 'fixed', digits: 3 };
  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [{ label: '분석에 포함', n: n }, { label: '결측으로 제외', n: dataset.n - n }, { _kind: 'total', label: '합계', n: dataset.n }],
    notes: ['목록별 결측 제거: 선택한 모든 변수가 유효한 사례만 썼다. 표준화: ' + (std === 'z' ? 'z점수' : std === 'range01' ? '범위 0~1' : '하지 않음') + ' (변수별).']
  });
  blocks.push({
    type: 'table', caption: '군집 결과 요약',
    columns: [{ key: 'label', label: '' }, { key: 'value', label: '값', format: F3 }],
    rows: [
      { label: '군집 수', value: K, _int: true }, { label: '시작 횟수 (가장 좋은 시작 번호)', value: starts, _txt: starts + ' (' + res.bestStart + '번째)' },
      { label: '반복 횟수', value: res.iterations, _int: true }, { label: '수렴 여부', value: null, _txt: res.converged ? '수렴' : '수렴하지 않음' },
      { label: '총 제곱합 (TSS)', value: tss }, { label: '군집 내 제곱합 (WSS)', value: res.wss },
      { _kind: 'total', label: '설명 비율 (1 − WSS/TSS)', value: tss > 0 ? 1 - res.wss / tss : null }
    ].map(function (r) { if (r._int) r.value = Math.round(r.value); if (r._txt !== undefined) { r.value = r._txt; } delete r._int; delete r._txt; return r; }),
    notes: ['제곱합은 군집한 공간(표준화했으면 표준화한 값)에서 계산했다. 난수 시드 ' + seed + '로 같은 설정에서는 항상 같은 결과가 나온다. 초기 중심은 k-means++로 고르고 군집 번호는 큰 군집부터 매겼다 (초기 중심을 다르게 고르는 프로그램과는 결과가 다를 수 있다).']
  });

  var cenCols = [{ key: 'label', label: '군집' }];
  names.forEach(function (nm, j) { cenCols.push({ key: 'v' + j, label: nm, format: F3 }); });
  var cenRows = res.centers.map(function (c, g) { var r = { label: '군집 ' + (g + 1) }; c.forEach(function (v, j) { r['v' + j] = v; }); return r; });
  blocks.push({
    type: 'table', caption: '최종 군집 중심' + (std === 'none' ? '' : ' (표준화한 값)'), columns: cenCols, rows: cenRows,
    notes: std === 'none' ? [] : ['표준화한 값의 군집 평균이다. 0은 전체 평균, 양수는 평균보다 높음을 뜻한다. 원래 단위는 다음 표.']
  });
  var origCenters = res.centers.map(function (c, g) {
    var m = new Array(k).fill(0);
    rowsOk.forEach(function (r, idx) { if (res.labels[idx] === g) for (t = 0; t < k; t++) m[t] += cols[t][idx] / res.sizes[g]; });
    return m;
  });
  if (std !== 'none') {
    blocks.push({
      type: 'table', caption: '최종 군집 중심 (원래 단위)', columns: cenCols,
      rows: origCenters.map(function (c, g) { var r = { label: '군집 ' + (g + 1) }; c.forEach(function (v, j) { r['v' + j] = v; }); return r; })
    });
  }

  // 군집별 사례 수, 군집 내 제곱합
  var wssBy = new Array(K).fill(0), distTo = new Array(n);
  X.forEach(function (x, idx) { var q = kmSq_(x, res.centers[res.labels[idx]]); wssBy[res.labels[idx]] += q; distTo[idx] = Math.sqrt(q); });
  blocks.push({
    type: 'table', caption: '군집별 사례 수',
    columns: [{ key: 'label', label: '군집' }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'pct', label: '%', format: { type: 'fixed', digits: 1 } },
      { key: 'wss', label: '군집 내 제곱합', format: F3 }, { key: 'avgd', label: '중심까지 평균 거리', format: F3 }],
    rows: res.sizes.map(function (sz, g) { return { label: '군집 ' + (g + 1), n: sz, pct: sz / n * 100, wss: wssBy[g], avgd: Math.sqrt(wssBy[g] / sz) }; })
      .concat([{ _kind: 'total', label: '합계', n: n, pct: 100, wss: res.wss, avgd: null }]),
    notes: ['"중심까지 평균 거리"는 군집 내 제곱합을 사례 수로 나눈 값의 제곱근(RMS 거리)이다.']
  });

  // 군집 중심 사이 거리
  var cdCols = [{ key: 'label', label: '' }];
  res.centers.forEach(function (c, g) { cdCols.push({ key: 'c' + g, label: '군집 ' + (g + 1), format: F3 }); });
  blocks.push({
    type: 'table', caption: '군집 중심 사이의 거리', columns: cdCols,
    rows: res.centers.map(function (c, g) { var r = { label: '군집 ' + (g + 1) }; res.centers.forEach(function (c2, h) { r['c' + h] = Math.sqrt(kmSq_(c, c2)); }); return r; })
  });

  // 변수별 F (참고)
  var fRows = vars.map(function (v, j) {
    var gm = 0, ss = 0, sb = 0, g, idx;
    for (idx = 0; idx < n; idx++) gm += cols[j][idx] / n;
    for (g = 0; g < K; g++) {
      var m = origCenters[g][j];
      sb += res.sizes[g] * (m - gm) * (m - gm);
    }
    for (idx = 0; idx < n; idx++) { var dm = cols[j][idx] - origCenters[res.labels[idx]][j]; ss += dm * dm; }
    var df1 = K - 1, df2 = n - K, F = ss > 0 ? (sb / df1) / (ss / df2) : null;
    return { name: names[j], F: F, df1: df1, df2: df2, p: F === null ? null : d.pf(F, df1, df2, false) };
  });
  blocks.push({
    type: 'table', caption: '변수별 군집 간 차이 (일원 분산분석, 참고)',
    columns: [{ key: 'name', label: '변수' }, { key: 'F', label: 'F', format: F3 }, { key: 'df1', label: '자유도1', format: { type: 'int' } }, { key: 'df2', label: '자유도2', format: { type: 'int' } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }],
    rows: fRows,
    notes: ['F가 큰 변수일수록 군집을 나누는 데 크게 기여했다. 군집이 이 자료에서 만들어졌으므로 유의확률은 검정으로 해석하지 말고 변수의 상대적 중요도를 보는 참고로만 쓴다.']
  });

  // 군집 수 판단 보조: k = 2~10 엘보와 평균 실루엣
  if (opt.elbow) {
    if (n > KM_ELBOW_MAX_N_) {
      warnings.push('사례가 ' + n + '건으로 많아(한도 ' + KM_ELBOW_MAX_N_ + '건) 군집 수 판단 보조(엘보)는 계산하지 않았습니다.');
    } else {
      var maxK = Math.min(KM_MAX_K_, n - 1), sil = n <= KM_SIL_MAX_N_, eRows = [], pts = [], kk;
      for (kk = 2; kk <= maxK; kk++) {
        var rr = kk === K ? res : kmeansRun(X, kk, seed, starts, maxIter);
        eRows.push({ k: kk, wss: rr.wss, explained: tss > 0 ? 1 - rr.wss / tss : null, sil: sil ? kmeansSilhouette(X, rr.labels, kk) : null, _bold: kk === K });
        pts.push([kk, rr.wss]);
      }
      blocks.push({
        type: 'table', caption: '군집 수 판단 보조 (k = 2~' + maxK + ')',
        columns: [{ key: 'k', label: '군집 수', format: { type: 'int' } }, { key: 'wss', label: '군집 내 제곱합', format: F3 }, { key: 'explained', label: '설명 비율', format: F3 }, { key: 'sil', label: '평균 실루엣', format: F3 }],
        rows: eRows,
        notes: ['굵은 줄이 선택한 군집 수이다. 군집 내 제곱합이 급격히 꺾이는 지점(엘보)과 평균 실루엣이 큰 군집 수를 참고한다 (실루엣 0.5 이상이면 구조가 뚜렷, 0.25 미만이면 약함). 정답이 있는 것은 아니고 해석 가능성도 함께 본다.' +
          (sil ? '' : ' 사례가 ' + KM_SIL_MAX_N_ + '건을 넘어 실루엣은 계산하지 않았다.')]
      });
      blocks.push(d.lineChartBlock('군집 수와 군집 내 제곱합 (엘보 그림)', '군집 수', '군집 내 제곱합', [{ label: '군집 내 제곱합', points: pts }]));
    }
  }

  var result = { analysis: 'kmeans', title: 'k-평균 군집분석', blocks: blocks, warnings: warnings };
  if (opt.save) {
    var lab = new Array(dataset.n).fill(null), dis = new Array(dataset.n).fill(null);
    rowsOk.forEach(function (r, idx) { lab[r] = res.labels[idx] + 1; dis[r] = distTo[idx]; });
    var vl = [];
    for (var g2 = 1; g2 <= K; g2++) vl.push(g2 + '=군집 ' + g2);
    result.saveColumns = [
      { name: '군집번호', label: '군집 번호 (k-평균, k=' + K + ')', scale: '명목', values: lab, valueLabels: vl.join(', ') },
      { name: '군집중심거리', label: '소속 군집 중심까지의 거리', scale: '연속', values: dis, fixed: 4 }
    ];
  }
  return result;
}

if (typeof module !== 'undefined') {
  module.exports = { kmeansPlusPlus: kmeansPlusPlus, kmeansLloyd: kmeansLloyd, kmeansRun: kmeansRun, kmeansSilhouette: kmeansSilhouette, runKmeans: runKmeans };
}
