/**
 * McNemar 검정 (같은 사람의 전·후 두 범주 변수, 2×2). 순수 함수: (dataset, request) → 결과 객체.
 *
 * - 두 변수(변수1 = 전, 변수2 = 후)의 값을 합쳐 정확히 2가지여야 한다. 작은 값을 0, 큰 값을 1(양성)로 본다.
 * - 표: a = (전 0, 후 0), b = (전 0, 후 1), c = (전 1, 후 0), d = (전 1, 후 1). 검정에는 값이 바뀐 쌍(불일치 쌍 b, c)만 쓴다.
 * - 정확 유의확률(양측): 불일치 쌍 n = b + c 에 대해 min(1, 2·P(X ≤ min(b, c))), X ~ Binomial(n, 0.5).
 * - 카이제곱(연속성 수정) = (|b − c| − 1)² / (b + c), 자유도 1. 불일치 쌍이 25건 이상이면 점근 근사도 쓸 만하다.
 * - 양성 비율의 변화 (b − c)/n 의 Wald 신뢰구간: ± z·√(b + c − (b − c)²/n)/n. 오즈비 b/c 의 정확 신뢰구간은
 *   불일치 쌍 중 b 의 비율에 Clopper–Pearson 구간을 구한 뒤 π/(1 − π) 로 바꾼다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function mcnDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return { getVars: data.getVars, compareValues_: data.compareValues_, pchisq: dist.pchisq, qnorm: dist.qnorm, ibetaPQ_: dist.ibetaPQ_ };
  }
  return { getVars: getVars, compareValues_: compareValues_, pchisq: pchisq, qnorm: qnorm, ibetaPQ_: ibetaPQ_ };
}

function mcnName_(v) { return v.label ? v.label : v.name; }

/** P(X ≤ k), X ~ Binomial(n, 0.5): 정규화된 불완전 베타 I_0.5(n − k, k + 1) */
function mcnemarBinomCdf_(k, n) {
  var d = mcnDeps_();
  if (k >= n) return 1;
  return d.ibetaPQ_(0.5, 0.5, n - k, k + 1)[0];
}

/** 정확 유의확률(양측). 불일치 쌍이 없으면 1 */
function mcnemarExactP(b, c) {
  var n = b + c;
  if (n === 0) return 1;
  return Math.min(1, 2 * mcnemarBinomCdf_(Math.min(b, c), n));
}

/** 베타(a, b)의 분위수 (이분법). 0 < p < 1 */
function mcnemarBetaQuantile_(p, a, b) {
  var d = mcnDeps_(), lo = 0, hi = 1, i;
  for (i = 0; i < 200; i++) {
    var mid = (lo + hi) / 2;
    if (d.ibetaPQ_(mid, 1 - mid, a, b)[0] < p) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}

/**
 * b, c, n(전체 쌍 수)와 신뢰수준으로 모든 통계량을 계산한다.
 * 반환: { n, nd, exactP, chi2, chiP, diff, diffLo, diffHi, odds, oddsLo, oddsHi } (계산할 수 없는 값은 null)
 */
function mcnemarTest(a, b, c, d, conf) {
  var dep = mcnDeps_();
  var n = a + b + c + d, nd = b + c;
  var res = { n: n, nd: nd, exactP: mcnemarExactP(b, c), chi2: null, chiP: null, diff: null, diffLo: null, diffHi: null, odds: null, oddsLo: null, oddsHi: null };
  if (nd > 0) {
    res.chi2 = Math.pow(Math.abs(b - c) - 1, 2) / nd;
    res.chiP = dep.pchisq(res.chi2, 1, false);
  }
  if (n > 0) {
    var z = dep.qnorm(1 - (1 - conf / 100) / 2);
    var se = Math.sqrt(Math.max(0, nd - (b - c) * (b - c) / n)) / n;
    res.diff = (b - c) / n;
    res.diffLo = res.diff - z * se; res.diffHi = res.diff + z * se;
  }
  if (nd > 0) {
    var alpha = (1 - conf / 100) / 2;
    var piLo = b > 0 ? mcnemarBetaQuantile_(alpha, b, c + 1) : 0;
    var piHi = c > 0 ? mcnemarBetaQuantile_(1 - alpha, b + 1, c) : 1;
    var odds = function (p) { return p >= 1 ? null : p / (1 - p); };
    res.odds = c > 0 ? b / c : null;
    res.oddsLo = odds(piLo); res.oddsHi = odds(piHi);
  }
  return res;
}

function runMcnemar(dataset, request) {
  var dp = mcnDeps_();
  var first = dp.getVars(dataset, request.first), second = dp.getVars(dataset, request.second);
  if (first.length !== 1 || second.length !== 1) throw new Error('변수1(전)과 변수2(후)를 하나씩 선택하세요.');
  var v1 = first[0], v2 = second[0];
  if (v1.name === v2.name) throw new Error('변수1과 변수2는 서로 다른 변수여야 합니다.');
  var conf = 95;
  var rowsOk = [], i;
  for (i = 0; i < dataset.n; i++) if (dataset.data[v1.name][i] !== null && dataset.data[v2.name][i] !== null) rowsOk.push(i);
  var n = rowsOk.length, warnings = [];
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건).');
  if (n < 2) throw new Error('사용할 수 있는 사례가 2건 이상 있어야 합니다 (현재 ' + n + '건).');
  var x = rowsOk.map(function (r) { return dataset.data[v1.name][r]; }), y = rowsOk.map(function (r) { return dataset.data[v2.name][r]; });
  var seen = {}, levels = [];
  x.concat(y).forEach(function (v) { var key = typeof v + ':' + v; if (!seen[key]) { seen[key] = true; levels.push({ value: v, key: key }); } });
  levels.sort(dp.compareValues_);
  if (levels.length !== 2) {
    throw new Error('McNemar 검정은 두 변수의 값이 합쳐서 정확히 2가지여야 합니다 (예: 0/1, 아니오/예). 지금은 ' + levels.length + '가지입니다' +
      (levels.length > 2 ? ' (' + levels.slice(0, 6).map(function (l) { return l.value; }).join(', ') + (levels.length > 6 ? ' …' : '') + '). 3가지 이상이면 값을 두 범주로 묶어 주세요.' : '. 값이 모두 같으면 변화를 볼 수 없습니다.'));
  }
  levels.forEach(function (l) {
    var lab = v1.valueLabels[String(l.value)];
    if (lab === undefined) lab = v2.valueLabels[String(l.value)];
    l.label = lab !== undefined ? lab : String(l.value);
  });
  var idx = {};
  levels.forEach(function (l, k) { idx[l.key] = k; });
  var cell = [[0, 0], [0, 0]];
  for (i = 0; i < n; i++) cell[idx[typeof x[i] + ':' + x[i]]][idx[typeof y[i] + ':' + y[i]]]++;
  var a = cell[0][0], b = cell[0][1], c = cell[1][0], d = cell[1][1];
  var t = mcnemarTest(a, b, c, d, conf);
  var F3 = { type: 'fixed', digits: 3 };
  var L0 = levels[0].label, L1 = levels[1].label;

  if (t.nd === 0) warnings.push('값이 바뀐 쌍이 하나도 없어 변화를 검정할 수 없습니다 (유의확률 1).');
  else if (t.nd < 10) warnings.push('값이 바뀐 쌍이 ' + t.nd + '건뿐이라 검정력이 매우 낮습니다. 정확 유의확률을 보세요.');

  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [{ label: '분석에 포함 (쌍)', n: n }, { label: '결측으로 제외', n: dataset.n - n }, { _kind: 'total', label: '합계', n: dataset.n }],
    notes: ['변수1 "' + mcnName_(v1) + '"(전)과 변수2 "' + mcnName_(v2) + '"(후)가 모두 유효한 쌍만 썼다. 값 "' + L1 + '"을(를) 양성(1), "' + L0 + '"을(를) 0으로 본다.']
  });
  blocks.push({
    type: 'table', caption: '교차표: 전(' + mcnName_(v1) + ') × 후(' + mcnName_(v2) + ')',
    columns: [{ key: 'row', label: '전 \\ 후' }, { key: 'c0', label: L0, format: { type: 'int' } }, { key: 'c1', label: L1, format: { type: 'int' } }, { key: 'total', label: '합계', format: { type: 'int' } }],
    rows: [
      { row: L0, c0: a, c1: b, total: a + b },
      { row: L1, c0: c, c1: d, total: c + d },
      { _kind: 'total', row: '합계', c0: a + c, c1: b + d, total: n }
    ],
    notes: ['검정에는 값이 바뀐 쌍(전 ' + L0 + ' → 후 ' + L1 + ' = ' + b + '건, 전 ' + L1 + ' → 후 ' + L0 + ' = ' + c + '건)만 쓴다. 안 바뀐 쌍은 검정에 영향을 주지 않는다.']
  });
  blocks.push({
    type: 'table', caption: 'McNemar 검정 통계량',
    columns: [
      { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'nd', label: '값이 바뀐 쌍', format: { type: 'int' } },
      { key: 'exactP', label: '정확 유의확률(양측)', format: { type: 'pvalue' } },
      { key: 'chi2', label: '카이제곱(연속성 수정)', format: F3 }, { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'chiP', label: '점근 유의확률', format: { type: 'pvalue' } }
    ],
    rows: [{ _bold: true, n: t.n, nd: t.nd, exactP: t.exactP, chi2: t.chi2, df: t.nd > 0 ? 1 : null, chiP: t.chiP }],
    notes: ['정확 유의확률은 값이 바뀐 쌍이 ' + L0 + '→' + L1 + ' 방향과 반대 방향으로 같은 확률(0.5)로 일어난다는 가정의 이항 검정이다. 값이 바뀐 쌍이 25건 미만이면 정확 유의확률을 보고, 많으면 점근 유의확률도 비슷하다.']
  });
  blocks.push({
    type: 'table', caption: '변화 요약',
    columns: [{ key: 'label', label: '' }, { key: 'value', label: '값', format: F3 }, { key: 'lo', label: conf + '% 하한', format: F3 }, { key: 'hi', label: conf + '% 상한', format: F3 }],
    rows: [
      { label: '"' + L1 + '" 비율 (전)', value: (c + d) / n, lo: null, hi: null },
      { label: '"' + L1 + '" 비율 (후)', value: (b + d) / n, lo: null, hi: null },
      { _bold: true, label: '비율 차이 (후 − 전)', value: t.diff, lo: t.diffLo, hi: t.diffHi },
      { label: '오즈비 (바뀐 쌍 중 ' + L0 + '→' + L1 + ' ÷ ' + L1 + '→' + L0 + ')', value: t.odds, lo: t.oddsLo, hi: t.oddsHi }
    ],
    notes: ['비율 차이의 구간은 Wald 근사, 오즈비의 구간은 정확(Clopper–Pearson) 구간이다. 오즈비가 1보다 크면 ' + L0 + '→' + L1 + ' 변화가 더 많다. 오즈비는 ' + L1 + '→' + L0 + ' 변화가 0건이면 정의되지 않는다.']
  });
  return { analysis: 'mcnemar', title: 'McNemar 검정', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = { mcnemarExactP: mcnemarExactP, mcnemarTest: mcnemarTest, runMcnemar: runMcnemar };
}
