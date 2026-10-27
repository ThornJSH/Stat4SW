/**
 * 이항 검정: 한 변수를 성공/실패로 나눠 성공 확률이 기준 비율(검정 비율)과 다른지 정확 검정한다.
 * 순수 함수: (dataset, request) → 결과 객체. 변수마다 결과 한 줄 (변수별로 결측을 따로 제외).
 *
 * - 성공: 변수 값이 옵션 "성공으로 볼 값"과 같은 사례, 나머지는 모두 실패 (숫자는 수치로, 글자는 문자열로 비교).
 * - 유의확률: X ~ Binomial(n, p0). 단측은 P(X ≥ x)(성공 비율이 크다) 또는 P(X ≤ x)(작다).
 *   양측은 관측값 x 보다 확률이 같거나 작은 모든 결과의 확률 합 (R binom.test, scipy binomtest 와 같다; 상대 허용 1e-7).
 *   p0 = 0.5 이면 단측 p 의 2배(최대 1)와 같다.
 * - 꼬리 확률은 정규화된 불완전 베타로 계산한다: P(X ≤ k) = I_{1−p}(n−k, k+1), P(X ≥ k) = I_p(k, n−k+1).
 * - 신뢰구간: Clopper–Pearson(정확, 베타 분위수)과 Wilson 점수 구간.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function bnDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return { getVars: data.getVars, lgamma: dist.lgamma, ibetaPQ_: dist.ibetaPQ_, qnorm: dist.qnorm };
  }
  return { getVars: getVars, lgamma: lgamma, ibetaPQ_: ibetaPQ_, qnorm: qnorm };
}

var BN_MAX_VARS_ = 20;
var BN_F3_ = { type: 'fixed', digits: 3 };
function bnName_(v) { return v.label ? v.label : v.name; }

/** 이항 확률 P(X = i) 의 로그 */
function bnLogPmf_(i, n, p) {
  var d = bnDeps_();
  return d.lgamma(n + 1) - d.lgamma(i + 1) - d.lgamma(n - i + 1) + i * Math.log(p) + (n - i) * Math.log(1 - p);
}

/** P(X ≤ k) */
function bnCdf(k, n, p) {
  if (k < 0) return 0;
  if (k >= n) return 1;
  return bnDeps_().ibetaPQ_(1 - p, p, n - k, k + 1)[0];
}

/** P(X ≥ k) */
function bnSf(k, n, p) {
  if (k <= 0) return 1;
  if (k > n) return 0;
  return bnDeps_().ibetaPQ_(p, 1 - p, k, n - k + 1)[0];
}

/**
 * 정확 유의확률. alt: 'two' | 'greater' | 'less'
 */
function binomTest(x, n, p0, alt) {
  if (alt === 'greater') return bnSf(x, n, p0);
  if (alt === 'less') return bnCdf(x, n, p0);
  var mean = p0 * n;
  if (x === mean) return 1;
  var logD = bnLogPmf_(x, n, p0) + Math.log1p(1e-7);   // 관측값보다 확률이 같거나 작은 결과의 경계 (상대 허용 1e-7)
  var i, lo, hi, mid;
  if (x < mean) {
    // 최빈값 오른쪽에서 확률이 경계 이하가 되는 가장 작은 값 b (확률은 오른쪽으로 갈수록 줄어든다)
    lo = Math.ceil(mean); hi = n + 1;
    while (lo < hi) {
      mid = Math.floor((lo + hi) / 2);
      if (bnLogPmf_(mid, n, p0) <= logD) hi = mid; else lo = mid + 1;
    }
    return Math.min(1, bnCdf(x, n, p0) + bnSf(lo, n, p0));
  }
  // 최빈값 왼쪽에서 확률이 경계 이하인 가장 큰 값 a (왼쪽으로 갈수록 줄어든다)
  lo = -1; hi = Math.floor(mean);
  while (lo < hi) {
    mid = Math.floor((lo + hi + 1) / 2);
    if (bnLogPmf_(mid, n, p0) <= logD) lo = mid; else hi = mid - 1;
  }
  i = lo;
  return Math.min(1, bnCdf(i, n, p0) + bnSf(x, n, p0));
}

/** 베타(a, b) 분포의 분위수 (이분법) */
function bnBetaQuantile_(prob, a, b) {
  var d = bnDeps_(), lo = 0, hi = 1;
  for (var it = 0; it < 100; it++) {
    var mid = (lo + hi) / 2;
    if (d.ibetaPQ_(mid, 1 - mid, a, b)[0] < prob) lo = mid; else hi = mid;
    if (hi - lo < 1e-16) break;
  }
  return (lo + hi) / 2;
}

/** Clopper–Pearson 정확 구간. conf: 신뢰수준(%) */
function binomExactCi(x, n, conf) {
  var alpha = 1 - conf / 100;
  return [x === 0 ? 0 : bnBetaQuantile_(alpha / 2, x, n - x + 1), x === n ? 1 : bnBetaQuantile_(1 - alpha / 2, x + 1, n - x)];
}

/** Wilson 점수 구간 */
function binomWilsonCi(x, n, conf) {
  var z = bnDeps_().qnorm(1 - (1 - conf / 100) / 2), ph = x / n, z2 = z * z;
  var center = (x + z2 / 2) / (n + z2);
  var half = z * Math.sqrt(n) / (n + z2) * Math.sqrt(ph * (1 - ph) + z2 / (4 * n));
  return [Math.max(0, center - half), Math.min(1, center + half)];
}

function bnIsSuccess_(v, text) {
  var num = Number(text);
  if (typeof v === 'number' && text !== '' && isFinite(num)) return v === num;
  return String(v).trim() === String(text).trim();
}

var BN_TAILS_ = {
  two: { label: '양측', col: '정확 유의확률(양측)' },
  greater: { label: '단측: 성공 비율이 검정 비율보다 크다', col: '정확 유의확률(단측: 크다)' },
  less: { label: '단측: 성공 비율이 검정 비율보다 작다', col: '정확 유의확률(단측: 작다)' }
};

function runBinomial(dataset, request) {
  var d = bnDeps_();
  var vars = d.getVars(dataset, request.vars);
  if (vars.length > BN_MAX_VARS_) throw new Error('변수는 ' + BN_MAX_VARS_ + '개까지 선택할 수 있습니다.');
  var opt = request.options || {};
  var succ = opt.success === undefined || opt.success === null || String(opt.success).trim() === '' ? '1' : String(opt.success).trim();
  var p0 = opt.p0 === undefined || opt.p0 === null || String(opt.p0).trim() === '' ? 0.5 : Number(opt.p0);
  if (!isFinite(p0) || p0 <= 0 || p0 >= 1) throw new Error('"검정 비율"은(는) 0보다 크고 1보다 작은 숫자여야 합니다 (예: 0.5).');
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  var tails = BN_TAILS_[opt.tails] ? opt.tails : 'two';
  var warnings = [], freqRows = [], testRows = [], ciRows = [];

  vars.forEach(function (v) {
    var col = dataset.data[v.name], x = 0, n = 0, seen = {}, order = [], i;
    for (i = 0; i < dataset.n; i++) {
      if (col[i] === null) continue;
      n++;
      var key = String(col[i]);
      if (!(key in seen)) { seen[key] = 0; order.push(col[i]); }
      seen[key]++;
      if (bnIsSuccess_(col[i], succ)) x++;
    }
    var name = bnName_(v);
    if (n === 0) throw new Error('변수 "' + name + '"에 값이 하나도 없습니다.');
    if (x === 0) {
      throw new Error('변수 "' + name + '"에 성공으로 볼 값 "' + succ + '"이(가) 없습니다. 변수에 있는 값: ' + order.slice(0, 10).map(String).join(', ') + (order.length > 10 ? ', …' : '') + '. 성공으로 볼 값을 다시 입력하세요.');
    }
    if (order.length > 2) warnings.push('"' + name + '": 값이 ' + order.length + '가지 있습니다. "' + succ + '"인 사례만 성공이고 나머지는 모두 실패로 처리했습니다.');
    var lab = v.valueLabels[succ];
    var ci = binomExactCi(x, n, conf), wi = binomWilsonCi(x, n, conf);
    freqRows.push({ name: name, value: succ + (lab !== undefined ? ' (' + lab + ')' : ''), x: x, y: n - x, n: n, miss: dataset.n - n });
    testRows.push({ name: name, n: n, x: x, ph: x / n, p0: p0, p: binomTest(x, n, p0, tails) });
    ciRows.push({ name: name, ph: x / n, lo: ci[0], hi: ci[1], wlo: wi[0], whi: wi[1] });
  });

  var blocks = [];
  blocks.push({
    type: 'table', caption: '성공·실패 빈도',
    columns: [
      { key: 'name', label: '변수' }, { key: 'value', label: '성공으로 본 값' },
      { key: 'x', label: '성공', format: { type: 'int' } }, { key: 'y', label: '실패', format: { type: 'int' } },
      { key: 'n', label: '유효 N', format: { type: 'int' } }, { key: 'miss', label: '결측(제외)', format: { type: 'int' } }
    ],
    rows: freqRows,
    notes: ['변수마다 값이 있는 사례만 쓴다. 성공으로 본 값과 같은 사례가 성공이고 나머지는 모두 실패이다.']
  });
  blocks.push({
    type: 'table', caption: '이항 검정 (귀무가설: 성공 확률 = 검정 비율, ' + BN_TAILS_[tails].label + ')',
    columns: [
      { key: 'name', label: '변수' }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'x', label: '성공 수', format: { type: 'int' } },
      { key: 'ph', label: '관측 비율', format: BN_F3_ }, { key: 'p0', label: '검정 비율', format: BN_F3_ },
      { key: 'p', label: BN_TAILS_[tails].col, format: { type: 'pvalue' } }
    ],
    rows: testRows,
    notes: [tails === 'two'
      ? '양측 유의확률은 관측한 성공 수보다 일어날 확률이 같거나 작은 모든 결과의 확률 합이다 (R의 binom.test와 같다). 검정 비율이 0.5이면 단측 유의확률의 2배와 같다.'
      : '단측 유의확률은 성공 수가 관측값 이상(크다) 또는 이하(작다)일 정확한 이항 확률이다. 방향은 자료를 보기 전에 정해야 한다.',
      '유의확률이 작으면 성공 확률이 검정 비율과 다르다고 본다. 정규 근사가 아닌 정확한 이항 분포로 계산하므로 표본이 작아도 쓸 수 있다.']
  });
  blocks.push({
    type: 'table', caption: '성공 비율의 ' + conf + '% 신뢰구간',
    columns: [
      { key: 'name', label: '변수' }, { key: 'ph', label: '관측 비율', format: BN_F3_ },
      { key: 'lo', label: '정확(Clopper–Pearson) 하한', format: BN_F3_ }, { key: 'hi', label: '정확(Clopper–Pearson) 상한', format: BN_F3_ },
      { key: 'wlo', label: 'Wilson 하한', format: BN_F3_ }, { key: 'whi', label: 'Wilson 상한', format: BN_F3_ }
    ],
    rows: ciRows,
    notes: ['Clopper–Pearson 구간은 정확하지만 다소 넓고(보수적), Wilson 구간은 표본이 작거나 비율이 0·1에 가까워도 안정적이다. 검정 비율이 구간 밖이면 같은 신뢰수준에서 유의하다는 뜻과 대체로 일치한다.']
  });
  return { analysis: 'binom', title: '이항 검정', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = { binomTest: binomTest, binomExactCi: binomExactCi, binomWilsonCi: binomWilsonCi, bnCdf: bnCdf, bnSf: bnSf, runBinomial: runBinomial };
}
