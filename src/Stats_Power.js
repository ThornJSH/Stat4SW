/**
 * 검정력·표본 크기 계산 (데이터를 쓰지 않는 분석). 순수 함수: (dataset 은 쓰지 않음, request) → 결과 객체.
 *
 * 검정 6가지의 검정력:
 *  - 독립표본 t: 비중심 t, ncp = d √(n1 n2/(n1 + n2)), 자유도 n1 + n2 − 2 (n2 = ⌈집단 크기 비 × n1⌉)
 *  - 대응표본·일표본 t: ncp = d √n, 자유도 n − 1
 *  - 일원 분산분석(균형): 비중심 F, 자유도 (k − 1, N − k), λ = f² N (N = k × 집단당 n)
 *  - 상관(피어슨): Fisher z 근사 (pingouin 과 같은 보정판). 임계 상관 r_c = √(t_c²/(t_c² + N − 2)) (t_c: 자유도 N − 2 의 임계값),
 *    z_r = atanh|r| + |r|/(2(N − 1)) (편향 보정), 검정력 = Φ((z_r − atanh r_c)√(N − 3)) (+ 양측이면 Φ((−z_r − atanh r_c)√(N − 3))). 정확 분포와 표본 크기가 1명 안팎으로 같다
 *  - 두 비율: Cohen 의 h = 2asin√p1 − 2asin√p2 의 정규 근사, δ = |h| √(n1 n2/(n1 + n2))
 *  - 선형 회귀(전체 R²): 비중심 F, 자유도 (u, N − u − 1), λ = f² N
 * 표본 크기는 목표 검정력 이상이 되는 가장 작은 정수, 탐지 가능한 효과크기는 이분법으로 구한다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function pwDeps_() {
  if (typeof module !== 'undefined') {
    var dist = require('./Distributions.js');
    return {
      pnorm: dist.pnorm, qnorm: dist.qnorm, pnct: dist.pnct, qt: dist.qt, qf: dist.qf, lgamma: dist.lgamma, ibetaPQ_: dist.ibetaPQ_,
      lineChartBlock: require('./Charts.js').lineChartBlock
    };
  }
  return { pnorm: pnorm, qnorm: qnorm, pnct: pnct, qt: qt, qf: qf, lgamma: lgamma, ibetaPQ_: ibetaPQ_, lineChartBlock: lineChartBlock };
}

var PW_MAX_N_ = 1000000;
var PW_TESTS_ = {
  ttest2: { label: '독립표본 t검정', es: 'd', unit: '집단당 N', conv: [0.2, 0.5, 0.8], medium: 0.5, hi: 10, minN: 2 },
  ttest1: { label: '대응표본·일표본 t검정', es: 'd', unit: '쌍(사례) 수 N', conv: [0.2, 0.5, 0.8], medium: 0.5, hi: 10, minN: 2 },
  anova: { label: '일원 분산분석', es: 'f', unit: '집단당 N', conv: [0.1, 0.25, 0.4], medium: 0.25, hi: 5, minN: 2 },
  corr: { label: '상관계수 (피어슨)', es: 'r', unit: '총 N', conv: [0.1, 0.3, 0.5], medium: 0.3, hi: 0.9999, minN: 4 },
  prop: { label: '두 비율 비교', es: 'h', unit: '집단당 N', conv: [0.2, 0.5, 0.8], medium: 0.5, hi: 3.14, minN: 2 },
  reg: { label: '선형 회귀 (전체 R²)', es: 'f²', unit: '총 N', conv: [0.02, 0.15, 0.35], medium: 0.15, hi: 100, minN: 3 }
};
var PW_SIZE_WORDS_ = ['작음', '중간', '큼'];

/* ---------- 비중심 F ---------- */

/** 비중심 F 분포의 누적확률 P(F' ≤ x): 포아송(λ/2) 가중의 중심 F(불완전 베타) 합. 모드에서 양쪽으로 퍼져 나가며 더한다. */
function pncf(x, d1, d2, lambda) {
  var d = pwDeps_();
  if (!(x > 0)) return 0;
  var y = d1 * x / (d1 * x + d2), yc = d2 / (d1 * x + d2);
  var half = lambda / 2;
  if (!(half > 0)) return d.ibetaPQ_(y, yc, d1 / 2, d2 / 2)[0];
  var logHalf = Math.log(half);
  var weight = function (j) { return Math.exp(-half + j * logHalf - d.lgamma(j + 1)); };
  var mode = Math.floor(half), sum = 0, j, w;
  for (j = mode; ; j++) {
    w = weight(j);
    sum += w * d.ibetaPQ_(y, yc, d1 / 2 + j, d2 / 2)[0];
    if (j > mode && w < 1e-18) break;
    if (j > mode + 100000) break;
  }
  for (j = mode - 1; j >= 0; j--) {
    w = weight(j);
    sum += w * d.ibetaPQ_(y, yc, d1 / 2 + j, d2 / 2)[0];
    if (w < 1e-18) break;
  }
  return Math.min(1, sum);
}

/* ---------- 검정력 ---------- */

/** opt: { alpha, tails (2|1), ratio, k, u }. n 은 검정마다 단위가 다르다 (PW_TESTS_[test].unit). 반환: 검정력 (0~1) */
function powerOf(test, n, es, opt) {
  var d = pwDeps_();
  var alpha = opt.alpha, two = opt.tails !== 1;
  var t, df, ncp, tc, n1, n2;
  if (test === 'ttest2' || test === 'prop') {
    n1 = n; n2 = Math.ceil(opt.ratio * n);
    var scale = Math.sqrt(n1 * n2 / (n1 + n2));
    if (test === 'prop') {
      var z = d.qnorm(two ? alpha / 2 : alpha, false), delta = Math.abs(es) * scale;
      return d.pnorm(delta - z, true) + (two ? d.pnorm(-delta - z, true) : 0);
    }
    df = n1 + n2 - 2; ncp = Math.abs(es) * scale;
  } else if (test === 'ttest1') {
    df = n - 1; ncp = Math.abs(es) * Math.sqrt(n);
  } else if (test === 'anova') {
    var N = opt.k * n, df1 = opt.k - 1, df2 = N - opt.k;
    var fc = d.qf(alpha, df1, df2, false);
    return 1 - pncf(fc, df1, df2, es * es * N);
  } else if (test === 'reg') {
    var df2r = n - opt.u - 1;
    var fcr = d.qf(alpha, opt.u, df2r, false);
    return 1 - pncf(fcr, opt.u, df2r, es * n);
  } else if (test === 'corr') {
    var rr = Math.min(Math.abs(es), 0.9999999), dof = n - 2;
    var ttt = d.qt(two ? alpha / 2 : alpha, dof, false);
    var zrc = Math.atanh(Math.sqrt(ttt * ttt / (ttt * ttt + dof))), zr = Math.atanh(rr) + rr / (2 * (n - 1)), sq = Math.sqrt(n - 3);
    return d.pnorm((zr - zrc) * sq, true) + (two ? d.pnorm((-zr - zrc) * sq, true) : 0);
  } else {
    throw new Error('알 수 없는 검정입니다: ' + test);
  }
  tc = d.qt(two ? alpha / 2 : alpha, df, false);
  var upper = 1 - d.pnct(tc, df, ncp);
  return two ? upper + d.pnct(-tc, df, ncp) : upper;
}

/** 목표 검정력 이상이 되는 가장 작은 N (단위는 검정마다). 한도를 넘으면 null */
function solveN(test, es, target, opt) {
  var lo = opt.minN || PW_TESTS_[test].minN, hi = lo;
  if (powerOf(test, lo, es, opt) >= target) return lo;
  while (powerOf(test, hi, es, opt) < target) {
    lo = hi;
    hi = hi * 2;
    if (hi > PW_MAX_N_) return null;
  }
  while (hi - lo > 1) {
    var mid = Math.floor((lo + hi) / 2);
    if (powerOf(test, mid, es, opt) >= target) hi = mid; else lo = mid;
  }
  return hi;
}

/** 주어진 N 과 목표 검정력에서 탐지할 수 있는 최소 효과크기 (이분법). 한도 안에서 달성할 수 없으면 null */
function solveEffect(test, n, target, opt) {
  var lo = 1e-6, hi = PW_TESTS_[test].hi;
  if (powerOf(test, n, hi, opt) < target) return null;
  for (var i = 0; i < 200; i++) {
    var mid = (lo + hi) / 2;
    if (powerOf(test, n, mid, opt) >= target) hi = mid; else lo = mid;
    if (hi - lo < 1e-12 * Math.max(1, hi)) break;
  }
  return hi;
}

/* ---------- 실행 ---------- */

function pwNum_(opt, key, label, def, min, max, inclusive) {
  var v = opt[key];
  if (v === undefined || v === null || String(v).trim() === '') return def;
  var x = Number(v);
  if (!isFinite(x)) throw new Error('"' + label + '"은(는) 숫자여야 합니다.');
  if (min !== null && (inclusive ? x < min : x <= min)) throw new Error('"' + label + '"은(는) ' + min + (inclusive ? ' 이상' : ' 초과') + '이어야 합니다.');
  if (max !== null && (inclusive ? x > max : x >= max)) throw new Error('"' + label + '"은(는) ' + max + (inclusive ? ' 이하' : ' 미만') + '이어야 합니다.');
  return x;
}

function pwF_(x, digits) { return Number(x).toFixed(digits === undefined ? 3 : digits); }
function pwCeil_(x) { return Math.ceil(x - 1e-9); }

function runPower(dataset, request) {
  var d = pwDeps_();
  var opt = request.options || {};
  var test = PW_TESTS_[opt.test] ? opt.test : 'ttest2';
  var info = PW_TESTS_[test];
  var solve = opt.solve === 'power' || opt.solve === 'effect' ? opt.solve : 'n';
  var alpha = pwNum_(opt, 'alpha', '유의수준', 0.05, 0, 1, false);
  var target = pwNum_(opt, 'power', '목표 검정력', 0.8, 0, 1, false);
  var tails = opt.tails === 'one' ? 1 : 2;
  var ratio = pwNum_(opt, 'ratio', '집단 크기 비', 1, 0, null, false);
  var dropout = pwNum_(opt, 'dropout', '예상 탈락률(%)', 0, 0, 100, true);
  if (dropout >= 100) throw new Error('"예상 탈락률(%)"은(는) 0 이상 100 미만이어야 합니다.');
  var k = pwNum_(opt, 'k', '집단 수', 3, 2, 50, true);
  var u = pwNum_(opt, 'u', '독립변수 수', 3, 1, 50, true);
  if (k !== Math.floor(k) || u !== Math.floor(u)) throw new Error('"집단 수"와 "독립변수 수"는 정수여야 합니다.');
  var nGiven = pwNum_(opt, 'n', '표본 크기', null, 0, null, false);
  if (nGiven !== null && nGiven !== Math.floor(nGiven)) throw new Error('"표본 크기"는 정수여야 합니다.');
  var popt = { alpha: alpha, tails: tails, ratio: ratio, k: k, u: u };
  var warnings = [];

  // 효과크기 (비율 검정은 두 비율에서 h 를 계산)
  var es = null, esNote = '', usedMedium = false;
  var p1 = pwNum_(opt, 'p1', '집단 1의 비율', null, 0, 1, false), p2 = pwNum_(opt, 'p2', '집단 2의 비율', null, 0, 1, false);
  var effIn = pwNum_(opt, 'effect', '효과크기', null, 0, null, false);
  if (test === 'prop' && p1 !== null && p2 !== null) {
    if (p1 === p2) throw new Error('두 비율이 같으면 차이가 없어 표본 크기를 구할 수 없습니다.');
    if (p1 >= 1 || p2 >= 1) throw new Error('비율은 0보다 크고 1보다 작아야 합니다.');
    es = Math.abs(2 * Math.asin(Math.sqrt(p1)) - 2 * Math.asin(Math.sqrt(p2)));
    esNote = '집단 1 ' + pwF_(p1, 3) + ', 집단 2 ' + pwF_(p2, 3) + ' → Cohen의 h = ' + pwF_(es);
  } else if (solve === 'effect') {
    if (test === 'prop') { /* h 로 풀고 두 비율은 입력하지 않는다 */ }
  } else if (effIn !== null) {
    if (test === 'corr' && effIn >= 1) throw new Error('상관계수 r 은 1보다 작아야 합니다.');
    es = effIn;
  } else if (test === 'prop') {
    throw new Error('두 비율 비교는 "집단 1의 비율"과 "집단 2의 비율"을 입력하거나 효과크기 h 를 입력하세요.');
  } else if (solve !== 'effect') {
    es = info.medium; usedMedium = true;
  }
  if (solve === 'effect' && nGiven === null) throw new Error('탐지 가능한 효과크기를 구하려면 "표본 크기"를 입력하세요.');
  if (solve === 'power' && nGiven === null) throw new Error('검정력을 구하려면 "표본 크기"를 입력하세요.');
  if (nGiven !== null && nGiven < info.minN) throw new Error('이 검정의 표본 크기는 ' + info.minN + ' 이상이어야 합니다.');
  if (test === 'reg' && nGiven !== null && nGiven < u + 2) throw new Error('회귀의 총 N 은 독립변수 수 + 2 이상이어야 합니다.');
  var minN = test === 'reg' ? u + 2 : info.minN;        // 회귀는 독립변수 수에 따라 달라진다
  popt.minN = minN;

  var unitName = info.unit;
  var twoGroup = test === 'ttest2' || test === 'prop';
  var totalOf = function (n) { return twoGroup ? n + Math.ceil(ratio * n) : test === 'anova' ? n * k : n; };
  var blocks = [];
  var inputRows = [
    { item: '검정', value: info.label },
    { item: '구하려는 값', value: solve === 'n' ? '필요한 표본 크기' : solve === 'power' ? '검정력' : '탐지 가능한 효과크기' },
    { item: '유의수준 (α)', value: pwF_(alpha, 3) + (tails === 1 ? ' (단측)' : ' (양측)') }
  ];
  if (solve !== 'power') inputRows.push({ item: '목표 검정력 (1 − β)', value: pwF_(target, 3) });
  if (es !== null) inputRows.push({ item: '효과크기 (' + info.es + ')', value: esNote || (pwF_(es) + (usedMedium ? ' (입력이 없어 중간 효과 사용)' : '')) });
  if (nGiven !== null && solve !== 'n') inputRows.push({ item: '표본 크기 (' + unitName + ')', value: String(nGiven) });
  if (test === 'anova') inputRows.push({ item: '집단 수', value: String(k) });
  if (test === 'reg') inputRows.push({ item: '독립변수 수', value: String(u) });
  if (twoGroup && ratio !== 1) inputRows.push({ item: '집단 크기 비 (집단 2 / 집단 1)', value: pwF_(ratio, 2) });
  if (dropout > 0) inputRows.push({ item: '예상 탈락률', value: pwF_(dropout, 1) + '%' });
  blocks.push({ type: 'table', caption: '입력 조건', columns: [{ key: 'item', label: '' }, { key: 'value', label: '값' }], rows: inputRows, notes: [] });

  var resRows = [];
  var nStar = null, shownN = null;
  if (solve === 'n') {
    nStar = solveN(test, es, target, popt);
    if (nStar === null) throw new Error('효과크기가 너무 작아 필요한 표본 크기가 ' + PW_MAX_N_ + '명을 넘습니다. 효과크기를 다시 확인하세요.');
    var ach = powerOf(test, nStar, es, popt);
    shownN = nStar;
    resRows.push({ item: '필요한 표본 크기 (' + unitName + ')', value: String(nStar) });
    if (twoGroup) resRows.push({ item: '집단 2의 N', value: String(Math.ceil(ratio * nStar)) });
    resRows.push({ item: '총 인원', value: String(totalOf(nStar)) });
    resRows.push({ item: '그 크기에서 실제 검정력', value: pwF_(ach, 4) });
    if (dropout > 0) {
      var rec = pwCeil_(nStar / (1 - dropout / 100));
      resRows.push({ item: '탈락을 감안한 모집 인원 (' + unitName + ')', value: String(rec) });
      resRows.push({ item: '탈락을 감안한 총 모집 인원', value: String(twoGroup ? rec + Math.ceil(ratio * rec) : test === 'anova' ? rec * k : rec) });
    }
    if (nStar > 2000) warnings.push('필요한 표본 크기가 ' + nStar + '로 매우 큽니다. 효과크기가 실제로 그만큼 작은지(작은 효과를 확인하려는 것인지) 다시 확인하세요.');
  } else if (solve === 'power') {
    shownN = nGiven;
    var pw = powerOf(test, nGiven, es, popt);
    resRows.push({ item: '검정력 (1 − β)', value: pwF_(pw, 4) });
    resRows.push({ item: '제2종 오류 확률 (β)', value: pwF_(1 - pw, 4) });
    resRows.push({ item: '총 인원', value: String(totalOf(nGiven)) });
    if (pw < 0.8) warnings.push('검정력이 ' + pwF_(pw, 3) + '로 흔한 기준(0.80)에 못 미칩니다. 효과가 실제로 있어도 놓칠 가능성이 큽니다.');
  } else {
    shownN = nGiven;
    es = solveEffect(test, nGiven, target, popt);
    if (es === null) throw new Error('이 표본 크기로는 목표 검정력 ' + target + '을 달성할 수 있는 효과크기가 없습니다. 표본 크기를 늘리세요.');
    var wi = es < (info.conv[0] + info.conv[1]) / 2 ? '작은 효과 이하' : es < (info.conv[1] + info.conv[2]) / 2 ? '중간 정도의 효과' : '큰 효과 이상';
    resRows.push({ item: '탐지 가능한 최소 효과크기 (' + info.es + ')', value: pwF_(es, 4) });
    resRows.push({ item: '크기 해석 (Cohen의 관행)', value: wi + ' (작음 ' + info.conv[0] + ', 중간 ' + info.conv[1] + ', 큼 ' + info.conv[2] + ')' });
    resRows.push({ item: '총 인원', value: String(totalOf(nGiven)) });
  }
  blocks.push({
    type: 'table', caption: '결과', columns: [{ key: 'item', label: '' }, { key: 'value', label: '값' }], rows: resRows,
    notes: [solve === 'n' ? '표본 크기는 목표 검정력 이상이 되는 가장 작은 정수이다 (그보다 1 작으면 검정력이 목표에 못 미친다).' : '검정력은 효과가 실제로 있을 때 그 효과를 유의하게 찾아낼 확률이다. 흔히 0.80 이상을 목표로 한다.']
      .concat(test === 'corr' ? ['상관계수는 Fisher z 변환의 근사(임계 상관과 편향 보정 포함)로 구했다. 정확한 분포를 쓰는 G*Power와는 필요한 표본 크기가 1명 안팎으로 다를 수 있다.'] : [])
      .concat(test === 'prop' ? ['두 비율은 Cohen의 h(arcsine 변환)의 정규 근사로 구했다. 비율이 0이나 1에 매우 가까우면 정확 검정과 차이가 날 수 있다.'] : [])
      .concat(test === 'anova' ? ['집단 크기가 모두 같은 균형 설계를 가정한다. 총 N = 집단 수 × 집단당 N.'] : [])
      .concat(test === 'reg' ? ['전체 R² 의 F검정(독립변수 ' + u + '개 모두 효과 없음)의 검정력이다. 효과크기 f² = R²/(1 − R²). 총 N 은 모든 변수를 포함한 사례 수이다.'] : [])
  });

  // 효과크기 관행별 표
  var cRows = [];
  if (solve === 'power') {
    info.conv.forEach(function (v, q) { cRows.push({ lab: PW_SIZE_WORDS_[q], es: v, val: pwF_(powerOf(test, nGiven, v, popt), 3) }); });
  } else {
    info.conv.forEach(function (v, q) {
      var nn = solveN(test, v, target, popt);
      cRows.push({ lab: PW_SIZE_WORDS_[q], es: v, val: nn === null ? '한도 초과' : String(nn) + (twoGroup ? ' (총 ' + totalOf(nn) + ')' : test === 'anova' ? ' (총 ' + totalOf(nn) + ')' : '') });
    });
  }
  blocks.push({
    type: 'table', caption: solve === 'power' ? '효과크기별 검정력 (표본 크기 ' + nGiven + ')' : '효과크기별 필요 표본 크기 (검정력 ' + pwF_(target, 2) + ')',
    columns: [{ key: 'lab', label: '효과 크기 (Cohen의 관행)' }, { key: 'es', label: info.es, format: { type: 'fixed', digits: 2 } }, { key: 'val', label: solve === 'power' ? '검정력' : '필요한 ' + unitName }],
    rows: cRows,
    notes: ['효과크기를 모를 때는 선행연구의 값이나 실제로 의미 있다고 볼 최소한의 차이를 쓰는 것이 좋고, 관행(작음·중간·큼)은 참고용이다.']
  });

  // 검정력 곡선
  var center = shownN, loN = Math.max(minN, Math.floor(center * 0.4)), hiN = Math.max(loN + 4, Math.ceil(center * 1.8));
  var step = Math.max(1, Math.ceil((hiN - loN) / 40));
  var curveEs = es !== null ? es : info.medium;
  var pts = [];
  for (var nn2 = loN; nn2 <= hiN; nn2 += step) pts.push([nn2, powerOf(test, nn2, curveEs, popt)]);
  blocks.push(d.lineChartBlock('표본 크기에 따른 검정력 (' + info.es + ' = ' + pwF_(curveEs, 3) + ', α = ' + alpha + ')', unitName, '검정력', [
    { label: '검정력', points: pts }, { label: '목표 검정력', points: [[pts[0][0], target], [pts[pts.length - 1][0], target]] }
  ]));

  return { analysis: 'power', title: '검정력·표본 크기 계산 (' + info.label + ')', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = { pncf: pncf, powerOf: powerOf, solveN: solveN, solveEffect: solveEffect, runPower: runPower };
}
