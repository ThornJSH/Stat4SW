/**
 * 효과크기: Cohen d, Hedges g. 순수 함수. t검정 세 종류(독립표본, 대응표본, 일표본)가 결과 아래에 효과크기 표로 보여 준다.
 *
 * - 독립표본: d = (평균1 − 평균2) / 합동 표준편차, 자유도 n1 + n2 − 2, 유효 표본 n1·n2/(n1 + n2)
 * - 대응표본: d_z = 평균 차이 / 차이의 표준편차, 유효 표본 n, 자유도 n − 1.
 *   d_av = 평균 차이 / 두 측정 표준편차의 평균 (다른 연구와 비교하기 좋음, 신뢰구간은 없음)
 * - 일표본: d = (평균 − 검정값) / 표본 표준편차, 유효 표본 n, 자유도 n − 1
 * - Hedges 보정 J = 1 − 3/(4·df − 1): g = d·J
 * - 신뢰구간(정확): t = d·√n_eff 를 비중심 t 분포(비중심 모수 δ = δ_pop·√n_eff)의 관측값으로 보고
 *   P(T ≤ t; δ_L) = 1 − α/2, P(T ≤ t; δ_U) = α/2 인 δ 를 구해 d = δ/√n_eff. Hedges g 의 구간은 d 의 구간 × J.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function effDeps_() {
  if (typeof module !== 'undefined') return { pnct: require('./Distributions.js').pnct };
  return { pnct: pnct };
}

/** Hedges 보정 계수 */
function hedgesJ(df) { return 1 - 3 / (4 * df - 1); }

/** pnct(t, df, δ) = p 인 δ (δ가 커질수록 누적확률은 줄어든다) */
function effSolveNcp_(t, df, p) {
  var pn = effDeps_().pnct;
  var f = function (dl) { return pn(t, df, dl) - p; };
  var lo = t - 8, hi = t + 8, i;
  for (i = 0; i < 60 && f(lo) < 0; i++) lo -= 8;
  for (i = 0; i < 60 && f(hi) > 0; i++) hi += 8;
  for (i = 0; i < 200; i++) {
    var mid = (lo + hi) / 2;
    if (f(mid) > 0) lo = mid; else hi = mid;
    if (hi - lo < 1e-13 * (1 + Math.abs(mid))) break;
  }
  return (lo + hi) / 2;
}

/**
 * Cohen d 의 정확 신뢰구간. d: 관측 d, nEff: 유효 표본, df, conf: 신뢰수준(%)
 * 반환: { lo, hi, approx } (approx: |t| 가 커서(약 35 초과) 비중심 t 의 정규 근사를 쓴 구간)
 */
function effectCI(d, nEff, df, conf) {
  var t = d * Math.sqrt(nEff), a = (1 - conf / 100) / 2, s = Math.sqrt(nEff);
  return { lo: effSolveNcp_(t, df, 1 - a) / s, hi: effSolveNcp_(t, df, a) / s, approx: Math.abs(t) > 35 };
}

function effBuild_(d, nEff, df, conf, extra) {
  var ci = effectCI(d, nEff, df, conf), J = hedgesJ(df);
  var r = { d: d, dLo: ci.lo, dHi: ci.hi, J: J, g: d * J, gLo: ci.lo * J, gHi: ci.hi * J, t: d * Math.sqrt(nEff), df: df, nEff: nEff, approx: ci.approx };
  Object.keys(extra || {}).forEach(function (k) { r[k] = extra[k]; });
  return r;
}

/** 독립표본. s1, s2: { n, mean, variance }. 계산할 수 없으면(사례 2개 미만, 합동 분산 0) null */
function effectIndependent(s1, s2, conf) {
  if (s1.n < 2 || s2.n < 2) return null;
  var df = s1.n + s2.n - 2;
  var sp = Math.sqrt(((s1.n - 1) * s1.variance + (s2.n - 1) * s2.variance) / df);
  if (!(sp > 0)) return null;
  return effBuild_((s1.mean - s2.mean) / sp, s1.n * s2.n / (s1.n + s2.n), df, conf, { std: sp });
}

/** 대응표본. diff: { n, mean, sd } (변수1 − 변수2), sdX, sdY: 각 측정의 표준편차 */
function effectPaired(diff, sdX, sdY, conf) {
  if (diff.n < 2 || !(diff.sd > 0)) return null;
  var r = effBuild_(diff.mean / diff.sd, diff.n, diff.n - 1, conf, { std: diff.sd });
  var avg = (sdX + sdY) / 2;
  r.dav = avg > 0 ? diff.mean / avg : null;
  r.stdAv = avg;
  return r;
}

/** 일표본. s: { n, mean, sd } */
function effectOneSample(s, mu, conf) {
  if (s.n < 2 || !(s.sd > 0)) return null;
  return effBuild_((s.mean - mu) / s.sd, s.n, s.n - 1, conf, { std: s.sd });
}

/** Cohen(1988)의 참고 기준: |d| < 0.2 거의 없음, < 0.5 작음, < 0.8 중간, 그 이상 큼 */
function effectSizeWord(d) {
  var a = Math.abs(d);
  return a < 0.2 ? '거의 없음' : a < 0.5 ? '작음' : a < 0.8 ? '중간' : '큼';
}

if (typeof module !== 'undefined') {
  module.exports = { hedgesJ: hedgesJ, effectCI: effectCI, effectIndependent: effectIndependent, effectPaired: effectPaired, effectOneSample: effectOneSample, effectSizeWord: effectSizeWord };
}
