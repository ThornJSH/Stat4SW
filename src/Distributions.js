/**
 * 확률분포: 정규, t, F, χ², Studentized range(Tukey)의 누적확률(p*)과 역함수(q*), 그리고 p값 도우미.
 * 순수 함수만 둔다. 외부 라이브러리 없이 직접 구현한다.
 *
 * 인자 규칙 (R/scipy와 같은 방식)
 *   p*(x, ..., lowerTail=true)  : lowerTail=false 이면 오른쪽(위) 꼬리 확률 P(X > x)
 *   q*(p, ..., lowerTail=true)  : lowerTail=false 이면 p를 오른쪽 꼬리 확률로 해석
 * 꼬리 확률은 1에서 빼서 구하지 않고 해당 꼬리를 직접 계산하므로 p값이 매우 작아도 상대 정확도를 유지한다.
 *
 * 정확도 (scipy 1.17 대비, tests/distributions.test.js)
 *   - 정규·χ²·t·F 의 p값: 상대오차 1e-10 이내. 꼬리 확률이 1e-300 근처여도 상대 정확도를 유지한다.
 *   - t, F 는 자유도가 아주 크면 오차가 df에 비례해 커진다 (≈ 3e-17 × df).
 *     df ≤ 1e6 에서 3e-11 이하, df = 1e9 에서 약 2e-8. 이 도구의 사례 수(수천 건)에서는 영향이 없다.
 *   - 역함수: 상대오차 1e-9 이내.
 *   - Studentized range (ptukey, qtukey): 누적확률 절대오차 1e-10 이내(scipy 대조), 위꼬리는 1 − P 라서
 *     p값이 1e-10보다 작으면 상대 정확도는 없다. 한 번 계산에 수 ms(qtukey는 수십 ms)가 든다.
 *   - 한계: |t| 가 1e154 근처를 넘는 극단값(꼬리 확률 1e-300 미만이면서 자유도 2 이하)은 다루지 않는다.
 *
 * 구현 요약
 *   - 불완전 감마(χ², 정규)·불완전 베타(t, F): 급수/연분수 + Loader의 saddle-point 접두 계수
 *     (자유도가 매우 커도 지수부 상쇄로 인한 정밀도 손실이 없다)
 *   - 정규 역함수: Acklam 근사 + Halley 보정
 *   - t, F, χ² 역함수: 꼬리 확률의 로그에 대한 Newton + 구간 이분법
 */

var LN_SQRT_2PI_ = 0.918938533204672741780329736406; // ln(sqrt(2π))
var LN_2PI_ = 1.837877066409345483560659472811;      // ln(2π)
var CF_EPS_ = 1e-15;
var CF_FPMIN_ = 1e-300;

/* ---------- 기초 함수 ---------- */

/** ln Γ(x), x > 0. 15 이상으로 올린 뒤 Stirling 급수를 쓴다. */
function lgamma(x) {
  if (isNaN(x) || x < 0) return NaN;
  if (x === 0) return Infinity;
  if (x === Infinity) return Infinity;
  var shift = 0;
  while (x < 15) { shift += Math.log(x); x += 1; }
  var x2 = 1 / (x * x);
  var series = (1 / 12 - x2 * (1 / 360 - x2 * (1 / 1260 - x2 * (1 / 1680 - x2 * (1 / 1188 - x2 * (691 / 360360 - x2 / 156)))))) / x;
  return (x - 0.5) * Math.log(x) - x + LN_SQRT_2PI_ + series - shift;
}

/** ln Γ(n+1) - [(n+½)ln n - n + ½ln(2π)] : Stirling 급수의 오차항 */
function stirlerr_(n) {
  if (n <= 15) return lgamma(n + 1) - (n + 0.5) * Math.log(n) + n - LN_SQRT_2PI_;
  var n2 = 1 / (n * n);
  return (1 / 12 - n2 * (1 / 360 - n2 * (1 / 1260 - n2 * (1 / 1680 - n2 * (1 / 1188 - n2 * (691 / 360360 - n2 / 156)))))) / n;
}

/** x ln(x/np) + np - x 를 상쇄 오차 없이 계산 */
function bd0_(x, np) {
  if (Math.abs(x - np) < 0.1 * (x + np)) {
    var v = (x - np) / (x + np);
    var s = (x - np) * v;
    if (Math.abs(s) < 2.2250738585072014e-308) return s;
    var ej = 2 * x * v;
    v = v * v;
    for (var j = 1; j < 1000; j++) {
      ej *= v;
      var s1 = s + ej / (2 * j + 1);
      if (s1 === s) return s1;
      s = s1;
    }
    return s;
  }
  return x * Math.log(x / np) + np - x;
}

/** λ^x e^-λ / Γ(x+1)  (x ≥ 0 실수) */
function dpoisRaw_(x, lambda) {
  if (x === 0) return Math.exp(-lambda);
  if (lambda === 0) return 0;
  return Math.exp(-stirlerr_(x) - bd0_(x, lambda)) / Math.sqrt(2 * Math.PI * x);
}

/** Γ(x+y+1)/(Γ(x+1)Γ(y+1)) p^x q^y  (x, y > 0 실수) */
function dbinomRaw_(x, y, p, q) {
  var n = x + y;
  // n·p + n·q = n 이 정확히 성립해야 한다. 어긋난 만큼(≈ n·1e-16)이 지수부 오차로 남는다.
  var np = n * p, nq = n * q;
  if (p < q) nq = n - np; else np = n - nq;
  var lc = stirlerr_(n) - stirlerr_(x) - stirlerr_(y) - bd0_(x, np) - bd0_(y, nq);
  var lf = LN_2PI_ + Math.log(x) + Math.log(y) - Math.log(n);
  return Math.exp(lc - 0.5 * lf);
}

/* ---------- 불완전 감마 / 베타 ---------- */

/** [P(a,x), Q(a,x)] 정규화된 불완전 감마. 작은 쪽 꼬리를 직접 계산한다. */
function gammaPQ_(a, x) {
  if (x <= 0) return [0, 1];
  if (x === Infinity) return [1, 0];
  // x^a e^-x / Γ(a+1). a = 0.5 (정규분포)는 Γ(1.5) = √π/2 로 닫힌 형태라 바로 계산한다.
  var pref = a === 0.5 ? 2 * Math.sqrt(x / Math.PI) * Math.exp(-x) : dpoisRaw_(a, x);
  if (x < a + 1) {
    var sum = 1, term = 1;
    for (var k = 1; k < 1000000; k++) {
      term *= x / (a + k);
      sum += term;
      if (term < sum * 1e-17) break;
    }
    var p = Math.min(1, pref * sum);
    return [p, 1 - p];
  }
  var b = x + 1 - a, c = 1 / CF_FPMIN_, d = 1 / b, h = d;
  for (var i = 1; i < 1000000; i++) {
    var an = -i * (i - a);
    b += 2;
    d = an * d + b; if (Math.abs(d) < CF_FPMIN_) d = CF_FPMIN_;
    c = b + an / c; if (Math.abs(c) < CF_FPMIN_) c = CF_FPMIN_;
    d = 1 / d;
    var del = d * c;
    h *= del;
    if (Math.abs(del - 1) < CF_EPS_) break;
  }
  var q = Math.min(1, a * pref * h);
  return [1 - q, q];
}

/**
 * 불완전 베타의 연분수 (Lentz). y = 1 - x 를 따로 받는다.
 * 초기값 d = 1 - (a+b)x/(a+1) 은 x가 1에 가까우면 큰 수끼리 빼서 정밀도를 잃으므로,
 * 그 경우 같은 값을 y로 표현한 ((1-b) + (a+b)y)/(a+1) 로 계산한다.
 */
function betacf_(x, y, a, b) {
  var qab = a + b, qap = a + 1, qam = a - 1;
  var c = 1, d = x > 0.5 ? ((1 - b) + qab * y) / qap : 1 - qab * x / qap;
  if (Math.abs(d) < CF_FPMIN_) d = CF_FPMIN_;
  d = 1 / d;
  var h = d;
  for (var m = 1; m < 100000; m++) {
    var m2 = 2 * m;
    var aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < CF_FPMIN_) d = CF_FPMIN_;
    c = 1 + aa / c; if (Math.abs(c) < CF_FPMIN_) c = CF_FPMIN_;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < CF_FPMIN_) d = CF_FPMIN_;
    c = 1 + aa / c; if (Math.abs(c) < CF_FPMIN_) c = CF_FPMIN_;
    d = 1 / d;
    var del = d * c;
    h *= del;
    if (Math.abs(del - 1) < CF_EPS_) break;
  }
  return h;
}

/**
 * [I_x(a,b), 1 - I_x(a,b)] 정규화된 불완전 베타. y = 1 - x 를 따로 받아 x가 1에 가까울 때의 상쇄를 피한다.
 */
function ibetaPQ_(x, y, a, b) {
  if (x <= 0) return [0, 1];
  if (y <= 0) return [1, 0];
  var pre = dbinomRaw_(a, b, x, y); // x^a y^b Γ(a+b+1)/(Γ(a+1)Γ(b+1))
  var n = a + b;
  if (x < (a + 1) / (n + 2)) {
    var p = Math.min(1, pre * (b / n) * betacf_(x, y, a, b));
    return [p, 1 - p];
  }
  var q = Math.min(1, pre * (a / n) * betacf_(y, x, b, a));
  return [1 - q, q];
}

/* ---------- 정규분포 ---------- */

function pnorm(x, lowerTail) {
  if (isNaN(x)) return NaN;
  if (lowerTail === false) x = -x; // P(Z > x) = P(Z < -x)
  if (x === 0) return 0.5;
  var pq = gammaPQ_(0.5, x * x / 2); // p = P(|Z| < |x|), q = P(|Z| > |x|)
  return x < 0 ? 0.5 * pq[1] : 0.5 + 0.5 * pq[0];
}

/** (p, lowerTail) → 작은 쪽 꼬리 확률 u(≤0.5), 그 꼬리가 아래쪽인지, 분위수의 부호 */
function tailSelect_(p, lowerTail) {
  var lower = lowerTail !== false;
  if (lower) return p <= 0.5 ? { u: p, useLower: true, sign: 1 } : { u: 1 - p, useLower: false, sign: -1 };
  return p <= 0.5 ? { u: p, useLower: false, sign: -1 } : { u: 1 - p, useLower: true, sign: 1 };
}

/** 아래꼬리 확률이 u (0 < u ≤ 0.5)인 표준정규 분위수 (≤ 0) */
function qnormLowerSmall_(u) {
  var a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  var b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  var c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  var d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  var x, q;
  if (u < 0.02425) {
    q = Math.sqrt(-2 * Math.log(u));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else {
    q = u - 0.5;
    var r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }
  for (var i = 0; i < 2; i++) { // Halley 보정
    var e = pnorm(x) - u;
    var w = e * Math.sqrt(2 * Math.PI) * Math.exp(x * x / 2);
    x = x - w / (1 + x * w / 2);
  }
  return x;
}

function qnorm(p, lowerTail) {
  if (isNaN(p) || p < 0 || p > 1) return NaN;
  var lower = lowerTail !== false;
  if (p === 0) return lower ? -Infinity : Infinity;
  if (p === 1) return lower ? Infinity : -Infinity;
  var s = tailSelect_(p, lowerTail);
  return s.u === 0.5 ? 0 : s.sign * qnormLowerSmall_(s.u);
}

/* ---------- 역함수 공통 풀이 ---------- */

/**
 * 양의 지지집합을 가진 분포에서 꼬리 확률이 u가 되는 x를 찾는다.
 * v = ln x 에 대한 ln(꼬리)의 Newton이며, 벗어나면 구간 이분법으로 대체한다.
 * (꼬리가 x의 거듭제곱처럼 줄어드는 영역에서 거의 한 번에 수렴한다)
 * tailFn: x → 꼬리 확률, pdfFn: x → 밀도, useLower: tailFn이 아래꼬리(증가)인지
 */
function solveTail_(u, useLower, tailFn, pdfFn, x0) {
  var sgn = useLower ? 1 : -1;
  var lu = Math.log(u);
  var VMAX = 709, VMIN = -745;
  // sgn·g(v) < 0 이면 x = e^v 가 해보다 작다
  function g(v) {
    var t = tailFn(Math.exp(v));
    return t > 0 ? Math.log(t) - lu : -Infinity;
  }
  var v0 = Math.log(x0 > 0 && isFinite(x0) ? x0 : 1);
  var lo = v0, hi = v0, step;
  if (sgn * g(v0) < 0) {
    for (step = 1; sgn * g(hi) < 0; step *= 2) {
      lo = hi; hi += step;
      if (hi > VMAX) { // 경계에서 한 번 더 확인한다 (구간을 건너뛰어 놓치지 않도록)
        hi = VMAX;
        if (sgn * g(hi) < 0) return Infinity;
        break;
      }
    }
  } else {
    for (step = 1; sgn * g(lo) > 0; step *= 2) {
      hi = lo; lo -= step;
      if (lo < VMIN) {
        lo = VMIN;
        if (sgn * g(lo) > 0) return 0;
        break;
      }
    }
  }
  var v = lo < v0 && v0 < hi ? v0 : (lo + hi) / 2;
  for (var i = 0; i < 200; i++) {
    var gv = g(v);
    var s = sgn * gv;
    if (s === 0) break;
    if (s < 0) lo = v; else hi = v;
    var next = NaN;
    if (isFinite(gv)) {
      var x = Math.exp(v);
      next = v - gv / (sgn * x * pdfFn(x) / tailFn(x));
    }
    if (!(next > lo && next < hi)) next = (lo + hi) / 2;
    var done = Math.abs(next - v) <= 4e-16 * Math.max(1, Math.abs(next));
    v = next;
    if (done) break;
  }
  return Math.exp(v);
}

/* ---------- t 분포 ---------- */

function pt(t, df, lowerTail) {
  if (isNaN(t) || isNaN(df) || df <= 0) return NaN;
  var lower = lowerTail !== false;
  if (!isFinite(df)) return pnorm(t, lowerTail);
  if (!isFinite(t)) return (t > 0) === lower ? 1 : 0;
  var t2 = t * t, x, y;
  if (isFinite(t2)) {
    x = df / (df + t2);
    y = t2 / (df + t2);
  } else { // |t| > 1e154: t²의 오버플로를 피한다 (그 이상의 극단값은 0 또는 1로 수렴)
    x = 1 / (1 + t * (t / df));
    y = 1 / (1 + (df / t) / t);
  }
  var pq = ibetaPQ_(x, y, df / 2, 0.5); // p = P(|T| > |t|)
  var far = 0.5 * pq[0];         // |t| 바깥쪽 꼬리 P(T > |t|) = P(T < -|t|)
  var near = 0.5 + 0.5 * pq[1];  // 나머지 쪽 P(T < |t|) = P(T > -|t|)
  if (t >= 0) return lower ? near : far;
  return lower ? far : near;
}

function dtPdf_(x, df) {
  return Math.exp(lgamma((df + 1) / 2) - lgamma(df / 2) - 0.5 * Math.log(df * Math.PI) - (df + 1) / 2 * Math.log(1 + x * x / df));
}

function qt(p, df, lowerTail) {
  if (isNaN(p) || isNaN(df) || p < 0 || p > 1 || df <= 0) return NaN;
  var lower = lowerTail !== false;
  if (p === 0) return lower ? -Infinity : Infinity;
  if (p === 1) return lower ? Infinity : -Infinity;
  if (!isFinite(df)) return qnorm(p, lowerTail);
  var s = tailSelect_(p, lowerTail);
  if (s.u === 0.5) return 0;
  var z = -qnormLowerSmall_(s.u);
  var x0 = z + (z * z * z + z) / (4 * df);
  var x = solveTail_(s.u, false, function (v) { return pt(v, df, false); }, function (v) { return dtPdf_(v, df); }, x0);
  return -s.sign * x;
}

/* ---------- 비중심 t 분포 ---------- */

/**
 * 비중심 t 분포의 누적확률 P(T ≤ t), 자유도 df, 비중심 모수 ncp (Lenth 1989, AS 243의 급수; R의 pnt 와 같은 방식).
 * 효과크기(Cohen d)의 정확 신뢰구간에 쓴다. 자유도가 4e5를 넘거나 |ncp|가 약 37.6을 넘으면 정규 근사로 바꾼다.
 */
function pnct(t, df, ncp) {
  if (isNaN(t) || isNaN(df) || isNaN(ncp) || df <= 0) return NaN;
  if (!isFinite(t)) return t > 0 ? 1 : 0;
  var negdel = false, tt = t, del = ncp;
  if (t < 0) { negdel = true; tt = -t; del = -ncp; }
  var tnc;
  if (df > 4e5 || del * del > 2 * Math.LN2 * 1021) {
    var s0 = 1 / (4 * df);
    var z = (tt * (1 - s0) - del) / Math.sqrt(1 + tt * tt * 2 * s0);
    tnc = pnorm(z, true);
    return negdel ? 1 - tnc : tnc;
  }
  var x = tt * tt;
  x = x / (x + df);
  if (x > 0) {
    var lambda = del * del;
    var p = 0.5 * Math.exp(-0.5 * lambda);
    var q = Math.sqrt(2 / Math.PI) * p * del;
    var s = 0.5 - p;
    if (s < 1e-7) s = -0.5 * (Math.exp(-0.5 * lambda) - 1);
    var a = 0.5, b = 0.5 * df;
    var rxb = Math.pow(1 - x, b);
    var albeta = 0.5 * Math.log(Math.PI) + lgamma(b) - lgamma(0.5 + b);
    var y = (tt * tt) / (tt * tt + df), y2 = df / (tt * tt + df);
    var xodd = ibetaPQ_(y, y2, a, b)[0];
    var godd = 2 * rxb * Math.exp(a * Math.log(x) - albeta);
    tnc = b * x;
    var xeven = tnc < 2.220446049250313e-16 ? tnc : 1 - rxb;
    var geven = tnc * rxb;
    tnc = p * xodd + q * xeven;
    for (var it = 1; it <= 1000; it++) {
      a += 1;
      xodd -= godd;
      xeven -= geven;
      godd *= x * (a + b - 1) / a;
      geven *= x * (a + b - 0.5) / (a + 0.5);
      p *= lambda / (2 * it);
      q *= lambda / (2 * it + 1);
      tnc += p * xodd + q * xeven;
      s -= p;
      if (s < -1e-10) break;
      if (s >= 0 && Math.abs(2 * s * (xodd - godd)) < 1e-12) break;
    }
  } else {
    tnc = 0;
  }
  tnc += pnorm(-del, true);
  tnc = Math.min(tnc, 1);
  return negdel ? 1 - tnc : tnc;
}

/* ---------- χ² 분포 ---------- */

function pchisq(x, df, lowerTail) {
  if (isNaN(x) || isNaN(df) || df <= 0) return NaN;
  var lower = lowerTail !== false;
  var pq = gammaPQ_(df / 2, x / 2);
  return lower ? pq[0] : pq[1];
}

function dchisqPdf_(x, df) {
  var k = df / 2;
  return Math.exp((k - 1) * Math.log(x) - x / 2 - k * Math.LN2 - lgamma(k));
}

function qchisq(p, df, lowerTail) {
  if (isNaN(p) || isNaN(df) || p < 0 || p > 1 || df <= 0) return NaN;
  var lower = lowerTail !== false;
  if (p === 0) return lower ? 0 : Infinity;
  if (p === 1) return lower ? Infinity : 0;
  var s = tailSelect_(p, lowerTail);
  var z = qnorm(p, lowerTail);
  var w = 2 / (9 * df);
  var x0 = df * Math.pow(1 - w + z * Math.sqrt(w), 3); // Wilson–Hilferty
  if (!(x0 > 0)) x0 = df * 0.1;
  return solveTail_(s.u, s.useLower, function (v) { return pchisq(v, df, s.useLower); }, function (v) { return dchisqPdf_(v, df); }, x0);
}

/* ---------- F 분포 ---------- */

function pf(f, df1, df2, lowerTail) {
  if (isNaN(f) || isNaN(df1) || isNaN(df2) || df1 <= 0 || df2 <= 0) return NaN;
  var lower = lowerTail !== false;
  if (f <= 0) return lower ? 0 : 1;
  if (f === Infinity) return lower ? 1 : 0;
  if (!isFinite(df1) && !isFinite(df2)) return NaN;
  if (!isFinite(df2)) return pchisq(f * df1, df1, lowerTail);
  if (!isFinite(df1)) return pchisq(df2 / f, df2, lower ? false : true);
  var u = df1 * f;
  var pq = ibetaPQ_(u / (u + df2), df2 / (u + df2), df1 / 2, df2 / 2);
  return lower ? pq[0] : pq[1];
}

function dfPdf_(x, df1, df2) {
  var a = df1 / 2, b = df2 / 2;
  return Math.exp(a * Math.log(df1 / df2) + (a - 1) * Math.log(x) - (a + b) * Math.log(1 + df1 * x / df2) + lgamma(a + b) - lgamma(a) - lgamma(b));
}

function qf(p, df1, df2, lowerTail) {
  if (isNaN(p) || isNaN(df1) || isNaN(df2) || p < 0 || p > 1 || df1 <= 0 || df2 <= 0) return NaN;
  var lower = lowerTail !== false;
  if (p === 0) return lower ? 0 : Infinity;
  if (p === 1) return lower ? Infinity : 0;
  var s = tailSelect_(p, lowerTail);
  return solveTail_(s.u, s.useLower, function (v) { return pf(v, df1, df2, s.useLower); }, function (v) { return dfPdf_(v, df1, df2); }, 1);
}

/* ---------- Studentized range 분포 (Tukey) ---------- */

/*
 * q = (k개 표준정규 표본의 범위) / S,  S = sqrt(χ²_df / df)  의 분포.
 *   P(Q ≤ q) = ∫ f_S(s) · W(q·s) ds,   W(w) = k ∫ φ(z) [Φ(z) − Φ(z−w)]^(k−1) dz
 * 안쪽(W)과 바깥쪽(S) 적분 모두 가우스–르장드르 구적이다. 노드와 가중치는 실행 시점에 계산한다.
 * 정확도: 누적확률의 절대오차 약 1e-10. 위꼬리는 1 − P 로 구하므로 p값이 1e-10보다 작으면 상대 정확도는 없다.
 */

var GL_CACHE_ = {};

/** n점 가우스–르장드르 노드와 가중치 ([-1, 1]) */
function gaussLegendre_(n) {
  if (GL_CACHE_[n]) return GL_CACHE_[n];
  var x = [], w = [], m = (n + 1) >> 1;
  for (var i = 0; i < m; i++) {
    var z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5)), pp = 1, z1;
    for (var it = 0; it < 100; it++) {
      var p1 = 1, p2 = 0;
      for (var j = 0; j < n; j++) {
        var p3 = p2;
        p2 = p1;
        p1 = ((2 * j + 1) * z * p2 - j * p3) / (j + 1);
      }
      pp = n * (z * p1 - p2) / (z * z - 1);
      z1 = z;
      z = z1 - p1 / pp;
      if (Math.abs(z - z1) < 1e-15) break;
    }
    x[i] = -z; x[n - 1 - i] = z;
    w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * pp * pp);
  }
  return (GL_CACHE_[n] = { x: x, w: w });
}

var TUKEY_Z_ = null; // 안쪽 적분의 z 노드 (정규 밀도와 누적확률을 미리 계산해 둔다)
var TUKEY_ZMAX_ = 8;

function tukeyZGrid_() {
  if (TUKEY_Z_) return TUKEY_Z_;
  var gl = gaussLegendre_(16);
  var g = { z: [], wt: [], phi: [], cdf: [] };
  for (var p = 0; p < 2 * TUKEY_ZMAX_; p++) { // 폭 1인 구간 16개
    var mid = -TUKEY_ZMAX_ + p + 0.5;
    for (var i = 0; i < gl.x.length; i++) {
      var z = mid + 0.5 * gl.x[i];
      g.z.push(z);
      g.wt.push(0.5 * gl.w[i]);
      g.phi.push(Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI));
      g.cdf.push(pnorm(z));
    }
  }
  return (TUKEY_Z_ = g);
}

/** k개 표준정규 표본의 범위가 w 이하일 확률 */
function rangeCdf_(w, k) {
  if (!(w > 0)) return 0;
  if (w === Infinity) return 1;
  var g = tukeyZGrid_(), km1 = k - 1, sum = 0;
  for (var i = 0; i < g.z.length; i++) {
    var d = g.cdf[i] - pnorm(g.z[i] - w);
    if (d > 0) sum += g.wt[i] * g.phi[i] * (km1 === 1 ? d : Math.pow(d, km1));
  }
  return Math.min(1, k * sum);
}

var TUKEY_BOUNDS_ = {};
var TUKEY_WMAX_ = 20; // 범위 w ≥ 20 이면 W(w) = 1 로 본다 (k ≤ 100 에서 오차 1e-13 이하)

/** S = sqrt(χ²/df) 의 질량이 거의 모두 들어 있는 구간 (양쪽 꼬리 1e-16) */
function tukeyBounds_(df) {
  if (TUKEY_BOUNDS_[df]) return TUKEY_BOUNDS_[df];
  var keys = Object.keys(TUKEY_BOUNDS_);
  if (keys.length > 256) delete TUKEY_BOUNDS_[keys[0]]; // 게임-하웰처럼 df가 계속 바뀌어도 메모리를 제한한다
  return (TUKEY_BOUNDS_[df] = {
    lo: Math.sqrt(qchisq(1e-16, df) / df),
    hi: Math.sqrt(qchisq(1e-16, df, false) / df)
  });
}

/** P(Q ≤ q) */
function ptukeyLower_(q, k, df) {
  if (!(q > 0)) return 0;
  if (q === Infinity) return 1;
  if (!isFinite(df)) return rangeCdf_(q, k);
  var b = tukeyBounds_(df);
  // q·s ≥ WMAX 인 구간은 W = 1 이므로 그 구간의 확률 P(S > WMAX/q) 는 χ² 위꼬리로 바로 더한다
  var top = b.hi, tail = 0, sc = TUKEY_WMAX_ / q;
  if (sc < top) {
    top = sc;
    tail = pchisq(df * sc * sc, df, false);
  }
  if (top <= b.lo) return Math.min(1, tail);
  // W(q·s) 는 s 에 대해 폭 ~ 1/q 로 변하므로 구간 폭을 4/q 이하로 제한한다
  var width = Math.min(3 / Math.sqrt(2 * df), 1, 4 / q);
  var panels = Math.max(1, Math.ceil((top - b.lo) / width));
  var h = (top - b.lo) / panels;
  var gl = gaussLegendre_(16), sum = 0;
  for (var p = 0; p < panels; p++) {
    var mid = b.lo + (p + 0.5) * h;
    for (var i = 0; i < gl.x.length; i++) {
      var s = mid + 0.5 * h * gl.x[i];
      // S 의 밀도: f_S(s) = df · dpoisRaw(df/2, df·s²/2) / s
      var f = df * dpoisRaw_(df / 2, df * s * s / 2) / s;
      sum += 0.5 * h * gl.w[i] * f * rangeCdf_(q * s, k);
    }
  }
  return Math.min(1, Math.max(0, sum + tail));
}

/** k: 비교하는 집단(평균)의 수 (k ≥ 2), df: 오차 자유도 */
function ptukey(q, k, df, lowerTail) {
  if (isNaN(q) || isNaN(k) || isNaN(df) || k < 2 || df <= 0) return NaN;
  var lower = ptukeyLower_(q, k, df);
  return lowerTail === false ? 1 - lower : lower;
}

function qtukey(p, k, df, lowerTail) {
  if (isNaN(p) || isNaN(k) || isNaN(df) || p < 0 || p > 1 || k < 2 || df <= 0) return NaN;
  var target = lowerTail === false ? 1 - p : p;
  if (target <= 0) return 0;
  if (target >= 1) return Infinity;
  // k = 2 일 때 Q = √2·|T| 이므로 그 분위수가 해의 하한이다 (같은 q에서 P는 k가 클수록 작다)
  var lo = Math.SQRT2 * qt((1 - target) / 2, df, false);
  if (!isFinite(df)) lo = Math.SQRT2 * qnorm((1 - target) / 2, false);
  var f = function (q) { return ptukeyLower_(q, k, df) - target; };
  var flo = f(lo);
  if (flo >= 0) return lo; // k = 2
  var hi = lo * 1.5 + 1, fhi = f(hi);
  for (var g = 0; fhi < 0 && g < 60; g++) { lo = hi; flo = fhi; hi *= 2; fhi = f(hi); }
  // Illinois 방식 가위치법 (구간을 항상 유지)
  var side = 0;
  for (var i = 0; i < 100; i++) {
    var q = (lo * fhi - hi * flo) / (fhi - flo);
    if (!(q > lo && q < hi)) q = (lo + hi) / 2;
    var fq = f(q);
    if (fq === 0 || Math.abs(hi - lo) < 1e-13 * q) return q;
    if (fq * fhi > 0) {
      hi = q; fhi = fq;
      if (side === -1) flo /= 2;
      side = -1;
    } else {
      lo = q; flo = fq;
      if (side === 1) fhi /= 2;
      side = 1;
    }
    if (Math.abs(fq) < 1e-14) return q;
  }
  return (lo + hi) / 2;
}

/* ---------- p값 도우미 ---------- */

/** 양측 p값 */
function pValueNormal(z) { return 2 * pnorm(-Math.abs(z)); }
function pValueT(t, df) { return 2 * pt(-Math.abs(t), df); }
/** 오른쪽 꼬리 p값 */
function pValueF(f, df1, df2) { return pf(f, df1, df2, false); }
function pValueChisq(x, df) { return pchisq(x, df, false); }

if (typeof module !== 'undefined') {
  module.exports = {
    lgamma: lgamma, ibetaPQ_: ibetaPQ_,
    pnorm: pnorm, qnorm: qnorm,
    pt: pt, qt: qt, pnct: pnct,
    pchisq: pchisq, qchisq: qchisq,
    pf: pf, qf: qf,
    ptukey: ptukey, qtukey: qtukey,
    pValueNormal: pValueNormal, pValueT: pValueT, pValueF: pValueF, pValueChisq: pValueChisq
  };
}
