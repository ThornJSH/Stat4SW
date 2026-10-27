/**
 * 일원 분산분석(ANOVA)과 사후검정. 순수 함수: (dataset, request) → 결과 객체.
 * 평균차이의 방향은 집단 I − 집단 J 이다.
 *
 * 사후검정: Tukey HSD(표본 크기가 다르면 Tukey-Kramer), Games-Howell(등분산 가정 안 함),
 *           Bonferroni, Scheffé, LSD(보정 없음), Šidák, Dunnett(통제집단 대비, 양측)
 * Dunnett: 통제집단과 각 집단의 합동분산 t 통계량 T_i 의 다변량 t 분포에서 P(max|T_i| < q)를 구한다.
 *   T_i = (a_i Z_i − b_i Z_0)/u (a_i = √(n_c/(n_i+n_c)), b_i = √(n_i/(n_i+n_c)), u = √(χ²_df/df))이므로 Z_0 와 u 에 대한 이중 적분이고,
 *   가우스–르장드르 구적(Z_0 는 [−8, 8], u 는 χ² 분위수 공간 (0,1))으로 계산한다. 유의확률 = 1 − P(|t_i|), 임계값은 P(q) = 신뢰수준을 풀어 구한다.
 * 강건 검정: Welch, Brown-Forsythe
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function anovaDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var desc = require('./Stats_Descriptive.js');
    var asm = require('./Stats_Assumptions.js');
    return {
      getVars: data.getVars, splitGroups: data.splitGroups,
      centralSums: desc.centralSums, leveneTest: asm.leveneTest,
      pf: dist.pf, qf: dist.qf, pt: dist.pt, qt: dist.qt, ptukey: dist.ptukey, qtukey: dist.qtukey, pnorm: dist.pnorm, qchisq: dist.qchisq, lgamma: dist.lgamma,
      boxplotBlock: require('./Charts.js').boxplotBlock,
      advAnova_: asm.advAnova_, advAdd_: asm.advAdd_, advAppend_: asm.advAppend_
    };
  }
  return {
    getVars: getVars, splitGroups: splitGroups, centralSums: centralSums, leveneTest: leveneTest,
    pf: pf, qf: qf, pt: pt, qt: qt, ptukey: ptukey, qtukey: qtukey, pnorm: pnorm, qchisq: qchisq, lgamma: lgamma,
    boxplotBlock: boxplotBlock,
    advAnova_: advAnova_, advAdd_: advAdd_, advAppend_: advAppend_
  };
}

var ANOVA_MAX_GROUPS_ = 20;

/* ---------- 계산 ---------- */

/** 집단별 N, 평균, 분산 (분산은 N ≥ 2 일 때) */
function groupMoments_(groups) {
  var d = anovaDeps_();
  return groups.map(function (g) {
    var cs = d.centralSums(g);
    return { n: g.length, mean: cs.mean, c2: cs.c2, variance: g.length > 1 ? cs.c2 / (g.length - 1) : null };
  });
}

/** 한 집단의 기술통계와 평균의 신뢰구간 */
function groupSummary_(values, conf) {
  var d = anovaDeps_();
  var n = values.length;
  var r = { n: n, mean: null, sd: null, se: null, lo: null, hi: null, min: null, max: null };
  if (n === 0) return r;
  var cs = d.centralSums(values);
  r.mean = cs.mean;
  // (Math.min.apply 는 큰 배열에서 호출 스택 한계로 실패하므로 반복문을 쓴다)
  r.min = values[0]; r.max = values[0];
  for (var q = 1; q < n; q++) {
    if (values[q] < r.min) r.min = values[q];
    else if (values[q] > r.max) r.max = values[q];
  }
  if (n > 1) {
    r.sd = Math.sqrt(cs.c2 / (n - 1));
    r.se = r.sd / Math.sqrt(n);
    var tc = d.qt((1 - conf / 100) / 2, n - 1, false);
    r.lo = r.mean - tc * r.se;
    r.hi = r.mean + tc * r.se;
  }
  return r;
}

/**
 * 일원 분산분석. groups: 집단별 숫자 배열.
 * 집단 내 변동이 0이거나 자유도가 0이면 F, p 는 null.
 */
function oneWayAnova(groups) {
  var d = anovaDeps_();
  var m = groupMoments_(groups);
  var k = groups.length, N = 0, i;
  m.forEach(function (g) { N += g.n; });
  var grand = 0;
  m.forEach(function (g) { grand += g.n * g.mean; });
  grand /= N;

  var ssb = 0, ssw = 0;
  m.forEach(function (g) {
    ssb += g.n * (g.mean - grand) * (g.mean - grand);
    ssw += g.c2;
  });
  var r = {
    k: k, N: N, ssb: ssb, ssw: ssw, sst: ssb + ssw,
    dfb: k - 1, dfw: N - k, dft: N - 1,
    msb: ssb / (k - 1), msw: null, F: null, p: null, eta2: null, omega2: null
  };
  if (r.dfw < 1) return r;
  r.msw = ssw / r.dfw;
  if (r.sst > 0) r.eta2 = ssb / r.sst;
  if (!(r.msw > 0)) return r;
  r.F = r.msb / r.msw;
  r.p = d.pf(r.F, r.dfb, r.dfw, false);
  r.omega2 = (ssb - r.dfb * r.msw) / (r.sst + r.msw);
  return r;
}

/** Welch 분산분석. 모든 집단의 사례가 2개 이상이고 분산이 0보다 커야 한다. */
function welchAnova(groups) {
  var d = anovaDeps_();
  var m = groupMoments_(groups), k = groups.length;
  var r = { F: null, df1: k - 1, df2: null, p: null };
  if (k < 2) return r;
  for (var i = 0; i < k; i++) if (m[i].n < 2 || !(m[i].variance > 0)) return r;
  var W = 0, wm = 0;
  var w = m.map(function (g) { return g.n / g.variance; });
  w.forEach(function (x, j) { W += x; wm += x * m[j].mean; });
  wm /= W;
  var num = 0, lam = 0;
  w.forEach(function (x, j) {
    num += x * (m[j].mean - wm) * (m[j].mean - wm);
    lam += (1 - x / W) * (1 - x / W) / (m[j].n - 1);
  });
  num /= (k - 1);
  r.F = num / (1 + 2 * (k - 2) / (k * k - 1) * lam);
  r.df2 = (k * k - 1) / (3 * lam);
  r.p = d.pf(r.F, r.df1, r.df2, false);
  return r;
}

/** Brown-Forsythe 검정 (F*). Welch와 같은 조건이 필요하다. */
function brownForsythe(groups) {
  var d = anovaDeps_();
  var m = groupMoments_(groups), k = groups.length, N = 0;
  m.forEach(function (g) { N += g.n; });
  var r = { F: null, df1: k - 1, df2: null, p: null };
  if (k < 2) return r;
  for (var i = 0; i < k; i++) if (m[i].n < 2 || !(m[i].variance > 0)) return r;
  var grand = 0;
  m.forEach(function (g) { grand += g.n * g.mean; });
  grand /= N;
  var num = 0, den = 0;
  m.forEach(function (g) {
    num += g.n * (g.mean - grand) * (g.mean - grand);
    den += (1 - g.n / N) * g.variance;
  });
  r.F = num / den;
  var inv = 0;
  m.forEach(function (g) {
    var c = (1 - g.n / N) * g.variance / den;
    inv += c * c / (g.n - 1);
  });
  r.df2 = 1 / inv;
  r.p = d.pf(r.F, r.df1, r.df2, false);
  return r;
}

var POSTHOC_METHODS_ = [
  { key: 'tukey', label: 'Tukey HSD', note: 'Tukey HSD (표본 크기가 다르면 Tukey-Kramer). 표준오차는 √(MSE(1/nᵢ+1/nⱼ)), 신뢰구간은 스튜던트화 범위 분포의 임계값을 쓴다.' },
  { key: 'gameshowell', label: 'Games-Howell', note: 'Games-Howell: 등분산을 가정하지 않는다 (Welch 자유도와 스튜던트화 범위 분포). 모든 집단의 사례가 2개 이상이어야 한다.' },
  { key: 'bonferroni', label: 'Bonferroni', note: 'Bonferroni: 합동분산 t검정의 유의확률에 비교 횟수 m = k(k−1)/2 를 곱한다 (최대 1).' },
  { key: 'scheffe', label: 'Scheffé', note: 'Scheffé: F 분포에 기반하며 가장 보수적이다.' },
  { key: 'lsd', label: 'LSD (최소유의차)', note: 'LSD: 합동분산 t검정을 다중비교 보정 없이 그대로 쓴다 (유의확률과 신뢰구간 모두 보정 없음). 비교가 늘면 1종 오류가 커지므로 전체 F가 유의할 때, 집단이 적을 때만 쓰는 것이 관례이다.' },
  { key: 'sidak', label: 'Šidák', note: 'Šidák: 합동분산 t검정의 유의확률을 1−(1−p)^m 으로 보정한다 (m = k(k−1)/2). Bonferroni보다 약간 덜 보수적이다.' },
  { key: 'dunnett', label: 'Dunnett (통제집단 대비)', note: 'Dunnett: 통제집단 하나와 나머지 집단을 각각 비교하는 양측 검정이다 (평균차이 = (I) 집단 − 통제집단, 합동분산, 다변량 t 분포로 보정). 통제집단과의 비교만 하므로 모든 쌍을 비교하는 방법보다 검정력이 높다.' }
];

/* ---------- Dunnett 다변량 t 확률 ---------- */

var ANOVA_GL_ = {};
/** n점 가우스–르장드르 구적의 [−1, 1] 위 절점과 가중치 (뉴턴법, 결과를 저장해 재사용) */
function anovaGL_(n) {
  if (ANOVA_GL_[n]) return ANOVA_GL_[n];
  var x = [], w = [], i, j;
  for (i = 0; i < n; i++) {
    var z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5)), pp = 0;
    for (var it = 0; it < 100; it++) {
      var p1 = 1, p2 = 0;
      for (j = 1; j <= n; j++) { var p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
      pp = n * (z * p1 - p2) / (z * z - 1);
      var dz = p1 / pp;
      z -= dz;
      if (Math.abs(dz) < 1e-15) break;
    }
    x.push(z); w.push(2 / ((1 - z * z) * pp * pp));
  }
  ANOVA_GL_[n] = { x: x, w: w };
  return ANOVA_GL_[n];
}

var ANOVA_DUNNETT_NODES_ = 8;

/** ns: 비교 집단의 크기들, nc: 통제집단 크기, df: 오차 자유도 → P(모든 |T_i| < q)를 구하는 함수 */
function anovaDunnettProb_(ns, nc, df) {
  var d = anovaDeps_();
  var cnt = {}, cats = [];
  ns.forEach(function (n) {
    if (!(n in cnt)) { cnt[n] = 0; cats.push(n); }
    cnt[n]++;
  });
  var parts = cats.map(function (n) { return { a: Math.sqrt(nc / (n + nc)), b: Math.sqrt(n / (n + nc)), c: cnt[n] }; });
  // Z_0 는 [−8, 8] 을 같은 폭의 구간으로 나눠 적분한다 (통제집단이 작으면 인자가 가팔라서 구간 수를 늘린다)
  var aMin = Math.min.apply(null, parts.map(function (pt) { return pt.a; }));
  var zPanels = aMin < 0.5 ? 32 : 16, gz = anovaGL_(ANOVA_DUNNETT_NODES_), zn = [];
  for (var zp = 0; zp < zPanels; zp++) {
    var zhalf = 8 / zPanels, zmid = -8 + zhalf * (2 * zp + 1);
    for (var zr = 0; zr < gz.x.length; zr++) {
      var z = zmid + zhalf * gz.x[zr];
      zn.push({ z: z, w: zhalf * gz.w[zr] * Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI) });
    }
  }
  // u = √(χ²/df)의 밀도를 분위수로 나눈 구간마다 적분한다 (구간 안에서는 부드러워 빨리 수렴)
  var bp = [1e-14, 1e-4, 1e-2, 0.1, 0.5].map(function (p) { return Math.sqrt(d.qchisq(p, df, true) / df); })
    .concat([0.1, 1e-2, 1e-4, 1e-14].map(function (p) { return Math.sqrt(d.qchisq(p, df, false) / df); }));
  var gp = anovaGL_(12), kk = df / 2;
  var lc = Math.LN2 + Math.log(df) + (kk - 1) * Math.log(df) - kk * Math.LN2 - d.lgamma(kk);
  var un = [];
  for (var s = 0; s < bp.length - 1; s++) {
    var half = (bp[s + 1] - bp[s]) / 2, mid = (bp[s + 1] + bp[s]) / 2;
    for (var r = 0; r < gp.x.length; r++) {
      var uu = mid + half * gp.x[r];
      un.push({ u: uu, w: half * gp.w[r] * Math.exp(lc + (2 * kk - 1) * Math.log(uu) - df * uu * uu / 2) });
    }
  }
  return function (q) {
    var total = 0;
    for (var a = 0; a < un.length; a++) {
      var qu = q * un[a].u, inner = 0;
      for (var b = 0; b < zn.length; b++) {
        var prod = 1;
        for (var c = 0; c < parts.length; c++) {
          var pt = parts[c], bz = pt.b * zn[b].z;
          var diff = d.pnorm((qu + bz) / pt.a) - d.pnorm((-qu + bz) / pt.a);
          prod *= pt.c === 1 ? diff : Math.pow(diff, pt.c);
        }
        inner += zn[b].w * prod;
      }
      total += un[a].w * inner;
    }
    return Math.min(1, Math.max(0, total));
  };
}

/** 단조 증가 함수 g 의 해 g(x) = target (구간 [lo, hi], 일리노이 거짓위치법) */
function anovaSolve_(g, target, lo, hi) {
  var flo = g(lo) - target, fhi = g(hi) - target, side = 0;
  for (var it = 0; it < 80; it++) {
    if (Math.abs(fhi) < 1e-14 || Math.abs(hi - lo) < 1e-13) break;
    var x = (lo * fhi - hi * flo) / (fhi - flo);
    var fx = g(x) - target;
    if (fx * fhi > 0) {
      hi = x; fhi = fx;
      if (side === -1) flo /= 2;
      side = -1;
    } else {
      lo = x; flo = fx;
      if (side === 1) fhi /= 2;
      side = 1;
    }
    if (Math.abs(fx) < 1e-14) return x;
  }
  return Math.abs(flo) < Math.abs(fhi) ? lo : hi;
}

/**
 * 사후검정의 쌍별 비교. 모든 쌍 (i < j) 을 한 번씩, 평균차이 = 집단 i − 집단 j.
 * 반환: [{ i, j, diff, se, p, lo, hi }]  (계산할 수 없으면 se, p, lo, hi 는 null)
 * conf: 신뢰수준(%)
 * Dunnett 은 통제집단(control, 집단 번호, 기본 0)과 나머지 집단의 쌍만 돌려준다 ({ i: 나머지 집단, j: 통제집단 }, 평균차이 = 집단 i − 통제집단).
 */
function postHoc(method, groups, conf, control) {
  var d = anovaDeps_();
  var m = groupMoments_(groups), k = groups.length;
  var a = oneWayAnova(groups);
  var alpha = 1 - conf / 100;
  var pairs = k * (k - 1) / 2;
  var rows = [];

  if (method === 'dunnett') {
    var ctrl = control === undefined || control === null ? 0 : control;
    if (ctrl < 0 || ctrl >= k) throw new Error('통제집단 번호가 범위를 벗어났습니다.');
    var okD = a.msw !== null && a.msw > 0;
    var others = [];
    for (var o = 0; o < k; o++) if (o !== ctrl) others.push(o);
    var prob = okD ? anovaDunnettProb_(others.map(function (g) { return m[g].n; }), m[ctrl].n, a.dfw) : null;
    var qd = null;
    if (okD) {
      var lo0 = d.qt(alpha / 2, a.dfw, false), hi0 = d.qt(alpha / (2 * others.length), a.dfw, false);
      qd = others.length === 1 ? lo0 : anovaSolve_(prob, 1 - alpha, lo0 * 0.999, hi0 * 1.001);
    }
    others.forEach(function (g) {
      var r = { i: g, j: ctrl, diff: m[g].mean - m[ctrl].mean, se: null, p: null, lo: null, hi: null };
      rows.push(r);
      if (!okD) return;
      r.se = Math.sqrt(a.msw * (1 / m[g].n + 1 / m[ctrl].n));
      r.p = others.length === 1 ? 2 * d.pt(-Math.abs(r.diff / r.se), a.dfw) : Math.min(1, Math.max(0, 1 - prob(Math.abs(r.diff / r.se))));
      r.lo = r.diff - qd * r.se;
      r.hi = r.diff + qd * r.se;
    });
    return rows;
  }

  // 합동분산을 쓰는 방법의 공통 값
  var pooledOk = a.msw !== null && a.msw > 0;
  var qcPooled = null, tcBonf = null, fcScheffe = null, tcLsd = null, tcSidak = null;
  if (pooledOk && method === 'tukey') qcPooled = d.qtukey(alpha, k, a.dfw, false);
  if (pooledOk && method === 'bonferroni') tcBonf = d.qt(alpha / (2 * pairs), a.dfw, false);
  if (pooledOk && method === 'lsd') tcLsd = d.qt(alpha / 2, a.dfw, false);
  if (pooledOk && method === 'sidak') tcSidak = d.qt((1 - Math.pow(1 - alpha, 1 / pairs)) / 2, a.dfw, false);
  if (pooledOk && method === 'scheffe') fcScheffe = Math.sqrt((k - 1) * d.qf(alpha, k - 1, a.dfw, false));

  for (var i = 0; i < k - 1; i++) {
    for (var j = i + 1; j < k; j++) {
      var row = { i: i, j: j, diff: m[i].mean - m[j].mean, se: null, p: null, lo: null, hi: null };
      rows.push(row);
      var t, hw;
      if (method === 'gameshowell') {
        if (m[i].n < 2 || m[j].n < 2 || !(m[i].variance > 0 || m[j].variance > 0)) continue;
        var vi = m[i].variance / m[i].n, vj = m[j].variance / m[j].n;
        row.se = Math.sqrt(vi + vj);
        var wdf = (vi + vj) * (vi + vj) / (vi * vi / (m[i].n - 1) + vj * vj / (m[j].n - 1));
        row.p = d.ptukey(Math.abs(row.diff) * Math.SQRT2 / row.se, k, wdf, false);
        hw = d.qtukey(alpha, k, wdf, false) * row.se / Math.SQRT2;
      } else {
        if (!pooledOk) continue;
        row.se = Math.sqrt(a.msw * (1 / m[i].n + 1 / m[j].n));
        t = row.diff / row.se;
        if (method === 'tukey') {
          row.p = d.ptukey(Math.abs(t) * Math.SQRT2, k, a.dfw, false);
          hw = qcPooled * row.se / Math.SQRT2;
        } else if (method === 'bonferroni') {
          row.p = Math.min(1, pairs * d.pt(-Math.abs(t), a.dfw) * 2);
          hw = tcBonf * row.se;
        } else if (method === 'lsd') {
          row.p = 2 * d.pt(-Math.abs(t), a.dfw);
          hw = tcLsd * row.se;
        } else if (method === 'sidak') {
          row.p = -Math.expm1(pairs * Math.log1p(-Math.min(1, 2 * d.pt(-Math.abs(t), a.dfw))));
          hw = tcSidak * row.se;
        } else if (method === 'scheffe') {
          row.p = d.pf(t * t / (k - 1), k - 1, a.dfw, false);
          hw = fcScheffe * row.se;
        } else {
          throw new Error('알 수 없는 사후검정입니다: ' + method);
        }
      }
      row.lo = row.diff - hw;
      row.hi = row.diff + hw;
    }
  }
  return rows;
}

/* ---------- 결과 ---------- */

function nameOf_anova_(v) { return v.label ? v.label : v.name; }

function runOneWayAnova(dataset, request) {
  var d = anovaDeps_();
  var depVars = d.getVars(dataset, request.dependent, { numeric: true });
  var factorVars = d.getVars(dataset, request.factor);
  if (factorVars.length !== 1) throw new Error('집단변수는 하나만 선택하세요.');
  var factor = factorVars[0];
  depVars.forEach(function (v) {
    if (v.name === factor.name) throw new Error('종속변수와 집단변수에 같은 변수 "' + v.name + '"를 쓸 수 없습니다.');
  });
  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  var alpha = 1 - conf / 100;
  var factorName = nameOf_anova_(factor);
  var methods = POSTHOC_METHODS_.filter(function (mt) { return opt[mt.key]; });
  var controlText = opt.control === undefined || opt.control === null ? '' : String(opt.control).trim();
  var controlLabel = '';
  var warnings = [];

  var descRows = [], anovaRows = [], levRows = [], effectRows = [], robustRows = [];
  var postRows = {};
  methods.forEach(function (mt) { postRows[mt.key] = []; });
  var groupLabels = [], boxBlocks = [];
  var adv = { rows: [], recs: [] };

  depVars.forEach(function (v) {
    var name = nameOf_anova_(v);
    var split = d.splitGroups(dataset, v.name, factor);
    var groups = split.groups;
    if (groups.length < 2) {
      throw new Error('집단변수 "' + factorName + '"의 집단이 ' + groups.length + '개입니다 (종속변수 "' + name + '"에 유효한 사례 기준). 집단이 2개 이상이어야 합니다.');
    }
    if (groups.length > ANOVA_MAX_GROUPS_) {
      throw new Error('집단변수 "' + factorName + '"의 집단이 ' + groups.length + '개입니다 (최대 ' + ANOVA_MAX_GROUPS_ + '개). 범주가 적은 변수를 선택하세요.');
    }
    var values = groups.map(function (g) { return g.values; });
    if (opt.assumptions) d.advAdd_(adv, d.advAnova_(name, groups.map(function (g) { return { label: g.label, values: g.values }; }), v.scale === 'ordinal'));
    if (opt.boxplot) {
      boxBlocks.push(d.boxplotBlock('상자그림: ' + name + ' (집단변수: ' + factorName + ')',
        groups.map(function (g) { return { label: g.label, values: g.values }; }), name));
    }
    var all = [];
    values.forEach(function (vs) { all = all.concat(vs); });

    // 기술통계
    groups.forEach(function (g, idx) {
      var s = groupSummary_(g.values, conf);
      descRows.push({ name: idx === 0 ? name : '', group: g.label, n: s.n, mean: s.mean, sd: s.sd, se: s.se, lo: s.lo, hi: s.hi, min: s.min, max: s.max });
      if (s.n < 2) warnings.push('"' + name + '": 집단 "' + g.label + '"의 사례가 1개뿐입니다.');
    });
    var tot = groupSummary_(all, conf);
    descRows.push({ _kind: 'total', name: '', group: '전체', n: tot.n, mean: tot.mean, sd: tot.sd, se: tot.se, lo: tot.lo, hi: tot.hi, min: tot.min, max: tot.max });

    // 분산분석
    var a = oneWayAnova(values);
    if (a.F === null) warnings.push('"' + name + '": 집단 내 변동이 없거나 자유도가 0이라 F를 계산할 수 없습니다.');
    anovaRows.push({ name: name, source: '집단-간', ss: a.ssb, df: a.dfb, ms: a.msb, F: a.F, p: a.p });
    anovaRows.push({ name: '', source: '집단-내', ss: a.ssw, df: a.dfw, ms: a.msw });
    anovaRows.push({ _kind: 'total', name: '', source: '합계', ss: a.sst, df: a.dft });
    effectRows.push({ name: name, eta2: a.eta2, omega2: a.omega2 });

    if (opt.levene) {
      var lev = d.leveneTest(values, false);
      levRows.push({ name: name, W: lev.W, df1: lev.df1, df2: lev.df2, p: lev.p });
    }
    if (opt.welch) {
      var w = welchAnova(values), bf = brownForsythe(values);
      robustRows.push({ name: name, test: 'Welch', F: w.F, df1: w.df1, df2: w.df2, p: w.p });
      robustRows.push({ name: '', test: 'Brown-Forsythe', F: bf.F, df1: bf.df1, df2: bf.df2, p: bf.p });
      if (w.F === null) warnings.push('"' + name + '": 모든 집단의 사례가 2개 이상이고 분산이 0보다 커야 Welch·Brown-Forsythe 검정을 할 수 있습니다.');
    }

    // 사후검정
    var ctrlIdx = 0;
    if (opt.dunnett && controlText !== '') {
      ctrlIdx = -1;
      groups.forEach(function (g, q) { if (ctrlIdx < 0 && (String(g.value) === controlText || String(g.label).trim() === controlText)) ctrlIdx = q; });
      if (ctrlIdx < 0) throw new Error('Dunnett의 통제집단 "' + controlText + '"이(가) 집단변수 "' + factorName + '"의 값이나 값 레이블에 없습니다 (집단: ' + groups.map(function (g) { return String(g.label); }).join(', ') + '). 비워 두면 첫 집단이 통제집단입니다.');
    }
    controlLabel = String(groups[ctrlIdx].label);
    methods.forEach(function (mt) {
      var res = postHoc(mt.key, values, conf, ctrlIdx);
      if (res.length && res.every(function (r) { return r.p === null; })) {
        warnings.push('"' + name + '": ' + mt.label + ' 검정을 할 수 없습니다 (집단 내 변동이 없거나 사례가 부족).');
      }
      res.forEach(function (r, idx) {
        postRows[mt.key].push({
          name: idx === 0 ? name : '', gi: groups[r.i].label, gj: groups[r.j].label,
          diff: r.diff, sig: r.p !== null && r.p < alpha ? '*' : '',
          se: r.se, p: r.p, lo: r.lo, hi: r.hi
        });
      });
    });
    groupLabels.push(groups.length);
  });

  var F3 = { type: 'fixed', digits: 3 };
  var blocks = [];

  blocks.push({
    type: 'table',
    caption: '기술통계 (집단변수: ' + factorName + ')',
    columns: [
      { key: 'name', label: '종속변수' }, { key: 'group', label: '집단' },
      { key: 'n', label: 'N', format: { type: 'int' } },
      { key: 'mean', label: '평균', format: F3 }, { key: 'sd', label: '표준편차', format: F3 },
      { key: 'se', label: '표준오차', format: F3 },
      { key: 'lo', label: '평균의 ' + conf + '% 신뢰구간 하한', format: F3 },
      { key: 'hi', label: '평균의 ' + conf + '% 신뢰구간 상한', format: F3 },
      { key: 'min', label: '최솟값', format: { type: 'trim', digits: 3 } },
      { key: 'max', label: '최댓값', format: { type: 'trim', digits: 3 } }
    ],
    rows: descRows,
    notes: ['종속변수와 집단변수가 모두 유효한 사례만 사용했다.']
  });

  if (opt.levene) {
    blocks.push({
      type: 'table',
      caption: '등분산 검정 (Levene, 평균 기준)',
      columns: [
        { key: 'name', label: '종속변수' },
        { key: 'W', label: 'Levene 통계량', format: F3 },
        { key: 'df1', label: '자유도1', format: { type: 'int' } },
        { key: 'df2', label: '자유도2', format: { type: 'int' } },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
      ],
      rows: levRows,
      notes: ['유의확률이 작으면(예: .05 미만) 등분산 가정을 기각한다. 이때는 Welch 검정과 Games-Howell 사후검정을 참고한다.']
    });
  }

  blocks.push({
    type: 'table',
    caption: '분산분석',
    columns: [
      { key: 'name', label: '종속변수' }, { key: 'source', label: '' },
      { key: 'ss', label: '제곱합', format: F3 },
      { key: 'df', label: '자유도', format: { type: 'int' } },
      { key: 'ms', label: '평균제곱', format: F3 },
      { key: 'F', label: 'F', format: F3 },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
    ],
    rows: anovaRows,
    notes: []
  });

  if (opt.effect) {
    blocks.push({
      type: 'table',
      caption: '효과크기',
      columns: [
        { key: 'name', label: '종속변수' },
        { key: 'eta2', label: 'η² (에타제곱)', format: F3 },
        { key: 'omega2', label: 'ω² (오메가제곱)', format: F3 }
      ],
      rows: effectRows,
      notes: ['η² = 집단-간 제곱합 / 합계 제곱합. ω²는 모집단 효과크기의 덜 편향된 추정이다.']
    });
  }

  if (opt.welch) {
    blocks.push({
      type: 'table',
      caption: '평균 동일성의 강건 검정',
      columns: [
        { key: 'name', label: '종속변수' }, { key: 'test', label: '' },
        { key: 'F', label: '통계량', format: F3 },
        { key: 'df1', label: '자유도1', format: { type: 'int' } },
        { key: 'df2', label: '자유도2', format: { type: 'trim', digits: 3 } },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
      ],
      rows: robustRows,
      notes: ['등분산을 가정하지 않는 검정이다.']
    });
  }

  methods.forEach(function (mt) {
    blocks.push({
      type: 'table',
      caption: '다중비교: ' + mt.label,
      columns: [
        { key: 'name', label: '종속변수' },
        { key: 'gi', label: '(I) 집단' }, { key: 'gj', label: '(J) 집단' },
        { key: 'diff', label: '평균차이(I−J)', format: F3 },
        { key: 'sig', label: '', align: 'left' },
        { key: 'se', label: '표준오차', format: F3 },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
        { key: 'lo', label: '차이의 ' + conf + '% 신뢰구간 하한', format: F3 },
        { key: 'hi', label: '차이의 ' + conf + '% 신뢰구간 상한', format: F3 }
      ],
      rows: postRows[mt.key],
      notes: [
        mt.note,
        (mt.key === 'dunnett' ? '통제집단은 "' + controlLabel + '"(J)이다' + (controlText === '' ? ' (옵션을 비워 첫 집단을 썼다)' : '') + '. 통제집단과 나머지 집단의 쌍만 표시했다.' : '모든 쌍을 한 번씩(집단변수 값 순으로 I < J) 표시했다.') + ' * 는 유의확률이 ' + alpha.toFixed(3).replace(/0+$/, '').replace(/\.$/, '') + ' 미만인 쌍이다.'
      ]
    });
  });

  if (adv.rows.length) blocks = d.advAppend_(blocks, adv, true);
  blocks = blocks.concat(boxBlocks);
  return { analysis: 'anova', title: '일원 분산분석', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = {
    oneWayAnova: oneWayAnova,
    welchAnova: welchAnova,
    brownForsythe: brownForsythe,
    postHoc: postHoc,
    runOneWayAnova: runOneWayAnova
  };
}
