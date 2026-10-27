/**
 * 분석 가정 검정: 정규성(Shapiro-Wilk, Jarque-Bera), 등분산(Levene). 순수 함수: (dataset, request) → 결과 객체.
 * p값은 Distributions.js(χ², F)로 계산한다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function assumeDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var desc = require('./Stats_Descriptive.js');
    return {
      getVars: data.getVars, splitGroups: data.splitGroups, pchisq: dist.pchisq, pf: dist.pf, pnorm: dist.pnorm, qnorm: dist.qnorm,
      centralSums: desc.centralSums, describeColumn: desc.describeColumn
    };
  }
  return { getVars: getVars, splitGroups: splitGroups, pchisq: pchisq, pf: pf, pnorm: pnorm, qnorm: qnorm, centralSums: centralSums, describeColumn: describeColumn };
}

/* ---------- Jarque-Bera ---------- */

/**
 * JB = n/6 · (g1² + g2²/4), g1 = m3/m2^1.5, g2 = m4/m2² - 3 (m_k = Σ(x-m)^k / n, 모집단식 적률).
 * 귀무가설(정규분포)에서 JB는 근사적으로 자유도 2인 χ² 분포를 따른다.
 * n < 3 이거나 분산이 0이면 계산할 수 없어 값은 null.
 */
function jarqueBera(values) {
  var d = assumeDeps_();
  var n = values.length;
  var r = { n: n, g1: null, g2: null, jb: null, df: 2, p: null };
  if (n < 3) return r;
  var cs = d.centralSums(values);
  if (!(cs.c2 > 0)) return r;
  var m2 = cs.c2 / n;
  r.g1 = (cs.c3 / n) / (m2 * Math.sqrt(m2));
  r.g2 = (cs.c4 / n) / (m2 * m2) - 3;
  r.jb = n / 6 * (r.g1 * r.g1 + r.g2 * r.g2 / 4);
  r.p = d.pchisq(r.jb, 2, false);
  return r;
}

/* ---------- Shapiro-Wilk ---------- */

var SW_MAX_N_ = 5000;

/** 다항식 cc[0] + cc[1]x + cc[2]x² + … (Royston의 poly) */
function swPoly_(cc, x) {
  var r = cc[0];
  if (cc.length > 1) {
    var p = x * cc[cc.length - 1];
    for (var j = cc.length - 2; j > 0; j--) p = (p + cc[j]) * x;
    r += p;
  }
  return r;
}

/**
 * Shapiro-Wilk 검정 (Royston 1992, 1995; Applied Statistics AS R94). 표본 크기 3~5000.
 * W = (Σ aᵢ x₍ᵢ₎)² / Σ(xᵢ - x̄)², 계수 a는 정규 순서통계량의 근사로 만들고,
 * p값은 n = 3은 정확식, 4~11은 -log(γ - log(1-W)), 12 이상은 log(1-W)가 근사적으로 정규분포를 따른다는 식으로 구한다.
 * 모든 값이 같거나(범위 0) 표본 크기가 범위를 벗어나면 w, p는 null. 이 근사는 scipy.stats.shapiro 와 같은 알고리즘이다.
 */
function shapiroWilk(values) {
  var d = assumeDeps_();
  var n = values.length;
  var r = { n: n, w: null, p: null };
  if (n < 3 || n > SW_MAX_N_) return r;
  var x = values.slice().sort(function (p, q) { return p - q; });
  var range = x[n - 1] - x[0];
  if (!(range >= 1e-19)) return r;

  var nn2 = Math.floor(n / 2);
  var a = [0];                                    // 1부터 쓴다
  var i;
  if (n === 3) {
    a[1] = Math.SQRT1_2;
  } else {
    var an25 = n + 0.25, summ2 = 0;
    for (i = 1; i <= nn2; i++) {
      a[i] = d.qnorm((i - 0.375) / an25);
      summ2 += a[i] * a[i];
    }
    summ2 *= 2;
    var ssumm2 = Math.sqrt(summ2), rsn = 1 / Math.sqrt(n);
    var a1 = swPoly_([0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056], rsn) - a[1] / ssumm2;
    var i1, fac;
    if (n > 5) {
      i1 = 3;
      var a2 = -a[2] / ssumm2 + swPoly_([0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633], rsn);
      fac = Math.sqrt((summ2 - 2 * a[1] * a[1] - 2 * a[2] * a[2]) / (1 - 2 * a1 * a1 - 2 * a2 * a2));
      a[2] = a2;
    } else {
      i1 = 2;
      fac = Math.sqrt((summ2 - 2 * a[1] * a[1]) / (1 - 2 * a1 * a1));
    }
    a[1] = a1;
    for (i = i1; i <= nn2; i++) a[i] /= -fac;
  }

  // 계수 c: 작은 쪽 절반은 -a, 큰 쪽 절반은 +a, (n이 홀수면) 가운데는 0. 범위로 나눈 x와 상관을 W로 쓴다.
  var c = [], sa = 0, sx = 0, k;
  for (k = 0; k < n; k++) {
    var j = n - 1 - k;
    c[k] = k === j ? 0 : (k < j ? -a[k + 1] : a[j + 1]);
    sa += c[k];
    sx += x[k] / range;
  }
  sa /= n; sx /= n;
  var ssa = 0, ssx = 0, sax = 0;
  for (k = 0; k < n; k++) {
    var asa = c[k] - sa, xsx = x[k] / range - sx;
    ssa += asa * asa; ssx += xsx * xsx; sax += asa * xsx;
  }
  // 1-W 를 직접 계산해 W가 1에 아주 가까울 때의 반올림 오차를 피한다
  var ssassx = Math.sqrt(ssa * ssx);
  var w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
  r.w = 1 - w1;

  if (n === 3) {
    var p3 = 1.90985931710274 * (Math.asin(Math.sqrt(r.w)) - 1.04719755119660);   // 6/π · (asin(√W) - π/3)
    r.p = p3 < 0 ? 0 : p3;
    return r;
  }
  var y = Math.log(w1), xx = Math.log(n), m, s;
  if (n <= 11) {
    var gamma = swPoly_([-2.273, 0.459], n);
    if (y >= gamma) { r.p = 1e-99; return r; }
    y = -Math.log(gamma - y);
    m = swPoly_([0.544, -0.39978, 0.025054, -6.714e-4], n);
    s = Math.exp(swPoly_([1.3822, -0.77857, 0.062767, -0.0020322], n));
  } else {
    m = swPoly_([-1.5861, -0.31082, -0.083751, 0.0038915], xx);
    s = Math.exp(swPoly_([-0.4803, -0.082676, 0.0030302], xx));
  }
  r.p = d.pnorm((y - m) / s, false);
  return r;
}

function runNormality(dataset, request) {
  var d = assumeDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  var warnings = [];

  var swRows = [];
  var rows = vars.map(function (v) {
    var values = dataset.data[v.name].filter(function (x) { return x !== null; });
    var r = jarqueBera(values);
    var sw = shapiroWilk(values);
    var name = v.label ? v.label : v.name;
    if (r.jb === null) warnings.push('"' + name + '": 유효한 값이 3개 미만이거나 모든 값이 같아 검정할 수 없습니다.');
    if (sw.w === null && values.length > SW_MAX_N_) {
      warnings.push('"' + name + '": 사례가 ' + SW_MAX_N_ + '개를 넘어 Shapiro-Wilk 검정은 계산하지 않았습니다. Jarque-Bera 결과를 보세요.');
    }
    swRows.push({ name: name, n: sw.n, w: sw.w, df: sw.w === null ? null : sw.n, p: sw.p });
    return { name: name, n: r.n, g1: r.g1, g2: r.g2, jb: r.jb, df: r.jb === null ? null : r.df, p: r.p };
  });

  var verdicts = advNormalityParagraphs_(vars.map(function (v) {
    return { name: v.label ? v.label : v.name, values: dataset.data[v.name].filter(function (x) { return x !== null; }) };
  }));

  return {
    analysis: 'normality',
    title: '정규성 검정',
    blocks: [{
      type: 'table',
      caption: 'Shapiro-Wilk 검정',
      columns: [
        { key: 'name', label: '변수' },
        { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'w', label: 'W 통계량', format: { type: 'fixed', digits: 3 } },
        { key: 'df', label: '자유도', format: { type: 'int' } },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
      ],
      rows: swRows,
      notes: [
        '귀무가설: 변수가 정규분포를 따른다. 유의확률이 작으면(예: .05 미만) 정규성을 기각한다. 표본 크기 3~' + SW_MAX_N_ + '에서 계산한다.',
        'Royston(1995)의 근사식(AS R94)이며 scipy.stats.shapiro, R의 shapiro.test와 같은 방법이다. 작은 표본에서도 쓸 수 있어 Jarque-Bera보다 먼저 본다.',
        '표본이 크면(대략 100 이상) 사소한 이탈에도 정규성이 기각되므로 왜도·첨도, 히스토그램도 함께 본다.'
      ]
    }, {
      type: 'table',
      caption: 'Jarque-Bera 검정',
      columns: [
        { key: 'name', label: '변수' },
        { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'g1', label: '왜도(g1)', format: { type: 'fixed', digits: 3 } },
        { key: 'g2', label: '초과첨도(g2)', format: { type: 'fixed', digits: 3 } },
        { key: 'jb', label: 'JB 통계량', format: { type: 'fixed', digits: 3 } },
        { key: 'df', label: '자유도', format: { type: 'int' } },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
      ],
      rows: rows,
      notes: [
        '귀무가설: 변수가 정규분포를 따른다. 유의확률이 작으면(예: .05 미만) 정규성을 기각한다.',
        'JB는 모집단식 왜도·첨도(g1, g2)를 쓰므로 기술통계의 표본 보정 값(G1, G2)과 조금 다르다.',
        '점근 근사(χ², 자유도 2)라서 표본이 작을 때(대략 30 미만)는 p값이 부정확하다.'
      ]
    }].concat(verdicts),
    warnings: warnings
  };
}

/* ---------- Levene ---------- */

function median_(sorted) {
  var n = sorted.length;
  return n % 2 ? sorted[(n - 1) / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2;
}

/**
 * Levene 등분산 검정. groups: 집단별 숫자 배열, useMedian: 중앙값 기준(Brown-Forsythe) 여부.
 * Z_ij = |x_ij - 집단 중심|, W = (N-k)/(k-1) · Σ n_i(Z̄_i - Z̄)² / ΣΣ(Z_ij - Z̄_i)²,
 * 귀무가설에서 W ~ F(k-1, N-k). 집단 내 Z의 변동이 0이면 W는 null.
 */
function leveneTest(groups, useMedian) {
  var d = assumeDeps_();
  var k = groups.length;
  var N = 0;
  groups.forEach(function (g) { N += g.length; });
  var r = { W: null, df1: k - 1, df2: N - k, p: null };
  if (k < 2 || N - k < 1) return r;

  var zs = groups.map(function (g) {
    var center;
    if (useMedian) center = median_(g.slice().sort(function (a, b) { return a - b; }));
    else center = d.centralSums(g).mean;
    return g.map(function (x) { return Math.abs(x - center); });
  });
  var zMeans = zs.map(function (z) { return d.centralSums(z).mean; });
  var total = 0;
  zs.forEach(function (z, i) { z.forEach(function (v) { total += v; }); });
  var grand = total / N;

  var between = 0, within = 0;
  zs.forEach(function (z, i) {
    between += z.length * (zMeans[i] - grand) * (zMeans[i] - grand);
    z.forEach(function (v) { within += (v - zMeans[i]) * (v - zMeans[i]); });
  });
  if (!(within > 0)) return r;
  r.W = (N - k) / (k - 1) * between / within;
  r.p = d.pf(r.W, r.df1, r.df2, false);
  return r;
}

var LEVENE_MAX_GROUPS_ = 50;

function runLevene(dataset, request) {
  var d = assumeDeps_();
  var depVars = d.getVars(dataset, request.dependent, { numeric: true });
  var factorVars = d.getVars(dataset, request.factor);
  if (factorVars.length !== 1) throw new Error('집단변수는 하나만 선택하세요.');
  var factor = factorVars[0];
  depVars.forEach(function (v) {
    if (v.name === factor.name) throw new Error('종속변수와 집단변수에 같은 변수 "' + v.name + '"를 쓸 수 없습니다.');
  });
  var useMedian = !!(request.options && request.options.median);
  var factorName = factor.label ? factor.label : factor.name;
  var warnings = [];
  var groupRows = [];
  var testRows = [];
  var levItems = [];

  depVars.forEach(function (v) {
    var name = v.label ? v.label : v.name;
    var split = d.splitGroups(dataset, v.name, factor);
    var groups = split.groups, used = split.used;
    if (groups.length > LEVENE_MAX_GROUPS_) {
      throw new Error('집단변수 "' + factorName + '"의 집단이 ' + groups.length + '개입니다. 범주가 적은 변수를 선택하세요.');
    }

    groups.forEach(function (g) {
      var s = d.describeColumn(g.values);
      groupRows.push({ name: name, group: g.label, n: s.n, mean: s.mean, sd: s.sd });
      if (s.n < 2) warnings.push('"' + name + '": 집단 "' + g.label + '"의 사례가 1개뿐입니다.');
    });

    var r = leveneTest(groups.map(function (g) { return g.values; }), useMedian);
    if (r.W === null) {
      warnings.push('"' + name + '": 집단이 2개 미만이거나 집단 내 변동이 없어 검정할 수 없습니다.');
    }
    testRows.push({ name: name, N: used, W: r.W, df1: r.df1, df2: r.df2, p: r.p });
    levItems.push({ name: name, p: r.p });
  });

  return {
    analysis: 'levene',
    title: '등분산 검정',
    blocks: [
      {
        type: 'table',
        caption: '집단별 기술통계 (집단변수: ' + factorName + ')',
        columns: [
          { key: 'name', label: '종속변수' },
          { key: 'group', label: '집단' },
          { key: 'n', label: 'N', format: { type: 'int' } },
          { key: 'mean', label: '평균', format: { type: 'fixed', digits: 3 } },
          { key: 'sd', label: '표준편차', format: { type: 'fixed', digits: 3 } }
        ],
        rows: groupRows,
        notes: []
      },
      {
        type: 'table',
        caption: 'Levene의 등분산 검정',
        columns: [
          { key: 'name', label: '종속변수' },
          { key: 'N', label: 'N', format: { type: 'int' } },
          { key: 'W', label: 'Levene 통계량', format: { type: 'fixed', digits: 3 } },
          { key: 'df1', label: '자유도1', format: { type: 'int' } },
          { key: 'df2', label: '자유도2', format: { type: 'int' } },
          { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
        ],
        rows: testRows,
        notes: [
          '귀무가설: 집단 간 분산이 같다. 유의확률이 작으면(예: .05 미만) 등분산 가정을 기각한다.',
          useMedian
            ? '집단 중심: 중앙값 (Brown-Forsythe 변형, 정규가 아닌 분포에 더 안정적).'
            : '집단 중심: 평균 (원래의 Levene 검정).',
          '종속변수와 집단변수가 모두 유효한 사례만 사용했다.'
        ]
      }
    ].concat(advLeveneParagraphs_(levItems)),
    warnings: warnings
  };
}

/* ---------- 가정 점검과 비모수 검정 추천 ---------- */

var ADV_ALPHA_ = 0.05;
var ADV_LARGE_N_ = 100;      // 이 이상이면 정규성 검정이 사소한 이탈에도 기각하므로 왜도·첨도가 큰 경우에만 비모수를 권한다
var ADV_SMALL_N_ = 10;       // 이 미만이면 정규성을 기각하지 못해도 확인이 어렵다고 알린다
var ADV_SKEW_ = 2;           // West, Finch & Curran (1995): |왜도| < 2, |첨도| < 7 이면 정규성 이탈이 크지 않다고 본다
var ADV_KURT_ = 7;
var ADV_MAX_PARAGRAPHS_ = 10;

var ADV_STATUS_TEXT_ = {
  ok: '기각하지 못함',
  small: '기각하지 못함 (표본이 작아 확인이 어려움)',
  reject: '정규성 기각',
  largeMild: '기각 (표본이 커서 영향은 작을 수 있음)',
  untestable: '검정할 수 없음'
};

/** 유의확률 문구: ".043", "< .001" */
function advP_(p) {
  if (p === null || p === undefined) return '';
  return p < 0.001 ? '< .001' : p.toFixed(3).replace(/^0/, '');
}

/**
 * 한 변수(또는 집단)의 정규성 판정. 3~5000건은 Shapiro-Wilk, 그 밖은 Jarque-Bera.
 * status: ok | small | reject | largeMild | untestable
 */
function advNormality_(values) {
  var n = values.length;
  var jb = jarqueBera(values);
  var r = { n: n, test: null, p: null, status: 'untestable' };
  if (jb.jb === null) return r;
  if (n <= SW_MAX_N_) { r.test = 'Shapiro-Wilk'; r.p = shapiroWilk(values).p; }
  else { r.test = 'Jarque-Bera'; r.p = jb.p; }
  if (r.p >= ADV_ALPHA_) r.status = n < ADV_SMALL_N_ ? 'small' : 'ok';
  else if (n >= ADV_LARGE_N_ && Math.abs(jb.g1) < ADV_SKEW_ && Math.abs(jb.g2) < ADV_KURT_) r.status = 'largeMild';
  else r.status = 'reject';
  return r;
}

/**
 * 평균 비교 계열 한 변수의 점검. o = { name, checks: [{label, values}], levene: {p, n} 또는 undefined, alt, welchHint, ordinal }
 * 반환 { rows, recs }: 점검표 행과 추천 문구(문단).
 */
function advMeanTest_(o) {
  var rows = [], bad = [], mild = [], levBad = false;
  o.checks.forEach(function (c) {
    var r = advNormality_(c.values);
    var target = o.name + (c.label ? ' — ' + c.label : '');
    rows.push({ target: target, n: r.n, test: r.test || '', p: r.p, verdict: ADV_STATUS_TEXT_[r.status] });
    if (r.status === 'reject') bad.push(c.label || o.name);
    if (r.status === 'largeMild') mild.push(c.label || o.name);
  });
  if (o.levene !== undefined) {
    var lp = o.levene.p;
    levBad = lp !== null && lp < ADV_ALPHA_;
    rows.push({ target: o.name + ' — 집단 간 분산', n: o.levene.n, test: 'Levene', p: lp,
      verdict: lp === null ? '검정할 수 없음' : (levBad ? '등분산 기각' : '기각하지 못함') });
  }
  var recs = [];
  if (bad.length) {
    recs.push('"' + o.name + '": ' + bad.join(', ') + '의 정규성이 기각되었습니다. 평균 비교 결과와 함께 ' + o.alt + ' 결과도 확인하세요.' +
      (levBad ? ' (집단 간 분산도 같지 않습니다.)' : ''));
  } else if (levBad) {
    recs.push('"' + o.name + '": 집단 간 분산이 같지 않습니다. ' + o.welchHint);
  }
  if (!bad.length && mild.length) {
    recs.push('"' + o.name + '": ' + mild.join(', ') + '에서 정규성이 기각되었지만 표본이 크고(N ≥ ' + ADV_LARGE_N_ + ') 왜도·첨도가 크지 않아 평균 비교 결과를 그대로 써도 큰 문제는 없습니다. 히스토그램도 확인하세요.');
  }
  if (o.ordinal) {
    recs.push('"' + o.name + '"은(는) 서열 척도로 지정되어 있어 평균보다 순위를 쓰는 ' + o.alt + '이(가) 적절합니다.');
  }
  return { rows: rows, recs: recs };
}

/** 점검표와 추천 문단을 결과 블록으로 만든다. rows·recs 는 advMeanTest_ 결과를 이어 붙인 것. */
function advBlocks_(rows, recs, usedLevene) {
  var anySmall = rows.some(function (r) { return r.verdict === ADV_STATUS_TEXT_.small; });
  var notes = [
    '귀무가설: 정규분포를 따른다' + (usedLevene ? ' / 집단 간 분산이 같다' : '') + '. 유의확률이 .05 미만이면 기각한다. 정규성은 사례가 ' + SW_MAX_N_ + '건 이하이면 Shapiro-Wilk, 그 이상이면 Jarque-Bera로 본다.',
    '사례가 ' + ADV_LARGE_N_ + '건 이상이면 사소한 이탈도 기각되므로 |왜도| ≥ ' + ADV_SKEW_ + ' 또는 |첨도| ≥ ' + ADV_KURT_ + '처럼 이탈이 클 때만 비모수 검정을 권한다.'
  ];
  if (anySmall) notes.push('사례가 ' + ADV_SMALL_N_ + '건 미만이면 검정력이 낮아 정규성을 기각하지 못해도 정규라고 단정할 수 없다. 히스토그램이나 상자그림도 확인한다.');
  notes.push('이 점검은 판단을 돕는 안내이며 분석 결과를 바꾸지 않는다. 옵션 "가정 점검"을 끄면 이 표와 추천이 나오지 않는다.');
  var blocks = [{
    type: 'table',
    caption: '가정 점검 (정규성' + (usedLevene ? ', 등분산' : '') + ')',
    columns: [
      { key: 'target', label: '대상' },
      { key: 'n', label: 'N', format: { type: 'int' } },
      { key: 'test', label: '검정' },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
      { key: 'verdict', label: '판정' }
    ],
    rows: rows,
    notes: notes
  }];
  if (recs.length) {
    recs.slice(0, ADV_MAX_PARAGRAPHS_).forEach(function (t) { blocks.push({ type: 'text', text: '▶ 추천: ' + t }); });
    if (recs.length > ADV_MAX_PARAGRAPHS_) blocks.push({ type: 'text', text: '▶ 그 밖에 ' + (recs.length - ADV_MAX_PARAGRAPHS_) + '개 변수에 대한 추천은 점검표의 판정을 보세요.' });
  } else {
    blocks.push({ type: 'text', text: '▶ 정규성' + (usedLevene ? '·등분산' : '') + ' 가정을 기각하지 못해 이 검정 결과를 그대로 해석해도 됩니다.' });
  }
  return blocks;
}

var ADV_ALT_MW_ = 'Mann-Whitney U 검정';
var ADV_ALT_WILCOXON_PAIRED_ = 'Wilcoxon 부호순위 검정 (대응표본)';
var ADV_ALT_WILCOXON_ONE_ = 'Wilcoxon 부호순위 검정 (일표본)';
var ADV_ALT_KW_ = 'Kruskal-Wallis 검정';

/** 독립표본 t검정 한 종속변수: groups = [{label, values}] (2개) */
function advIndependentT_(name, groups, ordinal) {
  var lev = leveneTest(groups.map(function (g) { return g.values; }), false);
  var n = 0;
  groups.forEach(function (g) { n += g.values.length; });
  return advMeanTest_({
    name: name, checks: groups, levene: { p: lev.p, n: n }, ordinal: ordinal, alt: ADV_ALT_MW_,
    welchHint: '표의 "등분산을 가정하지 않음"(Welch) 행의 결과를 보세요.'
  });
}

/** 대응표본 t검정: 차이(변수1 − 변수2)의 정규성 */
function advPaired_(name, diffs, ordinal) {
  return advMeanTest_({ name: name, checks: [{ label: '차이(변수1 − 변수2)', values: diffs }], ordinal: ordinal, alt: ADV_ALT_WILCOXON_PAIRED_ });
}

/** 일표본 t검정: 변수 자체의 정규성 */
function advOneSample_(name, values, ordinal) {
  return advMeanTest_({ name: name, checks: [{ label: '', values: values }], ordinal: ordinal, alt: ADV_ALT_WILCOXON_ONE_ });
}

/** 일원 분산분석 한 종속변수: groups = [{label, values}] */
function advAnova_(name, groups, ordinal) {
  var lev = leveneTest(groups.map(function (g) { return g.values; }), false);
  var n = 0;
  groups.forEach(function (g) { n += g.values.length; });
  return advMeanTest_({
    name: name, checks: groups, levene: { p: lev.p, n: n }, ordinal: ordinal, alt: ADV_ALT_KW_,
    welchHint: 'Welch 검정과 Games-Howell 사후검정(옵션)의 결과를 보세요.'
  });
}

/** 피어슨 상관: vars = [{name, values, ordinal}] 각 변수의 정규성. spearmanOn: 스피어만도 같이 계산했는지 */
function advCorrelation_(vars, spearmanOn) {
  var rows = [], bad = [], mild = [], ord = [];
  vars.forEach(function (v) {
    var r = advNormality_(v.values);
    rows.push({ target: v.name, n: r.n, test: r.test || '', p: r.p, verdict: ADV_STATUS_TEXT_[r.status] });
    if (r.status === 'reject') bad.push(v.name);
    if (r.status === 'largeMild') mild.push(v.name);
    if (v.ordinal) ord.push(v.name);
  });
  var recs = [];
  var how = spearmanOn ? '스피어만 순위상관 결과를 우선 보세요.' : '옵션에서 스피어만 순위상관을 켜서 함께 확인하세요.';
  if (bad.length) recs.push(bad.map(function (s) { return '"' + s + '"'; }).join(', ') + '의 정규성이 기각되었습니다. 피어슨 상관은 정규성 이탈과 이상치에 민감합니다. ' + how);
  if (ord.length) recs.push(ord.map(function (s) { return '"' + s + '"'; }).join(', ') + '은(는) 서열 척도로 지정되어 있어 피어슨보다 순위를 쓰는 스피어만이 적절합니다. ' + how);
  if (!bad.length && mild.length) recs.push(mild.map(function (s) { return '"' + s + '"'; }).join(', ') + ': 정규성이 기각되었지만 표본이 크고 왜도·첨도가 크지 않아 피어슨 상관을 그대로 써도 큰 문제는 없습니다. 산점도로 이상치와 비선형 관계도 확인하세요.');
  return { rows: rows, recs: recs };
}

/** 정규성 검정 결과 아래의 판정 문단 (변수별) */
function advNormalityParagraphs_(items) {
  var blocks = [];
  items.slice(0, ADV_MAX_PARAGRAPHS_).forEach(function (it) {
    var r = advNormality_(it.values);
    var t;
    if (r.status === 'untestable') return;
    if (r.status === 'reject') {
      t = '"' + it.name + '": 정규성이 기각되었습니다 (' + r.test + ' p = ' + advP_(r.p) + '). 이 변수로 t검정·분산분석·피어슨 상관을 쓸 때는 비모수 대안(Mann-Whitney U, Wilcoxon 부호순위, Kruskal-Wallis, 스피어만 상관)을 함께 고려하세요.';
    } else if (r.status === 'largeMild') {
      t = '"' + it.name + '": 정규성이 기각되었지만 (' + r.test + ' p = ' + advP_(r.p) + ') 표본이 크고(N ≥ ' + ADV_LARGE_N_ + ') 왜도·첨도가 크지 않아 영향은 작을 수 있습니다. 히스토그램도 확인하세요.';
    } else {
      t = '"' + it.name + '": 정규성을 기각하지 못했습니다 (' + r.test + ' p = ' + advP_(r.p) + ').' + (r.status === 'small' ? ' 다만 사례가 ' + ADV_SMALL_N_ + '건 미만이라 정규라고 단정할 수는 없습니다.' : '');
    }
    blocks.push({ type: 'text', text: '▶ ' + t });
  });
  if (items.length > ADV_MAX_PARAGRAPHS_) blocks.push({ type: 'text', text: '▶ 그 밖의 변수는 표의 유의확률을 보세요.' });
  return blocks;
}

/** 등분산 검정 결과 아래의 판정 문단: items = [{name, p}] */
function advLeveneParagraphs_(items) {
  var blocks = [];
  items.slice(0, ADV_MAX_PARAGRAPHS_).forEach(function (it) {
    if (it.p === null) return;
    var t;
    if (it.p < ADV_ALPHA_) {
      t = '"' + it.name + '": 집단 간 분산이 같지 않습니다 (p = ' + advP_(it.p) + '). 평균 비교에는 Welch 검정(t검정의 "등분산을 가정하지 않음" 행, 분산분석의 Welch 옵션)과 Games-Howell 사후검정을 쓰거나 Mann-Whitney U, Kruskal-Wallis 같은 비모수 검정을 고려하세요. 비모수 검정도 집단 분포의 모양이 비슷하다는 전제가 있습니다.';
    } else {
      t = '"' + it.name + '": 등분산 가정을 기각하지 못했습니다 (p = ' + advP_(it.p) + ').';
    }
    blocks.push({ type: 'text', text: '▶ ' + t });
  });
  if (items.length > ADV_MAX_PARAGRAPHS_) blocks.push({ type: 'text', text: '▶ 그 밖의 변수는 표의 유의확률을 보세요.' });
  return blocks;
}


/** 점검 결과 모으기: adv = { rows: [], recs: [] } 에 변수별 결과를 더하고, 끝에서 블록으로 만든다 */
function advAdd_(adv, part) {
  adv.rows = adv.rows.concat(part.rows);
  adv.recs = adv.recs.concat(part.recs);
}
function advAppend_(blocks, adv, usedLevene) {
  return blocks.concat(advBlocks_(adv.rows, adv.recs, usedLevene));
}

if (typeof module !== 'undefined') {
  module.exports = {
    advNormality_: advNormality_, advMeanTest_: advMeanTest_, advBlocks_: advBlocks_, advAdd_: advAdd_, advAppend_: advAppend_,
    advIndependentT_: advIndependentT_, advPaired_: advPaired_, advOneSample_: advOneSample_, advAnova_: advAnova_,
    advCorrelation_: advCorrelation_, advNormalityParagraphs_: advNormalityParagraphs_, advLeveneParagraphs_: advLeveneParagraphs_,
    jarqueBera: jarqueBera,
    shapiroWilk: shapiroWilk,
    leveneTest: leveneTest,
    runNormality: runNormality,
    runLevene: runLevene
  };
}
