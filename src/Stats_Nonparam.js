/**
 * 비모수 검정: Mann-Whitney U, Wilcoxon 부호순위(대응표본·일표본), Kruskal-Wallis + Dunn 사후검정, Friedman + Dunn–Bonferroni 사후비교.
 * 순수 함수: (dataset, request) → 결과 객체. 결과표 구성: 순위표 + 검정 통계량표.
 *
 * - 점근 p값은 동순위(ties) 보정을 하고 연속성 보정은 하지 않는다.
 * - 정확 p값: 동순위가 없을 때만 계산한다. Mann-Whitney는 각 집단 50건 이하, Wilcoxon은 짝 50건 이하.
 */

var NP_EXACT_MAX_ = 50;       // 정확 p값을 계산하는 표본 크기의 상한
var NP_MAX_GROUPS_ = 50;

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function npDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return {
      getVars: data.getVars, splitGroups: data.splitGroups, pnorm: dist.pnorm, pchisq: dist.pchisq,
      boxplotBlock: require('./Charts.js').boxplotBlock
    };
  }
  return { getVars: getVars, splitGroups: splitGroups, pnorm: pnorm, pchisq: pchisq, boxplotBlock: boxplotBlock };
}

var NP_F3_ = { type: 'fixed', digits: 3 };
function npName_(v) { return v.label ? v.label : v.name; }

/* ---------- 순위와 중앙값 ---------- */

/** 평균순위와 동순위 항 Σ(t³ − t) (t: 같은 값의 개수) */
function npRanks_(values) {
  var n = values.length;
  var idx = new Array(n), i;
  for (i = 0; i < n; i++) idx[i] = i;
  idx.sort(function (a, b) { return values[a] - values[b]; });
  var ranks = new Array(n), tieSum = 0;
  i = 0;
  while (i < n) {
    var j = i;
    while (j + 1 < n && values[idx[j + 1]] === values[idx[i]]) j++;
    var avg = (i + j) / 2 + 1;
    for (var k = i; k <= j; k++) ranks[idx[k]] = avg;
    var t = j - i + 1;
    tieSum += t * t * t - t;
    i = j + 1;
  }
  return { ranks: ranks, tieSum: tieSum };
}

function npMedian_(values) {
  var n = values.length;
  if (n === 0) return null;
  var s = values.slice().sort(function (a, b) { return a - b; });
  return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2;
}

/** 양측 p값: 2 · P(Z > |z|), 1을 넘지 않게 */
function npTwoSided_(z) {
  return Math.min(1, 2 * npDeps_().pnorm(Math.abs(z), false));
}

/* ---------- Mann-Whitney U ---------- */

/**
 * 정확 분포: 동순위가 없을 때 P(U ≤ u). 순위 0..N-1 중 m개를 고르는 방법의 수를 합(S)별로 센다 (U = S - m(m-1)/2).
 * 덧셈만 쓰므로 반올림 오차가 거의 없다.
 */
function mwExactCdf_(m, n, u) {
  var N = m + n;
  var cap = u + m * (m - 1) / 2;
  var dp = [], k, s, v;
  for (k = 0; k <= m; k++) { dp.push(new Array(cap + 1)); for (s = 0; s <= cap; s++) dp[k][s] = 0; }
  dp[0][0] = 1;
  for (v = 0; v < N; v++) {
    for (k = Math.min(m, v + 1); k >= 1; k--) {
      for (s = cap; s >= v; s--) dp[k][s] += dp[k - 1][s - v];
    }
  }
  var count = 0;
  for (s = 0; s <= cap; s++) count += dp[m][s];
  var total = 1;                                   // C(N, m)
  for (k = 1; k <= m; k++) total = total * (n + k) / k;
  return count / total;
}

/**
 * a, b: 두 집단의 값. U는 두 U 중 작은 값, Z는 그 U 기준(≤ 0), W는 작은 집단(같으면 집단 1)의 순위합.
 * 계산할 수 없으면 해당 값은 null.
 */
function mannWhitney(a, b) {
  var d = npDeps_();
  var n1 = a.length, n2 = b.length, N = n1 + n2;
  var r = { n1: n1, n2: n2, u: null, u1: null, u2: null, w: null, z: null, p: null, pExact: null, effect: null,
    rankSum1: null, rankSum2: null, meanRank1: null, meanRank2: null, ties: false };
  if (n1 < 1 || n2 < 1) return r;
  var rk = npRanks_(a.concat(b));
  var rs1 = 0, i;
  for (i = 0; i < n1; i++) rs1 += rk.ranks[i];
  var rs2 = N * (N + 1) / 2 - rs1;
  r.rankSum1 = rs1; r.rankSum2 = rs2;
  r.meanRank1 = rs1 / n1; r.meanRank2 = rs2 / n2;
  r.u1 = rs1 - n1 * (n1 + 1) / 2;
  r.u2 = n1 * n2 - r.u1;
  r.u = Math.min(r.u1, r.u2);
  r.w = n1 <= n2 ? rs1 : rs2;
  r.ties = rk.tieSum > 0;
  var varU = n1 * n2 / 12 * ((N + 1) - rk.tieSum / (N * (N - 1)));
  if (N > 1 && varU > 0) {
    r.z = (r.u - n1 * n2 / 2) / Math.sqrt(varU);
    r.p = npTwoSided_(r.z);
    r.effect = Math.abs(r.z) / Math.sqrt(N);
  }
  if (!r.ties && n1 <= NP_EXACT_MAX_ && n2 <= NP_EXACT_MAX_) {
    r.pExact = Math.min(1, 2 * mwExactCdf_(n1, n2, r.u));
  }
  return r;
}

/* ---------- Wilcoxon 부호순위 ---------- */

/** 정확 분포: 동순위·0이 없을 때 P(W ≤ w), 순위 1..n 부분집합의 합 */
function signedRankExactCdf_(n, w) {
  var max = n * (n + 1) / 2, dp = new Array(max + 1), s, k;
  for (s = 0; s <= max; s++) dp[s] = 0;
  dp[0] = 1;
  for (k = 1; k <= n; k++) for (s = max; s >= k; s--) dp[s] += dp[s - k];
  var count = 0;
  for (s = 0; s <= w && s <= max; s++) count += dp[s];
  return count / Math.pow(2, n);
}

/**
 * diffs: 차이 값(변수1 − 변수2, 또는 값 − 검정값). 0인 차이는 순위에서 뺀다 (zero_method 'wilcox').
 * Z = (W+ − 기댓값) / 표준편차, W+ 는 양의 차이의 순위합 (양수면 양의 차이가 우세).
 */
function signedRank(diffs) {
  var nonzero = diffs.filter(function (x) { return x !== 0; });
  var n = nonzero.length;
  var r = { total: diffs.length, n: n, nZero: diffs.length - n, nPos: 0, nNeg: 0, wPlus: null, wMinus: null,
    meanPos: null, meanNeg: null, z: null, p: null, pExact: null, effect: null, ties: false };
  if (n === 0) return r;
  var rk = npRanks_(nonzero.map(function (x) { return Math.abs(x); }));
  var wPlus = 0, wMinus = 0, i;
  for (i = 0; i < n; i++) {
    if (nonzero[i] > 0) { wPlus += rk.ranks[i]; r.nPos++; } else { wMinus += rk.ranks[i]; r.nNeg++; }
  }
  r.wPlus = wPlus; r.wMinus = wMinus;
  r.meanPos = r.nPos ? wPlus / r.nPos : null;
  r.meanNeg = r.nNeg ? wMinus / r.nNeg : null;
  r.ties = rk.tieSum > 0;
  var variance = n * (n + 1) * (2 * n + 1) / 24 - rk.tieSum / 48;
  if (variance > 0) {
    r.z = (wPlus - n * (n + 1) / 4) / Math.sqrt(variance);
    r.p = npTwoSided_(r.z);
    r.effect = Math.abs(r.z) / Math.sqrt(n);
  }
  if (!r.ties && n <= NP_EXACT_MAX_) {
    r.pExact = Math.min(1, 2 * signedRankExactCdf_(n, Math.min(wPlus, wMinus)));
  }
  return r;
}

/* ---------- Kruskal-Wallis와 Dunn ---------- */

/** groups: 집단별 값 배열. H는 동순위 보정, p는 χ²(k−1) 근사. η²_H = (H − k + 1)/(N − k) (음수면 0). */
function kruskalWallis(groups) {
  var d = npDeps_();
  var k = groups.length, all = [], ns = [], i, g;
  groups.forEach(function (gr) { ns.push(gr.length); all = all.concat(gr); });
  var N = all.length;
  var r = { k: k, N: N, ns: ns, rankSums: [], meanRanks: [], h: null, df: k - 1, p: null, eta2: null, tieSum: 0 };
  if (k < 2 || N < 2) return r;
  var rk = npRanks_(all);
  r.tieSum = rk.tieSum;
  var pos = 0, sumSq = 0;
  for (g = 0; g < k; g++) {
    var rs = 0;
    for (i = 0; i < ns[g]; i++) rs += rk.ranks[pos + i];
    pos += ns[g];
    r.rankSums.push(rs);
    r.meanRanks.push(ns[g] ? rs / ns[g] : null);
    if (ns[g]) sumSq += rs * rs / ns[g];
  }
  var h0 = 12 / (N * (N + 1)) * sumSq - 3 * (N + 1);
  var corr = 1 - rk.tieSum / (N * N * N - N);
  if (corr > 0) {
    r.h = Math.max(0, h0 / corr);
    r.p = d.pchisq(r.h, k - 1, false);
    if (N > k) r.eta2 = Math.max(0, (r.h - k + 1) / (N - k));
  }
  return r;
}

/**
 * Dunn 사후검정 (모든 쌍, Bonferroni 보정). z = (평균순위 I − 평균순위 J) / sqrt((N(N+1)/12 − Σ(t³−t)/(12(N−1))) · (1/nI + 1/nJ)).
 * 반환: [{ i, j, diff, z, p, pAdj }] (i < j). 계산할 수 없으면 z, p, pAdj 는 null.
 */
function dunnTest(groups) {
  var d = npDeps_();
  var kw = kruskalWallis(groups);
  var k = groups.length, N = kw.N;
  var pairs = [];
  var base = N * (N + 1) / 12 - kw.tieSum / (12 * (N - 1));
  var m = k * (k - 1) / 2;
  for (var i = 0; i < k; i++) {
    for (var j = i + 1; j < k; j++) {
      var row = { i: i, j: j, diff: null, z: null, p: null, pAdj: null };
      if (kw.ns[i] > 0 && kw.ns[j] > 0) {
        row.diff = kw.meanRanks[i] - kw.meanRanks[j];
        var se2 = base * (1 / kw.ns[i] + 1 / kw.ns[j]);
        if (N > 1 && se2 > 0) {
          row.z = row.diff / Math.sqrt(se2);
          row.p = Math.min(1, 2 * d.pnorm(Math.abs(row.z), false));
          row.pAdj = Math.min(1, row.p * m);
        }
      }
      pairs.push(row);
    }
  }
  return pairs;
}

/* ---------- 결과: 공통 ---------- */

function npPlainNumber_(options, key, label, def) {
  var v = options ? options[key] : undefined;
  if (v === undefined || v === null || v === '') return def;
  var x = Number(v);
  if (!isFinite(x)) throw new Error('"' + label + '"은(는) 숫자여야 합니다.');
  return x;
}

function npPvalueColumns_() {
  return [
    { key: 'p', label: '점근 유의확률(양측)', format: { type: 'pvalue' } },
    { key: 'pExact', label: '정확 유의확률(양측)', format: { type: 'pvalue' } }
  ];
}

var NP_ASYMPTOTIC_NOTE_ = '점근 유의확률은 동순위 보정을 하고 연속성 보정은 하지 않는다.';

/** 두 집단 검정(독립표본) 입력 검사: 종속변수 여러 개 + 집단변수 하나 */
function npDepAndFactor_(dataset, request, d) {
  var depVars = d.getVars(dataset, request.dependent, { numeric: true });
  var factorVars = d.getVars(dataset, request.factor);
  if (factorVars.length !== 1) throw new Error('집단변수는 하나만 선택하세요.');
  var factor = factorVars[0];
  depVars.forEach(function (v) {
    if (v.name === factor.name) throw new Error('종속변수와 집단변수에 같은 변수 "' + v.name + '"를 쓸 수 없습니다.');
  });
  return { depVars: depVars, factor: factor, factorName: npName_(factor) };
}

/* ---------- 결과: Mann-Whitney U ---------- */

function runMannWhitney(dataset, request) {
  var d = npDeps_();
  var inp = npDepAndFactor_(dataset, request, d);
  var opt = request.options || {};
  var warnings = [], rankRows = [], testRows = [], boxBlocks = [], groupNote = [];

  inp.depVars.forEach(function (v) {
    var name = npName_(v);
    var split = d.splitGroups(dataset, v.name, inp.factor);
    if (split.groups.length !== 2) {
      throw new Error('집단변수 "' + inp.factorName + '"의 집단이 ' + split.groups.length + '개입니다 (종속변수 "' + name + '"에 유효한 사례 기준). Mann-Whitney U 검정은 집단이 정확히 2개여야 합니다. 3개 이상이면 Kruskal-Wallis 검정을 쓰세요.');
    }
    var g1 = split.groups[0], g2 = split.groups[1];
    var r = mannWhitney(g1.values, g2.values);
    if (opt.boxplot) {
      boxBlocks.push(d.boxplotBlock('상자그림: ' + name + ' (집단변수: ' + inp.factorName + ')',
        [{ label: g1.label, values: g1.values }, { label: g2.label, values: g2.values }], name));
    }
    rankRows.push({ name: name, group: g1.label, n: r.n1, median: npMedian_(g1.values), meanRank: r.meanRank1, rankSum: r.rankSum1 });
    rankRows.push({ name: '', group: g2.label, n: r.n2, median: npMedian_(g2.values), meanRank: r.meanRank2, rankSum: r.rankSum2 });
    rankRows.push({ _kind: 'total', name: '', group: '합계', n: r.n1 + r.n2, median: null, meanRank: null, rankSum: null });
    if (r.z === null) warnings.push('"' + name + '": 모든 값이 같아 검정할 수 없습니다.');
    testRows.push({ name: name, u: r.u, w: r.w, z: r.z, p: r.p, pExact: r.pExact, effect: r.effect });
    var note = '집단1: ' + g1.label + ', 집단2: ' + g2.label;
    if (groupNote.indexOf(note) < 0) groupNote.push(note);
  });

  var testCols = [
    { key: 'name', label: '종속변수' },
    { key: 'u', label: 'Mann-Whitney U', format: { type: 'trim', digits: 3 } },
    { key: 'w', label: 'Wilcoxon W', format: { type: 'trim', digits: 3 } },
    { key: 'z', label: 'Z', format: NP_F3_ }
  ].concat(npPvalueColumns_());
  if (opt.effect) testCols.push({ key: 'effect', label: '효과크기 r', format: NP_F3_ });

  var blocks = [
    {
      type: 'table',
      caption: '순위 (집단변수: ' + inp.factorName + ')',
      columns: [
        { key: 'name', label: '종속변수' }, { key: 'group', label: '집단' },
        { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'median', label: '중앙값', format: { type: 'trim', digits: 3 } },
        { key: 'meanRank', label: '평균 순위', format: { type: 'fixed', digits: 2 } },
        { key: 'rankSum', label: '순위합', format: { type: 'trim', digits: 2 } }
      ],
      rows: rankRows,
      notes: ['종속변수와 집단변수가 모두 유효한 사례만 사용했다. 평균 순위가 큰 집단의 값이 전반적으로 크다.']
    },
    {
      type: 'table',
      caption: 'Mann-Whitney U 검정 통계량',
      columns: testCols,
      rows: testRows,
      notes: [
        '귀무가설: 두 집단의 분포가 같다 (두 집단의 값이 어느 쪽이 더 크다고 할 수 없다). 유의확률이 작으면(예: .05 미만) 기각한다.',
        'U는 두 U 중 작은 값이고 Z는 그 U 기준이라 0 이하이다. Wilcoxon W는 사례가 적은 집단(같으면 집단1)의 순위합이다. ' + groupNote.join(' / ') + '.',
        NP_ASYMPTOTIC_NOTE_,
        '정확 유의확률은 동순위가 없고 각 집단이 ' + NP_EXACT_MAX_ + '건 이하일 때만 계산한다.'
      ].concat(opt.effect ? ['효과크기 r = |Z| / √N (대략 0.1 작음, 0.3 중간, 0.5 큼).'] : [])
    }
  ];
  return { analysis: 'mannWhitney', title: 'Mann-Whitney U 검정', blocks: blocks.concat(boxBlocks), warnings: warnings };
}

/* ---------- 결과: Wilcoxon 부호순위 ---------- */

function npSignedRankBlocks_(opts) {
  // opts: { rows (순위표 행), testRows, labels: { neg, pos, tie }, firstCaption, notes, effect }
  var testCols = [
    { key: 'name', label: opts.nameLabel },
    { key: 'z', label: 'Z', format: NP_F3_ }
  ].concat(npPvalueColumns_());
  if (opts.effect) testCols.push({ key: 'effect', label: '효과크기 r', format: NP_F3_ });
  return [
    {
      type: 'table',
      caption: '순위',
      columns: [
        { key: 'name', label: opts.nameLabel }, { key: 'kind', label: '' },
        { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'meanRank', label: '평균 순위', format: { type: 'fixed', digits: 2 } },
        { key: 'rankSum', label: '순위합', format: { type: 'trim', digits: 2 } }
      ],
      rows: opts.rankRows,
      notes: opts.rankNotes
    },
    {
      type: 'table',
      caption: 'Wilcoxon 부호순위 검정 통계량',
      columns: testCols,
      rows: opts.testRows,
      notes: [
        '귀무가설: 차이의 분포가 0을 중심으로 대칭이다 (중앙값 차이가 0). 유의확률이 작으면(예: .05 미만) 기각한다.',
        'Z는 양의 순위합 기준이다. 양수이면 양의 차이가 우세하고 음수이면 음의 차이가 우세하다 (유의확률은 방향과 무관).',
        NP_ASYMPTOTIC_NOTE_,
        '차이가 0인 사례는 순위에서 뺀다. 정확 유의확률은 동순위·0이 없고 순위를 매긴 사례가 ' + NP_EXACT_MAX_ + '건 이하일 때만 계산한다.'
      ].concat(opts.effect ? ['효과크기 r = |Z| / √(순위를 매긴 사례 수) (대략 0.1 작음, 0.3 중간, 0.5 큼).'] : [])
    }
  ];
}

function npRankRows_(name, labels, r) {
  return [
    { name: name, kind: labels.neg, n: r.nNeg, meanRank: r.meanNeg, rankSum: r.wMinus },
    { name: '', kind: labels.pos, n: r.nPos, meanRank: r.meanPos, rankSum: r.wPlus },
    { name: '', kind: labels.tie, n: r.nZero, meanRank: null, rankSum: null },
    { _kind: 'total', name: '', kind: '합계', n: r.total, meanRank: null, rankSum: null }
  ];
}

function runWilcoxonPaired(dataset, request) {
  var d = npDeps_();
  var first = d.getVars(dataset, request.first, { numeric: true });
  var second = d.getVars(dataset, request.second, { numeric: true });
  if (first.length !== 1 || second.length !== 1) throw new Error('변수1과 변수2에 변수를 하나씩 선택하세요.');
  var v1 = first[0], v2 = second[0];
  if (v1.name === v2.name) throw new Error('변수1과 변수2에 같은 변수 "' + v1.name + '"를 쓸 수 없습니다.');
  var opt = request.options || {};
  var n1 = npName_(v1), n2 = npName_(v2);

  var x = [], y = [], diffs = [];
  for (var i = 0; i < dataset.n; i++) {
    var a = dataset.data[v1.name][i], b = dataset.data[v2.name][i];
    if (a !== null && b !== null) { x.push(a); y.push(b); diffs.push(a - b); }
  }
  var r = signedRank(diffs);
  var warnings = [];
  if (r.z === null) warnings.push('순위를 매길 수 있는 쌍이 없거나(모두 차이가 0) 변동이 없어 검정할 수 없습니다.');
  var dropped = dataset.n - x.length;
  var pairName = n1 + ' − ' + n2;

  var blocks = [
    {
      type: 'table',
      caption: '기술통계',
      columns: [
        { key: 'name', label: '변수' }, { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'median', label: '중앙값', format: { type: 'trim', digits: 3 } }
      ],
      rows: [{ name: n1, n: x.length, median: npMedian_(x) }, { name: n2, n: y.length, median: npMedian_(y) }],
      notes: ['두 변수가 모두 유효한 쌍만 사용했다' + (dropped > 0 ? ' (' + dropped + '건 제외).' : '.')]
    }
  ].concat(npSignedRankBlocks_({
    nameLabel: '쌍',
    rankRows: npRankRows_(pairName, { neg: '음의 순위 (' + n1 + ' < ' + n2 + ')', pos: '양의 순위 (' + n1 + ' > ' + n2 + ')', tie: '동순위 (' + n1 + ' = ' + n2 + ')' }, r),
    rankNotes: ['차이 = 변수1 − 변수2. 순위는 차이의 절댓값으로 매긴다.'],
    testRows: [{ name: pairName, z: r.z, p: r.p, pExact: r.pExact, effect: r.effect }],
    effect: !!opt.effect
  }));
  return { analysis: 'wilcoxonPaired', title: 'Wilcoxon 부호순위 검정 (대응표본)', blocks: blocks, warnings: warnings };
}

function runWilcoxonOne(dataset, request) {
  var d = npDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  var opt = request.options || {};
  var mu = npPlainNumber_(opt, 'testValue', '검정값', 0);
  var warnings = [], descRows = [], rankRows = [], testRows = [];

  vars.forEach(function (v) {
    var name = npName_(v);
    var values = dataset.data[v.name].filter(function (x) { return x !== null; });
    var r = signedRank(values.map(function (x) { return x - mu; }));
    if (r.z === null) warnings.push('"' + name + '": 검정값과 다른 값이 없거나 변동이 없어 검정할 수 없습니다.');
    descRows.push({ name: name, n: values.length, median: npMedian_(values) });
    rankRows = rankRows.concat(npRankRows_(name, { neg: '검정값보다 작음', pos: '검정값보다 큼', tie: '검정값과 같음' }, r));
    testRows.push({ name: name, z: r.z, p: r.p, pExact: r.pExact, effect: r.effect });
  });

  var blocks = [
    {
      type: 'table',
      caption: '기술통계',
      columns: [
        { key: 'name', label: '변수' }, { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'median', label: '중앙값', format: { type: 'trim', digits: 3 } }
      ],
      rows: descRows,
      notes: []
    }
  ].concat(npSignedRankBlocks_({
    nameLabel: '변수',
    rankRows: rankRows,
    rankNotes: ['차이 = 값 − 검정값 (' + mu + '). 순위는 차이의 절댓값으로 매긴다.'],
    testRows: testRows,
    effect: !!opt.effect
  }));
  blocks[2].caption = 'Wilcoxon 부호순위 검정 통계량 (검정값 = ' + mu + ')';
  blocks[2].notes[0] = '귀무가설: 모집단의 중앙값이 검정값과 같다 (분포가 대칭이라고 가정). 유의확률이 작으면(예: .05 미만) 기각한다.';
  return { analysis: 'wilcoxonOne', title: 'Wilcoxon 부호순위 검정 (일표본)', blocks: blocks, warnings: warnings };
}

/* ---------- 결과: Kruskal-Wallis ---------- */

function runKruskalWallis(dataset, request) {
  var d = npDeps_();
  var inp = npDepAndFactor_(dataset, request, d);
  var opt = request.options || {};
  var warnings = [], rankRows = [], testRows = [], dunnRows = [], boxBlocks = [];

  inp.depVars.forEach(function (v) {
    var name = npName_(v);
    var split = d.splitGroups(dataset, v.name, inp.factor);
    var groups = split.groups;
    if (groups.length < 2) {
      throw new Error('집단변수 "' + inp.factorName + '"의 집단이 ' + groups.length + '개입니다 (종속변수 "' + name + '"에 유효한 사례 기준). 집단이 2개 이상이어야 합니다.');
    }
    if (groups.length > NP_MAX_GROUPS_) {
      throw new Error('집단변수 "' + inp.factorName + '"의 집단이 ' + groups.length + '개입니다 (최대 ' + NP_MAX_GROUPS_ + '개). 범주가 적은 변수를 선택하세요.');
    }
    var values = groups.map(function (g) { return g.values; });
    var r = kruskalWallis(values);
    if (opt.boxplot) {
      boxBlocks.push(d.boxplotBlock('상자그림: ' + name + ' (집단변수: ' + inp.factorName + ')',
        groups.map(function (g) { return { label: g.label, values: g.values }; }), name));
    }
    groups.forEach(function (g, idx) {
      rankRows.push({ name: idx === 0 ? name : '', group: g.label, n: r.ns[idx], median: npMedian_(g.values), meanRank: r.meanRanks[idx] });
    });
    rankRows.push({ _kind: 'total', name: '', group: '합계', n: r.N, median: null, meanRank: null });
    if (r.h === null) warnings.push('"' + name + '": 모든 값이 같아 검정할 수 없습니다.');
    testRows.push({ name: name, h: r.h, df: r.df, p: r.p, eta2: r.eta2 });

    if (opt.dunn) {
      var pairs = dunnTest(values);
      pairs.forEach(function (pr, idx) {
        dunnRows.push({
          name: idx === 0 ? name : '', gi: groups[pr.i].label, gj: groups[pr.j].label,
          diff: pr.diff, z: pr.z, p: pr.p, pAdj: pr.pAdj, sig: pr.pAdj !== null && pr.pAdj < 0.05 ? '*' : ''
        });
      });
    }
  });

  var testCols = [
    { key: 'name', label: '종속변수' },
    { key: 'h', label: 'Kruskal-Wallis H', format: NP_F3_ },
    { key: 'df', label: '자유도', format: { type: 'int' } },
    { key: 'p', label: '점근 유의확률', format: { type: 'pvalue' } }
  ];
  if (opt.effect) testCols.push({ key: 'eta2', label: '효과크기 η²(H)', format: NP_F3_ });

  var blocks = [
    {
      type: 'table',
      caption: '순위 (집단변수: ' + inp.factorName + ')',
      columns: [
        { key: 'name', label: '종속변수' }, { key: 'group', label: '집단' },
        { key: 'n', label: 'N', format: { type: 'int' } },
        { key: 'median', label: '중앙값', format: { type: 'trim', digits: 3 } },
        { key: 'meanRank', label: '평균 순위', format: { type: 'fixed', digits: 2 } }
      ],
      rows: rankRows,
      notes: ['종속변수와 집단변수가 모두 유효한 사례만 사용했다.']
    },
    {
      type: 'table',
      caption: 'Kruskal-Wallis 검정 통계량',
      columns: testCols,
      rows: testRows,
      notes: [
        '귀무가설: 모든 집단의 분포가 같다. 유의확률이 작으면(예: .05 미만) 어느 한 집단 이상이 다르다는 뜻이며, 어느 집단 사이인지는 사후검정으로 본다.',
        'H는 동순위 보정을 했고, 유의확률은 χ² 근사(자유도 = 집단 수 − 1)이다. 집단이 작으면(대략 5건 미만) 근사가 부정확하다.'
      ].concat(opt.effect ? ['효과크기 η²(H) = (H − k + 1) / (N − k), 음수이면 0 (대략 0.01 작음, 0.06 중간, 0.14 큼).'] : [])
    }
  ];
  if (opt.dunn) {
    blocks.push({
      type: 'table',
      caption: '다중비교: Dunn 사후검정 (Bonferroni 보정)',
      columns: [
        { key: 'name', label: '종속변수' },
        { key: 'gi', label: '(I) 집단' }, { key: 'gj', label: '(J) 집단' },
        { key: 'diff', label: '평균 순위 차이(I−J)', format: { type: 'fixed', digits: 2 } },
        { key: 'sig', label: '', align: 'left' },
        { key: 'z', label: 'Z', format: NP_F3_ },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
        { key: 'pAdj', label: '보정 유의확률', format: { type: 'pvalue' } }
      ],
      rows: dunnRows,
      notes: [
        'Dunn(1964)의 순위합 기반 쌍별 비교이며 동순위 보정을 했다. 보정 유의확률 = 유의확률 × 비교 횟수(최대 1)이다.',
        '모든 쌍을 한 번씩(집단변수 값 순으로 I < J) 표시했다. * 는 보정 유의확률이 .05 미만인 쌍이다. 전체 검정(H)이 유의하지 않으면 사후검정은 참고만 한다.'
      ]
    });
  }
  return { analysis: 'kruskalWallis', title: 'Kruskal-Wallis 검정', blocks: blocks.concat(boxBlocks), warnings: warnings };
}

/* ---------- Friedman 검정 (같은 사람의 3회 이상 반복 측정, 비모수) ---------- */

var NP_FRIEDMAN_MAX_VARS_ = 20;

/**
 * Friedman 검정. columns: 측정(변수)별 값 배열, 모두 같은 길이 n (사람마다 하나씩, 결측 없음).
 * 사람(행)마다 값에 순위를 매기고(동순위는 평균 순위) 측정별 순위 합 R_j 로
 * χ² = 12·ΣR_j² / (n·k·(k + 1)) − 3·n·(k + 1) 을 동순위 보정 계수 C = 1 − ΣΣ(t³ − t) / (n·k·(k² − 1)) 로 나눈다 (자유도 k − 1, χ² 근사).
 * Kendall W = χ² / (n·(k − 1)). 0/1 자료에서는 Cochran Q 와 같다.
 * 모든 사람의 값이 전부 같아 C = 0 이면 chi2, p, w 는 null.
 */
function friedman(columns) {
  var d = npDeps_();
  var k = columns.length, n = k ? columns[0].length : 0, i, j;
  var r = { n: n, k: k, rankSums: [], meanRanks: [], chi2: null, df: k - 1, p: null, w: null, tieTotal: 0, corr: null };
  if (k < 2 || n < 1) return r;
  var sums = new Array(k).fill(0), tieTotal = 0;
  for (i = 0; i < n; i++) {
    var row = new Array(k);
    for (j = 0; j < k; j++) row[j] = columns[j][i];
    var rk = npRanks_(row);
    for (j = 0; j < k; j++) sums[j] += rk.ranks[j];
    tieTotal += rk.tieSum;
  }
  var sumSq = 0;
  for (j = 0; j < k; j++) { r.rankSums.push(sums[j]); r.meanRanks.push(sums[j] / n); sumSq += sums[j] * sums[j]; }
  r.tieTotal = tieTotal;
  r.corr = 1 - tieTotal / (n * k * (k * k - 1));
  if (r.corr > 0) {
    var chi0 = 12 * sumSq / (n * k * (k + 1)) - 3 * n * (k + 1);
    r.chi2 = Math.max(0, chi0 / r.corr);
    r.p = d.pchisq(r.chi2, k - 1, false);
    r.w = r.chi2 / (n * (k - 1));
  }
  return r;
}

/**
 * Friedman 사후비교: Dunn–Bonferroni. z = (평균 순위 I − 평균 순위 J) / √(k(k + 1)/(6n) · C), C 는 동순위 보정 계수.
 * 반환: [{ i, j, diff, z, p, pAdj }] (i < j). 검정할 수 없으면 z, p, pAdj 는 null.
 */
function friedmanPairs(fr) {
  var d = npDeps_(), k = fr.k, n = fr.n, pairs = [], m = k * (k - 1) / 2;
  for (var i = 0; i < k; i++) {
    for (var j = i + 1; j < k; j++) {
      var row = { i: i, j: j, diff: fr.meanRanks[i] - fr.meanRanks[j], z: null, p: null, pAdj: null };
      if (fr.corr !== null && fr.corr > 0) {
        row.z = row.diff / Math.sqrt(k * (k + 1) / (6 * n) * fr.corr);
        row.p = Math.min(1, 2 * d.pnorm(Math.abs(row.z), false));
        row.pAdj = Math.min(1, row.p * m);
      }
      pairs.push(row);
    }
  }
  return pairs;
}

function runFriedman(dataset, request) {
  var d = npDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < 3) throw new Error('Friedman 검정은 같은 사람을 3회 이상 측정한 변수를 3개 이상 선택하세요 (2회 측정이면 Wilcoxon 부호순위 검정을 쓰세요).');
  if (vars.length > NP_FRIEDMAN_MAX_VARS_) throw new Error('변수가 ' + vars.length + '개입니다 (최대 ' + NP_FRIEDMAN_MAX_VARS_ + '개).');
  var opt = request.options || {};
  var rowsOk = [], i, a;
  for (i = 0; i < dataset.n; i++) {
    var ok = true;
    for (a = 0; a < vars.length; a++) if (dataset.data[vars[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length, k = vars.length, warnings = [];
  if (n < 2) throw new Error('모든 변수가 유효한 사례가 2건 이상 있어야 합니다 (현재 ' + n + '건).');
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건). 한 변수라도 비어 있는 사람은 통째로 뺍니다.');
  var cols = vars.map(function (v) { return rowsOk.map(function (r) { return dataset.data[v.name][r]; }); });
  var names = vars.map(npName_);
  var fr = friedman(cols);
  if (fr.chi2 === null) warnings.push('모든 사람의 값이 전부 같아 검정할 수 없습니다.');
  else if (n < 10) warnings.push('사례가 ' + n + '건으로 적어 χ² 근사가 부정확할 수 있습니다.');

  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [{ label: '분석에 포함 (사람)', n: n }, { label: '결측으로 제외', n: dataset.n - n }, { label: '측정 횟수 (변수)', n: k }, { _kind: 'total', label: '합계', n: dataset.n }],
    notes: ['측정 순서는 변수를 보낸 순서이다. 한 사람이 모든 변수에 값이 있어야 하며, 비어 있는 사람은 통째로 제외했다.']
  });
  blocks.push({
    type: 'table', caption: '순위',
    columns: [
      { key: 'name', label: '측정(변수)' }, { key: 'n', label: 'N', format: { type: 'int' } },
      { key: 'median', label: '중앙값', format: { type: 'trim', digits: 3 } }, { key: 'meanRank', label: '평균 순위', format: { type: 'fixed', digits: 2 } }
    ],
    rows: names.map(function (nm, j) { return { name: nm, n: n, median: npMedian_(cols[j]), meanRank: fr.meanRanks[j] }; }),
    notes: ['사람마다 측정값에 순위를 매겨(동순위는 평균 순위) 측정별로 평균한 값이다. 평균 순위가 클수록 값이 큰 쪽이다.']
  });
  var testCols = [
    { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'chi2', label: 'Friedman χ²', format: NP_F3_ },
    { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'p', label: '점근 유의확률', format: { type: 'pvalue' } }
  ];
  if (opt.effect) testCols.push({ key: 'w', label: "Kendall's W", format: NP_F3_ });
  blocks.push({
    type: 'table', caption: 'Friedman 검정 통계량',
    columns: testCols,
    rows: [{ _bold: true, n: n, chi2: fr.chi2, df: fr.df, p: fr.p, w: fr.w }],
    notes: [
      '귀무가설: 모든 측정에서 값의 분포가 같다. 유의확률이 작으면(예: .05 미만) 어느 측정 이상이 다르다는 뜻이며, 어느 측정 사이인지는 사후비교로 본다.',
      'χ²는 동순위 보정을 했고 유의확률은 χ² 근사(자유도 = 측정 횟수 − 1)이다. 사례가 적으면(대략 10건 미만) 근사가 부정확하다. 0/1로 코딩된 자료에서는 이 값이 Cochran의 Q와 같다.'
    ].concat(opt.effect ? ["Kendall의 W = χ² / (N · (k − 1)) 는 0~1이며 측정 사이에 순위가 얼마나 일관되게 달라지는지를 보여 준다 (관행: 0.1 작음, 0.3 중간, 0.5 큼)."] : [])
  });
  if (opt.posthoc && fr.chi2 !== null) {
    blocks.push({
      type: 'table', caption: '사후비교: Dunn–Bonferroni (모든 쌍)',
      columns: [
        { key: 'gi', label: '(I) 측정' }, { key: 'gj', label: '(J) 측정' },
        { key: 'diff', label: '평균 순위 차이(I−J)', format: { type: 'fixed', digits: 2 } }, { key: 'sig', label: '', align: 'left' },
        { key: 'z', label: 'Z', format: NP_F3_ }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'pAdj', label: '보정 유의확률', format: { type: 'pvalue' } }
      ],
      rows: friedmanPairs(fr).map(function (pr) {
        return { gi: names[pr.i], gj: names[pr.j], diff: pr.diff, z: pr.z, p: pr.p, pAdj: pr.pAdj, sig: pr.pAdj !== null && pr.pAdj < 0.05 ? '*' : '' };
      }),
      notes: [
        '평균 순위의 차이를 표준오차 √(k(k + 1)/(6n) · C)로 나눈 Z(C: 동순위 보정 계수)로 쌍마다 비교하고, 보정 유의확률 = 유의확률 × 비교 횟수(최대 1)이다 (Dunn–Bonferroni).',
        '모든 쌍을 한 번씩(변수를 보낸 순서로 I < J) 보여 준다. * 는 보정 유의확률이 .05 미만인 쌍이다. 전체 검정이 유의하지 않으면 사후비교는 참고만 한다.'
      ]
    });
  }
  if (opt.boxplot) blocks.push(d.boxplotBlock('상자그림: 측정별 값', names.map(function (nm, j) { return { label: nm, values: cols[j] }; }), '값'));
  return { analysis: 'friedman', title: 'Friedman 검정', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = {
    npRanks_: npRanks_,
    mannWhitney: mannWhitney,
    mwExactCdf_: mwExactCdf_,
    signedRank: signedRank,
    kruskalWallis: kruskalWallis,
    dunnTest: dunnTest,
    runMannWhitney: runMannWhitney,
    runWilcoxonPaired: runWilcoxonPaired,
    runWilcoxonOne: runWilcoxonOne,
    runKruskalWallis: runKruskalWallis,
    friedman: friedman,
    friedmanPairs: friedmanPairs,
    runFriedman: runFriedman
  };
}
