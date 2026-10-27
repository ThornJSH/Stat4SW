/**
 * 반복측정 분산분석 (같은 사람을 반복 측정한 평균 비교). 순수 함수: (dataset, request) → 결과 객체.
 * - 집단변수가 없으면 일원 반복측정(3회 이상 측정), 있으면 **혼합 설계**(집단 사이 요인 1개 × 반복측정 1개, 2회 이상 측정)이다.
 *   혼합 설계의 계산은 파일 아래쪽 mixedAnova 의 머리말을 본다.
 *
 * - 분산을 측정(시간), 사람(개체), 오차로 나눈다: SS측정 = n·Σ(측정 평균 − 전체 평균)², SS사람 = k·Σ(사람 평균 − 전체 평균)²,
 *   SS오차 = SS전체 − SS측정 − SS사람. F = MS측정 / MS오차, 자유도 (k − 1, (n − 1)(k − 1)).
 * - 구형성(측정 사이 차이의 분산이 같다는 가정): Mauchly 검정. S 는 측정 k개의 표본 공분산(n − 1), C 는 (k − 1) × k 직교 대비(Helmert),
 *   Σ = C·S·Cᵀ, W = det(Σ) / (tr(Σ)/(k − 1))^(k − 1). χ² = −(n − 1)·d·ln W, d = 1 − (2(k − 1)² + (k − 1) + 2)/(6(k − 1)(n − 1)),
 *   자유도 k(k − 1)/2 − 1, 유의확률은 2차 근사 (R의 mauchly.test와 pingouin이 같은 값을 냄).
 * - 교정: Greenhouse–Geisser ε = tr(Σ)²/((k − 1)·tr(Σ²)), Huynh–Feldt ε = min(1, (n(k − 1)ε_GG − 2)/((k − 1)(n − 1 − (k − 1)ε_GG))), 하한 ε = 1/(k − 1).
 *   교정한 자유도 = ε × 자유도, F는 그대로.
 * - 효과크기: 부분 η² = SS측정/(SS측정 + SS오차), 일반화 η² = SS측정/SS전체 (사람 변동까지 포함한 분모).
 * - 사후비교: 측정 쌍마다 대응표본 t(자유도 n − 1), Bonferroni 보정 유의확률 = 유의확률 × 비교 횟수, 신뢰구간은 보정한 수준(α/비교 횟수).
 * - 한 변수라도 비어 있는 사람은 통째로 제외한다 (목록별 결측 제거).
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function repDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return {
      getVars: data.getVars, compareValues_: data.compareValues_, pf: dist.pf, qt: dist.qt, pt: dist.pt, pchisq: dist.pchisq,
      leveneTest: require('./Stats_Assumptions.js').leveneTest, lineChartBlock: require('./Charts.js').lineChartBlock
    };
  }
  return { getVars: getVars, compareValues_: compareValues_, pf: pf, qt: qt, pt: pt, pchisq: pchisq, leveneTest: leveneTest, lineChartBlock: lineChartBlock };
}

var REP_MAX_VARS_ = 12;
var REP_MAX_GROUPS_ = 10;
function repName_(v) { return v.label ? v.label : v.name; }
var REP_F3_ = { type: 'fixed', digits: 3 };

/** 촐레스키로 대칭 행렬의 로그 행렬식. 양의 정부호가 아니면 null */
function repLogDet_(A) {
  var p = A.length, L = [], i, j, k, logdet = 0;
  for (i = 0; i < p; i++) {
    L.push(new Array(p).fill(0));
    for (j = 0; j <= i; j++) {
      var s = A[i][j];
      for (k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) {
        if (!(s > 1e-300)) return null;
        L[i][i] = Math.sqrt(s);
        logdet += 2 * Math.log(L[i][i]);
      } else {
        L[i][j] = s / L[j][j];
      }
    }
  }
  return logdet;
}

/**
 * 일원 반복측정 ANOVA. columns: 측정(변수)별 값 배열, 같은 길이 n (결측 없음). conf: 신뢰수준(%)
 * 반환: { n, k, means, sds, ssTime, ssSubj, ssErr, ssTotal, dfTime, dfErr, msTime, msErr, F, p, mauchly, eps, corr, pEta2, gEta2, pairs }
 * F, p 를 계산할 수 없으면(오차 제곱합 0) null. mauchly 는 n < k 이거나 공분산이 특이하면 null.
 */
function repeatedAnova(columns, conf) {
  var d = repDeps_();
  var k = columns.length, n = k ? columns[0].length : 0, i, j, t;
  var r = { n: n, k: k, means: [], sds: [], ssTime: null, ssSubj: null, ssErr: null, ssTotal: null, dfTime: k - 1, dfErr: (n - 1) * (k - 1),
    msTime: null, msErr: null, F: null, p: null, mauchly: null, eps: null, corr: null, pEta2: null, gEta2: null, pairs: [] };
  if (k < 3 || n < 2) return r;
  var grand = 0, subjMeans = new Array(n).fill(0);
  for (j = 0; j < k; j++) {
    var sm = 0;
    for (i = 0; i < n; i++) { sm += columns[j][i]; subjMeans[i] += columns[j][i] / k; }
    r.means.push(sm / n); grand += sm / (n * k);
  }
  var ssTime = 0, ssSubj = 0, ssTotal = 0;
  for (j = 0; j < k; j++) ssTime += n * (r.means[j] - grand) * (r.means[j] - grand);
  for (i = 0; i < n; i++) ssSubj += k * (subjMeans[i] - grand) * (subjMeans[i] - grand);
  for (j = 0; j < k; j++) {
    var ss = 0;
    for (i = 0; i < n; i++) { ssTotal += (columns[j][i] - grand) * (columns[j][i] - grand); ss += (columns[j][i] - r.means[j]) * (columns[j][i] - r.means[j]); }
    r.sds.push(Math.sqrt(ss / (n - 1)));
  }
  var ssErr = ssTotal - ssTime - ssSubj;
  if (ssErr < 0) ssErr = 0;
  r.ssTime = ssTime; r.ssSubj = ssSubj; r.ssErr = ssErr; r.ssTotal = ssTotal;
  r.msTime = ssTime / r.dfTime; r.msErr = ssErr / r.dfErr;
  r.pEta2 = ssTime + ssErr > 0 ? ssTime / (ssTime + ssErr) : null;
  r.gEta2 = ssTotal > 0 ? ssTime / ssTotal : null;
  if (r.msErr > 0) { r.F = r.msTime / r.msErr; r.p = d.pf(r.F, r.dfTime, r.dfErr, false); }

  // 표본 공분산 S (k × k)와 직교 대비의 공분산 Σ = C·S·Cᵀ
  var S = [];
  for (i = 0; i < k; i++) {
    S.push(new Array(k));
    for (j = 0; j <= i; j++) {
      var acc = 0;
      for (t = 0; t < n; t++) acc += (columns[i][t] - r.means[i]) * (columns[j][t] - r.means[j]);
      S[i][j] = acc / (n - 1); S[j][i] = S[i][j];
    }
  }
  var C = [], p1 = k - 1, a, b;
  for (a = 1; a <= p1; a++) {
    var row = new Array(k).fill(0);
    for (b = 0; b < a; b++) row[b] = 1 / Math.sqrt(a * (a + 1));
    row[a] = -a / Math.sqrt(a * (a + 1));
    C.push(row);
  }
  var CS = C.map(function (cr) { return S[0].map(function (_, c2) { var s2 = 0; for (var q = 0; q < k; q++) s2 += cr[q] * S[q][c2]; return s2; }); });
  var Sig = CS.map(function (row2) { return C.map(function (cr) { var s3 = 0; for (var q = 0; q < k; q++) s3 += row2[q] * cr[q]; return s3; }); });
  var tr = 0, tr2 = 0;
  for (a = 0; a < p1; a++) { tr += Sig[a][a]; for (b = 0; b < p1; b++) tr2 += Sig[a][b] * Sig[b][a]; }
  if (tr > 0 && tr2 > 0) {
    var gg = tr * tr / (p1 * tr2);
    var hf = Math.min(1, (n * p1 * gg - 2) / (p1 * (n - 1 - p1 * gg)));
    if (!isFinite(hf) || hf < gg) hf = Math.max(gg, Math.min(1, isFinite(hf) ? hf : 1));
    r.eps = { gg: gg, hf: hf, lb: 1 / p1 };
    r.corr = {};
    ['gg', 'hf', 'lb'].forEach(function (key) {
      var e = r.eps[key];
      r.corr[key] = { eps: e, df1: e * r.dfTime, df2: e * r.dfErr, p: r.F === null ? null : d.pf(r.F, e * r.dfTime, e * r.dfErr, false) };
    });
    // Mauchly: n − 1 ≥ k − 1 (사람이 측정 횟수 이상)이고 Σ 가 양의 정부호일 때만
    var logdet = n >= k ? repLogDet_(Sig) : null;
    if (logdet !== null) {
      var logW = logdet - p1 * Math.log(tr / p1), dm = 1 - (2 * p1 * p1 + p1 + 2) / (6 * p1 * (n - 1));
      var chi2 = -(n - 1) * dm * logW, ddf = k * (k - 1) / 2 - 1;
      var w2 = (p1 + 2) * (p1 - 1) * (p1 - 2) * (2 * p1 * p1 * p1 + 6 * p1 * p1 + 3 * k + 2) / (288 * Math.pow((n - 1) * p1 * dm, 2));
      var q1 = d.pchisq(chi2, ddf, false), q2 = d.pchisq(chi2, ddf + 4, false);
      r.mauchly = { W: Math.exp(logW), chi2: chi2, df: ddf, p: Math.min(1, Math.max(0, q1 + w2 * (q2 - q1))) };
    }
  }

  // 사후비교: 쌍별 대응표본 t, Bonferroni
  r.pairs = repPairwise_(columns, conf);
  return r;
}

/**
 * 측정 쌍마다 대응표본 t(자유도 n − 1), Bonferroni 보정 유의확률(× 비교 횟수, 최대 1), 보정한 수준(α/비교 횟수)의 신뢰구간.
 * columns: 측정별 값 배열 (결측 없음, 2개 이상). 반환: [{ i, j, diff, se, t, df, p, pAdj, lo, hi }]  (i < j, 차이 = i − j)
 */
function repPairwise_(columns, conf) {
  var d = repDeps_();
  var k = columns.length, n = columns[0].length, i, j, t, pairs = [];
  var m = k * (k - 1) / 2, alpha = 1 - conf / 100;
  var tcrit = n > 1 ? d.qt(1 - alpha / (2 * m), n - 1) : null;
  for (i = 0; i < k; i++) {
    for (j = i + 1; j < k; j++) {
      var dm2 = 0, ssd = 0;
      for (t = 0; t < n; t++) dm2 += (columns[i][t] - columns[j][t]) / n;
      for (t = 0; t < n; t++) { var df0 = columns[i][t] - columns[j][t] - dm2; ssd += df0 * df0; }
      var se = n > 1 ? Math.sqrt(ssd / (n - 1)) / Math.sqrt(n) : 0;
      var pr = { i: i, j: j, diff: dm2, se: se, t: null, df: n - 1, p: null, pAdj: null, lo: null, hi: null };
      if (se > 0) {
        pr.t = dm2 / se; pr.p = Math.min(1, 2 * d.pt(Math.abs(pr.t), n - 1, false)); pr.pAdj = Math.min(1, pr.p * m);
        pr.lo = dm2 - tcrit * se; pr.hi = dm2 + tcrit * se;
      }
      pairs.push(pr);
    }
  }
  return pairs;
}

function runRepeatedAnova(dataset, request) {
  if (request.factor && request.factor.length) return runMixedAnova_(dataset, request);
  var d = repDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < 3) throw new Error('반복측정 분산분석은 같은 사람을 3회 이상 측정한 변수를 3개 이상 선택하세요 (2회 측정이면 대응표본 t검정을 쓰세요).');
  if (vars.length > REP_MAX_VARS_) throw new Error('변수가 ' + vars.length + '개입니다 (최대 ' + REP_MAX_VARS_ + '개).');
  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  var rowsOk = [], i, a;
  for (i = 0; i < dataset.n; i++) {
    var ok = true;
    for (a = 0; a < vars.length; a++) if (dataset.data[vars[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length, k = vars.length, warnings = [];
  if (n < 3) throw new Error('모든 변수가 유효한 사례가 3건 이상 있어야 합니다 (현재 ' + n + '건).');
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건). 한 변수라도 비어 있는 사람은 통째로 뺍니다.');
  var cols = vars.map(function (v) { return rowsOk.map(function (r) { return dataset.data[v.name][r]; }); });
  var names = vars.map(repName_);
  var res = repeatedAnova(cols, conf);
  if (res.F === null) warnings.push('측정 사이의 차이가 사람마다 완전히 같아 오차 제곱합이 0입니다. F를 계산할 수 없습니다.');
  if (n < 10) warnings.push('사례가 ' + n + '건으로 적어 구형성 검정과 F 검정의 정확도가 낮습니다.');

  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [{ label: '분석에 포함 (사람)', n: n }, { label: '결측으로 제외', n: dataset.n - n }, { label: '측정 횟수 (변수)', n: k }, { _kind: 'total', label: '합계', n: dataset.n }],
    notes: ['측정 순서는 변수를 보낸 순서이다 (' + names.join(' → ') + '). 한 사람이 모든 변수에 값이 있어야 하며 비어 있는 사람은 통째로 제외했다.']
  });
  blocks.push({
    type: 'table', caption: '기술통계',
    columns: [{ key: 'name', label: '측정(변수)' }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'mean', label: '평균', format: REP_F3_ }, { key: 'sd', label: '표준편차', format: REP_F3_ }, { key: 'se', label: '평균의 표준오차', format: REP_F3_ }],
    rows: names.map(function (nm, j) { return { name: nm, n: n, mean: res.means[j], sd: res.sds[j], se: res.sds[j] / Math.sqrt(n) }; })
  });
  var M = res.mauchly, E = res.eps;
  blocks.push({
    type: 'table', caption: 'Mauchly 구형성 검정',
    columns: [
      { key: 'W', label: "Mauchly's W", format: REP_F3_ }, { key: 'chi2', label: '근사 카이제곱', format: REP_F3_ }, { key: 'df', label: '자유도', format: { type: 'int' } },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'gg', label: 'ε Greenhouse–Geisser', format: REP_F3_ }, { key: 'hf', label: 'ε Huynh–Feldt', format: REP_F3_ }, { key: 'lb', label: 'ε 하한', format: REP_F3_ }
    ],
    rows: [{ W: M ? M.W : null, chi2: M ? M.chi2 : null, df: M ? M.df : null, p: M ? M.p : null, gg: E ? E.gg : null, hf: E ? E.hf : null, lb: E ? E.lb : null }],
    notes: [
      '구형성은 "측정 사이 차이들의 분산이 모두 같다"는 가정이다. 유의확률이 .05 미만이면 구형성이 깨진 것이고, 그때는 아래 검정에서 ε로 자유도를 줄여 교정한 행(Greenhouse–Geisser 또는 Huynh–Feldt)을 읽는다. ε은 1에 가까울수록 구형성에 가깝고 최소값은 1/(측정 횟수 − 1)이다.',
      '유의확률은 2차 근사이며 R의 mauchly.test와 같은 값이다 (측정이 3회이면 1차 근사와 같음).' + (M ? '' : ' 사람 수가 측정 횟수보다 적거나 공분산 행렬이 특이해 이 검정은 계산하지 못했다.')
    ]
  });
  var rec = 'sph';
  if (M && M.p < 0.05) rec = E.gg > 0.75 ? 'hf' : 'gg';
  var corrRows = [['sph', '구형성 가정', res.dfTime, res.dfErr, res.p], ['gg', 'Greenhouse–Geisser', res.corr && res.corr.gg.df1, res.corr && res.corr.gg.df2, res.corr && res.corr.gg.p],
    ['hf', 'Huynh–Feldt', res.corr && res.corr.hf.df1, res.corr && res.corr.hf.df2, res.corr && res.corr.hf.p], ['lb', '하한', res.corr && res.corr.lb.df1, res.corr && res.corr.lb.df2, res.corr && res.corr.lb.p]];
  blocks.push({
    type: 'table', caption: '개체-내 효과 검정 (측정)',
    columns: [
      { key: 'name', label: '교정' }, { key: 'ssT', label: 'SS(측정)', format: REP_F3_ }, { key: 'dfT', label: '자유도(측정)', format: { type: 'trim', digits: 3 } }, { key: 'msT', label: 'MS(측정)', format: REP_F3_ },
      { key: 'ssE', label: 'SS(오차)', format: REP_F3_ }, { key: 'dfE', label: '자유도(오차)', format: { type: 'trim', digits: 3 } }, { key: 'msE', label: 'MS(오차)', format: REP_F3_ },
      { key: 'F', label: 'F', format: REP_F3_ }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
      { key: 'peta', label: '부분 η²', format: REP_F3_ }, { key: 'geta', label: '일반화 η²', format: REP_F3_ }
    ],
    rows: corrRows.map(function (cr) {
      return { _bold: cr[0] === rec, name: cr[1], ssT: res.ssTime, dfT: cr[2], msT: cr[2] ? res.ssTime / cr[2] : null, ssE: res.ssErr, dfE: cr[3], msE: cr[3] ? res.ssErr / cr[3] : null,
        F: res.F, p: cr[4], peta: res.pEta2, geta: res.gEta2 };
    }),
    notes: [
      '읽을 행을 굵게 표시했다: 구형성 검정 유의확률이 .05 이상이면 "구형성 가정", 미만이면 ε(Greenhouse–Geisser)가 .75 이하일 때 Greenhouse–Geisser, 초과일 때 Huynh–Feldt (Girden의 관행). F와 제곱합은 모든 행에서 같고 자유도와 유의확률만 교정된다.',
      '부분 η² = SS(측정) / (SS(측정) + SS(오차)), 일반화 η² = SS(측정) / SS(전체) (연구 사이 비교에 권장). SS(사람) = ' + (res.ssSubj === null ? '.' : res.ssSubj.toFixed(3)) + ', SS(전체) = ' + (res.ssTotal === null ? '.' : res.ssTotal.toFixed(3)) + '.'
    ]
  });
  if (opt.posthoc) {
    blocks.push({
      type: 'table', caption: '쌍별 비교 (Bonferroni 보정)',
      columns: [
        { key: 'gi', label: '(I) 측정' }, { key: 'gj', label: '(J) 측정' }, { key: 'diff', label: '평균차이(I−J)', format: REP_F3_ }, { key: 'sig', label: '', align: 'left' },
        { key: 'se', label: '표준오차', format: REP_F3_ }, { key: 't', label: 't', format: REP_F3_ }, { key: 'df', label: '자유도', format: { type: 'int' } },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'pAdj', label: '보정 유의확률', format: { type: 'pvalue' } },
        { key: 'lo', label: conf + '% 신뢰구간 하한(보정)', format: REP_F3_ }, { key: 'hi', label: conf + '% 신뢰구간 상한(보정)', format: REP_F3_ }
      ],
      rows: res.pairs.map(function (pr) {
        return { gi: names[pr.i], gj: names[pr.j], diff: pr.diff, se: pr.se, t: pr.t, df: pr.df, p: pr.p, pAdj: pr.pAdj, lo: pr.lo, hi: pr.hi, sig: pr.pAdj !== null && pr.pAdj < 0.05 ? '*' : '' };
      }),
      notes: ['측정 쌍마다 대응표본 t검정을 하고 보정 유의확률 = 유의확률 × 비교 횟수(최대 1), 신뢰구간은 신뢰수준을 비교 횟수로 나눠 보정한 구간이다. * 는 보정 유의확률이 .05 미만인 쌍이다. 전체 검정이 유의하지 않으면 참고만 한다.',
        '모든 쌍을 한 번씩(변수를 보낸 순서로 I < J) 보여 준다. 구형성이 깨졌어도 쌍별 비교는 쌍마다 차이의 분산을 따로 쓰므로 그대로 쓸 수 있다.']
    });
  }
  if (opt.plot) {
    blocks.push(d.lineChartBlock('측정별 평균', '측정 순서 (' + names.join(' → ') + ')', '평균', [{ label: '평균', points: res.means.map(function (mv, j) { return [j + 1, mv]; }) }]));
  }
  if (M && M.p < 0.05) warnings.push('구형성이 깨졌습니다 (Mauchly p = ' + (M.p < 0.001 ? '<.001' : M.p.toFixed(3)) + ', ε(G-G) = ' + E.gg.toFixed(3) + '). 굵게 표시한 교정 행을 읽으세요.');
  if (res.F !== null && n >= 3) warnings.push('정규성이 의심되면 순위로 비교하는 Friedman 검정을 함께 보세요.');
  return { analysis: 'rmanova', title: '일원 반복측정 분산분석', blocks: blocks, warnings: warnings };
}

/* ---------- 혼합 설계 (집단 × 반복측정) ---------- */

/**
 * 혼합 설계 분산분석 (집단 g개 × 측정 k개, 집단 크기가 달라도 되며 제곱합은 유형 III).
 * columns: 측정별 값 배열(결측 없음), gidx: 사람별 집단 번호(0부터), g: 집단 수. 반환 구조는 아래 주석의 항목들.
 * - 개체-간: 사람별 평균에 √k 를 곱한 값 w 의 일원 ANOVA. SS집단 = Σ n_h (w̄_h − w̄)², SS오차 = Σ (w − w̄_h)², 절편 SS = g² (비가중 평균)² / Σ(1/n_h)
 * - 개체-내: 측정 k개를 직교(Helmert) 대비 k − 1개로 바꾼 점수 Z. 시간 = Z 의 절편(집단 평균의 비가중 평균 = 0 검정) 제곱합의 합,
 *   시간×집단 = Z 의 집단 효과 제곱합의 합, 오차 = 집단 안 잔차 제곱합의 합. 자유도 k − 1, (g − 1)(k − 1), (N − g)(k − 1).
 * - 구형성: 집단 안에서 합친 Z 의 공분산 Σ = (잔차'잔차)/(N − g)에 Mauchly 검정(자유도 N − g, 2차 근사), ε Greenhouse–Geisser = tr(Σ)²/((k − 1) tr(Σ²)),
 *   Huynh–Feldt = ((N − g + 1)(k − 1) ε − 2)/((k − 1)(N − g − (k − 1) ε)) (Lecoutre 보정; R car), 하한 1/(k − 1). 교정한 자유도 = ε × 자유도.
 */
function mixedAnova(columns, gidx, g) {
  var d = repDeps_();
  var k = columns.length, N = columns[0].length, p = k - 1, i, j, h, c, a;
  var nh = new Array(g).fill(0);
  for (i = 0; i < N; i++) nh[gidx[i]]++;
  var inv = 0;
  nh.forEach(function (v) { inv += 1 / v; });
  var r = { N: N, k: k, g: g, nh: nh };
  // 사람별 평균(√k 배)
  var w = new Array(N).fill(0), wh = new Array(g).fill(0), wAll = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < k; j++) w[i] += columns[j][i];
    w[i] /= Math.sqrt(k);
    wh[gidx[i]] += w[i];
    wAll += w[i];
  }
  wAll /= N;
  for (h = 0; h < g; h++) wh[h] /= nh[h];
  var ssGroup = 0, ssErrB = 0, muW = 0;
  for (h = 0; h < g; h++) { ssGroup += nh[h] * (wh[h] - wAll) * (wh[h] - wAll); muW += wh[h] / g; }
  for (i = 0; i < N; i++) ssErrB += (w[i] - wh[gidx[i]]) * (w[i] - wh[gidx[i]]);
  var dfEb = N - g, msErrB = ssErrB / dfEb;
  var ssIntB = g * g * muW * muW / inv;
  var Fb = function (ss, df) { return msErrB > 0 ? (ss / df) / msErrB : null; };
  r.between = {
    intercept: { ss: ssIntB, df: 1, ms: ssIntB, F: Fb(ssIntB, 1) },
    group: { ss: ssGroup, df: g - 1, ms: ssGroup / (g - 1), F: Fb(ssGroup, g - 1), peta: ssGroup + ssErrB > 0 ? ssGroup / (ssGroup + ssErrB) : null },
    err: { ss: ssErrB, df: dfEb, ms: msErrB }
  };
  r.between.intercept.p = r.between.intercept.F === null ? null : d.pf(r.between.intercept.F, 1, dfEb, false);
  r.between.group.p = r.between.group.F === null ? null : d.pf(r.between.group.F, g - 1, dfEb, false);

  // 직교 대비 점수 Z (N × p)
  var Z = [];
  for (i = 0; i < N; i++) {
    var row = [];
    for (a = 1; a <= p; a++) {
      var s = 0;
      for (j = 0; j < a; j++) s += columns[j][i] / Math.sqrt(a * (a + 1));
      s -= a * columns[a][i] / Math.sqrt(a * (a + 1));
      row.push(s);
    }
    Z.push(row);
  }
  var zbar = [];
  for (h = 0; h < g; h++) zbar.push(new Array(p).fill(0));
  var zw = new Array(p).fill(0);
  for (i = 0; i < N; i++) for (c = 0; c < p; c++) { zbar[gidx[i]][c] += Z[i][c] / nh[gidx[i]]; zw[c] += Z[i][c] / N; }
  var ssTime = 0, ssTG = 0, ssErrW = 0;
  for (c = 0; c < p; c++) {
    var mu = 0;
    for (h = 0; h < g; h++) { mu += zbar[h][c] / g; ssTG += nh[h] * (zbar[h][c] - zw[c]) * (zbar[h][c] - zw[c]); }
    ssTime += g * g * mu * mu / inv;
  }
  var res = Z.map(function (z, s2) { return z.map(function (v, c2) { return v - zbar[gidx[s2]][c2]; }); });
  res.forEach(function (rw) { rw.forEach(function (v) { ssErrW += v * v; }); });
  var dfEw = dfEb * p;
  var msErrW = ssErrW / dfEw;
  r.errw = { ss: ssErrW, df: dfEw, ms: msErrW };
  var effect = function (ss, df) {
    var F = msErrW > 0 ? (ss / df) / msErrW : null;
    return { ss: ss, df: df, ms: ss / df, F: F, p: F === null ? null : d.pf(F, df, dfEw, false), peta: ss + ssErrW > 0 ? ss / (ss + ssErrW) : null, corr: null };
  };
  r.time = effect(ssTime, p);
  r.inter = effect(ssTG, (g - 1) * p);

  // 집단 안에서 합친 공분산, 구형성, ε
  r.eps = null; r.mauchly = null;
  if (p >= 2 && dfEb >= 1) {
    var Sig = [];
    for (a = 0; a < p; a++) {
      Sig.push(new Array(p));
      for (var b2 = 0; b2 <= a; b2++) {
        var acc = 0;
        for (i = 0; i < N; i++) acc += res[i][a] * res[i][b2];
        Sig[a][b2] = acc / dfEb; Sig[b2][a] = Sig[a][b2];
      }
    }
    var tr = 0, tr2 = 0, b3;
    for (a = 0; a < p; a++) { tr += Sig[a][a]; for (b3 = 0; b3 < p; b3++) tr2 += Sig[a][b3] * Sig[b3][a]; }
    if (tr > 0 && tr2 > 0) {
      var gg = tr * tr / (p * tr2);
      var hf = ((dfEb + 1) * p * gg - 2) / (p * (dfEb - p * gg));
      hf = isFinite(hf) ? Math.min(1, hf) : 1;
      hf = Math.max(hf, gg);
      r.eps = { gg: gg, hf: hf, lb: 1 / p };
      [r.time, r.inter].forEach(function (ef) {
        ef.corr = {};
        ['gg', 'hf', 'lb'].forEach(function (key) {
          var e = r.eps[key];
          ef.corr[key] = { eps: e, df1: e * ef.df, df2: e * dfEw, p: ef.F === null ? null : d.pf(ef.F, e * ef.df, e * dfEw, false) };
        });
      });
      var logdet = dfEb >= p ? repLogDet_(Sig) : null;
      if (logdet !== null) {
        var logW = logdet - p * Math.log(tr / p), dm = 1 - (2 * p * p + p + 2) / (6 * p * dfEb);
        var chi2 = -dfEb * dm * logW, ddf = k * (k - 1) / 2 - 1;
        var w2 = (p + 2) * (p - 1) * (p - 2) * (2 * p * p * p + 6 * p * p + 3 * k + 2) / (288 * Math.pow(dfEb * p * dm, 2));
        var q1 = d.pchisq(chi2, ddf, false), q2 = d.pchisq(chi2, ddf + 4, false);
        r.mauchly = { W: Math.exp(logW), chi2: chi2, df: ddf, p: Math.min(1, Math.max(0, q1 + w2 * (q2 - q1))) };
      }
    }
  }

  // 집단 × 측정 평균·SD·N, 측정별 전체 평균·SD
  r.cells = [];
  for (h = 0; h < g; h++) {
    r.cells.push([]);
    for (j = 0; j < k; j++) {
      var vals = [];
      for (i = 0; i < N; i++) if (gidx[i] === h) vals.push(columns[j][i]);
      var m = 0, ss = 0;
      vals.forEach(function (v) { m += v; });
      m /= vals.length;
      vals.forEach(function (v) { ss += (v - m) * (v - m); });
      r.cells[h].push({ n: vals.length, mean: m, sd: vals.length > 1 ? Math.sqrt(ss / (vals.length - 1)) : null });
    }
  }
  return r;
}

function repLevels_(v, vals) {
  var d = repDeps_();
  var cats = [], seen = {};
  vals.forEach(function (x) {
    var key = typeof x + ':' + x;
    if (!seen[key]) { seen[key] = true; cats.push({ value: x }); }
  });
  cats.sort(d.compareValues_);
  var pos = {};
  cats.forEach(function (c, i) {
    pos[typeof c.value + ':' + c.value] = i;
    var label = v.valueLabels[String(c.value)];
    c.label = label !== undefined ? label : String(c.value);
  });
  return { levels: cats, index: vals.map(function (x) { return pos[typeof x + ':' + x]; }) };
}

function repP_(p) { return p === null ? '.' : (p < 0.001 ? '<.001' : p.toFixed(3).replace(/^0/, '')); }

function runMixedAnova_(dataset, request) {
  var d = repDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < 2) throw new Error('집단변수를 넣은 반복측정(혼합 설계)은 같은 사람을 2회 이상 측정한 변수를 2개 이상 선택하세요.');
  if (vars.length > REP_MAX_VARS_) throw new Error('변수가 ' + vars.length + '개입니다 (최대 ' + REP_MAX_VARS_ + '개).');
  var fv = d.getVars(dataset, request.factor);
  if (fv.length !== 1) throw new Error('집단변수는 하나만 선택하세요.');
  if (vars.some(function (v) { return v.name === fv[0].name; })) throw new Error('변수 "' + fv[0].name + '"를 측정 변수와 집단변수에 동시에 쓸 수 없습니다.');
  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  var rowsOk = [], i, a;
  for (i = 0; i < dataset.n; i++) {
    var ok = dataset.data[fv[0].name][i] !== null;
    for (a = 0; ok && a < vars.length; a++) if (dataset.data[vars[a].name][i] === null) ok = false;
    if (ok) rowsOk.push(i);
  }
  var N = rowsOk.length, k = vars.length, warnings = [];
  var lv = repLevels_(fv[0], rowsOk.map(function (r) { return dataset.data[fv[0].name][r]; }));
  var g = lv.levels.length, groupName = repName_(fv[0]);
  if (g < 2) throw new Error('집단변수 "' + groupName + '"의 집단이 ' + g + '개입니다. 2개 이상이어야 합니다.');
  if (g > REP_MAX_GROUPS_) throw new Error('집단변수 "' + groupName + '"의 집단이 ' + g + '개입니다 (최대 ' + REP_MAX_GROUPS_ + '개). 범주가 적은 변수를 선택하세요.');
  if (N <= g) throw new Error('유효한 사례(' + N + '건)가 집단 수(' + g + ')보다 많아야 합니다.');
  var cols = vars.map(function (v) { return rowsOk.map(function (r) { return dataset.data[v.name][r]; }); });
  var names = vars.map(repName_);
  var labelOf = function (h) { return lv.levels[h].label; };
  var res = mixedAnova(cols, lv.index, g);
  if (res.nh.some(function (n) { return n < 1; })) throw new Error('사례가 없는 집단이 있습니다.');
  if (dataset.n > N) warnings.push('목록별 결측 제거로 ' + (dataset.n - N) + '건을 제외했습니다 (사용한 사례 ' + N + '건). 측정 변수나 집단변수 중 하나라도 비어 있는 사람은 통째로 뺍니다.');
  var balanced = res.nh.every(function (n) { return n === res.nh[0]; });
  if (!balanced) warnings.push('집단마다 사람 수가 다릅니다 (' + res.nh.join(', ') + '명). 제곱합은 유형 III(비가중 평균 기준)으로 계산했습니다.');
  if (res.nh.some(function (n) { return n < 5; })) warnings.push('사람이 5명 미만인 집단이 있습니다. 구형성 검정과 F 검정의 정확도가 낮습니다.');
  if (!(res.errw.ms > 0)) warnings.push('개체-내 오차 제곱합이 0입니다. 측정 변화가 집단 안에서 사람마다 완전히 같아 F를 계산할 수 없습니다.');

  var F3 = REP_F3_;
  var blocks = [];
  var sumRows = [{ label: '분석에 포함 (사람)', n: N }, { label: '결측으로 제외', n: dataset.n - N }, { label: '측정 횟수 (변수)', n: k }, { label: '집단 수', n: g }];
  res.nh.forEach(function (n, h) { sumRows.push({ label: '집단: ' + labelOf(h), n: n }); });
  sumRows.push({ _kind: 'total', label: '합계', n: dataset.n });
  blocks.push({
    type: 'table', caption: '처리 요약', columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }], rows: sumRows,
    notes: ['혼합 설계: 집단변수 "' + groupName + '"(집단 사이) × 측정(같은 사람 안). 측정 순서는 변수를 보낸 순서이다 (' + names.join(' → ') + ').']
  });

  var descRows = [];
  for (var h = 0; h < g; h++) {
    for (var j = 0; j < k; j++) {
      var cell = res.cells[h][j];
      descRows.push({ grp: labelOf(h), name: names[j], n: cell.n, mean: cell.mean, sd: cell.sd, se: cell.sd === null ? null : cell.sd / Math.sqrt(cell.n) });
    }
  }
  for (j = 0; j < k; j++) {
    var all = cols[j], m = 0, ss = 0;
    all.forEach(function (v) { m += v; });
    m /= N;
    all.forEach(function (v) { ss += (v - m) * (v - m); });
    var sd = Math.sqrt(ss / (N - 1));
    descRows.push({ grp: '전체', name: names[j], n: N, mean: m, sd: sd, se: sd / Math.sqrt(N) });
  }
  blocks.push({
    type: 'table', caption: '기술통계 (집단 × 측정)',
    columns: [{ key: 'grp', label: groupName }, { key: 'name', label: '측정(변수)' }, { key: 'n', label: 'N', format: { type: 'int' } },
      { key: 'mean', label: '평균', format: F3 }, { key: 'sd', label: '표준편차', format: F3 }, { key: 'se', label: '평균의 표준오차', format: F3 }],
    rows: descRows
  });

  if (opt.levene !== false) {
    var lrows = [];
    for (j = 0; j < k; j++) {
      var grp = [];
      for (h = 0; h < g; h++) grp.push(cols[j].filter(function (_, s) { return lv.index[s] === h; }));
      var lt = d.leveneTest(grp, false);
      lrows.push({ name: names[j], W: lt.W, df1: lt.df1, df2: lt.df2, p: lt.p });
      if (lt.p !== null && lt.p < 0.05) warnings.push('"' + names[j] + '"에서 Levene 등분산 검정이 유의합니다 (p = ' + repP_(lt.p) + '). 집단 사이 분산이 같다는 가정이 의심됩니다.');
    }
    blocks.push({
      type: 'table', caption: 'Levene 등분산 검정 (측정별)',
      columns: [{ key: 'name', label: '측정(변수)' }, { key: 'W', label: 'Levene 통계량', format: F3 }, { key: 'df1', label: '자유도1', format: { type: 'int' } },
        { key: 'df2', label: '자유도2', format: { type: 'int' } }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }],
      rows: lrows,
      notes: ['각 측정에서 집단 사이 분산이 같다는 귀무가설을 평균 기준으로 검정한다. 유의하면(p < .05) 집단 사이 비교(개체-간 효과)에 주의한다.']
    });
  }

  var M = res.mauchly, E = res.eps;
  if (k >= 3) {
    blocks.push({
      type: 'table', caption: 'Mauchly 구형성 검정',
      columns: [
        { key: 'W', label: "Mauchly's W", format: F3 }, { key: 'chi2', label: '근사 카이제곱', format: F3 }, { key: 'df', label: '자유도', format: { type: 'int' } },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'gg', label: 'ε Greenhouse–Geisser', format: F3 }, { key: 'hf', label: 'ε Huynh–Feldt', format: F3 }, { key: 'lb', label: 'ε 하한', format: F3 }
      ],
      rows: [{ W: M ? M.W : null, chi2: M ? M.chi2 : null, df: M ? M.df : null, p: M ? M.p : null, gg: E ? E.gg : null, hf: E ? E.hf : null, lb: E ? E.lb : null }],
      notes: [
        '구형성은 "측정 사이 차이들의 분산이 모두 같다"는 가정이다. 집단 안에서 합친 공분산(자유도 ' + (N - g) + ')으로 검정한다. 유의확률이 .05 미만이면 아래 두 효과 검정에서 ε로 자유도를 줄여 교정한 행을 읽는다.',
        '유의확률은 2차 근사이며 R의 mauchly.test와 같은 방식이다.' + (M ? '' : ' 사람 수가 부족하거나 공분산 행렬이 특이해 이 검정은 계산하지 못했다.')
      ]
    });
  } else {
    blocks.push({ type: 'text', text: '▶ 측정이 2회라 구형성 검정과 ε 교정이 필요 없습니다 (아래 검정은 구형성 가정 행만 있습니다).' });
  }

  var rec = 'sph';
  if (M && M.p < 0.05) rec = E.gg > 0.75 ? 'hf' : 'gg';
  var effTable = function (caption, ef, label) {
    var rowsDef = [['sph', '구형성 가정', ef.df, res.errw.df, ef.p]];
    if (k >= 3) {
      ['gg', 'hf', 'lb'].forEach(function (key) {
        var nm = key === 'gg' ? 'Greenhouse–Geisser' : key === 'hf' ? 'Huynh–Feldt' : '하한';
        rowsDef.push([key, nm, ef.corr && ef.corr[key].df1, ef.corr && ef.corr[key].df2, ef.corr && ef.corr[key].p]);
      });
    }
    return {
      type: 'table', caption: caption,
      columns: [
        { key: 'name', label: '교정' }, { key: 'ssT', label: 'SS(' + label + ')', format: F3 }, { key: 'dfT', label: '자유도(' + label + ')', format: { type: 'trim', digits: 3 } }, { key: 'msT', label: 'MS(' + label + ')', format: F3 },
        { key: 'ssE', label: 'SS(오차)', format: F3 }, { key: 'dfE', label: '자유도(오차)', format: { type: 'trim', digits: 3 } }, { key: 'msE', label: 'MS(오차)', format: F3 },
        { key: 'F', label: 'F', format: F3 }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'peta', label: '부분 η²', format: F3 }
      ],
      rows: rowsDef.map(function (cr) {
        return { _bold: cr[0] === rec, name: cr[1], ssT: ef.ss, dfT: cr[2], msT: cr[2] ? ef.ss / cr[2] : null, ssE: res.errw.ss, dfE: cr[3], msE: cr[3] ? res.errw.ss / cr[3] : null, F: ef.F, p: cr[4], peta: ef.peta };
      }),
      notes: k >= 3 ? ['읽을 행을 굵게 표시했다: 구형성 검정 유의확률이 .05 이상이면 "구형성 가정", 미만이면 ε(Greenhouse–Geisser)가 .75 이하일 때 Greenhouse–Geisser, 초과일 때 Huynh–Feldt. F와 제곱합은 모든 행에서 같고 자유도와 유의확률만 교정된다. Huynh–Feldt는 Lecoutre 보정식이다 (R car 와 같음).'] : []
    };
  };
  blocks.push(effTable('개체-내 효과 검정 (측정)', res.time, '측정'));
  blocks.push(effTable('개체-내 효과 검정 (측정 × ' + groupName + ')', res.inter, '측정×집단'));
  var intP = res.inter.corr && rec !== 'sph' ? res.inter.corr[rec].p : res.inter.p;
  var bNotes = ['제곱합은 유형 III이다. 개체-간 효과는 사람별 평균(측정 전체의 평균)을 집단 사이에서 비교한 것이다. 부분 η² = 제곱합 /(제곱합 + 오차 제곱합).'];
  if (intP !== null && intP < 0.05) bNotes.push('측정 × 집단 상호작용이 유의합니다 (p = ' + repP_(intP) + '). 집단에 따라 측정 변화가 다르다는 뜻이니 집단 주효과나 측정 주효과를 그대로 해석하지 말고 상호작용 그림과 단순 효과를 보세요.');
  else if (intP !== null) bNotes.push('측정 × 집단 상호작용이 유의하지 않습니다. 측정 주효과와 집단 주효과를 각각 해석할 수 있습니다.');
  blocks.push({
    type: 'table', caption: '개체-간 효과 검정 (집단)',
    columns: [{ key: 'source', label: '' }, { key: 'ss', label: '제곱합', format: F3 }, { key: 'df', label: '자유도', format: { type: 'int' } }, { key: 'ms', label: '평균제곱', format: F3 },
      { key: 'F', label: 'F', format: F3 }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'peta', label: '부분 η²', format: F3 }],
    rows: [
      { source: '절편', ss: res.between.intercept.ss, df: 1, ms: res.between.intercept.ms, F: res.between.intercept.F, p: res.between.intercept.p },
      { source: groupName, ss: res.between.group.ss, df: res.between.group.df, ms: res.between.group.ms, F: res.between.group.F, p: res.between.group.p, peta: res.between.group.peta },
      { source: '오차', ss: res.between.err.ss, df: res.between.err.df, ms: res.between.err.ms }
    ],
    notes: bNotes
  });

  // 단순 효과: 집단별 측정 효과(집단 안 반복측정 ANOVA), 측정별 집단 효과(일원 ANOVA)
  if (opt.simple !== false) {
    var sRows = [];
    for (h = 0; h < g; h++) {
      var sub = cols.map(function (c) { return c.filter(function (_, s) { return lv.index[s] === h; }); });
      var row = { grp: labelOf(h), n: res.nh[h], df1: k - 1, df2: (res.nh[h] - 1) * (k - 1), F: null, p: null, pgg: null };
      if (res.nh[h] >= 2) {
        if (k >= 3) {
          var ra = repeatedAnova(sub, conf);
          row.F = ra.F; row.p = ra.p; row.pgg = ra.corr ? ra.corr.gg.p : null;
        } else {
          var pr = repPairwise_(sub, conf)[0];
          row.F = pr.t === null ? null : pr.t * pr.t; row.p = pr.p;
        }
      }
      sRows.push(row);
    }
    blocks.push({
      type: 'table', caption: '측정의 단순 효과 (집단별)',
      columns: [{ key: 'grp', label: groupName }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'df1', label: '자유도1', format: { type: 'int' } }, { key: 'df2', label: '자유도2', format: { type: 'int' } },
        { key: 'F', label: 'F', format: F3 }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }].concat(k >= 3 ? [{ key: 'pgg', label: '유의확률(Greenhouse–Geisser)', format: { type: 'pvalue' } }] : []),
      rows: sRows,
      notes: ['각 집단 안에서 측정 사이 평균이 달라졌는지 따로 반복측정 분산분석한 결과이다 (그 집단의 오차항만 사용' + (k === 2 ? '; 측정이 2회라 대응표본 t검정의 t²과 같음' : '; 구형성이 의심되면 Greenhouse–Geisser 유의확률을 본다') + '). 집단 수만큼 검정하므로 유의수준을 나눠서(예: .05/집단 수) 보는 것이 안전하다. 상호작용이 유의하지 않으면 참고만 한다.']
    });
    var gRows = [];
    for (j = 0; j < k; j++) {
      var gm = [], gn = [], tot = 0, sw = 0, within = 0;
      for (h = 0; h < g; h++) { gm.push(res.cells[h][j].mean); gn.push(res.cells[h][j].n); tot += gn[h]; sw += gn[h] * gm[h]; }
      var gmean = sw / tot, between = 0;
      for (h = 0; h < g; h++) between += gn[h] * (gm[h] - gmean) * (gm[h] - gmean);
      for (i = 0; i < N; i++) within += (cols[j][i] - gm[lv.index[i]]) * (cols[j][i] - gm[lv.index[i]]);
      var Fg = within > 0 ? (between / (g - 1)) / (within / (N - g)) : null;
      gRows.push({ name: names[j], df1: g - 1, df2: N - g, F: Fg, p: Fg === null ? null : d.pf(Fg, g - 1, N - g, false) });
    }
    blocks.push({
      type: 'table', caption: '집단의 단순 효과 (측정별)',
      columns: [{ key: 'name', label: '측정(변수)' }, { key: 'df1', label: '자유도1', format: { type: 'int' } }, { key: 'df2', label: '자유도2', format: { type: 'int' } },
        { key: 'F', label: 'F', format: F3 }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } }],
      rows: gRows,
      notes: ['각 측정에서 집단 사이 평균이 다른지 일원 분산분석한 결과이다 (등분산 가정; 집단이 2개면 독립표본 t검정의 t²과 같음). 측정 수만큼 검정하므로 유의수준을 나눠서 보는 것이 안전하다.']
    });
  }

  // 쌍별 비교: 집단 안에서 측정 쌍마다
  if (opt.posthoc) {
    var pRows = [];
    for (h = 0; h < g; h++) {
      var sub2 = cols.map(function (c) { return c.filter(function (_, s) { return lv.index[s] === h; }); });
      if (res.nh[h] < 2) continue;
      repPairwise_(sub2, conf).forEach(function (pr2) {
        pRows.push({ grp: labelOf(h), gi: names[pr2.i], gj: names[pr2.j], diff: pr2.diff, se: pr2.se, t: pr2.t, df: pr2.df, p: pr2.p, pAdj: pr2.pAdj, lo: pr2.lo, hi: pr2.hi, sig: pr2.pAdj !== null && pr2.pAdj < 0.05 ? '*' : '' });
      });
    }
    blocks.push({
      type: 'table', caption: '쌍별 비교 (집단별, Bonferroni 보정)',
      columns: [{ key: 'grp', label: groupName }, { key: 'gi', label: '(I) 측정' }, { key: 'gj', label: '(J) 측정' }, { key: 'diff', label: '평균차이(I−J)', format: F3 }, { key: 'sig', label: '', align: 'left' },
        { key: 'se', label: '표준오차', format: F3 }, { key: 't', label: 't', format: F3 }, { key: 'df', label: '자유도', format: { type: 'int' } },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'pAdj', label: '보정 유의확률', format: { type: 'pvalue' } },
        { key: 'lo', label: conf + '% 신뢰구간 하한(보정)', format: F3 }, { key: 'hi', label: conf + '% 신뢰구간 상한(보정)', format: F3 }],
      rows: pRows,
      notes: ['각 집단 안에서 측정 쌍마다 대응표본 t검정을 하고 보정 유의확률 = 유의확률 × 한 집단의 비교 횟수(최대 1)이다. * 는 보정 유의확률이 .05 미만인 쌍이다. 집단 수만큼 반복되므로 필요하면 더 엄격하게 본다.']
    });
  }

  if (opt.plot) {
    var series = [];
    for (h = 0; h < g; h++) series.push({ label: groupName + ' = ' + labelOf(h), points: res.cells[h].map(function (c, jj) { return [jj + 1, c.mean]; }) });
    blocks.push(d.lineChartBlock('상호작용 그림: 측정별 집단 평균', '측정 순서', '평균', series, names));
  }
  if (M && M.p < 0.05) warnings.push('구형성이 깨졌습니다 (Mauchly p = ' + repP_(M.p) + ', ε(G-G) = ' + E.gg.toFixed(3) + '). 굵게 표시한 교정 행을 읽으세요.');
  return { analysis: 'rmanova', title: '혼합 설계 분산분석 (집단 × 반복측정)', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = { repeatedAnova: repeatedAnova, mixedAnova: mixedAnova, repPairwise_: repPairwise_, runRepeatedAnova: runRepeatedAnova };
}
