/**
 * 성향점수 매칭(PSM). 순수 함수: (dataset, request) → 결과 객체.
 *
 * 1. 성향점수: 처치 ~ 공변량(숫자) + 범주형(더미) 로지스틱 회귀 (runLogistic 을 그대로 써서 모형 표와 경고를 함께 보여 준다).
 *    처치는 값이 2가지이고 기본은 큰 값이 처치(1)이다. 매칭은 logit 성향점수(선형 예측값)로 한다.
 * 2. 매칭 (탐욕적 최근접 이웃):
 *    - 처치 사례를 logit 성향점수 내림차순(동점은 데이터 순서)으로 한 번씩 처리하고, 1:k 는 k 라운드로 한다
 *      (첫 라운드에서 대조를 못 찾은 처치는 탈락, 이후 라운드에서 못 찾으면 그 처치는 다음 라운드에서 빠진다).
 *    - 거리: |logit 차이|(ps) 또는 공변량 마할라노비스 거리(maha, 공분산은 전체 사례 n−1, 더미 포함).
 *    - 캘리퍼 = 배수 × SD(logit, n−1): 안에 있는 후보만 쓴다 (0이면 없음, 마할라노비스 매칭도 캘리퍼 안에서만 고른다).
 *    - 동점은 데이터 순서가 앞선 사례. 복원(replace)이면 대조 사례를 다시 쓸 수 있다 (같은 처치에 같은 대조는 한 번만).
 * 3. 가중치: 짝지은 처치 1, 대조는 속한 세트마다 1/(세트의 대조 수)를 더한다. 매칭 ID: 짝지은 처치를 데이터 순서대로 1, 2, …
 *    (대조 사례는 속한 세트 중 가장 작은 ID).
 * 4. 균형: SMD = (처치 평균 − 대조 평균) / 처치군 SD(매칭 전, n−1). 매칭 후는 짝지은 처치의 평균과 가중 대조 평균 (MatchIt·cobalt의 ATT 방식).
 * 5. 효과(결과 변수가 있을 때): 세트별 차이 d = y처치 − mean(y대조)의 평균(ATT), SE = sd(d)/√m, t(m−1). 매칭 전은 Welch t.
 *    미관측 교란은 보정하지 못한다.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function psmDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var dd = require('./Stats_Distance.js');
    return {
      getVars: data.getVars, runLogistic: require('./Stats_Logistic.js').runLogistic, pt: dist.pt, qt: dist.qt,
      distCovariance_: dd.distCovariance_, distInverse_: dd.distInverse_, boxplotBlock: require('./Charts.js').boxplotBlock
    };
  }
  return { getVars: getVars, runLogistic: runLogistic, pt: pt, qt: qt, distCovariance_: distCovariance_, distInverse_: distInverse_, boxplotBlock: boxplotBlock };
}

var PSM_MAX_N_ = 5000;
var PSM_MAX_RATIO_ = 5;
var PSM_SMD_LIMIT_ = 0.1;
function psmName_(v) { return v.label ? v.label : v.name; }

function psmMean_(a) { var s = 0, i; for (i = 0; i < a.length; i++) s += a[i]; return a.length ? s / a.length : null; }
function psmSd_(a) {
  var m = psmMean_(a), s = 0, i;
  if (a.length < 2) return null;
  for (i = 0; i < a.length; i++) s += (a[i] - m) * (a[i] - m);
  return Math.sqrt(s / (a.length - 1));
}

/**
 * 매칭. eta: logit 성향점수, treat: 0/1, covCols: 공변량 열들(마할라노비스용), o: { caliper(SD의 배수, 0=없음), ratio, replace, method('ps'|'maha') }
 * 반환: { sets: [처치 인덱스 → 대조 인덱스 배열(없으면 빈 배열)], matched: 짝지은 처치 인덱스(오름차순), weights, ids(처치·대조 모두, 없으면 0), width(캘리퍼 폭 또는 null) }
 * 마할라노비스 매칭에서 공분산 행렬이 특이하면 { singular: true }
 */
function psmMatch(eta, treat, covCols, o) {
  var d = psmDeps_(), n = eta.length, i, c, t;
  var width = o.caliper > 0 ? o.caliper * psmSd_(eta) : null;
  var VI = null;
  if (o.method === 'maha') {
    VI = d.distInverse_(d.distCovariance_(covCols).C);
    if (!VI) return { singular: true };
  }
  var dist = function (a, b) {
    if (o.method !== 'maha') return Math.abs(eta[a] - eta[b]);
    var k = covCols.length, df = new Array(k), q = 0, r, s;
    for (r = 0; r < k; r++) df[r] = covCols[r][a] - covCols[r][b];
    for (r = 0; r < k; r++) for (s = 0; s < k; s++) q += df[r] * VI[r][s] * df[s];
    return Math.sqrt(q < 0 ? 0 : q);
  };
  var treated = [], controls = [];
  for (i = 0; i < n; i++) (treat[i] === 1 ? treated : controls).push(i);
  var order = treated.slice().sort(function (a, b) { return eta[b] - eta[a] || a - b; });
  var sets = new Array(n), used = new Array(n).fill(false);
  for (i = 0; i < n; i++) sets[i] = [];
  var active = order;
  for (var rnd = 0; rnd < o.ratio; rnd++) {
    var next = [];
    for (var a = 0; a < active.length; a++) {
      t = active[a];
      var best = -1, bd = null;
      for (var q = 0; q < controls.length; q++) {
        c = controls[q];
        if (sets[t].indexOf(c) >= 0) continue;
        if (!o.replace && used[c]) continue;
        if (width !== null && Math.abs(eta[t] - eta[c]) > width) continue;
        var dd = dist(t, c);
        if (bd === null || dd < bd) { best = c; bd = dd; }
      }
      if (best < 0) continue;
      sets[t].push(best);
      used[best] = true;
      next.push(t);
    }
    active = next;
  }
  var matched = treated.filter(function (x) { return sets[x].length > 0; });       // treated 는 오름차순
  var weights = new Array(n).fill(0), ids = new Array(n).fill(0);
  matched.forEach(function (tt, k) {
    weights[tt] += 1;
    ids[tt] = k + 1;
    sets[tt].forEach(function (cc) {
      weights[cc] += 1 / sets[tt].length;
      if (ids[cc] === 0 || k + 1 < ids[cc]) ids[cc] = k + 1;
    });
  });
  return { sets: sets, matched: matched, weights: weights, ids: ids, width: width };
}

/** 균형: cols/names, treat, weights, matched → [{ name, mt0, mc0, smd0, mt1, mc1, smd1 }] */
function psmBalance(cols, names, treat, weights, matched) {
  var ti = [], ci = [], i;
  for (i = 0; i < treat.length; i++) (treat[i] === 1 ? ti : ci).push(i);
  return cols.map(function (x, j) {
    var xt = ti.map(function (r) { return x[r]; }), xc = ci.map(function (r) { return x[r]; });
    var sd = psmSd_(xt);
    var mt0 = psmMean_(xt), mc0 = psmMean_(xc);
    var mt1 = psmMean_(matched.map(function (r) { return x[r]; }));
    var sw = 0, sx = 0;
    ci.forEach(function (r) { sw += weights[r]; sx += weights[r] * x[r]; });
    var mc1 = sw > 0 ? sx / sw : null;
    var ok = sd !== null && sd > 0;
    return { name: names[j], mt0: mt0, mc0: mc0, smd0: ok ? (mt0 - mc0) / sd : null, mt1: mt1, mc1: mc1, smd1: ok && mc1 !== null ? (mt1 - mc1) / sd : null };
  });
}

/** 세트별 차이로 ATT 와 t 검정. y: 결과(사례 순서), sets, matched → { att, se, t, df, p, lo, hi, mTreat, mCtrl } */
function psmEffect(y, sets, matched, conf) {
  var d = psmDeps_();
  var diffs = matched.map(function (t) { return y[t] - psmMean_(sets[t].map(function (c) { return y[c]; })); });
  var m = diffs.length, att = psmMean_(diffs), sd = psmSd_(diffs);
  var res = {
    att: att, se: null, t: null, df: m - 1, p: null, lo: null, hi: null, m: m,
    mTreat: psmMean_(matched.map(function (t) { return y[t]; })),
    mCtrl: psmMean_(matched.map(function (t) { return psmMean_(sets[t].map(function (c) { return y[c]; })); }))
  };
  if (m >= 2 && sd !== null) {
    res.se = sd / Math.sqrt(m);
    if (res.se > 0) {
      res.t = att / res.se;
      res.p = 2 * d.pt(Math.abs(res.t), m - 1, false);
    }
    var tc = d.qt(1 - (1 - conf / 100) / 2, m - 1);
    res.lo = att - tc * res.se; res.hi = att + tc * res.se;
  }
  return res;
}

/** 매칭 전 두 집단 평균 차이의 Welch t 검정 */
function psmWelch(a, b, conf) {
  var d = psmDeps_();
  var ma = psmMean_(a), mb = psmMean_(b), va = psmSd_(a), vb = psmSd_(b);
  var res = { ma: ma, mb: mb, diff: ma - mb, se: null, t: null, df: null, p: null, lo: null, hi: null };
  if (va === null || vb === null) return res;
  va = va * va; vb = vb * vb;
  var s2 = va / a.length + vb / b.length;
  if (!(s2 > 0)) return res;
  res.se = Math.sqrt(s2);
  res.t = res.diff / res.se;
  res.df = s2 * s2 / (va * va / (a.length * a.length * (a.length - 1)) + vb * vb / (b.length * b.length * (b.length - 1)));
  res.p = 2 * d.pt(Math.abs(res.t), res.df, false);
  var tc = d.qt(1 - (1 - conf / 100) / 2, res.df);
  res.lo = res.diff - tc * res.se; res.hi = res.diff + tc * res.se;
  return res;
}

/** 결측 사례를 뺀 데이터셋 사본 (결과 변수가 결측인 사례를 목록별로 제거할 때). map: 새 인덱스 → 원래 인덱스 */
function psmSubset_(dataset, keep) {
  var data = {};
  Object.keys(dataset.data).forEach(function (k) { data[k] = keep.map(function (i) { return dataset.data[k][i]; }); });
  return { n: keep.length, vars: dataset.vars, data: data, warnings: [], rowNumbers: keep.map(function (i) { return dataset.rowNumbers[i]; }) };
}

function runPsm(dataset, request) {
  var d = psmDeps_();
  var tv = d.getVars(dataset, request.treatment);
  if (tv.length !== 1) throw new Error('처치 변수는 하나만 선택하세요.');
  var covs = request.covariates && request.covariates.length ? d.getVars(dataset, request.covariates, { numeric: true }) : [];
  var facs = request.factors && request.factors.length ? d.getVars(dataset, request.factors) : [];
  if (covs.length + facs.length === 0) throw new Error('공변량(처치 선택에 영향을 주는 변수)을 하나 이상 선택하세요.');
  var out = request.outcome && request.outcome.length ? d.getVars(dataset, request.outcome, { numeric: true }) : [];
  if (out.length > 1) throw new Error('결과 변수는 하나만 선택하세요.');
  var seen = {};
  seen[tv[0].name] = '처치 변수';
  covs.concat(facs, out).forEach(function (v) {
    if (seen[v.name]) throw new Error('변수 "' + v.name + '"를 ' + seen[v.name] + '와(과) 다른 칸에 동시에 쓸 수 없습니다.');
    seen[v.name] = '다른 칸';
  });
  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  var caliper = opt.caliper === undefined || opt.caliper === null || opt.caliper === '' ? 0.2 : Number(opt.caliper);
  if (!isFinite(caliper) || caliper < 0) throw new Error('"캘리퍼"는 0 이상의 숫자여야 합니다 (0이면 캘리퍼를 쓰지 않습니다).');
  var ratio = opt.ratio === undefined || opt.ratio === null || opt.ratio === '' ? 1 : Number(opt.ratio);
  if (!isFinite(ratio) || ratio < 1 || ratio > PSM_MAX_RATIO_ || ratio !== Math.floor(ratio)) throw new Error('"매칭 비율(처치 1명당 대조 수)"은 1~' + PSM_MAX_RATIO_ + '의 정수여야 합니다.');
  var method = opt.method === 'maha' ? 'maha' : 'ps';
  var replace = !!opt.replace;

  // 결과 변수가 결측인 사례는 목록별로 뺀다
  var keep = null, base = dataset, warnings = [];
  if (out.length) {
    keep = [];
    for (var i = 0; i < dataset.n; i++) if (dataset.data[out[0].name][i] !== null) keep.push(i);
    if (keep.length < dataset.n) { base = psmSubset_(dataset, keep); }
  }
  var lres = d.runLogistic(base, {
    dependent: request.treatment, covariates: request.covariates || [], factors: request.factors || [],
    options: { hl: false, roc: false, flip: !!opt.flip, conf: conf, cutoff: 0.5 }, returnInternals: true
  });
  var I = lres.internals;
  delete lres.internals;
  var treat = I.y, n = treat.length, predictors = I.predictors, B = I.fit.B;
  if (n > PSM_MAX_N_) throw new Error('사례가 ' + n + '건입니다. 성향점수 매칭은 최대 ' + PSM_MAX_N_ + '건까지 계산합니다 (후보 쌍을 모두 비교하므로). 사례를 줄이거나 나눠 주세요.');
  var origOf = function (r) { return keep ? keep[I.rowsOk[r]] : I.rowsOk[r]; };    // 분석 사례 r → 원래 데이터셋 인덱스
  var eta = new Array(n), r, j;
  for (r = 0; r < n; r++) {
    var e = B[0];
    for (j = 0; j < predictors.length; j++) e += B[j + 1] * predictors[j].values[r];
    eta[r] = e;
  }
  var covCols = predictors.map(function (p) { return p.values; });
  var names = predictors.map(function (p) { return p.name; });
  var yv = out.length ? I.rowsOk.map(function (rr) { return base.data[out[0].name][rr]; }) : null;
  var nTreat = 0;
  for (r = 0; r < n; r++) nTreat += treat[r];
  var nCtrl = n - nTreat;
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건).');
  lres.warnings.forEach(function (w) { if (w.indexOf('목록별 결측 제거') !== 0) warnings.push('성향점수 모형: ' + w); });
  if (nTreat < 10) warnings.push('처치 사례가 ' + nTreat + '건으로 매우 적습니다. 매칭 결과가 불안정할 수 있습니다.');

  var M = psmMatch(eta, treat, covCols, { caliper: caliper, ratio: ratio, replace: replace, method: method });
  if (M.singular) throw new Error('공변량의 공분산 행렬을 계산할 수 없습니다 (공변량이 서로 완전히 겹침). 성향점수 최근접(ps) 방식을 쓰거나 중복된 공변량을 빼세요.');
  var m = M.matched.length;
  if (m === 0) throw new Error('매칭된 사례가 없습니다. 캘리퍼를 넓히거나(0이면 캘리퍼 없음) 매칭 방식을 바꿔 보세요.');
  if (m < nTreat) warnings.push('처치 ' + nTreat + '건 중 ' + (nTreat - m) + '건(' + ((nTreat - m) / nTreat * 100).toFixed(1) + '%)은 짝을 찾지 못해 버렸습니다. 결과는 짝지은 처치 사례에만 해당합니다 (ATT를 처치 전체로 일반화하기 어려움).');
  if (m < 10) warnings.push('짝지은 처치 사례가 ' + m + '건뿐이라 효과 추정이 불안정합니다.');

  var F3 = { type: 'fixed', digits: 3 };
  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [
      { label: '분석에 포함', n: n }, { label: '결측으로 제외', n: dataset.n - n },
      { label: '처치군', n: nTreat }, { label: '대조군', n: nCtrl }, { _kind: 'total', label: '합계', n: dataset.n }
    ],
    notes: ['목록별 결측 제거: 처치·공변량' + (out.length ? '·결과 변수' : '') + '가 모두 유효한 사례만 사용했다.']
  });
  var coding = lres.blocks.filter(function (b) { return b.type === 'table' && b.caption.indexOf('종속변수 코딩') === 0; })[0];
  if (coding) {
    blocks.push({
      type: 'table', caption: coding.caption.replace('종속변수 코딩', '처치 변수 코딩'), columns: coding.columns, rows: coding.rows,
      notes: ['내부 값 1이 처치(성향점수 모형의 "사건")이다. 기본은 큰 값이 처치이고, 옵션 "작은 값을 처치로"로 바꿀 수 있다.'].concat(coding.notes.slice(1))
    });
  }
  lres.blocks.forEach(function (b) {
    if (b.type === 'table' && (b.caption.indexOf('모형 검정') === 0 || b.caption.indexOf('모형 요약') === 0 || b.caption.indexOf('방정식의 변수') === 0)) {
      blocks.push({ type: 'table', caption: '성향점수 모형: ' + b.caption, columns: b.columns, rows: b.rows, notes: b.notes });
    }
  });

  // 공통 지지 영역
  var etaT = [], etaC = [];
  for (r = 0; r < n; r++) (treat[r] === 1 ? etaT : etaC).push(eta[r]);
  var minOf = function (a) { return Math.min.apply(null, a); }, maxOf = function (a) { return Math.max.apply(null, a); };
  var supLo = Math.max(minOf(etaT), minOf(etaC)), supHi = Math.min(maxOf(etaT), maxOf(etaC));
  var outside = etaT.filter(function (x) { return x < minOf(etaC) || x > maxOf(etaC); }).length;
  blocks.push({
    type: 'table', caption: '공통 지지 영역 (logit 성향점수)',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'min', label: '최솟값', format: F3 }, { key: 'mean', label: '평균', format: F3 }, { key: 'max', label: '최댓값', format: F3 }],
    rows: [
      { label: '처치군', n: nTreat, min: minOf(etaT), mean: psmMean_(etaT), max: maxOf(etaT) },
      { label: '대조군', n: nCtrl, min: minOf(etaC), mean: psmMean_(etaC), max: maxOf(etaC) }
    ],
    notes: [supLo <= supHi
      ? '두 집단이 겹치는 구간은 ' + supLo.toFixed(3) + ' ~ ' + supHi.toFixed(3) + '이다. 처치 사례 중 대조군 범위 밖은 ' + outside + '건이다.'
      : '두 집단의 성향점수 범위가 겹치지 않는다. 매칭할 수 없다.']
  });
  if (outside / nTreat > 0.1) warnings.push('처치 사례의 ' + (outside / nTreat * 100).toFixed(0) + '%(' + outside + '건)가 대조군의 성향점수 범위 밖입니다. 공통 지지가 부족해 매칭 결과를 일반화하기 어렵습니다.');

  var usedCtrl = {}, usedN = 0;
  M.matched.forEach(function (t) { M.sets[t].forEach(function (c) { if (!usedCtrl[c]) { usedCtrl[c] = true; usedN++; } }); });
  var methodText = (method === 'maha' ? '캘리퍼 안에서 공변량 마할라노비스 거리가 가장 가까운 대조' : '성향점수(logit)가 가장 가까운 대조') + ', 1:' + ratio + ', ' + (replace ? '복원' : '비복원') +
    (M.width !== null ? ', 캘리퍼 ' + caliper + ' × SD = ' + M.width.toFixed(4) : ', 캘리퍼 없음');
  blocks.push({
    type: 'table', caption: '매칭 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [
      { label: '처치군', n: nTreat }, { label: '짝지은 처치', n: m }, { label: '버린 처치', n: nTreat - m },
      { label: '대조군', n: nCtrl }, { label: '사용한 대조', n: usedN }, { label: '쓰지 않은 대조', n: nCtrl - usedN }
    ],
    notes: ['매칭: ' + methodText + '. 처치 사례를 성향점수가 높은 순으로 처리하고 동점은 데이터 순서가 앞선 사례를 고른다.',
      '매칭은 관찰한 공변량만 맞춘다. 관찰하지 못한 교란 변수는 보정하지 못한다.'].concat(replace ? ['복원 매칭이라 한 대조가 여러 처치와 짝지어질 수 있다. 데이터 탭에 저장되는 매칭 ID는 그중 가장 작은 번호이고, 가중치가 정확한 정보이다.'] : [])
  });

  // 균형 진단
  var bal = psmBalance(covCols.concat([eta]), names.concat(['성향점수(logit)']), treat, M.weights, M.matched);
  var bad = bal.filter(function (b) { return b.smd1 !== null && Math.abs(b.smd1) >= PSM_SMD_LIMIT_; });
  blocks.push({
    type: 'table', caption: '균형 진단 (표준화 평균 차이, SMD)',
    columns: [
      { key: 'name', label: '변수' }, { key: 'mt0', label: '처치 평균 (전)', format: F3 }, { key: 'mc0', label: '대조 평균 (전)', format: F3 }, { key: 'smd0', label: 'SMD (전)', format: F3 },
      { key: 'mt1', label: '처치 평균 (후)', format: F3 }, { key: 'mc1', label: '대조 평균 (후)', format: F3 }, { key: 'smd1', label: 'SMD (후)', format: F3 }, { key: 'verdict', label: '매칭 후 판정' }
    ],
    rows: bal.map(function (b) { return Object.assign({}, b, { verdict: b.smd1 === null ? '' : Math.abs(b.smd1) < PSM_SMD_LIMIT_ ? '균형' : '불균형', _bold: b.smd1 !== null && Math.abs(b.smd1) >= PSM_SMD_LIMIT_ }); }),
    notes: ['SMD = (처치 평균 − 대조 평균) / 처치군 표준편차(매칭 전). 매칭 후는 짝지은 처치의 평균과 가중 대조 평균이다. |SMD| < ' + PSM_SMD_LIMIT_ + '이면 균형이 잡힌 것으로 본다 (관행).']
  });
  if (bad.length) warnings.push('매칭 뒤에도 |SMD| ≥ ' + PSM_SMD_LIMIT_ + '인 변수가 있습니다: ' + bad.map(function (b) { return b.name; }).join(', ') + '. 캘리퍼를 줄이거나 공변량·매칭 방식을 바꿔 보세요.');

  // 성향점수 분포 그림 (매칭 전·후)
  var pick = function (idxs) { return idxs.map(function (x) { return eta[x]; }); };
  var tIdx = [], cIdx = [];
  for (r = 0; r < n; r++) (treat[r] === 1 ? tIdx : cIdx).push(r);
  var usedIdx = Object.keys(usedCtrl).map(Number).sort(function (a, b) { return a - b; });
  blocks.push(d.boxplotBlock('성향점수(logit) 분포: 매칭 전·후', [
    { label: '처치 (전)', values: pick(tIdx) }, { label: '대조 (전)', values: pick(cIdx) },
    { label: '처치 (후)', values: pick(M.matched) }, { label: '대조 (후)', values: pick(usedIdx) }
  ], 'logit(성향점수)'));

  // 효과 추정
  if (yv) {
    var wel = psmWelch(tIdx.map(function (x) { return yv[x]; }), cIdx.map(function (x) { return yv[x]; }), conf);
    var eff = psmEffect(yv, M.sets, M.matched, conf);
    var lab = psmName_(out[0]);
    blocks.push({
      type: 'table', caption: '처치 효과 추정 (결과 변수: ' + lab + ')',
      columns: [
        { key: 'label', label: '' }, { key: 'mt', label: '처치 평균', format: F3 }, { key: 'mc', label: '대조 평균', format: F3 }, { key: 'diff', label: '차이', format: F3 },
        { key: 'se', label: '표준오차', format: F3 }, { key: 't', label: 't', format: F3 }, { key: 'df', label: '자유도', format: { type: 'fixed', digits: 2 } },
        { key: 'p', label: '유의확률', format: { type: 'pvalue' } }, { key: 'lo', label: conf + '% 하한', format: F3 }, { key: 'hi', label: conf + '% 상한', format: F3 }
      ],
      rows: [
        { label: '매칭 전 (단순 비교, Welch)', mt: wel.ma, mc: wel.mb, diff: wel.diff, se: wel.se, t: wel.t, df: wel.df, p: wel.p, lo: wel.lo, hi: wel.hi },
        { _kind: 'group', _bold: true, label: '매칭 후 (ATT, 세트별 차이의 t 검정)', mt: eff.mTreat, mc: eff.mCtrl, diff: eff.att, se: eff.se, t: eff.t, df: eff.df, p: eff.p, lo: eff.lo, hi: eff.hi }
      ],
      notes: ['ATT(처치받은 사례에서의 평균 처치 효과) = 짝지은 처치 ' + m + '건의 (결과 − 짝 대조의 평균 결과)의 평균이다. 표준오차는 이 차이로 계산해 짝지은 자료를 반영한다.',
        '매칭 후 대조 평균은 처치 사례별 짝 대조 평균의 평균이다. 매칭 전은 두 집단 전체의 단순 비교이다.'].concat(replace ? ['복원 매칭에서는 같은 대조가 여러 세트에 쓰여 표준오차가 다소 낙관적일 수 있다.'] : [])
    });
  }

  var result = { analysis: 'psm', title: '성향점수 매칭', blocks: blocks, warnings: warnings };
  if (opt.save) {
    var prob = new Array(dataset.n).fill(null), idc = new Array(dataset.n).fill(null), wc = new Array(dataset.n).fill(null);
    for (r = 0; r < n; r++) {
      var o = origOf(r);
      prob[o] = 1 / (1 + Math.exp(-eta[r]));
      idc[o] = M.ids[r] > 0 ? M.ids[r] : null;
      wc[o] = M.weights[r];
    }
    result.saveColumns = [
      { name: '성향점수', label: '성향점수 (' + psmName_(tv[0]) + ' 처치 확률)', scale: '연속', values: prob, fixed: 4 },
      { name: '매칭ID', label: '매칭 ID (짝지은 세트 번호, 짝이 없으면 빈 값)', scale: '명목', values: idc },
      { name: '매칭가중치', label: '매칭 가중치 (짝이 없으면 0)', scale: '연속', values: wc, fixed: 3 }
    ];
  }
  return result;
}

if (typeof module !== 'undefined') {
  module.exports = { psmMatch: psmMatch, psmBalance: psmBalance, psmEffect: psmEffect, psmWelch: psmWelch, runPsm: runPsm };
}
