/**
 * 평균 비교: 일표본, 독립표본, 대응표본 t검정. 순수 함수: (dataset, request) → 결과 객체.
 * 평균차이의 방향: 집단1 − 집단2, 변수1 − 변수2.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function cmpDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var desc = require('./Stats_Descriptive.js');
    var asm = require('./Stats_Assumptions.js');
    var rel = require('./Stats_Relation.js');
    var eff = require('./Stats_Effect.js');
    return {
      effectIndependent: eff.effectIndependent, effectPaired: eff.effectPaired, effectOneSample: eff.effectOneSample, effectSizeWord: eff.effectSizeWord,
      getVars: data.getVars, splitGroups: data.splitGroups,
      qt: dist.qt, pt: dist.pt, pValueT: dist.pValueT,
      centralSums: desc.centralSums, leveneTest: asm.leveneTest, pearsonR: rel.pearsonR,
      boxplotBlock: require('./Charts.js').boxplotBlock,
      advIndependentT_: asm.advIndependentT_, advPaired_: asm.advPaired_, advOneSample_: asm.advOneSample_, advAnova_: asm.advAnova_,
      advCorrelation_: asm.advCorrelation_, advAdd_: asm.advAdd_, advAppend_: asm.advAppend_
    };
  }
  return {
    effectIndependent: effectIndependent, effectPaired: effectPaired, effectOneSample: effectOneSample, effectSizeWord: effectSizeWord,
    getVars: getVars, splitGroups: splitGroups, qt: qt, pt: pt, pValueT: pValueT,
    centralSums: centralSums, leveneTest: leveneTest, pearsonR: pearsonR,
    boxplotBlock: boxplotBlock,
    advIndependentT_: advIndependentT_, advPaired_: advPaired_, advOneSample_: advOneSample_, advAnova_: advAnova_,
    advCorrelation_: advCorrelation_, advAdd_: advAdd_, advAppend_: advAppend_
  };
}

function nameOf_(v) { return v.label ? v.label : v.name; }

/** 숫자 옵션을 읽는다. 값이 없으면 기본값, 숫자가 아니거나 범위 밖이면 사용자용 오류 메시지. */
function numberOption_(options, key, label, def, min, max) {
  var raw = options ? options[key] : undefined;
  if (raw === undefined || raw === null || raw === '') return def;
  var v = Number(raw);
  if (!isFinite(v)) throw new Error('"' + label + '"에는 숫자를 입력하세요.');
  if ((min !== undefined && v <= min) || (max !== undefined && v >= max)) {
    throw new Error('"' + label + '"은(는) ' + (min !== undefined ? min + ' 초과 ' : '') + (max !== undefined ? max + ' 미만' : '') + '이어야 합니다.');
  }
  return v;
}

function confidenceOption_(options) {
  return numberOption_(options, 'conf', '신뢰수준(%)', 95, 0, 100);
}

function ciLabel_(conf, side) {
  return '차이의 ' + conf + '% 신뢰구간 ' + side;
}

/* ---------- 계산 ---------- */

/** t분포 위꼬리 임계값: P(T > tc) = (1 - 신뢰수준)/2 */
function tCritical_(conf, df) {
  return cmpDeps_().qt((1 - conf / 100) / 2, df, false);
}

/** 한 집단의 평균, 표준편차, 표준오차 (계산할 수 없는 값은 null) */
function groupStats_(values) {
  var n = values.length;
  var r = { n: n, mean: null, variance: null, sd: null, se: null };
  if (n === 0) return r;
  var cs = cmpDeps_().centralSums(values);
  r.mean = cs.mean;
  if (n > 1) {
    r.variance = cs.c2 / (n - 1);
    r.sd = Math.sqrt(r.variance);
    r.se = r.sd / Math.sqrt(n);
  }
  return r;
}

/** 일표본 t검정. diff = 평균 − 검정값. 표준오차가 0이면 t와 p는 null. */
function oneSampleT(values, mu, conf) {
  var d = cmpDeps_();
  var s = groupStats_(values);
  var r = { n: s.n, mean: s.mean, sd: s.sd, se: s.se, diff: null, t: null, df: null, p1: null, p2: null, lo: null, hi: null };
  if (s.n < 2) return r;
  r.diff = s.mean - mu;
  r.df = s.n - 1;
  if (!(s.se > 0)) return r;
  r.t = r.diff / s.se;
  r.p2 = d.pValueT(r.t, r.df);
  r.p1 = r.p2 / 2;
  var tc = tCritical_(conf, r.df);
  r.lo = r.diff - tc * s.se;
  r.hi = r.diff + tc * s.se;
  return r;
}

function tFromSe_(diff, se, df, conf) {
  var d = cmpDeps_();
  var r = { diff: diff, se: se, df: df, t: null, p1: null, p2: null, lo: null, hi: null };
  if (!(se > 0) || !(df > 0)) return r;
  r.t = diff / se;
  r.p2 = d.pValueT(r.t, df);
  r.p1 = r.p2 / 2;
  var tc = tCritical_(conf, df);
  r.lo = diff - tc * se;
  r.hi = diff + tc * se;
  return r;
}

/**
 * 독립표본 t검정 (집단1 − 집단2).
 * equal: 합동분산(등분산 가정), welch: Welch–Satterthwaite 자유도(등분산 가정 안 함), levene: 평균 기준 Levene.
 * 어느 집단이든 사례가 2개 미만이면 검정은 null.
 */
function independentT(a, b, conf) {
  var d = cmpDeps_();
  var sa = groupStats_(a), sb = groupStats_(b);
  var out = { a: sa, b: sb, equal: null, welch: null, levene: null };
  if (sa.n < 1 || sb.n < 1) return out;
  var diff = sa.mean - sb.mean;
  if (sa.n < 2 || sb.n < 2) {
    out.equal = tFromSe_(diff, NaN, NaN, conf);
    out.welch = tFromSe_(diff, NaN, NaN, conf);
    return out;
  }
  var df = sa.n + sb.n - 2;
  var sp2 = ((sa.n - 1) * sa.variance + (sb.n - 1) * sb.variance) / df;
  out.equal = tFromSe_(diff, Math.sqrt(sp2 * (1 / sa.n + 1 / sb.n)), df, conf);

  var va = sa.variance / sa.n, vb = sb.variance / sb.n;
  var se = Math.sqrt(va + vb);
  var wdf = (va + vb) * (va + vb) / (va * va / (sa.n - 1) + vb * vb / (sb.n - 1));
  out.welch = tFromSe_(diff, se, wdf, conf);
  out.levene = d.leveneTest([a, b], false);
  return out;
}

/** 대응표본 t검정 (x − y). 두 배열은 같은 사례 순서의 유효한 쌍. */
function pairedT(x, y, conf) {
  var diffs = x.map(function (v, i) { return v - y[i]; });
  return {
    x: groupStats_(x),
    y: groupStats_(y),
    corr: cmpDeps_().pearsonR(x, y),
    test: oneSampleT(diffs, 0, conf)
  };
}

/* ---------- 결과 표 도우미 ---------- */

var F3_ = { type: 'fixed', digits: 3 };

/** 효과크기 표의 행: e 는 Stats_Effect 의 결과. paired 이면 d_z, g_z, d_av 세 줄 */
function effectRows_(name, e, paired) {
  var d = cmpDeps_();
  var word = d.effectSizeWord;
  var rows = [
    { name: name, kind: paired ? "Cohen's d_z" : "Cohen's d", std: e.std, est: e.d, lo: e.dLo, hi: e.dHi, size: word(e.d) },
    { name: '', kind: paired ? "Hedges' g_z" : "Hedges' g", std: e.std, est: e.g, lo: e.gLo, hi: e.gHi, size: word(e.g) }
  ];
  if (paired && e.dav !== null && e.dav !== undefined) rows.push({ name: '', kind: "Cohen's d_av", std: e.stdAv, est: e.dav, lo: null, hi: null, size: word(e.dav) });
  return rows;
}

/** 효과크기 표 블록. approx: 신뢰구간이 근사인 변수가 있는지 */
function effectBlock_(rows, conf, notes, approx) {
  return {
    type: 'table', caption: '효과크기',
    columns: [
      { key: 'name', label: '' }, { key: 'kind', label: '지표' }, { key: 'std', label: '표준화 기준값', format: F3_ },
      { key: 'est', label: '추정값', format: F3_ }, { key: 'lo', label: conf + '% 신뢰구간 하한', format: F3_ }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: F3_ },
      { key: 'size', label: '크기(참고)' }
    ],
    rows: rows,
    notes: notes.concat([
      "Hedges' g = d × 보정계수 (1 − 3/(4·자유도 − 1)). 표본이 작을 때 d의 편향을 줄인 값이다. 신뢰구간은 비중심 t 분포로 구한 정확한 구간이다" + (approx ? ' (|t|가 매우 커서 일부는 정규 근사)' : '') + '.',
      '크기(참고)는 Cohen(1988)의 관행(0.2 작음, 0.5 중간, 0.8 큼)이며 절대적인 기준이 아니다. 분야와 맥락에 맞게 해석하세요.'
    ])
  };
}

/** blocks 에서 caption 이 prefix 로 시작하는 블록 바로 뒤에 새 블록을 넣는다 (없으면 맨 뒤) */
function insertAfterCaption_(blocks, prefix, block) {
  for (var i = 0; i < blocks.length; i++) {
    if (blocks[i].type === 'table' && blocks[i].caption.indexOf(prefix) === 0) { blocks.splice(i + 1, 0, block); return; }
  }
  blocks.push(block);
}

function statsColumns_(first) {
  return first.concat([
    { key: 'n', label: 'N', format: { type: 'int' } },
    { key: 'mean', label: '평균', format: F3_ },
    { key: 'sd', label: '표준편차', format: F3_ },
    { key: 'se', label: '평균의 표준오차', format: F3_ }
  ]);
}

/* ---------- 일표본 t검정 ---------- */

function runTTestOne(dataset, request) {
  var d = cmpDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  var mu = numberOption_(request.options, 'testValue', '검정값', 0);
  var conf = confidenceOption_(request.options);
  var warnings = [];
  var statRows = [], testRows = [], effRows = [], effApprox = false;
  var adv = { rows: [], recs: [] };

  vars.forEach(function (v) {
    var name = nameOf_(v);
    var values = dataset.data[v.name].filter(function (x) { return x !== null; });
    var r = oneSampleT(values, mu, conf);
    if (request.options && request.options.effect) {
      var e1 = d.effectOneSample({ n: r.n, mean: r.mean, sd: r.sd }, mu, conf);
      if (e1) { effRows = effRows.concat(effectRows_(name, e1, false)); effApprox = effApprox || e1.approx; }
      else warnings.push('"' + name + '": 효과크기를 계산할 수 없습니다 (유효한 값이 2개 미만이거나 표준편차가 0).');
    }
    if (r.t === null) warnings.push('"' + name + '": 유효한 값이 2개 미만이거나 표준편차가 0이라 검정할 수 없습니다.');
    statRows.push({ name: name, n: r.n, mean: r.mean, sd: r.sd, se: r.se });
    testRows.push({ name: name, t: r.t, df: r.df, p1: r.p1, p2: r.p2, diff: r.diff, lo: r.lo, hi: r.hi });
    if (request.options && request.options.assumptions) d.advAdd_(adv, d.advOneSample_(name, values, v.scale === 'ordinal'));
  });

  var result1 = {
    analysis: 'ttest1',
    title: '일표본 t검정',
    blocks: [
      { type: 'table', caption: '일표본 통계량', columns: statsColumns_([{ key: 'name', label: '변수' }]), rows: statRows, notes: [] },
      {
        type: 'table',
        caption: '일표본 검정 (검정값 = ' + mu + ')',
        columns: [
          { key: 'name', label: '변수' },
          { key: 't', label: 't', format: F3_ },
          { key: 'df', label: '자유도', format: { type: 'int' } },
          { key: 'p1', label: '유의확률(단측)', format: { type: 'pvalue' } },
          { key: 'p2', label: '유의확률(양측)', format: { type: 'pvalue' } },
          { key: 'diff', label: '평균차이', format: F3_ },
          { key: 'lo', label: ciLabel_(conf, '하한'), format: F3_ },
          { key: 'hi', label: ciLabel_(conf, '상한'), format: F3_ }
        ],
        rows: testRows,
        notes: ['평균차이 = 평균 − 검정값. 단측 유의확률은 양측의 절반이다.']
      }
    ],
    warnings: warnings
  };
  if (effRows.length) insertAfterCaption_(result1.blocks, '일표본 검정', effectBlock_(effRows, conf, ['d = (평균 − 검정값) / 표본 표준편차 (표준화 기준값 = 표본 표준편차).'], effApprox));
  if (adv.rows.length) result1.blocks = d.advAppend_(result1.blocks, adv, false);
  return result1;
}

/* ---------- 독립표본 t검정 ---------- */

function runTTestIndependent(dataset, request) {
  var d = cmpDeps_();
  var depVars = d.getVars(dataset, request.dependent, { numeric: true });
  var factorVars = d.getVars(dataset, request.factor);
  if (factorVars.length !== 1) throw new Error('집단변수는 하나만 선택하세요.');
  var factor = factorVars[0];
  depVars.forEach(function (v) {
    if (v.name === factor.name) throw new Error('종속변수와 집단변수에 같은 변수 "' + v.name + '"를 쓸 수 없습니다.');
  });
  var conf = confidenceOption_(request.options);
  var factorName = nameOf_(factor);
  var warnings = [];
  var groupRows = [], testRows = [], groupNotes = [], boxBlocks = [], effRows = [], effApprox = false;
  var adv = { rows: [], recs: [] };

  depVars.forEach(function (v) {
    var name = nameOf_(v);
    var split = d.splitGroups(dataset, v.name, factor);
    if (split.groups.length !== 2) {
      throw new Error('집단변수 "' + factorName + '"의 집단이 ' + split.groups.length + '개입니다 (종속변수 "' + name + '"에 유효한 사례 기준). 독립표본 t검정은 집단이 정확히 2개여야 합니다.');
    }
    var g1 = split.groups[0], g2 = split.groups[1];
    var r = independentT(g1.values, g2.values, conf);
    if (request.options && request.options.effect) {
      var e2 = d.effectIndependent(r.a, r.b, conf);
      if (e2) { effRows = effRows.concat(effectRows_(name, e2, false)); effApprox = effApprox || e2.approx; }
      else warnings.push('"' + name + '": 효과크기를 계산할 수 없습니다 (각 집단의 사례가 2개 이상이고 변동이 있어야 합니다).');
    }
    if (request.options && request.options.assumptions) {
      d.advAdd_(adv, d.advIndependentT_(name, [{ label: g1.label, values: g1.values }, { label: g2.label, values: g2.values }], v.scale === 'ordinal'));
    }
    if (request.options && request.options.boxplot) {
      boxBlocks.push(d.boxplotBlock('상자그림: ' + name + ' (집단변수: ' + factorName + ')',
        [{ label: g1.label, values: g1.values }, { label: g2.label, values: g2.values }], name));
    }

    [[g1, r.a], [g2, r.b]].forEach(function (pair) {
      groupRows.push({ name: name, group: pair[0].label, n: pair[1].n, mean: pair[1].mean, sd: pair[1].sd, se: pair[1].se });
    });
    var note = '평균차이 = 집단1 − 집단2 (집단1: ' + g1.label + ', 집단2: ' + g2.label + ')';
    if (groupNotes.indexOf(note) < 0) groupNotes.push(note);

    var lev = r.levene || { W: null, p: null };
    var eq = r.equal || {}, we = r.welch || {};
    if (!r.equal || eq.t === null) {
      warnings.push('"' + name + '": 각 집단의 사례가 2개 이상이고 변동이 있어야 검정할 수 있습니다.');
    }
    // Levene 유의확률 ≥ .05: 등분산이라는 귀무가설을 기각할 근거가 부족하므로 "등분산을 가정함" 행을 읽는다.
    // < .05: 등분산을 기각하므로 "등분산을 가정하지 않음"(Welch) 행을 읽는다. 읽을 행을 굵게 표시한다.
    var readEqual = lev.p !== null && lev.p !== undefined ? lev.p >= 0.05 : null;
    if (readEqual === null && eq.t !== null && eq.t !== undefined) {
      warnings.push('"' + name + '": Levene 검정을 할 수 없어(집단 내 변동이 없음) 어느 행을 읽을지 표시하지 않았습니다.');
    }
    testRows.push({
      _bold: readEqual === true, name: name, assume: '등분산을 가정함', levF: lev.W, levP: lev.p,
      t: eq.t, df: eq.df, p1: eq.p1, p2: eq.p2, diff: eq.diff, se: eq.se, lo: eq.lo, hi: eq.hi
    });
    testRows.push({
      _bold: readEqual === false, name: '', assume: '등분산을 가정하지 않음', levF: null, levP: null,
      t: we.t, df: we.df, p1: we.p1, p2: we.p2, diff: we.diff, se: we.se, lo: we.lo, hi: we.hi
    });
  });

  var result = {
    analysis: 'ttestInd',
    title: '독립표본 t검정',
    blocks: [
      {
        type: 'table',
        caption: '집단통계량 (집단변수: ' + factorName + ')',
        columns: statsColumns_([{ key: 'name', label: '종속변수' }, { key: 'group', label: '집단' }]),
        rows: groupRows,
        notes: []
      },
      {
        type: 'table',
        caption: '독립표본 검정',
        columns: [
          { key: 'name', label: '종속변수' },
          { key: 'assume', label: '' },
          { key: 'levF', label: 'Levene F', format: F3_ },
          { key: 'levP', label: 'Levene 유의확률', format: { type: 'pvalue' } },
          { key: 't', label: 't', format: F3_ },
          { key: 'df', label: '자유도', format: { type: 'trim', digits: 3 } },
          { key: 'p1', label: '유의확률(단측)', format: { type: 'pvalue' } },
          { key: 'p2', label: '유의확률(양측)', format: { type: 'pvalue' } },
          { key: 'diff', label: '평균차이', format: F3_ },
          { key: 'se', label: '차이의 표준오차', format: F3_ },
          { key: 'lo', label: ciLabel_(conf, '하한'), format: F3_ },
          { key: 'hi', label: ciLabel_(conf, '상한'), format: F3_ }
        ],
        rows: testRows,
        notes: groupNotes.concat([
          'Levene 검정은 평균 기준이다. Levene 유의확률이 .05 이상이면 등분산이라는 귀무가설을 기각할 근거가 부족하므로 "등분산을 가정함" 행을, .05 미만이면 "등분산을 가정하지 않음"(Welch) 행을 읽는다. 읽을 행을 굵게 표시했다.',
          '단측 유의확률은 양측의 절반이다. 종속변수와 집단변수가 모두 유효한 사례만 사용했다.'
        ])
      }
    ],
    warnings: warnings
  };
  if (effRows.length) insertAfterCaption_(result.blocks, '독립표본 검정', effectBlock_(effRows, conf, ['d = (집단1 평균 − 집단2 평균) / 합동 표준편차 (표준화 기준값 = 합동 표준편차). 방향은 평균차이와 같다 (집단1 − 집단2). 등분산이 크게 어긋난 경우에도 합동 표준편차를 쓰는 값이다.'], effApprox));
  if (adv.rows.length) result.blocks = d.advAppend_(result.blocks, adv, true);
  result.blocks = result.blocks.concat(boxBlocks);
  return result;
}

/* ---------- 대응표본 t검정 ---------- */

function runTTestPaired(dataset, request) {
  var d = cmpDeps_();
  var first = d.getVars(dataset, request.first, { numeric: true });
  var second = d.getVars(dataset, request.second, { numeric: true });
  if (first.length !== 1 || second.length !== 1) throw new Error('변수1과 변수2에 변수를 하나씩 선택하세요.');
  var v1 = first[0], v2 = second[0];
  if (v1.name === v2.name) throw new Error('변수1과 변수2에 같은 변수 "' + v1.name + '"를 쓸 수 없습니다.');
  var conf = confidenceOption_(request.options);
  var n1 = nameOf_(v1), n2 = nameOf_(v2);

  var x = [], y = [];
  for (var i = 0; i < dataset.n; i++) {
    var a = dataset.data[v1.name][i], b = dataset.data[v2.name][i];
    if (a !== null && b !== null) { x.push(a); y.push(b); }
  }
  var r = pairedT(x, y, conf);
  var t = r.test;
  var warnings = [];
  var adv = { rows: [], recs: [] };
  var effRows = [], effApprox = false;
  if (request.options && request.options.effect) {
    var e3 = d.effectPaired({ n: t.n, mean: t.mean, sd: t.sd }, r.x.sd, r.y.sd, conf);
    if (e3) { effRows = effectRows_(n1 + ' − ' + n2, e3, true); effApprox = e3.approx; }
    else warnings.push('효과크기를 계산할 수 없습니다 (유효한 쌍이 2개 미만이거나 차이의 표준편차가 0).');
  }
  if (request.options && request.options.assumptions) {
    var diffs = x.map(function (xv, i) { return xv - y[i]; });
    d.advAdd_(adv, d.advPaired_(n1 + ' − ' + n2, diffs, v1.scale === 'ordinal' || v2.scale === 'ordinal'));
  }
  if (t.t === null) warnings.push('유효한 쌍이 2개 미만이거나 차이의 표준편차가 0이라 검정할 수 없습니다.');
  var dropped = dataset.n - x.length;
  var pairName = n1 + ' − ' + n2;

  var resultP = {
    analysis: 'ttestPaired',
    title: '대응표본 t검정',
    blocks: [
      {
        type: 'table',
        caption: '대응표본 통계량',
        columns: [
          { key: 'name', label: '변수' },
          { key: 'mean', label: '평균', format: F3_ },
          { key: 'n', label: 'N', format: { type: 'int' } },
          { key: 'sd', label: '표준편차', format: F3_ },
          { key: 'se', label: '평균의 표준오차', format: F3_ }
        ],
        rows: [
          { name: n1, mean: r.x.mean, n: r.x.n, sd: r.x.sd, se: r.x.se },
          { name: n2, mean: r.y.mean, n: r.y.n, sd: r.y.sd, se: r.y.se }
        ],
        notes: []
      },
      {
        type: 'table',
        caption: '대응표본 상관계수',
        columns: [
          { key: 'name', label: '쌍' },
          { key: 'n', label: 'N', format: { type: 'int' } },
          { key: 'r', label: '상관계수', format: F3_ },
          { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
        ],
        rows: [{ name: n1 + ' & ' + n2, n: r.corr.n, r: r.corr.r, p: r.corr.p }],
        notes: []
      },
      {
        type: 'table',
        caption: '대응표본 검정',
        columns: [
          { key: 'name', label: '쌍' },
          { key: 'diff', label: '평균', format: F3_ },
          { key: 'sd', label: '표준편차', format: F3_ },
          { key: 'se', label: '평균의 표준오차', format: F3_ },
          { key: 'lo', label: ciLabel_(conf, '하한'), format: F3_ },
          { key: 'hi', label: ciLabel_(conf, '상한'), format: F3_ },
          { key: 't', label: 't', format: F3_ },
          { key: 'df', label: '자유도', format: { type: 'int' } },
          { key: 'p1', label: '유의확률(단측)', format: { type: 'pvalue' } },
          { key: 'p2', label: '유의확률(양측)', format: { type: 'pvalue' } }
        ],
        rows: [{ name: pairName, diff: t.mean, sd: t.sd, se: t.se, lo: t.lo, hi: t.hi, t: t.t, df: t.df, p1: t.p1, p2: t.p2 }],
        notes: [
          '평균 = 변수1 − 변수2 의 평균. 단측 유의확률은 양측의 절반이다.',
          '두 변수가 모두 유효한 쌍만 사용했다' + (dropped > 0 ? ' (' + dropped + '건 제외).' : '.')
        ]
      }
    ],
    warnings: warnings
  };
  if (effRows.length) insertAfterCaption_(resultP.blocks, '대응표본 검정', effectBlock_(effRows, conf, ["d_z = 평균 차이 / 차이의 표준편차 (대응표본의 Cohen's d, 표준화 기준값 = 차이의 표준편차). d_av = 평균 차이 / 두 측정 표준편차의 평균 (표준화 기준값 = 그 평균, 다른 연구와 비교할 때 권장되며 신뢰구간은 제공하지 않는다). 방향은 변수1 − 변수2."], effApprox));
  if (adv.rows.length) resultP.blocks = d.advAppend_(resultP.blocks, adv, false);
  return resultP;
}

if (typeof module !== 'undefined') {
  module.exports = {
    oneSampleT: oneSampleT,
    independentT: independentT,
    pairedT: pairedT,
    runTTestOne: runTTestOne,
    runTTestIndependent: runTTestIndependent,
    runTTestPaired: runTTestPaired
  };
}
