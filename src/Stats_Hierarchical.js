/**
 * 위계적(블록) 회귀. 순수 함수: (dataset, request) → 결과 객체.
 * 독립변수를 1~4단계로 나눠 단계마다 변수를 누적해 OLS 를 적합하고, 단계 사이의 R² 변화와 F 변화 검정을 낸다.
 * 모든 적합은 선형 회귀의 olsFit 을 그대로 쓴다. 표 구성은 블록 회귀의 모형 요약(변화량 통계량)과 같다.
 *
 * 변수 처리: 명목 척도 변수는 범주형으로 보고 더미(기준 범주 제외)로 바꿔 한 단계에 함께 넣는다.
 * 숫자·서열 척도 변수는 숫자 점수로 넣는다.
 *
 * 변화 F: F = ((SSE_이전 − SSE_현재) / 추가한 계수 수) / (SSE_현재 / 현재 잔차 자유도), 자유도 (추가한 계수 수, 현재 잔차 자유도).
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function hierDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var reg = require('./Stats_Regression.js');
    return {
      getVars: data.getVars, pf: dist.pf, olsFit: reg.olsFit, regDummies_: reg.regDummies_,
      regCollinearityBlock_: reg.regCollinearityBlock_, regName_: reg.regName_
    };
  }
  return {
    getVars: getVars, pf: pf, olsFit: olsFit, regDummies_: regDummies_,
    regCollinearityBlock_: regCollinearityBlock_, regName_: regName_
  };
}

var HIER_MAX_BLOCKS_ = 4;

/**
 * 단계별 적합 결과에 변화량 통계량을 붙인다. fits: [olsFit 결과] (1단계부터). 반환: [{ dR2, dF, df1, df2, pF }]
 * 1단계의 변화량은 모형 자체의 R², F 이다.
 */
function hierChanges(fits) {
  var d = hierDeps_();
  return fits.map(function (fit, k) {
    if (k === 0) return { dR2: fit.r2, dF: fit.F, df1: fit.dfm, df2: fit.dfr, pF: fit.pF };
    var prev = fits[k - 1];
    var m = fit.p - prev.p;
    var dR2 = fit.r2 === null || prev.r2 === null ? null : Math.max(0, fit.r2 - prev.r2);
    var F = null, pF = null;
    if (fit.mse > 0) {
      F = (Math.max(0, prev.sse - fit.sse) / m) / fit.mse;
      pF = d.pf(F, m, fit.dfr, false);
    }
    return { dR2: dR2, dF: F, df1: m, df2: fit.dfr, pF: pF };
  });
}

function runHierarchical(dataset, request) {
  var d = hierDeps_();
  var dep = d.getVars(dataset, request.dependent, { numeric: true });
  if (dep.length !== 1) throw new Error('종속변수는 하나만 선택하세요.');

  // 단계별 변수 이름: 앞에서부터 채워야 한다
  var names = [], b, last = 0;
  for (b = 1; b <= HIER_MAX_BLOCKS_; b++) {
    var list = request['block' + b] || [];
    names.push(list);
    if (list.length) last = b;
  }
  if (!names[0].length) throw new Error('1단계에 독립변수를 하나 이상 선택하세요.');
  for (b = 2; b <= last; b++) {
    if (!names[b - 1].length) throw new Error(b + '단계가 비어 있는데 더 뒤 단계에 변수가 있습니다. 단계를 앞에서부터 채우세요.');
  }
  names = names.slice(0, last);

  var seen = {};
  seen[dep[0].name] = '종속변수';
  var blockVars = names.map(function (list, bi) {
    var vars = d.getVars(dataset, list);
    vars.forEach(function (v) {
      if (seen[v.name]) {
        throw new Error('변수 "' + v.name + '"를 ' + seen[v.name] + '와(과) ' + (bi + 1) + '단계에 동시에 쓸 수 없습니다.');
      }
      seen[v.name] = (bi + 1) + '단계 독립변수';
      if (v.scale !== 'nominal') d.getVars(dataset, [v.name], { numeric: true });
    });
    return vars;
  });

  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');

  // 모든 단계의 변수 기준으로 목록별 결측 제거
  var all = dep.slice();
  blockVars.forEach(function (vars) { all = all.concat(vars); });
  var rowsOk = [];
  for (var i = 0; i < dataset.n; i++) {
    var ok = true;
    for (var a = 0; a < all.length; a++) if (dataset.data[all[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length;
  var pick = function (name) { return rowsOk.map(function (r) { return dataset.data[name][r]; }); };
  var y = pick(dep[0].name);

  // 단계마다 새로 들어가는 열
  var refNotes = [], ordinals = [];
  var added = blockVars.map(function (vars) {
    var cols = [];
    vars.forEach(function (v) {
      if (v.scale === 'nominal') {
        var dm = d.regDummies_(v, pick(v.name), opt.lastRef);
        refNotes.push(dm.note);
        dm.predictors.forEach(function (pr) { cols.push(pr); });
      } else {
        if (v.scale === 'ordinal') ordinals.push(d.regName_(v));
        cols.push({ name: d.regName_(v), values: pick(v.name) });
      }
    });
    return cols;
  });

  var fits = [], cum = [], predsByStep = [];
  added.forEach(function (cols, bi) {
    cum = cum.concat(cols);
    predsByStep.push(cum);
    var fit = d.olsFit(cum, y, conf);
    var label = (bi + 1) + '단계';
    if (fit.tooFew) throw new Error('사례 수(' + fit.n + ')가 ' + label + ' 모형의 계수의 수(' + fit.p + ', 상수 포함)보다 많아야 합니다.');
    if (fit.constant) throw new Error('변수 "' + fit.constant + '"의 값이 모두 같아 회귀에 쓸 수 없습니다.');
    if (fit.collinear) throw new Error('변수 "' + fit.collinear + '"이(가) 앞선 독립변수들의 선형결합이라 ' + label + ' 모형을 계산할 수 없습니다 (완전 공선성). 이 변수를 빼 보세요.');
    fits.push(fit);
  });
  var ch = hierChanges(fits);
  var lastFit = fits[fits.length - 1];

  var warnings = [];
  if (fits.length === 1) warnings.push('1단계만 있어 일반 선형 회귀와 같습니다. 2단계 이상에 변수를 보내면 단계별 변화를 비교합니다.');
  fits.forEach(function (f, bi) {
    if (f.mse === 0) warnings.push((bi + 1) + '단계 모형의 잔차가 모두 0입니다 (완전한 적합). 표준오차, t, F를 계산할 수 없습니다.');
  });
  lastFit.coef.forEach(function (c) {
    if (c.vif !== null && c.vif > 10) {
      warnings.push('마지막 단계에서 "' + c.name + '"의 VIF가 ' + c.vif.toFixed(1) + '로 10을 넘습니다. 독립변수 사이의 공선성이 큽니다.');
    }
  });
  var outliers = lastFit.s > 0 ? lastFit.resid.filter(function (e) { return Math.abs(e / lastFit.s) > 3; }).length : 0;

  var F3 = { type: 'fixed', digits: 3 };
  var tables = [];

  var sumCols = [
    { key: 'step', label: '' },
    { key: 'r', label: 'R', format: F3 },
    { key: 'r2', label: 'R 제곱', format: F3 },
    { key: 'adj', label: '수정된 R 제곱', format: F3 },
    { key: 's', label: '추정값의 표준오차', format: F3 },
    { key: 'dR2', label: 'R 제곱 변화량', format: F3 },
    { key: 'dF', label: 'F 변화량', format: F3 },
    { key: 'df1', label: '자유도1', format: { type: 'int' } },
    { key: 'df2', label: '자유도2', format: { type: 'int' } },
    { key: 'pF', label: '유의확률 F 변화량', format: { type: 'pvalue' } }
  ];
  if (opt.dw) sumCols.push({ key: 'dw', label: '더빈-왓슨', format: F3 });
  var sumRows = fits.map(function (f, bi) {
    var c = ch[bi];
    var row = { step: (bi + 1) + '단계', r: f.r, r2: f.r2, adj: f.adjR2, s: f.s, dR2: c.dR2, dF: c.dF, df1: c.df1, df2: c.df2, pF: c.pF };
    if (opt.dw) row.dw = f.dw;
    return row;
  });
  var sumNotes = added.map(function (cols, bi) {
    return (bi === 0 ? '1단계 투입: ' : (bi + 1) + '단계 추가: ') + cols.map(function (c) { return c.name; }).join(', ');
  });
  sumNotes.push('각 단계는 앞 단계의 변수를 모두 포함한다. 변화량 통계량은 그 단계에서 추가한 변수들이 설명력을 유의하게 높였는지 본다 (1단계는 모형 전체의 값). 유의확률이 .05보다 작으면 유의한 증가이다.');
  sumNotes.push('독립변수: 단계마다 누적. 유효 사례 ' + n + '건' + (dataset.n > n ? ' (결측으로 ' + (dataset.n - n) + '건 제외).' : '.'));
  if (ordinals.length) sumNotes.push('서열 척도 변수(' + ordinals.join(', ') + ')는 숫자 점수로 투입했다. 범주로 다루려면 변수 척도를 명목으로 바꾸세요.');
  tables.push({ type: 'table', caption: '모형 요약 (종속변수: ' + d.regName_(dep[0]) + ')', columns: sumCols, rows: sumRows, notes: sumNotes });

  var anovaRows = [];
  fits.forEach(function (f, bi) {
    var label = (bi + 1) + '단계';
    anovaRows.push({ model: label, source: '회귀', ss: f.ssr, df: f.dfm, ms: f.msr, F: f.F, p: f.pF });
    anovaRows.push({ model: '', source: '잔차', ss: f.sse, df: f.dfr, ms: f.mse });
    anovaRows.push({ _kind: 'total', model: '', source: '합계', ss: f.sst, df: f.dft });
  });
  tables.push({
    type: 'table', caption: '분산분석 (단계별 모형 전체)',
    columns: [
      { key: 'model', label: '' },
      { key: 'source', label: '' },
      { key: 'ss', label: '제곱합', format: F3 },
      { key: 'df', label: '자유도', format: { type: 'int' } },
      { key: 'ms', label: '평균제곱', format: F3 },
      { key: 'F', label: 'F', format: F3 },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
    ],
    rows: anovaRows,
    notes: ['단계별 모형 전체가 종속변수를 설명하는지 보는 검정이다 (앞 단계와의 차이는 모형 요약의 변화량 통계량).']
  });

  var coefCols = [
    { key: 'model', label: '' },
    { key: 'name', label: '' },
    { key: 'B', label: 'B', format: F3 },
    { key: 'se', label: '표준오차', format: F3 },
    { key: 'beta', label: '베타', format: F3 },
    { key: 't', label: 't', format: F3 },
    { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
  ];
  if (opt.ci) {
    coefCols.push({ key: 'lo', label: 'B의 ' + conf + '% 신뢰구간 하한', format: F3 });
    coefCols.push({ key: 'hi', label: 'B의 ' + conf + '% 신뢰구간 상한', format: F3 });
  }
  coefCols.push({ key: 'tol', label: '공차', format: F3 });
  coefCols.push({ key: 'vif', label: 'VIF', format: F3 });
  var coefRows = [];
  fits.forEach(function (f, bi) {
    f.coef.forEach(function (c, ci) {
      coefRows.push(Object.assign({ model: ci === 0 ? (bi + 1) + '단계' : '' }, c));
    });
  });
  var coefNotes = ['베타는 표준화 계수이고, 공차 = 1/VIF 이다. 상수에는 베타·공차·VIF가 없다. 계수는 단계마다 그 단계의 모형에서 구한 값이다.'].concat(refNotes);
  if (outliers > 0) coefNotes.push('마지막 단계에서 표준화 잔차(잔차/추정값의 표준오차)의 절댓값이 3을 넘는 사례가 ' + outliers + '건 있다.');
  tables.push({ type: 'table', caption: '계수', columns: coefCols, rows: coefRows, notes: coefNotes });

  if (opt.collinearity) {
    tables.push(d.regCollinearityBlock_(predsByStep[predsByStep.length - 1], n, lastFit.coef, '공선성 진단 (마지막 단계 모형)'));
  }

  return { analysis: 'hierarchical', title: '위계적 회귀', blocks: tables, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = {
    hierChanges: hierChanges,
    runHierarchical: runHierarchical
  };
}
