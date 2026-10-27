/**
 * 조절효과 분석 (Hayes PROCESS 모델 1에 해당): Y = b0 + b1·X + b2·W + b3·X·W (+ 공변량).
 * 순수 함수: (dataset, request) → 결과 객체. 회귀 계산은 Stats_Regression.js의 olsFit(QR)을 쓴다.
 *
 * - 독립변수 X는 숫자, 조절변수 W는 숫자(연속) 또는 범주형(척도가 명목, 더미 변환). 공변량은 숫자.
 * - 평균 중심화(기본): X와 연속형 W를 각각 평균으로 중심화한 뒤 곱한다. 그러면 b1은 W가 평균일 때(범주형이면 기준 범주일 때) X의 효과이다.
 * - 상호작용 검정: 상호작용항을 넣었을 때의 R² 증가량(ΔR²)과 F 변화량.
 * - 단순기울기: 조절변수 수준별 X의 효과. 표준오차는 계수 공분산으로 √(V11 + 2w·V13 + w²·V33).
 * - Johnson-Neyman: X의 효과가 유의해지는 조절변수의 구간 (연속형 W).
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function modDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    var desc = require('./Stats_Descriptive.js');
    var reg = require('./Stats_Regression.js');
    return {
      getVars: data.getVars, compareValues_: data.compareValues_,
      pf: dist.pf, pt: dist.pt, qt: dist.qt, centralSums: desc.centralSums,
      olsFit: reg.olsFit, lineChartBlock: require('./Charts.js').lineChartBlock
    };
  }
  return {
    getVars: getVars, compareValues_: compareValues_, pf: pf, pt: pt, qt: qt, centralSums: centralSums,
    olsFit: olsFit, lineChartBlock: lineChartBlock
  };
}

var MOD_MAX_CATEGORIES_ = 10;
function modName_(v) { return v.label ? v.label : v.name; }

/**
 * Johnson-Neyman 구간. 모형 좌표의 계수 b1(X), b3(X·W), 분산·공분산 v11, v13, v33, 임계 t, 조절변수의 관측 범위 [lo, hi](모형 좌표).
 * 반환: { roots: [범위 안의 경계값], regions: [{ lo, hi, sig, sign }] } (sig: X의 효과가 유의한 구간인지)
 * 조건: (b1 + b3·w)² = t²·(v11 + 2w·v13 + w²·v33) 의 이차방정식 풀이.
 */
function johnsonNeyman(b1, b3, v11, v13, v33, tcrit, lo, hi) {
  var A = b3 * b3 - tcrit * tcrit * v33;
  var B = 2 * (b1 * b3 - tcrit * tcrit * v13);
  var C = b1 * b1 - tcrit * tcrit * v11;
  var roots = [];
  var scale = Math.max(Math.abs(A), Math.abs(B), Math.abs(C), 1e-300);
  if (Math.abs(A) <= 1e-12 * scale) {
    if (Math.abs(B) > 1e-12 * scale) roots.push(-C / B);
  } else {
    var disc = B * B - 4 * A * C;
    if (disc >= 0) {
      var sq = Math.sqrt(disc);
      roots.push((-B - sq) / (2 * A), (-B + sq) / (2 * A));
    }
  }
  roots = roots.filter(function (r) { return r > lo && r < hi; }).sort(function (a, b) { return a - b; });
  var cuts = [lo].concat(roots, [hi]);
  var regions = [];
  for (var i = 0; i + 1 < cuts.length; i++) {
    var mid = (cuts[i] + cuts[i + 1]) / 2;
    var slope = b1 + b3 * mid, v = v11 + 2 * mid * v13 + mid * mid * v33;
    var sig = v > 0 && Math.abs(slope) / Math.sqrt(v) > tcrit;
    regions.push({ lo: cuts[i], hi: cuts[i + 1], sig: sig, sign: slope >= 0 ? 1 : -1 });
  }
  return { roots: roots, regions: regions };
}

/**
 * 핵심 계산. input = {
 *   y: [..], x: [..], xName,
 *   w: { type: 'continuous', values, name } | { type: 'categorical', codes: [범주 번호], labels: [범주 레이블], refIdx, name },
 *   covs: [{ name, values }], center: boolean, conf: 95
 * } (모든 배열은 결측이 없는 같은 길이)
 * 반환: { fit, reduced, predictors, interactionIdx, deltaR2, F, df1, df2, pF, slopes, jn, plot, xMean, wMean, wSd } 또는 { error }
 */
function moderationFit(input) {
  var d = modDeps_();
  var n = input.y.length, i, j;
  var conf = input.conf;
  var xs = d.centralSums(input.x);
  var xOff = input.center ? xs.mean : 0;
  var xc = input.x.map(function (v) { return v - xOff; });
  var predictors = [{ name: input.xName, values: xc }];
  var interaction = [];
  var wInfo = { type: input.w.type };

  if (input.w.type === 'continuous') {
    var ws = d.centralSums(input.w.values);
    var wOff = input.center ? ws.mean : 0;
    var wc = input.w.values.map(function (v) { return v - wOff; });
    wInfo.mean = ws.mean; wInfo.sd = ws.n > 1 ? Math.sqrt(ws.c2 / (ws.n - 1)) : 0; wInfo.off = wOff;
    predictors.push({ name: input.w.name, values: wc });
    interaction.push({ name: input.xName + ' × ' + input.w.name, values: xc.map(function (v, k) { return v * wc[k]; }) });
  } else {
    wInfo.labels = input.w.labels; wInfo.refIdx = input.w.refIdx;
    wInfo.dummyOf = [];
    input.w.labels.forEach(function (lab, c) {
      if (c === input.w.refIdx) return;
      var dv = input.w.codes.map(function (code) { return code === c ? 1 : 0; });
      predictors.push({ name: input.w.name + ' = ' + lab, values: dv });
      interaction.push({ name: input.xName + ' × ' + input.w.name + ' = ' + lab, values: xc.map(function (v, k) { return v * dv[k]; }) });
      wInfo.dummyOf.push(c);
    });
  }
  var nMain = predictors.length;                               // X, W(또는 더미들)
  var covCols = input.covs.map(function (c) { return { name: c.name, values: c.values }; });
  var full = predictors.concat(interaction, covCols);
  var reducedP = predictors.concat(covCols);
  var fit = d.olsFit(full, input.y, conf);
  var reduced = d.olsFit(reducedP, input.y, conf);
  var bad = fit.tooFew || fit.constant || fit.collinear ? fit : (reduced.tooFew || reduced.constant || reduced.collinear ? reduced : null);
  if (bad) return { error: bad };

  var q = interaction.length;
  var deltaR2 = fit.r2 !== null && reduced.r2 !== null ? fit.r2 - reduced.r2 : null;
  var F = null, pF = null;
  if (deltaR2 !== null && fit.mse > 0) {
    F = (deltaR2 / q) / ((1 - fit.r2) / fit.dfr);
    pF = d.pf(F, q, fit.dfr, false);
  }
  var res = {
    fit: fit, reduced: reduced, names: full.map(function (p) { return p.name; }), nMain: nMain, q: q,
    deltaR2: deltaR2, F: F, df1: q, df2: fit.dfr, pF: pF, wInfo: wInfo,
    xMean: xs.mean, xOff: xOff, conf: conf
  };

  // 단순기울기: 조절변수 수준별 X의 효과. 계수 번호: B[0] 상수, B[1+j] 는 j 번째 설계 열 (X = 0)
  var B = fit.B, V = fit.slopeCov;
  var tc = d.qt((1 - conf / 100) / 2, fit.dfr, false);
  function slopeRow(label, wValue, b, v) {
    var se = v !== null && v >= 0 ? Math.sqrt(v) : null;
    var t = se !== null && se > 0 ? b / se : null;
    return {
      label: label, w: wValue, b: b, se: se, t: t,
      p: t !== null ? d.pt(-Math.abs(t), fit.dfr) * 2 : null,
      lo: se !== null ? b - tc * se : null, hi: se !== null ? b + tc * se : null
    };
  }
  var slopes = [];
  var ix = 0, iw = nMain, iInt = nMain;                        // 공분산 행렬에서 X, 첫 상호작용항의 위치 (설계 열 번호)
  if (input.w.type === 'continuous') {
    var iXW = nMain;                                            // X, W 다음이 X·W
    var levels = [['평균 − 1SD', -1], ['평균', 0], ['평균 + 1SD', 1]];
    levels.forEach(function (lv) {
      var off = lv[1] * wInfo.sd;                              // 평균에서의 거리
      var wm = (input.center ? 0 : wInfo.mean) + off;           // 모형 좌표의 W 값
      var b = B[1 + ix] + B[1 + iXW] * wm;
      var v = V ? V[ix][ix] + 2 * wm * V[ix][iXW] + wm * wm * V[iXW][iXW] : null;
      slopes.push(slopeRow(lv[0], wInfo.mean + off, b, v));
    });
    res.wCoordOffset = input.center ? wInfo.mean : 0;
    res.iXW = iXW;
  } else {
    var dummyPos = {};
    wInfo.dummyOf.forEach(function (c, k) { dummyPos[c] = k; });
    input.w.labels.forEach(function (lab, c) {
      if (c === input.w.refIdx) {
        slopes.push(slopeRow(lab + ' (기준)', null, B[1 + ix], V ? V[ix][ix] : null));
      } else {
        var pos = nMain + dummyPos[c];                          // 상호작용항의 설계 열 번호
        slopes.push(slopeRow(lab, null, B[1 + ix] + B[1 + pos], V ? V[ix][ix] + 2 * V[ix][pos] + V[pos][pos] : null));
      }
    });
  }
  res.slopes = slopes;

  // Johnson-Neyman (연속형 W)
  res.jn = null;
  if (input.w.type === 'continuous' && V) {
    var wvals = input.w.values;
    var wmin = Infinity, wmax = -Infinity;
    wvals.forEach(function (v) { if (v < wmin) wmin = v; if (v > wmax) wmax = v; });
    var shift = input.center ? wInfo.mean : 0;                  // 모형 좌표 = 원래 값 − shift
    var jn = johnsonNeyman(B[1 + ix], B[1 + res.iXW], V[ix][ix], V[ix][res.iXW], V[res.iXW][res.iXW], tc, wmin - shift, wmax - shift);
    var regions = jn.regions.map(function (r) {
      var lo = r.lo + shift, hi = r.hi + shift;
      var cnt = wvals.filter(function (v) { return v >= lo && v <= hi; }).length;
      return { lo: lo, hi: hi, sig: r.sig, sign: r.sign, pct: cnt / n * 100 };
    });
    // 경계에 걸린 사례가 두 구간에 겹쳐 세어지지 않게, 마지막 구간을 제외한 구간은 hi 미만만 센다
    regions.forEach(function (r, k) {
      if (k < regions.length - 1) r.pct = wvals.filter(function (v) { return v >= r.lo && v < r.hi; }).length / n * 100;
    });
    res.jn = { roots: jn.roots.map(function (r) { return r + shift; }), regions: regions, min: wmin, max: wmax };
  }

  // 그림 자료: 다른 변수는 평균(공변량은 원래 평균)에 둔다
  var xmin = Infinity, xmax = -Infinity;
  input.x.forEach(function (v) { if (v < xmin) xmin = v; if (v > xmax) xmax = v; });
  var covMeans = input.covs.map(function (c) { return d.centralSums(c.values).mean; });
  function predict(xRaw, dummyIdx, wModel) {
    var xm = xRaw - xOff;
    var y = B[1] * 0 + B[0] + B[1] * xm;
    if (input.w.type === 'continuous') y += B[2] * wModel + B[3] * xm * wModel;
    else if (dummyIdx !== null) y += B[1 + 1 + dummyIdx] + B[1 + nMain + dummyIdx] * xm;
    for (var c = 0; c < covMeans.length; c++) y += B[1 + nMain + q + c] * covMeans[c];
    return y;
  }
  var series = [];
  if (input.w.type === 'continuous') {
    slopes.forEach(function (s) {
      var wModel = s.w - (input.center ? wInfo.mean : 0);
      series.push({ label: s.label + ' (' + Number(s.w.toFixed(2)) + ')', points: [[xmin, predict(xmin, null, wModel)], [xmax, predict(xmax, null, wModel)]] });
    });
  } else {
    input.w.labels.forEach(function (lab, c) {
      var di = c === input.w.refIdx ? null : wInfo.dummyOf.indexOf(c);
      series.push({ label: lab, points: [[xmin, predict(xmin, di, 0)], [xmax, predict(xmax, di, 0)]] });
    });
  }
  res.plot = { series: series, xmin: xmin, xmax: xmax };
  return res;
}

function runModeration(dataset, request) {
  var d = modDeps_();
  var dep = d.getVars(dataset, request.dependent, { numeric: true });
  var ind = d.getVars(dataset, request.independent, { numeric: true });
  var modv = d.getVars(dataset, request.moderator);
  var covs = request.covariates && request.covariates.length ? d.getVars(dataset, request.covariates, { numeric: true }) : [];
  if (dep.length !== 1) throw new Error('종속변수는 하나만 선택하세요.');
  if (ind.length !== 1) throw new Error('독립변수는 하나만 선택하세요.');
  if (modv.length !== 1) throw new Error('조절변수는 하나만 선택하세요.');
  var seen = {};
  seen[dep[0].name] = '종속변수';
  [[ind[0], '독립변수'], [modv[0], '조절변수']].concat(covs.map(function (v) { return [v, '공변량']; })).forEach(function (p) {
    if (seen[p[0].name]) throw new Error('변수 "' + p[0].name + '"를 ' + seen[p[0].name] + '와(과) ' + p[1] + '에 동시에 쓸 수 없습니다.');
    seen[p[0].name] = p[1];
  });
  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');
  var center = opt.center !== false;                              // 기본 켬 (서버는 명시적으로 false 일 때만 끈다)

  var W = modv[0];
  var categorical = W.scale === 'nominal';
  if (!categorical) d.getVars(dataset, [W.name], { numeric: true });   // 연속형이면 숫자여야 한다

  var all = dep.concat(ind, [W], covs);
  var rowsOk = [];
  for (var i = 0; i < dataset.n; i++) {
    var ok = true;
    for (var a = 0; a < all.length; a++) if (dataset.data[all[a].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var n = rowsOk.length;
  var pick = function (name) { return rowsOk.map(function (r) { return dataset.data[name][r]; }); };

  var input = {
    y: pick(dep[0].name), x: pick(ind[0].name), xName: modName_(ind[0]),
    covs: covs.map(function (v) { return { name: modName_(v), values: pick(v.name) }; }),
    center: center, conf: conf
  };
  var cats = null;
  if (categorical) {
    var vals = pick(W.name), seenC = {};
    cats = [];
    vals.forEach(function (x) {
      var key = typeof x + ':' + x;
      if (!(key in seenC)) { seenC[key] = cats.length; cats.push({ value: x }); }
    });
    cats.sort(d.compareValues_);
    cats.forEach(function (c) {
      var label = W.valueLabels[String(c.value)];
      c.label = label !== undefined ? label : String(c.value);
    });
    if (cats.length < 2) throw new Error('범주형 조절변수 "' + modName_(W) + '"의 범주가 ' + cats.length + '개입니다. 2개 이상이어야 합니다.');
    if (cats.length > MOD_MAX_CATEGORIES_) throw new Error('범주형 조절변수 "' + modName_(W) + '"의 범주가 ' + cats.length + '개입니다 (최대 ' + MOD_MAX_CATEGORIES_ + '개).');
    var refIdx = opt.lastRef ? cats.length - 1 : 0;
    input.w = {
      type: 'categorical', name: modName_(W), labels: cats.map(function (c) { return c.label; }), refIdx: refIdx,
      codes: vals.map(function (x) { for (var c = 0; c < cats.length; c++) if (cats[c].value === x) return c; return -1; })
    };
  } else {
    input.w = { type: 'continuous', name: modName_(W), values: pick(W.name) };
  }

  var r = moderationFit(input);
  if (r.error) {
    var e = r.error;
    if (e.tooFew) throw new Error('사례 수(' + e.n + ')가 계수의 수(' + e.p + ', 상수 포함)보다 많아야 합니다.');
    if (e.constant) throw new Error('변수 "' + e.constant + '"의 값이 모두 같아 분석에 쓸 수 없습니다.');
    throw new Error('변수 "' + e.collinear + '"이(가) 다른 변수들의 선형결합이라 계산할 수 없습니다 (완전 공선성). 범주가 너무 작거나 독립변수와 조절변수가 같은 정보를 담고 있는지 확인하세요.');
  }

  var fit = r.fit, F3 = { type: 'fixed', digits: 3 };
  var warnings = [];
  if (dataset.n > n) warnings.push('목록별 결측 제거로 ' + (dataset.n - n) + '건을 제외했습니다 (사용한 사례 ' + n + '건).');
  if (fit.mse === 0) warnings.push('잔차가 모두 0입니다 (완전한 적합). 표준오차, t, F를 계산할 수 없습니다.');
  if (n < 30) warnings.push('사례가 ' + n + '건으로 적어 상호작용 검정의 검정력이 낮습니다.');
  if (!center) {
    fit.coef.forEach(function (c) {
      if (c.vif !== null && c.vif > 10) warnings.push('"' + c.name + '"의 VIF가 ' + c.vif.toFixed(1) + '로 큽니다. 평균 중심화 옵션을 켜면 줄어듭니다.');
    });
  }
  if (r.pF !== null && r.pF >= 0.05) warnings.push('상호작용(조절효과)이 유의하지 않습니다 (ΔR² 검정 p = ' + (r.pF < 0.001 ? '<.001' : r.pF.toFixed(3).replace(/^0/, '')) + '). 단순기울기는 참고로만 보세요.');

  var blocks = [];
  blocks.push({
    type: 'table',
    caption: '모형 요약 (종속변수: ' + modName_(dep[0]) + ')',
    columns: [
      { key: 'r', label: 'R', format: F3 }, { key: 'r2', label: 'R 제곱', format: F3 }, { key: 'adj', label: '수정된 R 제곱', format: F3 },
      { key: 's', label: '추정값의 표준오차', format: F3 }, { key: 'F', label: 'F', format: F3 },
      { key: 'df1', label: '자유도1', format: { type: 'int' } }, { key: 'df2', label: '자유도2', format: { type: 'int' } },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
    ],
    rows: [{ r: fit.r, r2: fit.r2, adj: fit.adjR2, s: fit.s, F: fit.F, df1: fit.dfm, df2: fit.dfr, p: fit.pF }],
    notes: ['독립변수 ' + modName_(ind[0]) + ', 조절변수 ' + modName_(W) + (covs.length ? ', 공변량 ' + covs.map(modName_).join(', ') : '') + '. 유효 사례 ' + n + '건.']
  });

  var centerNote = center
    ? '독립변수' + (categorical ? '' : '와 조절변수') + '를 평균으로 중심화했다. 그래서 "' + modName_(ind[0]) + '"의 계수는 ' + (categorical ? '기준 범주(' + input.w.labels[input.w.refIdx] + ')에서' : '"' + modName_(W) + '"가 평균일 때') + ' 독립변수의 효과이다.'
    : '중심화하지 않았다. "' + modName_(ind[0]) + '"의 계수는 조절변수가 0일 때(범주형이면 기준 범주일 때) 독립변수의 효과이다.';
  blocks.push({
    type: 'table',
    caption: '계수',
    columns: [
      { key: 'name', label: '' },
      { key: 'B', label: 'B', format: F3 }, { key: 'se', label: '표준오차', format: F3 },
      { key: 't', label: 't', format: F3 }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
      { key: 'lo', label: conf + '% 신뢰구간 하한', format: F3 }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: F3 }
    ],
    rows: fit.coef.map(function (c, idx) {
      var isInt = idx >= 1 + r.nMain && idx < 1 + r.nMain + r.q;
      return { _bold: isInt, name: c.name, B: c.B, se: c.se, t: c.t, p: c.p, lo: c.lo, hi: c.hi };
    }),
    notes: [centerNote, '굵게 표시한 행이 상호작용(조절효과) 항이다.']
  });

  blocks.push({
    type: 'table',
    caption: '상호작용 검정 (R² 변화량)',
    columns: [
      { key: 'dr2', label: 'R² 변화량(ΔR²)', format: { type: 'fixed', digits: 4 } },
      { key: 'F', label: 'F 변화량', format: F3 },
      { key: 'df1', label: '자유도1', format: { type: 'int' } }, { key: 'df2', label: '자유도2', format: { type: 'int' } },
      { key: 'p', label: '유의확률', format: { type: 'pvalue' } }
    ],
    rows: [{ dr2: r.deltaR2, F: r.F, df1: r.df1, df2: r.df2, p: r.pF }],
    notes: ['상호작용항을 넣기 전과 후의 R² 차이가 0인지 검정한다. 유의하면 조절효과가 있다 (범주형 조절변수는 모든 상호작용항을 한꺼번에 검정).']
  });

  var slopeCols = [{ key: 'label', label: '조절변수 수준' }];
  if (!categorical) slopeCols.push({ key: 'w', label: '조절변수 값', format: { type: 'trim', digits: 3 } });
  slopeCols = slopeCols.concat([
    { key: 'b', label: '독립변수의 효과(B)', format: F3 }, { key: 'se', label: '표준오차', format: F3 },
    { key: 't', label: 't', format: F3 }, { key: 'p', label: '유의확률', format: { type: 'pvalue' } },
    { key: 'lo', label: conf + '% 신뢰구간 하한', format: F3 }, { key: 'hi', label: conf + '% 신뢰구간 상한', format: F3 }
  ]);
  blocks.push({
    type: 'table',
    caption: '조건부 효과 (단순기울기): ' + modName_(ind[0]) + '의 효과',
    columns: slopeCols,
    rows: r.slopes,
    notes: [categorical
      ? '조절변수의 범주마다 독립변수의 효과(기울기)이다. 기준 범주의 효과는 B 그대로, 다른 범주는 B + 상호작용 계수이며 표준오차는 계수 공분산으로 구했다.'
      : '조절변수가 평균 − 1SD, 평균, 평균 + 1SD일 때 독립변수의 효과(기울기)이다. 표준오차는 √(Var(b1) + 2w·Cov(b1,b3) + w²·Var(b3)).']
  });

  if (r.jn) {
    var jnRows = r.jn.regions.map(function (g) {
      return {
        range: Number(g.lo.toFixed(3)) + ' ~ ' + Number(g.hi.toFixed(3)),
        result: g.sig ? '유의함 (' + (g.sign > 0 ? '양의 효과' : '음의 효과') + ')' : '유의하지 않음',
        pct: g.pct
      };
    });
    var jnNotes = ['조절변수의 관측 범위(' + Number(r.jn.min.toFixed(3)) + ' ~ ' + Number(r.jn.max.toFixed(3)) + ')를 독립변수의 효과가 유의한 구간과 아닌 구간으로 나눈 것이다 (유의수준 ' + Number(((100 - conf) / 100).toFixed(4)) + ').'];
    if (r.jn.roots.length) jnNotes.push('경계값: ' + r.jn.roots.map(function (v) { return Number(v.toFixed(3)); }).join(', '));
    else jnNotes.push('관측 범위 안에 경계값이 없어 전체 범위에서 결과가 같다.');
    if (opt.jn !== false) {
      blocks.push({
        type: 'table',
        caption: 'Johnson-Neyman: 조절변수 값에 따른 ' + modName_(ind[0]) + '의 효과의 유의성',
        columns: [
          { key: 'range', label: modName_(W) + ' 값의 구간' }, { key: 'result', label: '독립변수의 효과' },
          { key: 'pct', label: '사례 비율(%)', format: { type: 'fixed', digits: 1 } }
        ],
        rows: jnRows,
        notes: jnNotes
      });
    }
  }

  if (opt.plot !== false) {
    blocks.push(d.lineChartBlock('상호작용 그림: ' + modName_(dep[0]) + ' (' + modName_(W) + '별)', modName_(ind[0]), '예측된 ' + modName_(dep[0]), r.plot.series));
  }
  return { analysis: 'moderation', title: '조절효과 분석', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = {
    johnsonNeyman: johnsonNeyman,
    moderationFit: moderationFit,
    runModeration: runModeration
  };
}
