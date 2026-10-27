/**
 * 신뢰도 분석: 크론바흐 α. 순수 함수: (dataset, request) → 결과 객체.
 * 결과표 구성: 신뢰도 분석 보고 형식(사례 처리 요약, 신뢰도 통계량, 문항 통계량, 문항-전체 통계량, 문항 간 상관, 척도 통계량).
 *
 * - 결측: 목록별 제거 (선택한 모든 문항이 유효한 사례만 사용)
 * - α = k/(k−1) · (1 − Σ분산 / 합계 분산), 분산은 표본 분산(n−1). 이분 문항이면 KR-20과 같다.
 * - 표준화 α = k·r̄ / (1 + (k−1)·r̄), r̄ 는 문항 간 상관의 평균
 * - 신뢰구간: Feldt(1965). (1−α̂)/(1−α) 가 F(n−1, (n−1)(k−1))을 따른다는 성질
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function alphaDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return { getVars: data.getVars, qf: dist.qf };
  }
  return { getVars: getVars, qf: qf };
}

var ALPHA_MAX_ITEMS_ = 100;
function rlbName_(v) { return v.label ? v.label : v.name; }

/**
 * columns: 문항별 값 배열 (결측은 null, 길이는 모두 같음). 목록별 결측 제거 후 계산한다.
 * 반환: { k, n, nAll, alpha, alphaStd, scaleMean, scaleVar, items:[{mean, sd, n, scaleMeanIfDeleted, scaleVarIfDeleted, corrected, alphaIfDeleted}],
 *         r(문항 간 상관 행렬 또는 null), interitem({mean, min, max} 또는 null), ci({conf, lo, hi} 또는 null) }
 * 계산할 수 없는 값은 null (사례 2개 미만, 합계 분산 0, 분산이 0인 문항의 상관 등).
 */
function cronbachAlpha(columns, conf) {
  var d = alphaDeps_();
  var k = columns.length, nAll = k ? columns[0].length : 0, i, j, r;
  var keep = [];
  for (r = 0; r < nAll; r++) {
    var ok = true;
    for (i = 0; i < k; i++) if (columns[i][r] === null) { ok = false; break; }
    if (ok) keep.push(r);
  }
  var n = keep.length;
  var res = { k: k, n: n, nAll: nAll, alpha: null, alphaStd: null, scaleMean: null, scaleVar: null, items: [], r: null, interitem: null, ci: null };
  for (i = 0; i < k; i++) res.items.push({ mean: null, sd: null, n: n, scaleMeanIfDeleted: null, scaleVarIfDeleted: null, corrected: null, alphaIfDeleted: null });
  if (n < 2 || k < 1) return res;

  // 평균과 편차 (두 번째 패스로 평균의 반올림 오차 보정)
  var xs = [], means = [], devs = [], devSums = [];
  for (i = 0; i < k; i++) {
    var x = keep.map(function (row) { return columns[i][row]; });
    var sum = 0, t;
    for (t = 0; t < n; t++) sum += x[t];
    var m = sum / n, corr = 0;
    for (t = 0; t < n; t++) corr += x[t] - m;
    m += corr / n;
    var dv = new Array(n), ds = 0;
    for (t = 0; t < n; t++) { dv[t] = x[t] - m; ds += dv[t]; }
    xs.push(x); means.push(m); devs.push(dv); devSums.push(ds);
  }
  // 공분산 행렬 (표본, n−1)
  var C = [];
  for (i = 0; i < k; i++) {
    C.push(new Array(k));
  }
  for (i = 0; i < k; i++) {
    for (j = i; j < k; j++) {
      var s = 0;
      for (var q = 0; q < n; q++) s += devs[i][q] * devs[j][q];
      C[i][j] = C[j][i] = (s - devSums[i] * devSums[j] / n) / (n - 1);
    }
  }
  var trace = 0, S = 0;
  for (i = 0; i < k; i++) { trace += C[i][i]; for (j = 0; j < k; j++) S += C[i][j]; }
  var meanSum = 0;
  for (i = 0; i < k; i++) meanSum += means[i];
  res.scaleMean = meanSum;
  res.scaleVar = S;
  if (k >= 2 && S > 0) res.alpha = k / (k - 1) * (1 - trace / S);

  for (i = 0; i < k; i++) {
    var it = res.items[i], rowSum = 0;
    for (j = 0; j < k; j++) rowSum += C[i][j];
    var varRest = S - 2 * rowSum + C[i][i];                 // 문항 i 를 뺀 합계의 분산
    it.mean = means[i];
    it.sd = Math.sqrt(Math.max(0, C[i][i]));
    it.scaleMeanIfDeleted = meanSum - means[i];
    it.scaleVarIfDeleted = Math.max(0, varRest);
    if (C[i][i] > 0 && varRest > 0) it.corrected = Math.max(-1, Math.min(1, (rowSum - C[i][i]) / Math.sqrt(C[i][i] * varRest)));
    if (k >= 3 && varRest > 0) it.alphaIfDeleted = (k - 1) / (k - 2) * (1 - (trace - C[i][i]) / varRest);
  }

  // 문항 간 상관 (모든 문항의 분산이 0보다 클 때만)
  var allVar = true;
  for (i = 0; i < k; i++) if (!(C[i][i] > 0)) allVar = false;
  if (allVar && k >= 2) {
    var R = [], sum2 = 0, cnt = 0, mn = Infinity, mx = -Infinity;
    for (i = 0; i < k; i++) {
      R.push(new Array(k));
      for (j = 0; j < k; j++) {
        R[i][j] = i === j ? 1 : Math.max(-1, Math.min(1, C[i][j] / Math.sqrt(C[i][i] * C[j][j])));
        if (i < j) { sum2 += R[i][j]; cnt++; if (R[i][j] < mn) mn = R[i][j]; if (R[i][j] > mx) mx = R[i][j]; }
      }
    }
    var rbar = sum2 / cnt;
    res.r = R;
    res.interitem = { mean: rbar, min: mn, max: mx };
    if (1 + (k - 1) * rbar !== 0) res.alphaStd = k * rbar / (1 + (k - 1) * rbar);
  }

  // Feldt 신뢰구간
  var confLevel = conf === undefined || conf === null ? 95 : conf;
  if (res.alpha !== null && n >= 3 && k >= 2) {
    var g = 1 - confLevel / 100, df1 = n - 1, df2 = (n - 1) * (k - 1);
    var fLo = d.qf(g / 2, df1, df2, true), fHi = d.qf(g / 2, df1, df2, false);
    res.ci = { conf: confLevel, lo: 1 - (1 - res.alpha) / fLo, hi: 1 - (1 - res.alpha) / fHi };
  }
  return res;
}

function runReliability(dataset, request) {
  var d = alphaDeps_();
  var vars = d.getVars(dataset, request.vars, { numeric: true });
  if (vars.length < 2) throw new Error('신뢰도 분석에는 문항(변수)을 2개 이상 선택하세요.');
  if (vars.length > ALPHA_MAX_ITEMS_) throw new Error('문항은 ' + ALPHA_MAX_ITEMS_ + '개까지 선택할 수 있습니다.');
  var opt = request.options || {};
  var conf = opt.conf === undefined || opt.conf === null || opt.conf === '' ? 95 : Number(opt.conf);
  if (!isFinite(conf) || conf <= 0 || conf >= 100) throw new Error('"신뢰수준(%)"은(는) 0 초과 100 미만이어야 합니다.');

  var columns = vars.map(function (v) { return dataset.data[v.name]; });
  var r = cronbachAlpha(columns, conf);
  var names = vars.map(rlbName_);
  var warnings = [];
  var F3 = { type: 'fixed', digits: 3 };

  if (r.n < 2) throw new Error('모든 문항에 응답한 사례가 ' + r.n + '건이라 계산할 수 없습니다 (2건 이상 필요). 결측이 많은 문항을 빼 보세요.');
  if (r.nAll - r.n > 0) warnings.push('목록별 결측 제거로 ' + (r.nAll - r.n) + '건을 제외했습니다 (선택한 문항 중 하나라도 빠진 사례). 사용한 사례는 ' + r.n + '건입니다.');
  if (r.alpha === null) warnings.push('문항 합계의 분산이 0이라 크론바흐 α를 계산할 수 없습니다.');
  r.items.forEach(function (it, i) {
    if (it.sd === 0) warnings.push('"' + names[i] + '": 모든 값이 같아 문항-전체 상관을 계산할 수 없습니다. 이 문항을 빼는 것을 고려하세요.');
  });
  if (r.alpha !== null && r.alpha < 0) {
    warnings.push('크론바흐 α가 음수입니다. 문항 사이의 평균 공분산이 음수라는 뜻으로, 역문항(부정 문항)을 뒤집지 않았거나 문항이 같은 특성을 재지 않을 수 있습니다. 메뉴 Stat4SW 통계분석 > 변수 계산 > 역문항 처리로 역문항을 뒤집은 새 열을 만들어 다시 분석하세요.');
  } else {
    var neg = [];
    r.items.forEach(function (it, i) { if (it.corrected !== null && it.corrected < 0) neg.push('"' + names[i] + '"'); });
    if (neg.length) warnings.push('수정된 문항-전체 상관이 음수인 문항이 있습니다: ' + neg.join(', ') + '. 역문항이라면 먼저 역문항 처리(메뉴 Stat4SW 통계분석 > 변수 계산)로 뒤집은 열을 쓰세요.');
  }
  if (r.n < 30) warnings.push('사례가 ' + r.n + '건으로 적어 α의 신뢰구간이 매우 넓을 수 있습니다.');

  var blocks = [];
  blocks.push({
    type: 'table',
    caption: '사례 처리 요약',
    columns: [
      { key: 'label', label: '' },
      { key: 'n', label: 'N', format: { type: 'int' } },
      { key: 'pct', label: '%', format: { type: 'fixed', digits: 1 } }
    ],
    rows: [
      { label: '유효', n: r.n, pct: r.nAll ? r.n / r.nAll * 100 : null },
      { label: '제외됨', n: r.nAll - r.n, pct: r.nAll ? (r.nAll - r.n) / r.nAll * 100 : null },
      { _kind: 'total', label: '합계', n: r.nAll, pct: r.nAll ? 100 : null }
    ],
    notes: ['목록별 결측 제거: 선택한 모든 문항에 응답한 사례만 사용했다.']
  });

  var statCols = [
    { key: 'alpha', label: '크론바흐 α', format: F3 },
    { key: 'alphaStd', label: '표준화된 문항에 기초한 α', format: F3 },
    { key: 'k', label: '문항 수', format: { type: 'int' } }
  ];
  var statRow = { alpha: r.alpha, alphaStd: r.alphaStd, k: r.k };
  if (opt.ci) {
    statCols.splice(2, 0, { key: 'lo', label: 'α의 ' + conf + '% 신뢰구간 하한', format: F3 }, { key: 'hi', label: 'α의 ' + conf + '% 신뢰구간 상한', format: F3 });
    statRow.lo = r.ci ? r.ci.lo : null;
    statRow.hi = r.ci ? r.ci.hi : null;
  }
  var statNotes = [
    'α = k/(k−1) · (1 − 문항 분산의 합 / 합계 분산), 표본 분산(n−1). 문항이 모두 이분(0/1)이면 KR-20과 같다.',
    '표준화 α는 문항 간 상관의 평균(r̄)으로 k·r̄ / (1 + (k−1)·r̄). 문항의 척도(분산)가 크게 다를 때 참고한다.',
    '흔히 쓰는 기준은 .70 이상이지만 문항 수가 많으면 α가 커지고, 목적(탐색, 연구, 개인 판정)에 따라 요구 수준이 다르다. 높을수록 좋은 것만은 아니며 .95를 넘으면 문항이 중복일 수 있다.'
  ];
  if (opt.ci) statNotes.push('신뢰구간은 Feldt(1965) 방법이다 ((1−α̂)/(1−α) ~ F(n−1, (n−1)(k−1))).');
  blocks.push({ type: 'table', caption: '신뢰도 통계량', columns: statCols, rows: [statRow], notes: statNotes });

  if (opt.itemDesc) {
    blocks.push({
      type: 'table',
      caption: '문항 통계량',
      columns: [
        { key: 'name', label: '문항' },
        { key: 'mean', label: '평균', format: F3 },
        { key: 'sd', label: '표준편차', format: F3 },
        { key: 'n', label: 'N', format: { type: 'int' } }
      ],
      rows: r.items.map(function (it, i) { return { name: names[i], mean: it.mean, sd: it.sd, n: it.n }; }),
      notes: []
    });
  }

  if (opt.itemTotal) {
    var better = false;
    var rows = r.items.map(function (it, i) {
      var up = r.alpha !== null && it.alphaIfDeleted !== null && it.alphaIfDeleted > r.alpha + 1e-12;
      if (up) better = true;
      return {
        _bold: up, name: names[i], meanDel: it.scaleMeanIfDeleted, varDel: it.scaleVarIfDeleted,
        corrected: it.corrected, alphaDel: it.alphaIfDeleted
      };
    });
    var itNotes = [
      '수정된 문항-전체 상관은 그 문항과 나머지 문항 합계의 상관이다. 대략 .30 미만이면 척도와 잘 맞지 않는 문항으로 본다.',
      '"문항 제외 시 α"는 그 문항을 뺐을 때의 크론바흐 α이다 (문항이 3개 이상일 때 계산한다).'
    ];
    if (better) itNotes.push('굵게 표시한 문항은 제외하면 α가 올라간다. 내용상 문제가 없다면 빼는 것을 고려하되, 문항을 빼기 전에 의미(내용 타당도)를 먼저 검토한다.');
    blocks.push({
      type: 'table',
      caption: '문항-전체 통계량',
      columns: [
        { key: 'name', label: '문항' },
        { key: 'meanDel', label: '문항 제외 시 척도 평균', format: F3 },
        { key: 'varDel', label: '문항 제외 시 척도 분산', format: F3 },
        { key: 'corrected', label: '수정된 문항-전체 상관', format: F3 },
        { key: 'alphaDel', label: '문항 제외 시 α', format: F3 }
      ],
      rows: rows,
      notes: itNotes
    });
  }

  if (opt.interitem) {
    if (r.r) {
      var cols = [{ key: 'name', label: '' }];
      names.forEach(function (nm, idx) { cols.push({ key: 'c' + idx, label: nm, format: F3 }); });
      blocks.push({
        type: 'table',
        caption: '문항 간 상관행렬',
        columns: cols,
        rows: r.r.map(function (rowVals, i) {
          var row = { name: names[i] };
          rowVals.forEach(function (v, j) { row['c' + j] = v; });
          return row;
        }),
        notes: ['피어슨 상관계수. 사용한 사례만 계산했다.']
      });
      blocks.push({
        type: 'table',
        caption: '문항 간 상관 요약',
        columns: [
          { key: 'mean', label: '평균', format: F3 },
          { key: 'min', label: '최솟값', format: F3 },
          { key: 'max', label: '최댓값', format: F3 }
        ],
        rows: [{ mean: r.interitem.mean, min: r.interitem.min, max: r.interitem.max }],
        notes: ['평균 문항 간 상관은 대략 .15~.50이 바람직하다고 보는 경우가 많다.']
      });
    } else {
      warnings.push('분산이 0인 문항이 있어 문항 간 상관행렬을 만들 수 없습니다.');
    }
  }

  if (opt.scale) {
    blocks.push({
      type: 'table',
      caption: '척도 통계량 (문항 합계)',
      columns: [
        { key: 'mean', label: '평균', format: F3 },
        { key: 'var', label: '분산', format: F3 },
        { key: 'sd', label: '표준편차', format: F3 },
        { key: 'k', label: '문항 수', format: { type: 'int' } }
      ],
      rows: [{ mean: r.scaleMean, var: r.scaleVar, sd: r.scaleVar === null ? null : Math.sqrt(Math.max(0, r.scaleVar)), k: r.k }],
      notes: ['선택한 문항의 값을 모두 더한 점수의 통계량이다.']
    });
  }

  return { analysis: 'reliability', title: '신뢰도 분석 (크론바흐 α)', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = {
    cronbachAlpha: cronbachAlpha,
    runReliability: runReliability
  };
}
