/**
 * IPA(중요도–수행도 분석, Importance-Performance Analysis; Martilla & James 1977).
 * 순수 함수: (dataset, request) → 결과 객체. 중요도 문항들(request.importance)과 수행도(만족도) 문항들(request.performance)을
 * 같은 개수, 같은 순서로 짝지어 문항마다 두 평균을 구하고, 십자 기준선으로 4사분면에 나눈다.
 *
 * - 문항 평균·표준편차: 문항마다 유효한 값으로 (결측은 문항별로 제외, 옵션으로 목록별 제거).
 * - 기준선(중요도 = 세로, 수행도 = 가로): 방식 'mean' 문항 평균들의 평균, 'median' 문항 평균들의 중앙값, 'custom' 직접 지정한 값.
 *   mean·median 은 여기에 이동량(± 조정)을 더한다. 직접 지정은 이동량을 쓰지 않는다.
 * - 사분면: 중요도 ≥ 기준선이면 "높음", 수행도 ≥ 기준선이면 "높음" (같으면 높음).
 *   집중 개선(중요도 높음·수행도 낮음), 현상 유지(둘 다 높음), 낮은 우선순위(둘 다 낮음), 과잉 투자(중요도 낮음·수행도 높음).
 *   어느 한 축에서 기준선과의 거리가 0.1 이하인 문항은 "경계"(기준선을 조금만 옮겨도 사분면이 바뀔 수 있음)로 표시한다.
 * - 격차 검정(옵션): 문항마다 짝이 모두 있는 사례의 (수행도 − 중요도) 대응표본 t검정, 효과크기 d_z, Bonferroni 보정 유의확률.
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function ipDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return { getVars: data.getVars, pValueT: dist.pValueT, ipaChartBlock: require('./Charts.js').ipaChartBlock };
  }
  return { getVars: getVars, pValueT: pValueT, ipaChartBlock: ipaChartBlock };
}

var IP_MAX_ITEMS_ = 30;
var IP_EDGE_ = 0.1;
var IP_QUADRANTS_ = ['집중 개선', '현상 유지', '낮은 우선순위', '과잉 투자'];
var IP_F2_ = { type: 'fixed', digits: 2 };
var IP_F3_ = { type: 'fixed', digits: 3 };
function ipName_(v) { return v.label ? v.label : v.name; }

function ipNum_(opt, key, label, def) {
  var v = opt[key];
  if (v === undefined || v === null || String(v).trim() === '') return def;
  var x = Number(v);
  if (!isFinite(x)) throw new Error('"' + label + '"은(는) 숫자여야 합니다.');
  return x;
}

/** 값 배열(null 제외)의 { n, mean, sd } */
function ipMoments_(col, use) {
  var n = 0, s = 0, i, v;
  for (i = 0; i < col.length; i++) { if (col[i] !== null && (!use || use[i])) { n++; s += col[i]; } }
  var mean = n ? s / n : null, ss = 0;
  for (i = 0; i < col.length; i++) { if (col[i] !== null && (!use || use[i])) { v = col[i] - mean; ss += v * v; } }
  return { n: n, mean: mean, sd: n > 1 ? Math.sqrt(ss / (n - 1)) : null };
}

function ipMedian_(a) {
  var s = a.slice().sort(function (x, y) { return x - y; }), n = s.length;
  return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2;
}

/** 사분면 번호: 0 집중 개선, 1 현상 유지, 2 낮은 우선순위, 3 과잉 투자 */
function ipaQuadrant(imp, perf, impLine, perfLine) {
  var hi = imp >= impLine - 1e-12, hp = perf >= perfLine - 1e-12;
  return hi ? (hp ? 1 : 0) : (hp ? 3 : 2);
}

function runIpa(dataset, request) {
  var d = ipDeps_();
  var impVars = d.getVars(dataset, request.importance, { numeric: true });
  var perfVars = d.getVars(dataset, request.performance, { numeric: true });
  var k = impVars.length;
  if (k < 3) throw new Error('IPA는 문항을 3개 이상 짝지어 선택하세요 (중요도 문항 ' + k + '개).');
  if (perfVars.length !== k) throw new Error('중요도 문항(' + k + '개)과 수행도 문항(' + perfVars.length + '개)의 개수가 같아야 합니다. 같은 순서로 짝지어 보내세요.');
  if (k > IP_MAX_ITEMS_) throw new Error('문항은 ' + IP_MAX_ITEMS_ + '쌍까지 선택할 수 있습니다.');
  var seen = {}, i, j;
  impVars.forEach(function (v) { seen[v.name] = true; });
  perfVars.forEach(function (v) { if (seen[v.name]) throw new Error('변수 "' + v.name + '"를 중요도와 수행도에 동시에 쓸 수 없습니다.'); });

  var opt = request.options || {};
  var baseline = opt.baseline === 'median' || opt.baseline === 'custom' ? opt.baseline : 'mean';
  var impShift = ipNum_(opt, 'impShift', '중요도 기준선 이동', 0), perfShift = ipNum_(opt, 'perfShift', '수행도 기준선 이동', 0);
  var impLine0 = ipNum_(opt, 'impLine', '중요도 기준선 값', null), perfLine0 = ipNum_(opt, 'perfLine', '수행도 기준선 값', null);
  if (baseline === 'custom' && (impLine0 === null || perfLine0 === null)) throw new Error('기준선 방식을 "직접 지정"으로 하면 중요도 기준선 값과 수행도 기준선 값을 모두 입력하세요.');
  var gapTest = opt.gap !== false;

  // 사례 선택 (목록별 제거 옵션)
  var use = null, dropped = 0;
  if (opt.listwise) {
    use = new Array(dataset.n);
    for (i = 0; i < dataset.n; i++) {
      var ok = true;
      for (j = 0; j < k; j++) if (dataset.data[impVars[j].name][i] === null || dataset.data[perfVars[j].name][i] === null) { ok = false; break; }
      use[i] = ok;
      if (!ok) dropped++;
    }
  }
  var used = dataset.n - dropped;
  var items = [], warnings = [];
  for (j = 0; j < k; j++) {
    var ci = dataset.data[impVars[j].name], cp = dataset.data[perfVars[j].name];
    var mi = ipMoments_(ci, use), mp = ipMoments_(cp, use);
    if (mi.n < 2 || mp.n < 2) throw new Error('문항 "' + ipName_(impVars[j]) + '" / "' + ipName_(perfVars[j]) + '"에 값이 2개 이상 있어야 합니다.');
    var np = 0, sd = 0, sdd = 0, diffs = [];
    for (i = 0; i < dataset.n; i++) {
      if ((use && !use[i]) || ci[i] === null || cp[i] === null) continue;
      diffs.push(cp[i] - ci[i]); sd += cp[i] - ci[i];
    }
    np = diffs.length;
    var mdiff = np ? sd / np : null, sdiff = null;
    if (np > 1) { for (i = 0; i < np; i++) sdd += (diffs[i] - mdiff) * (diffs[i] - mdiff); sdiff = Math.sqrt(sdd / (np - 1)); }
    var test = { n: np, mdiff: mdiff, sd: sdiff, t: null, df: null, p: null, dz: null };
    if (np > 2 && sdiff > 0) {
      test.t = mdiff / (sdiff / Math.sqrt(np)); test.df = np - 1; test.p = d.pValueT(test.t, test.df); test.dz = mdiff / sdiff;
    }
    items.push({ name: ipName_(impVars[j]), perfName: ipName_(perfVars[j]), imp: mi, perf: mp, gap: mp.mean - mi.mean, test: test });
  }

  // 기준선
  var impMeans = items.map(function (x) { return x.imp.mean; }), perfMeans = items.map(function (x) { return x.perf.mean; });
  var avg = function (a) { var s = 0; a.forEach(function (v) { s += v; }); return s / a.length; };
  var baseImp, basePerf;
  if (baseline === 'mean') { baseImp = avg(impMeans); basePerf = avg(perfMeans); }
  else if (baseline === 'median') { baseImp = ipMedian_(impMeans); basePerf = ipMedian_(perfMeans); }
  else { baseImp = impLine0; basePerf = perfLine0; }
  var shiftI = baseline === 'custom' ? 0 : impShift, shiftP = baseline === 'custom' ? 0 : perfShift;
  var lineImp = baseImp + shiftI, linePerf = basePerf + shiftP;
  var baseLabel = baseline === 'mean' ? '문항 평균들의 평균' : baseline === 'median' ? '문항 평균들의 중앙값' : '직접 지정';

  var counts = [0, 0, 0, 0], byQuad = [[], [], [], []], edgeCount = 0;
  items.forEach(function (x) {
    x.q = ipaQuadrant(x.imp.mean, x.perf.mean, lineImp, linePerf);
    x.edge = Math.abs(x.imp.mean - lineImp) <= IP_EDGE_ + 1e-12 || Math.abs(x.perf.mean - linePerf) <= IP_EDGE_ + 1e-12;
    if (x.edge) edgeCount++;
    counts[x.q]++; byQuad[x.q].push(x.name);
  });

  // 경고
  var maxI = -Infinity, maxP = -Infinity;
  impVars.forEach(function (v) { dataset.data[v.name].forEach(function (x) { if (x !== null && x > maxI) maxI = x; }); });
  perfVars.forEach(function (v) { dataset.data[v.name].forEach(function (x) { if (x !== null && x > maxP) maxP = x; }); });
  if (Math.abs(maxI - maxP) >= 2) warnings.push('중요도 문항의 최댓값(' + maxI + ')과 수행도 문항의 최댓값(' + maxP + ')이 크게 달라 척도 범위가 다른 것 같습니다. 같은 척도(예: 1~5점)가 아니면 두 평균을 같은 그림에서 비교하기 어렵습니다.');
  if (dropped > 0) warnings.push('목록별 제거로 ' + dropped + '건을 제외했습니다 (선택한 모든 문항에 값이 있는 ' + used + '건 사용).');
  if (edgeCount > 0) warnings.push('기준선에서 ' + IP_EDGE_ + '점 이내인 "경계" 문항이 ' + edgeCount + '개 있습니다. 기준선을 조금만 옮겨도 사분면이 바뀔 수 있으니 해석에 주의하세요.');
  if (used < 30) warnings.push('사례가 ' + used + '건으로 적어 문항 평균이 불안정할 수 있습니다.');

  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: '값', format: { type: 'int' } }],
    rows: [{ label: '데이터의 사례 수', n: dataset.n }, { label: '분석에 쓴 사례 수 (목록별 제거일 때만 줄어듦)', n: used }, { label: '문항 쌍 (중요도 ↔ 수행도)', n: k }],
    notes: ['중요도 문항과 수행도 문항을 보낸 순서대로 짝지었다 (' + items.map(function (x) { return x.name + ' ↔ ' + x.perfName; }).slice(0, 6).join(', ') + (k > 6 ? ', …' : '') + '). 문항 평균은 문항마다 값이 있는 사례로 구했다.']
  });
  blocks.push({
    type: 'table', caption: '문항별 중요도·수행도와 사분면',
    columns: [
      { key: 'name', label: '문항 (중요도 변수)' }, { key: 'perfName', label: '수행도 변수' },
      { key: 'ni', label: 'N(중요도)', format: { type: 'int' } }, { key: 'mi', label: '중요도 평균', format: IP_F2_ }, { key: 'si', label: '표준편차', format: IP_F2_ },
      { key: 'np', label: 'N(수행도)', format: { type: 'int' } }, { key: 'mp', label: '수행도 평균', format: IP_F2_ }, { key: 'sp', label: '표준편차', format: IP_F2_ },
      { key: 'gap', label: '격차(수행도−중요도)', format: IP_F2_ },
      { key: 'di', label: '중요도−기준선', format: IP_F2_ }, { key: 'dp', label: '수행도−기준선', format: IP_F2_ },
      { key: 'quad', label: '사분면' }, { key: 'edge', label: '경계' }
    ],
    rows: items.map(function (x) {
      return { _bold: x.q === 0, name: x.name, perfName: x.perfName, ni: x.imp.n, mi: x.imp.mean, si: x.imp.sd, np: x.perf.n, mp: x.perf.mean, sp: x.perf.sd,
        gap: x.gap, di: x.imp.mean - lineImp, dp: x.perf.mean - linePerf, quad: IP_QUADRANTS_[x.q], edge: x.edge ? '★' : '' };
    }),
    notes: ['사분면은 중요도 평균과 수행도 평균을 기준선과 비교해 정한다 (같으면 높음으로 본다). 집중 개선 문항은 굵게 표시했다. ★ 경계는 기준선과의 거리가 한 축에서 ' + IP_EDGE_ + '점 이내인 문항이다.']
  });
  blocks.push({
    type: 'table', caption: '기준선',
    columns: [
      { key: 'axis', label: '축' }, { key: 'how', label: '방식' }, { key: 'base', label: '계산·지정값', format: IP_F3_ },
      { key: 'shift', label: '이동량', format: IP_F3_ }, { key: 'line', label: '최종 기준선', format: IP_F3_ }
    ],
    rows: [
      { axis: '중요도 (세로)', how: baseLabel, base: baseImp, shift: shiftI, line: lineImp },
      { axis: '수행도 (가로)', how: baseLabel, base: basePerf, shift: shiftP, line: linePerf }
    ],
    notes: ['기준선을 어디에 두느냐에 따라 사분면이 달라진다. 보고서에는 어떤 기준선을 썼는지(방식과 이동량)를 함께 적으세요.' + (baseline === 'custom' ? ' 직접 지정한 값은 이동량을 쓰지 않는다.' : ' 이동량은 계산된 값에 더한 조정값이다(옵션).')]
  });
  blocks.push({
    type: 'table', caption: '사분면별 문항',
    columns: [{ key: 'quad', label: '사분면' }, { key: 'desc', label: '의미' }, { key: 'n', label: '문항 수', format: { type: 'int' } }, { key: 'items', label: '문항' }],
    rows: [
      { quad: '집중 개선', desc: '중요도 높음, 수행도 낮음', n: counts[0], items: byQuad[0].join(', ') },
      { quad: '현상 유지', desc: '중요도 높음, 수행도 높음', n: counts[1], items: byQuad[1].join(', ') },
      { quad: '낮은 우선순위', desc: '중요도 낮음, 수행도 낮음', n: counts[2], items: byQuad[2].join(', ') },
      { quad: '과잉 투자', desc: '중요도 낮음, 수행도 높음', n: counts[3], items: byQuad[3].join(', ') }
    ],
    notes: ['일반적으로 집중 개선 문항에 자원을 먼저 투입하고, 현상 유지 문항은 지금 수준을 지키며, 과잉 투자 문항은 자원을 줄일 여지가 있는지 본다.']
  });
  if (gapTest) {
    var m = items.length;
    blocks.push({
      type: 'table', caption: '격차 검정 (대응표본 t검정, 수행도 − 중요도)',
      columns: [
        { key: 'name', label: '문항' }, { key: 'n', label: 'N(짝)', format: { type: 'int' } }, { key: 'md', label: '평균 차이', format: IP_F3_ },
        { key: 'sd', label: '표준편차', format: IP_F3_ }, { key: 't', label: 't', format: IP_F3_ }, { key: 'df', label: '자유도', format: { type: 'int' } },
        { key: 'p', label: '유의확률(양측)', format: { type: 'pvalue' } }, { key: 'pa', label: '보정 유의확률 (×' + m + ')', format: { type: 'pvalue' } },
        { key: 'dz', label: '효과크기 d_z', format: IP_F3_ }
      ],
      rows: items.map(function (x) {
        return { name: x.name, n: x.test.n, md: x.test.mdiff, sd: x.test.sd, t: x.test.t, df: x.test.df, p: x.test.p, pa: x.test.p === null ? null : Math.min(1, x.test.p * m), dz: x.test.dz };
      }),
      notes: ['음수 평균 차이는 수행도가 중요도보다 낮다는 뜻이다. 두 값이 모두 있는 사례만 쓴다(문항 평균의 차이와 짝지은 차이의 평균은 결측이 있으면 조금 다를 수 있다). 문항이 여러 개라 Bonferroni로 보정한 유의확률을 함께 보였다.']
    });
  }
  if (opt.plot !== false) {
    blocks.push(d.ipaChartBlock('IPA 사분면 그림', items.map(function (x) { return { label: x.name, x: x.perf.mean, y: x.imp.mean, q: x.q }; }), linePerf, lineImp, '수행도(만족도) 평균', '중요도 평균'));
  }
  return { analysis: 'ipa', title: 'IPA 중요도–수행도 분석', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = { ipaQuadrant: ipaQuadrant, runIpa: runIpa };
}
