/**
 * Cochran Q 검정: 같은 사람을 3회 이상 측정한 예/아니오(0/1) 반응의 성공 비율이 측정마다 달라졌는지 검정한다.
 * McNemar 검정(2회)의 3회 이상 확장이다. 순수 함수: (dataset, request) → 결과 객체.
 * 측정마다 변수 한 열, 사람이 행이다. 하나라도 비어 있는 사례는 통째로 제외한다 (목록별 제거).
 *
 * - 성공: 변수 값이 옵션 "성공으로 볼 값"과 같은 사례, 나머지는 모두 실패 (0/1 로 코딩된 자료는 값 1이 성공).
 * - Q = (k − 1)[k Σ C_j² − T²] / (k T − Σ R_i²), C_j = 측정 j 의 성공 수, R_i = 사람 i 의 성공 수, T = 전체 성공 수.
 *   H0 아래 Q ~ χ²(k − 1). 모든 측정에서 응답이 같은 사례(R_i = 0 또는 k)는 Q 에 기여하지 않는다.
 *   0/1 자료에 Friedman 검정을 하면 Q 와 같은 값이 나온다.
 * - 효과크기 Kendall의 W = Q / (N (k − 1)).
 * - 사후비교: 측정 쌍마다 정확 McNemar 검정(이항, 불일치 쌍만)이고 유의확률에 비교 횟수를 곱한다 (Bonferroni, 최대 1).
 */

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function cqDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    var dist = require('./Distributions.js');
    return { getVars: data.getVars, pchisq: dist.pchisq, mcnemarExactP: require('./Stats_Mcnemar.js').mcnemarExactP, lineChartBlock: require('./Charts.js').lineChartBlock };
  }
  return { getVars: getVars, pchisq: pchisq, mcnemarExactP: mcnemarExactP, lineChartBlock: lineChartBlock };
}

var CQ_MAX_VARS_ = 20;
var CQ_F3_ = { type: 'fixed', digits: 3 };
function cqName_(v) { return v.label ? v.label : v.name; }

/**
 * x: 사례 × 측정의 0/1 행렬. 반환: { n, k, q, df, p, w, colsum, informative } (분모가 0 이면 q, p, w 는 null)
 */
function cochranQ(x, pchisq) {
  var n = x.length, k = x[0].length, i, j;
  var col = new Array(k).fill(0), T = 0, sumR2 = 0, informative = 0;
  for (i = 0; i < n; i++) {
    var r = 0;
    for (j = 0; j < k; j++) { r += x[i][j]; col[j] += x[i][j]; }
    T += r; sumR2 += r * r;
    if (r > 0 && r < k) informative++;
  }
  var sumC2 = 0;
  for (j = 0; j < k; j++) sumC2 += col[j] * col[j];
  var den = k * T - sumR2;
  var res = { n: n, k: k, q: null, df: k - 1, p: null, w: null, colsum: col, informative: informative };
  if (!(den > 0)) return res;
  res.q = (k - 1) * (k * sumC2 - T * T) / den;
  res.p = pchisq(res.q, k - 1, false);
  res.w = res.q / (n * (k - 1));
  return res;
}

function runCochran(dataset, request) {
  var d = cqDeps_();
  var vars = d.getVars(dataset, request.vars);
  if (vars.length < 3) throw new Error('Cochran Q 검정은 측정 변수를 3개 이상 선택하세요 (2번 측정이면 McNemar 검정을 쓰세요).');
  if (vars.length > CQ_MAX_VARS_) throw new Error('측정 변수는 ' + CQ_MAX_VARS_ + '개까지 선택할 수 있습니다.');
  var opt = request.options || {};
  var succ = opt.success === undefined || opt.success === null || String(opt.success).trim() === '' ? '1' : String(opt.success).trim();
  var num = Number(succ);
  var isSucc = function (v) { return typeof v === 'number' && isFinite(num) ? v === num : String(v).trim() === succ; };
  var k = vars.length, i, j;
  var rowsOk = [];
  for (i = 0; i < dataset.n; i++) {
    var ok = true;
    for (j = 0; j < k; j++) if (dataset.data[vars[j].name][i] === null) { ok = false; break; }
    if (ok) rowsOk.push(i);
  }
  var N = rowsOk.length;
  if (N < 2) throw new Error('모든 측정 변수에 값이 있는 사례가 ' + N + '건입니다 (2건 이상 필요).');
  var seen = {}, order = [];
  var x = rowsOk.map(function (r) {
    return vars.map(function (v) {
      var val = dataset.data[v.name][r], key = String(val);
      if (!(key in seen)) { seen[key] = 0; order.push(val); }
      seen[key]++;
      return isSucc(val) ? 1 : 0;
    });
  });
  var totalSucc = 0;
  x.forEach(function (row) { row.forEach(function (v) { totalSucc += v; }); });
  if (totalSucc === 0) throw new Error('성공으로 볼 값 "' + succ + '"이(가) 측정 변수에 없습니다. 변수에 있는 값: ' + order.slice(0, 10).map(String).join(', ') + (order.length > 10 ? ', …' : '') + '. 성공으로 볼 값을 다시 입력하세요.');
  var warnings = [];
  if (order.length > 2) warnings.push('측정 변수에 값이 ' + order.length + '가지 있습니다. "' + succ + '"인 응답만 성공이고 나머지는 모두 실패로 처리했습니다.');
  if (dataset.n > N) warnings.push('목록별 결측 제거로 ' + (dataset.n - N) + '건을 제외했습니다 (측정 변수 중 하나라도 비어 있는 사람). 사용한 사례는 ' + N + '건입니다.');
  var res = cochranQ(x, d.pchisq);
  if (res.q === null) warnings.push('모든 사례가 모든 측정에서 같은 응답이라 Q 를 계산할 수 없습니다.');
  else if (res.informative < 10) warnings.push('측정마다 응답이 달라진 사례가 ' + res.informative + '건뿐이라 검정의 근사(카이제곱)가 부정확할 수 있습니다.');
  var names = vars.map(cqName_);

  var blocks = [];
  blocks.push({
    type: 'table', caption: '처리 요약',
    columns: [{ key: 'label', label: '' }, { key: 'n', label: 'N', format: { type: 'int' } }],
    rows: [
      { label: '분석에 포함 (사람)', n: N }, { label: '결측으로 제외', n: dataset.n - N }, { _kind: 'total', label: '합계', n: dataset.n },
      { label: '측정 횟수 (변수)', n: k }, { label: '측정마다 응답이 달라진 사례', n: res.informative }
    ],
    notes: ['측정 변수: ' + names.join(', ') + '. 성공으로 본 값은 "' + succ + '"이고 나머지는 실패이다. 모든 측정에서 응답이 같은 사례는 Q 에 영향을 주지 않는다.']
  });
  blocks.push({
    type: 'table', caption: '측정별 성공 빈도',
    columns: [{ key: 'name', label: '측정(변수)' }, { key: 'x', label: '성공', format: { type: 'int' } }, { key: 'y', label: '실패', format: { type: 'int' } }, { key: 'pct', label: '성공 비율', format: CQ_F3_ }],
    rows: names.map(function (nm, q) { return { name: nm, x: res.colsum[q], y: N - res.colsum[q], pct: res.colsum[q] / N }; })
  });
  blocks.push({
    type: 'table', caption: 'Cochran Q 검정',
    columns: [
      { key: 'n', label: 'N', format: { type: 'int' } }, { key: 'q', label: 'Cochran Q', format: CQ_F3_ }, { key: 'df', label: '자유도', format: { type: 'int' } },
      { key: 'p', label: '점근 유의확률', format: { type: 'pvalue' } }, { key: 'w', label: "Kendall의 W", format: CQ_F3_ }
    ],
    rows: [{ n: N, q: res.q, df: res.df, p: res.p, w: res.w }],
    notes: ['귀무가설: 성공 비율이 모든 측정에서 같다. 유의하면 어느 측정 사이가 다른지 아래 쌍별 비교로 본다. Kendall의 W = Q / (N(k−1))는 0~1의 효과크기이다 (같은 자료의 Friedman 검정과 Q 가 같다).']
  });
  if (opt.posthoc !== false && res.q !== null) {
    var pairRows = [], m = k * (k - 1) / 2;
    for (i = 0; i < k - 1; i++) {
      for (j = i + 1; j < k; j++) {
        var b = 0, c = 0;
        x.forEach(function (row) { if (row[i] === 1 && row[j] === 0) b++; else if (row[i] === 0 && row[j] === 1) c++; });
        var pExact = d.mcnemarExactP(b, c);
        pairRows.push({ a: names[i], b: names[j], bc: b, cb: c, p: pExact, pAdj: Math.min(1, pExact * m) });
      }
    }
    blocks.push({
      type: 'table', caption: '쌍별 비교 (정확 McNemar, Bonferroni 보정)',
      columns: [
        { key: 'a', label: '측정 A' }, { key: 'b', label: '측정 B' }, { key: 'bc', label: 'A 성공·B 실패', format: { type: 'int' } }, { key: 'cb', label: 'A 실패·B 성공', format: { type: 'int' } },
        { key: 'p', label: '정확 유의확률', format: { type: 'pvalue' } }, { key: 'pAdj', label: '보정 유의확률 (×' + m + ')', format: { type: 'pvalue' } }
      ],
      rows: pairRows,
      notes: ['두 측정 사이에서 응답이 바뀐 사람(불일치 쌍)만 쓰는 정확 이항 McNemar 검정이고, 유의확률에 비교 횟수(' + m + ')를 곱했다 (1을 넘으면 1).']
    });
  }
  if (opt.plot !== false) {
    blocks.push(d.lineChartBlock('측정별 성공 비율', '측정', '성공 비율', [{ label: '성공 비율', points: names.map(function (nm, q) { return [q + 1, res.colsum[q] / N]; }) }], names));
  }
  return { analysis: 'cochran', title: 'Cochran Q 검정 (반복측정 예/아니오)', blocks: blocks, warnings: warnings };
}

if (typeof module !== 'undefined') {
  module.exports = { cochranQ: cochranQ, runCochran: runCochran };
}
