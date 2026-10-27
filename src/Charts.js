/**
 * 차트: 막대그래프, 히스토그램, 산점도, 상자그림, 선 그래프, 경로도, 나무 그림, IPA 사분면 그림. 외부 라이브러리 없이 SVG를 직접 만든다.
 * 결과 객체의 차트 블록 { type: 'chart', kind, caption, ... } 을 만드는 순수 함수(*Block)와
 * 그것을 SVG 문자열로 바꾸는 renderChartHtml 이 들어 있다. 시트로는 분석 창이 SVG를 PNG 로 바꿔 그림으로 넣는다 (Export.js).
 * 색과 선은 속성으로 직접 지정해 표를 복사할 때처럼 그림도 그대로 복사·저장되게 한다.
 */

var CHART_BAR_COLOR_ = '#4e79a7';
var CHART_POINT_COLOR_ = '#4e79a7';
var CHART_LINE_COLOR_ = '#e15759';
var CHART_AXIS_COLOR_ = '#5f6368';
var CHART_GRID_COLOR_ = '#e0e3e7';
var CHART_FONT_ = "'Malgun Gothic', 'Apple SD Gothic Neo', sans-serif";
var CHART_MAX_BARS_ = 30;
var CHART_MAX_POINTS_ = 3000;
var CHART_MAX_PAIRS_ = 10;
var CHART_SERIES_COLORS_ = ['#4e79a7', '#e15759', '#59a14f', '#f28e2b', '#b07aa1', '#76b7b2', '#9c755f', '#bab0ac'];
var CHART_MAX_SERIES_ = 8;

function chartEsc_(s) {
  return String(s).replace(/[&<>"']/g, function (c) {
    return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
  });
}

/* ---------- 계산 ---------- */

/** 보기 좋은 눈금: lo..hi 를 덮는 1, 2, 5 × 10^k 간격 */
function niceTicks_(lo, hi, count) {
  if (!(hi > lo)) { lo -= 0.5; hi += 0.5; }
  var raw = (hi - lo) / count;
  var mag = Math.pow(10, Math.floor(Math.log10(raw)));
  var norm = raw / mag;
  var step = (norm <= 1 ? 1 : norm <= 2 ? 2 : norm <= 5 ? 5 : 10) * mag;
  var start = Math.floor(lo / step + 1e-9) * step;
  var end = Math.ceil(hi / step - 1e-9) * step;
  var ticks = [];
  for (var i = 0; start + i * step <= end + step * 1e-9; i++) ticks.push(Number((start + i * step).toPrecision(12)));
  return { ticks: ticks, min: ticks[0], max: ticks[ticks.length - 1], step: step };
}

function tickLabel_(v) {
  return String(Number(v.toPrecision(10)));
}

function medianSorted_(a, from, to) { // a[from..to] (포함)
  var n = to - from + 1;
  var mid = from + (n - 1) / 2;
  return n % 2 ? a[mid] : (a[mid - 0.5] + a[mid + 0.5]) / 2;
}

/**
 * 상자그림 통계. 상자는 튜키의 힌지(중앙값을 기준으로 나눈 아래·위 절반의 중앙값, n이 홀수면 중앙값을 양쪽에 포함),
 * 수염은 울타리(상자에서 1.5 × 힌지 간격) 안의 최솟값·최댓값, 울타리 밖은 이상치(○),
 * 3 × 힌지 간격 밖은 극단값(*).
 */
function boxStats(values) {
  var a = values.slice().sort(function (x, y) { return x - y; });
  var n = a.length;
  var r = { n: n, q1: null, median: null, q3: null, lo: null, hi: null, outliers: [], extremes: [] };
  if (n === 0) return r;
  r.median = medianSorted_(a, 0, n - 1);
  var half = Math.ceil(n / 2);            // 아래 절반의 끝 (n이 홀수면 중앙값 포함)
  r.q1 = medianSorted_(a, 0, half - 1);
  r.q3 = medianSorted_(a, n - half, n - 1);
  var iqr = r.q3 - r.q1;
  var inner = [r.q1 - 1.5 * iqr, r.q3 + 1.5 * iqr], outer = [r.q1 - 3 * iqr, r.q3 + 3 * iqr];
  var within = [];
  a.forEach(function (v) {
    if (v >= inner[0] && v <= inner[1]) within.push(v);
    else if (v < outer[0] || v > outer[1]) r.extremes.push(v);
    else r.outliers.push(v);
  });
  r.lo = within.length ? within[0] : r.q1;
  r.hi = within.length ? within[within.length - 1] : r.q3;
  return r;
}

/** 히스토그램 구간: Sturges 구간 수를 기준으로 보기 좋은 간격을 고른다. 구간은 [하한, 상한), 마지막만 상한 포함. */
function histogramBins(values) {
  var n = values.length;
  // Math.min.apply(null, 배열)은 배열이 약 12만 개를 넘으면 호출 스택 한계로 RangeError 가 나므로 반복문으로 구한다
  var min = values[0], max = values[0];
  for (var q = 1; q < n; q++) {
    if (values[q] < min) min = values[q];
    else if (values[q] > max) max = values[q];
  }
  if (!(max > min)) return [{ lo: min - 0.5, hi: max + 0.5, count: n }];
  var k = Math.ceil(Math.log(n) / Math.LN2) + 1;
  var raw = (max - min) / k;
  var mag = Math.pow(10, Math.floor(Math.log10(raw)));
  var mults = [1, 2, 2.5, 5, 10], step = 10 * mag;
  for (var m = 0; m < mults.length; m++) {
    if (mults[m] * mag >= raw - 1e-12) { step = mults[m] * mag; break; }
  }
  var start = Math.floor(min / step + 1e-9) * step;
  var bins = [];
  for (var lo = start; lo <= max + step * 1e-9 && bins.length < 100; lo += step) {
    bins.push({ lo: Number(lo.toPrecision(12)), hi: Number((lo + step).toPrecision(12)), count: 0 });
  }
  values.forEach(function (v) {
    var idx = Math.floor((v - start) / step + 1e-9);
    if (idx >= bins.length) idx = bins.length - 1;
    if (idx < 0) idx = 0;
    bins[idx].count++;
  });
  return bins;
}

function meanSd_(values) {
  var n = values.length, s = 0, i;
  for (i = 0; i < n; i++) s += values[i];
  var mean = s / n, ss = 0;
  for (i = 0; i < n; i++) ss += (values[i] - mean) * (values[i] - mean);
  return { mean: mean, sd: n > 1 ? Math.sqrt(ss / (n - 1)) : null };
}

/* ---------- 차트 블록 (스펙) ---------- */

function barChartBlock(caption, labels, counts, xLabel) {
  return { type: 'chart', kind: 'bar', caption: caption, xLabel: xLabel || '', yLabel: '빈도', labels: labels, values: counts };
}

function histogramBlock(caption, values, xLabel) {
  var ms = meanSd_(values);
  return {
    type: 'chart', kind: 'histogram', caption: caption, xLabel: xLabel || '', yLabel: '빈도',
    bins: histogramBins(values), n: values.length, mean: ms.mean, sd: ms.sd
  };
}

function scatterBlock(caption, x, y, xLabel, yLabel, r) {
  var n = x.length, i;
  var mx = meanSd_(x).mean, my = meanSd_(y).mean, sxy = 0, sxx = 0;
  for (i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) * (x[i] - mx); }
  var line = sxx > 0 ? { slope: sxy / sxx, intercept: my - sxy / sxx * mx } : null;
  var points = [];
  var stride = Math.max(1, Math.ceil(n / CHART_MAX_POINTS_));
  for (i = 0; i < n; i += stride) points.push([x[i], y[i]]);
  return { type: 'chart', kind: 'scatter', caption: caption, xLabel: xLabel, yLabel: yLabel, points: points, n: n, line: line, r: r };
}

/**
 * 선 그래프: series = [{ label, points: [[x, y], ...] }] (조절효과의 수준별 회귀선 등).
 * xTicks(선택): 범주 이름 배열. 주면 x 눈금을 1, 2, … 위치에 이 이름으로 표시한다 (점의 x 는 1부터 시작하는 번호) — 이원 분산분석의 상호작용 그림
 */
function lineChartBlock(caption, xLabel, yLabel, series, xTicks) {
  var block = {
    type: 'chart', kind: 'line', caption: caption, xLabel: xLabel || '', yLabel: yLabel || '',
    series: series.slice(0, CHART_MAX_SERIES_).map(function (s) { return { label: s.label, points: s.points }; })
  };
  if (xTicks && xTicks.length) block.xTicks = xTicks.map(String);
  return block;
}

/**
 * IPA 사분면 그림: points = [{ label, x(수행도), y(중요도), q(사분면 번호 0~3) }], xLine·yLine = 십자 기준선,
 * 사분면 번호: 0 집중 개선(중요도 높음·수행도 낮음, 왼쪽 위), 1 현상 유지(오른쪽 위), 2 낮은 우선순위(왼쪽 아래), 3 과잉 투자(오른쪽 아래)
 */
var IPA_QUADRANTS_ = ['집중 개선', '현상 유지', '낮은 우선순위', '과잉 투자'];
var IPA_COLORS_ = ['#e15759', '#59a14f', '#8d8d8d', '#f28e2b'];
function ipaChartBlock(caption, points, xLine, yLine, xLabel, yLabel) {
  return { type: 'chart', kind: 'ipa', caption: caption, xLabel: xLabel || '수행도', yLabel: yLabel || '중요도', points: points, xLine: xLine, yLine: yLine };
}

/**
 * 경로도(매개효과): spec = { x, y, mediators: [{ label, a, b }], direct, total } (a, b, direct, total 은 표시할 계수 문구, 예: "0.512***")
 */
function pathDiagramBlock(caption, spec) {
  return { type: 'chart', kind: 'path', caption: caption, spec: spec };
}

/**
 * 나무 그림(의사결정나무). nodes: 전위 순서 [{ id, left, right (자식 id, 잎이면 0), edge (부모에서 이 노드로 오는 조건 문구), lines: [문구...], cls (분류 나무의 예측 범주 번호, 회귀 -1) }]
 */
function treeChartBlock(caption, nodes) {
  return { type: 'chart', kind: 'tree', caption: caption, nodes: nodes };
}

/** groups: [{ label, values }] */
function boxplotBlock(caption, groups, yLabel) {
  return {
    type: 'chart', kind: 'box', caption: caption, yLabel: yLabel,
    groups: groups.map(function (g) {
      var s = boxStats(g.values);
      return { label: g.label, n: s.n, q1: s.q1, median: s.median, q3: s.q3, lo: s.lo, hi: s.hi, outliers: s.outliers, extremes: s.extremes };
    })
  };
}

/* ---------- SVG ---------- */

function svgText_(x, y, text, opts) {
  opts = opts || {};
  return '<text x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" font-size="' + (opts.size || 11) + '"' +
    ' font-family="' + CHART_FONT_ + '" fill="' + (opts.fill || '#202124') + '"' +
    (opts.anchor ? ' text-anchor="' + opts.anchor + '"' : '') +
    (opts.weight ? ' font-weight="' + opts.weight + '"' : '') +
    (opts.rotate ? ' transform="rotate(' + opts.rotate + ' ' + x.toFixed(1) + ' ' + y.toFixed(1) + ')"' : '') +
    '>' + chartEsc_(text) + '</text>';
}

function scale_(d0, d1, r0, r1) {
  return function (v) { return r0 + (v - d0) / (d1 - d0) * (r1 - r0); };
}

function fmt_(v) { return Number(v.toFixed(2)).toString(); }

/** 공통 틀: 제목, y축(눈금·격자), x축 선. 안쪽 그림은 body 로 넣는다. */
function svgFrame_(W, H, m, caption, yTicks, ysc, xAxisLabel, yAxisLabel, body) {
  var plotW = W - m.left - m.right;
  var out = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' + W + ' ' + H + '" width="' + W + '" height="' + H + '" role="img">';
  out += '<title>' + chartEsc_(caption) + '</title>';
  out += '<rect x="0" y="0" width="' + W + '" height="' + H + '" fill="#ffffff"/>';
  out += svgText_(m.left, 18, caption, { size: 13, weight: 'bold' });
  yTicks.forEach(function (t) {
    var y = ysc(t);
    out += '<line x1="' + m.left + '" y1="' + y.toFixed(1) + '" x2="' + (m.left + plotW) + '" y2="' + y.toFixed(1) + '" stroke="' + CHART_GRID_COLOR_ + '"/>';
    out += svgText_(m.left - 6, y + 4, tickLabel_(t), { anchor: 'end', fill: CHART_AXIS_COLOR_ });
  });
  out += body;
  out += '<line x1="' + m.left + '" y1="' + ysc(yTicks[0]).toFixed(1) + '" x2="' + (m.left + plotW) + '" y2="' + ysc(yTicks[0]).toFixed(1) + '" stroke="' + CHART_AXIS_COLOR_ + '"/>';
  out += '<line x1="' + m.left + '" y1="' + ysc(yTicks[yTicks.length - 1]).toFixed(1) + '" x2="' + m.left + '" y2="' + ysc(yTicks[0]).toFixed(1) + '" stroke="' + CHART_AXIS_COLOR_ + '"/>';
  if (xAxisLabel) out += svgText_(m.left + plotW / 2, H - 8, xAxisLabel, { anchor: 'middle', fill: CHART_AXIS_COLOR_ });
  if (yAxisLabel) out += svgText_(14, (m.top + (H - m.bottom)) / 2, yAxisLabel, { anchor: 'middle', fill: CHART_AXIS_COLOR_, rotate: -90 });
  return out;
}

function renderBar_(b) {
  var n = b.labels.length;
  var rotate = n > 8 || b.labels.some(function (l) { return String(l).length > 6; });
  var m = { left: 56, right: 16, top: 30, bottom: rotate ? 84 : 52 };
  var W = Math.max(380, Math.min(760, 70 + n * 44)), H = 300;
  var maxV = Math.max.apply(null, b.values);
  var yt = niceTicks_(0, maxV, 5);
  var ysc = scale_(yt.min, yt.max, H - m.bottom, m.top);
  var plotW = W - m.left - m.right, band = plotW / n, bw = band * 0.7;
  var body = '';
  b.values.forEach(function (v, i) {
    var x = m.left + i * band + (band - bw) / 2, y = ysc(v);
    body += '<rect x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' + bw.toFixed(1) + '" height="' + (ysc(0) - y).toFixed(1) + '" fill="' + CHART_BAR_COLOR_ + '"/>';
    if (n <= 20) body += svgText_(x + bw / 2, y - 4, String(v), { size: 10, anchor: 'middle' });
    var lab = String(b.labels[i]);
    if (lab.length > 12) lab = lab.slice(0, 11) + '…';
    var lx = x + bw / 2, ly = H - m.bottom + 14;
    body += rotate ? svgText_(lx, ly, lab, { size: 10, anchor: 'end', rotate: -35 }) : svgText_(lx, ly, lab, { size: 10, anchor: 'middle' });
  });
  return svgFrame_(W, H, m, b.caption, yt.ticks, ysc, b.xLabel, b.yLabel, body) + '</svg>';
}

function renderHistogram_(b) {
  var m = { left: 56, right: 16, top: 30, bottom: 52 };
  var W = 480, H = 300;
  var x0 = b.bins[0].lo, x1 = b.bins[b.bins.length - 1].hi;
  var width = b.bins[0].hi - b.bins[0].lo;
  var curve = [];
  var peak = Math.max.apply(null, b.bins.map(function (bin) { return bin.count; }));
  if (b.sd > 0) {
    var pdfMax = b.n * width / (b.sd * Math.sqrt(2 * Math.PI));
    for (var i = 0; i <= 80; i++) {
      var x = x0 + (x1 - x0) * i / 80;
      var z = (x - b.mean) / b.sd;
      curve.push([x, b.n * width / (b.sd * Math.sqrt(2 * Math.PI)) * Math.exp(-z * z / 2)]);
    }
    peak = Math.max(peak, pdfMax);
  }
  var yt = niceTicks_(0, peak, 5);
  var ysc = scale_(yt.min, yt.max, H - m.bottom, m.top);
  var xsc = scale_(x0, x1, m.left, W - m.right);
  var body = '';
  b.bins.forEach(function (bin, idx) {
    body += '<rect x="' + xsc(bin.lo).toFixed(1) + '" y="' + ysc(bin.count).toFixed(1) + '" width="' + (xsc(bin.hi) - xsc(bin.lo)).toFixed(1) + '" height="' + (ysc(0) - ysc(bin.count)).toFixed(1) + '" fill="' + CHART_BAR_COLOR_ + '" stroke="#ffffff"/>';
  });
  var every = b.bins.length > 12 ? 2 : 1;
  for (var k = 0; k <= b.bins.length; k += every) {
    var edge = k < b.bins.length ? b.bins[k].lo : b.bins[b.bins.length - 1].hi;
    body += svgText_(xsc(edge), H - m.bottom + 14, tickLabel_(edge), { size: 10, anchor: 'middle', fill: CHART_AXIS_COLOR_ });
  }
  if (curve.length) {
    body += '<polyline fill="none" stroke="' + CHART_LINE_COLOR_ + '" stroke-width="1.5" points="' + curve.map(function (p) { return xsc(p[0]).toFixed(1) + ',' + ysc(p[1]).toFixed(1); }).join(' ') + '"/>';
  }
  var tx = W - m.right - 4;
  body += svgText_(tx, m.top + 12, '평균 = ' + fmt_(b.mean), { size: 10, anchor: 'end' });
  if (b.sd !== null) body += svgText_(tx, m.top + 25, '표준편차 = ' + fmt_(b.sd), { size: 10, anchor: 'end' });
  body += svgText_(tx, m.top + 38, 'N = ' + b.n, { size: 10, anchor: 'end' });
  return svgFrame_(W, H, m, b.caption, yt.ticks, ysc, b.xLabel, b.yLabel, body) + '</svg>';
}

function renderScatter_(b) {
  var m = { left: 56, right: 16, top: 34, bottom: 52 };
  var W = 420, H = 320;
  var xs = b.points.map(function (p) { return p[0]; }), ys = b.points.map(function (p) { return p[1]; });
  var xt = niceTicks_(Math.min.apply(null, xs), Math.max.apply(null, xs), 6);
  var yt = niceTicks_(Math.min.apply(null, ys), Math.max.apply(null, ys), 5);
  var xsc = scale_(xt.min, xt.max, m.left, W - m.right);
  var ysc = scale_(yt.min, yt.max, H - m.bottom, m.top);
  var body = '';
  xt.ticks.forEach(function (t) {
    body += '<line x1="' + xsc(t).toFixed(1) + '" y1="' + ysc(yt.min).toFixed(1) + '" x2="' + xsc(t).toFixed(1) + '" y2="' + (ysc(yt.min) + 4).toFixed(1) + '" stroke="' + CHART_AXIS_COLOR_ + '"/>';
    body += svgText_(xsc(t), H - m.bottom + 16, tickLabel_(t), { size: 10, anchor: 'middle', fill: CHART_AXIS_COLOR_ });
  });
  b.points.forEach(function (p) {
    body += '<circle cx="' + xsc(p[0]).toFixed(1) + '" cy="' + ysc(p[1]).toFixed(1) + '" r="3" fill="' + CHART_POINT_COLOR_ + '" fill-opacity="0.55"/>';
  });
  if (b.line) {
    var ya = b.line.intercept + b.line.slope * xt.min, yb = b.line.intercept + b.line.slope * xt.max;
    body += '<line x1="' + xsc(xt.min).toFixed(1) + '" y1="' + ysc(Math.max(yt.min, Math.min(yt.max, ya))).toFixed(1) +
      '" x2="' + xsc(xt.max).toFixed(1) + '" y2="' + ysc(Math.max(yt.min, Math.min(yt.max, yb))).toFixed(1) +
      '" stroke="' + CHART_LINE_COLOR_ + '" stroke-width="1.5"/>';
  }
  var sub = 'N = ' + b.n + (b.r !== null && b.r !== undefined ? ', r = ' + b.r.toFixed(3) : '') + (b.points.length < b.n ? ' (점 ' + b.points.length + '개만 표시)' : '');
  body += svgText_(W - m.right - 4, m.top - 6, sub, { size: 10, anchor: 'end', fill: CHART_AXIS_COLOR_ });
  return svgFrame_(W, H, m, b.caption, yt.ticks, ysc, b.xLabel, b.yLabel, body) + '</svg>';
}

/** 범례 글자의 대략적인 폭(px): 한글은 넓고 영문·숫자는 좁다 */
function legendTextWidth_(text) {
  var w = 0;
  for (var i = 0; i < text.length; i++) w += text.charCodeAt(i) > 255 ? 10 : 5.8;
  return w;
}

function renderLine_(b) {
  var m = { left: 56, right: 16, top: 44, bottom: 52 };
  var legendW = 0;
  var labels = b.series.map(function (s) {
    var lab = String(s.label);
    return lab.length > 22 ? lab.slice(0, 21) + '…' : lab;
  });
  labels.forEach(function (lab) { legendW += 30 + legendTextWidth_(lab); });
  var W = Math.max(460, Math.ceil(m.left + legendW + m.right)), H = 330;
  var xs = [], ys = [];
  b.series.forEach(function (s) { s.points.forEach(function (p) { xs.push(p[0]); ys.push(p[1]); }); });
  var cat = b.xTicks && b.xTicks.length ? b.xTicks : null;
  var xt = cat
    ? { min: 0.5, max: cat.length + 0.5, ticks: cat.map(function (_, k) { return k + 1; }) }
    : niceTicks_(Math.min.apply(null, xs), Math.max.apply(null, xs), 6);
  var yt = niceTicks_(Math.min.apply(null, ys), Math.max.apply(null, ys), 5);
  var xsc = scale_(xt.min, xt.max, m.left, W - m.right);
  var ysc = scale_(yt.min, yt.max, H - m.bottom, m.top);
  var body = '';
  xt.ticks.forEach(function (t) {
    var tl = cat ? (cat[t - 1].length > 12 ? cat[t - 1].slice(0, 11) + '…' : cat[t - 1]) : tickLabel_(t);
    body += '<line x1="' + xsc(t).toFixed(1) + '" y1="' + ysc(yt.min).toFixed(1) + '" x2="' + xsc(t).toFixed(1) + '" y2="' + (ysc(yt.min) + 4).toFixed(1) + '" stroke="' + CHART_AXIS_COLOR_ + '"/>';
    body += svgText_(xsc(t), H - m.bottom + 16, tl, { size: 10, anchor: 'middle', fill: CHART_AXIS_COLOR_ });
  });
  var lx = m.left;
  b.series.forEach(function (s, i) {
    var color = CHART_SERIES_COLORS_[i % CHART_SERIES_COLORS_.length];
    body += '<polyline fill="none" stroke="' + color + '" stroke-width="2" points="' +
      s.points.map(function (p) { return xsc(p[0]).toFixed(1) + ',' + ysc(p[1]).toFixed(1); }).join(' ') + '"/>';
    var lab = labels[i];
    body += '<line x1="' + lx + '" y1="' + (m.top - 14) + '" x2="' + (lx + 14) + '" y2="' + (m.top - 14) + '" stroke="' + color + '" stroke-width="2"/>';
    body += svgText_(lx + 18, m.top - 10, lab, { size: 10, fill: CHART_AXIS_COLOR_ });
    lx += 30 + legendTextWidth_(lab);
  });
  return svgFrame_(W, H, m, b.caption, yt.ticks, ysc, b.xLabel, b.yLabel, body) + '</svg>';
}

/** IPA 사분면 그림: 십자 기준선, 사분면 이름, 문항 이름이 붙은 점 (이름이 겹치면 오른쪽, 왼쪽, 위, 아래 순으로 빈 곳을 찾는다) */
function renderIpa_(b) {
  var m = { left: 60, right: 22, top: 40, bottom: 54 };
  var W = 600, H = 480;
  var xs = b.points.map(function (p) { return p.x; }).concat([b.xLine]);
  var ys = b.points.map(function (p) { return p.y; }).concat([b.yLine]);
  var padX = (Math.max.apply(null, xs) - Math.min.apply(null, xs)) * 0.08 || 0.5;
  var padY = (Math.max.apply(null, ys) - Math.min.apply(null, ys)) * 0.08 || 0.5;
  var xt = niceTicks_(Math.min.apply(null, xs) - padX, Math.max.apply(null, xs) + padX, 6);
  var yt = niceTicks_(Math.min.apply(null, ys) - padY, Math.max.apply(null, ys) + padY, 6);
  var xsc = scale_(xt.min, xt.max, m.left, W - m.right);
  var ysc = scale_(yt.min, yt.max, H - m.bottom, m.top);
  var cx = xsc(b.xLine), cy = ysc(b.yLine), x0 = m.left, x1 = W - m.right, y0 = m.top, y1 = H - m.bottom;
  var body = '';
  // 사분면 바탕색 (왼쪽 위가 집중 개선)
  var quads = [[x0, y0, cx, cy, 0], [cx, y0, x1, cy, 1], [x0, cy, cx, y1, 2], [cx, cy, x1, y1, 3]];
  quads.forEach(function (q) {
    body += '<rect x="' + q[0].toFixed(1) + '" y="' + q[1].toFixed(1) + '" width="' + Math.max(0, q[2] - q[0]).toFixed(1) + '" height="' + Math.max(0, q[3] - q[1]).toFixed(1) +
      '" fill="' + IPA_COLORS_[q[4]] + '" fill-opacity="0.07"/>';
  });
  xt.ticks.forEach(function (t) {
    body += '<line x1="' + xsc(t).toFixed(1) + '" y1="' + y1.toFixed(1) + '" x2="' + xsc(t).toFixed(1) + '" y2="' + (y1 + 4).toFixed(1) + '" stroke="' + CHART_AXIS_COLOR_ + '"/>';
    body += svgText_(xsc(t), y1 + 17, tickLabel_(t), { size: 10, anchor: 'middle', fill: CHART_AXIS_COLOR_ });
  });
  body += '<line x1="' + cx.toFixed(1) + '" y1="' + y0 + '" x2="' + cx.toFixed(1) + '" y2="' + y1 + '" stroke="#202124" stroke-width="1.2" stroke-dasharray="5,4"/>';
  body += '<line x1="' + x0 + '" y1="' + cy.toFixed(1) + '" x2="' + x1 + '" y2="' + cy.toFixed(1) + '" stroke="#202124" stroke-width="1.2" stroke-dasharray="5,4"/>';
  body += svgText_(Math.min(x1 - 2, cx + 4), y0 + 26, '수행도 기준 ' + fmt_(b.xLine), { size: 10, fill: CHART_AXIS_COLOR_, anchor: cx + 90 > x1 ? 'end' : 'start' });
  body += svgText_(x1 - 4, cy - 4, '중요도 기준 ' + fmt_(b.yLine), { size: 10, fill: CHART_AXIS_COLOR_, anchor: 'end' });
  body += svgText_(x0 + 6, y0 + 14, '집중 개선', { size: 12, weight: 'bold', fill: IPA_COLORS_[0] });
  body += svgText_(x1 - 6, y0 + 14, '현상 유지', { size: 12, weight: 'bold', fill: IPA_COLORS_[1], anchor: 'end' });
  body += svgText_(x0 + 6, y1 - 8, '낮은 우선순위', { size: 12, weight: 'bold', fill: IPA_COLORS_[2] });
  body += svgText_(x1 - 6, y1 - 8, '과잉 투자', { size: 12, weight: 'bold', fill: IPA_COLORS_[3], anchor: 'end' });
  // 점과 이름
  var placed = [], px = [];
  var overlap = function (a, c) { return a.x0 < c.x1 && c.x0 < a.x1 && a.y0 < c.y1 && c.y0 < a.y1; };
  b.points.forEach(function (p) { var X = xsc(p.x), Y = ysc(p.y); px.push({ x0: X - 5, x1: X + 5, y0: Y - 5, y1: Y + 5 }); });
  b.points.forEach(function (p, i) {
    var X = xsc(p.x), Y = ysc(p.y);
    var lab = String(p.label);
    if (lab.length > 14) lab = lab.slice(0, 13) + '…';
    var w = legendTextWidth_(lab) + 4, h = 12;
    var cands = [
      { x: X + 8, y: Y + 4, anchor: 'start', box: { x0: X + 7, x1: X + 7 + w, y0: Y - 8, y1: Y + 5 } },
      { x: X - 8, y: Y + 4, anchor: 'end', box: { x0: X - 7 - w, x1: X - 7, y0: Y - 8, y1: Y + 5 } },
      { x: X, y: Y - 9, anchor: 'middle', box: { x0: X - w / 2, x1: X + w / 2, y0: Y - 9 - h, y1: Y - 7 } },
      { x: X, y: Y + 18, anchor: 'middle', box: { x0: X - w / 2, x1: X + w / 2, y0: Y + 7, y1: Y + 19 } }
    ];
    var best = null, bestScore = Infinity;
    cands.forEach(function (c) {
      var score = 0;
      if (c.box.x0 < x0 + 1 || c.box.x1 > x1 - 1) score += 1000;
      placed.forEach(function (o) { if (overlap(c.box, o)) score += 10; });
      px.forEach(function (o, k) { if (k !== i && overlap(c.box, o)) score += 5; });
      if (score < bestScore) { bestScore = score; best = c; }
    });
    placed.push(best.box);
    body += '<circle cx="' + X.toFixed(1) + '" cy="' + Y.toFixed(1) + '" r="5" fill="' + IPA_COLORS_[p.q] + '" stroke="#ffffff" stroke-width="1"/>';
    body += svgText_(best.x, best.y, lab, { size: 11, anchor: best.anchor, fill: '#202124' });
  });
  return svgFrame_(W, H, m, b.caption, yt.ticks, ysc, b.xLabel, b.yLabel, body) + '</svg>';
}

function renderPath_(b) {
  var s = b.spec, m = s.mediators.length;
  var boxW = 128, boxH = 38;
  var W = Math.max(560, 90 + m * (boxW + 24) + 2 * boxW), H = 270;
  var yMed = 66, yBase = 210;
  var xBox = { cx: 24 + boxW / 2, cy: yBase }, yBox = { cx: W - 24 - boxW / 2, cy: yBase };
  function box(cx, cy, label) {
    var lab = String(label);
    if (lab.length > 11) lab = lab.slice(0, 10) + '…';
    return '<rect x="' + (cx - boxW / 2).toFixed(1) + '" y="' + (cy - boxH / 2).toFixed(1) + '" width="' + boxW + '" height="' + boxH + '" rx="4" fill="#f1f3f4" stroke="' + CHART_AXIS_COLOR_ + '"/>' +
      svgText_(cx, cy + 4, lab, { size: 12, anchor: 'middle' });
  }
  function arrow(x1, y1, x2, y2) {
    return '<line x1="' + x1.toFixed(1) + '" y1="' + y1.toFixed(1) + '" x2="' + x2.toFixed(1) + '" y2="' + y2.toFixed(1) + '" stroke="#202124" stroke-width="1.5" marker-end="url(#pathArrow)"/>';
  }
  var body = '<defs><marker id="pathArrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#202124"/></marker></defs>';
  body += svgText_(W / 2, 40, '* p<.05  ** p<.01  *** p<.001 (비표준화 계수)', { size: 10, anchor: 'middle', fill: CHART_AXIS_COLOR_ });
  s.mediators.forEach(function (md, i) {
    var cx = W * (i + 1) / (m + 1);
    body += arrow(xBox.cx + boxW / 4, yBase - boxH / 2, cx - boxW / 4, yMed + boxH / 2);
    body += arrow(cx + boxW / 4, yMed + boxH / 2, yBox.cx - boxW / 4, yBase - boxH / 2);
    body += svgText_((xBox.cx + cx) / 2 - 6, (yBase + yMed) / 2 + 2, 'a' + (m > 1 ? i + 1 : '') + ' = ' + md.a, { size: 11, anchor: 'end' });
    body += svgText_((yBox.cx + cx) / 2 + 6, (yBase + yMed) / 2 + 2, 'b' + (m > 1 ? i + 1 : '') + ' = ' + md.b, { size: 11, anchor: 'start' });
    body += box(cx, yMed, md.label);
  });
  body += arrow(xBox.cx + boxW / 2, yBase, yBox.cx - boxW / 2, yBase);
  body += svgText_(W / 2, yBase - 6, 'c\' = ' + s.direct + ' (직접효과)', { size: 11, anchor: 'middle' });
  body += svgText_(W / 2, yBase + 18, 'c = ' + s.total + ' (총효과)', { size: 11, anchor: 'middle', fill: CHART_AXIS_COLOR_ });
  body += box(xBox.cx, xBox.cy, s.x) + box(yBox.cx, yBox.cy, s.y);
  return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' + W + ' ' + H + '" width="' + W + '" height="' + H + '" role="img">' +
    '<title>' + chartEsc_(b.caption) + '</title><rect x="0" y="0" width="' + W + '" height="' + H + '" fill="#ffffff"/>' +
    svgText_(16, 18, b.caption, { size: 13, weight: 'bold' }) + body + '</svg>';
}

function renderBox_(b) {
  var m = { left: 56, right: 16, top: 30, bottom: 56 };
  var k = b.groups.length;
  var W = Math.max(380, Math.min(760, 90 + k * 90)), H = 320;
  var lo = Infinity, hi = -Infinity;
  b.groups.forEach(function (g) {
    if (g.n === 0) return;
    var all = [g.lo, g.hi].concat(g.outliers, g.extremes);
    all.forEach(function (v) { if (v < lo) lo = v; if (v > hi) hi = v; });
  });
  var yt = niceTicks_(lo, hi, 5);
  var ysc = scale_(yt.min, yt.max, H - m.bottom, m.top);
  var band = (W - m.left - m.right) / k, bw = Math.min(60, band * 0.5);
  var body = '';
  b.groups.forEach(function (g, i) {
    var cx = m.left + (i + 0.5) * band;
    var lab = String(g.label);
    if (lab.length > 10) lab = lab.slice(0, 9) + '…';
    body += svgText_(cx, H - m.bottom + 16, lab, { size: 10, anchor: 'middle' });
    body += svgText_(cx, H - m.bottom + 30, 'N=' + g.n, { size: 9, anchor: 'middle', fill: CHART_AXIS_COLOR_ });
    if (g.n === 0) return;
    var s = CHART_AXIS_COLOR_;
    body += '<line x1="' + cx + '" y1="' + ysc(g.hi).toFixed(1) + '" x2="' + cx + '" y2="' + ysc(g.q3).toFixed(1) + '" stroke="' + s + '"/>';
    body += '<line x1="' + cx + '" y1="' + ysc(g.q1).toFixed(1) + '" x2="' + cx + '" y2="' + ysc(g.lo).toFixed(1) + '" stroke="' + s + '"/>';
    body += '<line x1="' + (cx - bw / 4).toFixed(1) + '" y1="' + ysc(g.hi).toFixed(1) + '" x2="' + (cx + bw / 4).toFixed(1) + '" y2="' + ysc(g.hi).toFixed(1) + '" stroke="' + s + '"/>';
    body += '<line x1="' + (cx - bw / 4).toFixed(1) + '" y1="' + ysc(g.lo).toFixed(1) + '" x2="' + (cx + bw / 4).toFixed(1) + '" y2="' + ysc(g.lo).toFixed(1) + '" stroke="' + s + '"/>';
    body += '<rect x="' + (cx - bw / 2).toFixed(1) + '" y="' + ysc(g.q3).toFixed(1) + '" width="' + bw.toFixed(1) + '" height="' + Math.max(1, ysc(g.q1) - ysc(g.q3)).toFixed(1) + '" fill="#cfe0f1" stroke="' + s + '"/>';
    body += '<line x1="' + (cx - bw / 2).toFixed(1) + '" y1="' + ysc(g.median).toFixed(1) + '" x2="' + (cx + bw / 2).toFixed(1) + '" y2="' + ysc(g.median).toFixed(1) + '" stroke="#202124" stroke-width="2"/>';
    g.outliers.forEach(function (v) {
      body += '<circle cx="' + cx + '" cy="' + ysc(v).toFixed(1) + '" r="3" fill="none" stroke="' + s + '"/>';
    });
    g.extremes.forEach(function (v) {
      body += svgText_(cx, ysc(v) + 4, '*', { size: 13, anchor: 'middle' });
    });
  });
  return svgFrame_(W, H, m, b.caption, yt.ticks, ysc, '', b.yLabel, body) + '</svg>';
}

var TREE_PALETTE_ = ['#cfe0f1', '#fde4cf', '#d9efd5', '#f3d9ef', '#fff2b8', '#d6f0f0', '#e6dcc8', '#e0e0e0', '#f6d2d2', '#d8d8f0'];

function renderTree_(b) {
  var nodes = b.nodes, byId = {}, boxW = 156, gap = 16, lineH = 14, topPad = 34, edgeGap = 36;
  nodes.forEach(function (nd) { byId[nd.id] = nd; });
  var maxDepth = 0, leafSlot = 0;
  var place = function (nd, depth) {
    nd.depth = depth;
    if (depth > maxDepth) maxDepth = depth;
    if (!nd.left) { nd.x = leafSlot++; return; }
    place(byId[nd.left], depth + 1);
    place(byId[nd.right], depth + 1);
    nd.x = (byId[nd.left].x + byId[nd.right].x) / 2;
  };
  place(nodes[0], 0);
  var levelH = new Array(maxDepth + 1).fill(0);
  nodes.forEach(function (nd) { nd.h = nd.lines.length * lineH + 10; if (nd.h > levelH[nd.depth]) levelH[nd.depth] = nd.h; });
  var levelY = [], y = topPad;
  levelH.forEach(function (h, i) { levelY[i] = y; y += h + edgeGap; });
  var W = Math.max(480, leafSlot * (boxW + gap) + gap), H = y - edgeGap + 16;
  var cx = function (nd) { return gap + boxW / 2 + nd.x * (boxW + gap); };
  var trunc = function (s, n) { s = String(s); return s.length > n ? s.slice(0, n - 1) + '…' : s; };
  var body = '';
  nodes.forEach(function (nd) {
    if (!nd.left) return;
    [byId[nd.left], byId[nd.right]].forEach(function (ch) {
      var x1 = cx(nd), y1 = levelY[nd.depth] + nd.h, x2 = cx(ch), y2 = levelY[ch.depth];
      body += '<line x1="' + x1.toFixed(1) + '" y1="' + y1.toFixed(1) + '" x2="' + x2.toFixed(1) + '" y2="' + y2.toFixed(1) + '" stroke="' + CHART_AXIS_COLOR_ + '"/>';
      body += svgText_((x1 + x2) / 2, (y1 + y2) / 2 + 3, trunc(ch.edge, 26), { size: 9.5, anchor: 'middle', fill: CHART_AXIS_COLOR_ });
    });
  });
  nodes.forEach(function (nd) {
    var x = cx(nd) - boxW / 2, y0 = levelY[nd.depth];
    body += '<rect x="' + x.toFixed(1) + '" y="' + y0.toFixed(1) + '" width="' + boxW + '" height="' + nd.h + '" rx="4" fill="' + (nd.cls >= 0 ? TREE_PALETTE_[nd.cls % TREE_PALETTE_.length] : '#e8f0fe') + '" stroke="' + CHART_AXIS_COLOR_ + '"' + (nd.left ? '' : ' stroke-width="1.6"') + '/>';
    nd.lines.forEach(function (ln, i) { body += svgText_(cx(nd), y0 + 15 + i * lineH, trunc(ln, 22), { size: 10.5, anchor: 'middle', weight: i === 0 && nd.left ? 'bold' : undefined }); });
  });
  return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' + W + ' ' + H + '" width="' + W + '" height="' + H + '" role="img">' +
    '<title>' + chartEsc_(b.caption) + '</title><rect x="0" y="0" width="' + W + '" height="' + H + '" fill="#ffffff"/>' +
    svgText_(16, 18, b.caption, { size: 13, weight: 'bold' }) + body + '</svg>';
}

/** 차트 블록 → SVG 를 담은 HTML. 데이터가 비어 있으면 안내 문구. */
function renderChartHtml(block) {
  var svg;
  switch (block.kind) {
    case 'bar': svg = block.labels.length ? renderBar_(block) : null; break;
    case 'histogram': svg = block.n ? renderHistogram_(block) : null; break;
    case 'scatter': svg = block.points.length ? renderScatter_(block) : null; break;
    case 'line': svg = block.series.length && block.series.some(function (s) { return s.points.length > 1; }) ? renderLine_(block) : null; break;
    case 'path': svg = block.spec && block.spec.mediators.length ? renderPath_(block) : null; break;
    case 'ipa': svg = block.points && block.points.length ? renderIpa_(block) : null; break;
    case 'tree': svg = block.nodes && block.nodes.length ? renderTree_(block) : null; break;
    case 'box': svg = block.groups.some(function (g) { return g.n > 0; }) ? renderBox_(block) : null; break;
    default: svg = null;
  }
  if (!svg) return '<p class="note">' + chartEsc_(block.caption) + ': 그릴 데이터가 없습니다.</p>';
  return '<div class="chart">' + svg + '</div>';
}

if (typeof module !== 'undefined') {
  module.exports = {
    niceTicks_: niceTicks_,
    boxStats: boxStats,
    histogramBins: histogramBins,
    barChartBlock: barChartBlock,
    histogramBlock: histogramBlock,
    scatterBlock: scatterBlock,
    boxplotBlock: boxplotBlock,
    lineChartBlock: lineChartBlock,
    ipaChartBlock: ipaChartBlock,
    pathDiagramBlock: pathDiagramBlock,
    treeChartBlock: treeChartBlock,
    renderChartHtml: renderChartHtml,
    CHART_MAX_BARS_: CHART_MAX_BARS_,
    CHART_MAX_PAIRS_: CHART_MAX_PAIRS_
  };
}
