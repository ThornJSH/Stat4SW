/**
 * 결과 객체 → HTML (순수 함수). 형식은 doc/설계-결과객체.md 참고.
 * 스타일은 style.html 의 .result / table.spss 클래스가 담당한다.
 */

function escapeHtml_(s) {
  return String(s).replace(/[&<>"']/g, function (c) {
    return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
  });
}

/** toFixed 결과에서 "-0.000" 같은 음의 0을 없앤다 */
function toFixedClean_(v, digits) {
  var s = v.toFixed(digits);
  return /^-0(\.0*)?$/.test(s) ? s.slice(1) : s;
}

function formatCell_(v, fmt) {
  if (v === null || v === undefined || v === '') return '';
  if (typeof v === 'number' && !isFinite(v)) return '.';
  var type = fmt && fmt.type ? fmt.type : 'text';
  if (type === 'text' || typeof v !== 'number') return String(v);
  var digits = fmt.digits === undefined ? 3 : fmt.digits;
  switch (type) {
    case 'int': return toFixedClean_(Math.round(v), 0);
    case 'fixed': return toFixedClean_(v, digits);
    case 'trim': return String(Number(v.toFixed(digits)));
    case 'pvalue':
      if (v < 0.001) return '<.001';
      if (v >= 1) return '1.000';
      return v.toFixed(3).replace(/^0/, '');
    default: return String(v);
  }
}

function renderTable_(block) {
  var cols = block.columns;
  var html = '<table class="spss">';
  if (block.caption) html += '<caption>' + escapeHtml_(block.caption) + '</caption>';
  html += '<thead><tr>' + cols.map(function (c) {
    return '<th class="' + alignOf_(c) + '">' + escapeHtml_(c.label) + '</th>';
  }).join('') + '</tr></thead><tbody>';
  block.rows.forEach(function (row) {
    var classes = [];
    if (row._kind) classes.push(escapeHtml_(row._kind));
    if (row._bold) classes.push('strong');
    var cls = classes.length ? ' class="' + classes.join(' ') + '"' : '';
    html += '<tr' + cls + '>' + cols.map(function (c) {
      return '<td class="' + alignOf_(c) + '">' + escapeHtml_(formatCell_(row[c.key], row._fmt || c.format)) + '</td>';
    }).join('') + '</tr>';
  });
  html += '</tbody></table>';
  (block.notes || []).forEach(function (n) {
    html += '<p class="note">' + escapeHtml_(n) + '</p>';
  });
  return html;
}

function alignOf_(col) {
  if (col.align) return col.align;
  return col.format && col.format.type && col.format.type !== 'text' ? 'right' : 'left';
}

/** 차트 블록 → SVG HTML (Charts.js). 다른 파일의 함수라서 호출 시점에 찾는다. */
function chartHtml_(block) {
  return typeof module !== 'undefined' ? require('./Charts.js').renderChartHtml(block) : renderChartHtml(block);
}

function renderResultHtml(result) {
  var html = '<div class="result" data-analysis="' + escapeHtml_(result.analysis || '') + '">';
  html += '<h2>' + escapeHtml_(result.title) + '</h2>';
  var chartNo = 0;   // 차트 블록의 순번. 화면의 SVG와 결과 객체의 차트를 짝지어 시트로 내보낼 때 쓴다
  result.blocks.forEach(function (b) {
    if (b.type === 'table') html += renderTable_(b);
    else if (b.type === 'chart') html += '<div class="chart-slot" data-chart="' + (chartNo++) + '">' + chartHtml_(b) + '</div>';
    else html += '<p>' + escapeHtml_(b.text) + '</p>';
  });
  (result.warnings || []).forEach(function (w) {
    html += '<p class="warn">' + escapeHtml_(w) + '</p>';
  });
  return html + '</div>';
}

if (typeof module !== 'undefined') {
  module.exports = { renderResultHtml: renderResultHtml, formatCell_: formatCell_ };
}
