/**
 * 결과를 시트의 "결과" 탭에 쓴다.
 * resultToGrid: 결과 객체 → 격자(값, 숫자 서식, 행 종류, 정렬). 순수 함수라 Node에서 테스트한다.
 * exportResults / writeGrid_: 격자를 실제 시트에 쓰고 서식을 입힌다 (SpreadsheetApp 사용, 시트에서만 확인 가능).
 *
 * 숫자는 숫자로 쓰고 서식(정수, 소수점)을 입힌다. p값은 HTML과 같은 문자열(<.001, .043)로 쓴다.
 * 차트는 그림(PNG)으로 넣는다. 서버(Apps Script)는 SVG를 이미지로 바꿀 수 없으므로 대화상자(브라우저)에서
 * SVG를 PNG로 바꿔 함께 보내고, 서버는 그림 높이만큼 행을 비워 두고 insertImage 로 삽입한다.
 * 변환하지 못한 차트는 안내 줄만 남긴다.
 */

var EXPORT_ROW_PX_ = 21;          // 시트의 기본 행 높이(px). 그림이 차지할 행 수를 계산할 때 쓴다
var EXPORT_MAX_PNG_CHARS_ = 2000000;
var EXPORT_MAX_IMAGES_ = 30;

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function exportDeps_() {
  return typeof module !== 'undefined'
    ? { formatCell_: require('./Render.js').formatCell_, asUserError_: require('./Safety.js').asUserError_ }
    : { formatCell_: formatCell_, asUserError_: asUserError_ };
}

/** 값 하나를 시트 값과 숫자 서식으로. 서식이 null 이면 General. */
function exportCell_(v, fmt) {
  if (v === null || v === undefined || v === '') return { value: '', numFmt: null };
  if (typeof v === 'number') {
    if (!isFinite(v)) return { value: '.', numFmt: null };
    var type = fmt && fmt.type ? fmt.type : 'text';
    var digits = fmt && fmt.digits !== undefined ? fmt.digits : 3;
    if (type === 'int') return { value: Math.round(v), numFmt: '0' };
    if (type === 'fixed') return { value: v, numFmt: digits > 0 ? '0.' + new Array(digits + 1).join('0') : '0' };
    if (type === 'trim') return { value: Number(v.toFixed(digits)), numFmt: null };
    if (type === 'pvalue') return { value: exportDeps_().formatCell_(v, fmt), numFmt: '@' };
    return { value: v, numFmt: null };
  }
  var s = String(v);
  // =, +, -, @ 로 시작하는 문자열이 수식으로 해석되지 않도록 텍스트 서식(@)을 쓴다
  return { value: s, numFmt: /^[=+\-@]/.test(s) ? '@' : null };
}

function exportAlign_(col) {
  if (col.align) return col.align;
  return col.format && col.format.type && col.format.type !== 'text' ? 'right' : 'left';
}

/**
 * 결과 객체 → 격자. rows[i] = { kind, cells, formats, aligns }
 * kind: title, caption, header, body, group, subtotal, total, note, warn, text, image, blank
 * stamp: 제목 옆에 붙일 시각 문자열 (없으면 생략)
 * chartDims: 차트 블록 순서대로 { width, height }(그림으로 넣을 크기) 또는 null(그림 없음)
 * 반환: { rows, width, images } — images[k] = { rowIndex, rows, chartIndex, width, height } (그림이 놓일 자리)
 */
function resultToGrid(result, stamp, chartDims) {
  var rows = [], images = [];
  var width = 1, chartIdx = 0;
  function add(kind, cells, formats, aligns, bold) {
    rows.push({ kind: kind, cells: cells, formats: formats || null, aligns: aligns || null, bold: !!bold });
    if (cells.length > width) width = cells.length;
  }
  add('title', [result.title + (stamp ? '  (' + stamp + ')' : '')]);
  result.blocks.forEach(function (b) {
    if (b.type === 'table') {
      add('caption', [b.caption || '']);
      var cols = b.columns;
      add('header', cols.map(function (c) { return c.label; }), null, cols.map(exportAlign_));
      b.rows.forEach(function (row) {
        var cells = [], formats = [];
        cols.forEach(function (c) {
          var e = exportCell_(row[c.key], row._fmt || c.format);
          cells.push(e.value);
          formats.push(e.numFmt);
        });
        add(row._kind || 'body', cells, formats, cols.map(exportAlign_), row._bold);
      });
      (b.notes || []).forEach(function (n) { add('note', [n]); });
      add('blank', []);
    } else if (b.type === 'chart') {
      var dim = chartDims && chartDims[chartIdx];
      if (dim) {
        // 그림이 놓일 행을 비워 두고, 그 아래에 내용이 있는 줄을 둔다 (getLastRow 가 그림 아래까지 이어지도록)
        var need = Math.ceil(dim.height / EXPORT_ROW_PX_) + 1;
        images.push({ rowIndex: rows.length, rows: need, chartIndex: chartIdx, width: dim.width, height: dim.height });
        for (var k = 0; k < need; k++) add('image', []);
        add('note', ['(그림) ' + b.caption]);
      } else {
        add('note', ['(차트 "' + b.caption + '"를 그림으로 넣지 못했습니다. 분석 창에서 확인하세요.)']);
      }
      add('blank', []);
      chartIdx++;
    } else {
      add('text', [b.text]);
      add('blank', []);
    }
  });
  (result.warnings || []).forEach(function (w) { add('warn', [w]); });
  add('blank', []);
  add('blank', []);
  return { rows: rows, width: width, images: images };
}

/** 클라이언트에서 받은 결과가 쓸 수 있는 모양인지 확인한다 (아니면 사용자용 메시지로 오류). */
function validateExportResult_(res) {
  var bad = function () { throw new Error('내보낼 결과의 형식이 올바르지 않습니다. 분석을 다시 실행해 주세요.'); };
  if (!res || typeof res !== 'object' || typeof res.title !== 'string' || !Array.isArray(res.blocks)) bad();
  res.blocks.forEach(function (b) {
    if (!b || typeof b !== 'object') bad();
    if (b.type === 'table') {
      if (!Array.isArray(b.columns) || !Array.isArray(b.rows)) bad();
      b.columns.forEach(function (c) { if (!c || typeof c.key !== 'string') bad(); });
    }
  });
}

/**
 * 클라이언트가 보낸 차트 그림을 확인하고 결과별·차트별 배열로 정리한다.
 * images[i][j] = { png(base64), width, height } 또는 null — i: 결과, j: 그 결과의 j번째 차트 블록.
 * images 가 없으면 그림 없이 내보낸다. 모양이 틀리면 사용자용 메시지로 오류.
 */
function validateExportImages_(images, results) {
  var bad = function () { throw new Error('차트 그림의 형식이 올바르지 않습니다. 다시 시도하거나 새로고침한 뒤 내보내 주세요.'); };
  var out = results.map(function () { return []; });
  if (images === undefined || images === null) return out;
  if (!Array.isArray(images) || images.length !== results.length) bad();
  var total = 0;
  images.forEach(function (list, i) {
    var charts = results[i].blocks.filter(function (b) { return b.type === 'chart'; }).length;
    if (!Array.isArray(list) || list.length !== charts) bad();
    list.forEach(function (it, j) {
      if (it === null || it === undefined) { out[i][j] = null; return; }
      if (typeof it !== 'object' || typeof it.png !== 'string' || it.png.length > EXPORT_MAX_PNG_CHARS_ ||
          !/^iVBORw0KGgo[A-Za-z0-9+/]*={0,2}$/.test(it.png)) bad();      // PNG 시그니처로 시작하는 base64 만 허용
      var w = Number(it.width), h = Number(it.height);
      if (!(w >= 40 && w <= 2000 && h >= 40 && h <= 2000)) bad();
      total++;
      if (total > EXPORT_MAX_IMAGES_) throw new Error('한 번에 내보낼 수 있는 차트는 ' + EXPORT_MAX_IMAGES_ + '개까지입니다. 결과를 나누어 내보내 주세요.');
      out[i][j] = { png: it.png, width: Math.round(w), height: Math.round(h) };
    });
  });
  return out;
}

/** JSON 으로 주고받을 때 NaN/Infinity 가 null 로 바뀌지 않도록 "." 문자열로 바꾼 복사본 */
function jsonSafeResult_(value) {
  if (typeof value === 'number') return isFinite(value) ? value : '.';
  if (Array.isArray(value)) return value.map(jsonSafeResult_);
  if (value && typeof value === 'object') {
    var out = {};
    Object.keys(value).forEach(function (k) { out[k] = jsonSafeResult_(value[k]); });
    return out;
  }
  return value;
}

/* ---------- 시트 쓰기 (SpreadsheetApp) ---------- */

var EXPORT_STYLE_ = {
  title: { weight: 'bold', size: 14 },
  caption: { weight: 'bold', size: 11 },
  header: { weight: 'bold', size: 10, background: '#f1f3f4' },
  note: { color: '#5f6368', style: 'italic', size: 9 },
  warn: { color: '#b06000', size: 9 }
};

/**
 * 격자를 sheet 의 startRow 부터 쓰고 서식을 입힌 뒤, 다음에 쓸 행 번호를 돌려준다.
 * blobs[k]: grid.images[k] 에 넣을 PNG Blob (없으면 그 자리는 비워 둔다)
 * report: { placed, errors } — 그림 삽입 결과를 모은다. 한 그림이 실패해도 나머지 내보내기는 계속한다.
 */
function writeGrid_(sheet, grid, startRow, blobs, report) {
  var n = grid.rows.length, w = grid.width, i, j;
  if (sheet.getMaxRows() < startRow + n) sheet.insertRowsAfter(sheet.getMaxRows(), startRow + n - sheet.getMaxRows());
  if (sheet.getMaxColumns() < w) sheet.insertColumnsAfter(sheet.getMaxColumns(), w - sheet.getMaxColumns());

  var values = [], numFmts = [], weights = [], sizes = [], colors = [], styles = [], backgrounds = [], aligns = [];
  for (i = 0; i < n; i++) {
    var row = grid.rows[i], st = EXPORT_STYLE_[row.kind] || {};
    var vr = [], nf = [], wt = [], sz = [], co = [], sy = [], bg = [], al = [];
    for (j = 0; j < w; j++) {
      vr.push(j < row.cells.length ? row.cells[j] : '');
      nf.push(row.formats && row.formats[j] ? row.formats[j] : 'General');
      wt.push(row.bold ? 'bold' : (st.weight || 'normal'));
      sz.push(st.size || 10);
      co.push(st.color || '#202124');
      sy.push(st.style || 'normal');
      bg.push(st.background || null);
      al.push(row.aligns && row.aligns[j] ? row.aligns[j] : 'left');
    }
    values.push(vr); numFmts.push(nf); weights.push(wt); sizes.push(sz); colors.push(co); styles.push(sy); backgrounds.push(bg); aligns.push(al);
  }
  var range = sheet.getRange(startRow, 1, n, w);
  range.setNumberFormats(numFmts);   // 서식을 먼저 정해야 '@' 칸의 문자열이 수식으로 해석되지 않는다
  range.setValues(values);
  range.setFontWeights(weights);
  range.setFontSizes(sizes);
  range.setFontColors(colors);
  range.setFontStyles(styles);
  range.setBackgrounds(backgrounds);
  range.setHorizontalAlignments(aligns);

  // 표의 선: 머리글 위·아래, 마지막 행 아래, 합계 행 위
  var solid = SpreadsheetApp.BorderStyle.SOLID, medium = SpreadsheetApp.BorderStyle.SOLID_MEDIUM;
  for (i = 0; i < n; i++) {
    var r = grid.rows[i], cols = r.cells.length;
    if (r.kind === 'header') {
      var tableEnd = i;
      while (tableEnd + 1 < n && ['body', 'group', 'subtotal', 'total'].indexOf(grid.rows[tableEnd + 1].kind) >= 0) tableEnd++;
      sheet.getRange(startRow + i, 1, 1, cols).setBorder(true, null, true, null, null, null, '#202124', medium);
      sheet.getRange(startRow + tableEnd, 1, 1, cols).setBorder(null, null, true, null, null, null, '#202124', medium);
    } else if (r.kind === 'total' || r.kind === 'subtotal') {
      sheet.getRange(startRow + i, 1, 1, cols).setBorder(true, null, null, null, null, null, '#202124', solid);
    } else if (r.kind === 'group') {
      sheet.getRange(startRow + i, 1, 1, cols).setBorder(true, null, null, null, null, null, '#c4c7c5', solid);
    }
  }
  // 차트 그림: 비워 둔 행의 첫 칸에 놓고 화면과 같은 크기로 맞춘다 (PNG는 2배 해상도)
  (grid.images || []).forEach(function (item, k) {
    if (!blobs || !blobs[k]) return;
    try {
      var img = sheet.insertImage(blobs[k], 1, startRow + item.rowIndex);
      img.setWidth(item.width);
      img.setHeight(item.height);
      if (report) report.placed++;
    } catch (err) {
      var msg = '차트 ' + (item.chartIndex + 1) + ' (' + (startRow + item.rowIndex) + '행): ' + (err && err.message ? err.message : err);
      if (report) report.errors.push(msg);
      console.error('[exportResults] 그림 삽입 실패 ' + msg);
    }
  });
  return startRow + n;
}

/**
 * 대화상자에서 호출: 결과(들)를 "결과" 탭 맨 아래에 이어 쓴다.
 * images: 결과별·차트별 PNG 그림 (validateExportImages_ 참고, 없으면 그림 없이)
 * 반환: { sheet, firstRow, lastRow, images(넣은 그림 수) }
 */
function exportResults(results, images) {
  try {
    return exportResultsUnsafe_(results, images);
  } catch (err) {
    throw exportDeps_().asUserError_(err, 'exportResults');
  }
}

function exportResultsUnsafe_(results, images) {
  if (!Array.isArray(results) || results.length === 0) throw new Error('내보낼 결과가 없습니다.');
  results.forEach(validateExportResult_);
  var imgs = validateExportImages_(images, results);
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var sheet = ss.getSheetByName(SHEET_RESULT);
  var created = false;
  if (!sheet) { sheet = ss.insertSheet(SHEET_RESULT); created = true; }
  if (created) {
    sheet.setColumnWidths(1, 1, 180);
    sheet.setColumnWidths(2, Math.max(1, Math.min(sheet.getMaxColumns() - 1, 20)), 110);
  }
  var row = sheet.getLastRow() === 0 ? 1 : sheet.getLastRow() + 2;
  var first = row;
  var stamp = Utilities.formatDate(new Date(), Session.getScriptTimeZone(), 'yyyy-MM-dd HH:mm:ss');
  var report = { placed: 0, errors: [] };
  var sizes = [];
  results.forEach(function (res, i) {
    var dims = imgs[i].map(function (it) { return it ? { width: it.width, height: it.height } : null; });
    var grid = resultToGrid(res, stamp, dims);
    var blobs = grid.images.map(function (item) {
      var it = imgs[i][item.chartIndex];
      var bytes = Utilities.base64Decode(it.png);
      sizes.push(bytes.length);
      return Utilities.newBlob(bytes, 'image/png', 'chart-' + (i + 1) + '-' + (item.chartIndex + 1) + '.png');
    });
    row = writeGrid_(sheet, grid, row, blobs, report);
  });
  sheet.activate();
  sheet.getRange(first, 1).activate();
  // 시트에 실제로 있는 그림 수 (insertImage 가 오류 없이 끝났는데도 그림이 안 보이는 경우를 가려내는 진단용)
  var onSheet = typeof sheet.getImages === 'function' ? sheet.getImages().length : null;
  console.log('[exportResults] 결과 ' + results.length + '개, 행 ' + first + '~' + (row - 1) + ', 그림 요청 ' + sizes.length +
    '개(PNG 바이트 ' + sizes.join(', ') + '), 삽입 ' + report.placed + '개, 오류 ' + report.errors.length + '개, 시트의 그림 수 ' + onSheet);
  return { sheet: SHEET_RESULT, firstRow: first, lastRow: row - 1, images: report.placed, imageErrors: report.errors, sheetImages: onSheet };
}

/**
 * "결과" 탭의 내용을 지운다: 값·서식·메모와 삽입한 그림. 탭 자체와 열 너비는 남긴다.
 * 반환: { rows, images } (지우기 전의 행 수와 그림 수)
 */
function resetResultSheet_(sheet) {
  var rows = sheet.getLastRow();
  var images = typeof sheet.getImages === 'function' ? sheet.getImages() : [];
  var count = images.length;
  images.forEach(function (img) { img.remove(); });
  sheet.clear();
  return { rows: rows, images: count };
}

/** "거리·유사도"의 "근접도 행렬을 새 탭에 저장"이 만드는 탭 이름 (거리행렬, 유사도행렬, 같은 이름이면 _2, _3 …) */
var MATRIX_TAB_PATTERN_ = /^(거리행렬|유사도행렬)(_\d+)?$/;

/** 이 도구가 만든 근접도 행렬 탭들 */
function matrixTabs_(ss) {
  return typeof ss.getSheets === 'function' ? ss.getSheets().filter(function (sh) { return MATRIX_TAB_PATTERN_.test(sh.getName()); }) : [];
}

/**
 * 메뉴 "결과 탭 지우기": 무엇을 지우는지 안내하고 확인을 받은 뒤 "결과" 탭의 내용을 모두 지우고,
 * 근접도 행렬 탭(거리행렬, 유사도행렬)은 탭째 삭제한다. 다른 탭은 건드리지 않는다.
 */
function clearResultSheet() {
  var ui = SpreadsheetApp.getUi();
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var sheet = ss.getSheetByName(SHEET_RESULT);
  var tabs = matrixTabs_(ss);
  var rows = sheet ? sheet.getLastRow() : 0;
  var images = sheet && typeof sheet.getImages === 'function' ? sheet.getImages().length : 0;
  var hasResult = !!sheet && (rows > 0 || images > 0);
  if (!hasResult && !tabs.length) {
    ui.alert('결과 탭 지우기', sheet ? '"' + SHEET_RESULT + '" 탭이 이미 비어 있습니다. 지울 거리·유사도 행렬 탭도 없습니다.' : '"' + SHEET_RESULT + '" 탭이 없습니다. 지울 내용이 없습니다.', ui.ButtonSet.OK);
    return;
  }
  var lines = [];
  if (hasResult) lines.push('· "' + SHEET_RESULT + '" 탭의 내용 전체 (' + rows + '행, 그림 ' + images + '개)');
  tabs.forEach(function (sh) {
    var r = sh.getLastRow(), c = typeof sh.getLastColumn === 'function' ? sh.getLastColumn() : 0;
    lines.push('· "' + sh.getName() + '" 탭 (거리·유사도 행렬' + (c ? ', ' + r + '행 × ' + c + '열' : '') + ') — 탭 자체를 삭제');
  });
  var answer = ui.alert('결과 탭 지우기',
    '다음을 지웁니다.\n' + lines.join('\n') + '\n\n' +
    (tabs.length ? '행렬 탭은 삭제하면 되돌릴 수 없습니다. 필요하면 먼저 다른 곳에 복사해 두세요.\n' : '') +
    '분석 창에 떠 있는 결과와 다른 탭은 그대로입니다. 계속할까요?', ui.ButtonSet.YES_NO);
  if (answer !== ui.Button.YES) return;
  var done = hasResult ? resetResultSheet_(sheet) : { rows: 0, images: 0 };
  var removed = [];
  tabs.forEach(function (sh) { removed.push(sh.getName()); ss.deleteSheet(sh); });
  var parts = [];
  if (hasResult) parts.push('"' + SHEET_RESULT + '" 탭을 비웠습니다 (' + done.rows + '행, 그림 ' + done.images + '개).');
  if (removed.length) parts.push('행렬 탭 ' + removed.length + '개를 삭제했습니다: ' + removed.join(', ') + '.');
  ui.alert('결과 탭 지우기', parts.join('\n'), ui.ButtonSet.OK);
}

if (typeof module !== 'undefined') {
  module.exports = {
    resetResultSheet_: resetResultSheet_,
    matrixTabs_: matrixTabs_,
    clearResultSheet: clearResultSheet,
    exportCell_: exportCell_,
    resultToGrid: resultToGrid,
    validateExportResult_: validateExportResult_,
    validateExportImages_: validateExportImages_,
    jsonSafeResult_: jsonSafeResult_,
    writeGrid_: writeGrid_,
    exportResults: exportResults
  };
}
