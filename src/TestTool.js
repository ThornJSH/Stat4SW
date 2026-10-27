/**
 * 시험 데이터 넣기: 기법별로 파이썬이 독립 계산한 시험 데이터(TestData.js)를 시트의 데이터·변수정의 탭에 넣고
 * 따라 할 안내와 기대 결과를 시험안내 탭에 적는다. 원래 데이터는 처음 한 번 백업 탭으로 복사해 두고 "되돌리기"로 복원한다.
 * 시험 데이터를 여러 번 바꿔 넣어도 백업은 처음 것 하나만 유지한다 (원래 데이터가 시험 데이터로 덮이지 않게).
 */

var TEST_BACKUP_PROP_ = 'STAT4SW_TEST_BACKUP';
var TEST_GUIDE_SHEET_ = '시험안내';

function testFormatNumber_(v) {
  if (typeof v !== 'number') return String(v);
  if (!isFinite(v)) return String(v);
  var s = Number(v.toPrecision(6));
  return String(s);
}

/** 기본 예제를 레지스트리 항목처럼 다룬다 (Sample.js) */
function testSampleEntry_() {
  var s = sampleDataset_();
  return {
    id: 'sample', analysis: '*', title: '기본 예제 (변수 ' + s.header.length + '개, 60명 — 대부분의 분석)',
    description: '모든 분석을 연습할 수 있는 기본 예제입니다. 변수마다 연습할 분석이 있습니다 (사용안내의 "예제 데이터의 변수" 표).',
    header: s.header, rows: s.rows, varRows: s.varRows,
    guide: ['사용 안내 문서의 "예제 데이터의 변수" 표에서 연습할 분석과 변수를 고른다.',
      '예: 독립표본 t검정(불안점수 × 성별), 매개효과(교육시간 → 자기효능감, 동기수준 → 직무성과), 로지스틱 회귀(이직여부, 변수 계산으로 만든 문항 평균과 근무연수)'],
    expected: []
  };
}

function testFindEntry_(id) {
  if (id === 'sample') return testSampleEntry_();
  for (var i = 0; i < TEST_DATASETS_.length; i++) if (TEST_DATASETS_[i].id === id) return TEST_DATASETS_[i];
  throw new Error('알 수 없는 시험 데이터입니다: ' + id);
}

function testBackupInfo_() {
  var raw = PropertiesService.getDocumentProperties().getProperty(TEST_BACKUP_PROP_);
  if (!raw) return null;
  try { return JSON.parse(raw); } catch (e) { return null; }
}

/** 대화상자: 시험 데이터 목록과 백업 상태 */
function listTestDatasets() {
  try {
    var items = [testSampleEntry_()].concat(TEST_DATASETS_).map(function (d) {
      return {
        id: d.id, analysis: d.analysis, analysisLabel: d.analysis === '*' ? '기본 예제' : (TEST_ANALYSIS_LABELS_[d.analysis] || d.analysis),
        title: d.title, description: d.description, n: d.rows.length, vars: d.header, guide: d.guide,
        expected: d.expected.map(function (e) { return { label: e.label, value: testFormatNumber_(e.value), info: !!e.info }; })
      };
    });
    return { items: items, backup: testBackupInfo_() };
  } catch (e) {
    throw asUserError_(e, 'listTestDatasets');
  }
}

/** 처음 한 번만 현재 데이터·변수정의 탭을 백업 탭으로 복사한다. 이미 백업이 있으면 그대로 둔다. */
function testBackupIfNeeded_(ss, props) {
  var info = testBackupInfo_();
  if (info) return info;
  var stamp = Utilities.formatDate(new Date(), Session.getScriptTimeZone(), 'yyyyMMdd_HHmm');
  info = { time: stamp, data: null, vars: null };
  var data = ss.getSheetByName(SHEET_DATA), vars = ss.getSheetByName(SHEET_VARS);
  if (data) {
    info.data = SHEET_DATA + '_백업_' + stamp;
    data.copyTo(ss).setName(info.data).setTabColor('#9aa0a6');
  }
  if (vars) {
    info.vars = SHEET_VARS + '_백업_' + stamp;
    vars.copyTo(ss).setName(info.vars).setTabColor('#9aa0a6');
  }
  props.setProperty(TEST_BACKUP_PROP_, JSON.stringify(info));
  return info;
}

function testBlank_(rows) {
  return rows.map(function (r) { return r.map(function (v) { return v === null || v === undefined ? '' : v; }); });
}

/** 데이터·변수정의 탭의 내용을 시험 데이터로 바꾼다 */
function testWriteData_(ss, entry) {
  createTemplateTabs_(ss);                                   // 탭이 없으면 만든다
  var data = ss.getSheetByName(SHEET_DATA), vars = ss.getSheetByName(SHEET_VARS);
  var rows = [entry.header].concat(testBlank_(entry.rows));
  if (data.getMaxColumns() < entry.header.length) data.insertColumnsAfter(data.getMaxColumns(), entry.header.length - data.getMaxColumns());
  data.clearContents();
  data.getRange(1, 1, rows.length, entry.header.length).setValues(rows);
  data.getRange(1, 1, 1, entry.header.length).setFontWeight('bold').setBackground('#e8eaed');
  data.setFrozenRows(1);
  vars.getRange(2, 1, Math.max(1, vars.getMaxRows() - 1), VARS_HEADER.length).clearContent();
  vars.getRange(2, 4, entry.varRows.length, 1).setNumberFormat('@');
  vars.getRange(2, 1, entry.varRows.length, VARS_HEADER.length).setValues(entry.varRows);
  ensureVarSheetValidation_(vars, data);
}

/** 시험안내 탭: 따라 할 순서와 기대 결과 */
function testWriteGuide_(ss, entry, backup) {
  var sheet = ss.getSheetByName(TEST_GUIDE_SHEET_) || ss.insertSheet(TEST_GUIDE_SHEET_);
  sheet.clearContents();
  var rows = [['시험 데이터: ' + entry.title, '', ''], [entry.description, '', ''], ['', '', ''], ['따라 할 순서', '', ''],
    ['1', '시트를 새로고침하고 분석 창을 연 뒤 변수 목록의 "새로고침"을 누른다.', '']];
  entry.guide.forEach(function (g, i) { rows.push([String(i + 2), g, '']); });
  if (entry.expected.length) {
    rows.push(['', '', ''], ['기대 결과 (파이썬이 독립 계산한 값. 이 도구의 결과와 소수 여섯째 자리까지 같아야 한다)', '', ''], ['항목', '기준값', '비고']);
    entry.expected.forEach(function (e) {
      rows.push([e.label, e.info ? String(e.value) : testFormatNumber_(e.value), e.info ? '참고(표의 문구와 비교)' : (e.abs ? '부호는 무시' : '')]);
    });
  }
  rows.push(['', '', ''], ['되돌리기', '메뉴 Stat4SW 통계분석 > 설정·관리 > 시험 데이터 > "원래 데이터로 되돌리기"를 누르면 데이터·변수정의 탭이 백업으로 복원되고 이 탭은 지워진다.', '']);
  if (backup && (backup.data || backup.vars)) rows.push(['백업 탭', [backup.data, backup.vars].filter(Boolean).join(', '), '']);
  var width = 3;
  sheet.getRange(1, 1, rows.length, width).setValues(rows);
  sheet.getRange(1, 1).setFontWeight('bold').setFontSize(13);
  sheet.setColumnWidth(1, 300);
  sheet.setColumnWidth(2, 560);
  sheet.setColumnWidth(3, 200);
  sheet.getRange(1, 1, rows.length, width).setWrap(true).setVerticalAlignment('top');
  sheet.setTabColor('#fbbc04');
  sheet.activate();
}

/** 대화상자: 시험 데이터를 데이터·변수정의 탭에 넣는다 (원래 데이터는 처음 한 번 백업) */
function insertTestDataset(id) {
  try {
    var entry = testFindEntry_(id);
    var ss = SpreadsheetApp.getActiveSpreadsheet();
    var props = PropertiesService.getDocumentProperties();
    var backup = testBackupIfNeeded_(ss, props);
    testWriteData_(ss, entry);
    testWriteGuide_(ss, entry, backup);
    return {
      id: entry.id, title: entry.title, rows: entry.rows.length, columns: entry.header.length, backup: backup,
      message: '"' + SHEET_DATA + '"와 "' + SHEET_VARS + '" 탭에 시험 데이터를 넣었습니다 (' + entry.rows.length + '행, 변수 ' + entry.header.length + '개). "' + TEST_GUIDE_SHEET_ + '" 탭에 따라 할 순서와 기대 결과가 있습니다.'
    };
  } catch (e) {
    throw asUserError_(e, 'insertTestDataset');
  }
}

/** 대화상자: 시험 데이터를 지우고 처음 백업한 원래 데이터·변수정의 탭으로 되돌린다 */
function restoreTestBackup() {
  try {
    var info = testBackupInfo_();
    if (!info) throw new Error('되돌릴 백업이 없습니다. (시험 데이터를 넣은 적이 없거나 이미 되돌렸습니다.)');
    var ss = SpreadsheetApp.getActiveSpreadsheet();
    var notes = [];
    [[SHEET_DATA, info.data], [SHEET_VARS, info.vars]].forEach(function (pair) {
      var current = ss.getSheetByName(pair[0]), backup = pair[1] ? ss.getSheetByName(pair[1]) : null;
      if (pair[1] && !backup) { notes.push('백업 탭 "' + pair[1] + '"을 찾을 수 없어 "' + pair[0] + '" 탭은 그대로 두었습니다.'); return; }
      if (current) ss.deleteSheet(current);
      if (backup) { backup.setName(pair[0]); backup.setTabColor(null); }
    });
    var guide = ss.getSheetByName(TEST_GUIDE_SHEET_);
    if (guide) ss.deleteSheet(guide);
    PropertiesService.getDocumentProperties().deleteProperty(TEST_BACKUP_PROP_);
    return { message: '원래 데이터로 되돌렸습니다.' + (notes.length ? ' ' + notes.join(' ') : ''), notes: notes };
  } catch (e) {
    throw asUserError_(e, 'restoreTestBackup');
  }
}

if (typeof module !== 'undefined') {
  module.exports = { testFormatNumber_: testFormatNumber_ };
}
