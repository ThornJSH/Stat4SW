/** "데이터" 탭을 (이미 있으면 모두 지우고) 템플릿 모양으로 만든다. 안내는 셀 값이 아니라 메모로 둔다(값으로 넣으면 변수 이름으로 읽힌다). */
function formatDataSheet_(data, reset) {
  if (reset) data.clear();
  data.setFrozenRows(1);
  data.getRange(1, 1).setNote('1행에 변수명을 입력하고 2행부터 응답을 붙여 넣으세요. 변수 설명은 "변수정의" 탭에 적습니다 (메뉴 Stat4SW 통계분석 > 변수 설정 편집).');
  data.getRange(1, 1, 1, 26).setFontWeight('bold').setBackground('#e8eaed');
  return data;
}

function buildDataSheet_(ss) { return formatDataSheet_(ss.insertSheet(SHEET_DATA)); }

/** "변수정의" 탭을 (이미 있으면 모두 지우고) 표준 머리글과 입력 보조(변수명·척도 드롭다운)로 만든다. */
function formatVarsSheet_(vars, dataSheet, reset) {
  if (reset) vars.clear();
  vars.setFrozenRows(1);
  vars.getRange(1, 1, 1, VARS_HEADER.length)
    .setValues([VARS_HEADER])
    .setFontWeight('bold')
    .setBackground('#e8eaed');
  vars.setColumnWidth(4, 220);
  ensureVarSheetValidation_(vars, dataSheet);
  return vars;
}

function buildVarsSheet_(ss) { return formatVarsSheet_(ss.insertSheet(SHEET_VARS), ss.getSheetByName(SHEET_DATA)); }

/**
 * 변수정의 탭의 입력 보조를 (다시) 건다: 변수명은 데이터 탭 1행을 목록으로 하는 드롭다운(목록에 없으면 경고 표시,
 * 아직 데이터를 붙이기 전에 적을 수도 있어 입력은 막지 않는다), 척도는 명목·서열·연속 드롭다운.
 */
function ensureVarSheetValidation_(varSheet, dataSheet) {
  var rows = 499;
  if (dataSheet) {
    varSheet.getRange(2, 1, rows, 1).setDataValidation(
      SpreadsheetApp.newDataValidation()
        .requireValueInRange(dataSheet.getRange(1, 1, 1, Math.max(1, dataSheet.getMaxColumns())), true)
        .setAllowInvalid(true)
        .setHelpText('데이터 탭 1행의 변수명 중에서 고르세요. 목록에 없는 이름은 경고 표시가 납니다.')
        .build()
    );
  }
  varSheet.getRange(2, 3, rows, 1).setDataValidation(
    SpreadsheetApp.newDataValidation().requireValueInList(['명목', '서열', '연속'], true).setAllowInvalid(false).build()
  );
}

/** 없는 탭(데이터, 변수정의)만 만든다. 이미 있는 탭은 건드리지 않는다. 만든 탭 이름을 돌려준다. */
function createTemplateTabs_(ss) {
  var created = [];
  if (!ss.getSheetByName(SHEET_DATA)) { buildDataSheet_(ss); created.push(SHEET_DATA); }
  if (!ss.getSheetByName(SHEET_VARS)) { buildVarsSheet_(ss); created.push(SHEET_VARS); }
  return created;
}

/**
 * 템플릿을 만들기 전에 지워질 내용을 알리는 확인 문구 (순수 함수).
 * dataInfo: { rows(머리글 제외 사례 행 수), cols(변수 수) } 또는 null(탭 없음), varCount: 변수정의의 정의 수 또는 null(탭 없음)
 */
function templateWarning_(dataInfo, varCount) {
  var lines = [];
  if (dataInfo) lines.push('· "' + SHEET_DATA + '" 탭: 변수 ' + dataInfo.cols + '개, 사례 ' + dataInfo.rows + '행');
  if (varCount !== null) lines.push('· "' + SHEET_VARS + '" 탭: 정의 ' + varCount + '개');
  return '아래 탭의 내용(값, 서식, 메모)을 모두 지우고 빈 템플릿으로 다시 만듭니다.\n' + lines.join('\n') + '\n\n' +
    '· 이 도구에는 되돌리기가 없습니다. 필요하면 먼저 시트를 복사해 두세요 (탭 이름을 우클릭 → 복사).\n' +
    '· 실수로 지웠다면 구글 시트의 "파일 > 버전 기록 > 버전 기록 보기"에서 이전 시점으로 되돌릴 수 있습니다 (시트 전체가 그 시점으로 돌아갑니다).\n' +
    '· 분석 결과로 데이터 탭에 저장한 열(예측값, 군집 번호 등)도 함께 지워집니다.\n' +
    '· "' + SHEET_RESULT + '" 탭과 다른 탭은 그대로입니다.\n\n' +
    '계속할까요?';
}

/** 탭에 내용이 있는지: 데이터는 1행(머리글)부터, 변수정의는 머리글 아래부터 */
function templateTabInfo_(ss) {
  var data = ss.getSheetByName(SHEET_DATA), vars = ss.getSheetByName(SHEET_VARS);
  var dataInfo = data ? { rows: Math.max(0, data.getLastRow() - 1), cols: data.getLastRow() > 0 && typeof data.getLastColumn === 'function' ? data.getLastColumn() : 0, has: data.getLastRow() > 0 } : null;
  var varCount = vars ? Math.max(0, vars.getLastRow() - 1) : null;
  return { data: data, vars: vars, dataInfo: dataInfo, varCount: varCount, hasContent: !!(dataInfo && dataInfo.has) || (varCount !== null && varCount > 0) };
}

/**
 * (메뉴에서는 startNewWork 를 쓴다. 이 함수는 예/아니오 확인창만 쓰는 단순판으로 시트 편집기에서 직접 실행하거나 시험에 쓴다.)
 * 데이터·변수정의 탭이 없으면 만들고, 이미 내용이 있으면 무엇이 지워지는지 안내한 뒤
 * 확인을 받아 두 탭의 내용을 모두 지우고 빈 템플릿으로 다시 만든다 (아니오를 누르면 아무것도 바꾸지 않는다).
 */
function createTemplate() {
  var ui = SpreadsheetApp.getUi();
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var info = templateTabInfo_(ss);
  if (info.hasContent) {
    var answer = ui.alert('템플릿 시트 만들기', templateWarning_(info.dataInfo, info.varCount), ui.ButtonSet.YES_NO);
    if (answer !== ui.Button.YES) return;
  }
  applyTemplate_(ss, info);
}

/** 새 작업을 시작할 때 사용자가 직접 입력해야 하는 확인 문구 */
var NEW_WORK_CONFIRM_ = '지우기';

/** 새 작업 시작 전에 지워질 내용을 알리는 안내 문구 (순수 함수). dataInfo·varCount 는 templateWarning_ 과 같다 */
function newWorkWarning_(dataInfo, varCount) {
  var lines = [];
  if (dataInfo) lines.push('· "' + SHEET_DATA + '" 탭: 변수 ' + dataInfo.cols + '개, 사례 ' + dataInfo.rows + '행');
  if (varCount !== null) lines.push('· "' + SHEET_VARS + '" 탭: 정의 ' + varCount + '개');
  return '새 작업을 시작하려고 아래 탭의 내용(값, 서식, 메모)을 모두 지우고 빈 템플릿으로 다시 만듭니다.\n' + lines.join('\n') + '\n\n' +
    '· 이 도구에는 되돌리기가 없습니다. 필요하면 먼저 시트를 복사해 두세요 (탭 이름을 우클릭 → 복사).\n' +
    '· 실수로 지웠다면 구글 시트의 "파일 > 버전 기록 > 버전 기록 보기"에서 이전 시점으로 되돌릴 수 있습니다 (시트 전체가 그 시점으로 돌아갑니다).\n' +
    '· 분석 결과로 데이터 탭에 저장한 열(예측값, 군집 번호 등)도 함께 지워집니다.\n' +
    '· "' + SHEET_RESULT + '" 탭은 이 단계에서 지우지 않고, 끝나고 따로 물어봅니다.\n\n' +
    '정말 지우려면 아래 칸에 ' + NEW_WORK_CONFIRM_ + ' 라고 입력하고 확인을 누르세요.';
}

/**
 * 메뉴 "새 작업 시작 (데이터 모두 지우기)": 데이터·변수정의 탭을 비우고 빈 템플릿으로 되돌린다. 탭이 없으면 (실수로 지웠거나 새 시트) 만든다.
 * 지워질 내용이 있을 때는 안내 뒤에 "지우기"를 직접 입력해야 진행하고 (취소하거나 다르게 입력하면 아무것도 바꾸지 않는다),
 * 끝나면 이전 분석 결과(결과 탭, 거리·유사도 행렬 탭)도 지울지 따로 묻는다. 고급 분석 모드와 상관없이 누구나 쓸 수 있다.
 */
function startNewWork() {
  var ui = SpreadsheetApp.getUi();
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var info = templateTabInfo_(ss);
  if (info.hasContent) {
    var res = ui.prompt('새 작업 시작 (데이터 모두 지우기)', newWorkWarning_(info.dataInfo, info.varCount), ui.ButtonSet.OK_CANCEL);
    if (res.getSelectedButton() !== ui.Button.OK) return;
    if (String(res.getResponseText()).trim() !== NEW_WORK_CONFIRM_) {
      ui.alert('새 작업 시작', '"' + NEW_WORK_CONFIRM_ + '"를 정확히 입력하지 않아 취소했습니다. 아무것도 바뀌지 않았습니다.', ui.ButtonSet.OK);
      return;
    }
  }
  applyTemplate_(ss, info, '새 작업 시작');
  newWorkClearResults_(ss, ui);
}

/** 새 작업을 시작한 뒤 이전 분석 결과도 지울지 묻고, 예이면 지운다 (결과 탭이 비어 있고 행렬 탭도 없으면 묻지 않는다) */
function newWorkClearResults_(ss, ui) {
  var sheet = ss.getSheetByName(SHEET_RESULT);
  var tabs = matrixTabs_(ss);
  var rows = sheet ? sheet.getLastRow() : 0;
  var images = sheet && typeof sheet.getImages === 'function' ? sheet.getImages().length : 0;
  var hasResult = !!sheet && (rows > 0 || images > 0);
  if (!hasResult && !tabs.length) return;
  var lines = [];
  if (hasResult) lines.push('· "' + SHEET_RESULT + '" 탭의 내용 전체 (' + rows + '행, 그림 ' + images + '개)');
  tabs.forEach(function (sh) { lines.push('· "' + sh.getName() + '" 탭 (거리·유사도 행렬) — 탭 자체를 삭제'); });
  var answer = ui.alert('새 작업 시작', '이전 분석 결과도 함께 지울까요?\n' + lines.join('\n') + '\n\n예 = 지움, 아니오 = 그대로 둠 (나중에 메뉴 "결과 탭 지우기"로 지울 수 있습니다).', ui.ButtonSet.YES_NO);
  if (answer !== ui.Button.YES) return;
  if (hasResult) resetResultSheet_(sheet);
  tabs.forEach(function (sh) { ss.deleteSheet(sh); });
  ui.alert('새 작업 시작', '이전 분석 결과를 지웠습니다. 이제 "' + SHEET_DATA + '" 탭 2행부터 새 데이터를 붙여 넣으세요.', ui.ButtonSet.OK);
}

/** 데이터·변수정의 탭을 (없으면 만들고, 있으면 비워) 템플릿으로 만들고 결과를 알린다 */
function applyTemplate_(ss, info, title) {
  var ui = SpreadsheetApp.getUi();
  var done = [];
  var data = info.data;
  if (data) { formatDataSheet_(data, true); done.push('"' + SHEET_DATA + '" 탭을 비우고 템플릿으로 다시 만듦'); }
  else { data = buildDataSheet_(ss); done.push('"' + SHEET_DATA + '" 탭을 만듦'); }
  if (info.vars) { formatVarsSheet_(info.vars, data, true); done.push('"' + SHEET_VARS + '" 탭을 비우고 템플릿으로 다시 만듦'); }
  else { buildVarsSheet_(ss); done.push('"' + SHEET_VARS + '" 탭을 만듦'); }
  ui.alert(title || '템플릿 시트 만들기', done.join('\n') + '.', ui.ButtonSet.OK);
}

/** 개발·시연용 예제 데이터. 데이터·변수정의 탭 내용을 덮어쓰므로 확인을 받는다. */
function fillSampleData() {
  var ui = SpreadsheetApp.getUi();
  var answer = ui.alert('예제 데이터', '"' + SHEET_DATA + '"와 "' + SHEET_VARS + '" 탭의 내용을 예제로 덮어씁니다. 이 도구에는 되돌리기가 없으니 필요하면 먼저 시트를 복사해 두세요. 실수로 덮어썼다면 "파일 > 버전 기록 > 버전 기록 보기"에서 이전 시점으로 되돌릴 수 있습니다. 계속할까요?', ui.ButtonSet.YES_NO);
  if (answer !== ui.Button.YES) return;

  var ss = SpreadsheetApp.getActiveSpreadsheet();
  createTemplateTabs_(ss);

  var sample = sampleDataset_();                  // Sample.js: 예제 데이터와 변수정의 (60명, 열마다 용도가 있다)
  var varRows = sample.varRows;
  var vars = ss.getSheetByName(SHEET_VARS);
  vars.getRange(2, 1, vars.getMaxRows() - 1, VARS_HEADER.length).clearContent();
  vars.getRange(2, 4, varRows.length, 1).setNumberFormat('@');       // "1=남" 같은 값 레이블이 수식·숫자로 해석되지 않도록
  vars.getRange(2, 1, varRows.length, VARS_HEADER.length).setValues(varRows);

  var rows = [sample.header].concat(sample.rows);
  var data = ss.getSheetByName(SHEET_DATA);
  if (data.getMaxColumns() < rows[0].length) data.insertColumnsAfter(data.getMaxColumns(), rows[0].length - data.getMaxColumns());   // 기본 26열보다 변수가 많다
  data.clearContents();
  data.getRange(1, 1, rows.length, rows[0].length).setValues(rows);
  data.getRange(1, 1, 1, rows[0].length).setFontWeight('bold').setBackground('#e8eaed');
  data.setFrozenRows(1);
  ensureVarSheetValidation_(ss.getSheetByName(SHEET_VARS), data);
}
