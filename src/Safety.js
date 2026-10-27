/**
 * 옵션 기억(사용자 속성)과 오류 안내.
 * sanitizeSavedOptions_, asUserError_ 는 순수 함수라 Node에서 테스트한다.
 */

var SETTINGS_KEY_ = 'stat4sw.options';
var SETTINGS_MAX_JSON_ = 8000;     // 속성 값 하나의 한도(9KB)보다 작게
var LARGE_DATA_ROWS_ = 20000;

/**
 * 저장·복원할 옵션을 걸러낸다. 모양: { 분석키: { 옵션키: true|false|숫자|짧은 문자열 } }
 * 모르는 모양이나 값은 버린다. (사용자 속성이 손상되었거나 다른 버전이 쓴 값이어도 안전하게)
 */
function sanitizeSavedOptions_(value) {
  var out = {};
  if (typeof value === 'string') {
    try { value = JSON.parse(value); } catch (e) { return out; }
  }
  if (!value || typeof value !== 'object' || Array.isArray(value)) return out;
  var keyOk = /^[A-Za-z][A-Za-z0-9]{0,30}$/;
  Object.keys(value).slice(0, 30).forEach(function (a) {
    var opts = value[a];
    if (!keyOk.test(a) || !opts || typeof opts !== 'object' || Array.isArray(opts)) return;
    var clean = {};
    Object.keys(opts).slice(0, 30).forEach(function (k) {
      var v = opts[k];
      if (!keyOk.test(k)) return;
      if (typeof v === 'boolean') clean[k] = v;
      else if (typeof v === 'number' && isFinite(v)) clean[k] = v;
      else if (typeof v === 'string' && v.length <= 20) clean[k] = v;
    });
    out[a] = clean;
  });
  return out;
}

/**
 * 프로그램 오류(TypeError 등)를 사용자용 안내로 바꾼다. 사용자가 고칠 수 있는 입력 오류(일반 Error)는 그대로 둔다.
 * 원인은 실행 로그(logger, 기본 console.error)에 남긴다.
 */
function asUserError_(e, where, logger) {
  var internal = e instanceof TypeError || e instanceof ReferenceError || e instanceof RangeError ||
    e instanceof SyntaxError || e instanceof EvalError;
  if (!internal) return e;
  (logger || console.error)('[' + where + '] ' + (e && e.stack ? e.stack : e));
  return new Error('예기치 않은 오류가 발생했습니다 (' + e.name + ': ' + e.message + '). ' +
    '선택한 변수와 옵션, 데이터를 확인하고, 계속되면 이 메시지를 알려 주세요.');
}

/**
 * 사례 수가 한도를 넘으면 사용자용 안내 문구, 아니면 null.
 * (Sheets는 입력 자체를 막을 수 없으므로 분석을 거부하고 어떻게 줄이는지 알려 준다)
 */
function rowLimitMessage_(n, max, tabName) {
  if (!(n > max)) return null;
  return '"' + tabName + '" 탭에 사례가 ' + n + '건 있습니다. 이 도구는 최대 ' + max + '건까지 분석할 수 있습니다. ' +
    '데이터를 나누거나 표본을 추출해 ' + max + '행(사례) 이하로 줄여 주세요.';
}

/** 사례가 많을 때의 경고 문구 (없으면 null) */
function largeDataWarning_(n) {
  return n > LARGE_DATA_ROWS_
    ? '사례가 ' + n + '건으로 많아 계산에 시간이 걸릴 수 있습니다. 이 도구는 수천 건 규모를 기준으로 만들었습니다.'
    : null;
}

/* ---------- 입력할 때의 안내 (단순 트리거) ---------- */

/**
 * 데이터 탭에서 한도를 넘는 행에 입력하면 그 자리에서 안내를 띄운다 (토스트).
 * 단순 트리거라 권한 승인 없이 동작하고, 입력을 막지는 못한다 (막는 것은 분석 시점의 거부).
 * 입력한 범위의 마지막 행이 MAX_ROWS + 1 행(머리글 1행 포함)을 넘으면 알린다.
 */
function onEdit(e) {
  try {
    if (!e || !e.range || !e.source) return;
    if (e.range.getSheet().getName() !== SHEET_DATA) return;
    var lastRow = e.range.getLastRow();
    if (lastRow - 1 > MAX_ROWS) {
      e.source.toast(
        lastRow + '행까지 입력했습니다. 이 도구는 "' + SHEET_DATA + '" 탭의 사례를 최대 ' + MAX_ROWS + '건(' + (MAX_ROWS + 1) + '행)까지 분석합니다. ' +
        '그 이상은 분석이 거부되니 데이터를 나누거나 줄여 주세요.',
        '행 수 제한 (' + MAX_ROWS + '건)', 20);
    }
  } catch (err) {
    // 안내가 실패해도 편집에는 영향이 없어야 한다
  }
}

/* ---------- 진단 정보 (문제를 알릴 때 붙여 넣는 요약) ---------- */

/**
 * 정보 화면의 "진단 정보 복사"가 호출: 버전, 모드, 시트 탭 크기 같은 요약 문자열.
 * 데이터 값이나 변수 이름, 시트 제목은 넣지 않는다 (개인정보가 섞이지 않게).
 */
function getDiagnostics() {
  var tz = 'Asia/Seoul';
  try { tz = Session.getScriptTimeZone() || tz; } catch (e) { /* 기본값 */ }
  var when = '';
  try { when = Utilities.formatDate(new Date(), tz, 'yyyy-MM-dd HH:mm'); } catch (e) { when = new Date().toISOString(); }
  var lines = ['Stat4SW 진단 정보', '버전: ' + APP_VERSION, '시각: ' + when + ' (' + tz + ')', '고급 분석: ' + (isAdvancedMode_() ? '켜짐' : '꺼짐')];
  try {
    var ss = SpreadsheetApp.getActiveSpreadsheet();
    var data = ss.getSheetByName(SHEET_DATA), vars = ss.getSheetByName(SHEET_VARS), res = ss.getSheetByName(SHEET_RESULT);
    var cols = data && typeof data.getLastColumn === 'function' ? data.getLastColumn() : '?';
    lines.push('데이터 탭: ' + (data ? '있음 (사례 ' + Math.max(0, data.getLastRow() - 1) + '행, 변수 ' + cols + '열)' : '없음'));
    lines.push('변수정의 탭: ' + (vars ? '있음 (정의 ' + Math.max(0, vars.getLastRow() - 1) + '개)' : '없음'));
    lines.push('결과 탭: ' + (res ? '있음 (' + res.getLastRow() + '행)' : '없음'));
  } catch (e) {
    lines.push('시트 정보를 읽지 못했습니다.');
  }
  return lines.join('\n');
}

/* ---------- 사용자 속성 (PropertiesService) ---------- */

/** 대화상자에서 호출: 저장된 옵션 (없거나 읽을 수 없으면 빈 객체) */
function getSavedOptions() {
  try {
    return sanitizeSavedOptions_(PropertiesService.getUserProperties().getProperty(SETTINGS_KEY_));
  } catch (e) {
    return {};
  }
}

/** 대화상자에서 호출: 옵션을 저장한다. 저장하지 못해도 분석에는 영향이 없으므로 오류를 던지지 않는다. */
function saveOptions(options) {
  try {
    var json = JSON.stringify(sanitizeSavedOptions_(options));
    if (json.length > SETTINGS_MAX_JSON_) return false;
    PropertiesService.getUserProperties().setProperty(SETTINGS_KEY_, json);
    return true;
  } catch (e) {
    return false;
  }
}

if (typeof module !== 'undefined') {
  module.exports = {
    rowLimitMessage_: rowLimitMessage_,
    sanitizeSavedOptions_: sanitizeSavedOptions_,
    asUserError_: asUserError_,
    largeDataWarning_: largeDataWarning_,
    getSavedOptions: getSavedOptions,
    saveOptions: saveOptions
  };
}
