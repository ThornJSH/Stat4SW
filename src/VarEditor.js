/**
 * 변수 설정 편집기: 변수정의 탭을 직접 타이핑하지 않고 대화상자에서 편집한다.
 * - 변수 목록은 데이터 탭 1행에서 가져오므로 변수명을 따로 입력·일치시킬 필요가 없다
 * - 값 레이블은 데이터에 실제로 있는 값(과 사례 수)을 표로 보여 주고 레이블만 입력받아 "1=남, 2=여" 형식으로 저장한다
 * 순수 함수(formatValueLabels_, formatMissingCodes_, buildEditorData_, findDefRow_, syncPlan_)는 Node에서 테스트한다.
 */

var VAR_SCALES_ = ['명목', '서열', '연속'];
var VAR_EDITOR_MAX_VALUES_ = 50;     // 값 레이블 표에 보여 줄 값의 최대 개수 (넘으면 범주형 변수가 아니라고 보고 표를 만들지 않는다)
var VAR_SCALE_KO_ = { nominal: '명목', ordinal: '서열', scale: '연속' };

/** 이름은 파일마다 고유해야 한다 (GAS는 모든 파일이 한 전역 범위). */
function veDeps_() {
  if (typeof module !== 'undefined') {
    var data = require('./Data.js');
    return { compareValues_: data.compareValues_, parseMissingCodes_: data.parseMissingCodes_, parseVarDefs: data.parseVarDefs };
  }
  return { compareValues_: compareValues_, parseMissingCodes_: parseMissingCodes_, parseVarDefs: parseVarDefs };
}

/* ---------- 순수 함수 ---------- */

/**
 * [{ key, label }] → 변수정의의 "값 레이블" 칸에 넣을 문자열 ("1=남, 2=여").
 * 레이블이 빈 값은 건너뛴다. 값 레이블 파서는 쉼표·세미콜론·줄바꿈으로 항목을, 첫 등호로 값과 레이블을 나누므로
 * 값에 그 문자가 있으면 저장할 수 없어 건너뛰고, 레이블의 쉼표·세미콜론·줄바꿈은 공백으로 바꾼다. 둘 다 warnings 로 알린다.
 */
function formatValueLabels_(pairs) {
  var parts = [], warnings = [];
  pairs.forEach(function (p) {
    var label = String(p.label === null || p.label === undefined ? '' : p.label).trim();
    if (label === '') return;
    var key = String(p.key);
    if (/[,;\r\n=]/.test(key)) {
      warnings.push('값 "' + key + '"에는 쉼표, 세미콜론, 줄바꿈, 등호가 있어 값 레이블로 저장할 수 없습니다.');
      return;
    }
    var clean = label.replace(/[,;\r\n]+/g, ' ').replace(/\s+/g, ' ').trim();
    if (clean !== label) warnings.push('값 "' + key + '"의 레이블에서 쉼표, 세미콜론, 줄바꿈을 공백으로 바꿨습니다.');
    parts.push(key + '=' + clean);
  });
  return { text: parts.join(', '), warnings: warnings };
}

/** 결측값 코드 입력 → 정리한 문자열 ("99 ;999,99" → "99, 999"). 중복은 한 번만. */
function formatMissingCodes_(text) {
  var seen = {}, out = [];
  veDeps_().parseMissingCodes_(text).forEach(function (c) {
    var s = String(c);
    if (!seen[s]) { seen[s] = true; out.push(s); }
  });
  return out.join(', ');
}

/**
 * 편집기에 보낼 한 변수의 정보. ds: loadDataset 결과, name: 변수명.
 * values: 데이터에 있는 값(사례 수)과, 이미 레이블이 있지만 데이터에는 없는 값(사례 수 0)을 합친 목록.
 */
function buildEditorData_(ds, name) {
  var v = null;
  for (var i = 0; i < ds.vars.length; i++) if (ds.vars[i].name === name) { v = ds.vars[i]; break; }
  if (!v) throw new Error('변수 "' + name + '"을(를) 찾을 수 없습니다. 데이터 탭 1행을 확인하세요.');
  var col = ds.data[name], tally = {}, entries = [], missing = 0;
  col.forEach(function (x) {
    if (x === null) { missing++; return; }
    var k = typeof x + ':' + x;
    if (!(k in tally)) { tally[k] = { value: x, count: 0 }; entries.push(tally[k]); }
    tally[k].count++;
  });
  entries.sort(veDeps_().compareValues_);
  var seenKey = {};
  var values = entries.map(function (e) {
    seenKey[String(e.value)] = true;
    var lab = v.valueLabels[String(e.value)];
    return { key: String(e.value), count: e.count, label: lab === undefined ? '' : lab, inData: true };
  });
  Object.keys(v.valueLabels).forEach(function (k) {
    if (!seenKey[k]) values.push({ key: k, count: 0, label: v.valueLabels[k], inData: false });
  });
  var truncated = values.length > VAR_EDITOR_MAX_VALUES_;
  return {
    name: name,
    label: v.label,
    scale: v.inferred ? '' : VAR_SCALE_KO_[v.scale] || '',          // 변수정의에 적힌 척도 ('' = 비워 둠 = 자동 판정)
    autoScale: VAR_SCALE_KO_[v.scale] || '',                        // 지금 분석에 쓰이는 척도 (자동 판정이면 판정 결과)
    inferred: !!v.inferred,
    missingCodes: v.missing.map(String).join(', '),
    n: col.length,
    missing: missing,
    valueCount: entries.length,
    truncated: truncated,
    values: truncated ? [] : values
  };
}

/** 변수정의 시트 값(머리글 포함 2차원 배열)에서 변수명의 행을 찾는다. 없으면 새로 쓸 행 번호. 행 번호는 시트 기준(1부터). */
function findDefRow_(values, name) {
  var lastFilled = 0;
  for (var i = 1; i < values.length; i++) {
    var cell = values[i][0] === null || values[i][0] === undefined ? '' : String(values[i][0]).trim();
    if (cell === '') continue;
    if (cell === name) return { row: i + 1, isNew: false };
    lastFilled = i;
  }
  return { row: lastFilled + 2, isNew: true };
}

/** 데이터 1행의 변수명과 변수정의의 변수명을 비교한다: 추가할 것, 데이터에 없는 정의 */
function syncPlan_(headerNames, defNames) {
  var inDefs = {}, inData = {};
  defNames.forEach(function (n) { inDefs[n] = true; });
  headerNames.forEach(function (n) { inData[n] = true; });
  return {
    toAdd: headerNames.filter(function (n, i) { return !inDefs[n] && headerNames.indexOf(n) === i; }),
    orphans: defNames.filter(function (n) { return !inData[n]; })
  };
}

/** 저장 요청을 검증하고 시트에 쓸 다섯 칸을 만든다: [변수명, 라벨, 척도, 값 레이블, 결측값 코드] + 경고 */
/** existingValueLabels: payload.keepValueLabels 가 true 일 때(값이 너무 많아 표를 만들지 않은 변수) 그대로 유지할 기존 값 레이블 문자열 */
function buildDefRow_(payload, headerNames, existingValueLabels) {
  if (!payload || typeof payload.name !== 'string' || headerNames.indexOf(payload.name) < 0) {
    throw new Error('"' + (payload && payload.name) + '"은(는) 데이터 탭 1행에 없는 변수입니다. 변수 목록을 새로고침하세요.');
  }
  var scale = payload.scale === undefined || payload.scale === null ? '' : String(payload.scale);
  if (scale !== '' && VAR_SCALES_.indexOf(scale) < 0) throw new Error('척도는 명목, 서열, 연속 중에서 고르세요.');
  var label = String(payload.label === undefined || payload.label === null ? '' : payload.label).trim();
  if (label.length > 200) throw new Error('라벨은 200자 이하로 입력하세요.');
  var pairs = Array.isArray(payload.values) ? payload.values : [];
  if (pairs.length > 500) throw new Error('값 레이블이 너무 많습니다.');
  var vl = payload.keepValueLabels === true ? { text: existingValueLabels || '', warnings: [] } : formatValueLabels_(pairs);
  return { cells: [payload.name, label, scale, vl.text, formatMissingCodes_(payload.missing)], warnings: vl.warnings };
}

/* ---------- 시트 읽기·쓰기 (SpreadsheetApp) ---------- */

function veSheets_() {
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  return { ss: ss, data: ss.getSheetByName(SHEET_DATA), vars: ss.getSheetByName(SHEET_VARS) };
}

/** 변수정의 탭이 없으면 표준 머리글로 만든다 */
function veEnsureVarSheet_(s) {
  if (s.vars) return s.vars;
  s.vars = buildVarsSheet_(s.ss);
  return s.vars;
}

function veDefNames_(varSheet) {
  if (!varSheet) return [];
  return veDeps_().parseVarDefs(varSheet.getDataRange().getValues()).defs.map(function (d) { return d.name; });
}

/** 대화상자에서 호출: 변수 목록(정의 여부)과 데이터에 없는 정의 */
function varEditorState() {
  try {
    var ds = loadDataset_();
    var s = veSheets_();
    var defNames = veDefNames_(s.vars);
    var names = ds.vars.map(function (v) { return v.name; });
    if (s.vars) ensureVarSheetValidation_(s.vars, s.data);      // 0.10.0 이전에 만든 시트에도 드롭다운을 건다
    return {
      n: ds.n,
      vars: ds.vars.map(function (v) {
        return { name: v.name, label: v.label, defined: defNames.indexOf(v.name) >= 0 };
      }),
      orphans: syncPlan_(names, defNames).orphans
    };
  } catch (e) {
    throw asUserError_(e, 'varEditorState');
  }
}

/** 대화상자에서 호출: 한 변수의 편집 정보 */
function varEditorLoad(name) {
  try {
    return buildEditorData_(loadDataset_(), name);
  } catch (e) {
    throw asUserError_(e, 'varEditorLoad');
  }
}

/** 대화상자에서 호출: 한 변수의 정의를 변수정의 탭에 저장한다 (그 변수의 행만 쓰고 다른 행·열은 건드리지 않는다) */
function varEditorSave(payload) {
  try {
    var ds = loadDataset_();
    var names = ds.vars.map(function (v) { return v.name; });
    var s = veSheets_();
    var varSheet = veEnsureVarSheet_(s);
    var existing = varSheet.getDataRange().getValues();
    var at = findDefRow_(existing, payload.name);
    var oldLabels = at.isNew ? '' : String(existing[at.row - 1][3] === undefined || existing[at.row - 1][3] === null ? '' : existing[at.row - 1][3]);
    var row = buildDefRow_(payload, names, oldLabels);
    var range = varSheet.getRange(at.row, 1, 1, 5);
    varSheet.getRange(at.row, 4).setNumberFormat('@');       // "1=남" 같은 문자열이 수식·숫자로 해석되지 않도록
    range.setValues([row.cells]);
    ensureVarSheetValidation_(varSheet, s.data);
    return { row: at.row, isNew: at.isNew, warnings: row.warnings, valueLabels: row.cells[3] };
  } catch (e) {
    throw asUserError_(e, 'varEditorSave');
  }
}

/** 대화상자에서 호출: 데이터에 있는데 변수정의에 없는 변수를 모두 추가한다 (변수명만 채우고 나머지는 비워 둔다) */
function varEditorSync() {
  try {
    var ds = loadDataset_();
    var s = veSheets_();
    var varSheet = veEnsureVarSheet_(s);
    var plan = syncPlan_(ds.vars.map(function (v) { return v.name; }), veDefNames_(varSheet));
    if (plan.toAdd.length) {
      var first = findDefRow_(varSheet.getDataRange().getValues(), '\u0000없는이름\u0000').row;
      varSheet.getRange(first, 1, plan.toAdd.length, 1).setValues(plan.toAdd.map(function (n) { return [n]; }));
    }
    ensureVarSheetValidation_(varSheet, s.data);
    return { added: plan.toAdd, orphans: plan.orphans };
  } catch (e) {
    throw asUserError_(e, 'varEditorSync');
  }
}

if (typeof module !== 'undefined') {
  module.exports = {
    formatValueLabels_: formatValueLabels_,
    formatMissingCodes_: formatMissingCodes_,
    buildEditorData_: buildEditorData_,
    findDefRow_: findDefRow_,
    syncPlan_: syncPlan_,
    buildDefRow_: buildDefRow_
  };
}
