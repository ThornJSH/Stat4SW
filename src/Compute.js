/**
 * 변수 계산: 역문항 처리, 문항 합계·평균, z점수. 결과는 데이터 탭의 새 열로 쓰고 변수정의에도 등록한다.
 * 계산은 순수 함수(reverseCode, rowAggregate, zScores, planCompute_)이고 Node에서 테스트한다.
 * 시트 읽기·쓰기(calcState, calcPlan, calcApply)는 대화상자(compute.html)에서 호출한다.
 *
 * 원칙: 원본 열은 바꾸지 않는다. 새 열 이름이 이미 있으면 거부한다. 결측(빈 칸, 결측값 코드, 계산할 수 없는 값)은 새 열에서 빈 칸이다.
 */

var CALC_MAX_VARS_ = 50;       // 한 번에 만들 수 있는 새 열의 최대 개수
var CALC_NAME_MAX_ = 60;       // 새 변수 이름의 최대 길이
var CALC_SCALE_TEXT_ = { ordinal: '서열', scale: '연속' };

/* ---------- 순수 함수 ---------- */

/** 숫자(또는 숫자 문자열)면 그 수, 아니면 null */
function calcNum_(v) {
  if (typeof v === 'number') return isFinite(v) ? v : null;
  if (typeof v === 'string' && v.trim() !== '') {
    var n = Number(v.trim());
    return isNaN(n) ? null : n;
  }
  return null;
}

/**
 * 역코딩: 새 값 = 최솟값 + 최댓값 − 값 (예: 1~5점 척도에서 1↔5, 2↔4, 3은 그대로).
 * 결측(null)은 결측으로 두고, 척도 범위(min~max) 밖의 값은 뒤집을 수 없어 결측으로 바꾸고 센다.
 * 반환: { values, outOfRange } (outOfRange: 범위 밖이라 결측이 된 건수)
 */
function reverseCode(values, min, max) {
  if (!isFinite(min) || !isFinite(max)) throw new Error('척도의 최솟값과 최댓값을 숫자로 입력하세요.');
  if (!(min < max)) throw new Error('척도의 최솟값은 최댓값보다 작아야 합니다.');
  var out = [], outOfRange = 0;
  values.forEach(function (v) {
    var x = v === null ? null : calcNum_(v);
    if (x === null) { out.push(null); return; }
    if (x < min || x > max) { outOfRange++; out.push(null); return; }
    out.push(min + max - x);
  });
  return { values: out, outOfRange: outOfRange };
}

/**
 * 행별 합계 또는 평균. columns: 같은 길이의 열 배열(결측은 null), kind: 'sum' | 'mean',
 * minValid: 유효한 문항이 이만큼 이상일 때만 계산한다 (모자라면 결측). 반환: { values, computed, missing }
 */
function rowAggregate(columns, kind, minValid) {
  if (kind !== 'sum' && kind !== 'mean') throw new Error('알 수 없는 계산입니다: ' + kind);
  var n = columns.length ? columns[0].length : 0;
  var out = [], computed = 0;
  for (var i = 0; i < n; i++) {
    var s = 0, k = 0;
    for (var j = 0; j < columns.length; j++) {
      var x = columns[j][i] === null ? null : calcNum_(columns[j][i]);
      if (x !== null) { s += x; k++; }
    }
    if (k >= minValid && k > 0) { out.push(kind === 'sum' ? s : s / k); computed++; } else out.push(null);
  }
  return { values: out, computed: computed, missing: n - computed };
}

/** z점수 = (값 − 평균) / 표준편차 (표본 표준편차, n−1). 유효한 값이 2개 미만이거나 표준편차가 0이면 모두 결측. 반환: { values, mean, sd, n } */
function zScores(values) {
  var nums = [];
  values.forEach(function (v) { var x = v === null ? null : calcNum_(v); if (x !== null) nums.push(x); });
  var n = nums.length;
  var res = { values: values.map(function () { return null; }), mean: null, sd: null, n: n };
  if (n < 2) return res;
  var sum = 0, i;
  for (i = 0; i < n; i++) sum += nums[i];
  var mean = sum / n, ss = 0;
  // 평균의 반올림 오차를 줄이기 위해 편차의 합으로 한 번 보정한다 (두 번째 패스)
  var corr = 0;
  for (i = 0; i < n; i++) corr += nums[i] - mean;
  mean += corr / n;
  for (i = 0; i < n; i++) ss += (nums[i] - mean) * (nums[i] - mean);
  var sd = Math.sqrt(ss / (n - 1));
  res.mean = mean; res.sd = sd;
  if (!(sd > 0)) return res;
  res.values = values.map(function (v) {
    var x = v === null ? null : calcNum_(v);
    return x === null ? null : (x - mean) / sd;
  });
  return res;
}

/** 새 변수 이름 검사: 앞뒤 공백을 지운 이름을 돌려준다. taken: 이미 쓰는 이름의 집합 */
function calcCheckName_(name, taken) {
  var s = name === null || name === undefined ? '' : String(name).trim();
  if (s === '') throw new Error('새 변수의 이름을 입력하세요.');
  if (s.length > CALC_NAME_MAX_) throw new Error('새 변수의 이름은 ' + CALC_NAME_MAX_ + '자 이하로 입력하세요: "' + s.slice(0, 20) + '…"');
  if (taken[s]) throw new Error('이미 있는 변수 이름입니다: "' + s + '". 다른 이름이나 접미사를 쓰세요.');
  taken[s] = true;
  return s;
}

function calcVarByName_(ds, name) {
  for (var i = 0; i < ds.vars.length; i++) if (ds.vars[i].name === name) return ds.vars[i];
  throw new Error('변수 "' + name + '"을(를) 찾을 수 없습니다. 데이터 탭 1행을 확인하세요.');
}

/** 변수의 열이 모두 숫자인지 확인하고, 아니면 사용자용 오류 */
function calcRequireNumeric_(ds, v) {
  var bad = 0;
  ds.data[v.name].forEach(function (x) { if (x !== null && calcNum_(x) === null) bad++; });
  if (bad > 0) {
    throw new Error('"' + (v.label || v.name) + '"에 숫자가 아닌 값이 ' + bad + '건 있습니다. 결측을 뜻하는 글자(예: ., NA)는 변수정의의 결측값 코드에 적고, 그 밖의 글자는 고치세요.');
  }
}

function calcSummary_(values) {
  var n = 0, sum = 0, min = null, max = null;
  values.forEach(function (v) {
    if (v === null) return;
    n++; sum += v;
    if (min === null || v < min) min = v;
    if (max === null || v > max) max = v;
  });
  return { valid: n, missing: values.length - n, mean: n ? sum / n : null, min: min, max: max };
}

/**
 * 계산 계획: 데이터셋과 요청에서 새 열들을 만든다 (시트는 건드리지 않는다).
 * request: { kind: 'reverse' | 'sum' | 'mean' | 'z', vars: [이름], min, max, suffix, name, minValid }
 * 반환: { columns: [{ name, label, scale, values, ...요약 }], warnings }
 */
function planCompute_(ds, request) {
  var kind = request && request.kind;
  if (['reverse', 'sum', 'mean', 'z'].indexOf(kind) < 0) throw new Error('알 수 없는 계산입니다: ' + kind);
  var names = (request.vars || []).filter(function (n, i, a) { return a.indexOf(n) === i; });
  if (names.length === 0) throw new Error('계산할 변수를 선택하세요.');
  if (names.length > CALC_MAX_VARS_) throw new Error('한 번에 ' + CALC_MAX_VARS_ + '개까지 계산할 수 있습니다.');
  if (ds.n === 0) throw new Error('데이터가 없습니다. 2행부터 데이터를 입력하세요.');
  var vars = names.map(function (n) { return calcVarByName_(ds, n); });
  vars.forEach(function (v) { calcRequireNumeric_(ds, v); });

  var taken = {}, warnings = [], columns = [];
  ds.vars.forEach(function (v) { taken[v.name] = true; });

  if (kind === 'reverse' || kind === 'z') {
    var suffix = request.suffix === undefined || request.suffix === null ? '' : String(request.suffix).trim();
    if (suffix === '') throw new Error('새 변수 이름을 만들 접미사를 입력하세요 (예: ' + (kind === 'reverse' ? '_r' : '_z') + ').');
    var min = Number(request.min), max = Number(request.max);
    vars.forEach(function (v) {
      var label = v.label || v.name, cell;
      var col = ds.data[v.name];
      if (kind === 'reverse') {
        if (v.scale === 'nominal') throw new Error('"' + label + '"은(는) 명목 척도라 역코딩할 수 없습니다. 척도를 서열이나 연속으로 바꾸세요.');
        if (request.min === '' || request.max === '' || request.min === null || request.max === null || request.min === undefined || request.max === undefined) {
          throw new Error('척도의 최솟값과 최댓값을 입력하세요.');
        }
        var r = reverseCode(col, min, max);
        if (r.outOfRange > 0) warnings.push('"' + label + '": 척도 범위(' + min + '~' + max + ') 밖의 값 ' + r.outOfRange + '건은 뒤집을 수 없어 결측으로 두었습니다.');
        cell = { name: calcCheckName_(v.name + suffix, taken), label: label + ' (역)', scale: v.inferred ? '' : (CALC_SCALE_TEXT_[v.scale] || ''), values: r.values };
      } else {
        var z = zScores(col);
        if (z.mean === null) warnings.push('"' + label + '": 유효한 값이 2개 미만이라 z점수를 계산할 수 없습니다.');
        else if (!(z.sd > 0)) warnings.push('"' + label + '": 모든 값이 같아(표준편차 0) z점수를 계산할 수 없습니다.');
        cell = { name: calcCheckName_(v.name + suffix, taken), label: label + ' (z점수)', scale: '연속', values: z.values, fixed: 3 };
      }
      columns.push(cell);
    });
  } else {
    if (vars.length < 2) throw new Error('합계·평균은 변수를 2개 이상 선택하세요.');
    var minValid = request.minValid === undefined || request.minValid === null || request.minValid === '' ? vars.length : Number(request.minValid);
    if (!isFinite(minValid) || minValid < 1 || minValid > vars.length || Math.floor(minValid) !== minValid) {
      throw new Error('"최소 유효 문항 수"는 1 이상 ' + vars.length + ' 이하의 정수여야 합니다.');
    }
    var agg = rowAggregate(names.map(function (n) { return ds.data[n]; }), kind, minValid);
    if (agg.computed === 0) warnings.push('계산된 값이 하나도 없습니다 (유효한 문항이 ' + minValid + '개 이상인 사례가 없습니다).');
    if (kind === 'sum' && minValid < vars.length) warnings.push('결측을 허용한 합계는 응답한 문항 수가 사례마다 달라 비교하기 어렵습니다. 평균을 쓰는 것을 고려하세요.');
    var nm = calcCheckName_(request.name, taken);
    columns.push({
      name: nm,
      label: vars.length + '문항 ' + (kind === 'sum' ? '합계' : '평균'),
      scale: '연속', values: agg.values, fixed: kind === 'mean' ? 3 : null
    });
    if (agg.missing > 0) warnings.push('결측 ' + agg.missing + '건: 유효한 문항이 ' + minValid + '개 미만인 사례는 비워 둡니다.');
  }
  columns.forEach(function (c) {
    var s = calcSummary_(c.values);
    c.valid = s.valid; c.missing = s.missing; c.mean = s.mean; c.min = s.min; c.max = s.max;
  });
  return { columns: columns, warnings: warnings };
}

/** 열 번호(1부터) → A, B, …, AA */
function calcColLetter_(n) {
  var s = '';
  while (n > 0) { var r = (n - 1) % 26; s = String.fromCharCode(65 + r) + s; n = Math.floor((n - 1) / 26); }
  return s;
}

/* ---------- 시트 읽기·쓰기 (SpreadsheetApp) ---------- */

/** 대화상자: 변수 목록과 역코딩 범위에 쓸 정보 */
function calcState() {
  try {
    var ds = loadDataset_();
    return {
      n: ds.n,
      vars: ds.vars.map(function (v) {
        var col = ds.data[v.name], numeric = true, valid = 0, min = null, max = null;
        col.forEach(function (x) {
          if (x === null) return;
          var num = calcNum_(x);
          if (num === null) { numeric = false; return; }
          valid++;
          if (min === null || num < min) min = num;
          if (max === null || num > max) max = num;
        });
        var lmin = null, lmax = null;
        Object.keys(v.valueLabels || {}).forEach(function (k) {
          var num = calcNum_(k);
          if (num === null) return;
          if (lmin === null || num < lmin) lmin = num;
          if (lmax === null || num > lmax) lmax = num;
        });
        return { name: v.name, label: v.label, scale: v.scale, inferred: v.inferred, numeric: numeric && valid > 0, valid: valid, dataMin: min, dataMax: max, labelMin: lmin, labelMax: lmax };
      })
    };
  } catch (e) {
    throw asUserError_(e, 'calcState');
  }
}

function calcPlanSummary_(plan, startColumn) {
  return {
    warnings: plan.warnings,
    startColumn: startColumn,
    columns: plan.columns.map(function (c, i) {
      return {
        name: c.name, label: c.label, scale: c.scale, column: calcColLetter_(startColumn + i),
        valid: c.valid, missing: c.missing, mean: c.mean, min: c.min, max: c.max,
        sample: c.values.slice(0, 5)
      };
    })
  };
}

/** 대화상자: 미리보기. 시트는 건드리지 않고 만들 열의 이름, 위치, 요약, 앞의 몇 값을 돌려준다. */
function calcPlan(request) {
  try {
    var ds = loadDataset_();
    var plan = planCompute_(ds, request);
    var sheet = SpreadsheetApp.getActiveSpreadsheet().getSheetByName(SHEET_DATA);
    var width = sheet.getDataRange().getValues()[0].length;
    return calcPlanSummary_(plan, width + 1);
  } catch (e) {
    throw asUserError_(e, 'calcPlan');
  }
}

/**
 * 새 열들을 데이터 탭 맨 오른쪽에 쓰고 변수정의에 등록한다 (변수 계산과 분석 결과 저장이 함께 쓴다).
 * columns: [{ name, label, scale, values (데이터셋의 사례 순서, 결측은 null), fixed (소수 자릿수, 선택), valueLabels (변수정의에 넣을 "값=레이블, ..." 문자열, 선택) }]
 * 이름 검사는 호출하는 쪽이 한다. 반환: { width (쓰기 전 열 수), rows, firstColumn, lastColumn }
 */
function writeNewColumns_(ds, columns) {
  var s = veSheets_();
  var all = s.data.getDataRange().getValues();
  var width = all[0].length, rowCount = all.length;
  var maxRow = 1;
  ds.rowNumbers.forEach(function (rn) { if (rn > maxRow) maxRow = rn; });
  if (maxRow > rowCount) throw new Error('데이터가 읽는 중에 바뀐 것 같습니다. 다시 시도하세요.');

  var k = columns.length;
  var matrix = [];
  for (var r = 0; r < rowCount - 1; r++) matrix.push(new Array(k).fill(''));
  columns.forEach(function (c, j) {
    c.values.forEach(function (v, i) { matrix[ds.rowNumbers[i] - 2][j] = v === null ? '' : v; });
  });

  if (s.data.getMaxColumns() < width + k) s.data.insertColumnsAfter(s.data.getMaxColumns(), width + k - s.data.getMaxColumns());
  s.data.getRange(1, width + 1, 1, k).setValues([columns.map(function (c) { return c.name; })]).setFontWeight('bold').setBackground('#e8eaed');
  s.data.getRange(2, width + 1, rowCount - 1, k).setValues(matrix);
  columns.forEach(function (c, j) {
    if (c.fixed) s.data.getRange(2, width + 1 + j, rowCount - 1, 1).setNumberFormat('0.' + new Array(c.fixed + 1).join('0'));
  });

  // 변수정의에 등록 (마지막 정의 다음 행부터)
  var varSheet = veEnsureVarSheet_(s);
  var first = findDefRow_(varSheet.getDataRange().getValues(), '\u0000없는이름\u0000').row;
  varSheet.getRange(first, 4, k, 1).setNumberFormat('@');            // "0=불합격" 같은 값 레이블이 수식·숫자로 해석되지 않도록
  varSheet.getRange(first, 1, k, 5).setValues(columns.map(function (c) { return [c.name, c.label, c.scale, c.valueLabels || '', '']; }));
  ensureVarSheetValidation_(varSheet, s.data);
  return { width: width, rows: rowCount - 1, firstColumn: calcColLetter_(width + 1), lastColumn: calcColLetter_(width + k) };
}

/** 이미 있는 이름이면 _2, _3 … 을 붙여 겹치지 않는 이름을 만든다. taken: 이름 → true */
function uniqueNewName_(base, taken) {
  var name = base, i = 2;
  while (taken[name]) { name = base + '_' + i; i++; }
  taken[name] = true;
  return name;
}

/** 대화상자: 데이터 탭 맨 오른쪽에 새 열을 쓰고 변수정의에 등록한다 */
function calcApply(request) {
  try {
    var ds = loadDataset_();
    var plan = planCompute_(ds, request);
    var w = writeNewColumns_(ds, plan.columns);
    var summary = calcPlanSummary_(plan, w.width + 1);
    summary.added = plan.columns.map(function (c) { return c.name; });
    summary.rows = w.rows;
    summary.firstColumn = w.firstColumn;
    summary.lastColumn = w.lastColumn;
    return summary;
  } catch (e) {
    throw asUserError_(e, 'calcApply');
  }
}

if (typeof module !== 'undefined') {
  module.exports = {
    calcNum_: calcNum_,
    reverseCode: reverseCode,
    rowAggregate: rowAggregate,
    zScores: zScores,
    planCompute_: planCompute_,
    calcColLetter_: calcColLetter_,
    uniqueNewName_: uniqueNewName_,
    calcCheckName_: calcCheckName_
  };
}
