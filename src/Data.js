/**
 * 데이터·변수정의 처리 (순수 함수만 둔다 — SpreadsheetApp 등 GAS API 사용 금지)
 * Node 테스트와 GAS에서 같은 코드를 쓰기 위해 파일 끝에 module.exports 가드를 둔다.
 */

var SCALE_MAP = { '명목': 'nominal', '서열': 'ordinal', '연속': 'scale' };

/** "1=남, 2=여" → { '1': '남', '2': '여' } */
function parseValueLabels_(text) {
  var labels = {};
  if (text === null || text === undefined) return labels;
  String(text).split(/[,;\n]+/).forEach(function (part) {
    var i = part.indexOf('=');
    if (i < 0) return;
    var key = part.slice(0, i).trim();
    var label = part.slice(i + 1).trim();
    if (key !== '') labels[key] = label;
  });
  return labels;
}

/** "99, 999" → [99, 999] (숫자로 읽히면 숫자, 아니면 문자열) */
function parseMissingCodes_(text) {
  if (text === null || text === undefined || text === '') return [];
  if (typeof text === 'number') return [text];
  return String(text).split(/[,;\s]+/).filter(function (s) { return s !== ''; }).map(function (s) {
    var n = Number(s);
    return isNaN(n) ? s : n;
  });
}

/**
 * 변수정의 탭 값(헤더 포함 2차원 배열) → 변수 정의 목록
 * 열 순서: 변수명, 라벨, 척도, 값 레이블, 결측값 코드
 */
function parseVarDefs(rows) {
  var defs = [];
  var warnings = [];
  var seen = {};
  for (var r = 1; r < rows.length; r++) {
    var row = rows[r];
    var name = row[0] === null || row[0] === undefined ? '' : String(row[0]).trim();
    if (name === '') continue;
    if (seen[name]) {
      warnings.push('변수정의: 변수명 "' + name + '"이 중복되어 첫 번째 정의만 사용합니다.');
      continue;
    }
    seen[name] = true;
    var scaleText = row[2] === null || row[2] === undefined ? '' : String(row[2]).trim();
    var scale = SCALE_MAP[scaleText] || null;
    if (scaleText !== '' && !scale) {
      warnings.push('변수정의: "' + name + '"의 척도 "' + scaleText + '"을(를) 알 수 없어 자동 판정합니다. (명목/서열/연속 중 선택)');
    }
    defs.push({
      name: name,
      label: row[1] === null || row[1] === undefined ? '' : String(row[1]).trim(),
      scale: scale,
      valueLabels: parseValueLabels_(row[3]),
      missing: parseMissingCodes_(row[4])
    });
  }
  return { defs: defs, warnings: warnings };
}

function isBlank_(v) {
  return v === null || v === undefined || (typeof v === 'string' && v.trim() === '');
}

function isMissingCode_(v, codes) {
  var s = typeof v === 'string' ? v.trim() : String(v);
  for (var i = 0; i < codes.length; i++) {
    if (String(codes[i]) === s) return true;
  }
  return false;
}

function toNumber_(v) {
  if (typeof v === 'number') return isFinite(v) ? v : null;
  if (typeof v === 'string') {
    var n = Number(v.trim());
    return isNaN(n) ? null : n;
  }
  return null; // Date, boolean 등
}

/**
 * 데이터 탭 값(1행=변수명) + 변수 정의 → 열 단위 데이터셋
 * 결측(빈 셀, 결측값 코드, 연속 변수의 숫자 아님)은 null로 통일한다.
 * 반환: { n, vars:[{name,label,scale,valueLabels,missing,inferred}], data:{name:[...]}, warnings, rowNumbers }
 */
function buildDataset(values, defs) {
  var warnings = [];
  var defByName = {};
  (defs || []).forEach(function (d) { defByName[d.name] = d; });

  var header = values.length ? values[0] : [];
  var cols = []; // { name, index }
  var seen = {};
  header.forEach(function (h, idx) {
    var name = isBlank_(h) ? '' : String(h).trim();
    if (name === '') return;
    if (seen[name]) {
      warnings.push('데이터: 변수명 "' + name + '"이 중복되어 첫 번째 열만 사용합니다.');
      return;
    }
    seen[name] = true;
    cols.push({ name: name, index: idx });
  });

  // 완전히 빈 행은 제외 (rowNumbers: 남은 사례가 시트의 몇 번째 행인지, 1행=변수명이라 2부터)
  var rows = [], rowNumbers = [];
  values.slice(1).forEach(function (row, i) {
    if (cols.some(function (c) { return !isBlank_(row[c.index]); })) { rows.push(row); rowNumbers.push(i + 2); }
  });

  var vars = [];
  var data = {};
  cols.forEach(function (c) {
    var def = defByName[c.name];
    var codes = def ? def.missing : [];
    var raw = rows.map(function (row) {
      var v = row[c.index];
      if (isBlank_(v)) return null;
      if (codes.length && isMissingCode_(v, codes)) return null;
      return typeof v === 'string' ? v.trim() : v;
    });

    var scale = def && def.scale ? def.scale : null;
    var inferred = false;
    if (!scale) {
      inferred = true;
      var allNumeric = raw.every(function (v) { return v === null || toNumber_(v) !== null; });
      scale = allNumeric ? 'scale' : 'nominal';
    }

    var column = raw;
    if (scale === 'scale') {
      var bad = 0;
      column = raw.map(function (v) {
        if (v === null) return null;
        var n = toNumber_(v);
        if (n === null) bad++;
        return n;
      });
      if (bad > 0) {
        warnings.push('데이터: 연속 변수 "' + c.name + '"에서 숫자가 아닌 값 ' + bad + '건을 결측으로 처리했습니다.');
      }
    }

    data[c.name] = column;
    vars.push({
      name: c.name,
      label: def ? def.label : '',
      scale: scale,
      valueLabels: def ? def.valueLabels : {},
      missing: codes,
      inferred: inferred
    });
  });

  (defs || []).forEach(function (d) {
    if (!seen[d.name]) warnings.push('변수정의: "' + d.name + '"에 해당하는 데이터 열이 없습니다.');
  });

  return { n: rows.length, vars: vars, data: data, warnings: warnings, rowNumbers: rowNumbers };
}

/** 지정한 변수가 모두 유효한 행의 인덱스 (목록별 결측 제거) */
function listwiseIndices(dataset, names) {
  var idx = [];
  for (var i = 0; i < dataset.n; i++) {
    var ok = true;
    for (var k = 0; k < names.length; k++) {
      if (dataset.data[names[k]][i] === null) { ok = false; break; }
    }
    if (ok) idx.push(i);
  }
  return idx;
}

/** 목록별 결측 제거 후 지정한 변수의 열만 뽑는다. { n, columns:{name:[...]}, dropped } */
function selectComplete(dataset, names) {
  var idx = listwiseIndices(dataset, names);
  var columns = {};
  names.forEach(function (nm) {
    columns[nm] = idx.map(function (i) { return dataset.data[nm][i]; });
  });
  return { n: idx.length, columns: columns, dropped: dataset.n - idx.length };
}

/**
 * 분석 요청의 변수 이름 목록을 데이터셋의 변수 정보로 바꾼다. 중복은 한 번만 쓴다.
 * 문제가 있으면 사용자에게 그대로 보일 메시지로 Error를 던진다.
 * opts.numeric: true 이면 숫자가 아닌 값이 있는 변수를 거부한다.
 */
function getVars(dataset, names, opts) {
  if (!names || names.length === 0) throw new Error('분석할 변수를 선택하세요.');
  var numeric = !!(opts && opts.numeric);
  var seen = {};
  var out = [];
  names.forEach(function (name) {
    if (seen[name]) return;
    seen[name] = true;
    var found = null;
    for (var i = 0; i < dataset.vars.length; i++) {
      if (dataset.vars[i].name === name) { found = dataset.vars[i]; break; }
    }
    if (!found) throw new Error('변수 "' + name + '"을(를) 찾을 수 없습니다. 새로고침한 뒤 다시 선택하세요.');
    if (numeric && dataset.data[name].some(function (x) { return x !== null && typeof x !== 'number'; })) {
      throw new Error('변수 "' + name + '"에 숫자가 아닌 값이 있어 이 분석에 쓸 수 없습니다.');
    }
    out.push(found);
  });
  return out;
}

/** 숫자는 값 순, 문자는 사전 순, 숫자가 문자보다 앞. {value} 객체를 비교한다. */
function compareValues_(a, b) {
  var ta = typeof a.value, tb = typeof b.value;
  if (ta !== tb) return ta === 'number' ? -1 : 1;
  return a.value < b.value ? -1 : a.value > b.value ? 1 : 0;
}

/**
 * 종속변수와 집단변수가 모두 유효한 사례만 모아 집단(집단변수 값 순)별로 나눈다.
 * 반환: { groups: [{ value, label, values }], used }  (label: 값 레이블이 있으면 레이블, 없으면 값 문자열)
 */
function splitGroups(dataset, depName, factor) {
  var values = dataset.data[depName];
  var keys = dataset.data[factor.name];
  var tally = {}, order = [], used = 0;
  for (var i = 0; i < dataset.n; i++) {
    if (values[i] === null || keys[i] === null) continue;
    used++;
    var key = typeof keys[i] + ':' + keys[i];
    if (!(key in tally)) { tally[key] = { value: keys[i], values: [] }; order.push(key); }
    tally[key].values.push(values[i]);
  }
  var groups = order.map(function (k) { return tally[k]; }).sort(compareValues_);
  groups.forEach(function (g) {
    var label = factor.valueLabels[String(g.value)];
    g.label = label !== undefined ? label : String(g.value);
  });
  return { groups: groups, used: used };
}

if (typeof module !== 'undefined') {
  module.exports = {
    compareValues_: compareValues_,
    splitGroups: splitGroups,
    parseValueLabels_: parseValueLabels_,
    parseMissingCodes_: parseMissingCodes_,
    parseVarDefs: parseVarDefs,
    buildDataset: buildDataset,
    listwiseIndices: listwiseIndices,
    selectComplete: selectComplete,
    getVars: getVars
  };
}
