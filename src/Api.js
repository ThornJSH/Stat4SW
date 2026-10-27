/** 시트에서 데이터셋을 읽는다. 탭이 없으면 사용자용 오류 메시지를 던진다. */
function loadDataset_() {
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var dataSheet = ss.getSheetByName(SHEET_DATA);
  if (!dataSheet) throw new Error('"' + SHEET_DATA + '" 탭이 없습니다. 메뉴 Stat4SW 통계분석 > 새 작업 시작 (데이터 모두 지우기)을 먼저 실행하세요.');
  var varSheet = ss.getSheetByName(SHEET_VARS);

  var warnings = [];
  var defs = [];
  if (varSheet) {
    var parsed = parseVarDefs(varSheet.getDataRange().getValues());
    defs = parsed.defs;
    warnings = warnings.concat(parsed.warnings);
  } else {
    warnings.push('"' + SHEET_VARS + '" 탭이 없어 척도를 자동 판정합니다.');
  }

  // 너무 큰 시트는 읽기 전에 막는다 (읽기만으로 시간이 오래 걸리므로). 빈 행이 섞여 있을 수 있어 넉넉하게 4배까지 읽는다.
  var lastRows = dataSheet.getLastRow() - 1;
  if (lastRows > MAX_ROWS * 4) throw new Error(rowLimitMessage_(lastRows, MAX_ROWS, SHEET_DATA));   // 안내에는 실제 한도(MAX_ROWS)를 쓴다
  var ds = buildDataset(dataSheet.getDataRange().getValues(), defs);
  var overLimit = rowLimitMessage_(ds.n, MAX_ROWS, SHEET_DATA);
  if (overLimit) throw new Error(overLimit);
  if (ds.vars.length === 0) {
    throw new Error('"' + SHEET_DATA + '" 탭 1행에 변수명을 입력하고 2행부터 데이터를 넣으세요. (메뉴 Stat4SW 통계분석 > 사용 안내에서 방법을 볼 수 있습니다.)');
  }
  ds.warnings = warnings.concat(ds.warnings);
  return ds;
}

/** 분석 키 → 분석 함수. 파일 로드 순서와 무관하도록 호출 시점에 찾는다. */
function analysisFor_(key) {
  switch (key) {
    case 'frequencies': return runFrequencies;
    case 'descriptives': return runDescriptives;
    case 'normality': return runNormality;
    case 'levene': return runLevene;
    case 'ttest1': return runTTestOne;
    case 'ttestInd': return runTTestIndependent;
    case 'ttestPaired': return runTTestPaired;
    case 'anova': return runOneWayAnova;
    case 'mannWhitney': return runMannWhitney;
    case 'wilcoxonPaired': return runWilcoxonPaired;
    case 'wilcoxonOne': return runWilcoxonOne;
    case 'kruskalWallis': return runKruskalWallis;
    case 'correlation': return runCorrelation;
    case 'crosstab': return runCrosstab;
    case 'reliability': return runReliability;
    case 'regression': return runRegression;
    case 'hierarchical': return runHierarchical;
    case 'efa': return runEfa;
    case 'anova2': return runAnova2;
    case 'kappa': return runKappa;
    case 'ologit': return runOrdinalLogit;
    case 'power': return runPower;
    case 'km': return runKm;
    case 'binom': return runBinomial;
    case 'cochran': return runCochran;
    case 'ipa': return runIpa;
    case 'icc': return runIcc;
    case 'logistic': return runLogistic;
    case 'moderation': return runModeration;
    case 'mediation': return runMediation;
    case 'distance': return runDistance;
    case 'mcnemar': return runMcnemar;
    case 'friedman': return runFriedman;
    case 'rmanova': return runRepeatedAnova;
    case 'psm': return runPsm;
    case 'kmeans': return runKmeans;
    case 'cart': return runCart;
    default: return null;
  }
}

/** 시트의 데이터가 필요 없는 분석 키 (client.html 의 noData 와 같아야 한다) */
function analysisNeedsData_(key) {
  return key !== 'power';
}

/**
 * 대화상자에서 호출: 분석을 실행하고 결과를 HTML로 돌려준다.
 * request = { analysis, vars: [...], options: {...} }
 * 사용자 입력 오류는 Error 메시지로 던져지고 화면에 그대로 보인다.
 */
function runAnalysis(request) {
  try {
    var key = request && request.analysis;
    var run = analysisFor_(key);
    if (!run) throw new Error('알 수 없는 분석입니다: ' + key);
    // 데이터를 쓰지 않는 분석(검정력·표본 크기 계산)은 시트를 읽지 않는다 (데이터 탭이 없어도 실행된다)
    var dataset = analysisNeedsData_(key) ? loadDataset_() : { n: 0, vars: [], data: {}, rowNumbers: [], warnings: [] };
    var result = run(dataset, request);
    // 분석이 새 열(예측값 등)을 돌려주면 데이터 탭에 저장한다 (옵션을 켠 경우에만 분석 함수가 만든다)
    if (result.saveColumns && result.saveColumns.length) {
      var taken = {};
      dataset.vars.forEach(function (v) { taken[v.name] = true; });
      result.saveColumns.forEach(function (c) { c.name = uniqueNewName_(c.name, taken); });
      var written = writeNewColumns_(dataset, result.saveColumns);
      result.blocks.push({ type: 'text', text: '▶ 데이터 탭 ' + written.firstColumn + (written.lastColumn !== written.firstColumn ? '~' + written.lastColumn : '') + '열에 ' +
        result.saveColumns.map(function (c) { return c.name; }).join(', ') + '을(를) 저장하고 변수정의에도 등록했습니다 (분석 창의 새로고침을 누르면 변수 목록에 나옵니다).' });
    }
    delete result.saveColumns;
    // 큰 표(근접도 행렬 등)는 화면에 다 싣지 않고 새 탭에 쓴다 (분석 함수가 옵션을 켠 경우에만 만든다)
    if (result.saveSheet) {
      var tab = writeResultTab_(result.saveSheet.name, result.saveSheet.rows);
      result.blocks.push({ type: 'text', text: '▶ 새 탭 "' + tab.name + '"에 ' + tab.rows + '행 × ' + tab.cols + '열을 저장했습니다 (첫 행과 첫 열은 이름).' });
      delete result.saveSheet;
    }
    var big = largeDataWarning_(dataset.n);
    if (big) result.warnings.unshift(big);
    // result 는 "시트로 내보내기"에 쓰도록 함께 돌려준다 (NaN, Infinity 는 JSON 으로 못 보내서 "." 로 바꾼 복사본)
    return { title: result.title, html: renderResultHtml(result), result: jsonSafeResult_(result) };
  } catch (e) {
    throw asUserError_(e, 'runAnalysis');
  }
}

/** 대화상자용: 변수 목록과 유효/결측 건수 요약 */
function getDatasetInfo() {
  try {
    var info = datasetInfo_(loadDataset_());
    var mode = getModeInfo_();
    info.advanced = mode.advanced;
    info.advancedKeys = mode.advancedKeys;
    return info;
  } catch (e) {
    throw asUserError_(e, 'getDatasetInfo');
  }
}

/**
 * 새 탭을 만들어 2차원 배열을 쓴다. 같은 이름이 있으면 _2, _3 … 을 붙인다. 첫 행·첫 열은 굵게 하고 고정한다.
 * 만든 탭으로 화면이 옮겨 가지 않도록 원래 보던 탭으로 되돌린다.
 */
function writeResultTab_(baseName, rows) {
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var name = baseName, i = 2;
  while (ss.getSheetByName(name)) { name = baseName + '_' + i; i++; }
  var nr = rows.length, nc = rows[0].length;
  var before = null;
  try { before = ss.getActiveSheet(); } catch (e) { before = null; }
  var sh = ss.insertSheet(name);
  if (sh.getMaxRows() < nr) sh.insertRowsAfter(sh.getMaxRows(), nr - sh.getMaxRows());
  if (sh.getMaxColumns() < nc) sh.insertColumnsAfter(sh.getMaxColumns(), nc - sh.getMaxColumns());
  sh.getRange(1, 1, nr, nc).setValues(rows);
  sh.getRange(1, 1, 1, nc).setFontWeight('bold').setBackground('#e8eaed');
  sh.getRange(1, 1, nr, 1).setFontWeight('bold').setBackground('#e8eaed');
  sh.setFrozenRows(1);
  sh.setFrozenColumns(1);
  if (before) { try { ss.setActiveSheet(before); } catch (e) { /* 되돌리지 못해도 결과에는 영향 없음 */ } }
  return { name: name, rows: nr, cols: nc };
}

/** 분석 창에서 "사건으로 볼 값", "기준 범주"를 고를 때 쓰는 변수의 값 목록 (한 변수당 값이 이 수를 넘으면 목록을 주지 않는다) */
var VALUE_SUMMARY_MAX_ = 60;

/**
 * 변수의 값 목록: { name, label, n(유효 사례), numeric(모든 값이 숫자), valueCount, truncated,
 * values: [{ key, label, count }] (값 순서, 데이터에 있는 값만) }. 값 레이블이 없으면 label 은 ''.
 */
function buildValueSummary_(ds, name) {
  var v = null;
  for (var i = 0; i < ds.vars.length; i++) if (ds.vars[i].name === name) { v = ds.vars[i]; break; }
  if (!v) throw new Error('변수 "' + name + '"을(를) 찾을 수 없습니다. 분석 창의 새로고침을 누르세요.');
  var tally = {}, entries = [], valid = 0, numeric = true;
  ds.data[name].forEach(function (x) {
    if (x === null) return;
    valid++;
    if (typeof x !== 'number') numeric = false;
    var k = typeof x + ':' + x;
    if (!(k in tally)) { tally[k] = { value: x, count: 0 }; entries.push(tally[k]); }
    tally[k].count++;
  });
  entries.sort(compareValues_);
  var truncated = entries.length > VALUE_SUMMARY_MAX_;
  return {
    name: name,
    label: v.label,
    n: valid,
    numeric: numeric && valid > 0,
    valueCount: entries.length,
    truncated: truncated,
    values: truncated ? [] : entries.map(function (e) {
      var lab = v.valueLabels[String(e.value)];
      return { key: String(e.value), label: lab === undefined ? '' : lab, count: e.count };
    })
  };
}

/** 대화상자용: 여러 변수의 값 목록 { 변수명: 요약 } */
function getValueSummaries(names) {
  try {
    var ds = loadDataset_(), out = {};
    (names || []).forEach(function (n) { out[n] = buildValueSummary_(ds, n); });
    return out;
  } catch (e) {
    throw asUserError_(e, 'getValueSummaries');
  }
}

function datasetInfo_(ds) {
  return {
    n: ds.n,
    warnings: ds.warnings,
    vars: ds.vars.map(function (v) {
      var missing = ds.data[v.name].filter(function (x) { return x === null; }).length;
      return {
        name: v.name,
        label: v.label,
        scale: v.scale,
        inferred: v.inferred,
        valid: ds.n - missing,
        missing: missing
      };
    })
  };
}
