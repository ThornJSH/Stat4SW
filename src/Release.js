/**
 * 배포용 템플릿 점검. 개발 중인 시트를 복사해 배포용 템플릿을 만든 뒤, Apps Script 편집기에서
 * checkTemplateReady 를 실행하면 "비워 두어야 할 것"이 남아 있는지 알려 준다 (메뉴에는 넣지 않는다).
 * evaluateTemplateReadiness_ 는 순수 함수라 Node에서 테스트한다.
 */

var TEMPLATE_VARS_HEADER_ = ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'];

function nonBlank_(v) {
  return !(v === null || v === undefined || (typeof v === 'string' && v.trim() === ''));
}

/**
 * 시트 요약 → 점검 결과.
 * s = { hasData, hasVars, hasResult, dataHeaderCells, dataRows, varHeader[], varRows, otherTabs[] }
 */
function evaluateTemplateReadiness_(s) {
  var problems = [];
  if (!s.hasData) problems.push('"데이터" 탭이 없습니다. (메뉴 Stat4SW 통계분석 > 새 작업 시작으로 만들 수 있습니다)');
  if (!s.hasVars) problems.push('"변수정의" 탭이 없습니다. (메뉴 Stat4SW 통계분석 > 새 작업 시작으로 만들 수 있습니다)');
  if (s.hasVars) {
    var ok = TEMPLATE_VARS_HEADER_.every(function (h, i) { return String(s.varHeader[i] || '').trim() === h; });
    if (!ok) problems.push('"변수정의" 탭 1행이 표준 머리글(' + TEMPLATE_VARS_HEADER_.join(', ') + ')과 다릅니다.');
    if (s.varRows > 0) problems.push('"변수정의" 탭에 변수 정의가 ' + s.varRows + '개 남아 있습니다. 배포용 템플릿은 2행부터 비워 두세요.');
  }
  if (s.hasData) {
    if (s.dataHeaderCells > 0) problems.push('"데이터" 탭 1행에 값이 ' + s.dataHeaderCells + '칸 남아 있습니다. 배포용 템플릿은 비워 두세요.');
    if (s.dataRows > 0) problems.push('"데이터" 탭에 데이터가 ' + s.dataRows + '행 남아 있습니다. 배포용 템플릿은 비워 두세요.');
  }
  if (s.hasResult) problems.push('"결과" 탭이 있습니다. 배포용 템플릿에서는 지우세요 (내보내면 새로 만들어집니다).');
  if (s.otherTabs.length) problems.push('그 밖의 탭이 있습니다: ' + s.otherTabs.join(', ') + '. 의도한 것이 아니면 지우세요.');
  return { ok: problems.length === 0, problems: problems };
}

/** 현재 스프레드시트를 읽어 요약한다 (SpreadsheetApp 사용) */
function collectTemplateSummary_() {
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var data = ss.getSheetByName(SHEET_DATA), vars = ss.getSheetByName(SHEET_VARS);
  var s = {
    hasData: !!data, hasVars: !!vars, hasResult: !!ss.getSheetByName(SHEET_RESULT),
    dataHeaderCells: 0, dataRows: 0, varHeader: [], varRows: 0,
    otherTabs: ss.getSheets().map(function (sh) { return sh.getName(); }).filter(function (n) {
      return n !== SHEET_DATA && n !== SHEET_VARS && n !== SHEET_RESULT;
    })
  };
  if (data) {
    var dv = data.getDataRange().getValues();
    s.dataHeaderCells = (dv[0] || []).filter(nonBlank_).length;
    s.dataRows = dv.slice(1).filter(function (row) { return row.some(nonBlank_); }).length;
  }
  if (vars) {
    var vv = vars.getDataRange().getValues();
    s.varHeader = (vv[0] || []).slice(0, TEMPLATE_VARS_HEADER_.length).map(String);
    s.varRows = vv.slice(1).filter(function (row) { return nonBlank_(row[0]); }).length;
  }
  return s;
}

/** Apps Script 편집기에서 실행: 배포용 템플릿으로 내보내기 전에 비워 둘 것이 남았는지 확인한다. */
function checkTemplateReady() {
  var r = evaluateTemplateReadiness_(collectTemplateSummary_());
  var text = r.ok
    ? '문제가 없습니다. 이 시트를 배포용 템플릿으로 쓸 수 있습니다.'
    : '다음을 확인하세요:\n\n- ' + r.problems.join('\n- ');
  Logger.log(text);
  SpreadsheetApp.getUi().alert('배포용 템플릿 점검 (버전 ' + APP_VERSION + ')', text, SpreadsheetApp.getUi().ButtonSet.OK);
}

/**
 * Apps Script 편집기에서 실행: 실제 Apps Script에서 분석이 얼마나 걸리는지 잰다 (MAX_ROWS 를 정하는 근거).
 * 합성 데이터를 메모리에서 만들어 분석 시간만 재고(시트는 건드리지 않음), 현재 "데이터" 탭을 읽는 시간도 잰다.
 * 결과는 실행 로그와 알림으로 보여 준다. 각 행 수에서 6분 한도의 몇 %를 쓰는지도 함께 적는다.
 */
function benchmarkAnalyses() {
  var seed = 12345;
  function rnd() { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; }
  function gauss() { return (rnd() + rnd() + rnd() + rnd() - 2) * 1.7; }
  function makeDataset(n) {
    var rows = [['g3', 'g2', 'x1', 'x2', 'x3', 'x4', 'x5', 'y']];
    for (var i = 0; i < n; i++) {
      var g3 = 1 + Math.floor(rnd() * 3), g2 = 1 + Math.floor(rnd() * 2);
      var x = [gauss() * 10 + 50, gauss() * 5 + 20, gauss() * 3, gauss() * 8, gauss()];
      rows.push([g3, g2, x[0], x[1], x[2], x[3], x[4], 0.3 * x[0] - 0.2 * x[1] + g3 * 2 + gauss() * 5]);
    }
    return buildDataset(rows, []);
  }
  function time(f) { var t = Date.now(); f(); return Date.now() - t; }

  var lines = ['Apps Script 실제 속도 (밀리초). 6분(360000ms) 한도 대비는 가장 무거운 항목 기준.'];
  [10000, 50000, 100000].forEach(function (n) {
    var ds;
    var build = time(function () { ds = makeDataset(n); });
    var t = {};
    t.desc = time(function () { runDescriptives(ds, { vars: ['x1', 'x2', 'x3'], options: { skewness: true, kurtosis: true, histogram: true } }); });
    t.corr = time(function () { runCorrelation(ds, { vars: ['x1', 'x2', 'x3', 'x4'], options: { pearson: true, spearman: true, scatter: true } }); });
    t.reg = time(function () { runRegression(ds, { dependent: ['y'], covariates: ['x1', 'x2', 'x3', 'x4', 'x5'], factors: ['g3'], options: { collinearity: true } }); });
    t.anova = time(function () { runOneWayAnova(ds, { dependent: ['y'], factor: ['g3'], options: { levene: true, welch: true, tukey: true, boxplot: true } }); });
    var worst = Math.max(t.desc, t.corr, t.reg, t.anova);
    lines.push(n + '행: 구성 ' + build + ' / 기술통계 ' + t.desc + ' / 상관(스피어만 포함) ' + t.corr + ' / 회귀 ' + t.reg + ' / ANOVA+Tukey ' + t.anova +
      '  → 가장 무거운 항목이 한도의 ' + (worst / 3600).toFixed(1) + '%');
  });

  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var data = ss.getSheetByName(SHEET_DATA);
  if (data) {
    var cells = data.getLastRow() * Math.max(1, data.getLastColumn());
    var read = time(function () { data.getDataRange().getValues(); });
    lines.push('현재 "' + SHEET_DATA + '" 탭 읽기: ' + data.getLastRow() + '행 × ' + data.getLastColumn() + '열(' + cells + '칸) ' + read + 'ms');
  }
  lines.push('현재 MAX_ROWS = ' + MAX_ROWS + ' (읽기 시간은 칸 수에 비례하므로 열이 많으면 같은 행 수에서 더 오래 걸린다)');
  var text = lines.join('\n');
  Logger.log(text);
  SpreadsheetApp.getUi().alert('분석 속도 측정 (버전 ' + APP_VERSION + ')', text, SpreadsheetApp.getUi().ButtonSet.OK);
}

if (typeof module !== 'undefined') {
  module.exports = { evaluateTemplateReadiness_: evaluateTemplateReadiness_, nonBlank_: nonBlank_ };
}
