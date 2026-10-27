/**
 * 확률분포 구현을 시트 내장 함수(NORM.S.DIST, T.DIST.RT, CHISQ.INV.RT, F.DIST.RT …)와 비교한다.
 * 편집기에서 compareDistributionsWithSheet 를 선택해 실행한다. (메뉴에는 넣지 않는다)
 * 임시 탭 "_분포검증"을 만들었다가 끝나면 지운다. 자유도는 시트 함수가 정수로 자르므로 정수만 쓴다.
 * 정밀 기준은 tests/ 의 scipy 대조 테스트이며, 이 함수는 GAS 환경에서 같은 코드가 돌아가는지 확인하는 용도다.
 */
function compareDistributionsWithSheet() {
  var specs = [
    { n: 'pnorm', a: [[-3], [-1.96], [0], [1], [2.5]], js: function (a) { return pnorm(a[0]); }, f: function (a) { return 'NORM.S.DIST(' + a[0] + ',TRUE)'; } },
    { n: 'qnorm', a: [[0.001], [0.025], [0.5], [0.975]], js: function (a) { return qnorm(a[0]); }, f: function (a) { return 'NORM.S.INV(' + a[0] + ')'; } },
    { n: 't 아래꼬리', a: [[-2, 5], [0.5, 10], [2, 30], [3, 100]], js: function (a) { return pt(a[0], a[1]); }, f: function (a) { return 'T.DIST(' + a[0] + ',' + a[1] + ',TRUE)'; } },
    { n: 't 위꼬리', a: [[0.5, 5], [2, 10], [3, 30], [5, 100]], js: function (a) { return pt(a[0], a[1], false); }, f: function (a) { return 'T.DIST.RT(' + a[0] + ',' + a[1] + ')'; } },
    { n: 't 양측 p값', a: [[2, 5], [2.5, 20], [4, 100]], js: function (a) { return pValueT(a[0], a[1]); }, f: function (a) { return 'T.DIST.2T(' + a[0] + ',' + a[1] + ')'; } },
    { n: 'qt 아래꼬리', a: [[0.025, 5], [0.5, 10], [0.975, 30]], js: function (a) { return qt(a[0], a[1]); }, f: function (a) { return 'T.INV(' + a[0] + ',' + a[1] + ')'; } },
    { n: 'qt 양측', a: [[0.05, 5], [0.01, 20], [0.05, 100]], js: function (a) { return qt(a[0] / 2, a[1], false); }, f: function (a) { return 'T.INV.2T(' + a[0] + ',' + a[1] + ')'; } },
    { n: 'χ² 아래꼬리', a: [[1, 1], [5, 3], [20, 10], [120, 100]], js: function (a) { return pchisq(a[0], a[1]); }, f: function (a) { return 'CHISQ.DIST(' + a[0] + ',' + a[1] + ',TRUE)'; } },
    { n: 'χ² 위꼬리', a: [[3.84, 1], [7.8, 3], [30, 10], [150, 100]], js: function (a) { return pValueChisq(a[0], a[1]); }, f: function (a) { return 'CHISQ.DIST.RT(' + a[0] + ',' + a[1] + ')'; } },
    { n: 'qχ² 아래꼬리', a: [[0.05, 1], [0.5, 5], [0.95, 20]], js: function (a) { return qchisq(a[0], a[1]); }, f: function (a) { return 'CHISQ.INV(' + a[0] + ',' + a[1] + ')'; } },
    { n: 'qχ² 위꼬리', a: [[0.05, 1], [0.01, 5], [0.05, 20]], js: function (a) { return qchisq(a[0], a[1], false); }, f: function (a) { return 'CHISQ.INV.RT(' + a[0] + ',' + a[1] + ')'; } },
    { n: 'F 아래꼬리', a: [[0.5, 2, 5], [1, 5, 10], [3, 10, 30]], js: function (a) { return pf(a[0], a[1], a[2]); }, f: function (a) { return 'F.DIST(' + a[0] + ',' + a[1] + ',' + a[2] + ',TRUE)'; } },
    { n: 'F 위꼬리', a: [[3, 2, 5], [4, 5, 10], [2.5, 10, 100]], js: function (a) { return pValueF(a[0], a[1], a[2]); }, f: function (a) { return 'F.DIST.RT(' + a[0] + ',' + a[1] + ',' + a[2] + ')'; } },
    { n: 'qF 아래꼬리', a: [[0.5, 2, 5], [0.95, 5, 10], [0.99, 10, 30]], js: function (a) { return qf(a[0], a[1], a[2]); }, f: function (a) { return 'F.INV(' + a[0] + ',' + a[1] + ',' + a[2] + ')'; } },
    { n: 'qF 위꼬리', a: [[0.05, 2, 5], [0.05, 5, 10], [0.01, 10, 100]], js: function (a) { return qf(a[0], a[1], a[2], false); }, f: function (a) { return 'F.INV.RT(' + a[0] + ',' + a[1] + ',' + a[2] + ')'; } }
  ];

  var ss = SpreadsheetApp.getActiveSpreadsheet();
  var tmpName = '_분포검증';
  var old = ss.getSheetByName(tmpName);
  if (old) ss.deleteSheet(old);
  var sh = ss.insertSheet(tmpName);

  try {
    var items = [];
    specs.forEach(function (s) {
      s.a.forEach(function (args) { items.push({ spec: s, args: args }); });
    });
    sh.getRange(1, 1, items.length, 1).setFormulas(items.map(function (it) { return ['=' + it.spec.f(it.args)]; }));
    SpreadsheetApp.flush();
    var sheetValues = sh.getRange(1, 1, items.length, 1).getValues();

    var worst = 0, fails = [];
    items.forEach(function (it, i) {
      var expected = sheetValues[i][0];
      var actual = it.spec.js(it.args);
      var label = it.spec.n + '(' + it.args.join(', ') + ')';
      if (typeof expected !== 'number') {
        fails.push(label + ' → 시트 오류: ' + expected);
        return;
      }
      var err = Math.abs(actual - expected) / Math.max(Math.abs(expected), 1e-300);
      worst = Math.max(worst, err);
      if (!(err < 1e-8)) fails.push(label + ' = ' + actual + ', 시트 ' + expected + ', 상대오차 ' + err.toExponential(2));
    });

    var summary = items.length + '건 비교, 최대 상대오차 ' + worst.toExponential(2) +
      (fails.length ? '\n\n불일치 ' + fails.length + '건:\n' + fails.join('\n') : '\n\n모두 일치합니다 (허용 1e-8).');
    Logger.log(summary);
    SpreadsheetApp.getUi().alert('확률분포 ↔ 시트 함수 비교', summary, SpreadsheetApp.getUi().ButtonSet.OK);
  } finally {
    ss.deleteSheet(sh);
  }
}
