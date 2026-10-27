const test = require('node:test');
const assert = require('node:assert/strict');
const Data = require('../src/Data.js');
const Desc = require('../src/Stats_Descriptive.js');
const Rel = require('../src/Stats_Relation.js');
const Reg = require('../src/Stats_Regression.js');
const An = require('../src/Stats_Anova.js');
const Ex = require('../src/Export.js');
const Ch = require('../src/Charts.js');

const ds = (values, defs) => Data.buildDataset(values, defs || []);

/* ---------- 셀 변환 ---------- */

test('exportCell_: 숫자 서식, p값 문자열, 비어 있는 값, 계산 불가, 수식 방지', () => {
  assert.deepEqual(Ex.exportCell_(7.6, { type: 'int' }), { value: 8, numFmt: '0' });
  assert.deepEqual(Ex.exportCell_(3.14159, { type: 'fixed', digits: 3 }), { value: 3.14159, numFmt: '0.000' });
  assert.deepEqual(Ex.exportCell_(58.3333, { type: 'fixed', digits: 1 }), { value: 58.3333, numFmt: '0.0' });
  assert.deepEqual(Ex.exportCell_(2.5004, { type: 'trim', digits: 3 }), { value: 2.5, numFmt: null });
  assert.deepEqual(Ex.exportCell_(0.0432, { type: 'pvalue' }), { value: '.043', numFmt: '@' });
  assert.deepEqual(Ex.exportCell_(0.0001, { type: 'pvalue' }), { value: '<.001', numFmt: '@' });
  assert.deepEqual(Ex.exportCell_(null, { type: 'fixed' }), { value: '', numFmt: null });
  assert.deepEqual(Ex.exportCell_(undefined, null), { value: '', numFmt: null });
  assert.deepEqual(Ex.exportCell_(NaN, { type: 'fixed' }), { value: '.', numFmt: null });
  assert.deepEqual(Ex.exportCell_(Infinity, { type: 'int' }), { value: '.', numFmt: null });
  assert.deepEqual(Ex.exportCell_('남', null), { value: '남', numFmt: null });
  for (const s of ['=SUM(A1:A9)', '+1+1', '-cmd', '@x']) {
    assert.deepEqual(Ex.exportCell_(s, null), { value: s, numFmt: '@' }, s);   // 텍스트 서식으로 수식 해석을 막는다
  }
  assert.deepEqual(Ex.exportCell_('정상 텍스트', null).numFmt, null);
});

/* ---------- 격자 ---------- */

function sampleResult() {
  const defs = Data.parseVarDefs([['h'], ['g', '집단', '명목', '1=대조, 2=처치', '']]).defs;
  const rows = [['x', 'y', 'g']];
  for (let i = 0; i < 24; i++) rows.push([i + (i % 5) * 0.5, (i * 7) % 13 + i * 0.2, 1 + (i % 2)]);
  return { dataset: ds(rows, defs) };
}

test('resultToGrid: 제목, 표, 각주, 경고, 빈 줄 구성과 행 종류', () => {
  const { dataset } = sampleResult();
  const res = Desc.runDescriptives(dataset, { vars: ['x', 'y'], options: { histogram: true } });
  res.warnings.push('테스트 경고');
  const g = Ex.resultToGrid(res, '2026-10-21 12:00:00');
  assert.equal(g.rows[0].kind, 'title');
  assert.equal(g.rows[0].cells[0], '기술통계  (2026-10-21 12:00:00)');
  assert.equal(g.rows[1].kind, 'caption');
  assert.equal(g.rows[2].kind, 'header');
  assert.deepEqual(g.rows[2].cells, ['변수', 'N', '최솟값', '최댓값', '평균', '표준편차']);
  assert.deepEqual(g.rows[2].aligns, ['left', 'right', 'right', 'right', 'right', 'right']);
  assert.equal(g.rows[3].kind, 'body');
  assert.equal(g.rows[3].cells[0], 'x');
  assert.equal(g.rows[3].formats[1], '0');         // N: 정수
  assert.equal(g.rows[3].formats[4], '0.000');     // 평균
  const total = g.rows.find(r => r.kind === 'total');
  assert.equal(total.cells[0], '유효 N(목록별)');
  assert.ok(g.rows.some(r => r.kind === 'note' && r.cells[0].includes('변수별로 유효한 사례')));
  assert.ok(g.rows.some(r => r.kind === 'note' && r.cells[0].includes('그림으로 넣지 못했습니다')));   // 히스토그램 (그림이 없으면 안내 줄)
  assert.ok(g.rows.some(r => r.kind === 'warn' && r.cells[0] === '테스트 경고'));
  assert.deepEqual(g.rows.slice(-2).map(r => r.kind), ['blank', 'blank']);
  assert.equal(g.width, 6);
});

test('resultToGrid: 행 단위 서식(_fmt), p값 문자열, 상관 행렬의 계산 불가(.)', () => {
  const d = ds([['a', 'b', 'c'], [1, 2, 5], [2, 4, 5], [3, 5, 5], [4, 9, 5]]);
  const res = Rel.runCorrelation(d, { vars: ['a', 'b', 'c'], options: { pearson: true } });
  const g = Ex.resultToGrid(res);
  const body = g.rows.filter(r => r.kind === 'body' || r.kind === 'group');
  assert.equal(body[0].cells[1], 'Pearson 상관계수');
  assert.equal(body[0].formats[2], '0.000');
  assert.equal(body[1].cells[1], '유의확률(양측)');
  assert.equal(typeof body[1].cells[3], 'string');     // p값은 문자열
  assert.equal(body[1].formats[3], '@');
  assert.equal(body[2].formats[2], '0');               // N 행은 정수 서식
  assert.equal(body[0].cells[4], '.');                 // b와 c(상수)의 상관은 계산 불가
});

test('resultToGrid: 값 레이블에 수식처럼 보이는 문자열이 있어도 텍스트 서식', () => {
  const defs = Data.parseVarDefs([['h'], ['g', '집단', '명목', '1==SUM(1), 2=-5', '']]).defs;
  const d = ds([['g'], [1], [2], [2]], defs);
  const g = Ex.resultToGrid(Desc.runFrequencies(d, { vars: ['g'], options: {} }));
  const labelRows = g.rows.filter(r => r.kind === 'body' && ['=SUM(1)', '-5'].includes(r.cells[1]));
  assert.equal(labelRows.length, 2);
  labelRows.forEach(r => assert.equal(r.formats[1], '@'));
});

test('jsonSafeResult_: NaN, Infinity 는 "."로 바뀌고 JSON 왕복 후에도 같다', () => {
  const safe = Ex.jsonSafeResult_({ a: NaN, b: [1, Infinity, { c: -Infinity, d: 'x' }], e: null, f: 2.5 });
  assert.deepEqual(safe, { a: '.', b: [1, '.', { c: '.', d: 'x' }], e: null, f: 2.5 });
  assert.deepEqual(JSON.parse(JSON.stringify(safe)), safe);
  const { dataset } = sampleResult();
  const res = Rel.runCorrelation(ds([['a', 'b'], [1, 5], [2, 5], [3, 5]]), { vars: ['a', 'b'], options: { pearson: true } });
  const rt = JSON.parse(JSON.stringify(Ex.jsonSafeResult_(res)));
  const g = Ex.resultToGrid(rt);
  assert.ok(g.rows.some(r => r.cells.includes('.')));
  assert.ok(dataset.n === 24);
});

test('validateExportResult_: 모양이 틀린 결과는 사용자용 메시지로 거부', () => {
  const ok = { title: 't', blocks: [{ type: 'table', columns: [{ key: 'a', label: 'A' }], rows: [{ a: 1 }] }, { type: 'text', text: 'x' }] };
  assert.doesNotThrow(() => Ex.validateExportResult_(ok));
  for (const bad of [null, 'x', {}, { title: 1, blocks: [] }, { title: 't', blocks: 'x' }, { title: 't', blocks: [null] },
    { title: 't', blocks: [{ type: 'table', columns: 'x', rows: [] }] }, { title: 't', blocks: [{ type: 'table', columns: [{ label: 'x' }], rows: [] }] }]) {
    assert.throws(() => Ex.validateExportResult_(bad), /형식이 올바르지 않습니다/);
  }
});

/* ---------- 시트 쓰기: 실제 API의 메서드 이름과 범위 크기 규칙을 따르는 기록용 mock ---------- */

function makeSheetMock(initialRows) {
  const calls = [];
  const state = { lastRow: initialRows || 0, maxRows: 1000, maxCols: 26, values: {}, activated: null };
  function range(r, c, nr, nc) {
    const dims = (arr, what) => {
      if (!Array.isArray(arr) || arr.length !== nr || arr.some(row => !Array.isArray(row) || row.length !== nc)) {
        throw new Error(`${what}: 범위(${nr}x${nc})와 배열 크기가 다릅니다`);
      }
    };
    const rec = (name, check) => (arg) => { if (check) dims(arg, name); calls.push({ name, r, c, nr, nc }); return rng; };
    const rng = {
      setNumberFormats: rec('setNumberFormats', true),
      setValues: (arr) => { dims(arr, 'setValues'); state.values[r] = arr; calls.push({ name: 'setValues', r, c, nr, nc }); return rng; },
      setFontWeights: rec('setFontWeights', true), setFontSizes: rec('setFontSizes', true), setFontColors: rec('setFontColors', true),
      setFontStyles: rec('setFontStyles', true), setBackgrounds: rec('setBackgrounds', true), setHorizontalAlignments: rec('setHorizontalAlignments', true),
      setBorder: (top, left, bottom, right, vertical, horizontal, color, style) => {
        if (style !== 'SOLID' && style !== 'SOLID_MEDIUM') throw new Error('setBorder: 잘못된 선 종류 ' + style);
        calls.push({ name: 'setBorder', r, c, nr, nc }); return rng;
      },
      activate: () => { state.activated = { r, c }; return rng; }
    };
    return rng;
  }
  const sheet = {
    getLastRow: () => state.lastRow, getMaxRows: () => state.maxRows, getMaxColumns: () => state.maxCols,
    insertRowsAfter: (pos, n) => { state.maxRows += n; calls.push({ name: 'insertRowsAfter', pos, n }); },
    insertColumnsAfter: (pos, n) => { state.maxCols += n; calls.push({ name: 'insertColumnsAfter', pos, n }); },
    getRange: (r, c, nr, nc) => {
      if (r < 1 || c < 1 || (nr !== undefined && nr < 1) || (nc !== undefined && nc < 1)) throw new Error(`getRange 인자 오류 ${r},${c},${nr},${nc}`);
      if (nr !== undefined && r + nr - 1 > state.maxRows) throw new Error('범위가 시트 행 수를 넘습니다');
      if (nc !== undefined && c + nc - 1 > state.maxCols) throw new Error('범위가 시트 열 수를 넘습니다');
      return range(r, c, nr === undefined ? 1 : nr, nc === undefined ? 1 : nc);
    },
    setColumnWidths: (c, n, w) => { calls.push({ name: 'setColumnWidths', c, n, w }); },
    activate: () => { calls.push({ name: 'sheet.activate' }); },
    getImages: () => calls.filter(c => c.name === 'insertImage' && !c.failed),
    insertImage: (blob, col, row) => {
      if (state.failImages) throw new Error(state.failImages);
      if (!blob || blob.contentType !== 'image/png') throw new Error('insertImage: PNG Blob 이 아닙니다');
      if (col < 1 || row < 1) throw new Error('insertImage 위치 오류');
      const call = { name: 'insertImage', col, row, blob };
      calls.push(call);
      const img = { setWidth: (w) => { call.width = w; return img; }, setHeight: (h) => { call.height = h; return img; } };
      return img;
    }
  };
  return { sheet, calls, state };
}

function installSheetApp(mock, existing) {
  const created = [];
  globalThis.SHEET_RESULT = '결과';
  globalThis.SpreadsheetApp = {
    BorderStyle: { SOLID: 'SOLID', SOLID_MEDIUM: 'SOLID_MEDIUM' },
    getActiveSpreadsheet: () => ({
      getSheetByName: (n) => (existing && n === '결과' ? mock.sheet : null),
      insertSheet: (n) => { created.push(n); return mock.sheet; }
    })
  };
  globalThis.Utilities = {
    formatDate: () => '2026-10-21 12:00:00',
    base64Decode: (s) => Buffer.from(s, 'base64'),
    newBlob: (data, contentType, name) => ({ data, contentType, name })
  };
  globalThis.Session = { getScriptTimeZone: () => 'Asia/Seoul' };
  return created;
}

test('exportResults: 새 결과 탭 만들기, 값·서식·선이 범위 크기에 맞게 기록, 결과 사이 간격', () => {
  const { dataset } = sampleResult();
  const a = Desc.runDescriptives(dataset, { vars: ['x', 'y'], options: {} });
  const b = An.runOneWayAnova(dataset, { dependent: ['x'], factor: ['g'], options: { levene: true, effect: true, tukey: true, boxplot: true } });
  const mock = makeSheetMock(0);
  const created = installSheetApp(mock, false);
  const out = Ex.exportResults([a, b].map(Ex.jsonSafeResult_));
  assert.deepEqual(created, ['결과']);
  assert.equal(out.sheet, '결과');
  assert.equal(out.firstRow, 1);
  assert.ok(out.lastRow > 20);
  assert.ok(mock.calls.some(c => c.name === 'setColumnWidths'));          // 새 탭은 열 너비를 정한다
  assert.ok(mock.calls.filter(c => c.name === 'setValues').length === 2);   // 결과마다 한 번
  assert.ok(mock.calls.filter(c => c.name === 'setBorder').length >= 6);
  assert.deepEqual(mock.state.activated, { r: 1, c: 1 });
  // 첫 값 범위는 1행부터, 두 번째 결과는 첫 결과 바로 아래에서 시작 (빈 줄 2개 포함)
  const sv = mock.calls.filter(c => c.name === 'setValues');
  assert.equal(sv[0].r, 1);
  assert.equal(sv[1].r, sv[0].r + sv[0].nr);
  assert.equal(out.lastRow, sv[1].r + sv[1].nr - 1);
  // 제목 행의 값
  assert.ok(String(mock.state.values[1][0][0]).startsWith('기술통계  (2026-10-21'));
});

test('exportResults: 기존 결과 탭이 있으면 마지막 행 아래 한 줄 띄우고 이어 쓴다', () => {
  const { dataset } = sampleResult();
  const a = Desc.runFrequencies(dataset, { vars: ['g'], options: {} });
  const mock = makeSheetMock(30);
  const created = installSheetApp(mock, true);
  const out = Ex.exportResults([a]);
  assert.deepEqual(created, []);
  assert.equal(out.firstRow, 32);
  assert.ok(!mock.calls.some(c => c.name === 'setColumnWidths'));        // 기존 탭의 너비는 건드리지 않는다
  assert.deepEqual(mock.state.activated, { r: 32, c: 1 });
});

test('exportResults: 시트의 행·열이 모자라면 늘린다', () => {
  const { dataset } = sampleResult();
  const a = An.runOneWayAnova(dataset, { dependent: ['x'], factor: ['g'], options: { levene: true, tukey: true } });
  const mock = makeSheetMock(990);                     // 최대 1000행 중 990행까지 사용 중
  mock.state.maxCols = 3;                              // 열도 모자라게
  installSheetApp(mock, true);
  assert.doesNotThrow(() => Ex.exportResults([a]));
  assert.ok(mock.calls.some(c => c.name === 'insertRowsAfter'));
  assert.ok(mock.calls.some(c => c.name === 'insertColumnsAfter'));
});

test('exportResults: 결과가 없거나 형식이 틀리면 오류, 시트는 건드리지 않는다', () => {
  const mock = makeSheetMock(0);
  installSheetApp(mock, false);
  assert.throws(() => Ex.exportResults([]), /내보낼 결과가 없습니다/);
  assert.throws(() => Ex.exportResults('x'), /내보낼 결과가 없습니다/);
  assert.throws(() => Ex.exportResults([{ title: 't' }]), /형식이 올바르지 않습니다/);
  assert.equal(mock.calls.length, 0);
});

test('회귀 결과(열 많음, 행 단위 서식 없음)도 격자로 변환되고 시트 쓰기가 오류 없이 끝난다', () => {
  const rows = [['y', 'a', 'b']];
  for (let i = 0; i < 20; i++) rows.push([3 + 2 * i + (i % 4), i, (i * 5) % 9]);
  const res = Reg.runRegression(ds(rows), { dependent: ['y'], covariates: ['a', 'b'], options: { ci: true, collinearity: true, dw: true } });
  const mock = makeSheetMock(0);
  installSheetApp(mock, false);
  const out = Ex.exportResults([Ex.jsonSafeResult_(res)]);
  assert.ok(out.lastRow > 15);
  const g = Ex.resultToGrid(res);
  assert.ok(g.width >= 7);                               // 공선성 진단 표의 열
});

/* ---------- 차트 그림 ---------- */

// 1x1 PNG (유효한 PNG 시그니처로 시작하는 base64)
const PNG = 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg==';

function chartResult() {
  return {
    title: 't', warnings: [],
    blocks: [
      { type: 'table', caption: '표', columns: [{ key: 'a', label: 'A' }], rows: [{ a: 1 }], notes: [] },
      Ch.barChartBlock('막대', ['a', 'b'], [1, 2], ''),
      Ch.barChartBlock('두 번째', ['a'], [1], '')
    ]
  };
}

test('resultToGrid: 차트 그림 자리(행 확보, 그림 아래 안내 줄), 그림 없는 차트는 안내 줄', () => {
  const g = Ex.resultToGrid(chartResult(), '', [{ width: 480, height: 300 }, null]);
  assert.equal(g.images.length, 1);
  const im = g.images[0];
  assert.equal(im.rows, 16);                         // ceil(300 / 21) + 1
  assert.equal(im.chartIndex, 0);
  assert.deepEqual([im.width, im.height], [480, 300]);
  for (let k = 0; k < im.rows; k++) assert.equal(g.rows[im.rowIndex + k].kind, 'image');
  assert.equal(g.rows[im.rowIndex + im.rows].kind, 'note');
  assert.equal(g.rows[im.rowIndex + im.rows].cells[0], '(그림) 막대');       // 그림 바로 아래의 내용 있는 줄
  assert.ok(g.rows.some(r => r.kind === 'note' && r.cells[0].includes('"두 번째"') && r.cells[0].includes('그림으로 넣지 못했습니다')));
  // 그림 자리 뒤에 내용 있는 행이 반드시 있다: getLastRow 가 그림 아래까지 이어져 다음 내보내기가 겹치지 않는다
  const lastContent = g.rows.map((r, i) => (r.cells.length ? i : -1)).reduce((a, b) => Math.max(a, b), -1);
  assert.ok(lastContent > im.rowIndex + im.rows - 1);
  // chartDims 를 주지 않으면 (이전 방식) 그림 없음
  assert.equal(Ex.resultToGrid(chartResult(), '').images.length, 0);
});

test('validateExportImages_: 정상, 없음, null, 그리고 잘못된 입력', () => {
  const results = [chartResult(), { title: 'x', blocks: [{ type: 'text', text: 't' }] }];
  const ok = [[{ png: PNG, width: 480, height: 300 }, null], []];
  const out = Ex.validateExportImages_(ok, results);
  assert.equal(out[0][0].png, PNG);
  assert.equal(out[0][1], null);
  assert.deepEqual(out[1], []);
  assert.deepEqual(Ex.validateExportImages_(undefined, results), [[], []]);   // 그림 없이 내보내기
  const bad = /차트 그림의 형식이 올바르지 않습니다/;
  const mk = (it, count = 1) => [[it].concat(new Array(count).fill(null)).slice(0, 2), []];
  assert.throws(() => Ex.validateExportImages_('x', results), bad);
  assert.throws(() => Ex.validateExportImages_([[]], results), bad);                               // 결과 수와 다름
  assert.throws(() => Ex.validateExportImages_([[null], []], results), bad);                       // 차트 수와 다름
  assert.throws(() => Ex.validateExportImages_(mk({ png: 'AAAA', width: 480, height: 300 }), results), bad);   // PNG 시그니처 아님
  assert.throws(() => Ex.validateExportImages_(mk({ png: PNG + '<script>', width: 480, height: 300 }), results), bad);
  assert.throws(() => Ex.validateExportImages_(mk({ png: PNG, width: 5, height: 300 }), results), bad);
  assert.throws(() => Ex.validateExportImages_(mk({ png: PNG, width: 480, height: 99999 }), results), bad);
  assert.throws(() => Ex.validateExportImages_(mk({ png: PNG, width: 'a', height: 300 }), results), bad);
  assert.throws(() => Ex.validateExportImages_(mk({ png: 'iVBORw0KGgo' + 'A'.repeat(2000001), width: 480, height: 300 }), results), bad);
  const many = [[]];
  const manyResult = { title: 'm', blocks: Array.from({ length: 31 }, () => Ch.barChartBlock('c', ['a'], [1], '')) };
  many[0] = Array.from({ length: 31 }, () => ({ png: PNG, width: 100, height: 100 }));
  assert.throws(() => Ex.validateExportImages_(many, [manyResult]), /30개까지/);
});

test('exportResults: 그림을 해당 위치에 화면 크기로 삽입하고, 변환하지 못한 차트는 안내 줄만 남긴다', () => {
  const { dataset } = sampleResult();
  const a = An.runOneWayAnova(dataset, { dependent: ['x'], factor: ['g'], options: { boxplot: true } });     // 차트 1개
  const b = Desc.runDescriptives(dataset, { vars: ['x'], options: { histogram: true } });                      // 차트 1개
  const mock = makeSheetMock(0);
  installSheetApp(mock, false);
  const results = [a, b].map(Ex.jsonSafeResult_);
  const out = Ex.exportResults(results, [[{ png: PNG, width: 480, height: 320 }], [null]]);
  const ins = mock.calls.filter(c => c.name === 'insertImage');
  assert.equal(ins.length, 1);
  assert.equal(out.images, 1);
  assert.equal(ins[0].col, 1);
  assert.deepEqual([ins[0].width, ins[0].height], [480, 320]);
  assert.equal(ins[0].blob.contentType, 'image/png');
  assert.equal(ins[0].blob.name, 'chart-1-1.png');
  assert.ok(Buffer.from(ins[0].blob.data).slice(0, 4).toString('hex') === '89504e47');      // 디코딩한 PNG 바이트
  // 그림은 첫 결과의 값 범위 안에서, 표 아래에 놓인다
  const sv = mock.calls.filter(c => c.name === 'setValues');
  assert.ok(ins[0].row > sv[0].r && ins[0].row < sv[1].r);
  const grid = Ex.resultToGrid(results[0], '', [{ width: 480, height: 320 }]);
  assert.equal(ins[0].row, sv[0].r + grid.images[0].rowIndex);
  assert.equal(sv[0].nr, grid.rows.length);
  // 두 번째 결과는 첫 결과 바로 아래(그림 행을 포함해서)에서 시작한다
  assert.equal(sv[1].r, sv[0].r + sv[0].nr);
});

test('exportResults: 그림 형식이 틀리면 시트를 건드리지 않고 오류', () => {
  const { dataset } = sampleResult();
  const a = An.runOneWayAnova(dataset, { dependent: ['x'], factor: ['g'], options: { boxplot: true } });
  const mock = makeSheetMock(0);
  installSheetApp(mock, false);
  assert.throws(() => Ex.exportResults([Ex.jsonSafeResult_(a)], [[{ png: 'bad', width: 480, height: 320 }]]), /차트 그림의 형식/);
  assert.equal(mock.calls.length, 0);
});

test('Render: 차트 슬롯의 순번이 결과 객체의 차트 블록 순서와 같다', () => {
  const R = require('../src/Render.js');
  const html = R.renderResultHtml(chartResult());
  const nums = [...html.matchAll(/class="chart-slot" data-chart="(\d+)"/g)].map(m => Number(m[1]));
  assert.deepEqual(nums, [0, 1]);
  assert.ok(html.indexOf('data-chart="0"') < html.indexOf('data-chart="1"'));
});

test('exportResults: 그림 삽입 결과를 돌려준다 (넣은 수, 오류, 시트의 그림 수)', () => {
  const { dataset } = sampleResult();
  const a = An.runOneWayAnova(dataset, { dependent: ['x'], factor: ['g'], options: { boxplot: true } });
  const mock = makeSheetMock(0);
  installSheetApp(mock, false);
  const out = Ex.exportResults([Ex.jsonSafeResult_(a)], [[{ png: PNG, width: 480, height: 320 }]]);
  assert.equal(out.images, 1);
  assert.deepEqual(out.imageErrors, []);
  assert.equal(out.sheetImages, 1);
});

test('exportResults: 그림 삽입이 실패해도 표는 모두 쓰고, 실패 사유와 위치를 알려 준다', () => {
  const { dataset } = sampleResult();
  const a = An.runOneWayAnova(dataset, { dependent: ['x'], factor: ['g'], options: { boxplot: true } });
  const mock = makeSheetMock(0);
  mock.state.failImages = '이미지를 삽입할 수 없습니다';
  installSheetApp(mock, false);
  const errors = [];
  const original = console.error;
  console.error = (m) => errors.push(String(m));
  let out;
  try { out = Ex.exportResults([Ex.jsonSafeResult_(a)], [[{ png: PNG, width: 480, height: 320 }]]); } finally { console.error = original; }
  assert.equal(out.images, 0);
  assert.equal(out.imageErrors.length, 1);
  assert.ok(out.imageErrors[0].includes('차트 1') && out.imageErrors[0].includes('이미지를 삽입할 수 없습니다') && /\d+행/.test(out.imageErrors[0]));
  assert.equal(out.sheetImages, 0);
  assert.ok(mock.calls.some(c => c.name === 'setValues'));          // 표는 이미 쓴 뒤
  assert.equal(errors.length, 1);                                     // 실행 로그에도 남는다
});

/* ---------- 결과 탭 지우기 (순수 도우미) ---------- */

test('resetResultSheet_: 그림을 모두 제거하고 시트를 비운 뒤, 지우기 전의 행 수와 그림 수를 돌려준다', () => {
  const removed = [];
  let cleared = 0;
  const sheet = {
    getLastRow: () => 33,
    getImages: () => [1, 2, 3].map(i => ({ remove: () => removed.push(i) })),
    clear: () => { cleared++; }
  };
  assert.deepEqual(Ex.resetResultSheet_(sheet), { rows: 33, images: 3 });
  assert.deepEqual(removed, [1, 2, 3]);
  assert.equal(cleared, 1);
  // getImages 가 없는 시트(구버전 환경)도 오류 없이 비운다
  let c2 = 0;
  assert.deepEqual(Ex.resetResultSheet_({ getLastRow: () => 5, clear: () => { c2++; } }), { rows: 5, images: 0 });
  assert.equal(c2, 1);
});

test('결과 탭을 지운 뒤 다시 내보내면 1행부터 쓴다 (그림 행이 남아 겹치지 않는다)', () => {
  const { dataset } = sampleResult();
  const a = Desc.runFrequencies(dataset, { vars: ['g'], options: { barChart: true } });
  // 지우기 전: 마지막 행이 있으면 그 아래로 이어 쓴다
  const before = makeSheetMock(30);
  installSheetApp(before, true);
  assert.equal(Ex.exportResults([Ex.jsonSafeResult_(a)]).firstRow, 32);
  // 지운 뒤(getLastRow = 0): 1행부터
  const after = makeSheetMock(0);
  installSheetApp(after, true);
  assert.equal(Ex.exportResults([Ex.jsonSafeResult_(a)]).firstRow, 1);
});
