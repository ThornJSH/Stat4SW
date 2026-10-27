const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const D = require('../src/Stats_Distance.js');
const dist = require('../src/Distributions.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'distance.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-10, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));

function dataset(X, names) {
  const nm = names || X[0].map((_, j) => 'v' + j);
  return Data.buildDataset([nm].concat(X), nm.map(n => ({ name: n, label: '', scale: 'scale', valueLabels: {}, missing: [] })));
}
function vectorsOf(X, target, std) {
  const cols = X[0].map((_, j) => X.map(r => r[j]));
  const sc = D.standardizeColumns(cols, std, cols.map((_, j) => 'v' + j));
  return target === 'case' ? X.map((_, i) => sc.map(c => c[i])) : sc;
}

/* ---------- scipy 기준값과 대조 ---------- */

test('proximityMatrix·standardizeColumns: scipy.spatial.distance 기준 행렬 39개와 1e-10 이내로 일치', () => {
  assert.ok(ref.pairs.length >= 35);
  for (const c of ref.pairs) {
    const vectors = vectorsOf(c.data, c.target, c.standardize);
    const opts = { p: c.p };
    if (c.measure === 'mahalanobis') {
      const cv = D.standardizeColumns(c.data[0].map((_, j) => c.data.map(r => r[j])), c.standardize, c.data[0].map((_, j) => 'v' + j));
      const k = cv.length, n = cv[0].length, mean = cv.map(x => x.reduce((s, v) => s + v, 0) / n);
      const S = Array.from({ length: k }, (_, i) => Array.from({ length: k }, (_, j) => { let s = 0; for (let t = 0; t < n; t++) s += (cv[i][t] - mean[i]) * (cv[j][t] - mean[j]); return s / (n - 1); }));
      const res = D.runDistance(dataset(c.data), { vars: c.data[0].map((_, j) => 'v' + j), options: { measure: 'mahalanobis', standardize: c.standardize } });
      const shown = table(res, '마할라노비스 거리 행렬');
      c.expected.forEach((row, i) => row.forEach((e, j) => assert.ok(near(shown.rows[i]['c' + j], e, 1e-9), `${c.name} [${i}][${j}] ${shown.rows[i]['c' + j]} vs ${e}`)));
      continue;
    }
    const M = D.proximityMatrix(vectors, c.measure, opts);
    assert.equal(M.length, c.expected.length, c.name);
    c.expected.forEach((row, i) => row.forEach((e, j) => assert.ok(near(M[i][j], e, 1e-10), `${c.name} [${i}][${j}] ${M[i][j]} vs ${e}`)));
  }
});

test('mahalanobisD2: numpy 기준값(D², χ² 유의확률, 임계값, 이상치 후보)과 같다', () => {
  const m = ref.maha;
  const cols = m.data[0].map((_, j) => m.data.map(r => r[j]));
  const r = D.mahalanobisD2(cols);
  assert.equal(r.k, m.k);
  r.d2.forEach((v, i) => assert.ok(near(v, m.d2[i], 1e-9), `D² ${i}: ${v} vs ${m.d2[i]}`));
  const crit = dist.qchisq(0.001, m.k, false);
  assert.ok(near(crit, m.critical, 1e-8));
  r.d2.forEach((v, i) => assert.ok(near(dist.pchisq(v, m.k, false), m.p[i], 1e-7)));
  assert.deepEqual(r.d2.map((v, i) => [v, i]).filter(x => x[0] > crit).map(x => x[1]), m.flagged);
});

/* ---------- 결과 객체 ---------- */

test('runDistance: 표 구성, 행렬은 50×50까지만 보이고 요약은 전체로 계산, 단위 경고', () => {
  const m = ref.maha, X = m.data;                                                   // 40사례 × 3변수
  const names = ['a', 'b', 'c'];
  const res = D.runDistance(dataset(X, names), { vars: names, options: { measure: 'euclid' } });
  assert.equal(res.title, '거리·유사도');
  assert.deepEqual(table(res, '처리 요약').rows.map(r => r.n), [40, 0, 3, 780]);
  const mat = table(res, '유클리드 거리 행렬 (사례 간)');
  assert.equal(mat.rows.length, 40);
  assert.equal(mat.columns.length, 41);
  assert.equal(mat.rows[0].label, '행 2');
  assert.equal(mat.rows[0].c0, 0);
  const sm = table(res, '유클리드 거리 요약');
  assert.ok(sm.rows[0].value <= sm.rows[1].value && sm.rows[0].value <= sm.rows[2].value);
  assert.ok(R.renderResultHtml(res).includes('거리 행렬'));
  assert.ok(Ex.resultToGrid(res).rows.length > 40);
  // 단위 경고: 한 변수만 단위가 크다
  const big = X.map(r => [r[0], r[1] * 100, r[2]]);
  const w = D.runDistance(dataset(big, names), { vars: names, options: { measure: 'euclid' } });
  assert.ok(w.warnings.some(x => x.includes('표준편차가 크게 다릅니다')));
  const z = D.runDistance(dataset(big, names), { vars: names, options: { measure: 'euclid', standardize: 'z' } });
  assert.ok(!z.warnings.some(x => x.includes('표준편차가 크게 다릅니다')));
  // 마할라노비스·상관 거리는 단위에 불변이라 경고하지 않는다
  assert.ok(!D.runDistance(dataset(big, names), { vars: names, options: { measure: 'mahalanobis' } }).warnings.length);
});

test('runDistance: 화면 행렬 한도(50)와 새 탭 저장 한도(500), 한도를 넘으면 안내', () => {
  const rows = [];
  for (let i = 0; i < 60; i++) rows.push([i % 7, (i * 3) % 11, (i * 5) % 13]);
  const names = ['a', 'b', 'c'];
  const res = D.runDistance(dataset(rows, names), { vars: names, options: { measure: 'cityblock', saveMatrix: true } });
  const mat = table(res, '맨해튼(블록) 거리 행렬');
  assert.equal(mat.rows.length, 50);
  assert.ok(mat.notes.some(n => n.includes('앞 50×50')));
  assert.equal(res.saveSheet.rows.length, 61);                                       // 머리글 + 60행
  assert.equal(res.saveSheet.rows[0].length, 61);
  assert.equal(res.saveSheet.name, '거리행렬');
  // 요약은 전체 쌍 1770개
  assert.ok(table(res, '맨해튼(블록) 거리 요약 (쌍 1770개)'));
  const many = [];
  for (let i = 0; i < 600; i++) many.push([i % 7, (i * 3) % 11, (i * 5) % 13 + (i % 2)]);
  assert.throws(() => D.runDistance(dataset(many, names), { vars: names, options: { measure: 'euclid', saveMatrix: true } }), /최대 500개/);
  const huge = [];
  for (let i = 0; i < 1200; i++) huge.push([i % 7, (i * 3) % 11, (i * 5) % 13 + (i % 2)]);
  assert.throws(() => D.runDistance(dataset(huge, names), { vars: names, options: { measure: 'euclid' } }), /1000건 이하/);
  const o = D.runDistance(dataset(huge, names), { vars: names, options: { measure: 'euclid', outlier: true } });      // 이상치 점검만
  assert.ok(o.warnings.some(w => w.includes('이상치 점검만')));
  assert.ok(table(o, '마할라노비스 이상치 점검'));
  assert.ok(!table(o, '유클리드 거리 행렬'));
});

test('runDistance: 마할라노비스 이상치 표는 기준값의 이상치 후보를 가리키고 D² 저장 열은 사례 순서에 맞는다', () => {
  const m = ref.maha, names = ['a', 'b', 'c'];
  const rows = [names].concat(m.data.map(r => r.slice()));
  rows[4][1] = '';                                                                     // 결측 한 건 → 그 사례 제외 (기준값은 40건 기준이라 표만 확인)
  const ds = Data.buildDataset(rows, []);
  const res = D.runDistance(ds, { vars: names, options: { measure: 'euclid', standardize: 'z', outlier: true, saveD2: true } });
  assert.ok(res.warnings.some(w => w.includes('1건을 제외')));
  const col = res.saveColumns[0];
  assert.equal(col.name, '마할라노비스_D2');
  assert.equal(col.values.length, 40);
  assert.equal(col.values[3], null);                                                   // 결측 사례(사례 번호 3)는 빈 값
  assert.ok(col.values.filter(v => v !== null).every(v => v >= 0));
  // 결측 없는 자료로 기준값과 같은 후보
  const ok = D.runDistance(dataset(m.data, names), { vars: names, options: { outlier: true, measure: 'euclid' } });
  const summary = table(ok, '마할라노비스 이상치 점검');
  assert.ok(near(summary.rows[2].value, m.critical, 1e-6));
  assert.equal(summary.rows[4].value, m.flagged.length);
  const list = table(ok, '이상치 후보');
  assert.equal(list.rows[0].label, '행 ' + (m.top + 2));
  assert.ok(near(list.rows[0].d2, m.d2[m.top], 1e-9) && near(list.rows[0].p, m.p[m.top], 1e-7));
});

test('runDistance: 입력 검증 (변수 수, 이분 자료 전용, 변수 간 마할라노비스, 민코프스키 차수, 표준화 불가, 특이 공분산)', () => {
  const X = ref.maha.data, names = ['a', 'b', 'c'];
  const run = (opts, ds, vars) => D.runDistance(ds || dataset(X, names), { vars: vars || names, options: opts });
  assert.throws(() => run({}, null, ['a']), /2개 이상/);
  assert.throws(() => run({ measure: 'jaccard' }), /0 또는 1/);
  assert.throws(() => run({ measure: 'mahalanobis', target: 'var' }), /사례 간 거리에만/);
  assert.throws(() => run({ measure: 'minkowski', p: '0.5' }), /1 이상/);
  assert.throws(() => run({ measure: 'nope' }), /알 수 없는 측도/);
  assert.throws(() => run({ standardize: 'z' }, dataset([[1, 5, 2], [2, 5, 4], [3, 5, 1]], names)), /"b"의 값이 모두 같아/);
  const dup = [[1, 2, 2], [2, 4, 4], [3, 6, 6], [4, 9, 9], [5, 10, 10]];                // b = c: 공분산 특이
  assert.throws(() => run({ measure: 'mahalanobis' }, dataset(dup, names)), /공분산 행렬/);
  assert.throws(() => run({ outlier: true }, dataset(dup, names)), /공분산 행렬/);
  assert.throws(() => run({ outlier: true }, dataset([[1, 2, 3], [2, 4, 1], [5, 1, 4]], names)), /사례 수가 변수 수보다/);
  const txt = Data.buildDataset([['a', 'b'], [1, 'x'], [2, 'y']], []);
  assert.throws(() => D.runDistance(txt, { vars: ['a', 'b'], options: {} }), /숫자가 아닌 값/);
});

test('runDistance: 이분 자료 유사도와 변수 간 거리, 계산할 수 없는 쌍은 제외하고 알린다', () => {
  const c = ref.pairs.find(x => x.binary && x.measure === 'jaccard');
  const names = c.data[0].map((_, j) => 'v' + j);
  const res = D.runDistance(dataset(c.data, names), { vars: names, options: { measure: 'jaccard' } });
  const mat = table(res, 'Jaccard 계수 (유사도) 행렬');
  c.expected.forEach((row, i) => row.forEach((e, j) => assert.ok(near(mat.rows[i]['c' + j], e, 1e-12))));
  assert.ok(mat.notes[0].includes('비슷하다'));
  assert.equal(table(res, 'Jaccard 계수 (유사도) 요약 (쌍 45개)').rows[0].label, '가장 비슷한 쌍');
  // 모두 0인 두 사례는 Jaccard 를 계산할 수 없다
  const z = D.runDistance(dataset([[0, 0, 0], [0, 0, 0], [1, 0, 1], [1, 1, 0]], ['a', 'b', 'c']), { vars: ['a', 'b', 'c'], options: { measure: 'jaccard' } });
  assert.ok(table(z, 'Jaccard 계수 (유사도) 요약').notes[0].includes('계산할 수 없는 쌍'));
  // 변수 간 거리: 변수 이름이 행렬의 이름이 된다
  const v = D.runDistance(dataset(ref.pairs[0].data, ['연령', '소득', '점수', '만족']), { vars: ['연령', '소득', '점수', '만족'], options: { target: 'var', measure: 'corrDist', standardize: 'z' } });
  const vm = table(v, '상관 거리 (1 − r) 행렬 (변수 간)');
  assert.deepEqual(vm.rows.map(r => r.label), ['연령', '소득', '점수', '만족']);
  assert.equal(vm.rows[0].c0, 0);
});

/* ---------- 시트 시뮬레이션: 새 탭 저장과 D² 열 저장 ---------- */

const vm = require('node:vm');
const SRC = path.join(__dirname, '..', 'src');
function mkSheet(name, rows) {
  const sh = {
    name, rows: rows.map(r => r.slice()), maxRows: 1000, maxCols: 26, frozen: [0, 0], bold: 0,
    getName: () => sh.name,
    getLastRow: () => sh.rows.length, getMaxRows: () => sh.maxRows, getMaxColumns: () => sh.maxCols,
    insertRowsAfter(pos, n) { sh.maxRows += n; }, insertColumnsAfter(pos, n) { sh.maxCols += n; },
    setFrozenRows(n) { sh.frozen[0] = n; }, setFrozenColumns(n) { sh.frozen[1] = n; },
    getDataRange: () => ({ getValues: () => { const w = Math.max(1, ...sh.rows.map(r => r.length)); return sh.rows.map(r => Array.from({ length: w }, (_, j) => (j < r.length ? r[j] : ''))); } }),
    getRange(r, c, nr, nc) {
      const rng = {
        setValues(arr) {
          assert.ok(r - 1 + arr.length <= sh.maxRows && c - 1 + arr[0].length <= sh.maxCols, '그리드를 넘어 쓰려 함 (' + sh.name + ')');
          arr.forEach((a, i) => { while (sh.rows.length < r + i) sh.rows.push([]); a.forEach((v, j) => { const row = sh.rows[r - 1 + i]; while (row.length < c - 1 + j) row.push(''); row[c - 1 + j] = v; }); });
          return rng;
        },
        setFontWeight() { sh.bold++; return rng; }, setBackground: () => rng, setNumberFormat: () => rng, setDataValidation: () => rng
      };
      return rng;
    }
  };
  return sh;
}
function makeBook(dataRows, varRows) {
  const sheets = [mkSheet('데이터', dataRows), mkSheet('변수정의', varRows)];
  let active = sheets[0];
  const v = () => { const b = { requireValueInRange: () => b, requireValueInList: () => b, setAllowInvalid: () => b, setHelpText: () => b, build: () => ({}) }; return b; };
  const ss = {
    getSheetByName: n => sheets.find(s => s.name === n) || null,
    insertSheet: n => { const s = mkSheet(n, []); sheets.push(s); active = s; return s; },
    getActiveSheet: () => active, setActiveSheet: s => { active = s; }
  };
  return { sheets, ss, app: { newDataValidation: v, getActiveSpreadsheet: () => ss }, getActive: () => active };
}

test('시트 시뮬레이션: 거리 행렬을 새 탭에 저장하고 D² 열을 데이터 탭에 쓴다 (같은 이름이면 _2, 보던 탭 유지)', () => {
  const m = ref.maha;
  const dataRows = [['x1', 'x2', 'x3']].concat(m.data);
  const varRows = [['변수명', '라벨', '척도', '값 레이블', '결측값 코드'], ['x1', '', '연속', '', ''], ['x2', '', '연속', '', ''], ['x3', '', '연속', '', '']];
  const book = makeBook(dataRows, varRows);
  const ctx = vm.createContext({ console: { log() {}, error() {} }, SpreadsheetApp: book.app, Utilities: {}, Session: {}, PropertiesService: {} });
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  ctx.__req = { analysis: 'distance', vars: ['x1', 'x2', 'x3'], options: { measure: 'euclid', standardize: 'z', outlier: true, saveD2: true, saveMatrix: true } };
  const out = JSON.parse(JSON.stringify(vm.runInContext('runAnalysis(__req)', ctx)));
  assert.ok(out.html.includes('거리행렬') && out.html.includes('마할라노비스_D2'));
  const tab = book.ss.getSheetByName('거리행렬');
  assert.equal(tab.rows.length, 41);
  assert.equal(tab.rows[0].length, 41);
  assert.equal(tab.rows[1][1], 0);                                                     // 대각선
  assert.deepEqual(tab.frozen, [1, 1]);
  assert.equal(book.getActive().name, '데이터');                                         // 새 탭으로 옮겨 가지 않음
  assert.equal(out.result.saveSheet, undefined);
  assert.equal(book.ss.getSheetByName('데이터').rows[0][3], '마할라노비스_D2');
  assert.ok(book.ss.getSheetByName('데이터').rows.slice(1).every(r => typeof r[3] === 'number' && r[3] >= 0));
  assert.deepEqual(book.ss.getSheetByName('변수정의').rows.slice(-1)[0].slice(0, 3), ['마할라노비스_D2', '마할라노비스 D² (x1, x2, x3)', '연속']);
  vm.runInContext('runAnalysis(__req)', ctx);
  assert.ok(book.ss.getSheetByName('거리행렬_2'));
  assert.equal(book.ss.getSheetByName('데이터').rows[0][4], '마할라노비스_D2_2');
});
