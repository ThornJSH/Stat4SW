const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const Data = require('../src/Data.js');
const K = require('../src/Compute.js');

const SRC = path.join(__dirname, '..', 'src');
const json = x => JSON.parse(JSON.stringify(x));
const close = (a, b, tol = 1e-12) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b));

/* ---------- 순수 함수 ---------- */

test('reverseCode: 새 값 = 최솟값 + 최댓값 − 값, 결측 유지, 범위 밖은 결측으로 세기', () => {
  const r = K.reverseCode([1, 2, 3, 4, 5, null, 6, 0, '2'], 1, 5);
  assert.deepEqual(r.values, [5, 4, 3, 2, 1, null, null, null, 4]);
  assert.equal(r.outOfRange, 2);
  assert.deepEqual(K.reverseCode([0, 10, 4], 0, 10).values, [10, 0, 6]);
  assert.deepEqual(K.reverseCode([1.5, 2.5], 1, 4).values, [3.5, 2.5]);
  assert.throws(() => K.reverseCode([1], 5, 1), /최솟값은 최댓값보다/);
  assert.throws(() => K.reverseCode([1], 3, 3), /최솟값은 최댓값보다/);
  assert.throws(() => K.reverseCode([1], NaN, 5), /숫자로 입력/);
  // 두 번 뒤집으면 원래 값
  const back = K.reverseCode(K.reverseCode([1, 2, 3, 4, 5], 1, 5).values, 1, 5).values;
  assert.deepEqual(back, [1, 2, 3, 4, 5]);
});

test('rowAggregate: 합계·평균, 최소 유효 문항 수, 결측 처리', () => {
  const cols = [[1, 2, null, 4], [3, null, null, 5], [5, 6, null, null]];
  const mean = K.rowAggregate(cols, 'mean', 2);
  assert.deepEqual(mean.values, [3, 4, null, 4.5]);
  assert.equal(mean.computed, 3);
  assert.equal(mean.missing, 1);
  const sum = K.rowAggregate(cols, 'sum', 3);
  assert.deepEqual(sum.values, [9, null, null, null]);
  assert.deepEqual(K.rowAggregate(cols, 'sum', 1).values, [9, 8, null, 9]);
  assert.deepEqual(K.rowAggregate([['1', 2]].map(r => [r[0]]).concat([[2]]), 'sum', 2).values, [3]);   // 숫자 문자열도 숫자로
  assert.throws(() => K.rowAggregate(cols, 'median', 1), /알 수 없는 계산/);
});

test('zScores: 평균 0, 표본 표준편차 1, 결측 유지, 계산 불가는 모두 결측', () => {
  const z = K.zScores([2, 4, 4, 4, 5, 5, 7, 9, null]);
  assert.equal(z.n, 8);
  assert.ok(close(z.mean, 5));
  assert.ok(close(z.sd, Math.sqrt(32 / 7)));
  assert.equal(z.values[8], null);
  const valid = z.values.filter(v => v !== null);
  assert.ok(Math.abs(valid.reduce((s, v) => s + v, 0)) < 1e-12);
  assert.ok(close(valid.reduce((s, v) => s + v * v, 0), 7));          // Σz² = n − 1
  assert.ok(close(z.values[0], (2 - 5) / Math.sqrt(32 / 7)));
  assert.deepEqual(K.zScores([5, 5, 5]).values, [null, null, null]);   // 표준편차 0
  assert.deepEqual(K.zScores([3, null]).values, [null, null]);         // 유효한 값 1개
  // 큰 평균에서도 안정적
  const big = K.zScores([1e9 + 1, 1e9 + 2, 1e9 + 3]);
  assert.ok(close(big.values[0], -1) && close(big.values[2], 1));
});

test('calcColLetter_, calcCheckName_', () => {
  assert.deepEqual([1, 5, 26, 27, 52, 53, 702, 703].map(K.calcColLetter_), ['A', 'E', 'Z', 'AA', 'AZ', 'BA', 'ZZ', 'AAA']);
  const taken = { a: true };
  assert.equal(K.calcCheckName_('  b  ', taken), 'b');
  assert.throws(() => K.calcCheckName_('b', taken), /이미 있는 변수 이름/);
  assert.throws(() => K.calcCheckName_('a', taken), /이미 있는 변수 이름/);
  assert.throws(() => K.calcCheckName_('  ', taken), /이름을 입력/);
  assert.throws(() => K.calcCheckName_('x'.repeat(61), {}), /60자 이하/);
});

/* ---------- planCompute_ ---------- */

const defs = [
  { name: 'q1', label: '문항1', scale: 'ordinal', valueLabels: { '1': '전혀', '5': '매우' }, missing: [] },
  { name: 'q2', label: '', scale: 'ordinal', valueLabels: {}, missing: [] },
  { name: 'q3', label: '', scale: 'ordinal', valueLabels: {}, missing: [99] },
  { name: 'g', label: '', scale: 'nominal', valueLabels: {}, missing: [] }
];
const rows = [['q1', 'q2', 'q3', 'g'], [1, 2, 3, 'a'], [5, 4, 99, 'b'], [3, '', 2, 'a'], ['', '', '', ''], [2, 2, 2, 'b'], [4, 6, 1, 'a']];
const ds = () => Data.buildDataset(rows, defs);

test('buildDataset: rowNumbers 는 빈 행을 뺀 사례의 시트 행 번호', () => {
  assert.deepEqual(ds().rowNumbers, [2, 3, 4, 6, 7]);
});

test('planCompute_: 역코딩 — 이름, 라벨, 척도, 범위 밖 경고', () => {
  const p = K.planCompute_(ds(), { kind: 'reverse', vars: ['q1', 'q2'], min: '1', max: '5', suffix: '_r' });
  assert.deepEqual(p.columns.map(c => c.name), ['q1_r', 'q2_r']);
  assert.deepEqual(p.columns[0].values, [5, 1, 3, 4, 2]);
  assert.equal(p.columns[0].label, '문항1 (역)');
  assert.equal(p.columns[1].label, 'q2 (역)');
  assert.equal(p.columns[0].scale, '서열');
  assert.deepEqual(p.columns[1].values, [4, 2, null, 4, null]);       // 6은 범위 밖 → 결측
  assert.ok(p.warnings.some(w => w.includes('범위(1~5) 밖의 값 1건')));
  assert.deepEqual([p.columns[0].valid, p.columns[0].missing], [5, 0]);
  assert.deepEqual([p.columns[1].valid, p.columns[1].missing], [3, 2]);
});

test('planCompute_: 결측값 코드(99)는 합계·평균에서 빠지고, 최소 유효 문항 수 적용', () => {
  const p = K.planCompute_(ds(), { kind: 'mean', vars: ['q1', 'q2', 'q3'], name: '평균점수', minValid: '2' });
  assert.equal(p.columns.length, 1);
  assert.equal(p.columns[0].name, '평균점수');
  assert.equal(p.columns[0].label, '3문항 평균');
  assert.deepEqual(p.columns[0].values, [2, 4.5, 2.5, 2, 11 / 3]);
  const all = K.planCompute_(ds(), { kind: 'sum', vars: ['q1', 'q2', 'q3'], name: '합계' });
  assert.deepEqual(all.columns[0].values, [6, null, null, 6, 11]);     // 기본: 3문항 모두 유효한 사례만
  assert.ok(all.warnings.some(w => w.includes('결측 2건')));
  const loose = K.planCompute_(ds(), { kind: 'sum', vars: ['q1', 'q2', 'q3'], name: '합계', minValid: 1 });
  assert.ok(loose.warnings.some(w => w.includes('평균을 쓰는 것을 고려')));
});

test('planCompute_: z점수', () => {
  const p = K.planCompute_(ds(), { kind: 'z', vars: ['q1'], suffix: '_z' });
  assert.equal(p.columns[0].name, 'q1_z');
  assert.equal(p.columns[0].label, '문항1 (z점수)');
  assert.equal(p.columns[0].scale, '연속');
  assert.ok(Math.abs(p.columns[0].mean) < 1e-12);
});

test('planCompute_: 입력 검증', () => {
  const run = req => K.planCompute_(ds(), req);
  assert.throws(() => run({ kind: 'reverse', vars: [], min: 1, max: 5, suffix: '_r' }), /변수를 선택/);
  assert.throws(() => run({ kind: 'reverse', vars: ['g'], min: 1, max: 5, suffix: '_r' }), /숫자가 아닌 값/);
  assert.throws(() => run({ kind: 'reverse', vars: ['q1'], min: '', max: '', suffix: '_r' }), /최솟값과 최댓값을 입력/);
  assert.throws(() => run({ kind: 'reverse', vars: ['q1'], min: 1, max: 5, suffix: '' }), /접미사/);
  assert.throws(() => run({ kind: 'reverse', vars: ['q1'], min: 5, max: 1, suffix: '_r' }), /최솟값은 최댓값보다/);
  assert.throws(() => run({ kind: 'sum', vars: ['q1'], name: 's' }), /2개 이상/);
  assert.throws(() => run({ kind: 'sum', vars: ['q1', 'q2'], name: 'q1' }), /이미 있는 변수 이름/);
  assert.throws(() => run({ kind: 'sum', vars: ['q1', 'q2'], name: 's', minValid: 3 }), /1 이상 2 이하/);
  assert.throws(() => run({ kind: 'sum', vars: ['q1', 'q2'], name: 's', minValid: 1.5 }), /정수/);
  assert.throws(() => run({ kind: 'reverse', vars: ['nope'], min: 1, max: 5, suffix: '_r' }), /찾을 수 없습니다/);
  assert.throws(() => run({ kind: 'x', vars: ['q1'] }), /알 수 없는 계산/);
  assert.throws(() => K.planCompute_(Data.buildDataset([['a']], []), { kind: 'z', vars: ['a'], suffix: '_z' }), /데이터가 없습니다/);
  // 명목 척도로 정의된 숫자 변수는 역코딩할 수 없다
  const nominalDs = Data.buildDataset([['n'], [1], [2], [3]], [{ name: 'n', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]);
  assert.throws(() => K.planCompute_(nominalDs, { kind: 'reverse', vars: ['n'], min: 1, max: 3, suffix: '_r' }), /명목 척도라 역코딩할 수 없습니다/);
  // 새 이름이 서로 겹침 (같은 접미사로 두 변수가 같은 이름이 될 수는 없지만 기존 열과 겹칠 수는 있다)
  const dup = Data.buildDataset([['a', 'a_r'], [1, 2], [2, 3]], []);
  assert.throws(() => K.planCompute_(dup, { kind: 'reverse', vars: ['a'], min: 1, max: 2, suffix: '_r' }), /이미 있는 변수 이름/);
});

/* ========== Apps Script 전역 범위 시뮬레이션: 시트 읽기·쓰기 ========== */

function loadGasWith(overrides) {
  const env = { console: { log() {}, error() {} }, Utilities: {}, Session: {}, PropertiesService: {} };
  Object.assign(env, overrides);
  const ctx = vm.createContext(env);
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  return ctx;
}

/** 값을 기억하는 시트 mock (getDataRange 는 실제처럼 직사각형) */
function makeSheet(rows, maxCols) {
  const log = { formats: [], validations: [], bold: [], inserted: 0 };
  let cols = maxCols || 26;
  const width = () => Math.max(1, ...rows.map(r => r.length));
  const sheet = {
    log, rows,
    getLastRow: () => rows.length,
    getMaxColumns: () => cols,
    insertColumnsAfter: (pos, n) => { cols += n; log.inserted += n; },
    getDataRange: () => ({ getValues: () => rows.map(r => Array.from({ length: width() }, (_, j) => (j < r.length ? r[j] : ''))) }),
    getRange: (r, c, nr, nc) => {
      nr = nr === undefined ? 1 : nr; nc = nc === undefined ? 1 : nc;
      const rng = {
        tag: 'range', r, c, nr, nc,
        setValues: (arr) => {
          if (arr.length !== nr || arr.some(a => a.length !== nc)) throw new Error('setValues 크기 불일치');
          arr.forEach((a, i) => {
            while (rows.length < r + i) rows.push([]);
            a.forEach((v, j) => { const row = rows[r - 1 + i]; while (row.length < c - 1 + j) row.push(''); row[c - 1 + j] = v; });
          });
          return rng;
        },
        setNumberFormat: (f) => { log.formats.push({ r, c, nr, nc, f }); return rng; },
        setFontWeight: (w) => { log.bold.push({ r, c, nc, w }); return rng; },
        setBackground: () => rng,
        setDataValidation: (v) => { log.validations.push({ r, c, nr, nc, v }); return rng; }
      };
      return rng;
    }
  };
  return sheet;
}

function validationBuilder() {
  const b = { _: {}, requireValueInRange() { return b; }, requireValueInList() { return b; }, setAllowInvalid() { return b; }, setHelpText() { return b; }, build() { return b._; } };
  return b;
}

const VAR_HEADER = ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'];

function calcApp(dataRows, varRows) {
  const data = makeSheet(dataRows), vars = makeSheet(varRows);
  const app = {
    newDataValidation: validationBuilder,
    getActiveSpreadsheet: () => ({ getSheetByName: (n) => (n === '데이터' ? data : n === '변수정의' ? vars : null) })
  };
  return { app, data, vars };
}

const SIM_DATA = () => [['q1', 'q2', 'q3', 'g'], [1, 2, 3, 'a'], [5, 4, 99, 'b'], [3, '', 2, 'a'], ['', '', '', ''], [2, 2, 2, 'b'], [4, 6, 1, 'a']];
const SIM_VARS = () => [VAR_HEADER, ['q1', '문항1', '서열', '1=전혀, 5=매우', ''], ['q2', '', '서열', '', ''], ['q3', '', '서열', '', '99'], ['g', '', '명목', '', '']];

test('시트 시뮬레이션: calcState 는 숫자 여부와 역코딩 범위 힌트를 준다', () => {
  const { app } = calcApp(SIM_DATA(), SIM_VARS());
  const ctx = loadGasWith({ SpreadsheetApp: app });
  const st = json(vm.runInContext('calcState()', ctx));
  assert.equal(st.n, 5);
  const by = Object.fromEntries(st.vars.map(v => [v.name, v]));
  assert.equal(by.g.numeric, false);
  assert.deepEqual([by.q1.numeric, by.q1.dataMin, by.q1.dataMax, by.q1.labelMin, by.q1.labelMax], [true, 1, 5, 1, 5]);
  assert.deepEqual([by.q3.dataMin, by.q3.dataMax, by.q3.labelMin], [1, 3, null]);        // 99는 결측 코드라 제외
});

test('시트 시뮬레이션: calcPlan 은 시트를 바꾸지 않고 열 위치와 요약을 돌려준다', () => {
  const { app, data, vars } = calcApp(SIM_DATA(), SIM_VARS());
  const ctx = loadGasWith({ SpreadsheetApp: app });
  const before = JSON.stringify([data.rows, vars.rows]);
  ctx.__req = { kind: 'reverse', vars: ['q1', 'q2'], min: '1', max: '5', suffix: '_r' };
  const p = json(vm.runInContext('calcPlan(__req)', ctx));
  assert.deepEqual(p.columns.map(c => [c.name, c.column]), [['q1_r', 'E'], ['q2_r', 'F']]);
  assert.deepEqual(p.columns[0].sample, [5, 1, 3, 4, 2]);
  assert.equal(JSON.stringify([data.rows, vars.rows]), before);
});

test('시트 시뮬레이션: calcApply 는 맨 오른쪽에 새 열을 쓰고(빈 행 정렬 유지) 변수정의에 등록한다', () => {
  const { app, data, vars } = calcApp(SIM_DATA(), SIM_VARS());
  const ctx = loadGasWith({ SpreadsheetApp: app });
  ctx.__req = { kind: 'reverse', vars: ['q1', 'q2'], min: '1', max: '5', suffix: '_r' };
  const r = json(vm.runInContext('calcApply(__req)', ctx));
  assert.deepEqual([r.firstColumn, r.lastColumn, r.rows, r.added], ['E', 'F', 6, ['q1_r', 'q2_r']]);
  // 머리글과 값: 원래 데이터 4열 오른쪽에, 시트의 행 위치 그대로 (5행은 빈 행)
  assert.deepEqual(data.rows[0].slice(4, 6), ['q1_r', 'q2_r']);
  assert.deepEqual(data.rows.slice(1).map(row => row.slice(4, 6)), [[5, 4], [1, 2], [3, ''], ['', ''], [4, 4], [2, '']]);
  assert.deepEqual(data.rows[1].slice(0, 4), [1, 2, 3, 'a']);                 // 원본 열은 그대로
  assert.ok(data.log.bold.some(b => b.r === 1 && b.c === 5 && b.w === 'bold'));
  // 변수정의에 마지막 정의 다음 행부터 등록
  assert.deepEqual(vars.rows.slice(5), [['q1_r', '문항1 (역)', '서열', '', ''], ['q2_r', 'q2 (역)', '서열', '', '']]);
  assert.deepEqual(vars.rows[1], ['q1', '문항1', '서열', '1=전혀, 5=매우', '']);   // 기존 정의 그대로
  // 새 열은 분석에도 바로 쓰인다 (끝에서 끝으로)
  const info = json(vm.runInContext('getDatasetInfo()', ctx));
  assert.deepEqual(info.vars.map(v => v.name), ['q1', 'q2', 'q3', 'g', 'q1_r', 'q2_r']);
  assert.equal(info.vars[4].valid, 5);
});

test('시트 시뮬레이션: 합계·평균 열, 열이 모자라면 늘리고, z점수는 숫자 서식', () => {
  const { app, data, vars } = calcApp(SIM_DATA(), SIM_VARS());
  const ctx = loadGasWith({ SpreadsheetApp: app });
  ctx.__req = { kind: 'mean', vars: ['q1', 'q2', 'q3'], name: '평균점수', minValid: 2 };
  const r = json(vm.runInContext('calcApply(__req)', ctx));
  assert.equal(r.firstColumn, 'E');
  assert.deepEqual(data.rows.slice(1).map(row => row[4]), [2, 4.5, 2.5, '', 2, 11 / 3]);
  assert.ok(data.log.formats.some(f => f.c === 5 && f.f === '0.000'));
  assert.deepEqual(vars.rows[5], ['평균점수', '3문항 평균', '연속', '', '']);
  ctx.__req2 = { kind: 'z', vars: ['q1'], suffix: '_z' };
  const r2 = json(vm.runInContext('calcApply(__req2)', ctx));
  assert.equal(r2.firstColumn, 'F');                                           // 방금 쓴 열 다음
  assert.deepEqual(vars.rows[6].slice(0, 3), ['q1_z', '문항1 (z점수)', '연속']);
  // 열 수가 모자라면 늘린다 (최대 열 수 4+1+1 → 처음 26이라 늘지 않음)
  assert.equal(data.log.inserted, 0);
  const small = calcApp(SIM_DATA(), SIM_VARS());
  small.data.getMaxColumns = () => 4;
  const ctx2 = loadGasWith({ SpreadsheetApp: small.app });
  ctx2.__r = { kind: 'reverse', vars: ['q1'], min: 1, max: 5, suffix: '_r' };
  vm.runInContext('calcApply(__r)', ctx2);
  assert.equal(small.data.log.inserted, 1);
});

test('시트 시뮬레이션: 이미 있는 이름과 숫자가 아닌 변수는 거부하고 시트를 바꾸지 않는다', () => {
  const { app, data, vars } = calcApp(SIM_DATA(), SIM_VARS());
  const ctx = loadGasWith({ SpreadsheetApp: app });
  const before = JSON.stringify([data.rows, vars.rows]);
  ctx.__a = { kind: 'sum', vars: ['q1', 'q2'], name: 'q3' };
  assert.throws(() => vm.runInContext('calcApply(__a)', ctx), /이미 있는 변수 이름/);
  ctx.__b = { kind: 'z', vars: ['g'], suffix: '_z' };
  assert.throws(() => vm.runInContext('calcApply(__b)', ctx), /숫자가 아닌 값/);
  assert.equal(JSON.stringify([data.rows, vars.rows]), before);
  // 두 번 같은 계산을 하면 두 번째는 이름이 겹쳐 거부된다
  ctx.__c = { kind: 'reverse', vars: ['q1'], min: 1, max: 5, suffix: '_r' };
  vm.runInContext('calcApply(__c)', ctx);
  assert.throws(() => vm.runInContext('calcApply(__c)', ctx), /이미 있는 변수 이름/);
});
