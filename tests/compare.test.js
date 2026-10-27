const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const C = require('../src/Stats_Compare.js');
const R = require('../src/Render.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'compare.json'), 'utf8'));

function relErr(a, e) {
  if (a === e) return 0;
  return Math.abs(a - e) / Math.max(Math.abs(e), 1e-300);
}
/** 값은 상대오차, 0에 가까운 값(평균차이의 신뢰구간 하한 등)은 절대오차 1e-10을 함께 허용 */
function close(a, e, tol, msg) {
  assert.ok(relErr(a, e) < tol || Math.abs(a - e) < 1e-10, `${msg}: ${a} vs ${e}`);
}

const ds = (values, defs) => Data.buildDataset(values, defs || []);
const tbl = (res, caption) => res.blocks.find(b => b.caption && b.caption.startsWith(caption));

/* ---------- 계산: scipy 대조 ---------- */

test('oneSampleT: scipy 대조 (통계량, 자유도, p, 평균차이, 신뢰구간)', () => {
  for (const c of ref.one) {
    const r = C.oneSampleT(c.values, c.mu, c.conf);
    const m = `${c.name} conf=${c.conf}`;
    close(r.t, c.t, 1e-9, m + ' t');
    assert.equal(r.df, c.df, m);
    close(r.p2, c.p, 1e-9, m + ' p');
    close(r.diff, c.diff, 1e-9, m + ' diff');
    close(r.lo, c.lo, 1e-8, m + ' lo');
    close(r.hi, c.hi, 1e-8, m + ' hi');
    close(r.p1, c.p / 2, 1e-9, m + ' p1');
  }
});

test('independentT: scipy 대조 (등분산 가정, Welch, Levene)', () => {
  for (const c of ref.ind) {
    const r = C.independentT(c.a, c.b, c.conf);
    const m = `${c.name} conf=${c.conf}`;
    for (const k of ['equal', 'welch']) {
      close(r[k].diff, c.diff, 1e-9, `${m} ${k} diff`);
      close(r[k].t, c[k].t, 1e-9, `${m} ${k} t`);
      close(r[k].df, c[k].df, 1e-9, `${m} ${k} df`);
      close(r[k].p2, c[k].p, 1e-9, `${m} ${k} p`);
      close(r[k].lo, c[k].lo, 1e-8, `${m} ${k} lo`);
      close(r[k].hi, c[k].hi, 1e-8, `${m} ${k} hi`);
    }
    close(r.levene.W, c.levene.W, 1e-9, `${m} levene W`);
    close(r.levene.p, c.levene.p, 1e-8, `${m} levene p`);
  }
});

test('pairedT: scipy 대조 (t검정과 상관계수)', () => {
  for (const c of ref.paired) {
    const r = C.pairedT(c.x, c.y, c.conf);
    const m = c.name;
    close(r.test.diff, c.diff, 1e-9, m + ' diff');
    close(r.test.t, c.t, 1e-9, m + ' t');
    assert.equal(r.test.df, c.df, m);
    close(r.test.p2, c.p, 1e-9, m + ' p');
    close(r.test.lo, c.lo, 1e-8, m + ' lo');
    close(r.test.hi, c.hi, 1e-8, m + ' hi');
    close(r.corr.r, c.r, 1e-10, m + ' r');
    close(r.corr.p, c.rp, 1e-8, m + ' r p');
  }
});

test('계산할 수 없는 경우는 null (표본 1개, 표준편차 0)', () => {
  assert.equal(C.oneSampleT([5], 0, 95).t, null);
  assert.equal(C.oneSampleT([5, 5, 5], 0, 95).t, null);
  assert.equal(C.oneSampleT([5, 5, 5], 0, 95).diff, 5);
  const r = C.independentT([1], [2, 3, 4], 95);
  assert.equal(r.equal.t, null);
  assert.equal(r.welch.t, null);
  assert.equal(C.independentT([2, 2], [3, 3], 95).equal.t, null); // 합동분산 0
});

/* ---------- 실행 함수 ---------- */

const defs = Data.parseVarDefs([['h'],
  ['g', '집단', '명목', '1=대조, 2=처치', ''],
  ['h3', '3집단', '명목', '', '']
]).defs;

function sample() {
  return ds([
    ['y', 'z', 'g', 'h3'],
    [10, 12, 1, 1], [12, 15, 1, 1], [11, 11, 1, 2], [13, 14, 1, 2], [9, 10, 1, 3],
    [15, 14, 2, 1], [18, 20, 2, 2], [14, 17, 2, 3], [20, 19, 2, 3], [17, 18, 2, 1],
    ['', 30, 2, 1], [16, '', 1, 2], [5, 5, '', 1]
  ], defs);
}

test('runTTestOne: 표 구조, 검정값·신뢰수준 옵션, 결측 제외', () => {
  const res = C.runTTestOne(sample(), { vars: ['y'], options: { testValue: 12, conf: 99 } });
  const s = tbl(res, '일표본 통계량');
  assert.equal(s.rows[0].n, 12);
  const t = tbl(res, '일표본 검정');
  assert.ok(t.caption.includes('검정값 = 12'));
  assert.equal(t.rows[0].df, 11);
  assert.ok(t.columns.some(c => c.label === '차이의 99% 신뢰구간 하한'));
  assert.ok(Math.abs(t.rows[0].diff - (s.rows[0].mean - 12)) < 1e-12);
  assert.equal(t.rows[0].p1, t.rows[0].p2 / 2);
  assert.ok(R.renderResultHtml(res).includes('일표본 검정'));
});

test('runTTestOne: 옵션 오류와 검정 불가 경고', () => {
  assert.throws(() => C.runTTestOne(sample(), { vars: ['y'], options: { testValue: 'abc' } }), /검정값/);
  assert.throws(() => C.runTTestOne(sample(), { vars: ['y'], options: { conf: 100 } }), /신뢰수준/);
  assert.throws(() => C.runTTestOne(sample(), { vars: ['y'], options: { conf: 0 } }), /신뢰수준/);
  const res = C.runTTestOne(ds([['c'], [4], [4], [4]]), { vars: ['c'] });
  assert.equal(res.warnings.length, 1);
});

test('runTTestIndependent: 집단 통계, 두 가정의 행, 값 레이블, 차이 방향', () => {
  const res = C.runTTestIndependent(sample(), { dependent: ['y'], factor: ['g'] });
  const g = tbl(res, '집단통계량');
  assert.deepEqual(g.rows.map(r => r.group), ['대조', '처치']);
  assert.deepEqual(g.rows.map(r => r.n), [6, 5]);   // 대조: y가 유효한 6건(행 12 포함), 처치: 5건(행 11은 y 결측)
  const t = tbl(res, '독립표본 검정');
  assert.deepEqual(t.rows.map(r => r.assume), ['등분산을 가정함', '등분산을 가정하지 않음']);
  assert.equal(t.rows[1].name, '');
  assert.equal(t.rows[1].levF, null);
  assert.ok(t.rows[0].diff < 0);                        // 대조 평균 < 처치 평균
  assert.equal(t.rows[0].df, 9);
  assert.ok(t.rows[1].df > 0 && t.rows[1].df <= 9);
  assert.ok(t.notes.some(n => n.includes('집단1: 대조') && n.includes('집단2: 처치')));
  assert.ok(R.renderResultHtml(res).includes('등분산을 가정하지 않음'));
});

test('runTTestIndependent: 집단이 2개가 아니면 오류, 종속·집단 같은 변수는 오류', () => {
  assert.throws(() => C.runTTestIndependent(sample(), { dependent: ['y'], factor: ['h3'] }), /3개/);
  assert.throws(() => C.runTTestIndependent(ds([['y', 'g'], [1, 1], [2, 1]]), { dependent: ['y'], factor: ['g'] }), /1개/);
  assert.throws(() => C.runTTestIndependent(sample(), { dependent: ['g'], factor: ['g'] }), /같은 변수/);
  assert.throws(() => C.runTTestIndependent(sample(), { dependent: ['y'], factor: [] }), /변수를 선택/);
});

test('runTTestIndependent: 종속변수 여러 개', () => {
  const res = C.runTTestIndependent(sample(), { dependent: ['y', 'z'], factor: ['g'] });
  assert.equal(tbl(res, '집단통계량').rows.length, 4);
  assert.equal(tbl(res, '독립표본 검정').rows.length, 4);
});

test('runTTestPaired: 세 표, 평균은 변수1 − 변수2, 결측 쌍 제외', () => {
  const res = C.runTTestPaired(sample(), { first: ['y'], second: ['z'] });
  const stats = tbl(res, '대응표본 통계량');
  assert.deepEqual(stats.rows.map(r => r.n), [11, 11]);          // y, z 모두 유효한 11쌍 (행 11은 y, 행 12는 z 결측)
  const t = tbl(res, '대응표본 검정');
  assert.equal(t.rows[0].name, 'y − z');
  assert.equal(t.rows[0].df, 10);
  assert.ok(Math.abs(t.rows[0].diff - (stats.rows[0].mean - stats.rows[1].mean)) < 1e-12);
  const corr = tbl(res, '대응표본 상관계수');
  assert.equal(corr.rows[0].n, 11);
  assert.ok(corr.rows[0].r > 0.8);
  assert.ok(t.notes.some(n => n.includes('2건 제외')));
});

test('runTTestPaired: 입력 오류', () => {
  assert.throws(() => C.runTTestPaired(sample(), { first: ['y'], second: ['y'] }), /같은 변수/);
  assert.throws(() => C.runTTestPaired(sample(), { first: ['y', 'z'], second: ['z'] }), /하나씩/);
  assert.throws(() => C.runTTestPaired(sample(), { first: [], second: ['z'] }), /변수를 선택/);
});

/* ---------- 독립표본 t검정: Levene 결과로 읽을 행 강조 ---------- */

test('runTTestIndependent: Levene p ≥ .05 이면 "등분산을 가정함" 행, p < .05 이면 "가정하지 않음" 행을 굵게', () => {
  const Data = require('../src/Data.js');
  const R = require('../src/Render.js');
  const Ex = require('../src/Export.js');
  const make = (a, b) => Data.buildDataset([['g', 'y']].concat(a.map(v => [1, v]), b.map(v => [2, v])), []);
  const run = d => C.runTTestIndependent(d, { dependent: ['y'], factor: ['g'], options: {} });
  const testTable = res => res.blocks.find(b => b.caption === '독립표본 검정');

  // 분산이 비슷함 → Levene p 큼 → 첫 행(등분산을 가정함)
  const eq = run(make([1, 2, 3, 4, 5, 6], [2, 3, 4, 5, 6, 7]));
  let rows = testTable(eq).rows;
  assert.ok(rows[0].levP >= 0.05);
  assert.equal(rows[0].assume, '등분산을 가정함');
  assert.equal(rows[0]._bold, true);
  assert.equal(rows[1]._bold, false);

  // 분산이 매우 다름 → Levene p 작음 → 둘째 행(등분산을 가정하지 않음)
  const uneq = run(make([10, 10.1, 9.9, 10.2, 9.8, 10.05, 10.15, 9.95], [1, 20, 5, 30, -10, 15, 40, -20]));
  rows = testTable(uneq).rows;
  assert.ok(rows[0].levP < 0.05);
  assert.equal(rows[0]._bold, false);
  assert.equal(rows[1].assume, '등분산을 가정하지 않음');
  assert.equal(rows[1]._bold, true);

  // 화면 HTML: 굵게 표시할 행에만 strong 클래스, 시트 내보내기: 굵은 글씨 표시
  const html = R.renderResultHtml(uneq);
  assert.equal((html.match(/<tr class="strong">/g) || []).length, 1);
  const grid = Ex.resultToGrid(uneq);
  const bolds = grid.rows.filter(r => r.bold);
  assert.equal(bolds.length, 1);
  assert.ok(bolds[0].cells.includes('등분산을 가정하지 않음'));
  assert.ok(testTable(uneq).notes.some(n => n.includes('굵게 표시')));
});

test('runTTestIndependent: 종속변수가 여러 개면 변수마다 따로 판정하고, Levene을 할 수 없으면 강조하지 않는다', () => {
  const Data = require('../src/Data.js');
  const rowsData = [['g', 'a', 'b']];
  const a = [1, 2, 3, 4, 5, 6, 7, 8], b = [10, 10.1, 9.9, 10.2, 9.8, 10.05, 10.15, 9.95];
  const a2 = [2, 3, 4, 5, 6, 7, 8, 9], b2 = [1, 20, 5, 30, -10, 15, 40, -20];
  for (let i = 0; i < 8; i++) rowsData.push([1, a[i], b[i]]);
  for (let i = 0; i < 8; i++) rowsData.push([2, a2[i], b2[i]]);
  const res = C.runTTestIndependent(Data.buildDataset(rowsData, []), { dependent: ['a', 'b'], factor: ['g'], options: {} });
  const rows = res.blocks.find(bl => bl.caption === '독립표본 검정').rows;
  assert.deepEqual(rows.map(r => r._bold), [true, false, false, true]);

  // 모든 값이 집단 안에서 같으면 Levene을 할 수 없다 → 어느 행도 강조하지 않음
  const flat = Data.buildDataset([['g', 'y'], [1, 5], [1, 5], [1, 5], [2, 7], [2, 7], [2, 7]], []);
  const r2 = C.runTTestIndependent(flat, { dependent: ['y'], factor: ['g'], options: {} });
  const rows2 = r2.blocks.find(bl => bl.caption === '독립표본 검정').rows;
  assert.ok(rows2.every(r => !r._bold));
});
