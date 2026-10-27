const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const A = require('../src/Stats_Assumptions.js');
const R = require('../src/Render.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'assumptions.json'), 'utf8'));

function relErr(a, e) {
  if (a === e) return 0;
  return Math.abs(a - e) / Math.max(Math.abs(e), 1e-300);
}

const ds = (values, defs) => Data.buildDataset(values, defs || []);
const table = (res, caption) => res.blocks.find(b => b.caption && b.caption.startsWith(caption));

/* ---------- Jarque-Bera ---------- */

test('jarqueBera: scipy 기준값 대조', () => {
  for (const c of ref.jb) {
    const r = A.jarqueBera(c.values);
    assert.equal(r.df, 2);
    assert.ok(relErr(r.g1, c.g1) < 1e-9 || Math.abs(r.g1 - c.g1) < 1e-12, `${c.name} g1 ${r.g1} vs ${c.g1}`);
    assert.ok(relErr(r.g2, c.g2) < 1e-9 || Math.abs(r.g2 - c.g2) < 1e-12, `${c.name} g2 ${r.g2} vs ${c.g2}`);
    assert.ok(relErr(r.jb, c.jb) < 1e-8, `${c.name} JB ${r.jb} vs ${c.jb}`);
    assert.ok(relErr(r.p, c.p) < 1e-8, `${c.name} p ${r.p} vs ${c.p}`);
  }
});

test('jarqueBera: 계산할 수 없는 경우는 null', () => {
  assert.equal(A.jarqueBera([1, 2]).jb, null);
  assert.equal(A.jarqueBera([]).jb, null);
  assert.equal(A.jarqueBera([4, 4, 4, 4]).jb, null);
  assert.equal(A.jarqueBera([4, 4, 4, 4]).p, null);
});

test('jarqueBera: 정규성에 가까운 대칭 자료는 p가 크고, 치우친 자료는 작다', () => {
  const sym = [-3, -2, -2, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3];
  assert.ok(A.jarqueBera(sym).p > 0.5);
  const skewed = Array.from({ length: 100 }, (_, i) => Math.exp(i / 20));
  assert.ok(A.jarqueBera(skewed).p < 0.01);
});

test('runNormality: 결과 표 구조, 결측 제외, 경고', () => {
  const d = ds([['a', 'b'], [1, 5], [2, 5], [3, 5], [4, ''], [10, 5], [2, 5], [3, 5]]);
  const res = A.runNormality(d, { vars: ['a', 'b'] });
  const t = table(res, 'Jarque-Bera');
  assert.deepEqual(t.columns.map(c => c.key), ['name', 'n', 'g1', 'g2', 'jb', 'df', 'p']);
  assert.equal(t.rows[0].n, 7);
  assert.equal(t.rows[1].n, 6);
  assert.equal(t.rows[1].jb, null);           // b는 모두 5라서 검정 불가
  assert.equal(t.rows[1].df, null);
  assert.equal(res.warnings.length, 1);
  assert.ok(R.renderResultHtml(res).includes('Jarque-Bera'));
});

test('runNormality: 숫자가 아닌 변수와 빈 선택은 오류', () => {
  assert.throws(() => A.runNormality(ds([['t'], ['x'], ['y'], ['z']]), { vars: ['t'] }), /숫자가 아닌 값/);
  assert.throws(() => A.runNormality(ds([['a'], [1]]), { vars: [] }), /변수를 선택/);
});

/* ---------- Shapiro-Wilk ---------- */

test('shapiroWilk: scipy 기준값 대조 (n=3 ~ 5000)', () => {
  assert.ok(ref.sw.length >= 15);
  for (const c of ref.sw) {
    const r = A.shapiroWilk(c.values);
    assert.equal(r.n, c.values.length);
    // scipy는 계수표 상수와 정규 분위수 근사가 달라 W는 1e-9, p는 W가 1에 가까운 n=5000에서 1e-5 안에서 일치한다
    assert.ok(relErr(r.w, c.W) < 1e-8, `${c.name} W ${r.w} vs ${c.W}`);
    assert.ok(relErr(r.p, c.p) < 1e-5, `${c.name} p ${r.p} vs ${c.p}`);
  }
});

test('shapiroWilk: 범위 밖과 상수는 null, 입력 순서와 무관', () => {
  assert.equal(A.shapiroWilk([1, 2]).w, null);
  assert.equal(A.shapiroWilk([]).p, null);
  assert.equal(A.shapiroWilk([5, 5, 5, 5, 5]).w, null);
  assert.equal(A.shapiroWilk(Array.from({ length: 5001 }, (_, i) => i)).w, null);
  assert.ok(A.shapiroWilk(Array.from({ length: 5000 }, (_, i) => Math.sin(i))).w !== null);
  const x = [3.1, 1.2, 9.9, 4.4, 2.2, 7.1, 5.0, 6.3];
  const y = x.slice().reverse();
  assert.equal(A.shapiroWilk(x).w, A.shapiroWilk(y).w);
});

test('shapiroWilk: 위치·척도 변환에 불변 (큰 평균에서도 W가 같다)', () => {
  const base = [2.1, 3.4, 3.9, 4.0, 4.4, 5.2, 6.8, 9.5, 10.1, 15.0];
  const w0 = A.shapiroWilk(base).w;
  const w1 = A.shapiroWilk(base.map(v => 1e6 + 3 * v)).w;
  assert.ok(relErr(w1, w0) < 1e-9, `${w1} vs ${w0}`);
});

test('runNormality: Shapiro-Wilk 표가 먼저 나오고 5000건 초과는 경고', () => {
  const d = ds([['a'], [1], [2], [4], [8], [9], [3], [5]]);
  const res = A.runNormality(d, { vars: ['a'] });
  const t = table(res, 'Shapiro-Wilk');
  assert.equal(res.blocks[0], t);
  assert.deepEqual(t.columns.map(c => c.key), ['name', 'n', 'w', 'df', 'p']);
  assert.equal(t.rows[0].n, 7);
  assert.equal(t.rows[0].df, 7);
  assert.ok(t.rows[0].w > 0 && t.rows[0].w <= 1);
  assert.ok(R.renderResultHtml(res).includes('Shapiro-Wilk'));

  const big = [['a']].concat(Array.from({ length: 5001 }, (_, i) => [Math.sin(i) * 10]));
  const r2 = A.runNormality(ds(big), { vars: ['a'] });
  assert.equal(table(r2, 'Shapiro-Wilk').rows[0].w, null);
  assert.ok(r2.warnings.some(w => w.includes('5000')));
  assert.ok(table(r2, 'Jarque-Bera').rows[0].jb !== null);
});

/* ---------- Levene ---------- */

test('leveneTest: scipy 기준값 대조 (평균 기준 Levene, 중앙값 기준 Brown-Forsythe)', () => {
  for (const c of ref.levene) {
    const r = A.leveneTest(c.groups, c.center === 'median');
    const N = c.groups.reduce((s, g) => s + g.length, 0);
    assert.equal(r.df1, c.groups.length - 1, c.name);
    assert.equal(r.df2, N - c.groups.length, c.name);
    assert.ok(relErr(r.W, c.W) < 1e-9, `${c.name} (${c.center}) W ${r.W} vs ${c.W}`);
    assert.ok(relErr(r.p, c.p) < 1e-8, `${c.name} (${c.center}) p ${r.p} vs ${c.p}`);
  }
});

test('leveneTest: 집단이 하나이거나 집단 내 변동이 없으면 null', () => {
  assert.equal(A.leveneTest([[1, 2, 3]], false).W, null);
  assert.equal(A.leveneTest([[2, 2, 2], [5, 5, 5]], false).W, null);
});

function twoGroupData() {
  // 집단 1: 퍼짐이 작음, 집단 2: 퍼짐이 큼, 집단 3: 사례가 1개, 일부 결측
  return ds([
    ['y', 'g'],
    [10, 1], [11, 1], [10, 1], [11, 1], [10, 1],
    [0, 2], [20, 2], [5, 2], [30, 2], [12, 2],
    [7, ''], ['', 1], [50, 3]
  ], Data.parseVarDefs([['h'], ['g', '집단', '명목', '1=A, 2=B, 3=C', '']]).defs);
}

test('runLevene: 집단별 표와 검정 표, 값 레이블, 결측 제외', () => {
  const res = A.runLevene(twoGroupData(), { dependent: ['y'], factor: ['g'] });
  const groups = table(res, '집단별 기술통계');
  assert.ok(groups.caption.includes('집단'));
  assert.deepEqual(groups.rows.map(r => r.group), ['A', 'B', 'C']);
  assert.deepEqual(groups.rows.map(r => r.n), [5, 5, 1]);
  const t = table(res, 'Levene');
  assert.equal(t.rows[0].N, 11);                   // 종속·집단이 모두 유효한 사례
  assert.equal(t.rows[0].df1, 2);
  assert.equal(t.rows[0].df2, 8);
  assert.ok(t.rows[0].p < 0.05);                   // 퍼짐이 크게 다르다
  assert.ok(res.warnings.some(w => w.includes('1개뿐')));
});

test('runLevene: 중앙값 옵션은 다른 통계량을 내고 각주에 표시', () => {
  const mean = A.runLevene(twoGroupData(), { dependent: ['y'], factor: ['g'] });
  const med = A.runLevene(twoGroupData(), { dependent: ['y'], factor: ['g'], options: { median: true } });
  const wMean = table(mean, 'Levene').rows[0].W;
  const wMed = table(med, 'Levene').rows[0].W;
  assert.notEqual(wMean, wMed);
  assert.ok(table(med, 'Levene').notes.some(n => n.includes('중앙값')));
  assert.ok(table(mean, 'Levene').notes.some(n => n.includes('평균')));
});

test('runLevene: 입력 오류', () => {
  const d = twoGroupData();
  assert.throws(() => A.runLevene(d, { dependent: ['y'], factor: [] }), /변수를 선택/);
  assert.throws(() => A.runLevene(d, { dependent: ['y'], factor: ['y'] }), /같은 변수/);
  assert.throws(() => A.runLevene(d, { dependent: [], factor: ['g'] }), /변수를 선택/);
  assert.throws(() => A.runLevene(d, { dependent: ['y'], factor: ['g', 'y'] }), /하나만/);
  assert.throws(() => A.runLevene(ds([['y', 'g'], ['x', 1]]), { dependent: ['y'], factor: ['g'] }), /숫자가 아닌 값/);
});

test('runLevene: 집단이 너무 많으면 오류', () => {
  const rows = [['y', 'g']];
  for (let i = 0; i < 60; i++) rows.push([i, i], [i + 0.5, i]);
  assert.throws(() => A.runLevene(ds(rows), { dependent: ['y'], factor: ['g'] }), /60개/);
});

test('runLevene: 종속변수 여러 개, 끝까지 렌더링', () => {
  const d = ds([['y1', 'y2', 'g'], [1, 5, 'a'], [2, 7, 'a'], [3, 4, 'a'], [2, 9, 'b'], [8, 1, 'b'], [4, 6, 'b'], [9, 2, 'b']]);
  const res = A.runLevene(d, { dependent: ['y1', 'y2'], factor: ['g'] });
  assert.equal(table(res, 'Levene').rows.length, 2);
  assert.equal(table(res, '집단별').rows.length, 4);
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('Levene의 등분산 검정') && html.includes('유의확률'));
});
