const test = require('node:test');
const assert = require('node:assert/strict');
const Data = require('../src/Data.js');
const A = require('../src/Stats_Assumptions.js');
const C = require('../src/Stats_Compare.js');
const V = require('../src/Stats_Anova.js');
const Rel = require('../src/Stats_Relation.js');
const R = require('../src/Render.js');

/* 재현 가능한 난수 */
function makeRng(seed) {
  let s = seed;
  const u = () => { s = (s * 1664525 + 1013904223) % 4294967296; return (s + 0.5) / 4294967296; };
  const normal = () => Math.sqrt(-2 * Math.log(u())) * Math.cos(2 * Math.PI * u());
  return { u, normal };
}
const rng = makeRng(20261024);
const normalSample = (n, mean = 0, sd = 1) => Array.from({ length: n }, () => mean + sd * rng.normal());
const skewedSample = n => Array.from({ length: n }, () => Math.exp(1.2 * rng.normal()));

const ds = (values, defs) => Data.buildDataset(values, defs || []);
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const texts = res => res.blocks.filter(b => b.type === 'text').map(b => b.text);

/** 열 형태 → 행 형태. cols: { 이름: 값 배열 } */
function toRows(cols) {
  const names = Object.keys(cols);
  const n = cols[names[0]].length;
  const rows = [names];
  for (let i = 0; i < n; i++) rows.push(names.map(k => cols[k][i]));
  return rows;
}

/* ---------- 정규성 판정 ---------- */

test('advNormality_: 상태 분류 (ok, small, reject, largeMild, untestable)', () => {
  assert.equal(A.advNormality_(normalSample(40)).status, 'ok');
  assert.equal(A.advNormality_(normalSample(6)).status, 'small');
  assert.equal(A.advNormality_(skewedSample(30)).status, 'reject');
  assert.equal(A.advNormality_([4, 4, 4, 4]).status, 'untestable');
  assert.equal(A.advNormality_([1, 2]).status, 'untestable');
  assert.equal(A.advNormality_(normalSample(40)).test, 'Shapiro-Wilk');
  // 큰 표본: 약간 이탈해 기각되지만 왜도·첨도가 작으면 largeMild, 크면 reject
  const mild = Array.from({ length: 400 }, (_, i) => { const z = rng.normal(); return z + 0.12 * z * z; });
  const r = A.advNormality_(mild);
  assert.equal(r.test, 'Shapiro-Wilk');
  assert.ok(r.p < 0.05, `p ${r.p}`);
  assert.equal(r.status, 'largeMild');
  assert.equal(A.advNormality_(skewedSample(400)).status, 'reject');
  // 5000건 초과는 Jarque-Bera
  assert.equal(A.advNormality_(normalSample(5200)).test, 'Jarque-Bera');
});

/* ---------- 독립표본 t검정 ---------- */

function twoGroupData(a, b, extra) {
  const g = a.map(() => 1).concat(b.map(() => 2));
  const y = a.concat(b);
  const cols = { 집단: g, 점수: y };
  if (extra) Object.assign(cols, extra);
  return ds(toRows(cols));
}

test('runTTestIndependent: 가정 점검 옵션을 켜면 정규성 기각 시 Mann-Whitney 추천', () => {
  const d = twoGroupData(skewedSample(25), skewedSample(25));
  const res = C.runTTestIndependent(d, { dependent: ['점수'], factor: ['집단'], options: { assumptions: true } });
  const t = table(res, '가정 점검');
  assert.ok(t, '점검표');
  assert.ok(t.rows.some(r => r.test === 'Shapiro-Wilk' && r.verdict === '정규성 기각'));
  assert.ok(t.rows.some(r => r.test === 'Levene'));
  const tx = texts(res);
  assert.ok(tx.some(s => s.startsWith('▶ 추천') && s.includes('Mann-Whitney U 검정')), tx.join('\n'));
  assert.ok(R.renderResultHtml(res).includes('Mann-Whitney U 검정'));
});

test('runTTestIndependent: 옵션을 끄거나 없으면 점검 블록이 없다', () => {
  const d = twoGroupData(skewedSample(25), skewedSample(25));
  const off = C.runTTestIndependent(d, { dependent: ['점수'], factor: ['집단'], options: {} });
  assert.ok(!table(off, '가정 점검'));
  assert.equal(texts(off).length, 0);
});

test('runTTestIndependent: 정규이고 등분산이면 그대로 써도 된다는 문구', () => {
  const d = twoGroupData(normalSample(40, 0, 1), normalSample(40, 0.5, 1));
  const res = C.runTTestIndependent(d, { dependent: ['점수'], factor: ['집단'], options: { assumptions: true } });
  const tx = texts(res);
  assert.equal(tx.length, 1);
  assert.ok(tx[0].includes('그대로 해석해도 됩니다'), tx[0]);
});

test('runTTestIndependent: 정규인데 분산만 다르면 Welch 안내 (비모수는 권하지 않음)', () => {
  const d = twoGroupData(normalSample(60, 0, 1), normalSample(60, 0, 4));
  const res = C.runTTestIndependent(d, { dependent: ['점수'], factor: ['집단'], options: { assumptions: true } });
  const tx = texts(res).join('\n');
  assert.ok(tx.includes('Welch'), tx);
  assert.ok(!tx.includes('Mann-Whitney'), tx);
});

test('runTTestIndependent: 서열 척도 변수는 순위 검정을 추천', () => {
  const y = normalSample(40, 3, 1).map(v => Math.round(v));
  const d = ds(toRows({ 집단: y.map((_, i) => i % 2 + 1), 만족: y }), [{ name: '만족', scale: 'ordinal', label: '', valueLabels: {}, missing: [] }]);
  const res = C.runTTestIndependent(d, { dependent: ['만족'], factor: ['집단'], options: { assumptions: true } });
  assert.ok(texts(res).some(s => s.includes('서열 척도') && s.includes('Mann-Whitney U 검정')), texts(res).join('\n'));
});

test('runTTestIndependent: 큰 표본의 경미한 이탈은 비모수를 권하지 않고 안내만', () => {
  const mild = n => Array.from({ length: n }, () => { const z = rng.normal(); return z + 0.12 * z * z; });
  const d = twoGroupData(mild(300), mild(300));
  const res = C.runTTestIndependent(d, { dependent: ['점수'], factor: ['집단'], options: { assumptions: true } });
  const tx = texts(res).join('\n');
  assert.ok(tx.includes('큰 문제는 없습니다'), tx);
  assert.ok(!tx.includes('Mann-Whitney'), tx);
});

/* ---------- 대응표본·일표본·ANOVA ---------- */

test('runTTestPaired: 차이의 정규성을 보고 Wilcoxon(대응표본)을 추천', () => {
  const x = normalSample(30, 10, 1);
  const y = x.map(v => v - skewedSample(1)[0]);
  const d = ds(toRows({ a: x, b: y }));
  const res = C.runTTestPaired(d, { first: ['a'], second: ['b'], options: { assumptions: true } });
  const t = table(res, '가정 점검');
  assert.equal(t.rows.length, 1);
  assert.ok(t.rows[0].target.includes('차이'));
  assert.ok(!t.rows.some(r => r.test === 'Levene'));
  assert.ok(texts(res).some(s => s.includes('Wilcoxon 부호순위 검정 (대응표본)')), texts(res).join('\n'));
});

test('runTTestOne: 변수의 정규성을 보고 Wilcoxon(일표본)을 추천', () => {
  const d = ds(toRows({ s: skewedSample(30) }));
  const res = C.runTTestOne(d, { vars: ['s'], options: { testValue: '1', assumptions: true } });
  assert.ok(texts(res).some(s => s.includes('Wilcoxon 부호순위 검정 (일표본)')), texts(res).join('\n'));
});

test('runOneWayAnova: 집단별 정규성과 Levene, Kruskal-Wallis 추천', () => {
  const y = skewedSample(20).concat(skewedSample(20), skewedSample(20));
  const g = y.map((_, i) => Math.floor(i / 20) + 1);
  const d = ds(toRows({ g, y }));
  const res = V.runOneWayAnova(d, { dependent: ['y'], factor: ['g'], options: { assumptions: true } });
  const t = table(res, '가정 점검');
  assert.equal(t.rows.filter(r => r.test === 'Shapiro-Wilk').length, 3);
  assert.equal(t.rows.filter(r => r.test === 'Levene').length, 1);
  assert.ok(texts(res).some(s => s.includes('Kruskal-Wallis 검정')), texts(res).join('\n'));
  const off = V.runOneWayAnova(d, { dependent: ['y'], factor: ['g'], options: {} });
  assert.ok(!table(off, '가정 점검'));
});

/* ---------- 상관분석 ---------- */

test('runCorrelation: 피어슨일 때 정규성이 깨지면 스피어만 추천, 스피어만도 켰으면 우선 보라는 문구', () => {
  const x = skewedSample(40);
  const y = x.map(v => v + normalSample(1)[0]);
  const d = ds(toRows({ x, y }));
  const on = Rel.runCorrelation(d, { vars: ['x', 'y'], options: { pearson: true, assumptions: true } });
  const t1 = texts(on).join('\n');
  assert.ok(t1.includes('스피어만 순위상관을 켜서'), t1);
  const both = Rel.runCorrelation(d, { vars: ['x', 'y'], options: { pearson: true, spearman: true, assumptions: true } });
  assert.ok(texts(both).join('\n').includes('스피어만 순위상관 결과를 우선'), texts(both).join('\n'));
  // 피어슨을 쓰지 않으면 점검하지 않는다
  const sp = Rel.runCorrelation(d, { vars: ['x', 'y'], options: { spearman: true, assumptions: true } });
  assert.ok(!table(sp, '가정 점검'));
  const off = Rel.runCorrelation(d, { vars: ['x', 'y'], options: { pearson: true } });
  assert.ok(!table(off, '가정 점검'));
});

/* ---------- 정규성·등분산 검정 결과의 판정 문단 ---------- */

test('runNormality: 변수별 판정 문단 (기각/기각하지 못함)', () => {
  const d = ds(toRows({ 치우침: skewedSample(40), 정규: normalSample(40) }));
  const res = A.runNormality(d, { vars: ['치우침', '정규'] });
  const tx = texts(res);
  assert.equal(tx.length, 2);
  assert.ok(tx[0].includes('"치우침"') && tx[0].includes('정규성이 기각되었습니다') && tx[0].includes('Kruskal-Wallis'), tx[0]);
  assert.ok(tx[1].includes('"정규"') && tx[1].includes('기각하지 못했습니다'), tx[1]);
});

test('runLevene: 등분산 기각 시 Welch와 비모수 안내, 기각하지 못하면 그대로', () => {
  const d1 = twoGroupData(normalSample(50, 0, 1), normalSample(50, 0, 5));
  const r1 = A.runLevene(d1, { dependent: ['점수'], factor: ['집단'], options: {} });
  assert.ok(texts(r1)[0].includes('Welch') && texts(r1)[0].includes('Kruskal-Wallis'), texts(r1)[0]);
  const d2 = twoGroupData(normalSample(50, 0, 1), normalSample(50, 3, 1));
  const r2 = A.runLevene(d2, { dependent: ['점수'], factor: ['집단'], options: {} });
  assert.ok(texts(r2)[0].includes('기각하지 못했습니다'), texts(r2)[0]);
});

test('판정 문단과 추천은 많은 변수에서 10개로 제한된다', () => {
  const cols = {};
  for (let i = 0; i < 14; i++) cols['v' + i] = skewedSample(30);
  const d = ds(toRows(cols));
  const res = A.runNormality(d, { vars: Object.keys(cols) });
  assert.equal(texts(res).length, 11);                       // 10개 + "그 밖의 변수는 표의 유의확률을 보세요"
});
