const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const Rel = require('../src/Stats_Relation.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'kendall.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-300) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const def = name => ({ name, label: '', scale: 'scale', valueLabels: {}, missing: [] });

test('kendallTauB: scipy kendalltau(타우-b)와 같다 — 타우, 자동 규칙의 유의확률, 점근 유의확률, 정확 유의확률(동순위 없음, n=300까지; 정확 계산은 O(n·불일치쌍)이라 큰 n 은 시험하지 않는다)', () => {
  assert.equal(ref.cases.length, 17);
  ref.cases.forEach(c => {
    const a = Rel.kendallTauB(c.x, c.y), b = Rel.kendallTauB(c.x, c.y, 'asymptotic');
    assert.ok(Math.abs(a.r - c.tau) < 1e-12, `${c.name} tau ${a.r} vs ${c.tau}`);
    assert.ok(near(a.p, c.p_auto, 1e-9), `${c.name} auto p ${a.p} vs ${c.p_auto}`);
    assert.ok(near(b.p, c.p_asym, 1e-9), `${c.name} asym p ${b.p} vs ${c.p_asym}`);
    assert.equal(a.method, !c.ties && c.n <= 33 ? 'exact' : 'asymptotic', c.name);
    if (c.p_exact !== undefined && c.n <= 300) assert.ok(near(Rel.kendallTauB(c.x, c.y, 'exact').p, c.p_exact, 1e-9), `${c.name} exact p`);
  });
});

test('kendallTauB: 성질 — 완전 일치 1, 완전 반대 −1, 한 변수가 상수이면 null, 대칭, 단조 변환에 불변, n < 3 은 null', () => {
  const x = [1, 2, 3, 4, 5, 6, 7];
  assert.equal(Rel.kendallTauB(x, x.map(v => v * 3 + 1)).r, 1);
  assert.equal(Rel.kendallTauB(x, x.map(v => -v)).r, -1);
  assert.equal(Rel.kendallTauB(x, [4, 4, 4, 4, 4, 4, 4]).r, null);
  const y = [3, 1, 4, 1, 5, 9, 2];
  assert.equal(Rel.kendallTauB(x, y).r, Rel.kendallTauB(y, x).r);
  assert.equal(Rel.kendallTauB(x, y).r, Rel.kendallTauB(x, y.map(v => v ** 3)).r);
  assert.equal(Rel.kendallTauB([1, 2], [2, 1]).r, null);
});

function dsOf(c) {
  const rows = [['x', 'y']].concat(c.x.map((v, i) => [v, c.y[i]]));
  return Data.buildDataset(rows, [def('x'), def('y')]);
}

test('runCorrelation: 켄달 타우-b 행렬(계수, 유의확률, N)과 설명이 기준값과 같고 피어슨·스피어만과 함께 나온다', () => {
  const exact = ref.cases[3], ties = ref.cases[10];   // n=20 정확, 정수 점수 n=120 동순위
  [[exact, '정확'], [ties, '정규 근사']].forEach(([c, word]) => {
    const res = Rel.runCorrelation(dsOf(c), { vars: ['x', 'y'], options: { kendall: true } });
    const m = table(res, '켄달 타우-b 상관 행렬');
    assert.ok(m, c.name);
    assert.deepEqual(m.rows.slice(0, 3).map(r => r.stat), ['Kendall 타우-b', '유의확률(양측)', 'N']);
    assert.ok(Math.abs(m.rows[0].c1 - c.tau) < 1e-12 && m.rows[0].c0 === 1);
    assert.ok(near(m.rows[1].c1, c.p_auto, 1e-9) && m.rows[2].c1 === c.n);
    assert.ok(m.notes[0].includes(word === '정확' ? '정확 분포로 구했다' : '정확 분포를 쓴다'), m.notes[0]);
    assert.ok(!table(res, '피어슨'));
  });
  const all = Rel.runCorrelation(dsOf(exact), { vars: ['x', 'y'], options: { pearson: true, spearman: true, kendall: true } });
  assert.deepEqual(all.blocks.filter(b => b.type === 'table').map(b => b.caption.split(' ')[0]), ['피어슨', '스피어만', '켄달']);
  assert.ok(R.renderResultHtml(all).includes('Kendall') && Ex.resultToGrid(all).rows.length > 8);
});

test('runCorrelation: 켄달 — 쌍별 결측, 상수 변수 경고, 통제변수와 함께 쓰면 오류, 방법을 하나도 안 고르면 오류', () => {
  const rows = [['a', 'b', 'c', 'd'], [1, 2, 5, 3], [2, 1, 5, 1], [3, 4, 5, 4], [4, 3, 5, 2], [5, null, 5, 6], [null, 6, 5, 5], [6, 5, 5, 7], [7, 7, 5, 9], [8, 6, 5, 8]];
  const ds = Data.buildDataset(rows, [def('a'), def('b'), def('c'), def('d')]);
  const res = Rel.runCorrelation(ds, { vars: ['a', 'b', 'c'], options: { kendall: true } });
  const m = table(res, '켄달');
  assert.equal(m.rows[2].c1, 7);            // a–b: 두 변수가 모두 있는 사례
  assert.ok(res.warnings.some(w => /값이 모두 같아/.test(w)));
  assert.throws(() => Rel.runCorrelation(ds, { vars: ['a', 'b'], controls: ['d'], options: { kendall: true } }), /켄달 타우-b는 통제변수/);
  assert.throws(() => Rel.runCorrelation(ds, { vars: ['a', 'b'], options: {} }), /켄달 타우-b\)를 하나 이상/);
  assert.doesNotThrow(() => Rel.runCorrelation(ds, { vars: ['a', 'b'], controls: ['d'], options: { pearson: true } }));
});

test('큰 문제: 켄달 타우-b 변수 10개(45쌍) × 5만 건(동순위 많음)이 Node 에서 10초 안에 끝난다', () => {
  let seed = 3;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  const names = Array.from({ length: 10 }, (_, i) => 'v' + i);
  const rows = [names];
  for (let i = 0; i < 50000; i++) { const z = rnd() + rnd(); rows.push(names.map((_, j) => Math.round((z * 0.5 + rnd() + rnd() * (j + 1) * 0.1) * 4))); }
  const ds = Data.buildDataset(rows, names.map(def));
  const t0 = Date.now();
  const res = Rel.runCorrelation(ds, { vars: names, options: { kendall: true } });
  const sec = (Date.now() - t0) / 1000;
  assert.equal(table(res, '켄달').rows.length, 30);
  assert.ok(sec < 10, `켄달 ${sec}초`);
});
