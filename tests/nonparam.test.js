const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const N = require('../src/Stats_Nonparam.js');
const R = require('../src/Render.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'nonparam.json'), 'utf8'));

function relErr(a, e) {
  if (a === e) return 0;
  return Math.abs(a - e) / Math.max(Math.abs(e), 1e-300);
}
const ds = (values, defs) => Data.buildDataset(values, defs || []);
const table = (res, caption) => res.blocks.find(b => b.caption && b.caption.startsWith(caption));

/* ---------- 순위 ---------- */

test('npRanks_: 평균순위와 동순위 항', () => {
  const r = N.npRanks_([10, 20, 20, 30, 20, 5]);
  assert.deepEqual(r.ranks, [2, 4, 4, 6, 4, 1]);
  assert.equal(r.tieSum, 3 * 3 * 3 - 3);          // 20이 3개
  assert.equal(N.npRanks_([3, 1, 2]).tieSum, 0);
});

/* ---------- Mann-Whitney ---------- */

test('mannWhitney: scipy 기준값 대조 (U, 점근 p, 정확 p)', () => {
  for (const c of ref.mw) {
    const r = N.mannWhitney(c.a, c.b);
    assert.equal(r.u, Math.min(c.u1, c.a.length * c.b.length - c.u1), c.name);
    assert.equal(r.u1, c.u1, c.name);
    assert.ok(relErr(r.p, c.pAsym) < 1e-10, `${c.name} p ${r.p} vs ${c.pAsym}`);
    assert.equal(r.ties, c.ties, c.name);
    if (c.pExact === null) assert.equal(r.pExact, null, c.name);
    else assert.ok(relErr(r.pExact, c.pExact) < 1e-10, `${c.name} exact ${r.pExact} vs ${c.pExact}`);
  }
});

test('mannWhitney: 순위합, W, Z의 부호와 효과크기, 불가능한 경우', () => {
  const r = N.mannWhitney([1, 2, 3], [4, 5, 6, 7]);
  assert.equal(r.rankSum1, 6);
  assert.equal(r.rankSum2, 22);
  assert.equal(r.u, 0);
  assert.equal(r.w, 6);                            // 사례가 적은 집단(집단1)의 순위합
  assert.ok(r.z < 0);
  assert.ok(Math.abs(r.effect - Math.abs(r.z) / Math.sqrt(7)) < 1e-15);
  assert.equal(N.mannWhitney([5, 5], [5, 5, 5]).z, null);     // 모두 같은 값
  assert.equal(N.mannWhitney([], [1]).u, null);
  // 집단 순서를 바꿔도 U(작은 값)와 p는 같다
  const a = [1.1, 5.2, 3.3, 8.4], b = [2.5, 9.1, 7.7, 6.6, 4.0];
  const x = N.mannWhitney(a, b), y = N.mannWhitney(b, a);
  assert.equal(x.u, y.u);
  assert.ok(relErr(x.p, y.p) < 1e-14);
});

test('mwExactCdf_: 전수 계산과 일치', () => {
  // m=3, n=4: 가능한 U 분포 (총 35가지)
  const counts = new Array(13).fill(0);
  const items = [1, 2, 3, 4, 5, 6, 7];
  for (let a = 0; a < 7; a++) for (let b = a + 1; b < 7; b++) for (let c = b + 1; c < 7; c++) {
    const sum = items[a] + items[b] + items[c];
    counts[sum - 6]++;
  }
  let cum = 0;
  for (let u = 0; u <= 12; u++) {
    cum += counts[u];
    assert.ok(relErr(N.mwExactCdf_(3, 4, u), cum / 35) < 1e-14, `u=${u}`);
  }
});

/* ---------- Wilcoxon ---------- */

test('signedRank: scipy 기준값 대조 (통계량, Z, 점근 p, 정확 p)', () => {
  for (const c of ref.wilcoxon) {
    const r = N.signedRank(c.diffs);
    assert.equal(Math.min(r.wPlus, r.wMinus), c.stat, c.name);
    assert.ok(Math.abs(Math.abs(r.z) - Math.abs(c.z)) < 1e-12, `${c.name} z ${r.z} vs ${c.z}`);
    assert.ok(relErr(r.p, c.pAsym) < 1e-10, `${c.name} p ${r.p} vs ${c.pAsym}`);
    if (c.pExact === null) assert.equal(r.pExact, null, c.name);
    else assert.ok(relErr(r.pExact, c.pExact) < 1e-10, `${c.name} exact ${r.pExact} vs ${c.pExact}`);
  }
});

test('signedRank: 0 차이는 순위에서 빠지고, Z의 부호는 양의 순위합 기준', () => {
  const r = N.signedRank([0, 0, 1, 2, 3, -1.5]);
  assert.equal(r.nZero, 2);
  assert.equal(r.n, 4);
  assert.equal(r.nPos, 3);
  assert.equal(r.nNeg, 1);
  assert.equal(r.wPlus + r.wMinus, 10);            // 순위 1..4의 합
  assert.ok(r.z > 0);
  assert.ok(N.signedRank([-1, -2, -3, -4, -5]).z < 0);
  assert.equal(N.signedRank([0, 0, 0]).z, null);
  assert.equal(N.signedRank([]).n, 0);
});

/* ---------- Kruskal-Wallis, Dunn ---------- */

test('kruskalWallis: scipy 기준값 대조 (H, p)', () => {
  for (const c of ref.kw) {
    const r = N.kruskalWallis(c.groups);
    assert.equal(r.df, c.groups.length - 1, c.name);
    assert.ok(relErr(r.h, c.H) < 1e-10, `${c.name} H ${r.h} vs ${c.H}`);
    assert.ok(relErr(r.p, c.p) < 1e-10, `${c.name} p ${r.p} vs ${c.p}`);
  }
});

test('dunnTest: 독립 구현(numpy) 기준값 대조 (평균 순위 차이, Z, 보정 p)', () => {
  for (const c of ref.kw) {
    const got = N.dunnTest(c.groups);
    assert.equal(got.length, c.dunn.length, c.name);
    got.forEach((g, i) => {
      const e = c.dunn[i];
      assert.equal(g.i, e.i); assert.equal(g.j, e.j);
      assert.ok(Math.abs(g.diff - e.diff) < 1e-10, `${c.name} diff`);
      assert.ok(relErr(g.z, e.z) < 1e-10, `${c.name} z ${g.z} vs ${e.z}`);
      assert.ok(relErr(g.p, e.p) < 1e-9, `${c.name} p`);
      assert.ok(relErr(g.pAdj, e.pAdj) < 1e-9, `${c.name} pAdj`);
    });
  }
});

test('kruskalWallis: 두 집단이고 동순위가 없으면 H = Z²이므로 Mann-Whitney 점근 p와 같다, η²는 0 이상', () => {
  const a = [1.2, 3.4, 2.2, 8.1, 5.5], b = [4.4, 6.1, 9.3, 7.7, 10.2, 12.0];
  const kw = N.kruskalWallis([a, b]);
  const mw = N.mannWhitney(a, b);
  assert.ok(relErr(kw.p, mw.p) < 1e-9, `${kw.p} vs ${mw.p}`);
  assert.ok(kw.eta2 >= 0);
  assert.equal(N.kruskalWallis([[3, 3], [3, 3, 3]]).h, null);        // 모두 같은 값
});

/* ---------- 결과 객체 ---------- */

const sample = () => {
  const rows = [['g', 'h', 'x', 'y']];
  const x = [5, 7, 8, 9, 12, 14, 15, 17, 20, 22, 6, 9, 11, 16, 19];
  const y = [6, 9, 7, 12, 11, 15, 19, 16, 25, 23, 8, 10, 13, 14, 18];
  for (let i = 0; i < 15; i++) rows.push([i % 3 + 1, i % 2 + 1, x[i], y[i]]);
  return ds(rows);
};

test('runMannWhitney: 표 구조, 집단 수 오류, 결측 제외', () => {
  const d = sample();
  const res = N.runMannWhitney(d, { dependent: ['x', 'y'], factor: ['h'], options: { effect: true } });
  const rank = table(res, '순위');
  const stat = table(res, 'Mann-Whitney U 검정 통계량');
  assert.deepEqual(stat.columns.map(c => c.key), ['name', 'u', 'w', 'z', 'p', 'pExact', 'effect']);
  assert.equal(stat.rows.length, 2);
  assert.equal(rank.rows.filter(r => r._kind === 'total').length, 2);
  assert.ok(R.renderResultHtml(res).includes('Mann-Whitney U'));
  assert.throws(() => N.runMannWhitney(d, { dependent: ['x'], factor: ['g'], options: {} }), /정확히 2개/);
  assert.throws(() => N.runMannWhitney(d, { dependent: ['x'], factor: ['x'], options: {} }), /같은 변수/);
  const noEffect = N.runMannWhitney(d, { dependent: ['x'], factor: ['h'], options: {} });
  assert.ok(!table(noEffect, 'Mann-Whitney U 검정 통계량').columns.some(c => c.key === 'effect'));
});

test('runWilcoxonPaired: 순위표(음·양·동순위)와 Z, 같은 변수와 변수 미선택은 오류', () => {
  const d = ds([['a', 'b'], [1, 2], [5, 3], [4, 4], [7, 9], [8, 3], [6, 1], [2, 8], ['', 5]]);
  const res = N.runWilcoxonPaired(d, { first: ['a'], second: ['b'], options: { effect: true } });
  const rank = table(res, '순위');
  const kinds = rank.rows.map(r => r.kind);
  assert.ok(kinds[0].startsWith('음의 순위') && kinds[1].startsWith('양의 순위') && kinds[2].startsWith('동순위') && kinds[3] === '합계');
  assert.equal(rank.rows[0].n, 3);                  // 1<2, 7<9, 2<8
  assert.equal(rank.rows[1].n, 3);                  // 5>3, 8>3, 6>1
  assert.equal(rank.rows[2].n, 1);                  // 4=4
  assert.equal(rank.rows[3].n, 7);                  // b가 비어 있는 쌍은 제외
  assert.ok(table(res, 'Wilcoxon').rows[0].z !== null);
  assert.throws(() => N.runWilcoxonPaired(d, { first: ['a'], second: ['a'], options: {} }), /같은 변수/);
  assert.throws(() => N.runWilcoxonPaired(d, { first: [], second: ['a'], options: {} }), /변수/);
});

test('runWilcoxonOne: 검정값 옵션, 검정값과 같은 값은 동순위 행', () => {
  const d = ds([['s'], [10], [12], [15], [20], [9], [30], [18]]);
  const res = N.runWilcoxonOne(d, { vars: ['s'], options: { testValue: '15' } });
  const rank = table(res, '순위');
  assert.equal(rank.rows[2].n, 1);                  // 15
  assert.equal(rank.rows[0].n, 3);                  // 10, 12, 9
  assert.equal(rank.rows[1].n, 3);                  // 20, 30, 18
  assert.ok(table(res, 'Wilcoxon').caption.includes('검정값 = 15'));
  assert.throws(() => N.runWilcoxonOne(d, { vars: ['s'], options: { testValue: 'abc' } }), /숫자/);
});

test('runKruskalWallis: 표 구조, Dunn 옵션, 집단 1개는 오류', () => {
  const d = sample();
  const res = N.runKruskalWallis(d, { dependent: ['x'], factor: ['g'], options: { dunn: true, effect: true, boxplot: true } });
  const stat = table(res, 'Kruskal-Wallis 검정 통계량');
  assert.deepEqual(stat.columns.map(c => c.key), ['name', 'h', 'df', 'p', 'eta2']);
  const dunn = table(res, '다중비교: Dunn');
  assert.equal(dunn.rows.length, 3);                // 3집단 → 3쌍
  assert.ok(res.blocks.some(b => b.type === 'chart'));
  const off = N.runKruskalWallis(d, { dependent: ['x'], factor: ['g'], options: {} });
  assert.ok(!table(off, '다중비교: Dunn'));
  assert.throws(() => N.runKruskalWallis(ds([['g', 'x'], [1, 1], [1, 2]]), { dependent: ['x'], factor: ['g'], options: {} }), /2개 이상/);
});
