const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const N = require('../src/Stats_Nonparam.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'friedman.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-10, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const colsOf = c => c.data[0].map((_, j) => c.data.map(r => r[j]));

test('friedman: scipy 기준값(χ², 유의확률, Kendall W, 평균 순위, 순위 합, 동순위 보정) 6가지', () => {
  assert.equal(ref.cases.length, 6);
  ref.cases.forEach(c => {
    const f = N.friedman(colsOf(c));
    assert.equal(f.n, c.n); assert.equal(f.k, c.k); assert.equal(f.df, c.k - 1);
    f.meanRanks.forEach((v, j) => assert.ok(near(v, c.mean_ranks[j], 1e-12), `${c.name} 평균 순위 ${j}`));
    f.rankSums.forEach((v, j) => assert.ok(near(v, c.rank_sums[j], 1e-12)));
    assert.equal(f.tieTotal, c.tie_total, c.name + ' 동순위 항');
    if (c.chi2 === null) { assert.equal(f.chi2, null); assert.equal(f.p, null); assert.equal(f.w, null); return; }
    assert.ok(near(f.chi2, c.chi2, 1e-10), `${c.name} χ² ${f.chi2} vs ${c.chi2}`);
    assert.ok(near(f.p, c.p, 1e-8), `${c.name} p ${f.p} vs ${c.p}`);
    assert.ok(near(f.w, c.w, 1e-10), c.name + ' W');
  });
});

test('friedmanPairs: Dunn–Bonferroni 사후비교가 독립 계산(z, 유의확률, 보정 유의확률)과 같다', () => {
  ref.cases.filter(c => c.chi2 !== null).forEach(c => {
    const pairs = N.friedmanPairs(N.friedman(colsOf(c)));
    assert.equal(pairs.length, c.k * (c.k - 1) / 2);
    pairs.forEach((p, q) => {
      const e = c.pairs[q];
      assert.deepEqual([p.i, p.j], [e.i, e.j]);
      assert.ok(near(p.diff, e.diff, 1e-10) && near(p.z, e.z, 1e-9), `${c.name} 쌍 ${q} z ${p.z} vs ${e.z}`);
      assert.ok(near(p.p, e.p, 1e-8) && near(p.pAdj, e.p_adj, 1e-8), `${c.name} 쌍 ${q} p`);
    });
  });
  const none = N.friedmanPairs(N.friedman(colsOf(ref.cases[5])));
  assert.ok(none.every(p => p.z === null && p.pAdj === null));
});

test('0/1 자료의 Friedman χ² 는 Cochran Q 와 같다 (도움말과 도우미가 이 사실에 기대고 있다)', () => {
  const c = ref.cases[4];
  assert.ok(c.cochran_q !== null);
  assert.ok(near(N.friedman(colsOf(c)).chi2, c.cochran_q, 1e-10));
});

function dsOf(c, extra) {
  const names = c.data[0].map((_, j) => '측정' + (j + 1));
  const rows = [names].concat(c.data.map(r => r.slice()));
  if (extra) extra(rows);
  return { names, ds: Data.buildDataset(rows, names.map(n => ({ name: n, label: '', scale: 'scale', valueLabels: {}, missing: [] }))) };
}

test('runFriedman: 표 구성(순위, 검정 통계량, 사후비교, 상자그림)과 값이 기준값과 같다', () => {
  const c = ref.cases[1];                                              // 리커트, 동순위 많음
  const { names, ds } = dsOf(c);
  const res = N.runFriedman(ds, { vars: names, options: { posthoc: true, effect: true, boxplot: true } });
  assert.equal(res.title, 'Friedman 검정');
  assert.deepEqual(table(res, '처리 요약').rows.map(r => r.n), [30, 0, 3, 30]);
  const rk = table(res, '순위').rows;
  rk.forEach((r, j) => assert.ok(near(r.meanRank, c.mean_ranks[j], 1e-12) && near(r.median, c.medians[j], 1e-12) && r.n === 30));
  const st = table(res, 'Friedman 검정 통계량').rows[0];
  assert.ok(near(st.chi2, c.chi2, 1e-10) && near(st.p, c.p, 1e-8) && near(st.w, c.w, 1e-10) && st.df === 2);
  const ph = table(res, '사후비교').rows;
  assert.equal(ph.length, 3);
  ph.forEach((r, q) => { assert.ok(near(r.z, c.pairs[q].z, 1e-9) && near(r.pAdj, c.pairs[q].p_adj, 1e-8)); assert.equal(r.sig, c.pairs[q].p_adj < 0.05 ? '*' : ''); });
  assert.equal(ph[0].gi, '측정1');
  assert.equal(res.blocks.find(b => b.type === 'chart').kind, 'box');
  assert.ok(R.renderResultHtml(res).includes('Friedman') && R.renderResultHtml(res).includes('Dunn'));
  assert.ok(Ex.resultToGrid(res).rows.length > 15);
  // 옵션을 끄면 사후비교·W·그림이 없다
  const off = N.runFriedman(ds, { vars: names, options: {} });
  assert.ok(!table(off, '사후비교') && !off.blocks.some(b => b.type === 'chart') && table(off, 'Friedman 검정 통계량').columns.length === 4);
});

test('runFriedman: 결측 사람은 통째로 제외, 경고, 입력 검증', () => {
  const c = ref.cases[0];
  const { names, ds } = dsOf(c, rows => { rows[3][1] = ''; rows[8][3] = ''; });
  const res = N.runFriedman(ds, { vars: names, options: {} });
  assert.deepEqual(table(res, '처리 요약').rows.map(r => r.n), [18, 2, 4, 20]);
  assert.ok(res.warnings.some(w => w.includes('2건을 제외') && w.includes('통째로')));
  assert.throws(() => N.runFriedman(ds, { vars: names.slice(0, 2), options: {} }), /3회 이상/);
  const same = dsOf(ref.cases[5]);
  const r2 = N.runFriedman(same.ds, { vars: same.names, options: { posthoc: true } });
  assert.ok(r2.warnings.some(w => w.includes('전부 같아')) && !table(r2, '사후비교'));
  const few = dsOf({ data: [[1, 2, 3], [2, 3, 1], [3, 1, 2], [1, 3, 2]] });
  assert.ok(N.runFriedman(few.ds, { vars: few.names, options: {} }).warnings.some(w => w.includes('근사가 부정확')));
  const one = Data.buildDataset([['a', 'b', 'c'], [1, 2, 3]], []);
  assert.throws(() => N.runFriedman(one, { vars: ['a', 'b', 'c'], options: {} }), /2건 이상/);
  const txt = Data.buildDataset([['a', 'b', 'c'], [1, 2, 'x'], [2, 3, 'y']], []);
  assert.throws(() => N.runFriedman(txt, { vars: ['a', 'b', 'c'], options: {} }), /숫자가 아닌 값/);
});
