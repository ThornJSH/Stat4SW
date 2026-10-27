const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const B = require('../src/Stats_Binomial.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'binom.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-14) => Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const def = (name, scale, vl) => ({ name, label: '', scale: scale || 'nominal', valueLabels: vl || {}, missing: [] });

test('binomTest: scipy binomtest(양측·단측)와 같다 — 20가지 (k, n, p0) 조합', () => {
  assert.equal(ref.cases.length, 20);
  ref.cases.forEach(c => {
    const m = `${c.k}/${c.n} p0=${c.p0}`;
    ['two', 'greater', 'less'].forEach(alt => {
      const p = B.binomTest(c.k, c.n, c.p0, alt);
      assert.ok(near(p, c[alt], 1e-9, 1e-14), `${m} ${alt}: ${p} vs ${c[alt]}`);
    });
  });
});

test('신뢰구간: Clopper–Pearson(정확)과 Wilson이 scipy proportion_ci·statsmodels proportion_confint와 같다', () => {
  ref.cases.forEach(c => {
    const ex = B.binomExactCi(c.k, c.n, c.conf), wi = B.binomWilsonCi(c.k, c.n, c.conf);
    const m = `${c.k}/${c.n} conf=${c.conf}`;
    assert.ok(near(ex[0], c.exact[0], 1e-9, 1e-12) && near(ex[1], c.exact[1], 1e-9, 1e-12), `${m} exact ${ex} vs ${c.exact}`);
    assert.ok(near(wi[0], c.wilson[0], 1e-10, 1e-12) && near(wi[1], c.wilson[1], 1e-10, 1e-12), `${m} wilson ${wi} vs ${c.wilson}`);
  });
});

test('binomTest: 성질 — p0 = 0.5 의 양측은 단측의 2배, 관측이 기대값이면 1, 꼬리 확률의 합은 1', () => {
  assert.ok(near(B.binomTest(27, 40, 0.5, 'two'), 2 * B.binomTest(27, 40, 0.5, 'greater'), 1e-12));
  assert.equal(B.binomTest(5, 10, 0.5, 'two'), 1);
  [[10, 0.3, 4], [200, 0.15, 33], [25, 0.8, 3]].forEach(([n, p, k]) => {
    assert.ok(near(B.bnCdf(k - 1, n, p) + B.bnSf(k, n, p), 1, 1e-12));
    assert.ok(B.binomTest(k, n, p, 'two') <= 1);
  });
});

function dsOfCounts(x, n, extra) {
  const rows = [['y']];
  for (let i = 0; i < n; i++) rows.push([i < x ? 1 : 0]);
  if (extra) extra(rows);
  return Data.buildDataset(rows, [def('y', 'nominal', { 1: '합격', 0: '불합격' })]);
}

test('runBinomial: 표 세 개의 값(빈도, 검정, 신뢰구간)이 기준값과 같다 (양측·단측, 99% 구간)', () => {
  const c = ref.cases[6];   // 45/48, p0 .8, 99%
  ['two', 'greater', 'less'].forEach(alt => {
    const res = B.runBinomial(dsOfCounts(c.k, c.n), { vars: ['y'], options: { success: '1', p0: String(c.p0), tails: alt, conf: String(c.conf) } });
    const f = table(res, '성공·실패'), t = table(res, '이항 검정'), ci = table(res, '성공 비율의');
    assert.deepEqual([f.rows[0].x, f.rows[0].y, f.rows[0].n, f.rows[0].miss], [45, 3, 48, 0]);
    assert.equal(f.rows[0].value, '1 (합격)');
    assert.ok(near(t.rows[0].p, c[alt], 1e-9, 1e-14) && t.rows[0].ph === 45 / 48 && t.rows[0].p0 === 0.8);
    assert.ok(near(ci.rows[0].lo, c.exact[0], 1e-9) && near(ci.rows[0].whi, c.wilson[1], 1e-10));
    assert.ok(ci.caption.includes('99%'));
  });
  const two = B.runBinomial(dsOfCounts(c.k, c.n), { vars: ['y'], options: { p0: String(c.p0), conf: '99' } });
  assert.ok(table(two, '이항 검정').caption.includes('양측') && table(two, '이항 검정').columns[5].label.includes('양측'));
  assert.ok(R.renderResultHtml(two).includes('Clopper') && Ex.resultToGrid(two).rows.length > 5);
});

test('runBinomial: 여러 변수(변수별 결측 제외), 값이 3가지 이상이면 경고, 글자 값과 성공값 지정', () => {
  const rows = [['a', 'c'], [1, 1], [0, 2], [1, 3], [1, 1], [null, 1], [0, 2], [1, 3]];
  const ds = Data.buildDataset(rows, [def('a'), def('c')]);
  const res = B.runBinomial(ds, { vars: ['a', 'c'], options: { success: '1', p0: '0.5' } });
  const f = table(res, '성공·실패');
  assert.deepEqual(f.rows.map(r => [r.x, r.y, r.n, r.miss]), [[4, 2, 6, 1], [3, 4, 7, 0]]);
  assert.ok(res.warnings.some(w => /"c".*3가지/.test(w)) && !res.warnings.some(w => /"a".*가지/.test(w)));
  assert.ok(near(table(res, '이항 검정').rows[0].p, B.binomTest(4, 6, 0.5, 'two'), 1e-14));
  const txt = Data.buildDataset([['b'], ['예'], ['아니오'], ['예'], ['예'], ['아니오']], [def('b')]);
  const r2 = B.runBinomial(txt, { vars: ['b'], options: { success: '예', p0: '0.5', tails: 'greater' } });
  assert.deepEqual([table(r2, '성공·실패').rows[0].x, table(r2, '성공·실패').rows[0].n], [3, 5]);
  assert.ok(near(table(r2, '이항 검정').rows[0].p, B.binomTest(3, 5, 0.5, 'greater'), 1e-14));
  assert.ok(table(r2, '이항 검정').columns[5].label.includes('크다'));
});

test('runBinomial: 오류 안내 — 성공값이 없음, 값이 없음, 검정 비율·신뢰수준 범위, 변수 수', () => {
  const ds = dsOfCounts(3, 10);
  assert.throws(() => B.runBinomial(ds, { vars: ['y'], options: { success: '9' } }), /성공으로 볼 값 "9".*변수에 있는 값: 1, 0/);
  assert.throws(() => B.runBinomial(ds, { vars: ['y'], options: { p0: '1' } }), /검정 비율/);
  assert.throws(() => B.runBinomial(ds, { vars: ['y'], options: { p0: '0' } }), /검정 비율/);
  assert.throws(() => B.runBinomial(ds, { vars: ['y'], options: { p0: 'abc' } }), /검정 비율/);
  assert.throws(() => B.runBinomial(ds, { vars: ['y'], options: { conf: '100' } }), /신뢰수준/);
  assert.throws(() => B.runBinomial(Data.buildDataset([['y'], [null], [null]], [def('y')]), { vars: ['y'], options: {} }), /값이 하나도 없습니다/);
  assert.throws(() => B.runBinomial(ds, { vars: [], options: {} }), /변수를 선택/);
});

test('큰 문제: 5만 건 × 변수 20개 이항 검정이 Node 에서 3초 안에 끝난다', () => {
  const names = Array.from({ length: 20 }, (_, i) => 'v' + i);
  let seed = 11;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  const rows = [names];
  for (let i = 0; i < 50000; i++) rows.push(names.map((_, j) => (rnd() < 0.4 + j * 0.01 ? 1 : 0)));
  const ds = Data.buildDataset(rows, names.map(n => def(n)));
  const t0 = Date.now();
  const res = B.runBinomial(ds, { vars: names, options: { success: '1', p0: '0.37', tails: 'two' } });
  const sec = (Date.now() - t0) / 1000;
  assert.equal(table(res, '이항 검정').rows.length, 20);
  assert.ok(sec < 3, `이항 검정 ${sec}초`);
});
