const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const A = require('../src/Stats_RepAnova.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'rmanova.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-10, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const colsOf = c => c.data[0].map((_, j) => c.data.map(r => r[j]));

test('repeatedAnova: pingouin·statsmodels·numpy가 서로 일치한 기준값과 같다 — SS, F, p, Mauchly, ε 3종과 교정 p, η², 쌍별 비교 6가지', () => {
  assert.equal(ref.cases.length, 6);
  ref.cases.forEach(c => {
    const r = A.repeatedAnova(colsOf(c), c.conf);
    const w = c.name;
    r.means.forEach((v, j) => assert.ok(near(v, c.means[j], 1e-12) && near(r.sds[j], c.sds[j], 1e-10), w + ' 평균·SD ' + j));
    assert.ok(near(r.ssTime, c.ss_time, 1e-10) && near(r.ssSubj, c.ss_subj, 1e-10) && near(r.ssErr, c.ss_err, 1e-9) && near(r.ssTotal, c.ss_total, 1e-10), w + ' SS');
    assert.deepEqual([r.dfTime, r.dfErr], [c.df_time, c.df_err]);
    assert.ok(near(r.F, c.F, 1e-9) && near(r.p, c.p, 1e-8), `${w} F ${r.F} vs ${c.F}`);
    assert.ok(near(r.mauchly.W, c.mauchly.W, 1e-8) && near(r.mauchly.chi2, c.mauchly.chi2, 1e-8) && r.mauchly.df === c.mauchly.df, w + ' Mauchly W, χ²');
    assert.ok(near(r.mauchly.p, c.mauchly.p, 1e-7), `${w} Mauchly p ${r.mauchly.p} vs ${c.mauchly.p}`);
    ['gg', 'hf', 'lb'].forEach(key => {
      assert.ok(near(r.eps[key], c.corr[key].eps, 1e-9), `${w} ε ${key}`);
      assert.ok(near(r.corr[key].df1, c.corr[key].df1, 1e-9) && near(r.corr[key].df2, c.corr[key].df2, 1e-9) && near(r.corr[key].p, c.corr[key].p, 1e-7), `${w} 교정 ${key} p ${r.corr[key].p} vs ${c.corr[key].p}`);
    });
    assert.ok(near(r.pEta2, c.peta2, 1e-10) && near(r.gEta2, c.geta2, 1e-10), w + ' η²');
    r.pairs.forEach((p, q) => {
      const e = c.pairs[q];
      assert.deepEqual([p.i, p.j, p.df], [e.i, e.j, e.df]);
      assert.ok(near(p.diff, e.diff, 1e-9) && near(p.se, e.se, 1e-9) && near(p.t, e.t, 1e-9) && near(p.p, e.p, 1e-8) && near(p.pAdj, e.p_bonf, 1e-8), `${w} 쌍 ${q}`);
      assert.ok(near(p.lo, e.lo, 1e-8) && near(p.hi, e.hi, 1e-8), `${w} 쌍 ${q} 구간`);
    });
  });
});

function dsOf(c, extra) {
  const names = c.data[0].map((_, j) => '시점' + (j + 1));
  const rows = [names].concat(c.data.map(r => r.slice()));
  if (extra) extra(rows);
  return { names, ds: Data.buildDataset(rows, names.map(n => ({ name: n, label: '', scale: 'scale', valueLabels: {}, missing: [] }))) };
}

test('runRepeatedAnova: 표 구성, 교정 행 선택(굵게), 쌍별 비교, 평균 선 그림', () => {
  const sph = ref.cases[0], vio = ref.cases[1];
  let { names, ds } = dsOf(sph);
  let res = A.runRepeatedAnova(ds, { vars: names, options: { posthoc: true, plot: true } });
  assert.equal(res.title, '일원 반복측정 분산분석');
  assert.deepEqual(table(res, '처리 요약').rows.map(r => r.n), [20, 0, 4, 20]);
  assert.ok(table(res, '기술통계').rows.every((r, j) => near(r.mean, sph.means[j], 1e-12)));
  const m = table(res, 'Mauchly').rows[0];
  assert.ok(near(m.p, sph.mauchly.p, 1e-7) && near(m.gg, sph.corr.gg.eps, 1e-9) && near(m.hf, sph.corr.hf.eps, 1e-9) && near(m.lb, 1 / 3, 1e-12));
  let eff = table(res, '개체-내 효과 검정').rows;
  assert.deepEqual(eff.map(r => r.name), ['구형성 가정', 'Greenhouse–Geisser', 'Huynh–Feldt', '하한']);
  assert.deepEqual(eff.map(r => !!r._bold), [true, false, false, false]);                       // 구형성 성립 → 첫 행
  assert.ok(near(eff[0].F, sph.F, 1e-9) && near(eff[0].peta, sph.peta2, 1e-10) && near(eff[0].geta, sph.geta2, 1e-10) && near(eff[1].p, sph.corr.gg.p, 1e-7));
  assert.ok(near(eff[1].dfT, sph.corr.gg.df1, 1e-9) && eff[0].dfT === 3 && eff[0].dfE === 57);
  assert.equal(table(res, '쌍별 비교').rows.length, 6);
  assert.equal(res.blocks.find(b => b.type === 'chart').kind, 'line');
  assert.ok(R.renderResultHtml(res).includes('Mauchly') && R.renderResultHtml(res).includes('Bonferroni'));
  assert.ok(Ex.resultToGrid(res).rows.length > 20);
  // 구형성 위반: ε(G-G) ≤ .75 → Greenhouse–Geisser 행
  ({ names, ds } = dsOf(vio));
  res = A.runRepeatedAnova(ds, { vars: names, options: {} });
  eff = table(res, '개체-내 효과 검정').rows;
  assert.deepEqual(eff.map(r => !!r._bold), [false, true, false, false]);
  assert.ok(res.warnings.some(w => w.includes('구형성이 깨졌습니다')));
  assert.ok(!table(res, '쌍별 비교') && !res.blocks.some(b => b.type === 'chart'));
  // 구형성이 깨졌지만 ε(G-G) > .75 인 자료(60명 × 4회)는 Huynh–Feldt 행
  const hfc = ref.cases[5];
  ({ names, ds } = dsOf(hfc));
  const r5 = A.runRepeatedAnova(ds, { vars: names, options: {} });
  assert.ok(hfc.mauchly.p < 0.05 && hfc.corr.gg.eps > 0.75);
  assert.deepEqual(table(r5, '개체-내 효과 검정').rows.map(r => !!r._bold), [false, false, true, false]);
});

test('runRepeatedAnova: 결측 사람은 통째로 제외, 경고, 입력 검증, 계산 불가 경우', () => {
  const c = ref.cases[0];
  const { names, ds } = dsOf(c, rows => { rows[3][1] = ''; rows[8][3] = ''; });
  const res = A.runRepeatedAnova(ds, { vars: names, options: {} });
  assert.deepEqual(table(res, '처리 요약').rows.map(r => r.n), [18, 2, 4, 20]);
  assert.ok(res.warnings.some(w => w.includes('2건을 제외') && w.includes('통째로')));
  assert.throws(() => A.runRepeatedAnova(ds, { vars: names.slice(0, 2), options: {} }), /3회 이상/);
  assert.throws(() => A.runRepeatedAnova(ds, { vars: names, options: { conf: '100' } }), /신뢰수준/);
  const few = Data.buildDataset([['a', 'b', 'c'], [1, 2, 3], [2, 3, 1]], []);
  assert.throws(() => A.runRepeatedAnova(few, { vars: ['a', 'b', 'c'], options: {} }), /3건 이상/);
  const txt = Data.buildDataset([['a', 'b', 'c'], [1, 2, 'x'], [2, 3, 'y'], [3, 1, 'z']], []);
  assert.throws(() => A.runRepeatedAnova(txt, { vars: ['a', 'b', 'c'], options: {} }), /숫자가 아닌 값/);
  // 사람마다 측정 사이 차이가 같음 → 오차 제곱합 0
  const flat = Data.buildDataset([['a', 'b', 'c'], [1, 2, 3], [2, 3, 4], [5, 6, 7], [9, 10, 11]], []);
  const rf = A.runRepeatedAnova(flat, { vars: ['a', 'b', 'c'], options: { posthoc: true } });
  assert.ok(rf.warnings.some(w => w.includes('오차 제곱합이 0')));
  assert.equal(table(rf, '개체-내 효과 검정').rows[0].F, null);
  // 사람 수가 측정 횟수와 같으면 Mauchly 는 계산하되(n ≥ k), 더 적으면 null 안내
  const tiny = Data.buildDataset([['a', 'b', 'c', 'd'], [1, 3, 2, 5], [2, 2, 4, 3], [4, 3, 6, 5]], []);
  const rt = A.runRepeatedAnova(tiny, { vars: ['a', 'b', 'c', 'd'], options: {} });
  assert.equal(table(rt, 'Mauchly').rows[0].W, null);
  assert.ok(table(rt, 'Mauchly').notes.some(n => n.includes('계산하지 못했다')));
});
