const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const A = require('../src/Stats_RepAnova.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'mixed.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-12) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const colsOf = c => c.Y[0].map((_, j) => c.Y.map(r => r[j]));

function dsOf(c, valueLabels, extra) {
  const names = c.Y[0].map((_, j) => 't' + (j + 1));
  const rows = [names.concat(['grp'])].concat(c.Y.map((r, i) => r.concat([c.G[i]])));
  if (extra) extra(rows);
  const defs = names.map(n => ({ name: n, label: '', scale: 'scale', valueLabels: {}, missing: [] }))
    .concat([{ name: 'grp', label: '집단', scale: 'nominal', valueLabels: valueLabels || {}, missing: [] }]);
  return { names, ds: Data.buildDataset(rows, defs) };
}
const g0 = c => c.G.map(v => v - 1);

test('mixedAnova: statsmodels(유형 III)·pingouin·numpy가 서로 일치한 기준값과 같다 — 개체-간, 시간, 시간×집단, 오차, ε, Mauchly, 교정 p (5가지 사례)', () => {
  assert.equal(ref.cases.length, 5);
  ref.cases.forEach(c => {
    const r = A.mixedAnova(colsOf(c), g0(c), c.g);
    const w = c.name;
    assert.deepEqual(r.nh, c.nh);
    const b = c.between;
    assert.ok(near(r.between.intercept.ss, b.intercept.ss, 1e-9) && near(r.between.intercept.F, b.intercept.F, 1e-8) && near(r.between.intercept.p, b.intercept.p, 1e-6, 1e-12), w + ' 절편');
    assert.ok(near(r.between.group.ss, b.group.ss, 1e-9) && r.between.group.df === b.group.df && near(r.between.group.F, b.group.F, 1e-8) && near(r.between.group.p, b.group.p, 1e-6, 1e-12) && near(r.between.group.peta, b.group.peta, 1e-9), w + ' 집단');
    assert.ok(near(r.between.err.ss, b.err.ss, 1e-9) && r.between.err.df === b.err.df && near(r.between.err.ms, b.err.ms, 1e-9), w + ' 개체-간 오차');
    assert.ok(near(r.errw.ss, c.errw.ss, 1e-9) && r.errw.df === c.errw.df, w + ' 개체-내 오차');
    ['time', 'inter'].forEach(key => {
      const e = c[key], x = r[key];
      assert.ok(near(x.ss, e.ss, 1e-8, 1e-10) && x.df === e.df && near(x.F, e.F, 1e-8) && near(x.p, e.p, 1e-6, 1e-12) && near(x.peta, e.peta, 1e-9), `${w} ${key} 효과`);
      if (c.k >= 3) {
        ['gg', 'hf', 'lb'].forEach(corr => {
          assert.ok(near(x.corr[corr].eps, e.corr[corr].eps, 1e-8) && near(x.corr[corr].df1, e.corr[corr].df1, 1e-8) && near(x.corr[corr].df2, e.corr[corr].df2, 1e-8) && near(x.corr[corr].p, e.corr[corr].p, 1e-6, 1e-12), `${w} ${key} ${corr}`);
        });
      }
    });
    if (c.k >= 3) {
      assert.ok(near(r.mauchly.W, c.mauchly.W, 1e-8) && near(r.mauchly.chi2, c.mauchly.chi2, 1e-8) && r.mauchly.df === c.mauchly.df && near(r.mauchly.p, c.mauchly.p, 1e-7, 1e-12), w + ' Mauchly');
    } else {
      assert.equal(r.mauchly, null);
      assert.equal(r.eps, null);
    }
    r.cells.forEach((row, h) => row.forEach((cell, j) => {
      const e = c.cells[h][j];
      assert.ok(cell.n === e.n && near(cell.mean, e.mean, 1e-12) && near(cell.sd, e.sd, 1e-10), `${w} 셀 ${h},${j}`);
    }));
  });
});

test('mixedAnova: 집단이 하나이면 일원 반복측정 분산분석(repeatedAnova)과 시간 효과, ε, Mauchly 가 같다', () => {
  const c = ref.cases[0];
  const cols = colsOf(c);
  const one = A.mixedAnova(cols, new Array(c.N).fill(0), 1);
  const rm = A.repeatedAnova(cols, 95);
  assert.ok(near(one.time.ss, rm.ssTime, 1e-9) && near(one.time.F, rm.F, 1e-9) && near(one.time.p, rm.p, 1e-9, 1e-12), '시간');
  assert.ok(near(one.eps.gg, rm.eps.gg, 1e-9) && near(one.eps.hf, rm.eps.hf, 1e-9) && near(one.eps.lb, rm.eps.lb, 1e-12));
  assert.ok(near(one.mauchly.W, rm.mauchly.W, 1e-9) && near(one.mauchly.p, rm.mauchly.p, 1e-9, 1e-12));
  assert.ok(near(one.errw.ss, rm.ssErr, 1e-9));
});

test('runRepeatedAnova(혼합 설계): 표 구성, 읽을 행 굵게, 개체-간 효과, 단순 효과, 쌍별 비교, 상호작용 그림, 렌더·내보내기', () => {
  const c = ref.cases[1];                                   // 3집단 × 4시점 불균형, 구형성 위반
  const { names, ds } = dsOf(c, { 1: '통제', 2: '프로그램A', 3: '프로그램B' });
  const res = A.runRepeatedAnova(ds, { vars: names, factor: ['grp'], options: { posthoc: true, plot: true } });
  assert.equal(res.title, '혼합 설계 분산분석 (집단 × 반복측정)');
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption.split(' (')[0]),
    ['처리 요약', '기술통계', 'Levene 등분산 검정', 'Mauchly 구형성 검정', '개체-내 효과 검정', '개체-내 효과 검정', '개체-간 효과 검정', '측정의 단순 효과', '집단의 단순 효과', '쌍별 비교']);
  const sum = table(res, '처리 요약').rows;
  assert.deepEqual(sum.map(r => r.n), [30, 0, 4, 3, 8, 12, 10, 30]);
  const desc = table(res, '기술통계').rows;
  assert.equal(desc.length, 3 * 4 + 4);
  assert.ok(desc[0].grp === '통제' && near(desc[0].mean, c.cells[0][0].mean, 1e-12) && desc[0].n === 8);
  assert.ok(desc[12].grp === '전체' && near(desc[12].mean, c.col_mean[0], 1e-12) && near(desc[12].sd, c.col_sd[0], 1e-10));
  table(res, 'Levene').rows.forEach((r, j) => assert.ok(near(r.W, c.levene[j].W, 1e-8) && near(r.p, c.levene[j].p, 1e-6, 1e-12) && r.df1 === 2 && r.df2 === 27));
  const mau = table(res, 'Mauchly').rows[0];
  assert.ok(near(mau.p, c.mauchly.p, 1e-7, 1e-12) && near(mau.gg, c.eps.gg, 1e-8) && near(mau.hf, c.eps.hf, 1e-8));
  const tt = res.blocks.filter(b => b.type === 'table' && b.caption.startsWith('개체-내 효과 검정'));
  assert.ok(tt[0].caption.includes('(측정)') && tt[1].caption.includes('(측정 × 집단)'));
  assert.deepEqual(tt[0].rows.map(r => r.name), ['구형성 가정', 'Greenhouse–Geisser', 'Huynh–Feldt', '하한']);
  assert.ok(c.mauchly.p < 0.05 && c.eps.gg <= 0.75);
  assert.deepEqual(tt[0].rows.map(r => !!r._bold), [false, true, false, false]);                // 구형성 위반, ε(G-G) ≤ .75 → Greenhouse–Geisser
  assert.ok(near(tt[0].rows[0].F, c.time.F, 1e-8) && near(tt[1].rows[0].F, c.inter.F, 1e-8) && near(tt[1].rows[1].p, c.inter.corr.gg.p, 1e-6, 1e-12));
  assert.ok(near(tt[1].rows[0].peta, c.inter.peta, 1e-9) && near(tt[0].rows[2].dfT, c.time.corr.hf.df1, 1e-8));
  const bt = table(res, '개체-간 효과 검정').rows;
  assert.deepEqual(bt.map(r => r.source), ['절편', '집단', '오차']);
  assert.ok(near(bt[1].F, c.between.group.F, 1e-8) && near(bt[1].peta, c.between.group.peta, 1e-9) && near(bt[0].ss, c.between.intercept.ss, 1e-9) && bt[2].df === 27);
  assert.ok(res.warnings.some(w => w.includes('집단마다 사람 수가 다릅니다')) && res.warnings.some(w => w.includes('구형성이 깨졌습니다')));
  const st = table(res, '측정의 단순 효과').rows;
  c.time_in_group.forEach((e, h) => {
    assert.ok(st[h].n === e.n && near(st[h].F, e.F, 1e-7) && near(st[h].p, e.p, 1e-6, 1e-12) && near(st[h].pgg, e.p_gg, 1e-6, 1e-12) && st[h].df1 === e.df1 && st[h].df2 === e.df2, `집단별 시간 효과 ${h}`);
  });
  const sg = table(res, '집단의 단순 효과').rows;
  c.group_at_time.forEach((e, j) => assert.ok(near(sg[j].F, e.F, 1e-7) && near(sg[j].p, e.p, 1e-6, 1e-12) && sg[j].df1 === 2 && sg[j].df2 === 27, `시점별 집단 효과 ${j}`));
  const pw = table(res, '쌍별 비교').rows;
  assert.equal(pw.length, 3 * 6);
  c.time_in_group.forEach((e, h) => e.pairs.forEach((pr, q) => {
    const row = pw[h * 6 + q];
    assert.ok(near(row.diff, pr.diff, 1e-9) && near(row.t, pr.t, 1e-8) && near(row.p, pr.p, 1e-6, 1e-12) && near(row.pAdj, pr.p_bonf, 1e-6, 1e-12) && near(row.lo, pr.lo, 1e-8) && near(row.hi, pr.hi, 1e-8), `쌍 ${h},${q}`);
  }));
  const ch = res.blocks.find(b => b.type === 'chart');
  assert.equal(ch.kind, 'line');
  assert.deepEqual(ch.xTicks, names);
  assert.equal(ch.series.length, 3);
  assert.ok(ch.series[1].label.includes('프로그램A') && near(ch.series[1].points[2][1], c.cells[1][2].mean, 1e-12));
  assert.ok(R.renderResultHtml(res).includes('혼합 설계') && R.renderResultHtml(res).includes('<svg'));
  assert.ok(Ex.resultToGrid(res).rows.length > 60);
});

test('runRepeatedAnova(혼합 설계): 측정이 2회이면 구형성 표 없이 구형성 가정 행만, 균형 자료는 pingouin 값, 옵션 끄기', () => {
  const c = ref.cases[2];                                   // 2집단 × 2시점
  const { names, ds } = dsOf(c);
  const res = A.runRepeatedAnova(ds, { vars: names, factor: ['grp'], options: { simple: true } });
  assert.ok(!table(res, 'Mauchly'));
  assert.ok(res.blocks.some(b => b.type === 'text' && b.text.includes('측정이 2회')));
  const tt = res.blocks.filter(b => b.type === 'table' && b.caption.startsWith('개체-내 효과 검정'));
  assert.deepEqual(tt.map(t => t.rows.length), [1, 1]);
  assert.ok(near(tt[1].rows[0].F, c.inter.F, 1e-8) && near(tt[0].rows[0].p, c.time.p, 1e-6, 1e-12) && tt[0].rows[0]._bold === true);
  const st = table(res, '측정의 단순 효과');
  assert.ok(!st.columns.some(x => x.key === 'pgg'));
  c.time_in_group.forEach((e, h) => assert.ok(near(st.rows[h].F, e.F, 1e-7) && near(st.rows[h].p, e.p, 1e-6, 1e-12) && st.rows[h].df1 === 1, `2회 집단별 ${h}`));
  assert.ok(!res.blocks.some(b => b.type === 'chart') && !table(res, '쌍별 비교'));
  // 옵션: 단순 효과·Levene 끄기
  const off = A.runRepeatedAnova(ds, { vars: names, factor: ['grp'], options: { simple: false, levene: false } });
  assert.ok(!table(off, '측정의 단순') && !table(off, '집단의 단순') && !table(off, 'Levene'));
  // 균형 자료(3집단 × 3시점)는 구형성 가정 행을 읽는 경우의 굵게 표시
  const b = ref.cases[3];
  const d2 = dsOf(b);
  const r2 = A.runRepeatedAnova(d2.ds, { vars: d2.names, factor: ['grp'], options: {} });
  assert.ok(b.mauchly.p >= 0.05);
  assert.deepEqual(r2.blocks.filter(x => x.type === 'table' && x.caption.startsWith('개체-내 효과 검정'))[1].rows.map(r => !!r._bold), [true, false, false, false]);
  assert.ok(!r2.warnings.some(w => w.includes('사람 수가 다릅니다')));
  assert.ok(table(r2, '개체-간 효과 검정').notes.some(n => n.includes('상호작용이 유의합니다')));
});

test('runRepeatedAnova(혼합 설계): 결측 제외, 입력 검증, 집단변수 없는 경로는 그대로', () => {
  const c = ref.cases[0];
  const { names, ds } = dsOf(c, null, rows => { rows[3][1] = ''; rows[9][3] = ''; });
  const res = A.runRepeatedAnova(ds, { vars: names, factor: ['grp'], options: {} });
  assert.ok(res.warnings.some(w => w.includes('2건을 제외') && w.includes('22건')));
  assert.deepEqual(table(res, '처리 요약').rows.map(r => r.n).slice(0, 2), [22, 2]);
  const { ds: ds0 } = dsOf(c);
  const run = (req, d) => A.runRepeatedAnova(d || ds0, Object.assign({ vars: ['t1', 't2', 't3'], factor: ['grp'], options: {} }, req));
  assert.throws(() => run({ vars: ['t1'] }), /2개 이상/);
  assert.throws(() => run({ factor: ['grp', 't1'] }), /하나만/);
  assert.throws(() => run({ factor: ['t1'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => run({ options: { conf: 100 } }), /신뢰수준/);
  const one = Data.buildDataset([['a', 'b', 'g'], [1, 2, 1], [2, 3, 1], [3, 1, 1]], [{ name: 'g', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]);
  assert.throws(() => A.runRepeatedAnova(one, { vars: ['a', 'b'], factor: ['g'], options: {} }), /집단이 1개/);
  const many = [['a', 'b', 'g']];
  for (let i = 0; i < 24; i++) many.push([i, i * 2 % 7, i % 12 + 1]);
  assert.throws(() => A.runRepeatedAnova(Data.buildDataset(many, [{ name: 'g', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]), { vars: ['a', 'b'], factor: ['g'], options: {} }), /집단이 12개/);
  const few = Data.buildDataset([['a', 'b', 'g'], [1, 2, 1], [2, 3, 2]], [{ name: 'g', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]);
  assert.throws(() => A.runRepeatedAnova(few, { vars: ['a', 'b'], factor: ['g'], options: {} }), /집단 수\(2\)보다 많아야/);
  // 집단변수를 비우면 예전처럼 3개 이상 안내
  assert.throws(() => A.runRepeatedAnova(ds0, { vars: ['t1', 't2'], options: {} }), /3개 이상/);
  assert.throws(() => A.runRepeatedAnova(ds0, { vars: ['t1', 't2'], factor: [], options: {} }), /3개 이상/);
  assert.equal(A.runRepeatedAnova(ds0, { vars: ['t1', 't2', 't3'], factor: [], options: {} }).title, '일원 반복측정 분산분석');
});
