const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const A2 = require('../src/Stats_Anova2.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');
const Ch = require('../src/Charts.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'anova2.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-12) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));

function dsOf(c, valueLabels) {
  const rows = [['y', 'A', 'B']].concat(c.y.map((v, i) => [v, c.A[i], c.B[i]]));
  const defs = [
    { name: 'y', label: '점수', scale: 'scale', valueLabels: {}, missing: [] },
    { name: 'A', label: '요인A', scale: 'nominal', valueLabels: (valueLabels && valueLabels.A) || {}, missing: [] },
    { name: 'B', label: '요인B', scale: 'nominal', valueLabels: (valueLabels && valueLabels.B) || {}, missing: [] }
  ];
  return Data.buildDataset(rows, defs);
}
const idx0 = arr => arr.map(v => v - 1);

test('anova2Analyze: statsmodels·pingouin·numpy가 서로 일치한 기준값과 같다 — 유형 III·II 제곱합, F, p, 부분 η², 절편, 오차, R² (5가지 사례)', () => {
  assert.equal(ref.cases.length, 5);
  ref.cases.forEach(c => {
    ['III', 'II'].forEach(type => {
      const r = A2.anova2Analyze(c.y, idx0(c.A), idx0(c.B), c.a, c.b, type);
      const w = `${c.name} 유형 ${type}`;
      const ex = c[type];
      r.terms.forEach(tm => {
        const e = ex.terms[tm.key];
        assert.ok(near(tm.ss, e.ss, 1e-8, 1e-10) && tm.df === e.df && near(tm.ms, e.ms, 1e-8, 1e-10), `${w} ${tm.key} 제곱합`);
        assert.ok(near(tm.F, e.F, 1e-7) && near(tm.p, e.p, 1e-6, 1e-12) && near(tm.peta, e.peta, 1e-8, 1e-12), `${w} ${tm.key} F, p, η² ${tm.F} vs ${e.F}`);
      });
      assert.ok(near(r.sse, c.sse, 1e-9) && r.dfe === c.dfe && near(r.mse, c.mse, 1e-9) && near(r.sst, c.sst, 1e-10) && near(r.sumSq, c.sumsq, 1e-10), w + ' 오차·합계');
      assert.ok(near(r.r2, c.r2, 1e-9) && near(r.adjR2, c.adj, 1e-9), w + ' R²');
      if (type === 'III') {
        assert.ok(near(r.intercept.ss, ex.intercept.ss, 1e-8) && near(r.intercept.F, ex.intercept.F, 1e-7) && near(r.intercept.p, ex.intercept.p, 1e-6, 1e-12), w + ' 절편');
        assert.ok(near(r.model.ss, ex.model.ss, 1e-8) && near(r.model.F, ex.model.F, 1e-7) && near(r.model.peta, ex.model.peta, 1e-8) && r.model.df === ex.model.df, w + ' 수정된 모형');
      } else {
        assert.equal(r.intercept, undefined);
      }
    });
    // 셀 통계와 관측 주변평균
    const r = A2.anova2Analyze(c.y, idx0(c.A), idx0(c.B), c.a, c.b, 'III');
    r.cells.forEach((row, i) => row.forEach((cell, j) => {
      const e = c.cells[i][j];
      assert.ok(cell.n === e.n && near(cell.mean, e.mean, 1e-12) && near(cell.sd, e.sd, 1e-10), `${c.name} 셀 ${i},${j}`);
    }));
    assert.deepEqual(r.rowN, c.row_n); assert.deepEqual(r.colN, c.col_n);
    r.rowMean.forEach((v, k) => assert.ok(near(v, c.row_mean[k], 1e-12)));
    r.colMean.forEach((v, k) => assert.ok(near(v, c.col_mean[k], 1e-12)));
  });
  // 균형 자료에서는 두 유형이 같고, 불균형 자료에서는 다르다
  const bal = ref.cases[0], unb = ref.cases[1];
  const s3 = A2.anova2Analyze(bal.y, idx0(bal.A), idx0(bal.B), bal.a, bal.b, 'III').terms, s2 = A2.anova2Analyze(bal.y, idx0(bal.A), idx0(bal.B), bal.a, bal.b, 'II').terms;
  s3.forEach((tm, k) => assert.ok(near(tm.ss, s2[k].ss, 1e-9, 1e-9)));
  const u3 = A2.anova2Analyze(unb.y, idx0(unb.A), idx0(unb.B), unb.a, unb.b, 'III').terms, u2 = A2.anova2Analyze(unb.y, idx0(unb.A), idx0(unb.B), unb.a, unb.b, 'II').terms;
  assert.ok(Math.abs(u3[0].ss - u2[0].ss) > 1e-3);
  assert.ok(near(u3[2].ss, u2[2].ss, 1e-9, 1e-9));              // 상호작용은 유형 II와 III 가 같다
});

test('anova2Emm·anova2PostHoc·anova2Simple: statsmodels(t_test, f_test)·scipy 기준값과 같다', () => {
  ref.cases.forEach(c => {
    const r = A2.anova2Analyze(c.y, idx0(c.A), idx0(c.B), c.a, c.b, 'III');
    const w = c.name;
    const em = A2.anova2Emm(r, c.conf);
    const cmp = (g, e, what) => assert.ok(near(g.mean, e.mean, 1e-9) && near(g.se, e.se, 1e-8) && near(g.lo, e.lo, 1e-8) && near(g.hi, e.hi, 1e-8), `${w} EMM ${what}`);
    em.rows.forEach((g, i) => cmp(g, c.emm.rows[i], '행 ' + i));
    em.cols.forEach((g, j) => cmp(g, c.emm.cols[j], '열 ' + j));
    em.cells.forEach((row, i) => row.forEach((g, j) => cmp(g, c.emm.cells[i][j], `셀 ${i},${j}`)));
    [['A', r.rowMean, r.rowN, c.a], ['B', r.colMean, r.colN, c.b]].forEach(([key, means, ns, k]) => {
      const exp = c.posthoc[key];
      if (!exp) { assert.ok(k < 3); return; }
      const tk = A2.anova2PostHoc('tukey', means, ns, r.mse, r.dfe, c.conf), bf = A2.anova2PostHoc('bonferroni', means, ns, r.mse, r.dfe, c.conf);
      assert.equal(tk.length, exp.length);
      exp.forEach((e, q) => {
        assert.ok(tk[q].i === e.i && tk[q].j === e.j && near(tk[q].diff, e.diff, 1e-10) && near(tk[q].se, e.se, 1e-9), `${w} ${key} 쌍 ${q}`);
        assert.ok(near(tk[q].p, e.p_tukey, 1e-5, 1e-9) && near(tk[q].lo, e.lo_t, 1e-5, 1e-6) && near(tk[q].hi, e.hi_t, 1e-5, 1e-6), `${w} ${key} Tukey ${q}: ${tk[q].p} vs ${e.p_tukey}`);
        assert.ok(near(bf[q].p, e.p_bonf, 1e-7, 1e-12) && near(bf[q].lo, e.lo_b, 1e-7, 1e-8) && near(bf[q].hi, e.hi_b, 1e-7, 1e-8), `${w} ${key} Bonferroni ${q}`);
      });
    });
    const sm = A2.anova2Simple(r);
    ['aAtB', 'bAtA'].forEach(key => sm[key].forEach((g, k) => {
      const e = c.simple[key][k];
      assert.ok(near(g.ss, e.ss, 1e-9, 1e-10) && g.df === e.df && near(g.F, e.F, 1e-7) && near(g.p, e.p, 1e-6, 1e-12), `${w} 단순 주효과 ${key} ${k}`);
    }));
  });
});

test('runAnova2: 표 구성, 기술통계(주변 행), Levene, 개체-간 효과 검정, 단순 주효과, 상호작용 그림, 렌더·내보내기', () => {
  const c = ref.cases[1];                                    // 불균형 3×3
  const ds = dsOf(c, { A: { 1: '낮음', 2: '중간', 3: '높음' }, B: { 1: '가', 2: '나', 3: '다' } });
  const res = A2.runAnova2(ds, { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: { emmeans: true } });
  assert.equal(res.title, '이원 분산분석');
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption.split(' (')[0]), ['기술통계', 'Levene 등분산 검정', '개체-간 효과 검정', '추정 주변평균', '사후검정: 요인A', '사후검정: 요인B', '단순 주효과']);
  const desc = table(res, '기술통계').rows;
  assert.equal(desc.length, 9 + 3 + 3 + 1);                   // 셀 9 + 행별 전체 3 + 열별 전체 3 + 전체
  assert.deepEqual(desc.slice(0, 4).map(r => [r.a, r.b]), [['낮음', '가'], ['낮음', '나'], ['낮음', '다'], ['낮음', '전체']]);
  assert.ok(near(desc[0].mean, c.cells[0][0].mean, 1e-12) && desc[0].n === c.cells[0][0].n && near(desc[0].sd, c.cells[0][0].sd, 1e-10));
  assert.ok(near(desc[3].mean, c.row_mean[0], 1e-12) && desc[3].n === c.row_n[0]);
  assert.deepEqual([desc[desc.length - 1].a, desc[desc.length - 1].b, desc[desc.length - 1].n], ['전체', '전체', c.n]);
  const lv = table(res, 'Levene').rows[0];
  assert.ok(near(lv.W, c.levene.W, 1e-8) && lv.df1 === c.levene.df1 && lv.df2 === c.levene.df2 && near(lv.p, c.levene.p, 1e-6, 1e-12));
  const tt = table(res, '개체-간 효과 검정');
  assert.deepEqual(tt.rows.map(r => r.source), ['수정된 모형', '절편', '요인A', '요인B', '요인A × 요인B', '오차', '합계', '수정 합계']);
  assert.ok(near(tt.rows[2].F, c.III.terms.A.F, 1e-7) && near(tt.rows[2].peta, c.III.terms.A.peta, 1e-8) && near(tt.rows[4].p, c.III.terms.AB.p, 1e-6, 1e-12));
  assert.ok(tt.notes[0].includes('유형 III') && tt.notes.some(n => n.includes('R 제곱 =')));
  assert.ok(res.warnings.some(w => w.includes('불균형')));
  const sim = table(res, '단순 주효과').rows;
  assert.equal(sim.length, 6);
  assert.ok(sim[0].eff.includes('요인A의 효과 (요인B = 가)') && sim[3].eff.includes('요인B의 효과 (요인A = 낮음)'));
  assert.ok(near(sim[0].F, c.simple.aAtB[0].F, 1e-7) && near(sim[3].p, c.simple.bAtA[0].p, 1e-6, 1e-12));
  const emm = table(res, '추정 주변평균').rows;
  assert.equal(emm.length, 3 + 3 + 9);
  assert.ok(near(emm[0].mean, c.emm.rows[0].mean, 1e-9) && near(emm[3].se, c.emm.cols[0].se, 1e-8) && emm[6].lev === '낮음, 가');
  const ph = table(res, '사후검정: 요인A');
  assert.equal(ph.rows.length, 3);
  assert.equal(ph.rows[0].i, '낮음');
  assert.ok(near(ph.rows[0].p, c.posthoc.A[0].p_tukey, 1e-5, 1e-9) && ph.caption.includes('Tukey HSD'));
  const ch = res.blocks.find(b => b.type === 'chart');
  assert.equal(ch.kind, 'line');
  assert.deepEqual(ch.xTicks, ['낮음', '중간', '높음']);
  assert.equal(ch.series.length, 3);
  assert.ok(ch.series[0].label.includes('요인B = 가') && near(ch.series[1].points[2][1], c.cells[2][1].mean, 1e-12));
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('개체-간 효과 검정') && html.includes('<svg') && html.includes('중간') && html.includes('높음'));
  assert.ok(Ex.resultToGrid(res).rows.length > 40);
});

test('runAnova2: 옵션 (유형 II, Bonferroni, 사후검정·단순 주효과·그림·Levene 끄기, 2수준 요인은 사후검정 없음)', () => {
  const c = ref.cases[3];                                    // 불균형 2×2
  const ds = dsOf(c);
  const r2 = A2.runAnova2(ds, { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: { ss: 'II', posthoc: 'bonferroni', levene: false, simple: false, plot: false } });
  assert.deepEqual(r2.blocks.filter(b => b.type === 'table').map(b => b.caption.split(' (')[0]), ['기술통계', '개체-간 효과 검정']);
  assert.ok(!r2.blocks.some(b => b.type === 'chart'));
  const tt = table(r2, '개체-간 효과 검정');
  assert.deepEqual(tt.rows.map(r => r.source), ['수정된 모형', '요인A', '요인B', '요인A × 요인B', '오차', '합계', '수정 합계']);
  assert.ok(near(tt.rows[1].ss, c.II.terms.A.ss, 1e-8, 1e-10) && near(tt.rows[2].ss, c.II.terms.B.ss, 1e-8, 1e-10));
  assert.ok(tt.notes[0].includes('유형 II'));
  const c5 = ref.cases[4];                                   // 4×3: 두 요인 모두 사후검정
  const d5 = dsOf(c5);
  const r5 = A2.runAnova2(d5, { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: { posthoc: 'bonferroni' } });
  const ph = table(r5, '사후검정: 요인A');
  assert.equal(ph.rows.length, 6);
  assert.ok(ph.caption.includes('Bonferroni') && near(ph.rows[0].p, c5.posthoc.A[0].p_bonf, 1e-7, 1e-12));
  assert.equal(table(r5, '사후검정: 요인B').rows.length, 3);
  const none = A2.runAnova2(d5, { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: { posthoc: 'none' } });
  assert.ok(!table(none, '사후검정'));
  // 2×2 에서는 기본 옵션이어도 사후검정 표가 없다
  assert.ok(!table(A2.runAnova2(ds, { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} }), '사후검정'));
});

test('runAnova2: 경고 (불균형, 작은 셀, Levene 유의, 상호작용 안내)와 결측 제외', () => {
  const small = ref.cases[3];
  const res = A2.runAnova2(dsOf(small), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} });
  assert.ok(res.warnings.some(w => w.includes('불균형')) && res.warnings.some(w => w.includes('5건 미만인 셀이 2개')));
  const bal = A2.runAnova2(dsOf(ref.cases[0]), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} });
  assert.ok(!bal.warnings.some(w => w.includes('불균형')));
  assert.ok(table(bal, '개체-간').notes.some(n => n.includes('상호작용이 유의합니다')));
  const noInt = A2.runAnova2(dsOf(ref.cases[2]), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} });
  assert.ok(table(noInt, '개체-간').notes.some(n => n.includes('상호작용이 유의하지 않습니다')) && table(noInt, '단순 주효과').notes[0].includes('참고만'));
  // 분산이 매우 다른 셀 → Levene 경고
  const rows = [['y', 'A', 'B']];
  for (let i = 0; i < 12; i++) for (const [a, b, sd] of [[1, 1, 0.1], [1, 2, 0.1], [2, 1, 6], [2, 2, 8]]) rows.push([10 + (i % 6 - 2.5) * sd + a + b, a, b]);
  const lv = A2.runAnova2(Data.buildDataset(rows, []), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} });
  assert.ok(lv.warnings.some(w => w.includes('Levene')));
  // 결측 사례는 제외
  const c = ref.cases[0];
  const miss = [['y', 'A', 'B']].concat(c.y.map((v, i) => [v, c.A[i], c.B[i]]));
  miss[3][0] = ''; miss[7][1] = '';
  const rm = A2.runAnova2(Data.buildDataset(miss, [{ name: 'A', scale: 'nominal', label: '', valueLabels: {}, missing: [] }, { name: 'B', scale: 'nominal', label: '', valueLabels: {}, missing: [] }]), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} });
  assert.ok(rm.warnings.some(w => w.includes('2건을 제외') && w.includes('58건')));
  assert.equal(table(rm, '개체-간').rows.find(r => r.source === '오차').df, 58 - 6);
});

test('runAnova2: 입력 검증과 계산 불가 경우', () => {
  const ds = dsOf(ref.cases[0]);
  const run = (req, d) => A2.runAnova2(d || ds, Object.assign({ dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} }, req));
  assert.throws(() => run({ dependent: ['y', 'A'] }), /하나만/);
  assert.throws(() => run({ factor1: ['A', 'B'] }), /변수 하나씩/);
  assert.throws(() => run({ factor2: ['A'] }), /"A"를 요인 1와\(과\) 요인 2에 동시에/);
  assert.throws(() => run({ factor1: ['y'] }), /"y"를 종속변수와\(과\) 요인 1에/);
  assert.throws(() => run({ options: { conf: 100 } }), /신뢰수준/);
  // 비어 있는 셀
  const rows = [['y', 'A', 'B']];
  [[1, 1], [1, 1], [1, 2], [1, 2], [2, 1], [2, 1], [2, 1]].forEach(([a, b], i) => rows.push([i + 1, a, b]));
  assert.throws(() => A2.runAnova2(Data.buildDataset(rows, []), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} }), /조합.*사례가 없습니다/);
  // 셀 수 이상이 아닌 사례 수
  const few = [['y', 'A', 'B'], [1, 1, 1], [2, 1, 2], [3, 2, 1], [4, 2, 2]];
  assert.throws(() => A2.runAnova2(Data.buildDataset(few, []), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} }), /셀 수\(4\)보다 많아야/);
  // 수준이 하나뿐인 요인, 수준이 너무 많은 요인
  const one = [['y', 'A', 'B']];
  for (let i = 0; i < 8; i++) one.push([i, 1, i % 2 + 1]);
  assert.throws(() => A2.runAnova2(Data.buildDataset(one, []), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} }), /수준이 1개/);
  const many = [['y', 'A', 'B']];
  for (let i = 0; i < 24; i++) many.push([i, i % 12 + 1, i % 2 + 1]);
  assert.throws(() => A2.runAnova2(Data.buildDataset(many, []), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} }), /수준이 12개/);
  // 숫자가 아닌 종속변수
  const txt = Data.buildDataset([['y', 'A', 'B'], ['x', 1, 1], ['y', 1, 2], ['z', 2, 1], ['w', 2, 2]], [{ name: 'y', scale: 'ordinal', label: '', valueLabels: {}, missing: [] }]);
  assert.throws(() => A2.runAnova2(txt, { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} }), /숫자가 아닌 값/);
  // 모든 셀 값이 같아 오차 0
  const flat = [['y', 'A', 'B']];
  for (let i = 0; i < 12; i++) flat.push([[5, 6, 7, 8][i % 4], i % 2 + 1, Math.floor(i / 2) % 2 + 1]);
  const rf = A2.runAnova2(Data.buildDataset(flat, []), { dependent: ['y'], factor1: ['A'], factor2: ['B'], options: {} });
  assert.ok(rf.warnings.some(w => w.includes('오차가 0')));
  assert.equal(table(rf, '개체-간').rows.find(r => r.source === 'A').F, null);
});

test('lineChartBlock: xTicks 는 범주 이름을 x 눈금으로 그린다 (없으면 기존과 같다)', () => {
  const plain = Ch.lineChartBlock('c', 'x', 'y', [{ label: 's', points: [[1, 2], [2, 3]] }]);
  assert.equal(plain.xTicks, undefined);
  const cat = Ch.lineChartBlock('c', 'x', 'y', [{ label: 's', points: [[1, 2], [2, 3], [3, 1]] }], ['가', '나', '다']);
  assert.deepEqual(cat.xTicks, ['가', '나', '다']);
  const svg = R.renderResultHtml({ analysis: 'x', title: 't', warnings: [], blocks: [cat] });
  assert.ok(svg.includes('>가<') && svg.includes('>나<') && svg.includes('>다<'));
  const svgPlain = R.renderResultHtml({ analysis: 'x', title: 't', warnings: [], blocks: [plain] });
  assert.ok(!svgPlain.includes('>가<'));
});
