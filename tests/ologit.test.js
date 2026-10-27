const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const O = require('../src/Stats_Ologit.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'ologit.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-12) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const def = (name, scale, vl) => ({ name, label: '', scale: scale || 'scale', valueLabels: vl || {}, missing: [] });

function dsOf(c, vl, extra) {
  const header = ['y'].concat(c.names);
  const rows = [header].concat(c.y.map((v, i) => [v].concat(c.cols.map(col => col[i]))));
  if (extra) extra(rows);
  return Data.buildDataset(rows, [def('y', 'ordinal', vl)].concat(c.names.map(nm => def(nm))));
}
const predsOf = c => c.names.map((nm, i) => ({ name: nm, values: c.cols[i] }));
const yOf = c => c.y.map(v => c.cats.indexOf(v));

test('ordinalFit: statsmodels OrderedModel 과 numpy 로그가능도·수치 헤세가 일치한 기준값과 같다 — 계수, 문턱값, 표준오차, 공분산, 로그가능도, 예측 확률 (5가지 자료)', () => {
  assert.equal(ref.cases.length, 5);
  ref.cases.forEach(c => {
    const f = O.ordinalFit(predsOf(c), yOf(c), c.J);
    const m = c.name;
    assert.ok(f.converged && !f.separation, m + ' 수렴');
    f.beta.forEach((v, i) => assert.ok(near(v, c.beta[i], 1e-6, 1e-8), `${m} 계수 ${i}: ${v} vs ${c.beta[i]}`));
    f.theta.forEach((v, j) => assert.ok(near(v, c.theta[j], 1e-6, 1e-8), `${m} 문턱값 ${j}`));
    f.se.forEach((v, i) => assert.ok(near(v, c.se[i], 1e-5, 1e-8), `${m} 표준오차 ${i}: ${v} vs ${c.se[i]}`));
    f.cov.forEach((row, i) => row.forEach((v, j) => assert.ok(Math.abs(v - c.cov[i][j]) <= 2e-4 * Math.sqrt(c.cov[i][i] * c.cov[j][j]) + 1e-9, `${m} 공분산 ${i},${j}: ${v} vs ${c.cov[i][j]}`)));
    assert.ok(near(f.ll, c.ll, 1e-10) && near(f.ll0, c.ll0, 1e-12), m + ' 로그가능도');
    assert.deepEqual(f.counts, c.counts);
    c.probs.forEach((row, i) => row.forEach((v, j) => assert.ok(near(f.probs[i][j], v, 1e-6, 1e-9), `${m} 확률 ${i},${j}`)));
    f.probs.forEach(row => assert.ok(near(row.reduce((s, v) => s + v, 0), 1, 1e-12)));
  });
});

test('brantTest: 문턱값별 이분형 로지스틱 계수는 statsmodels 와 같고 Wald 통계량은 numpy 식과 같다 (분리가 있으면 계산하지 않음)', () => {
  ref.cases.forEach((c, q) => {
    const b = O.brantTest(predsOf(c), yOf(c), c.J);
    if (q === 1) {                                           // 마지막 문턱값에서 완전 분리 → 계산 불가
      assert.equal(b.ok, false);
      assert.ok(b.reason.includes('분리'));
      return;
    }
    assert.ok(b.ok, c.name);
    b.splitB.forEach((row, m) => row.forEach((v, i) => assert.ok(near(v, c.brant.splitB[m][i], 1e-6, 1e-8), `${c.name} 문턱값 ${m} 계수 ${i}`)));
    const e = c.brant;
    assert.ok(near(b.omnibus.chi2, e.omnibus.chi2, 1e-6) && b.omnibus.df === e.omnibus.df && near(b.omnibus.p, e.omnibus.p, 1e-5, 1e-12), `${c.name} 전체 ${b.omnibus.chi2} vs ${e.omnibus.chi2}`);
    b.perVar.forEach((v, i) => assert.ok(near(v.chi2, e.perVar[i].chi2, 1e-6) && v.df === e.perVar[i].df && near(v.p, e.perVar[i].p, 1e-5, 1e-12), `${c.name} 변수 ${i}`));
  });
  // 위반 자료는 Brant 가 유의하고 성립 자료는 유의하지 않다
  assert.ok(ref.cases[2].brant.omnibus.p < 0.05 && ref.cases[0].brant.omnibus.p > 0.05);
});

test('runOrdinalLogit: 표 구성, 문턱값·위치 모수 값, Exp(B)와 구간, 모형 적합 정보, 유사 R², Brant 표, 렌더·내보내기', () => {
  const c = ref.cases[0];
  const ds = dsOf(c, { 1: '낮음', 2: '다소 낮음', 3: '다소 높음', 4: '높음' });
  const res = O.runOrdinalLogit(ds, { dependent: ['y'], covariates: c.names, options: {} });
  assert.equal(res.title, '서열 로지스틱 회귀 (비례오즈)');
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption.split(' (')[0]), ['케이스 처리 요약', '모형 적합 정보', '유사 R 제곱', '모수 추정값', '비례오즈 가정 점검']);
  const sum = table(res, '케이스 처리 요약').rows;
  assert.deepEqual(sum.slice(0, 4).map(r => r.n), c.counts);
  assert.ok(sum[0].label.includes('y = 낮음'));
  const fi = table(res, '모형 적합 정보').rows;
  assert.ok(near(fi[0].m2ll, -2 * c.ll0, 1e-12) && near(fi[1].m2ll, -2 * c.ll, 1e-10) && near(fi[1].chi2, c.chi2, 1e-9) && fi[1].df === c.df && near(fi[1].p, c.p, 1e-6, 1e-12));
  const pr = table(res, '유사 R 제곱').rows[0];
  assert.ok(near(pr.cs, c.cs, 1e-9) && near(pr.nk, c.nk, 1e-9) && near(pr.mf, c.mf, 1e-9));
  const pt = table(res, '모수 추정값').rows;
  assert.equal(pt.length, 3 + 2);
  assert.ok(pt[0].name.includes('문턱값 [y ≤ 낮음]') && pt[2].name.includes('문턱값 [y ≤ 다소 높음]'));
  pt.slice(0, 3).forEach((r, j) => assert.ok(near(r.B, c.theta[j], 1e-6, 1e-8) && near(r.se, c.se[c.k + j], 1e-5) && r.or === null));
  pt.slice(3).forEach((r, i) => {
    assert.equal(r.name, c.names[i]);
    const b = c.beta[i], se = c.se[i], z = 1.959963984540054;
    assert.ok(near(r.B, b, 1e-6, 1e-8) && near(r.se, se, 1e-5) && near(r.wald, (b / se) ** 2, 1e-5) && near(r.or, Math.exp(b), 1e-6) && near(r.orLo, Math.exp(b - z * se), 1e-5) && near(r.orHi, Math.exp(b + z * se), 1e-5), `위치 모수 ${i}`);
    assert.ok(near(r.lo, b - z * se, 1e-5, 1e-7) && near(r.hi, b + z * se, 1e-5, 1e-7));
  });
  const bt = table(res, '비례오즈 가정 점검');
  assert.equal(bt.rows.length, 3);
  assert.deepEqual(bt.rows[2].name, '전체 (Brant)');
  assert.ok(near(bt.rows[2].chi2, c.brant.omnibus.chi2, 1e-6) && bt.rows[2].df === 2 * 2 && near(bt.rows[0].b1, c.brant.splitB[1][0], 1e-6, 1e-8));
  assert.ok(!res.warnings.some(w => w.includes('비례오즈 가정이 의심')));
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('모수 추정값') && html.includes('문턱값'));
  assert.ok(Ex.resultToGrid(res).rows.length > 15);
});

test('runOrdinalLogit: 범주형 독립변수(더미), 기준 범주 마지막 옵션, 90% 신뢰수준', () => {
  const c = ref.cases[1];
  // 요인 g 와 숫자 a, b 로 다시 만든다 (기준값은 g = 2, g = 3 더미 = 첫 범주 기준)
  const header = ['y', 'a', 'b', 'g'];
  const rows = [header].concat(c.y.map((v, i) => [v, c.cols[0][i], c.cols[1][i], c.g[i]]));
  const ds = Data.buildDataset(rows, [def('y', 'ordinal'), def('a'), def('b'), def('g', 'nominal', { 1: '가', 2: '나', 3: '다' })]);
  const res = O.runOrdinalLogit(ds, { dependent: ['y'], covariates: ['a', 'b'], factors: ['g'], options: { parallel: false } });
  const pt = table(res, '모수 추정값');
  assert.deepEqual(pt.rows.slice(-4).map(r => r.name), ['a', 'b', 'g = 나', 'g = 다']);
  pt.rows.slice(-4).forEach((r, i) => assert.ok(near(r.B, c.beta[i], 1e-6, 1e-8), `더미 ${i}`));
  assert.ok(pt.notes.some(n => n.includes('"g"의 기준 범주: 가')));
  assert.ok(!table(res, '비례오즈'));
  const last = O.runOrdinalLogit(ds, { dependent: ['y'], covariates: ['a', 'b'], factors: ['g'], options: { parallel: false, lastRef: true } });
  assert.deepEqual(table(last, '모수 추정값').rows.slice(-2).map(r => r.name), ['g = 가', 'g = 나']);
  assert.ok(table(last, '모수 추정값').notes.some(n => n.includes('기준 범주: 다')));
  // Brant 를 켜면 분리 때문에 계산하지 못했다는 경고 (이 자료는 마지막 문턱값에서 분리)
  const br = O.runOrdinalLogit(ds, { dependent: ['y'], covariates: ['a', 'b'], factors: ['g'], options: {} });
  assert.ok(br.warnings.some(w => w.includes('비례오즈 가정 점검(Brant)을 계산하지 못했습니다')));
  // 90%
  const c5 = ref.cases[4];
  const r5 = O.runOrdinalLogit(dsOf(c5), { dependent: ['y'], covariates: c5.names, options: { conf: '90' } });
  const row = table(r5, '모수 추정값').rows[c5.J - 1];
  const z90 = 1.6448536269514722;
  assert.ok(near(row.lo, c5.beta[0] - z90 * c5.se[0], 1e-5, 1e-7) && table(r5, '모수 추정값').columns.some(x => x.label === 'B의 90% 신뢰구간 하한'));
});

test('runOrdinalLogit: 비례오즈 위반 경고, 분류표, 예측 범주·확률 저장, 결측 제외', () => {
  const v = ref.cases[2];                                   // 위반
  const res = O.runOrdinalLogit(dsOf(v), { dependent: ['y'], covariates: v.names, options: { classify: true } });
  assert.ok(res.warnings.some(w => w.includes('비례오즈 가정이 의심됩니다') && w.includes('x')));
  const ct = table(res, '분류표');
  assert.deepEqual(ct.rows.slice(0, 3).map(r => [r.c0, r.c1, r.c2]), v.table);
  assert.ok(near(ct.rows[3].pct, v.acc * 100, 1e-9) && ct.rows[3].sum === v.n);
  assert.ok(!table(O.runOrdinalLogit(dsOf(v), { dependent: ['y'], covariates: v.names, options: {} }), '분류표'));
  // 저장: 예측 범주와 범주별 확률 (제외된 사례는 null)
  const c = ref.cases[0];
  const ds = dsOf(c, { 1: '낮음', 2: '다소 낮음', 3: '다소 높음', 4: '높음' }, rows => { rows[3][1] = ''; });
  const sv = O.runOrdinalLogit(ds, { dependent: ['y'], covariates: c.names, options: { save: true, parallel: false } });
  assert.deepEqual(sv.saveColumns.map(x => x.name), ['예측범주', '확률_1', '확률_2', '확률_3', '확률_4']);
  assert.equal(sv.saveColumns[0].scale, '서열');
  assert.ok(sv.saveColumns[0].valueLabels.includes('1=낮음') && sv.saveColumns[0].valueLabels.includes('4=높음'));
  assert.equal(sv.saveColumns[0].values.length, c.n);
  assert.equal(sv.saveColumns[0].values[2], null);          // rows[3] = 사례 인덱스 2 (결측 → 제외)
  assert.ok(sv.saveColumns[1].fixed === 4 && sv.saveColumns[1].values[2] === null);
  const okRows = sv.saveColumns[0].values.map((_, i) => i).filter(i => sv.saveColumns[0].values[i] !== null);
  assert.equal(okRows.length, c.n - 1);
  okRows.forEach(i => {
    const ps = [1, 2, 3, 4].map(q => sv.saveColumns[q].values[i]);
    assert.ok(near(ps.reduce((s, x) => s + x, 0), 1, 1e-12));
    assert.equal(sv.saveColumns[0].values[i], ps.indexOf(Math.max(...ps)) + 1);
  });
  assert.ok(sv.warnings.some(w => w.includes('1건을 제외')));
  assert.equal(O.runOrdinalLogit(ds, { dependent: ['y'], covariates: c.names, options: {} }).saveColumns, undefined);
});

test('runOrdinalLogit: 입력 검증, 경고(작은 범주, 사례 수), 분리 경고', () => {
  const c = ref.cases[0];
  const ds = dsOf(c);
  const run = (req, d) => O.runOrdinalLogit(d || ds, Object.assign({ dependent: ['y'], covariates: ['x1'], options: {} }, req));
  assert.throws(() => run({ dependent: ['y', 'x1'] }), /하나만/);
  assert.throws(() => run({ covariates: [] }), /독립변수/);
  assert.throws(() => run({ covariates: ['y'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => run({ options: { conf: 100 } }), /신뢰수준/);
  const two = Data.buildDataset([['y', 'x'], [1, 1], [2, 2], [1, 3], [2, 4], [1, 2], [2, 1]], [def('y', 'ordinal')]);
  assert.throws(() => O.runOrdinalLogit(two, { dependent: ['y'], covariates: ['x'], options: {} }), /2가지.*이분형 로지스틱/);
  const many = [['y', 'x']];
  for (let i = 0; i < 44; i++) many.push([i % 11 + 1, Math.sin(i) + i * 0.1]);
  assert.throws(() => O.runOrdinalLogit(Data.buildDataset(many, [def('y', 'ordinal')]), { dependent: ['y'], covariates: ['x'], options: {} }), /11가지.*3~10가지/);
  const flat = [['y', 'x']];
  for (let i = 0; i < 30; i++) flat.push([i % 3 + 1, 5]);
  assert.throws(() => O.runOrdinalLogit(Data.buildDataset(flat, [def('y', 'ordinal')]), { dependent: ['y'], covariates: ['x'], options: {} }), /값이 모두 같아/);
  const few = [['y', 'a', 'b', 'c'], [1, 1, 5, 2], [2, 2, 3, 7], [3, 3, 1, 4], [1, 4, 2, 9], [2, 5, 8, 1]];
  assert.throws(() => O.runOrdinalLogit(Data.buildDataset(few, [def('y', 'ordinal')]), { dependent: ['y'], covariates: ['a', 'b', 'c'], options: {} }), /모수의 수/);
  // 작은 범주(5건 미만)와 사례 수 경고
  const smallRows = [['y', 'x', 'z']];
  for (let i = 0; i < 40; i++) smallRows.push([i < 18 ? 1 : i < 36 ? 2 : i < 38 ? 3 : 4, Math.cos(i) + (i < 18 ? 0 : 1.2) + (i >= 36 ? 2 : 0), Math.sin(i * 1.7)]);
  const sm = O.runOrdinalLogit(Data.buildDataset(smallRows, [def('y', 'ordinal')]), { dependent: ['y'], covariates: ['x', 'z'], options: { parallel: false } });
  assert.ok(sm.warnings.some(w => w.includes('5건 미만인 범주')) && sm.warnings.some(w => w.includes('10배에 못 미칩니다')));
  // 완전 분리: x 가 y 를 완벽히 가른다
  const sep = [['y', 'x']];
  for (let i = 0; i < 60; i++) sep.push([i < 20 ? 1 : i < 40 ? 2 : 3, i + (i < 20 ? 0 : i < 40 ? 100 : 200)]);
  const sr = O.runOrdinalLogit(Data.buildDataset(sep, [def('y', 'ordinal')]), { dependent: ['y'], covariates: ['x'], options: { parallel: false } });
  assert.ok(sr.warnings.some(w => w.includes('분리') || w.includes('수렴')));
  const txt = Data.buildDataset([['y', 'x'], [1, 'a'], [2, 'b'], [3, 'c'], [1, 'd'], [2, 'e'], [3, 'f']], [def('y', 'ordinal'), def('x', 'ordinal')]);
  assert.throws(() => O.runOrdinalLogit(txt, { dependent: ['y'], covariates: ['x'], options: {} }), /숫자가 아닌 값/);
});
