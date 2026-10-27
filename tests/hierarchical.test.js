const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const H = require('../src/Stats_Hierarchical.js');
const Reg = require('../src/Stats_Regression.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'hierarchical.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-10) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));

/** 기준값 사례에서 데이터셋과 요청을 만든다 (범주형은 명목 척도로 정의) */
function build(c, options, scaleOf) {
  const vars = [];
  c.blocks.forEach(b => b.forEach(v => vars.push(v)));
  const header = ['y'].concat(vars.map(v => v.name));
  const rows = [header];
  for (let i = 0; i < c.n; i++) rows.push([c.y[i]].concat(vars.map(v => v.values[i])));
  const defs = header.map(nm => {
    const v = vars.find(q => q.name === nm);
    const scale = nm === 'y' ? 'scale' : (scaleOf && scaleOf[nm]) || (v.kind === 'factor' ? 'nominal' : 'scale');
    return { name: nm, label: '', scale, valueLabels: {}, missing: [] };
  });
  const request = { dependent: ['y'], options: Object.assign({ conf: c.conf }, options) };
  c.blocks.forEach((b, bi) => { request['block' + (bi + 1)] = b.map(v => v.name); });
  return { dataset: Data.buildDataset(rows, defs), request };
}

test('runHierarchical: statsmodels·numpy 기준값과 같다 — 단계별 요약, 변화 F, 분산분석, 계수 (4가지 사례)', () => {
  assert.equal(ref.cases.length, 4);
  ref.cases.forEach(c => {
    const { dataset, request } = build(c, { ci: true, dw: true });
    const res = H.runHierarchical(dataset, request);
    const w = c.name;
    const sum = table(res, '모형 요약').rows, anova = table(res, '분산분석').rows, coef = table(res, '계수').rows;
    assert.equal(sum.length, c.steps.length, w);
    let ci = 0;
    c.steps.forEach((s, k) => {
      const t = `${w} ${k + 1}단계`;
      assert.equal(sum[k].step, (k + 1) + '단계');
      assert.ok(near(sum[k].r, s.r, 1e-9) && near(sum[k].r2, s.r2, 1e-9) && near(sum[k].adj, s.adj, 1e-9) && near(sum[k].s, s.s, 1e-9), t + ' 요약');
      assert.ok(near(sum[k].dR2, s.dR2, 1e-8, 1e-12) && near(sum[k].dF, s.dF, 1e-8) && sum[k].df1 === s.df1 && sum[k].df2 === s.df2 && near(sum[k].pF, s.pDF, 1e-6, 1e-12), `${t} 변화량 ${sum[k].dF} vs ${s.dF}`);
      assert.ok(near(sum[k].dw, s.dw, 1e-9), t + ' DW');
      const a = anova.slice(3 * k, 3 * k + 3);
      assert.deepEqual(a.map(r => r.source), ['회귀', '잔차', '합계']);
      assert.ok(near(a[0].ss, s.ssr, 1e-9) && near(a[1].ss, s.sse, 1e-9) && near(a[2].ss, s.sst, 1e-10), t + ' 제곱합');
      assert.ok(a[0].df === s.dfm && a[1].df === s.dfr && a[2].df === c.n - 1, t + ' 자유도');
      assert.ok(near(a[0].F, s.F, 1e-9) && near(a[0].p, s.pF, 1e-6, 1e-12), t + ' F');
      s.names.forEach((nm, j) => {
        const r = coef[ci++];
        assert.equal(r.name, nm, t);
        assert.equal(r.model, j === 0 ? (k + 1) + '단계' : '', t);
        assert.ok(near(r.B, s.B[j], 1e-8, 1e-9) && near(r.se, s.se[j], 1e-8) && near(r.t, s.t[j], 1e-7) && near(r.p, s.pval[j], 1e-6, 1e-12), `${t} ${nm}`);
        assert.ok(near(r.lo, s.lo[j], 1e-8, 1e-9) && near(r.hi, s.hi[j], 1e-8, 1e-9), `${t} ${nm} 구간`);
        if (j === 0) { assert.equal(r.beta, null); assert.equal(r.vif, null); }
        else assert.ok(near(r.beta, s.beta[j], 1e-8) && near(r.vif, s.vif[j], 1e-7) && near(r.tol, s.tol[j], 1e-7), `${t} ${nm} 베타·VIF`);
      });
    });
    assert.equal(ci, coef.length, w);
  });
});

test('runHierarchical: 1단계만 쓰면 선형 회귀와 같은 결과이고 안내 경고가 나온다', () => {
  const c = ref.cases[0];
  const { dataset, request } = build(c, {});
  const one = H.runHierarchical(dataset, { dependent: ['y'], block1: request.block1, options: {} });
  const lin = Reg.runRegression(dataset, { dependent: ['y'], covariates: request.block1, options: {} });
  const a = table(one, '계수').rows, b = table(lin, '계수').rows;
  assert.deepEqual(a.map(r => r.name), b.map(r => r.name));
  a.forEach((r, i) => assert.ok(near(r.B, b[i].B, 1e-12) && near(r.p, b[i].p, 1e-12), r.name));
  const s = table(one, '모형 요약').rows[0], l = table(lin, '모형 요약').rows[0];
  assert.ok(near(s.r2, l.r2, 1e-12) && near(s.dR2, l.r2, 1e-12) && near(s.dF, table(lin, '분산분석').rows[0].F, 1e-12));
  assert.ok(one.warnings.some(w => w.includes('1단계만 있어')));
});

test('runHierarchical: 표 구성, 투입 변수 안내, 서열 변수 안내, 범주형 기준 범주, 공선성 진단은 마지막 단계, 렌더·내보내기', () => {
  const c = ref.cases[1];
  const { dataset, request } = build(c, { collinearity: true });
  const res = H.runHierarchical(dataset, request);
  assert.equal(res.title, '위계적 회귀');
  assert.deepEqual(res.blocks.map(b => b.caption.split(' (')[0]), ['모형 요약', '분산분석', '계수', '공선성 진단']);
  const notes = table(res, '모형 요약').notes;
  assert.ok(notes[0].startsWith('1단계 투입: grp = 2, grp = 3, age') && notes[1].startsWith('2단계 추가: z1, z2'));
  assert.ok(table(res, '계수').notes.some(n => n.includes('"grp"의 기준 범주: 1')));
  assert.ok(!notes.some(n => n.includes('서열')));
  const diag = table(res, '공선성 진단');
  assert.equal(diag.rows.length, 6);                        // 차원 수 = 마지막 단계의 계수 수 (상수, grp = 2, grp = 3, age, z1, z2)
  assert.equal(diag.caption, '공선성 진단 (마지막 단계 모형)');
  // 서열로 정의하면 숫자 점수로 투입하고 안내한다
  const ord = build(ref.cases[3], {}, { l1: 'ordinal', l3: 'ordinal' });
  const r2 = H.runHierarchical(ord.dataset, ord.request);
  assert.ok(table(r2, '모형 요약').notes.some(n => n.includes('서열 척도 변수(l1, l3)는 숫자 점수로 투입')));
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('R 제곱 변화량') && html.includes('1단계 투입'));
  assert.ok(Ex.resultToGrid(res).rows.length > 20);
});

test('runHierarchical: 결측 사례는 모든 단계에서 같이 제외한다', () => {
  const c = ref.cases[0];
  const vars = c.blocks.flat();
  const rows = [['y'].concat(vars.map(v => v.name))];
  for (let i = 0; i < c.n; i++) rows.push([c.y[i]].concat(vars.map(v => v.values[i])));
  rows[3][vars.length] = '';            // 마지막 단계 변수 하나에 결측
  rows[10][1] = '';                     // 종속변수 결측
  const ds = Data.buildDataset(rows, []);
  const req = { dependent: ['y'], block1: ['age', 'ten'], block2: ['x1'], block3: ['x2', 'x3'], options: {} };
  const res = H.runHierarchical(ds, req);
  assert.ok(table(res, '모형 요약').notes.some(n => n.includes('유효 사례 58건') && n.includes('2건 제외')));
  assert.equal(table(res, '분산분석').rows[2].df, 57);        // 모든 단계가 같은 58사례
  assert.equal(table(res, '분산분석').rows[8].df, 57);
});

test('runHierarchical: 입력 검증과 계산 불가 경우', () => {
  const c = ref.cases[0];
  const { dataset } = build(c, {});
  const run = req => H.runHierarchical(dataset, Object.assign({ dependent: ['y'], options: {} }, req));
  assert.throws(() => run({}), /1단계에 독립변수/);
  assert.throws(() => run({ block2: ['x1'] }), /1단계에 독립변수/);
  assert.throws(() => run({ block1: ['age'], block3: ['x1'] }), /2단계가 비어 있는데/);
  assert.throws(() => run({ block1: ['age'], block2: ['age'] }), /"age"를 1단계 독립변수와\(과\) 2단계/);
  assert.throws(() => run({ block1: ['y'] }), /"y"를 종속변수와\(과\) 1단계/);
  assert.throws(() => run({ dependent: ['y', 'age'], block1: ['ten'] }), /하나만/);
  assert.throws(() => run({ block1: ['age'], options: { conf: 100 } }), /신뢰수준/);
  // 완전 공선성: 2단계에서 앞 변수의 배수를 추가
  const rows = [['y', 'a', 'b']];
  for (let i = 0; i < 12; i++) rows.push([Math.sin(i) * 3 + i, i % 5, 2 * (i % 5)]);
  assert.throws(() => H.runHierarchical(Data.buildDataset(rows, []), { dependent: ['y'], block1: ['a'], block2: ['b'] }), /"b".*2단계 모형.*완전 공선성/);
  // 사례 부족: 1단계는 되지만 2단계 계수가 너무 많음
  const few = [['y', 'a', 'b', 'c'], [1, 2, 3, 5], [2, 1, 4, 3], [3, 4, 2, 6], [5, 3, 1, 2]];
  assert.throws(() => H.runHierarchical(Data.buildDataset(few, []), { dependent: ['y'], block1: ['a'], block2: ['b', 'c'] }), /사례 수\(4\)가 2단계 모형의 계수의 수\(4/);
  // 숫자가 아닌 값(명목이 아닌 변수)
  const txt = [['y', 'a', 'g'], [1, 2, 'x'], [2, 1, 'y'], [3, 4, 'x'], [5, 3, 'y'], [6, 5, 'x']];
  const dsTxt = Data.buildDataset(txt, [{ name: 'g', label: '', scale: 'ordinal', valueLabels: {}, missing: [] }]);
  assert.throws(() => H.runHierarchical(dsTxt, { dependent: ['y'], block1: ['a'], block2: ['g'] }), /숫자가 아닌 값/);
  // 값이 모두 같은 변수
  const flat = [['y', 'a', 'b'], [1, 2, 7], [2, 1, 7], [3, 4, 7], [5, 3, 7], [6, 5, 7]];
  assert.throws(() => H.runHierarchical(Data.buildDataset(flat, []), { dependent: ['y'], block1: ['a'], block2: ['b'] }), /값이 모두 같아/);
});

test('runHierarchical: 완전한 적합이면 경고하고 변화 F 는 비운다', () => {
  const rows = [['y', 'a', 'b']];
  const a = [1, 2, 3, 4, 5, 6, 7, 8], b = [2, 1, 4, 3, 6, 5, 8, 9];
  for (let i = 0; i < 8; i++) rows.push([2 * a[i] + 3 * b[i] + 1, a[i], b[i]]);
  const res = H.runHierarchical(Data.buildDataset(rows, []), { dependent: ['y'], block1: ['a'], block2: ['b'], options: {} });
  assert.ok(res.warnings.some(w => w.includes('2단계 모형의 잔차가 모두 0')));
  const last = table(res, '모형 요약').rows[1];
  assert.equal(last.dF, null);
  assert.equal(last.pF, null);
});
