const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const C = require('../src/Stats_Cart.js');
const Ch = require('../src/Charts.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'cart.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-10, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));

const SETUP = [
  { key: 'A', names: ['x1', 'x2', 'noise', 'c'] }, { key: 'A', names: ['x1', 'x2', 'noise', 'c'] }, { key: 'A', names: ['x1', 'x2', 'noise', 'c'] },
  { key: 'B', names: ['x1', 'x2', 'c'] }, { key: 'C', names: ['x1', 'x2', 'c'] }
];
function pieces(ci) {
  const { key, names } = SETUP[ci], d = ref.data[key];
  const P = names.map(nm => ({ name: nm, type: nm === 'c' ? 'cat' : 'num', values: d[nm] }));
  return { P, y: d.y, names, d };
}

test('cartBuild·cartPruneSequence·cartCrossValidate: 파이썬 독립 구현과 같은 나무, 가지치기열, 교차검증 오차, 선택, 최종 나무 — 5가지 설정', () => {
  assert.equal(ref.cases.length, 5);
  ref.cases.forEach((c, ci) => {
    const { P, y, names } = pieces(ci);
    const idx = y.map((_, i) => i);
    const T = C.cartBuild(P, y, c.mode, c.C, c.params, idx);
    // 키운 나무
    assert.equal(T.nodes.length, c.full.length, c.name + ' 노드 수');
    c.full.forEach((e, k) => {
      const nd = T.nodes[k];
      assert.equal(nd.n, e.n, `${c.name} 노드 ${e.id} N`);
      assert.equal(nd.left, e.left); assert.equal(nd.right, e.right); assert.equal(nd.depth, e.depth);
      assert.ok(near(nd.pred, e.pred, 1e-9), `${c.name} 노드 ${e.id} 예측`);
      assert.ok(near(nd.err, e.err, 1e-9) || Math.abs(nd.err - e.err) < 1e-9, `${c.name} 노드 ${e.id} 위험`);
      if (e.split === null) assert.equal(nd.split, null, `${c.name} 노드 ${e.id} 잎`);
      else {
        assert.equal(nd.split.j, e.split.j, `${c.name} 노드 ${e.id} 변수`);
        assert.equal(nd.split.kind, e.split.kind);
        if (e.split.kind === 'num') assert.ok(near(nd.split.arg, e.split.arg, 1e-12), `${c.name} 노드 ${e.id} 임계값`);
        else assert.deepEqual(nd.split.arg, e.split.arg, `${c.name} 노드 ${e.id} 범주 집합`);
        assert.ok(near(nd.split.gain, e.split.gain, 1e-9), `${c.name} 노드 ${e.id} 이득 ${nd.split.gain} vs ${e.split.gain}`);
      }
    });
    // 가지치기열
    const seq = C.cartPruneSequence(T);
    assert.deepEqual(seq.leaves, c.leaves, c.name + ' 잎 수열');
    assert.equal(seq.alphas.length, c.alphas.length);
    c.alphas.forEach((a, k) => assert.ok(near(seq.alphas[k], a, 1e-9) || Math.abs(seq.alphas[k] - a) < 1e-12, `${c.name} α${k + 1}: ${seq.alphas[k]} vs ${a}`));
    // 교차검증과 선택
    let k = 0;
    if (c.rule !== 'none') {
      const cv = C.cartCrossValidate(P, y, c.mode, c.C, c.params, seq.alphas, c.seed, Math.min(c.nfold, y.length));
      c.xerror.forEach((v, q) => { assert.ok(near(cv.xerror[q], v, 1e-9) || Math.abs(cv.xerror[q] - v) < 1e-12, `${c.name} xerror ${q}: ${cv.xerror[q]} vs ${v}`); assert.ok(near(cv.se[q], c.se[q], 1e-8) || Math.abs(cv.se[q] - c.se[q]) < 1e-12, `${c.name} se ${q}`); });
      k = C.cartChoose(cv.xerror, cv.se, c.rule);
    }
    assert.equal(k, c.selected, c.name + ' 선택');
    // 최종 나무
    const fin = C.cartExport(T, seq.sets[k]);
    assert.equal(fin.length, c.final.length, c.name + ' 최종 노드 수');
    fin.forEach((nd, q) => {
      const e = c.final[q];
      assert.deepEqual([nd.id, nd.n, nd.left, nd.right, nd.depth], [e.id, e.n, e.left, e.right, e.depth], `${c.name} 최종 노드 ${e.id}`);
      assert.ok(near(nd.pred, e.pred, 1e-9));
      if (e.split) assert.ok(nd.split && nd.split.j === e.split.j && near(nd.split.gain, e.split.gain, 1e-9));
      else assert.equal(nd.split, null);
    });
  });
});

function makeDs(ci, extra) {
  const { names, d } = pieces(ci);
  const rows = [['y'].concat(names)];
  d.y.forEach((v, i) => rows.push([v].concat(names.map(nm => d[nm][i]))));
  if (extra) extra(rows);
  const defs = [{ name: 'c', label: '구분', scale: 'nominal', valueLabels: ci === 3 ? { '0': 'a', '1': 'b', '2': 'c', '3': 'd' } : { '0': 'A', '1': 'B', '2': 'C' }, missing: [] }];
  defs.push({ name: 'y', label: ci === 4 ? '점수' : '결과', scale: ci === 4 ? 'scale' : 'nominal', valueLabels: {}, missing: [] });
  return Data.buildDataset(rows, defs);
}
const reqFor = (ci, extra) => {
  const c = ref.cases[ci], { names } = pieces(ci);
  return Object.assign({ dependent: ['y'], covariates: names.filter(n => n !== 'c'), factors: ['c'], options: { maxDepth: String(c.params.maxDepth), minSplit: String(c.params.minSplit), minLeaf: String(c.params.minLeaf), prune: c.rule, folds: String(c.nfold), seed: String(c.seed) } }, extra);
};

test('runCart(분류): 표 구성, 가지치기 곡선, 나무 구조, 그림, 변수 중요도, 분류표가 기준값과 같다', () => {
  const c = ref.cases[0];
  const res = C.runCart(makeDs(0), reqFor(0));
  assert.equal(res.title, '의사결정나무');
  assert.deepEqual(table(res, '처리 요약').rows.map(r => r.n), [200, 0, 200]);
  const curve = table(res, '가지치기 곡선').rows;
  assert.equal(curve.length, c.leaves.length);
  curve.forEach((r, k) => { assert.equal(r.leaves, c.leaves[k]); assert.ok(near(r.cv, c.xerror[k], 1e-9)); assert.equal(!!r.sel, k === c.selected); });
  const tr = table(res, '나무 구조').rows;
  assert.equal(tr.length, c.final.length);
  assert.equal(tr[0].rule, '(전체)');
  assert.ok(tr.every((r, k) => r.n === c.final[k].n && (r.split === '잎') === (c.final[k].split === null)));
  const imp = table(res, '변수 중요도').rows;
  Object.keys(c.importance).forEach(nm => assert.ok(near(imp.find(r => r.name === nm).gain, c.importance[nm], 1e-9), '중요도 ' + nm));
  const conf = table(res, '분류표');
  c.confusion.forEach((row, a) => row.forEach((v, b) => assert.equal(conf.rows[a]['p' + b], v)));
  assert.ok(near(conf.rows[2].pct, c.accuracy * 100, 1e-9));
  assert.ok(conf.notes.some(n => n.includes('교차검증 정확도')) && conf.notes.some(n => n.includes('민감도')));
  const chart = res.blocks.find(b => b.type === 'chart' && b.kind === 'tree');
  assert.equal(chart.nodes.length, c.final.length);
  assert.ok(Ch.renderChartHtml(chart).includes('<svg') && Ch.renderChartHtml(chart).includes('예측'));
  assert.ok(res.blocks.some(b => b.type === 'chart' && b.kind === 'line'));
  assert.ok(R.renderResultHtml(res).includes('나무 구조'));
  assert.ok(Ex.resultToGrid(res).rows.length > 20);
});

test('runCart(3분류, 회귀): 노드 규칙 문구, 범주 표시, 회귀 성능이 기준값과 같다', () => {
  const b = ref.cases[3];
  const res = C.runCart(makeDs(3), reqFor(3));
  assert.deepEqual(table(res, '종속변수 분포').rows.map(r => r.orig), ['0', '1', '2']);
  const rules = table(res, '나무 구조').rows.map(r => r.rule + ' | ' + r.split).join('\n');
  assert.ok(/구분 ∈ \{[a-d, ]+\}/.test(rules) || /x1 ≤/.test(rules), rules);
  assert.equal(table(res, '분류표').rows[3].pct !== null, true);
  const r = ref.cases[4];
  const reg = C.runCart(makeDs(4), reqFor(4));
  const perf = table(reg, '모형 성능').rows;
  assert.ok(near(perf[0].value, r.r2, 1e-9) && near(perf[1].value, r.rmse, 1e-9));
  assert.ok(near(perf[2].value, Math.sqrt(r.xerror[r.selected]), 1e-9));
  assert.ok(table(reg, '나무 구조').columns.some(c => c.label === '평균'));
  assert.ok(near(table(res, '분류표').rows[3].pct / 100, b.accuracy, 1e-9));
});

test('runCart: 저장 열(나무예측, 나무노드), 결측 제외, 시드 재현', () => {
  const rows = ds => ds;
  const d1 = makeDs(0, r => { r[6][1] = ''; });
  const res = C.runCart(d1, reqFor(0, { options: { maxDepth: '3', prune: 'cv1se', seed: '3', save: true } }));
  const [pred, node] = res.saveColumns;
  assert.deepEqual([pred.name, node.name], ['나무예측', '나무노드']);
  assert.equal(pred.values.length, 200);
  assert.equal(pred.values[5], null);
  assert.ok(pred.values.filter(v => v !== null).every(v => v === 0 || v === 1));
  assert.equal(pred.valueLabels, '0=0, 1=1');
  const nodeIds = new Set(table(res, '나무 구조').rows.map(r => r.id));
  assert.ok(node.values.filter(v => v !== null).every(v => nodeIds.has(v)));
  assert.ok(res.warnings.some(w => w.includes('1건을 제외')));
  const again = C.runCart(d1, reqFor(0, { options: { maxDepth: '3', prune: 'cv1se', seed: '3', save: true } }));
  assert.deepEqual(again.saveColumns[1].values, node.values);
  const reg = C.runCart(makeDs(4), reqFor(4, { options: { prune: 'none', save: true } }));
  assert.ok(reg.saveColumns[0].values.every(v => typeof v === 'number') && reg.saveColumns[0].scale === '연속');
  assert.equal(C.runCart(makeDs(0), reqFor(0)).saveColumns, undefined);
});

test('runCart: 입력 검증과 경고', () => {
  const run = (ci, extra, o) => C.runCart(makeDs(ci), reqFor(ci, Object.assign({}, extra, o ? { options: Object.assign({}, reqFor(ci).options, o) } : {})));
  assert.throws(() => run(0, { covariates: [], factors: [] }), /예측변수/);
  assert.throws(() => run(0, { covariates: ['y'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => run(0, {}, { maxDepth: '0' }), /최대 깊이/);
  assert.throws(() => run(0, {}, { maxDepth: '9' }), /1~8/);
  assert.throws(() => run(0, {}, { minSplit: '1' }), /분할 최소/);
  assert.throws(() => run(0, {}, { minLeaf: '0' }), /잎 최소/);
  assert.throws(() => run(0, {}, { folds: '1' }), /폴드/);
  assert.throws(() => run(0, {}, { seed: 'a' }), /시드/);
  const many = Data.buildDataset([['y', 'x']].concat(Array.from({ length: 30 }, (_, i) => [i % 12, i])), [{ name: 'y', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]);
  assert.throws(() => C.runCart(many, { dependent: ['y'], covariates: ['x'], options: {} }), /12가지[\s\S]*최대 10가지|최대 10가지/);
  const one = Data.buildDataset([['y', 'x']].concat(Array.from({ length: 20 }, (_, i) => [1, i])), [{ name: 'y', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]);
  assert.throws(() => C.runCart(one, { dependent: ['y'], covariates: ['x'], options: {} }), /한 가지뿐/);
  const small = Data.buildDataset([['y', 'x']].concat(Array.from({ length: 8 }, (_, i) => [i % 2, i])), [{ name: 'y', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]);
  assert.throws(() => C.runCart(small, { dependent: ['y'], covariates: ['x'], options: {} }), /10건 이상/);
  const cats = Data.buildDataset([['y', 'g']].concat(Array.from({ length: 40 }, (_, i) => [i % 2, i % 14])), [{ name: 'y', label: '', scale: 'nominal', valueLabels: {}, missing: [] }, { name: 'g', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]);
  assert.throws(() => C.runCart(cats, { dependent: ['y'], factors: ['g'], options: {} }), /범주가 14개/);
  // 예측변수가 의미 없으면 (잡음만) 루트만 남아 경고가 날 수 있고, 항상 불안정성 경고는 있다
  const res = run(0, {}, {});
  assert.ok(res.warnings.some(w => w.includes('불안정성')));
  const stump = C.runCart(makeDs(0), reqFor(0, { options: { maxDepth: '1', prune: 'none' } }));
  assert.equal(table(stump, '나무 요약').rows[3].value, 2);
  // 사례가 5,000건을 넘으면 거부
  const big = Data.buildDataset([['y', 'x']].concat(Array.from({ length: 5100 }, (_, i) => [i % 2, i])), [{ name: 'y', label: '', scale: 'nominal', valueLabels: {}, missing: [] }]);
  assert.throws(() => C.runCart(big, { dependent: ['y'], covariates: ['x'], options: {} }), /5000건/);
});

test('treeChartBlock: 렌더 결과에 모든 노드와 가지 문구가 들어간다', () => {
  const nodes = [
    { id: 1, left: 2, right: 3, edge: '', lines: ['x ≤ 1?', 'N=10', '예측: A'], cls: 0 },
    { id: 2, left: 0, right: 0, edge: 'x ≤ 1', lines: ['N=6', '예측: A'], cls: 0 },
    { id: 3, left: 0, right: 0, edge: 'x > 1', lines: ['N=4', '예측: B'], cls: 1 }
  ];
  const html = Ch.renderChartHtml(Ch.treeChartBlock('테스트 나무', nodes));
  assert.ok(html.includes('x ≤ 1') && html.includes('x &gt; 1') && html.includes('N=6') && html.includes('N=4'));
  assert.equal((html.match(/<rect /g) || []).length, 4);                              // 배경 + 노드 3개
  assert.ok(Ch.renderChartHtml({ type: 'chart', kind: 'tree', caption: '빈', nodes: [] }).includes('그릴 데이터가 없습니다'));
});
