const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const Rel = require('../src/Stats_Relation.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'partial.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-12) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const tables = (res, caption) => res.blocks.filter(b => b.type === 'table' && b.caption.startsWith(caption));
const def = name => ({ name, label: '', scale: 'scale', valueLabels: {}, missing: [] });

function dsOf(c, extra) {
  const names = c.vnames.concat(c.cnames);
  const cols = c.V.concat(c.C);
  const rows = [names].concat(cols[0].map((_, i) => cols.map(col => col[i])));
  if (extra) extra(rows);
  return Data.buildDataset(rows, names.map(def));
}

test('partialCorrelation: pingouin·numpy 두 경로가 서로 일치한 기준값과 같다 — 피어슨·스피어만의 편상관, 유의확률, 신뢰구간, 통제 전 상관 (6가지 자료)', () => {
  assert.equal(ref.cases.length, 6);
  ref.cases.forEach(c => {
    ['pearson', 'spearman'].forEach(method => {
      const rank = method === 'spearman';
      const vc = c.V.map(col => rank ? Rel.rankAverage_(col) : col), cc = c.C.map(col => rank ? Rel.rankAverage_(col) : col);
      const res = Rel.partialCorrelation(vc, cc, c.conf);
      const e = c[method], w = `${c.name} ${method}`;
      assert.equal(res.df, e.df);
      assert.equal(res.n, c.n);
      for (let i = 0; i < c.k; i++) {
        for (let j = 0; j < c.k; j++) {
          assert.ok(near(res.r0[i][j], e.zero[i][j], 1e-9, 1e-12), `${w} 통제 전 ${i},${j}`);
          assert.ok(near(res.r[i][j], e.r[i][j], 1e-8, 1e-10), `${w} 편상관 ${i},${j}: ${res.r[i][j]} vs ${e.r[i][j]}`);
          if (i !== j) {
            assert.ok(near(res.p[i][j], e.p[i][j], 1e-6, 1e-13), `${w} p ${i},${j}`);
            assert.ok(near(res.lo[i][j], e.lo[i][j], 1e-7, 1e-9) && near(res.hi[i][j], e.hi[i][j], 1e-7, 1e-9), `${w} 구간 ${i},${j}`);
          }
        }
      }
    });
  });
  // 억제 사례는 통제 전후 부호가 다르다
  const sup = ref.cases[4];
  assert.ok(sup.pearson.zero[0][1] * sup.pearson.r[0][1] < 0);
});

test('runCorrelation(통제변수): 표 구성, 피어슨·스피어만, 통제 전 → 통제 후 비교표, 신뢰구간 옵션, 렌더·내보내기', () => {
  const c = ref.cases[1];                                  // 변수 4개, 통제변수 2개
  const ds = dsOf(c);
  const res = Rel.runCorrelation(ds, { vars: c.vnames, controls: c.cnames, options: { pearson: true, spearman: true, ci: true, conf: '95' } });
  assert.equal(res.title, '편상관분석 (통제변수: z1, z2)');
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption.split(' (')[0]),
    ['피어슨 편상관 행렬', '피어슨 편상관: 통제 전 → 통제 후 비교', '스피어만 순위 편상관 행렬', '스피어만 순위 편상관: 통제 전 → 통제 후 비교']);
  const mat = tables(res, '피어슨 편상관 행렬')[0];
  assert.equal(mat.rows.length, 4 * 3);                     // 변수마다 (계수, 유의확률, 자유도)
  assert.deepEqual(mat.rows.slice(0, 3).map(r => r.stat), ['피어슨 편상관 계수', '유의확률(양측)', '자유도']);
  const e = c.pearson;
  assert.ok(near(mat.rows[0].c1, e.r[0][1], 1e-8) && near(mat.rows[1].c1, e.p[0][1], 1e-6, 1e-13) && mat.rows[2].c1 === e.df && mat.rows[0].c0 === 1 && mat.rows[1].c0 === null);
  assert.ok(near(mat.rows[3 * 3].c2, e.r[3][2], 1e-8));
  const cmp = table(res, '피어슨 편상관: 통제 전');
  assert.equal(cmp.rows.length, 6);
  assert.deepEqual([cmp.rows[0].a, cmp.rows[0].b], ['a', 'b']);
  assert.ok(near(cmp.rows[0].r0, e.zero[0][1], 1e-9) && near(cmp.rows[0].r, e.r[0][1], 1e-8) && near(cmp.rows[0].diff, e.r[0][1] - e.zero[0][1], 1e-8) && near(cmp.rows[0].lo, e.lo[0][1], 1e-7, 1e-9) && near(cmp.rows[0].hi, e.hi[0][1], 1e-7, 1e-9));
  assert.ok(cmp.columns.some(x => x.key === 'lo') && cmp.notes.some(n => n.includes('Fisher z')));
  const sp = table(res, '스피어만 순위 편상관: 통제 전');
  assert.ok(near(sp.rows[1].r, c.spearman.r[0][2], 1e-8) && near(sp.rows[1].p, c.spearman.p[0][2], 1e-6, 1e-13));
  assert.ok(!table(Rel.runCorrelation(ds, { vars: c.vnames, controls: c.cnames, options: { pearson: true } }), '피어슨 편상관: 통제 전').columns.some(x => x.key === 'lo'));
  assert.ok(R.renderResultHtml(res).includes('편상관 행렬') && R.renderResultHtml(res).includes('통제 전'));
  assert.ok(Ex.resultToGrid(res).rows.length > 30);
});

test('runCorrelation(통제변수): 잔차 산점도, 결측 목록별 제외, 쌍별 옵션이어도 목록별 제거, 통제변수가 없으면 기존 결과 그대로', () => {
  const c = ref.cases[0];
  const ds = dsOf(c);
  const res = Rel.runCorrelation(ds, { vars: c.vnames, controls: c.cnames, options: { pearson: true, scatter: true } });
  const sc = res.blocks.filter(b => b.type === 'chart');
  assert.equal(sc.length, 3);                               // 변수 3개 → 3쌍
  assert.ok(sc[0].caption.includes('통제 후 잔차 산점도: x × y'));
  // 잔차 산점도의 상관은 편상관과 같다
  const pts = sc[0].points, mx = pts.reduce((s, p) => s + p[0], 0) / pts.length, my = pts.reduce((s, p) => s + p[1], 0) / pts.length;
  let sxy = 0, sxx = 0, syy = 0;
  pts.forEach(p => { sxy += (p[0] - mx) * (p[1] - my); sxx += (p[0] - mx) ** 2; syy += (p[1] - my) ** 2; });
  assert.ok(near(sxy / Math.sqrt(sxx * syy), c.pearson.r[0][1], 1e-8) && Math.abs(mx) < 1e-10 && Math.abs(my) < 1e-10);
  assert.ok(near(sc[0].r, c.pearson.r[0][1], 1e-8) && pts.length === c.n);
  const rk = Rel.runCorrelation(ds, { vars: c.vnames, controls: c.cnames, options: { spearman: true, scatter: true } });
  assert.ok(rk.blocks.find(b => b.type === 'chart').caption.startsWith('순위의 통제 후 잔차 산점도'));
  // 결측은 목록별로 제외
  const dm = dsOf(c, rows => { rows[3][0] = ''; rows[10][3] = ''; });
  const rm = Rel.runCorrelation(dm, { vars: c.vnames, controls: c.cnames, options: { pearson: true, listwise: false } });
  assert.ok(rm.warnings.some(w => w.includes('2건을 제외') && w.includes('58건')));
  assert.ok(table(rm, '피어슨 편상관 행렬').notes.some(n => n.includes('목록별 결측 제거')));
  assert.equal(table(rm, '피어슨 편상관 행렬').rows[2].c1, 58 - 2 - 1);
  // 통제변수 칸이 비어 있으면 기존 상관분석 (제목과 표 이름 그대로)
  const plain = Rel.runCorrelation(ds, { vars: c.vnames, controls: [], options: { pearson: true } });
  assert.equal(plain.title, '상관분석');
  assert.ok(table(plain, '피어슨 상관 행렬'));
});

test('runCorrelation(통제변수): 입력 검증과 계산 불가 경우', () => {
  const c = ref.cases[0];
  const ds = dsOf(c);
  const run = (req, d) => Rel.runCorrelation(d || ds, Object.assign({ vars: ['x', 'y'], controls: ['z'], options: { pearson: true } }, req));
  assert.throws(() => run({ vars: ['x'] }), /2개 이상/);
  assert.throws(() => run({ controls: ['x'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => run({ options: {} }), /상관계수 종류/);
  assert.throws(() => run({ options: { pearson: true, conf: 100 } }), /신뢰수준/);
  const wide = [['a', 'b'].concat(Array.from({ length: 11 }, (_, i) => 'z' + i))];
  for (let i = 0; i < 40; i++) wide.push([i, Math.sin(i), ...Array.from({ length: 11 }, (_, j) => Math.cos(i * (j + 2)) + i * 0.01 * j)]);
  assert.throws(() => Rel.runCorrelation(Data.buildDataset(wide, []), { vars: ['a', 'b'], controls: Array.from({ length: 11 }, (_, i) => 'z' + i), options: { pearson: true } }), /10개까지/);
  const few = Data.buildDataset([['x', 'y', 'z1', 'z2'], [1, 2, 3, 4], [2, 1, 5, 2], [3, 5, 2, 1], [4, 3, 1, 6]], []);
  assert.throws(() => Rel.runCorrelation(few, { vars: ['x', 'y'], controls: ['z1', 'z2'], options: { pearson: true } }), /5건 이상 필요/);
  const flat = Data.buildDataset([['x', 'y', 'z'], [1, 2, 7], [2, 1, 7], [3, 5, 7], [4, 3, 7], [5, 6, 7], [6, 4, 7]], []);
  assert.throws(() => Rel.runCorrelation(flat, { vars: ['x', 'y'], controls: ['z'], options: { pearson: true } }), /"z"의 값이 모두 같아/);
  const flatVar = Data.buildDataset([['x', 'y', 'z'], [1, 7, 7.5], [2, 7, 1], [3, 7, 2], [4, 7, 8], [5, 7, 3], [6, 7, 5]], []);
  assert.throws(() => Rel.runCorrelation(flatVar, { vars: ['x', 'y'], controls: ['z'], options: { pearson: true } }), /"y"의 값이 모두 같아/);
  const coll = [['x', 'y', 'z1', 'z2']];
  for (let i = 0; i < 12; i++) coll.push([Math.sin(i) * 3 + i, Math.cos(i * 2) * 2 + i, i % 4 + 0.1 * i, 2 * (i % 4 + 0.1 * i) + 1]);
  assert.throws(() => Rel.runCorrelation(Data.buildDataset(coll, []), { vars: ['x', 'y'], controls: ['z1', 'z2'], options: { pearson: true } }), /통제변수들이 서로 선형결합/);
  // 통제변수로 완전히 설명되는 변수는 편상관이 계산되지 않는다
  const exact = [['x', 'y', 'z']];
  for (let i = 0; i < 15; i++) exact.push([Math.sin(i * 1.3) * 4 + i * 0.2, 2 * (Math.cos(i * 0.7) * 3 + i * 0.1) + 1, Math.cos(i * 0.7) * 3 + i * 0.1]);
  const rex = Rel.runCorrelation(Data.buildDataset(exact, []), { vars: ['x', 'y'], controls: ['z'], options: { pearson: true } });
  assert.ok(Number.isNaN(table(rex, '피어슨 편상관 행렬').rows[0].c1));
  const txt = Data.buildDataset([['x', 'y', 'z'], [1, 2, 'a'], [2, 1, 'b'], [3, 5, 'c'], [4, 3, 'd'], [5, 6, 'e']], [{ name: 'z', label: '', scale: 'ordinal', valueLabels: {}, missing: [] }]);
  assert.throws(() => Rel.runCorrelation(txt, { vars: ['x', 'y'], controls: ['z'], options: { pearson: true } }), /숫자가 아닌 값/);
});
