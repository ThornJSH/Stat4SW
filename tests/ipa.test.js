const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const I = require('../src/Stats_Ipa.js');
const Ch = require('../src/Charts.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'ipa.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-12) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const def = name => ({ name, label: '', scale: 'ordinal', valueLabels: {}, missing: [] });

function dsOf(c) {
  const k = c.k, inames = Array.from({ length: k }, (_, j) => '중요' + (j + 1)), pnames = Array.from({ length: k }, (_, j) => '수행' + (j + 1));
  const rows = [inames.concat(pnames)].concat(c.imp.map((r, i) => r.concat(c.perf[i])));
  return { ds: Data.buildDataset(rows, inames.concat(pnames).map(def)), inames, pnames };
}
const optsOf = (s, listwise) => {
  const o = { baseline: s.baseline, impShift: String(s.impShift || 0), perfShift: String(s.perfShift || 0), gap: true, plot: true, listwise: !!listwise };
  if (s.baseline === 'custom') { o.impLine = String(s.impLine); o.perfLine = String(s.perfLine); }
  return o;
};

test('runIpa: numpy·scipy 기준값 — 문항 평균·표준편차, 기준선(평균·중앙값·이동·직접 지정), 사분면, 경계, 격차 검정 (3가지 자료 × 시나리오 5가지)', () => {
  assert.equal(ref.cases.length, 3);
  ref.cases.forEach(c => {
    const { ds, inames, pnames } = dsOf(c);
    c.results.forEach((e, si) => {
      const m = `${c.name} #${si}`;
      const res = I.runIpa(ds, { importance: inames, performance: pnames, options: optsOf(e.scenario, e.listwise) });
      const it = table(res, '문항별'), ln = table(res, '기준선'), gt = table(res, '격차 검정');
      assert.equal(it.rows.length, c.k);
      assert.ok(near(ln.rows[0].line, e.line_imp, 1e-12) && near(ln.rows[1].line, e.line_perf, 1e-12), `${m} 기준선`);
      assert.ok(near(ln.rows[0].base, e.base_imp, 1e-12) && near(ln.rows[1].base, e.base_perf, 1e-12), `${m} 기준 계산값`);
      e.items.forEach((x, j) => {
        const r = it.rows[j], g = gt.rows[j];
        assert.equal(r.ni, x.imp.n); assert.equal(r.np, x.perf.n);
        assert.ok(near(r.mi, x.imp.mean, 1e-12) && near(r.si, x.imp.sd, 1e-12) && near(r.mp, x.perf.mean, 1e-12) && near(r.sp, x.perf.sd, 1e-12), `${m} 문항 ${j} 평균`);
        assert.ok(near(r.gap, x.gap, 1e-12));
        assert.equal(r.quad, x.quadrant, `${m} 문항 ${j} 사분면`);
        assert.equal(r.edge === '★', x.edge, `${m} 문항 ${j} 경계`);
        assert.equal(g.n, x.pairs);
        ['mdiff:md', 'sddiff:sd', 't:t', 'df:df', 'dz:dz'].forEach(s => { const [a, b] = s.split(':'); assert.ok(near(g[b], x[a], 1e-9), `${m} 문항 ${j} ${a}: ${g[b]} vs ${x[a]}`); });
        assert.ok(near(g.p, x.p, 1e-8, 1e-14), `${m} 문항 ${j} p`);
        assert.ok(near(g.pa, x.p === null ? null : Math.min(1, x.p * c.k), 1e-8, 1e-14));
      });
      const quad = table(res, '사분면별');
      assert.deepEqual(quad.rows.map(r => r.n), [0, 1, 2, 3].map(q => e.items.filter(x => x.quadrant === ['집중 개선', '현상 유지', '낮은 우선순위', '과잉 투자'][q]).length));
    });
  });
});

test('runIpa: 그림 — 점 수, 기준선 위치, 사분면 이름, 이름, 내보내기', () => {
  const c = ref.cases[0], e = c.results[1];
  const { ds, inames, pnames } = dsOf(c);
  const res = I.runIpa(ds, { importance: inames, performance: pnames, options: optsOf(e.scenario) });
  const chart = res.blocks.find(b => b.type === 'chart');
  assert.equal(chart.kind, 'ipa');
  assert.equal(chart.points.length, c.k);
  assert.ok(near(chart.xLine, e.line_perf, 1e-12) && near(chart.yLine, e.line_imp, 1e-12));
  chart.points.forEach((p, j) => { assert.ok(near(p.x, e.items[j].perf.mean, 1e-12) && near(p.y, e.items[j].imp.mean, 1e-12)); assert.equal(['집중 개선', '현상 유지', '낮은 우선순위', '과잉 투자'][p.q], e.items[j].quadrant); });
  const html = R.renderResultHtml(res);
  assert.equal((html.match(/<circle/g) || []).length, c.k);
  ['집중 개선', '현상 유지', '낮은 우선순위', '과잉 투자', '중요도 기준', '수행도 기준'].forEach(s => assert.ok(html.includes(s), s));
  assert.ok(chart.points.every(p => html.includes('>' + p.label + '<')));
  assert.ok(Ex.resultToGrid(res, undefined, [{ width: 600, height: 480 }]).images.length === 1);
  const off = I.runIpa(ds, { importance: inames, performance: pnames, options: { gap: false, plot: false } });
  assert.ok(!table(off, '격차') && !off.blocks.some(b => b.type === 'chart'));
});

test('ipaChartBlock: 점이 많을 때 이름이 겹치지 않도록 배치하고 데이터가 없으면 안내 문구', () => {
  const pts = Array.from({ length: 20 }, (_, i) => ({ label: '문항' + (i + 1), x: 3 + (i % 5) * 0.05, y: 3.5 + Math.floor(i / 5) * 0.05, q: i % 4 }));
  const h = Ch.renderChartHtml(Ch.ipaChartBlock('빽빽한 그림', pts, 3.1, 3.55));
  assert.equal((h.match(/<circle/g) || []).length, 20);
  assert.ok(Ch.renderChartHtml(Ch.ipaChartBlock('빈 그림', [], 3, 3)).includes('그릴 데이터가 없습니다'));
});

test('runIpa: 오류와 경고 — 문항 수, 개수 불일치, 같은 변수, 직접 지정 값 누락, 숫자 아님, 척도 범위 경고, 목록별 제거 건수', () => {
  const c = ref.cases[2], { ds, inames, pnames } = dsOf(c);
  const run = (imp, perf, options) => I.runIpa(ds, { importance: imp, performance: perf, options: options || {} });
  assert.throws(() => run(inames.slice(0, 2), pnames.slice(0, 2)), /3개 이상/);
  assert.throws(() => run(inames, pnames.slice(0, 3)), /개수가 같아야/);
  assert.throws(() => run(inames, [inames[0]].concat(pnames.slice(1))), /동시에 쓸 수 없습니다/);
  assert.throws(() => run(inames, pnames, { baseline: 'custom', impLine: '3' }), /직접 지정.*모두 입력/);
  assert.throws(() => run(inames, pnames, { impShift: 'abc' }), /숫자여야/);
  const c0 = ref.cases[1], d0 = dsOf(c0);
  const lw = I.runIpa(d0.ds, { importance: d0.inames, performance: d0.pnames, options: { listwise: true } });
  assert.ok(lw.warnings.some(w => /목록별 제거로 \d+건을 제외/.test(w)));
  const rows = [['a', 'b', 'c', 'x', 'y', 'z']].concat([[5, 4, 3, 9, 8, 7], [4, 4, 3, 8, 9, 7], [5, 3, 2, 9, 7, 6], [4, 5, 3, 7, 8, 9]]);
  const wide = I.runIpa(Data.buildDataset(rows, ['a', 'b', 'c', 'x', 'y', 'z'].map(def)), { importance: ['a', 'b', 'c'], performance: ['x', 'y', 'z'], options: {} });
  assert.ok(wide.warnings.some(w => /척도 범위가 다른/.test(w)));
  const txt = Data.buildDataset([['a', 'b', 'c', 'x', 'y', 'z'], ['가', 1, 2, 3, 4, 5]], ['a', 'b', 'c', 'x', 'y', 'z'].map(def));
  assert.throws(() => I.runIpa(txt, { importance: ['a', 'b', 'c'], performance: ['x', 'y', 'z'], options: {} }), /숫자가 아닌 값/);
});

test('ipaQuadrant: 기준선과 같으면 높음으로 본다', () => {
  assert.equal(I.ipaQuadrant(4, 2, 3, 3), 0);
  assert.equal(I.ipaQuadrant(4, 4, 3, 3), 1);
  assert.equal(I.ipaQuadrant(2, 2, 3, 3), 2);
  assert.equal(I.ipaQuadrant(2, 4, 3, 3), 3);
  assert.equal(I.ipaQuadrant(3, 3, 3, 3), 1);
});

test('큰 문제: IPA 30쌍 × 5만 건(격차 검정·그림 포함)이 Node 에서 5초 안에 끝난다', () => {
  let seed = 13;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  const k = 30, names = [], rows = [];
  for (let j = 0; j < k; j++) names.push('i' + j);
  for (let j = 0; j < k; j++) names.push('p' + j);
  rows.push(names);
  for (let i = 0; i < 50000; i++) rows.push(names.map((_, j) => Math.max(1, Math.min(5, Math.round(3 + (j % k) * 0.03 + (rnd() + rnd() - 1) * 2)))));
  const ds = Data.buildDataset(rows, names.map(def));
  const t0 = Date.now();
  const res = I.runIpa(ds, { importance: names.slice(0, k), performance: names.slice(k), options: {} });
  const sec = (Date.now() - t0) / 1000;
  assert.equal(table(res, '문항별').rows.length, 30);
  assert.ok(sec < 5, `IPA ${sec}초`);
});
