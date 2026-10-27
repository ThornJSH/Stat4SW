const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const C = require('../src/Stats_Cochran.js');
const { pchisq } = require('../src/Distributions.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'cochran.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-14) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const def = name => ({ name, label: '', scale: 'nominal', valueLabels: {}, missing: [] });

function dsOf(x, names, extra) {
  const nm = names || x[0].map((_, j) => 't' + (j + 1));
  const rows = [nm].concat(x.map(r => r.slice()));
  if (extra) extra(rows);
  return Data.buildDataset(rows, nm.map(def));
}

test('cochranQ: statsmodels cochrans_q (= scipy friedmanchisquare) 와 같다 — Q, 자유도, 유의확률, W, 응답이 달라진 사례 수', () => {
  assert.equal(ref.cases.length, 5);
  ref.cases.forEach(c => {
    const r = C.cochranQ(c.x, pchisq);
    assert.ok(near(r.q, c.q, 1e-12, 1e-12), `${c.name} Q ${r.q} vs ${c.q}`);
    assert.equal(r.df, c.df);
    assert.ok(near(r.p, c.p, 1e-9, 1e-14), `${c.name} p ${r.p} vs ${c.p}`);
    assert.ok(near(r.w, c.w, 1e-12));
    assert.equal(r.informative, c.informative);
    assert.deepEqual(r.colsum, c.colsum);
  });
});

test('runCochran: 표 4개(처리 요약, 측정별 빈도, Q 검정, 쌍별 McNemar)와 그림이 기준값과 같다', () => {
  ref.cases.forEach(c => {
    const names = c.x[0].map((_, j) => '측정' + (j + 1));
    const res = C.runCochran(dsOf(c.x, names), { vars: names, options: { success: '1', posthoc: true, plot: true } });
    const sum = table(res, '처리 요약');
    assert.equal(sum.rows[0].n, c.n); assert.equal(sum.rows[4].n, c.informative);
    const fr = table(res, '측정별 성공 빈도');
    fr.rows.forEach((r, q) => { assert.equal(r.x, c.colsum[q]); assert.equal(r.y, c.n - c.colsum[q]); assert.ok(near(r.pct, c.colsum[q] / c.n, 1e-12)); });
    const t = table(res, 'Cochran Q');
    assert.ok(near(t.rows[0].q, c.q, 1e-12) && near(t.rows[0].p, c.p, 1e-9, 1e-14) && t.rows[0].df === c.k - 1 && near(t.rows[0].w, c.w, 1e-12));
    const pw = table(res, '쌍별 비교');
    const m = c.k * (c.k - 1) / 2;
    assert.equal(pw.rows.length, m);
    pw.rows.forEach((r, q) => {
      const e = c.pairs[q];
      assert.equal(r.bc, e.ps_fail); assert.equal(r.cb, e.pf_succ);
      assert.ok(near(r.p, e.p, 1e-9, 1e-14), `${c.name} 쌍 ${q} p ${r.p} vs ${e.p}`);
      assert.ok(near(r.pAdj, Math.min(1, e.p * m), 1e-9, 1e-14));
    });
    const chart = res.blocks.find(b => b.type === 'chart');
    assert.deepEqual(chart.xTicks, names);
    assert.ok(chart.series[0].points.every((p, q) => near(p[1], c.colsum[q] / c.n, 1e-12)));
  });
  const c0 = ref.cases[0], nm0 = c0.x[0].map((_, j) => 't' + j);
  const res = C.runCochran(dsOf(c0.x, nm0), { vars: nm0, options: {} });
  assert.ok(R.renderResultHtml(res).includes('Cochran') && Ex.resultToGrid(res).rows.length > 10);
});

test('runCochran: 옵션 — 쌍별·그림 끄기, 성공값 지정(0/글자), 결측 목록별 제거', () => {
  const x = ref.cases[0].x;
  const nm = ['a', 'b', 'c'];
  const off = C.runCochran(dsOf(x, nm), { vars: nm, options: { posthoc: false, plot: false } });
  assert.ok(!table(off, '쌍별') && !off.blocks.some(b => b.type === 'chart'));
  // 성공값을 0으로 하면 성공과 실패가 뒤바뀌지만 Q 는 같다 (대칭)
  const flip = C.runCochran(dsOf(x, nm), { vars: nm, options: { success: '0', posthoc: false } });
  assert.ok(near(table(flip, 'Cochran Q').rows[0].q, table(off, 'Cochran Q').rows[0].q, 1e-12));
  assert.equal(table(flip, '측정별').rows[0].x, 30 - table(off, '측정별').rows[0].x);
  // 글자 값
  const txt = Data.buildDataset([nm].concat(x.map(r => r.map(v => (v ? '예' : '아니오')))), nm.map(def));
  const rt = C.runCochran(txt, { vars: nm, options: { success: '예', posthoc: false } });
  assert.ok(near(table(rt, 'Cochran Q').rows[0].q, table(off, 'Cochran Q').rows[0].q, 1e-12));
  // 결측: 하나라도 비면 그 사람 제외
  const miss = C.runCochran(dsOf(x, nm, rows => { rows[1][1] = null; rows[2][3 - 1] = null; }), { vars: nm, options: { posthoc: false } });
  assert.equal(table(miss, '처리 요약').rows[0].n, 28);
  assert.ok(miss.warnings.some(w => /2건을 제외/.test(w)));
});

test('runCochran: 오류와 경고 — 변수 2개, 성공값 없음, 사례 부족, 모두 같은 응답, 응답이 달라진 사례가 적음, 값이 3가지 이상', () => {
  const x = ref.cases[0].x, nm = ['a', 'b', 'c'];
  assert.throws(() => C.runCochran(dsOf(x.map(r => r.slice(0, 2)), ['a', 'b']), { vars: ['a', 'b'], options: {} }), /3개 이상.*McNemar/);
  assert.throws(() => C.runCochran(dsOf(x, nm), { vars: nm, options: { success: '9' } }), /성공으로 볼 값 "9".*변수에 있는 값: /);
  assert.throws(() => C.runCochran(dsOf([[1, 0, 1]], nm), { vars: nm, options: {} }), /1건입니다/);
  const same = C.runCochran(dsOf([[1, 1, 1], [0, 0, 0], [1, 1, 1], [0, 0, 0]], nm), { vars: nm, options: {} });
  assert.ok(same.warnings.some(w => /같은 응답이라 Q 를 계산할 수 없습니다/.test(w)));
  assert.equal(table(same, 'Cochran Q').rows[0].q, null);
  assert.ok(!table(same, '쌍별'));
  const few = C.runCochran(dsOf(ref.cases[3].x, nm), { vars: nm, options: { posthoc: false } });
  assert.ok(few.warnings.some(w => /2건뿐/.test(w)));
  const three = C.runCochran(dsOf([[1, 2, 0], [0, 1, 1], [1, 0, 2], [1, 1, 1]], nm), { vars: nm, options: {} });
  assert.ok(three.warnings.some(w => /값이 3가지/.test(w)));
});

test('큰 문제: Cochran Q 20회 측정 × 5만 명(쌍별 190쌍·그림 포함)이 Node 에서 3초 안에 끝난다', () => {
  let seed = 5;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  const k = 20, n = 50000, names = Array.from({ length: k }, (_, j) => 't' + j);
  const x = [];
  for (let i = 0; i < n; i++) { const z = rnd() - 0.5; x.push(names.map((_, j) => (rnd() + z + j * 0.004 > 0.75 ? 1 : 0))); }
  const ds = dsOf(x, names);
  const t0 = Date.now();
  const res = C.runCochran(ds, { vars: names, options: { posthoc: true, plot: true } });
  const sec = (Date.now() - t0) / 1000;
  assert.equal(table(res, '쌍별').rows.length, 190);
  assert.ok(sec < 3, `Cochran Q ${sec}초`);
});
