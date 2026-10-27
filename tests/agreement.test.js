const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const A = require('../src/Stats_Agreement.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'agreement.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-12) => (a === null || e === null ? a === e : Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const colsOf = Y => Y[0].map((_, j) => Y.map(r => r[j]));
const scaleDef = (name, scale, vl) => ({ name, label: '', scale: scale || 'scale', valueLabels: vl || {}, missing: [] });

function dsOfCols(columns, names, scale, vl, extra) {
  const rows = [names].concat(columns[0].map((_, i) => columns.map(c => c[i])));
  if (extra) extra(rows);
  return Data.buildDataset(rows, names.map(n => scaleDef(n, scale, vl)));
}

test('cohenKappa: sklearn·statsmodels·numpy가 서로 일치한 기준값과 같다 — kappa, 표준오차, 영가설 z, 유의확률, 신뢰구간, 가중치 3가지 (3가지 자료)', () => {
  assert.equal(ref.kappa.length, 3);
  ref.kappa.forEach(c => {
    ['none', 'linear', 'quadratic'].forEach(w => {
      const e = c[w], r = A.cohenKappa(c.table, w, 95);
      const m = `${c.name} ${w}`;
      assert.equal(r.n, e.n);
      ['po', 'pe', 'kappa', 'se', 'z', 'lo', 'hi', 'kmax'].forEach(key => assert.ok(near(r[key], e[key], 1e-9, 1e-12), `${m} ${key}: ${r[key]} vs ${e[key]}`));
      assert.ok(near(r.p, e.p, 1e-7, 1e-14), `${m} p`);
    });
  });
  // 신뢰수준을 바꾸면 구간이 달라진다
  const t = ref.kappa[0].table;
  const k90 = A.cohenKappa(t, 'none', 90), k99 = A.cohenKappa(t, 'none', 99);
  assert.ok(k99.hi - k99.lo > k90.hi - k90.lo);
  // 우연 일치가 100% (한 범주만 사용)이면 kappa 는 null
  const degenerate = A.cohenKappa([[10, 0], [0, 0]], 'none', 95);
  assert.equal(degenerate.kappa, null);
  // 완전 일치
  const perfect = A.cohenKappa([[5, 0], [0, 5]], 'none', 95);
  assert.ok(near(perfect.kappa, 1, 1e-12));
});

test('fleissKappa: statsmodels fleiss_kappa 와 같고 (평가자 2명이면 Scott π), 범주별 kappa 의 가중평균 항등식이 성립한다', () => {
  assert.equal(ref.fleiss.length, 3);
  ref.fleiss.forEach(c => {
    const X = c.R.map(row => { const x = new Array(c.cats.length).fill(0); row.forEach(v => x[c.cats.indexOf(v)]++); return x; });
    const r = A.fleissKappa(X);
    const m = c.name;
    assert.ok(r.N === c.n && r.n === c.raters && near(r.pbar, c.pbar, 1e-12) && near(r.pe, c.pe, 1e-12) && near(r.kappa, c.kappa, 1e-10), m + ' kappa');
    assert.ok(near(r.se, c.se, 1e-10) && near(r.z, c.z, 1e-9) && near(r.pval, c.pval, 1e-7, 1e-14), m + ' 표준오차·z');
    r.kappaJ.forEach((v, j) => assert.ok(near(v, c.kappa_j[j], 1e-10) && near(r.zJ[j], c.z_j[j], 1e-9) && near(r.pJ[j], c.p_j[j], 1e-7, 1e-14), `${m} 범주 ${j}`));
    const num = r.p.reduce((s, p, j) => s + p * (1 - p) * r.kappaJ[j], 0), den = r.p.reduce((s, p) => s + p * (1 - p), 0);
    assert.ok(near(num / den, r.kappa, 1e-9));
  });
});

test('iccCompute: pingouin·numpy가 서로 일치한 기준값과 같다 — ICC 6유형의 값, F, 자유도, 유의확률, 신뢰구간 (6가지 자료, 90% 포함)', () => {
  assert.equal(ref.icc.length, 6);
  ref.icc.forEach(c => {
    const r = A.iccCompute(colsOf(c.Y), c.conf);
    const m = c.name;
    assert.ok(r.n === c.n && r.k === c.k && near(r.ssb, c.ssb, 1e-10) && near(r.ssj, c.ssj, 1e-10) && near(r.sse, c.sse, 1e-9) && near(r.msw, c.msw, 1e-10), m + ' 제곱합');
    r.rows.forEach((x, q) => {
      const e = c.rows[q];
      assert.equal(x.key, e.key);
      assert.ok(near(x.icc, e.icc, 1e-9, 1e-12) && near(x.F, e.F, 1e-9) && x.df1 === e.df1 && near(x.df2, e.df2, 1e-12) && near(x.p, e.p, 1e-7, 1e-14), `${m} ${x.key}`);
      assert.ok(near(x.lo, e.lo, 1e-6, 1e-9) && near(x.hi, e.hi, 1e-6, 1e-9), `${m} ${x.key} 신뢰구간 ${x.lo},${x.hi} vs ${e.lo},${e.hi}`);
    });
  });
  // 변동이 전혀 없으면 null
  const flat = A.iccCompute([[5, 5, 5, 5], [5, 5, 5, 5]], 95);
  assert.ok(flat.rows.every(x => x.icc === null));
});

test('runKappa(평가자 2명): 표 구성(교차표, kappa), 가중 옵션, 범주 레이블, 경고(kappa 역설), 렌더·내보내기', () => {
  const c = ref.kappa[0];                                  // 이진, 한쪽 범주가 많음
  const ds = dsOfCols([c.a, c.b], ['평가A', '평가B'], 'nominal', { 0: '아니오', 1: '예' });
  const res = A.runKappa(ds, { vars: ['평가A', '평가B'], options: {} });
  assert.equal(res.title, '평가자 일치도 (Cohen kappa)');
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption.split(' (')[0]), ['처리 요약', '평가 일치 교차표', 'Cohen kappa']);
  const cross = table(res, '평가 일치 교차표');
  assert.deepEqual(cross.rows.map(r => r.name), ['아니오', '예', '합계']);
  assert.deepEqual([cross.rows[0].c0, cross.rows[0].c1, cross.rows[1].c0, cross.rows[1].c1], [c.table[0][0], c.table[0][1], c.table[1][0], c.table[1][1]]);
  assert.equal(cross.rows[2].sum, c.none.n);
  const kp = table(res, 'Cohen kappa').rows[0];
  assert.ok(near(kp.kappa, c.none.kappa, 1e-9) && near(kp.se, c.none.se, 1e-9) && near(kp.lo, c.none.lo, 1e-9) && near(kp.z, c.none.z, 1e-9) && near(kp.po, c.none.po, 1e-12) && near(kp.kmax, c.none.kmax, 1e-9));
  assert.ok(table(res, 'Cohen kappa').notes[0].includes('Landis & Koch'));
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('평가 일치 교차표') && html.includes('아니오'));
  assert.ok(Ex.resultToGrid(res).rows.length > 10);
  // 가중 kappa (5점 서열): 선형·이차 옵션이 표 제목과 값에 반영된다
  const o = ref.kappa[1];
  const dso = dsOfCols([o.a, o.b], ['A', 'B'], 'ordinal');
  ['linear', 'quadratic'].forEach(w => {
    const r = A.runKappa(dso, { vars: ['A', 'B'], options: { weights: w, conf: '95' } });
    const t = table(r, w === 'linear' ? '선형 가중 kappa' : '이차 가중 kappa');
    assert.ok(t, w);
    assert.ok(near(t.rows[0].kappa, o[w].kappa, 1e-9) && near(t.rows[0].z, o[w].z, 1e-9));
    assert.ok(t.columns.every(x => x.key !== 'kmax'));
  });
  // 일치율은 높은데 kappa 가 낮은 경우의 경고 (한쪽 범주에 쏠림)
  const a = [], b = [];
  for (let i = 0; i < 100; i++) { a.push(i < 90 ? 1 : (i < 95 ? 0 : 1)); b.push(i < 90 ? 1 : (i < 95 ? 1 : 0)); }
  const paradox = A.runKappa(dsOfCols([a, b], ['A', 'B'], 'nominal'), { vars: ['A', 'B'], options: {} });
  assert.ok(paradox.warnings.some(w => w.includes('kappa 역설')));
});

test('runKappa(평가자 3명 이상): Fleiss kappa, 범주별 kappa, 쌍별 kappa 행렬(Light kappa), 가중 옵션은 쌍별에만', () => {
  const c = ref.fleiss[0];                                 // 3명 × 3범주
  const cols = colsOf(c.R);
  const names = ['A', 'B', 'C'];
  const ds = dsOfCols(cols, names, 'nominal');
  const res = A.runKappa(ds, { vars: names, options: {} });
  assert.equal(res.title, '평가자 일치도 (Fleiss kappa)');
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption.split(' (')[0]), ['처리 요약', 'Fleiss kappa', '평가자 쌍별 Cohen kappa']);
  const ft = table(res, 'Fleiss kappa').rows;
  assert.equal(ft.length, 4);
  assert.ok(near(ft[3].kappa, c.kappa, 1e-10) && near(ft[3].se, c.se, 1e-10) && near(ft[3].z, c.z, 1e-9) && near(ft[3].pv, c.pval, 1e-7, 1e-14));
  ft.slice(0, 3).forEach((r, j) => assert.ok(near(r.kappa, c.kappa_j[j], 1e-10) && r.cnt === c.counts[j] && near(r.p, c.p[j], 1e-12) && near(r.z, c.z_j[j], 1e-9)));
  const pr = table(res, '평가자 쌍별');
  pr.rows.forEach((row, i) => ['c0', 'c1', 'c2'].forEach((key, j) => assert.ok(near(row[key], c.pair[i][j], 1e-10), `쌍 ${i},${j}`)));
  assert.ok(pr.notes[0].includes('Light kappa') && pr.notes[0].includes(c.light.toFixed(3)));
  assert.ok(table(res, 'Fleiss kappa').notes[1].includes('평가자가 2명일 때만'));
  // 4명
  const c4 = ref.fleiss[1];
  const r4 = A.runKappa(dsOfCols(colsOf(c4.R), ['A', 'B', 'C', 'D'], 'nominal'), { vars: ['A', 'B', 'C', 'D'], options: { weights: 'quadratic' } });
  assert.ok(near(table(r4, 'Fleiss kappa').rows[4].kappa, c4.kappa, 1e-10));
  assert.ok(table(r4, '평가자 쌍별').caption.includes('이차 가중'));
});

test('runKappa: 결측 제외, 범주 레이블 우선순위, 입력 검증', () => {
  const c = ref.kappa[2];
  const ds = dsOfCols([c.a, c.b], ['A', 'B'], 'nominal', { 1: '하', 2: '중', 3: '상' }, rows => { rows[4][0] = ''; rows[9][1] = ''; });
  const res = A.runKappa(ds, { vars: ['A', 'B'], options: {} });
  assert.ok(res.warnings.some(w => w.includes('2건을 제외') && w.includes('48건')));
  assert.deepEqual(table(res, '평가 일치 교차표').rows.map(r => r.name), ['하', '중', '상', '합계']);
  // 평가자 B 가 쓰지 않은 범주(주변합 0)도 교차표에 남고 kappa 가 계산된다 (기준값 3번)
  const ds2 = dsOfCols([c.a, c.b], ['A', 'B'], 'nominal');
  const r2 = A.runKappa(ds2, { vars: ['A', 'B'], options: {} });
  assert.ok(near(table(r2, 'Cohen kappa').rows[0].kappa, c.none.kappa, 1e-9));
  const run = (vars, d) => A.runKappa(d || ds2, { vars: vars, options: {} });
  assert.throws(() => run(['A']), /2개 이상/);
  assert.throws(() => A.runKappa(ds2, { vars: ['A', 'B'], options: { conf: 100 } }), /신뢰수준/);
  const one = dsOfCols([[1, 1, 1], [1, 1, 1]], ['A', 'B'], 'nominal');
  assert.throws(() => run(['A', 'B'], one), /범주가 1개/);
  const many = dsOfCols([Array.from({ length: 25 }, (_, i) => i), Array.from({ length: 25 }, (_, i) => i)], ['A', 'B'], 'scale');
  assert.throws(() => run(['A', 'B'], many), /범주가 25개.*최대 20개/);
  const few = dsOfCols([[1], [2]], ['A', 'B'], 'nominal');
  assert.throws(() => run(['A', 'B'], few), /2건 이상/);
  const degenerate = dsOfCols([[1, 1, 1, 1], [1, 1, 1, 1]], ['A', 'B'], 'nominal', {}, null);
  assert.throws(() => run(['A', 'B'], degenerate), /범주가 1개/);
});

test('runIcc: 표 구성(기술통계, 분산분석표, ICC 6유형), 굵은 행, 신뢰수준, 결측 제외, 경고, 입력 검증', () => {
  const c = ref.icc[2];                                    // 3명 평가자 × 12대상, 편향 큼
  const names = ['평가1', '평가2', '평가3'];
  const ds = dsOfCols(colsOf(c.Y), names, 'scale');
  const res = A.runIcc(ds, { vars: names, options: {} });
  assert.equal(res.title, '급내상관계수 (ICC)');
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption.split(' (')[0]), ['처리 요약', '기술통계', '이원 분산분석표', '급내상관계수']);
  const an = table(res, '이원 분산분석표').rows;
  assert.ok(near(an[0].ss, c.ssb, 1e-10) && near(an[1].ss, c.ssj, 1e-10) && near(an[2].ss, c.sse, 1e-9) && an[3].df === 12 * 3 - 1);
  const it = table(res, '급내상관계수').rows;
  assert.deepEqual(it.map(r => r.label), ['ICC(1,1)', 'ICC(2,1)', 'ICC(3,1)', 'ICC(1,k)', 'ICC(2,k)', 'ICC(3,k)']);
  assert.deepEqual(it.map(r => !!r._bold), [false, true, true, false, false, false]);
  it.forEach((r, q) => assert.ok(near(r.icc, c.rows[q].icc, 1e-9, 1e-12) && near(r.lo, c.rows[q].lo, 1e-6, 1e-9) && near(r.hi, c.rows[q].hi, 1e-6, 1e-9) && near(r.p, c.rows[q].p, 1e-7, 1e-14)));
  assert.ok(c.rows[1].icc < c.rows[2].icc);                // 평가자 편향이 크면 절대 일치도 < 일관성
  assert.ok(it[0].desc.includes('일원 확률') && it[1].desc.includes('절대 일치도') && it[2].desc.includes('이원 혼합'));
  assert.ok(table(res, '급내상관계수').notes.some(n => n.includes('Koo & Li')));
  assert.ok(res.warnings.some(w => w.includes('대상이 12건')));
  assert.ok(R.renderResultHtml(res).includes('ICC(2,1)') && Ex.resultToGrid(res).rows.length > 12);
  // 90% 신뢰수준 (기준값에 90% 사례가 있다)
  const c90 = ref.icc[5];
  const r90 = A.runIcc(dsOfCols(colsOf(c90.Y), ['a', 'b', 'c'], 'scale'), { vars: ['a', 'b', 'c'], options: { conf: '90' } });
  const t90 = table(r90, '급내상관계수');
  assert.ok(t90.columns.some(x => x.label === '90% 신뢰구간 하한'));
  assert.ok(near(t90.rows[1].lo, c90.rows[1].lo, 1e-6, 1e-9) && near(t90.rows[1].hi, c90.rows[1].hi, 1e-6, 1e-9));
  // 결측 대상은 통째로 제외
  const dm = dsOfCols(colsOf(c.Y), names, 'scale', {}, rows => { rows[3][1] = ''; });
  const rm = A.runIcc(dm, { vars: names, options: {} });
  assert.ok(rm.warnings.some(w => w.includes('1건을 제외') && w.includes('11건')));
  assert.deepEqual(table(rm, '처리 요약').rows.slice(0, 2).map(r => r.n), [11, 1]);
  // 입력 검증
  assert.throws(() => A.runIcc(ds, { vars: ['평가1'], options: {} }), /2개 이상/);
  assert.throws(() => A.runIcc(ds, { vars: names, options: { conf: 0 } }), /신뢰수준/);
  const two = dsOfCols([[1, 2], [2, 3]], ['a', 'b'], 'scale');
  assert.throws(() => A.runIcc(two, { vars: ['a', 'b'], options: {} }), /3건 이상/);
  const flat = dsOfCols([[5, 5, 5, 5], [5, 5, 5, 5]], ['a', 'b'], 'scale');
  assert.ok(A.runIcc(flat, { vars: ['a', 'b'], options: {} }).warnings.some(w => w.includes('계산할 수 없습니다')));
  const txt = Data.buildDataset([['a', 'b'], [1, 'x'], [2, 'y'], [3, 'z'], [4, 'w']], [scaleDef('b', 'ordinal')]);
  assert.throws(() => A.runIcc(txt, { vars: ['a', 'b'], options: {} }), /숫자가 아닌 값/);
});

test('runIcc: 옵션 anova 를 끄면 이원 분산분석표가 없다', () => {
  const c = ref.icc[0];
  const names = ['a', 'b', 'c', 'd'];
  const ds = dsOfCols(colsOf(c.Y), names, 'scale');
  assert.ok(table(A.runIcc(ds, { vars: names, options: {} }), '이원 분산분석표'));
  assert.ok(!table(A.runIcc(ds, { vars: names, options: { anova: false } }), '이원 분산분석표'));
});
