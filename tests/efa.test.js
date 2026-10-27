const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const E = require('../src/Stats_Efa.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'efa.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-12) => Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const colsOf = c => c.data[0].map((_, j) => c.data.map(r => r[j]));
const nameOf = j => 'q' + (j + 1);

function dsOf(c, extra) {
  const names = c.data[0].map((_, j) => nameOf(j));
  const rows = [names].concat(c.data.map(r => r.slice()));
  if (extra) extra(rows);
  return { names, ds: Data.buildDataset(rows, names.map(n => ({ name: n, label: '', scale: 'scale', valueLabels: {}, missing: [] }))) };
}
const maxDiff = (A, B) => Math.max(...A.flat().map((v, i) => Math.abs(v - B.flat()[i])));

test('efaAnalyze: numpy·statsmodels·factor_analyzer가 서로 일치한 기준값과 같다 — KMO, Bartlett, 고유값, 적재값, 공통성, 회전 (7가지 사례)', () => {
  assert.equal(ref.cases.length, 7);
  ref.cases.forEach(c => {
    const r = E.efaAnalyze(colsOf(c), { method: c.method, rotation: c.rotation, nfactors: c.nfactors });
    const w = c.name;
    assert.equal(r.m, c.m, w + ' 요인 수');
    assert.equal(r.kaiser, c.kaiser);
    r.eig.forEach((v, k) => assert.ok(near(v, c.eig[k], 1e-9, 1e-12), `${w} 고유값 ${k}`));
    assert.ok(near(r.kmo, c.kmo, 1e-10) && near(r.bartlett.chi2, c.bartlett.chi2, 1e-9) && r.bartlett.df === c.bartlett.df && near(r.bartlett.p, c.bartlett.p, 1e-7, 1e-14), w + ' KMO·Bartlett');
    r.msa.forEach((v, k) => assert.ok(near(v, c.msa[k], 1e-10), `${w} MSA ${k}`));
    const iter = c.method === 'paf' ? 1e-5 : 1e-9;               // 주축분해는 수렴 기준 1e-6 에서 멈춘다
    assert.ok(maxDiff(r.L, c.L) < iter, `${w} 적재값 ${maxDiff(r.L, c.L)}`);
    r.h2.forEach((v, k) => assert.ok(Math.abs(v - c.h2[k]) < iter, `${w} 공통성 ${k}`));
    r.h2Initial.forEach((v, k) => assert.ok(Math.abs(v - c.h2_init[k]) < 1e-9, `${w} 초기 공통성 ${k}`));
    assert.ok(maxDiff([r.ssExt], [c.ss_ext]) < iter * 10, w + ' 추출 제곱합');
    if (c.rotated) {
      assert.ok(maxDiff(r.rotated, c.rotated) < (c.method === 'paf' ? 2e-5 : 1e-6), `${w} 회전 ${maxDiff(r.rotated, c.rotated)}`);
      assert.ok(maxDiff([r.ssRot], [c.ss_rot]) < 1e-4, w + ' 회전 제곱합');
      // 직각 회전은 공통성을 바꾸지 않는다
      r.rotated.forEach((row, k) => assert.ok(near(row.reduce((s, v) => s + v * v, 0), r.h2[k], 1e-9, 1e-9), `${w} 회전 뒤 공통성 ${k}`));
    } else {
      assert.equal(r.rotated, null, w);
    }
  });
});

test('runEfa: 표 구성(KMO·Bartlett, 공통성, 설명된 총분산, 요인 행렬, 회전된 요인 행렬), 스크리 그림, 정렬과 숨김, 렌더·내보내기', () => {
  const c = ref.cases[0];                                 // 300명 9변수 3요인 PAF + 베리맥스
  const { names, ds } = dsOf(c);
  const res = E.runEfa(ds, { vars: names, options: {} });
  assert.equal(res.title, '탐색적 요인분석');
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption), ['KMO와 Bartlett 검정', '공통성', '설명된 총분산', '요인 행렬 (회전 전)', '회전된 요인 행렬']);
  const k = table(res, 'KMO').rows[0];
  assert.ok(near(k.kmo, c.kmo, 1e-10) && near(k.chi2, c.bartlett.chi2, 1e-9) && k.df === 36);
  assert.ok(table(res, 'KMO').notes[0].includes('주축분해') && table(res, 'KMO').notes[0].includes('베리맥스') && table(res, 'KMO').notes[0].includes('고유값 1 이상'));
  const comm = table(res, '공통성').rows;
  assert.equal(comm.length, 9);
  assert.ok(near(comm[0].ext, c.h2[0], 1e-5) && near(comm[0].init, c.h2_init[0], 1e-9) && near(comm[0].msa, c.msa[0], 1e-10));
  const tv = table(res, '설명된 총분산').rows;
  assert.equal(tv.length, 9);
  assert.ok(near(tv[0].iTot, c.eig[0], 1e-9) && near(tv[0].iPct, 100 * c.eig[0] / 9, 1e-9) && near(tv[8].iCum, 100, 1e-9));
  assert.ok(tv[3].eTot === undefined && near(tv[2].eCum, 100 * c.ss_ext.reduce((s, v) => s + v, 0) / 9, 1e-5) && tv[2].rTot !== undefined);
  const rot = table(res, '회전된 요인 행렬');
  assert.deepEqual(rot.columns.map(x => x.label), ['', '요인 1', '요인 2', '요인 3']);
  assert.equal(rot.rows.length, 9);
  // 크기순 정렬: 가장 크게 적재된 요인 번호가 오름차순, 숨김: 0.3 미만은 비움
  const main = rot.rows.map(r => [1, 2, 3].map(f => r['f' + (f - 1)]).findIndex(v => v !== null && v >= 0.5));
  assert.deepEqual(main, [0, 0, 0, 1, 1, 1, 2, 2, 2]);
  assert.ok(rot.rows.every(r => [0, 1, 2].filter(f => r['f' + f] !== null).length >= 1));
  assert.ok(rot.rows.some(r => [0, 1, 2].some(f => r['f' + f] === null)));
  const ch = res.blocks.find(b => b.type === 'chart');
  assert.equal(ch.kind, 'line');
  assert.equal(ch.series[0].points.length, 9);
  assert.deepEqual(ch.series[1].points, [[1, 1], [9, 1]]);
  assert.ok(R.renderResultHtml(res).includes('Bartlett') && R.renderResultHtml(res).includes('회전된 요인 행렬'));
  assert.ok(Ex.resultToGrid(res).rows.length > 30);
  // 옵션: 스크리 끄기, 정렬 끄기, 숨김 0, 회전 없음, PCA, 요인 수 지정
  const r2 = E.runEfa(ds, { vars: names, options: { scree: false, sort: false, suppress: '0', rotation: 'none', method: 'pca', nfactors: '2' } });
  assert.ok(!r2.blocks.some(b => b.type === 'chart'));
  assert.deepEqual(r2.blocks.filter(b => b.type === 'table').map(b => b.caption), ['KMO와 Bartlett 검정', '공통성', '설명된 총분산', '성분 행렬']);
  const lt = table(r2, '성분 행렬');
  assert.deepEqual(lt.columns.map(x => x.label), ['', '성분 1', '성분 2']);
  assert.ok(lt.rows.every(r => r.f0 !== null && r.f1 !== null));
  assert.equal(lt.rows[0].name, nameOf(0));                 // 정렬 끄면 원래 순서
  assert.ok(table(r2, '설명된 총분산').rows[0].rTot === undefined);
  assert.ok(table(r2, 'KMO').notes[0].includes('직접 지정'));
});

test('runEfa: 요인 1개는 회전하지 않고 안내, 주성분의 공통성, 반복 횟수 표시', () => {
  const c = ref.cases[5];                                 // 1요인
  const { names, ds } = dsOf(c);
  const res = E.runEfa(ds, { vars: names, options: {} });
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption), ['KMO와 Bartlett 검정', '공통성', '설명된 총분산', '요인 행렬']);
  assert.ok(res.warnings.some(w => w.includes('요인이 1개라 회전하지 않았습니다')));
  assert.ok(table(res, '요인 행렬').notes[0].includes('반복 후 수렴'));
  const pc = ref.cases[2];                                // 주성분
  const d2 = dsOf(pc);
  const r2 = E.runEfa(d2.ds, { vars: d2.names, options: { method: 'pca' } });
  assert.ok(table(r2, '공통성').rows.every(r => r.init === 1));
  assert.ok(table(r2, '공통성').notes[0].includes('초기 공통성이 1'));
  assert.ok(table(r2, '설명된 총분산').caption === '설명된 총분산' && table(r2, '성분 행렬 (회전 전)'));
});

test('runEfa: 경고(교차적재, 낮은 공통성, KMO, 사례 수), 결측 사례 제외', () => {
  // 교차적재: 상관 2요인 + 교차적재 변수 포함 자료에서 요인 수 2
  const b = ref.cases[1];
  const { names, ds } = dsOf(b);
  const res = E.runEfa(ds, { vars: names, options: {} });
  assert.ok(res.warnings.every(w => !w.includes('사례 수')));
  // 무관한 잡음 변수를 하나 더 → 낮은 공통성과 사례 수 경고 (소표본)
  const small = { data: b.data.slice(0, 60).map((r, i) => r.concat([Math.sin(i * 12.9898) * 43758.5453 % 1])) };
  const d3 = dsOf(small);
  const r3 = E.runEfa(d3.ds, { vars: d3.names, options: { nfactors: '2' } });
  assert.ok(r3.warnings.some(w => w.includes('사례 수(60건)가 적습니다')));
  assert.ok(r3.warnings.some(w => w.includes('공통성이') && w.includes('q9')));
  // 결측 사례는 통째로 제외
  const d4 = dsOf(b, rows => { rows[3][1] = ''; rows[9][5] = ''; });
  const r4 = E.runEfa(d4.ds, { vars: d4.names, options: {} });
  assert.ok(r4.warnings.some(w => w.includes('2건을 제외') && w.includes('198건')));
  assert.ok(table(r4, 'KMO').notes[0].includes('사례 198건'));
});

test('runEfa: 입력 검증과 계산 불가 경우', () => {
  const c = ref.cases[3];
  const { names, ds } = dsOf(c);
  const run = (vars, options) => E.runEfa(ds, { vars: vars, options: options || {} });
  assert.throws(() => run(names.slice(0, 2)), /3개 이상/);
  assert.throws(() => run(names, { nfactors: '0' }), /요인 수/);
  assert.throws(() => run(names, { nfactors: '2.5' }), /요인 수/);
  assert.throws(() => run(names, { nfactors: '6' }), /1 이상 5 이하/);               // 주축분해는 변수 수 − 1 까지
  assert.doesNotThrow(() => run(names, { nfactors: '6', method: 'pca' }));
  assert.throws(() => run(names, { suppress: '1' }), /작은 적재값/);
  // 사례 수가 변수 수 이하
  const few = Data.buildDataset([['a', 'b', 'c'], [1, 2, 3], [2, 1, 4], [3, 5, 1]], []);
  assert.throws(() => E.runEfa(few, { vars: ['a', 'b', 'c'], options: {} }), /변수 수 3개보다 많아야/);
  // 값이 모두 같은 변수
  const rows = [['a', 'b', 'c', 'd']];
  for (let i = 0; i < 12; i++) rows.push([Math.sin(i) * 3 + i, Math.cos(i * 2) * 2, i % 4, 7]);
  assert.throws(() => E.runEfa(Data.buildDataset(rows, []), { vars: ['a', 'b', 'c', 'd'], options: {} }), /"d"의 값이 모두 같아/);
  // 같은 정보의 변수(완전 공선성) → 상관행렬이 특이
  const dup = [['a', 'b', 'c', 'd']];
  for (let i = 0; i < 12; i++) dup.push([Math.sin(i) * 3 + i, Math.cos(i * 2) * 2, i % 4, 2 * (Math.sin(i) * 3 + i) + 1]);
  assert.throws(() => E.runEfa(Data.buildDataset(dup, []), { vars: ['a', 'b', 'c', 'd'], options: {} }), /특이/);
  // 서로 직교인 변수: 고유값이 모두 1 → 고유값 1 이상 요인 없음
  const orth = [['a', 'b', 'c']];
  [[1, 1, 1], [1, -1, 1], [-1, 1, 1], [-1, -1, 1], [1, 1, -1], [1, -1, -1], [-1, 1, -1], [-1, -1, -1]].forEach(r => orth.push(r));
  assert.throws(() => E.runEfa(Data.buildDataset(orth, []), { vars: ['a', 'b', 'c'], options: {} }), /고유값이 1 이상인 요인이 없습니다/);
  // 숫자가 아닌 값
  const txt = Data.buildDataset([['a', 'b', 'c'], [1, 2, 'x'], [2, 3, 'y'], [3, 1, 'z'], [4, 5, 'w']], [{ name: 'c', label: '', scale: 'ordinal', valueLabels: {}, missing: [] }]);
  assert.throws(() => E.runEfa(txt, { vars: ['a', 'b', 'c'], options: {} }), /숫자가 아닌 값/);
});

test('efaVarimax·efaSortOrder: 직교 회전은 공통성을 보존하고 돌려 놓은 단순 구조를 되찾는다', () => {
  const base = [[0.8, 0.1], [0.7, -0.05], [0.75, 0], [0.05, 0.8], [0, 0.7], [-0.1, 0.75]];
  const th = 25 * Math.PI / 180, c = Math.cos(th), s = Math.sin(th);
  const mixed = base.map(([a, b]) => [a * c - b * s, a * s + b * c]);
  const crit = L => [0, 1].reduce((sum, f) => {
    const col = L.map((row, i) => row[f] * row[f] / L[i].reduce((q, v) => q + v * v, 0));   // Kaiser 정규화한 제곱 적재값
    const mean = col.reduce((q, v) => q + v, 0) / col.length;
    return sum + col.reduce((q, v) => q + (v - mean) * (v - mean), 0);
  }, 0);
  const r = E.efaCanonical(E.efaVarimax(mixed).L);
  assert.ok(crit(r) > crit(mixed));
  r.forEach((row, i) => {
    assert.ok(Math.abs(row[0] * row[0] + row[1] * row[1] - (base[i][0] ** 2 + base[i][1] ** 2)) < 1e-12, '공통성 보존 ' + i);
    const big = Math.abs(row[0]) > Math.abs(row[1]) ? 0 : 1;
    assert.ok(Math.abs(row[big]) > 0.65 && Math.abs(row[1 - big]) < 0.15, `행 ${i}: ${row}`);
  });
  // 크기순: 가장 크게 적재된 요인 번호 오름차순, 같은 요인 안에서는 적재값 큰 순
  assert.deepEqual(E.efaSortOrder([[0.1, 0.9], [0.8, 0.1], [0.7, 0.2], [0.2, 0.6]]), [1, 2, 0, 3]);
});

/* ---------- 2단계: 프로맥스, 병렬분석, 요인 점수 ---------- */

test('프로맥스: factor_analyzer(Kaiser 정규화)와 같은 정의로 numpy 가 만든 패턴·구조 행렬, 요인 상관과 같다 (6가지 사례)', () => {
  const withPromax = ref.cases.filter(c => c.promax);
  assert.equal(withPromax.length, 6);
  withPromax.forEach(c => {
    const w = c.name;
    const pm = E.efaPromax(c.L, 4);
    assert.ok(maxDiff(pm.P, c.promax.P) < 1e-6 && maxDiff(pm.phi, c.promax.phi) < 1e-6, `${w} 패턴·요인 상관`);
    assert.ok(maxDiff(E.efaStructure(pm.P, pm.phi), c.promax.S) < 1e-6, w + ' 구조');
    pm.phi.forEach((row, a) => { assert.ok(near(row[a], 1, 1e-9)); row.forEach((v, b) => assert.ok(near(v, pm.phi[b][a], 1e-9, 1e-12))); });
    // 사각 회전도 공통성은 diag(P Φ P') 로 보존된다
    pm.P.forEach((row, i) => {
      let s = 0;
      row.forEach((vi, a) => row.forEach((vj, b) => { s += vi * pm.phi[a][b] * vj; }));
      assert.ok(near(s, c.h2[i], 1e-5, 1e-6), `${w} 공통성 ${i}`);
    });
    const r = E.efaAnalyze(colsOf(c), { method: c.method, rotation: 'promax', nfactors: c.nfactors });
    assert.equal(r.oblique, true);
    assert.ok(maxDiff(r.rotated, c.promax.P) < (c.method === 'paf' ? 1e-4 : 1e-6), w + ' efaAnalyze 패턴');
    assert.equal(r.ssRot, undefined);
  });
  // 요인 1개는 회전하지 않는다
  const one = ref.cases[5];
  const r1 = E.efaAnalyze(colsOf(one), { method: 'paf', rotation: 'promax', nfactors: null });
  assert.equal(r1.rotated, null);
  assert.equal(r1.oblique, false);
});

test('runEfa(프로맥스): 패턴·구조·요인 상관 표, 회전 뒤 제곱합 열 없음, 안내', () => {
  const c = ref.cases[1];
  const { names, ds } = dsOf(c);
  const res = E.runEfa(ds, { vars: names, options: { rotation: 'promax', nfactors: '2' } });
  assert.deepEqual(res.blocks.filter(b => b.type === 'table').map(b => b.caption), ['KMO와 Bartlett 검정', '공통성', '설명된 총분산', '요인 행렬 (회전 전)', '패턴 행렬', '구조 행렬', '요인 상관 행렬']);
  assert.ok(table(res, 'KMO').notes[0].includes('프로맥스 (사각, Kaiser 정규화, κ=4)'));
  assert.ok(!table(res, '설명된 총분산').columns.some(x => x.key === 'rTot'));
  assert.ok(table(res, '설명된 총분산').notes.some(n => n.includes('사각 회전')));
  const phi = table(res, '요인 상관 행렬').rows;
  assert.ok(near(phi[0].f0, 1, 1e-9) && near(phi[1].f1, 1, 1e-9) && near(phi[0].f1, c.promax.phi[0][1], 1e-4));
  assert.ok(phi[0].f1 > 0.2 && phi[0].f1 < 0.8);                  // 상관 .4 로 만든 자료
  const pat = table(res, '패턴 행렬').rows, str = table(res, '구조 행렬').rows;
  assert.equal(pat.length, 8);
  const row = k => pat.find(r => r.name === nameOf(k));
  assert.ok(near(row(0).f0, c.promax.P[0][0], 1e-3));
  const srow = str.find(r => r.name === nameOf(2));
  assert.ok(near(srow.f0, c.promax.S[2][0], 1e-3) && near(srow.f1, c.promax.S[2][1], 1e-3));
  // 요인 1개이면 회전하지 않는다는 안내
  const d1 = dsOf(ref.cases[5]);
  const r1 = E.runEfa(d1.ds, { vars: d1.names, options: { rotation: 'promax' } });
  assert.ok(r1.warnings.some(w => w.includes('회전하지 않았습니다')));
  assert.ok(!r1.blocks.some(b => b.caption === '패턴 행렬'));
});

test('병렬분석: numpy 모의(4000회)의 고유값 평균·95백분위와 허용오차 안에서 같다 (난수 기반), 합 = 변수 수, 같은 시드는 같은 결과', () => {
  assert.equal(ref.parallel.length, 6);
  ref.parallel.forEach(rp => {
    const sims = E.efaParallelSims(rp.n, rp.p);
    const got = E.efaParallel(rp.n, rp.p, 2026, sims);
    assert.equal(got.sims, sims);
    got.mean.forEach((v, k) => assert.ok(Math.abs(v - rp.mean[k]) < 0.03, `n=${rp.n} p=${rp.p} 평균 ${k}: ${v} vs ${rp.mean[k]}`));
    got.p95.forEach((v, k) => assert.ok(Math.abs(v - rp.p95[k]) < 0.06, `n=${rp.n} p=${rp.p} 95백분위 ${k}: ${v} vs ${rp.p95[k]}`));
    assert.ok(near(got.mean.reduce((s, v) => s + v, 0), rp.p, 1e-9));
    got.mean.forEach((v, k) => { if (k) assert.ok(v <= got.mean[k - 1] + 1e-12); });
    got.p95.forEach((v, k) => assert.ok(v >= got.mean[k] - 1e-9));
  });
  const a = E.efaParallel(100, 6, 7, 20), b = E.efaParallel(100, 6, 7, 20), c = E.efaParallel(100, 6, 8, 20);
  assert.deepEqual(a, b);
  assert.notDeepEqual(a.mean, c.mean);
  assert.equal(E.efaParallelSims(300, 9), 100);
  assert.equal(E.efaParallelSims(5000, 50), 20);
});

test('runEfa(병렬분석): 표·그림·권장 요인 수, 시드 검증', () => {
  const c = ref.cases[0];                                  // 직교 3요인
  const { names, ds } = dsOf(c);
  const res = E.runEfa(ds, { vars: names, options: { parallel: true, seed: '11' } });
  const pt = table(res, '병렬분석');
  assert.ok(pt.rows.length >= 4 && pt.rows.length <= 9);
  assert.ok(near(pt.rows[0].obs, c.eig[0], 1e-9));
  assert.deepEqual(pt.rows.slice(0, 4).map(r => r.judge), ['유지', '유지', '유지', '제외']);
  assert.ok(pt.notes[1].includes('권하는 요인 수: 3개') && pt.notes[0].includes('시드 11'));
  const ch = res.blocks.find(b => b.type === 'chart');
  assert.deepEqual(ch.series.map(s => s.label), ['고유값', '기준 1', '난수 자료의 95백분위']);
  assert.equal(ch.series[2].points.length, 9);
  assert.ok(!res.warnings.some(w => w.includes('병렬분석은')));
  // 옵션을 끄면 표가 없다
  assert.ok(!table(E.runEfa(ds, { vars: names, options: {} }), '병렬분석'));
  // 요인이 없는 잡음 자료: 병렬분석이 권하는 수는 고유값 기준보다 적거나 같다
  let seed = 5;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  const rows = [['a', 'b', 'c', 'd', 'e', 'f']];
  for (let i = 0; i < 60; i++) rows.push(Array.from({ length: 6 }, () => Math.round((rnd() + rnd() + rnd() - 1.5) * 1000) / 1000));
  const noise = Data.buildDataset(rows, []);
  const rn = E.runEfa(noise, { vars: rows[0], options: { parallel: true } });
  const recN = Number(/권하는 요인 수: (\d+)개/.exec(table(rn, '병렬분석').notes[1])[1]);
  assert.ok(recN <= 1);
  assert.throws(() => E.runEfa(ds, { vars: names, options: { parallel: true, seed: '1.5' } }), /시드/);
});

test('runEfa(요인 점수 저장): 회귀 방법 점수가 기준값과 같고 결측 사례는 빈 값, 옵션을 끄면 저장하지 않는다', () => {
  const c = ref.cases[0];                                  // PAF + 베리맥스, 요인 3개
  const { names, ds } = dsOf(c);
  const res = E.runEfa(ds, { vars: names, options: { save: true } });
  assert.equal(res.saveColumns.length, 3);
  assert.deepEqual(res.saveColumns.map(x => x.name), ['요인점수1', '요인점수2', '요인점수3']);
  assert.ok(res.saveColumns[0].label.includes('주축분해·베리맥스') && res.saveColumns[0].scale === '연속' && res.saveColumns[0].fixed === 4);
  res.saveColumns.forEach((col, f) => {
    assert.equal(col.values.length, 300);
    for (let i = 0; i < 6; i++) assert.ok(near(col.values[i], c.scores.first[i][f], 1e-3, 1e-5), `점수 ${i},${f}`);
    const mean = col.values.reduce((s, v) => s + v, 0) / 300;
    const v = col.values.reduce((s, x) => s + (x - mean) * (x - mean), 0) / 299;
    assert.ok(Math.abs(mean) < 1e-9 && near(v, c.scores.var[f], 1e-3), `평균·분산 ${f}`);
  });
  assert.ok(res.warnings.some(w => w.includes('회귀 방법(Thurstone)')));
  assert.equal(E.runEfa(ds, { vars: names, options: {} }).saveColumns, undefined);
  // 프로맥스: 구조 행렬 기준 점수
  const rp = E.runEfa(ds, { vars: names, options: { save: true, rotation: 'promax' } });
  for (let i = 0; i < 6; i++) rp.saveColumns.forEach((col, f) => assert.ok(near(col.values[i], c.scores_promax.first[i][f], 1e-3, 1e-5), `프로맥스 점수 ${i},${f}`));
  assert.ok(rp.saveColumns[0].label.includes('프로맥스'));
  // 주성분은 성분점수, 결측 사례는 null (rows[5] = 사례 인덱스 4)
  const d2 = dsOf(ref.cases[2], rows => { rows[5][1] = ''; });
  const r2 = E.runEfa(d2.ds, { vars: d2.names, options: { save: true, method: 'pca' } });
  assert.equal(r2.saveColumns[0].name, '성분점수1');
  assert.equal(r2.saveColumns[0].values[4], null);
  assert.ok(r2.saveColumns[0].values[0] !== null && r2.saveColumns[0].values.length === 150);
  // 주성분 점수(회귀 방법)는 표준화 성분 점수와 같아 분산이 1이다
  const vals = r2.saveColumns[1].values.filter(v => v !== null);
  const mean = vals.reduce((s, v) => s + v, 0) / vals.length;
  assert.ok(Math.abs(vals.reduce((s, x) => s + (x - mean) ** 2, 0) / (vals.length - 1) - 1) < 1e-9);
});
