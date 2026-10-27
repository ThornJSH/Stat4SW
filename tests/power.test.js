const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const P = require('../src/Stats_Power.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'power.json'), 'utf8'));
const near = (a, e, tol, abs = 1e-12) => Math.abs(a - e) <= Math.max(abs, tol * Math.abs(e));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const popt = r => ({ alpha: r.alpha, tails: r.tails, ratio: r.ratio || 1, k: r.k || 3, u: r.u || 3, minN: r.test === 'reg' ? (r.u || 3) + 2 : undefined });
const val = (res, caption, item) => table(res, caption).rows.find(x => x.item.startsWith(item)).value;

test('pncf: scipy ncf.cdf 와 같다 (자유도·비중심도 8가지 × 분위 4개, 큰 비중심도 포함)', () => {
  assert.equal(ref.ncf.length, 32);
  ref.ncf.forEach(g => {
    const got = P.pncf(g.x, g.d1, g.d2, g.lam);
    assert.ok(near(got, g.cdf, 1e-8, 1e-12), `ncf x=${g.x} d=(${g.d1},${g.d2}) λ=${g.lam}: ${got} vs ${g.cdf}`);
  });
  assert.equal(P.pncf(0, 3, 20, 5), 0);
});

test('powerOf: statsmodels·pingouin·scipy 가 서로 일치한 기준값과 같다 — 6가지 검정의 검정력 (격자 수백 개)', () => {
  const rows = ref.rows.filter(r => r.power !== undefined);
  assert.ok(rows.length > 300);
  const byTest = {};
  rows.forEach(r => {
    const got = P.powerOf(r.test, r.n, r.es, popt(r));
    assert.ok(near(got, r.power, 1e-6, 1e-9), `${r.test} es=${r.es} n=${r.n} α=${r.alpha} tails=${r.tails} ratio=${r.ratio}: ${got} vs ${r.power}`);
    byTest[r.test] = (byTest[r.test] || 0) + 1;
  });
  assert.deepEqual(Object.keys(byTest).sort(), ['anova', 'corr', 'prop', 'reg', 'ttest1', 'ttest2']);
});

test('solveN: 목표 검정력 이상이 되는 가장 작은 정수와 같고 (n 에서는 이상, n−1 에서는 미만), 탐지 가능한 효과크기가 statsmodels 와 같다', () => {
  const nreq = ref.rows.filter(r => r.kind === 'nreq');
  assert.ok(nreq.length > 80);
  nreq.forEach(r => {
    const got = P.solveN(r.test, r.es, r.target, popt(r));
    assert.equal(got, r.nreq, `${r.test} es=${r.es} α=${r.alpha} tails=${r.tails} ratio=${r.ratio} k=${r.k} u=${r.u}`);
    assert.ok(P.powerOf(r.test, got, r.es, popt(r)) >= r.target);
    if (got > popt(r).minN) assert.ok(P.powerOf(r.test, got - 1, r.es, popt(r)) < r.target);
  });
  const eff = ref.rows.filter(r => r.kind === 'effect');
  assert.ok(eff.length >= 10);
  eff.forEach(r => {
    const got = P.solveEffect(r.test, r.n, r.target, popt(r));
    assert.ok(near(got, r.effect, 1e-5, 1e-8), `${r.test} 효과크기 n=${r.n} α=${r.alpha}: ${got} vs ${r.effect}`);
  });
  // 알려진 값 (G*Power 와 같은 관행): d=.5 → 집단당 64, f=.25(k=3) → 집단당 53, f²=.15(u=3) → 총 77
  assert.equal(P.solveN('ttest2', 0.5, 0.8, { alpha: 0.05, tails: 2, ratio: 1 }), 64);
  assert.equal(P.solveN('anova', 0.25, 0.8, { alpha: 0.05, tails: 2, ratio: 1, k: 3 }), 53);
  assert.equal(P.solveN('reg', 0.15, 0.8, { alpha: 0.05, tails: 2, ratio: 1, u: 3, minN: 5 }), 77);
  // 한도 안에서 달성할 수 없는 경우
  assert.equal(P.solveEffect('ttest1', 2, 0.9999, { alpha: 0.001, tails: 2, ratio: 1 }), null);
});

test('runPower(표본 크기): 표 구성, 총 인원, 실제 검정력, 효과크기 관행별 표, 곡선, 탈락 감안 모집 인원, 렌더·내보내기', () => {
  const res = P.runPower(null, { options: { test: 'ttest2', solve: 'n', effect: '0.5', alpha: '0.05', power: '0.8', dropout: '15' } });
  assert.equal(res.title, '검정력·표본 크기 계산 (독립표본 t검정)');
  assert.deepEqual(res.blocks.map(b => b.caption.split(' (')[0]), ['입력 조건', '결과', '효과크기별 필요 표본 크기', '표본 크기에 따른 검정력']);
  assert.equal(val(res, '결과', '필요한 표본 크기'), '64');
  assert.equal(val(res, '결과', '집단 2의 N'), '64');
  assert.equal(val(res, '결과', '총 인원'), '128');
  assert.ok(near(Number(val(res, '결과', '그 크기에서 실제 검정력')), ref.rows.find(r => r.test === 'ttest2' && r.es === 0.5 && r.n === 64 && r.ratio === 1 && r.alpha === 0.05 && r.tails === 2).power, 1e-3));
  assert.equal(val(res, '결과', '탈락을 감안한 모집 인원'), '76');                // ⌈64 / 0.85⌉
  assert.equal(val(res, '결과', '탈락을 감안한 총 모집 인원'), '152');
  const conv = table(res, '효과크기별 필요 표본 크기').rows;
  assert.deepEqual(conv.map(r => [r.lab, r.es]), [['작음', 0.2], ['중간', 0.5], ['큼', 0.8]]);
  assert.deepEqual(conv.map(r => r.val.split(' ')[0]), ['394', '64', '26']);       // statsmodels 의 같은 조건 표본 크기
  const ch = res.blocks.find(b => b.type === 'chart');
  assert.deepEqual(ch.series.map(s => s.label), ['검정력', '목표 검정력']);
  assert.ok(ch.series[0].points.length > 10 && ch.series[0].points.every((p, i, a) => i === 0 || (p[1] >= a[i - 1][1] - 1e-12)));
  assert.ok(R.renderResultHtml(res).includes('필요한 표본 크기') && R.renderResultHtml(res).includes('<svg'));
  assert.ok(Ex.resultToGrid(res).rows.length > 12);
  // 효과크기를 비우면 중간 효과, 집단 크기 비 2:1 이면 집단 2 가 두 배
  const med = P.runPower(null, { options: { test: 'anova', solve: 'n', k: '4' } });
  assert.ok(table(med, '입력 조건').rows.find(r => r.item.startsWith('효과크기')).value.includes('중간 효과 사용'));
  assert.equal(val(med, '결과', '총 인원'), String(4 * Number(val(med, '결과', '필요한 표본 크기'))));
  const r2 = P.runPower(null, { options: { test: 'ttest2', solve: 'n', effect: '0.5', ratio: '2' } });
  assert.equal(Number(val(r2, '결과', '집단 2의 N')), 2 * Number(val(r2, '결과', '필요한 표본 크기')));
});

test('runPower(검정력·효과크기 구하기, 두 비율, 단측, 회귀)', () => {
  const pw = P.runPower(null, { options: { test: 'ttest1', solve: 'power', effect: '0.4', n: '30' } });
  const e = ref.rows.find(r => r.test === 'ttest1' && r.es === 0.5 && r.n === 34);
  assert.ok(Number(val(pw, '결과', '검정력')) > 0.5 && Number(val(pw, '결과', '검정력')) < 0.9);
  assert.ok(table(pw, '효과크기별 검정력').rows.length === 3);
  const pw2 = P.runPower(null, { options: { test: 'ttest1', solve: 'power', effect: '0.5', n: '34', alpha: String(e.alpha), tails: 'two' } });
  assert.ok(near(Number(val(pw2, '결과', '검정력')), P.powerOf('ttest1', 34, 0.5, { alpha: e.alpha, tails: 2, ratio: 1 }), 1e-4));
  const low = P.runPower(null, { options: { test: 'ttest1', solve: 'power', effect: '0.2', n: '20' } });
  assert.ok(low.warnings.some(w => w.includes('(0.80)에 못 미칩니다')));
  const ef = P.runPower(null, { options: { test: 'ttest2', solve: 'effect', n: '40' } });
  const refEff = ref.rows.find(r => r.test === 'ttest2' && r.kind === 'effect' && r.alpha === 0.05 && r.tails === 2);
  assert.ok(near(Number(val(ef, '결과', '탐지 가능한 최소 효과크기')), refEff.effect, 1e-4));
  assert.ok(val(ef, '결과', '크기 해석').includes('중간'));
  // 두 비율: h 계산과 표본 크기
  const hp = ref.rows.find(r => r.kind === 'h' && r.p1 === 0.5);
  const pr = P.runPower(null, { options: { test: 'prop', solve: 'n', p1: '0.5', p2: '0.3' } });
  assert.ok(table(pr, '입력 조건').rows.find(r => r.item.startsWith('효과크기')).value.includes(hp.h.toFixed(3)));
  assert.equal(Number(val(pr, '결과', '필요한 표본 크기')), P.solveN('prop', hp.h, 0.8, { alpha: 0.05, tails: 2, ratio: 1 }));
  const one = P.runPower(null, { options: { test: 'corr', solve: 'n', effect: '0.3', tails: 'one' } });
  assert.equal(Number(val(one, '결과', '필요한 표본 크기')), ref.rows.find(r => r.test === 'corr' && r.kind === 'nreq' && r.es === 0.3 && r.alpha === 0.05 && r.tails === 1).nreq);
  assert.ok(table(one, '결과').notes.some(n => n.includes('Fisher z')));
  const rg = P.runPower(null, { options: { test: 'reg', solve: 'n', effect: '0.15', u: '3' } });
  assert.equal(val(rg, '결과', '필요한 표본 크기'), '77');
  assert.ok(table(rg, '입력 조건').rows.some(r => r.item === '독립변수 수' && r.value === '3'));
});

test('runPower: 입력 검증', () => {
  const run = o => P.runPower(null, { options: o });
  assert.throws(() => run({ test: 'ttest2', solve: 'power' }), /"표본 크기"를 입력/);
  assert.throws(() => run({ test: 'ttest2', solve: 'effect' }), /"표본 크기"를 입력/);
  assert.throws(() => run({ alpha: '0' }), /유의수준/);
  assert.throws(() => run({ alpha: '1' }), /유의수준/);
  assert.throws(() => run({ power: '1' }), /목표 검정력/);
  assert.throws(() => run({ effect: 'abc' }), /효과크기.*숫자/);
  assert.throws(() => run({ effect: '0' }), /효과크기/);
  assert.throws(() => run({ test: 'corr', effect: '1' }), /1보다 작아야/);
  assert.throws(() => run({ test: 'prop' }), /두 비율/);
  assert.throws(() => run({ test: 'prop', p1: '0.4', p2: '0.4' }), /두 비율이 같으면/);
  assert.throws(() => run({ test: 'anova', k: '1' }), /집단 수/);
  assert.throws(() => run({ test: 'anova', k: '2.5' }), /정수/);
  assert.throws(() => run({ test: 'ttest2', solve: 'power', n: '1.5' }), /정수/);
  assert.throws(() => run({ test: 'ttest2', solve: 'power', n: '1' }), /2 이상/);
  assert.throws(() => run({ test: 'reg', solve: 'power', n: '4', u: '5', effect: '0.15' }), /독립변수 수 \+ 2/);
  assert.throws(() => run({ dropout: '100' }), /탈락률/);
  assert.throws(() => run({ ratio: '0' }), /집단 크기 비/);
  assert.throws(() => run({ test: 'ttest2', effect: '0.0001' }), /1000000명을 넘습니다/);
  assert.throws(() => run({ test: 'ttest1', solve: 'effect', n: '2', power: '0.9999', alpha: '0.001' }), /달성할 수 있는 효과크기가 없습니다/);
  // 아주 작은 효과의 큰 표본 경고
  const big = run({ test: 'corr', effect: '0.05' });
  assert.ok(big.warnings.some(w => w.includes('매우 큽니다')));
});
