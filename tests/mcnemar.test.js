const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const M = require('../src/Stats_Mcnemar.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'mcnemar.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-10, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));

function dsOf(c, labels) {
  const rows = [['전', '후']];
  c.before.forEach((v, i) => rows.push([v, c.after[i]]));
  return Data.buildDataset(rows, labels ? [{ name: '전', label: '사전 취업', scale: 'nominal', valueLabels: labels, missing: [] }, { name: '후', label: '사후 취업', scale: 'nominal', valueLabels: labels, missing: [] }] : []);
}
const run = (c, extra, labels) => M.runMcnemar(dsOf(c, labels), Object.assign({ first: ['전'], second: ['후'], options: {} }, extra));

test('mcnemarTest: scipy 기준값(정확 p, 카이제곱 p, 비율 차이와 구간, 오즈비의 정확 구간) 5가지', () => {
  assert.equal(ref.cases.length, 5);
  ref.cases.forEach(c => {
    const t = M.mcnemarTest(c.a, c.b, c.c, c.d, 95);
    assert.ok(near(t.exactP, c.exact, 1e-9), `${c.name} 정확 p ${t.exactP} vs ${c.exact}`);
    assert.ok(near(t.chi2, c.chi2, 1e-12) && near(t.chiP, c.chi_p, 1e-8), c.name + ' 카이제곱');
    assert.ok(near(t.diff, c.diff, 1e-12) && near(t.diffLo, c.diff_lo, 1e-8) && near(t.diffHi, c.diff_hi, 1e-8), c.name + ' 비율 차이');
    assert.ok(near(t.odds, c.odds, 1e-12), c.name + ' 오즈비');
    assert.ok(near(t.oddsLo, c.odds_lo, 1e-7) && near(t.oddsHi, c.odds_hi, 1e-7), `${c.name} 오즈비 구간 ${t.oddsLo}..${t.oddsHi} vs ${c.odds_lo}..${c.odds_hi}`);
  });
  assert.equal(M.mcnemarExactP(0, 0), 1);
  assert.equal(M.mcnemarExactP(5, 5), 1);
});

test('runMcnemar: 행 단위 자료에서 만든 표와 검정 표가 기준값과 같다 (값이 1/2 코딩이어도, 값 레이블이 있으면 레이블로)', () => {
  ref.cases.forEach(c => {
    const res = run(c);
    assert.equal(res.title, 'McNemar 검정');
    const x = table(res, '교차표').rows;
    assert.deepEqual([x[0].c0, x[0].c1, x[1].c0, x[1].c1], [c.a, c.b, c.c, c.d], c.name);
    assert.equal(x[2].total, c.n);
    const s = table(res, 'McNemar 검정 통계량').rows[0];
    assert.equal(s.nd, c.nd);
    assert.ok(near(s.exactP, c.exact, 1e-9) && near(s.chi2, c.chi2, 1e-12) && near(s.chiP, c.chi_p, 1e-8), c.name);
    const ch = table(res, '변화 요약').rows;
    assert.ok(near(ch[0].value, c.p_before, 1e-12) && near(ch[1].value, c.p_after, 1e-12) && near(ch[2].value, c.diff, 1e-12) && near(ch[2].lo, c.diff_lo, 1e-8), c.name);
    assert.ok(near(ch[3].value, c.odds, 1e-12) && near(ch[3].lo, c.odds_lo, 1e-7) && near(ch[3].hi, c.odds_hi, 1e-7), c.name + ' 오즈비');
  });
  const lab = run(ref.cases[0], {}, { '0': '미취업', '1': '취업' });
  assert.deepEqual(table(lab, '교차표').columns.map(c => c.label), ['전 \\ 후', '미취업', '취업', '합계']);
  assert.ok(table(lab, '교차표').caption.includes('사전 취업') && table(lab, '교차표').notes[0].includes('미취업 → 후 취업'));
  const r = R.renderResultHtml(run(ref.cases[0]));
  assert.ok(r.includes('McNemar') && r.includes('정확 유의확률'));
  assert.ok(Ex.resultToGrid(run(ref.cases[0])).rows.length > 10);
});

test('runMcnemar: 결측 제외, 경고(불일치 쌍 적음·없음), 입력 검증', () => {
  const rows = [['전', '후'], [0, 1], [1, 0], ['', 1], [0, ''], [1, 1], [0, 0], [0, 0]];
  const ds = Data.buildDataset(rows, []);
  const res = M.runMcnemar(ds, { first: ['전'], second: ['후'], options: {} });
  assert.deepEqual(table(res, '처리 요약').rows.map(r => r.n), [5, 2, 7]);
  assert.ok(res.warnings.some(w => w.includes('2건을 제외')) && res.warnings.some(w => w.includes('2건뿐이라 검정력')));
  const same = Data.buildDataset([['전', '후'], [0, 0], [1, 1], [0, 0], [1, 1]], []);
  const nochange = M.runMcnemar(same, { first: ['전'], second: ['후'], options: {} });
  assert.ok(nochange.warnings.some(w => w.includes('하나도 없어')));
  assert.equal(table(nochange, 'McNemar 검정 통계량').rows[0].exactP, 1);
  const three = Data.buildDataset([['전', '후'], [0, 1], [1, 2], [2, 0]], []);
  assert.throws(() => M.runMcnemar(three, { first: ['전'], second: ['후'], options: {} }), /정확히 2가지[\s\S]*3가지/);
  const one = Data.buildDataset([['전', '후'], [1, 1], [1, 1], [1, 1]], []);
  assert.throws(() => M.runMcnemar(one, { first: ['전'], second: ['후'], options: {} }), /1가지/);
  assert.throws(() => M.runMcnemar(ds, { first: ['전'], second: ['전'], options: {} }), /서로 다른/);
  assert.throws(() => M.runMcnemar(ds, { first: [], second: ['후'], options: {} }), /변수를 선택|하나씩/);
  assert.throws(() => M.runMcnemar(Data.buildDataset([['전', '후'], [0, 1]], []), { first: ['전'], second: ['후'], options: {} }), /2건 이상/);
});
