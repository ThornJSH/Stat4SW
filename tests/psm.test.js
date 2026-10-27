const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const P = require('../src/Stats_Psm.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'psm.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-10, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));

const D = ref.data, N = D.y.length;
function makeDs(extra) {
  const rows = [['t', 'x1', 'x2', 'g', 'y']];
  for (let i = 0; i < N; i++) rows.push([D.treat[i], D.x1[i], D.x2[i], D.g[i], D.y[i]]);
  if (extra) extra(rows);
  return Data.buildDataset(rows, [
    { name: 't', label: '처치', scale: 'nominal', valueLabels: { '0': '미참여', '1': '참여' }, missing: [] },
    { name: 'g', label: '구분', scale: 'nominal', valueLabels: {}, missing: [] }
  ]);
}
const dsFull = makeDs();
const req = (c, extra) => Object.assign({
  treatment: ['t'], covariates: ['x1', 'x2'], factors: ['g'], outcome: ['y'],
  options: { method: c.options.method, caliper: String(c.options.caliper), ratio: String(c.options.ratio), replace: c.options.replace }
}, extra);

test('psmMatch: 파이썬 직접 구현과 같은 짝(세트), 가중치, ID, 캘리퍼 폭 — 4가지 설정', () => {
  const eta = ref.fit.eta, treat = D.treat;
  const G1 = D.g.map(v => (v === 1 ? 1 : 0)), G2 = D.g.map(v => (v === 2 ? 1 : 0));
  const cols = [D.x1, D.x2, G1, G2];
  assert.equal(ref.cases.length, 4);
  for (const c of ref.cases) {
    const o = c.options;
    const M = P.psmMatch(eta, treat, cols, { caliper: o.caliper, ratio: o.ratio, replace: o.replace, method: o.method });
    assert.deepEqual(M.matched, c.matched, c.name);
    Object.keys(c.sets).forEach(t => assert.deepEqual(M.sets[Number(t)], c.sets[t], `${c.name} 처치 ${t}`));
    assert.equal(M.matched.length, Object.keys(c.sets).length, c.name);
    M.weights.forEach((w, i) => assert.ok(near(w, c.weights[i], 1e-12) || Math.abs(w - c.weights[i]) < 1e-12, `${c.name} 가중치 ${i}`));
    M.matched.forEach(t => assert.equal(M.ids[t], c.ids[String(t)], `${c.name} ID ${t}`));
    if (c.width === null) assert.equal(M.width, null); else assert.ok(near(M.width, c.width, 1e-9), c.name + ' 폭');
  }
});

test('runPsm: 성향점수 모형, 매칭 요약, 균형(SMD 전·후), 효과(ATT, Welch)가 기준값과 같다', () => {
  for (const c of ref.cases) {
    const res = P.runPsm(dsFull, req(c));
    const e = c.effect;
    const sum = table(res, '매칭 요약').rows;
    assert.equal(sum[1].n, c.matched.length, c.name);
    assert.equal(sum[2].n, c.n_treated - c.matched.length, c.name);
    assert.equal(sum[4].n, c.n_used_controls, c.name);
    const bal = table(res, '균형 진단').rows;
    assert.equal(bal.length, 5);
    c.smd.forEach((s, k) => {
      for (const key of ['mt0', 'mc0', 'smd0', 'mt1', 'mc1', 'smd1']) assert.ok(near(bal[k][key], s[key], 1e-8), `${c.name} ${s.name}.${key}: ${bal[k][key]} vs ${s[key]}`);
    });
    const eff = table(res, '처치 효과 추정').rows;
    assert.ok(near(eff[1].diff, e.att, 1e-9) && near(eff[1].se, e.se, 1e-8) && near(eff[1].t, e.t, 1e-8) && eff[1].df === e.df, c.name + ' ATT');
    assert.ok(near(eff[1].p, e.p, 1e-7) && near(eff[1].lo, e.lo, 1e-8) && near(eff[1].hi, e.hi, 1e-8), c.name + ' ATT 추론');
    assert.ok(near(eff[1].mt, e.m_treat, 1e-9) && near(eff[1].mc, e.m_ctrl, 1e-9), c.name + ' 평균');
    assert.ok(near(eff[0].diff, e.naive_diff, 1e-9) && near(eff[0].t, e.naive_t, 1e-8) && near(eff[0].df, e.naive_df, 1e-8) && near(eff[0].p, e.naive_p, 1e-7), c.name + ' Welch');
  }
  // 성향점수 모형은 파이썬의 계수와 같고 공통 지지 영역도 같다
  const res = P.runPsm(dsFull, req(ref.cases[0]));
  const coef = table(res, '성향점수 모형: 방정식의 변수').rows;
  ref.fit.B.forEach((b, k) => assert.ok(near(coef[k].B, b, 1e-8), `B[${k}]`));
  const sup = table(res, '공통 지지 영역').rows;
  assert.ok(near(sup[0].min, ref.fit.eta_min_t, 1e-9) && near(sup[0].max, ref.fit.eta_max_t, 1e-9) && near(sup[1].min, ref.fit.eta_min_c, 1e-9) && near(sup[1].max, ref.fit.eta_max_c, 1e-9));
  assert.ok(table(res, '처치 변수 코딩').rows.some(r => r.orig === '참여' && r.code === 1));
  assert.ok(res.blocks.some(b => b.type === 'chart' && b.kind === 'box' && b.groups.length === 4));
  assert.ok(R.renderResultHtml(res).includes('SMD'));
  assert.ok(Ex.resultToGrid(res).rows.length > 30);
});

test('runPsm: 저장 열(성향점수, 매칭ID, 매칭가중치)은 사례 순서에 맞고 결측 사례는 빈 값, 결과 변수 결측은 목록별로 제외', () => {
  const c = ref.cases[0];
  const ds = makeDs(rows => { rows[6][2] = ''; rows[10][4] = ''; });                // 공변량 결측(사례 5), 결과 결측(사례 9)
  const res = P.runPsm(ds, req(c, { options: { method: 'ps', caliper: '0.2', ratio: '1', replace: false, save: true } }));
  const [prob, id, w] = res.saveColumns;
  assert.deepEqual([prob.name, id.name, w.name], ['성향점수', '매칭ID', '매칭가중치']);
  assert.equal(prob.values.length, N);
  assert.equal(prob.values[5], null);
  assert.equal(prob.values[9], null);
  assert.ok(prob.values.filter(v => v !== null).every(v => v > 0 && v < 1));
  assert.equal(w.values.filter(v => v !== null && v > 0).length >= 2 * 40, true);   // 처치와 대조가 같은 수만큼 짝지어짐 (1:1)
  const matchedIds = id.values.filter(v => v !== null);
  assert.equal(new Set(matchedIds).size * 2, matchedIds.length);                       // 1:1 비복원: 세트마다 처치 1 + 대조 1
  assert.ok(res.warnings.some(x => x.includes('2건을 제외')));
  // 저장하지 않으면 열이 없다
  assert.equal(P.runPsm(dsFull, req(c)).saveColumns, undefined);
});

test('runPsm: 입력 검증과 경고', () => {
  const c = ref.cases[0];
  const r = (extra, ds) => P.runPsm(ds || dsFull, Object.assign(req(c), extra));
  assert.throws(() => r({ covariates: [], factors: [] }), /공변량/);
  assert.throws(() => r({ treatment: ['g'] }), /정확히 2가지|동시에 쓸 수 없습니다/);
  assert.throws(() => r({ outcome: ['x1'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => r({ options: { caliper: '-1' } }), /캘리퍼/);
  assert.throws(() => r({ options: { ratio: '0' } }), /1~5/);
  assert.throws(() => r({ options: { ratio: '1.5' } }), /정수/);
  assert.throws(() => r({ options: { conf: '100' } }), /신뢰수준/);
  assert.throws(() => r({ options: { caliper: '0.0001' } }), /매칭된 사례가 없습니다/);
  // 결과 변수 없이도 실행되고 효과 표는 없다
  const noOut = P.runPsm(dsFull, { treatment: ['t'], covariates: ['x1'], factors: ['g'], options: { method: 'maha', caliper: '0.3' } });
  assert.ok(!table(noOut, '처치 효과') && table(noOut, '균형 진단'));
  // 못 짝지은 처치는 경고, 캘리퍼가 좁으면 매칭 뒤 균형 경고가 나올 수 있다
  const tight = r({ options: { caliper: '0.05', ratio: '1' } });
  assert.ok(tight.warnings.some(w => w.includes('짝을 찾지 못해 버렸습니다')));
  // 처치 사례가 적으면 경고
  const few = makeDs(rows => { let k = 0; for (let i = 1; i < rows.length; i++) { if (rows[i][0] === 1 && ++k > 8) rows[i][0] = 0; } });
  assert.ok(P.runPsm(few, { treatment: ['t'], covariates: ['x1', 'x2'], options: { caliper: '0' } }).warnings.some(w => w.includes('매우 적습니다')));
  // 작은 값을 처치로 바꾸면 처치 수가 달라진다
  const flip = r({ options: { method: 'ps', caliper: '0', flip: true } });
  assert.equal(table(flip, '처리 요약').rows[2].n, c.n_control);
});

test('psmBalance·psmEffect·psmWelch: 손으로 계산한 값', () => {
  // 처치 평균 6, 대조(전) 평균 4, 처치 SD(n−1) = 2 → SMD 1. 매칭 후 대조 가중 평균 5.5 → SMD 0.25
  const x = [4, 6, 8, 3, 4, 5, 6, 7];                   // 처치 0~2(4, 6, 8), 대조 3~7(3, 4, 5, 6, 7)
  const treat = [1, 1, 1, 0, 0, 0, 0, 0];
  const w = [1, 1, 1, 0, 0, 0.5, 1, 1.5];                // 매칭된 대조 가중 합 3
  const b = P.psmBalance([x], ['x'], treat, w, [0, 1, 2])[0];
  assert.equal(b.mt0, 6); assert.equal(b.mc0, 5); assert.ok(Math.abs(b.smd0 - 0.5) < 1e-12);
  assert.ok(Math.abs(b.mc1 - (0.5 * 5 + 1 * 6 + 1.5 * 7) / 3) < 1e-12);
  assert.ok(Math.abs(b.smd1 - (6 - 19 / 3) / 2) < 1e-12);
  // 효과: 세트 0 → 대조 [2, 4] 평균 3, 세트 1 → 대조 [5] 평균 5. y처치 = 10, 9 → 차이 7, 4 → ATT 5.5, SD 2.1213, SE 1.5
  const y = [10, 9, 0, 2, 4, 5];
  const e = P.psmEffect(y, [[3, 4], [5]], [0, 1], 95);
  assert.ok(Math.abs(e.att - 5.5) < 1e-12 && Math.abs(e.se - 1.5) < 1e-12 && e.df === 1);
  assert.ok(Math.abs(e.t - 5.5 / 1.5) < 1e-12);
  assert.ok(Math.abs(e.lo - (5.5 - 12.706204736432095 * 1.5)) < 1e-6);
  // Welch: a = [1, 2, 3], b = [2, 4, 6] → diff −2, se sqrt(1/3 + 4/3) = 1.2910, df 2.5, t −1.5492
  const wl = P.psmWelch([1, 2, 3], [2, 4, 6], 95);
  assert.ok(Math.abs(wl.diff + 2) < 1e-12 && Math.abs(wl.se - Math.sqrt(5 / 3)) < 1e-12 && Math.abs(wl.t + 2 / Math.sqrt(5 / 3)) < 1e-12);
  assert.ok(Math.abs(wl.df - (5 / 3) ** 2 / ((1 / 3) ** 2 / 2 + (4 / 3) ** 2 / 2)) < 1e-12);
});
