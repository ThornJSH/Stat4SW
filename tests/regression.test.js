const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const Reg = require('../src/Stats_Regression.js');
const R = require('../src/Render.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'regression.json'), 'utf8'));

function relErr(a, e) {
  if (a === e) return 0;
  return Math.abs(a - e) / Math.max(Math.abs(e), 1e-300);
}
function close(a, e, tol, msg, absTol = 1e-12) {
  assert.ok(relErr(a, e) < tol || Math.abs(a - e) < absTol, `${msg}: ${a} vs ${e}`);
}
const ds = (values, defs) => Data.buildDataset(values, defs || []);
const tbl = (res, caption) => res.blocks.find(b => b.caption && b.caption.startsWith(caption));

/** 기준값 사례에서 데이터셋과 요청을 만든다 */
function build(c, options) {
  const names = ['y'].concat(Object.keys(c.covariates));
  const n = c.y.length;
  const header = names.slice();
  if (c.factor) header.push(c.factor.name);
  const rows = [header];
  for (let i = 0; i < n; i++) {
    const row = [c.y[i]].concat(Object.keys(c.covariates).map(k => c.covariates[k][i]));
    if (c.factor) row.push(c.factor.values[i]);
    rows.push(row);
  }
  const request = { dependent: ['y'], covariates: Object.keys(c.covariates), options: Object.assign({ conf: c.conf || 95 }, options) };
  if (c.factor) request.factors = [c.factor.name];
  if (c.lastRef) request.options.lastRef = true;
  return { dataset: ds(rows), request };
}

/* ---------- 계수와 분산분석: 정확 산술 기준 ---------- */

test('runRegression: 정확 산술 기준값 대조 (계수, 표준오차, t, p, 신뢰구간, 베타, 공차, VIF)', () => {
  for (const c of ref.cases) {
    const { dataset, request } = build(c, { ci: true, dw: true });
    const res = Reg.runRegression(dataset, request);
    const e = c.expected;
    const coef = tbl(res, '계수').rows;
    assert.deepEqual(coef.map(r => r.name), e.names, c.name);
    coef.forEach((r, i) => {
      const m = `${c.name} ${e.names[i]}`;
      close(r.B, e.B[i], 1e-9, m + ' B', 1e-10);
      close(r.se, e.se[i], 1e-9, m + ' se');
      close(r.t, e.t[i], 1e-8, m + ' t');
      close(r.p, e.p[i], 1e-7, m + ' p');
      close(r.lo, e.lo[i], 1e-8, m + ' lo', 1e-9);
      close(r.hi, e.hi[i], 1e-8, m + ' hi', 1e-9);
      if (e.beta[i] === null) {
        assert.equal(r.beta, null); assert.equal(r.vif, null); assert.equal(r.tol, null);
      } else {
        close(r.beta, e.beta[i], 1e-9, m + ' beta');
        close(r.vif, e.vif[i], 1e-8, m + ' VIF');
        close(r.tol, e.tol[i], 1e-8, m + ' tol');
      }
    });
  }
});

test('runRegression: 모형 요약과 분산분석 (R², 수정 R², F, 제곱합, 더빈-왓슨)', () => {
  for (const c of ref.cases) {
    const { dataset, request } = build(c, { dw: true });
    const res = Reg.runRegression(dataset, request);
    const e = c.expected;
    const sm = tbl(res, '모형 요약').rows[0];
    close(sm.r, e.r, 1e-9, c.name + ' R');
    close(sm.r2, e.r2, 1e-9, c.name + ' R2');
    close(sm.adj, e.adjR2, 1e-9, c.name + ' adj');
    close(sm.s, e.s, 1e-9, c.name + ' s');
    close(sm.dw, e.dw, 1e-9, c.name + ' DW');
    const a = tbl(res, '분산분석').rows;
    assert.deepEqual(a.map(r => r.df), [e.dfm, e.dfr, e.dfm + e.dfr]);
    close(a[0].ss, e.ssr, 1e-8, c.name + ' SSR');
    close(a[1].ss, e.sse, 1e-8, c.name + ' SSE');
    close(a[2].ss, e.sst, 1e-9, c.name + ' SST');
    close(a[0].F, e.F, 1e-8, c.name + ' F');
    close(a[0].p, e.pF, 1e-7, c.name + ' p(F)');
  }
});

test('runRegression: 단순회귀 R²는 상관계수의 제곱 (scipy linregress와 교차 확인된 기준)', () => {
  const c = ref.cases[0];
  const e = c.expected;
  close(e.r, Math.sqrt(e.r2), 1e-12, 'r');
  const res = Reg.runRegression(...Object.values(build(c, {})));
  close(tbl(res, '모형 요약').rows[0].r2, e.r2, 1e-9, 'R2');
});

/* ---------- 공선성 진단: numpy eigh 기준 ---------- */

test('collinearityDiagnostics: numpy 고유분해 기준값 대조 (고유값, 조건지수, 분산 비율)', () => {
  for (const c of ref.cases.filter(x => x.diag)) {
    const names = Object.keys(c.covariates);
    const preds = names.map(k => ({ name: k, values: c.covariates[k] }));
    const got = Reg.collinearityDiagnostics(preds, c.y.length);
    got.eig.forEach((v, j) => {
      close(v, c.diag.eig[j], 1e-9, `${c.name} 고유값 ${j}`, 1e-12);
      close(got.cond[j], c.diag.cond[j], 1e-8, `${c.name} 조건지수 ${j}`);
      got.prop[j].forEach((p, i) => close(p, c.diag.prop[j][i], 1e-7, `${c.name} 분산 비율 ${j},${i}`, 1e-9));
    });
    // 계수마다 분산 비율의 합은 1
    for (let i = 0; i < names.length + 1; i++) {
      assert.ok(Math.abs(got.prop.reduce((s, row) => s + row[i], 0) - 1) < 1e-12);
    }
  }
});

test('runRegression: 공선성 진단 옵션이 표를 추가하고 강한 공선성은 VIF 경고', () => {
  const c = ref.cases.find(x => x.name === '강한 공선성');
  const { dataset, request } = build(c, { collinearity: true });
  const res = Reg.runRegression(dataset, request);
  const t = tbl(res, '공선성 진단');
  assert.equal(t.rows.length, 4);                         // 차원 = 상수 + 독립변수 3
  assert.equal(t.columns.length, 3 + 4);                  // 차원, 고유값, 조건지수 + 계수 4개의 분산 비율
  assert.ok(t.rows[3].cond > 30);                         // 마지막 차원의 조건지수가 큼
  assert.ok(res.warnings.some(w => w.includes('VIF')));
  const none = Reg.runRegression(build(c, {}).dataset, build(c, {}).request);
  assert.equal(tbl(none, '공선성 진단'), undefined);
});

/* ---------- 선형대수 단위 ---------- */

test('jacobiEigen_: 알려진 고유값과 직교성', () => {
  const e = Reg.jacobiEigen_([[2, 1, 0], [1, 2, 0], [0, 0, 5]]);
  close(e.values[0], 5, 1e-12, 'λ1'); close(e.values[1], 3, 1e-12, 'λ2'); close(e.values[2], 1, 1e-12, 'λ3');
  for (let a = 0; a < 3; a++) {
    for (let b = 0; b < 3; b++) {
      const dot = e.vectors[a].reduce((s, v, i) => s + v * e.vectors[b][i], 0);
      assert.ok(Math.abs(dot - (a === b ? 1 : 0)) < 1e-12);
    }
  }
});

test('qrLeastSquares_: 정확히 맞는 선형 관계와 공선성 감지', () => {
  const x1 = [1, 2, 3, 4, 5, 6], x2 = [2, 1, 4, 3, 6, 5];
  const y = x1.map((v, i) => 3 * v - 2 * x2[i]);
  const fit = Reg.qrLeastSquares_([x1, x2], y);
  assert.ok(fit.ok);
  close(fit.beta[0], 3, 1e-12, 'b1'); close(fit.beta[1], -2, 1e-12, 'b2');
  const bad = Reg.qrLeastSquares_([x1, x1.map(v => 2 * v)], y);
  assert.equal(bad.ok, false);
  assert.equal(bad.badCol, 1);
});

/* ---------- 실행 함수: 구조와 입력 오류 ---------- */

const defs = Data.parseVarDefs([['h'], ['g', '집단', '명목', '1=대조, 2=처치A, 3=처치B', '']]).defs;

function sample() {
  const rows = [['y', 'x', 'z', 'g']];
  for (let i = 0; i < 20; i++) rows.push([5 + 2 * i + (i % 3) * 0.7, i, (i * 7) % 11, 1 + (i % 3)]);
  rows.push(['', 3, 4, 1], [10, '', 4, 2], [12, 5, 6, '']);   // 결측 3건
  return ds(rows, defs);
}

test('runRegression: 표 구성(모형 요약, 분산분석, 계수), 더미 이름과 기준 범주, 결측 제외', () => {
  const res = Reg.runRegression(sample(), { dependent: ['y'], covariates: ['x'], factors: ['g'], options: {} });
  assert.deepEqual(res.blocks.map(b => b.caption.split(' ')[0]), ['모형', '분산분석', '계수']);
  const c = tbl(res, '계수');
  assert.deepEqual(c.rows.map(r => r.name), ['(상수)', 'x', '집단 = 처치A', '집단 = 처치B']);
  assert.ok(c.notes.some(n => n.includes('기준 범주: 대조')));
  assert.ok(tbl(res, '모형 요약').notes[0].includes('3건 제외'));
  assert.equal(tbl(res, '분산분석').rows[0].df, 3);
  assert.equal(tbl(res, '분산분석').rows[1].df, 16);        // 유효 20건 − 계수 4개
  assert.ok(!c.columns.some(k => k.key === 'lo'));          // 신뢰구간 옵션을 켜지 않으면 열 없음
  assert.ok(!tbl(res, '모형 요약').columns.some(k => k.key === 'dw'));
  const last = Reg.runRegression(sample(), { dependent: ['y'], factors: ['g'], options: { lastRef: true } });
  assert.deepEqual(tbl(last, '계수').rows.map(r => r.name), ['(상수)', '집단 = 대조', '집단 = 처치A']);
  assert.ok(tbl(last, '계수').notes.some(n => n.includes('기준 범주: 처치B')));
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('(상수)') && html.includes('VIF') && html.includes('집단 = 처치A'));
});

test('runRegression: 입력 오류', () => {
  const d = sample();
  assert.throws(() => Reg.runRegression(d, { dependent: [], covariates: ['x'] }), /변수를 선택/);
  assert.throws(() => Reg.runRegression(d, { dependent: ['y', 'x'], covariates: ['z'] }), /하나만/);
  assert.throws(() => Reg.runRegression(d, { dependent: ['y'] }), /독립변수/);
  assert.throws(() => Reg.runRegression(d, { dependent: ['y'], covariates: ['y'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => Reg.runRegression(d, { dependent: ['y'], covariates: ['x'], factors: ['x'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => Reg.runRegression(d, { dependent: ['y'], covariates: ['x'], options: { conf: 100 } }), /신뢰수준/);
  assert.throws(() => Reg.runRegression(ds([['y', 't'], [1, 'a'], [2, 'b']]), { dependent: ['y'], covariates: ['t'] }), /숫자가 아닌 값/);
});

test('runRegression: 상수 변수, 완전 공선성, 사례 부족, 범주 수 제한', () => {
  const flat = ds([['y', 'x'], [1, 5], [2, 5], [4, 5], [3, 5]]);
  assert.throws(() => Reg.runRegression(flat, { dependent: ['y'], covariates: ['x'] }), /값이 모두 같아/);
  const col = ds([['y', 'a', 'b'], [1, 1, 2], [3, 2, 4], [2, 3, 6], [6, 4, 8], [5, 5, 10]]);
  assert.throws(() => Reg.runRegression(col, { dependent: ['y'], covariates: ['a', 'b'] }), /"b".*완전 공선성/);
  const few = ds([['y', 'a', 'b'], [1, 1, 2], [3, 2, 1], [2, 3, 5]]);
  assert.throws(() => Reg.runRegression(few, { dependent: ['y'], covariates: ['a', 'b'] }), /사례 수\(3\)가 계수의 수\(3/);
  const rows = [['y', 'g']];
  for (let i = 0; i < 24; i++) rows.push([i, i % 12]);
  assert.throws(() => Reg.runRegression(ds(rows), { dependent: ['y'], factors: ['g'] }), /12개/);
  assert.throws(() => Reg.runRegression(ds([['y', 'g'], [1, 1], [2, 1], [3, 1], [4, 1]]), { dependent: ['y'], factors: ['g'] }), /1개입니다/);
});

test('runRegression: 완전한 적합이면 경고하고 검정 통계량은 비운다', () => {
  const d = ds([['y', 'x'], [3, 1], [5, 2], [7, 3], [9, 4], [11, 5]]);
  const res = Reg.runRegression(d, { dependent: ['y'], covariates: ['x'] });
  assert.ok(res.warnings.some(w => w.includes('완전한 적합')));
  const c = tbl(res, '계수').rows;
  close(c[0].B, 1, 1e-12, 'intercept'); close(c[1].B, 2, 1e-12, 'slope');
  assert.equal(c[1].se, null);
  assert.equal(tbl(res, '분산분석').rows[0].F, null);
  assert.ok(R.renderResultHtml(res).includes('분산분석'));
});
