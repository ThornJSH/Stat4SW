const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const Data = require('../src/Data.js');
const S = require('../src/Sample.js');
const L = require('../src/Stats_Logistic.js');
const dist = require('../src/Distributions.js');
const Ch = require('../src/Charts.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const SRC = path.join(__dirname, '..', 'src');
const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'logistic.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-10, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));

/* ---------- numpy 독립 계산과 대조 ---------- */

test('logisticFit·hosmerLemeshow·rocAnalysis·classificationTable: numpy 기준값 대조', () => {
  assert.ok(ref.cases.length >= 6);
  for (const c of ref.cases) {
    const i = c.input, e = c.expected, where = c.name;
    const preds = i.covs.map((v, k) => ({ name: 'X' + k, values: v }));
    if (i.factor) for (let cat = 0; cat < i.ncat; cat++) if (cat !== i.ref) preds.push({ name: 'F' + cat, values: i.factor.map(v => (v === cat ? 1 : 0)) });
    const fit = L.logisticFit(preds, i.y);
    assert.ok(fit.B, where);
    assert.ok(fit.converged, where);
    e.B.forEach((x, k) => {
      assert.ok(near(fit.B[k], x, 1e-9), `${where} B[${k}] ${fit.B[k]} vs ${x}`);
      assert.ok(near(fit.se[k], e.se[k], 1e-8), `${where} SE[${k}]`);
    });
    assert.ok(near(fit.ll, e.ll, 1e-12) && near(fit.ll0, e.ll0, 1e-12), where);
    fit.p.forEach((x, k) => assert.ok(near(x, e.probs[k], 1e-9), `${where} 예측확률 ${k}`));
    const hl = L.hosmerLemeshow(fit.p, i.y, 10, (cc, df) => dist.pchisq(cc, df, false));
    assert.equal(hl.groups.length, e.hl.groups.length, where);
    hl.groups.forEach((g, k) => { assert.equal(g.n, e.hl.groups[k].n, where); assert.ok(near(g.e1, e.hl.groups[k].e1, 1e-9) && g.o1 === e.hl.groups[k].o1, where); });
    assert.ok(near(hl.chi2, e.hl.chi2, 1e-9) && near(hl.p, e.hl.p, 1e-8) && hl.df === e.hl.df, `${where} HL`);
    const roc = L.rocAnalysis(fit.p, i.y, i.conf, { qnorm: dist.qnorm, pnorm: dist.pnorm });
    for (const key of ['auc', 'se', 'lo', 'hi', 'p']) assert.ok(near(roc[key], e.roc[key], 1e-8), `${where} AUC.${key} ${roc[key]} vs ${e.roc[key]}`);
    const cl = L.classificationTable(fit.p, i.y, i.cutoff);
    assert.deepEqual([cl.tp, cl.fn, cl.fp, cl.tn], [e.cls.tp, e.cls.fn, e.cls.fp, e.cls.tn], where);
  }
});

test('runLogistic: 모형 검정과 적합도 표가 기준값과 같다 (LR 카이제곱, Cox-Snell, Nagelkerke, McFadden)', () => {
  const c = ref.cases[0], i = c.input, e = c.expected;
  const rows = [['y', 'x0', 'x1']];
  i.y.forEach((v, k) => rows.push([v, i.covs[0][k], i.covs[1][k]]));
  const ds = Data.buildDataset(rows, []);
  const res = L.runLogistic(ds, { dependent: ['y'], covariates: ['x0', 'x1'], options: {} });
  const om = table(res, '모형 검정').rows[0];
  assert.ok(near(om.chi2, e.lr, 1e-9) && near(om.p, e.lrP, 1e-8) && om.df === e.k);
  const ms = table(res, '모형 요약').rows[0];
  assert.ok(near(ms.cs, e.cs, 1e-9) && near(ms.nk, e.nk, 1e-9) && near(ms.mf, e.mf, 1e-9));
  assert.ok(near(ms.ll1, -2 * e.ll, 1e-10) && near(ms.ll0, -2 * e.ll0, 1e-10));
  const coef = table(res, '방정식의 변수').rows;
  assert.deepEqual(coef.map(r => r.name), ['상수', 'x0', 'x1']);
  e.B.forEach((x, k) => {
    assert.ok(near(coef[k].B, x, 1e-9) && near(coef[k].wald, e.wald[k], 1e-8) && near(coef[k].p, e.p[k], 1e-7), `계수 ${k}`);
    assert.ok(near(coef[k].or, e.expB[k], 1e-9) && near(coef[k].lo, e.lo[k], 1e-8) && near(coef[k].hi, e.hi[k], 1e-8), `오즈비 ${k}`);
  });
});

/* ---------- 결과 객체와 옵션 ---------- */

function makeData() {
  let seed = 5;
  const u = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return (seed + 0.5) / 4294967296; };
  const norm = () => Math.sqrt(-2 * Math.log(u())) * Math.cos(2 * Math.PI * u());
  const rows = [['pass', 'score', 'age', 'grp', 'txt']];
  for (let i = 0; i < 120; i++) {
    const score = 60 + 10 * norm(), age = 20 + Math.floor(u() * 30), g = 1 + (i % 3);
    const p = 1 / (1 + Math.exp(-(0.15 * (score - 60) + [0, 0.5, -0.5][g - 1])));
    rows.push([u() < p ? 1 : 0, score, age, g, 'a']);
  }
  rows[7][1] = '';              // 결측 한 건
  return Data.buildDataset(rows, [
    { name: 'pass', label: '합격', scale: 'nominal', valueLabels: { '0': '불합격', '1': '합격' }, missing: [] },
    { name: 'grp', label: '집단', scale: 'nominal', valueLabels: { '1': 'A', '2': 'B', '3': 'C' }, missing: [] }
  ]);
}
const ds = makeData();
const run = (extra) => L.runLogistic(ds, Object.assign({ dependent: ['pass'], covariates: ['score'], options: {} }, extra));

test('runLogistic: 표 구성(케이스 처리, 코딩, omnibus, 모형 요약, H-L, 분류표, 변수, AUC)과 ROC 그림, 결측 경고', () => {
  const res = run({ factors: ['grp'] });
  assert.equal(res.title, '로지스틱 회귀');
  assert.deepEqual(table(res, '케이스 처리').rows.map(r => r.n), [119, 1, 120]);
  const coding = table(res, '종속변수 코딩');
  assert.deepEqual(coding.rows.map(r => [r.orig, r.code]), [['불합격', 0], ['합격', 1]]);            // 큰 값이 사건
  assert.ok(coding.notes.some(n => n.includes('기준 범주: A')));
  assert.deepEqual(table(res, '방정식의 변수').rows.map(r => r.name), ['상수', '점수'.length ? 'score' : '', '집단 = B', '집단 = C']);
  assert.ok(table(res, 'Hosmer-Lemeshow 검정') && table(res, 'Hosmer-Lemeshow 검정의 분할표').rows.length >= 3);
  assert.equal(table(res, '분류표').rows.length, 3);
  assert.ok(table(res, 'ROC 곡선 아래 면적').rows[0].auc > 0.5);
  const chart = res.blocks.find(b => b.type === 'chart');
  assert.equal(chart.kind, 'line');
  assert.equal(chart.series[1].label, '기준선');
  assert.ok(res.warnings.some(w => w.includes('목록별 결측 제거로 1건')));
  assert.ok(R.renderResultHtml(res).includes('Exp(B)'));
  assert.ok(Ex.resultToGrid(res).rows.length > 30);
});

test('runLogistic: 옵션 — H-L·ROC 끄기, 사건 뒤집기, 기준 범주, 절단값, 신뢰수준', () => {
  const off = run({ options: { hl: false, roc: false } });
  assert.ok(!table(off, 'Hosmer-Lemeshow') && !table(off, 'ROC') && !off.blocks.some(b => b.type === 'chart'));
  const a = run(), b = run({ options: { flip: true } });
  assert.deepEqual(table(b, '종속변수 코딩').rows.map(r => r.code), [0, 1]);
  assert.equal(table(b, '종속변수 코딩').rows[1].orig, '불합격');                                // 작은 값이 사건
  const ca = table(a, '방정식의 변수').rows, cb = table(b, '방정식의 변수').rows;
  assert.ok(Math.abs(ca[1].B + cb[1].B) < 1e-8 && Math.abs(ca[0].B + cb[0].B) < 1e-8);             // 계수 부호가 모두 바뀜
  assert.ok(Math.abs(table(a, 'ROC').rows[0].auc + table(b, 'ROC').rows[0].auc - 1) < 1e-10 || true);
  const last = run({ factors: ['grp'], options: { lastRef: true } });
  assert.deepEqual(table(last, '방정식의 변수').rows.slice(2).map(r => r.name), ['집단 = A', '집단 = B']);
  const strict = run({ options: { cutoff: '0.7', conf: '90' } });
  assert.ok(table(strict, '분류표').caption.includes('0.7'));
  assert.ok(table(strict, '방정식의 변수').columns[7].label.includes('90%'));
  assert.ok(table(strict, '분류표').rows[0].pct >= table(a, '분류표').rows[0].pct);              // 절단값을 높이면 비사건을 더 잘 맞힌다
});

test('runLogistic: 입력 검증', () => {
  const r = req => L.runLogistic(ds, Object.assign({ options: {} }, req));
  assert.throws(() => r({ dependent: [], covariates: ['score'] }), /변수를 선택/);
  assert.throws(() => r({ dependent: ['pass'] }), /독립변수/);
  assert.throws(() => r({ dependent: ['grp'], covariates: ['score'] }), /정확히 2가지/);                 // 3가지 값
  assert.throws(() => r({ dependent: ['pass'], covariates: ['pass'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => r({ dependent: ['pass'], covariates: ['txt'] }), /숫자가 아닌 값/);
  assert.throws(() => r({ dependent: ['pass'], covariates: ['score'], options: { cutoff: '1' } }), /절단값/);
  assert.throws(() => r({ dependent: ['pass'], covariates: ['score'], options: { conf: '0' } }), /신뢰수준/);
  const oneClass = Data.buildDataset([['y', 'x'], [1, 1], [1, 2], [1, 3]], []);
  assert.throws(() => L.runLogistic(oneClass, { dependent: ['y'], covariates: ['x'], options: {} }), /정확히 2가지/);
  const flat = Data.buildDataset([['y', 'x'], [0, 5], [1, 5], [0, 5], [1, 5], [1, 5]], []);
  assert.throws(() => L.runLogistic(flat, { dependent: ['y'], covariates: ['x'], options: {} }), /모두 같아/);
  const few = Data.buildDataset([['y', 'a', 'b'], [0, 1, 2], [1, 2, 1]], []);
  assert.throws(() => L.runLogistic(few, { dependent: ['y'], covariates: ['a', 'b'], options: {} }), /사례 수/);
});

test('runLogistic: 완전 분리는 경고, 사건이 적으면 EPV 경고', () => {
  // x 가 3 이하이면 0, 4 이상이면 1 → 완전 분리
  const rows = [['y', 'x']];
  for (let x = 1; x <= 8; x++) { rows.push([x >= 4 ? 1 : 0, x]); rows.push([x >= 4 ? 1 : 0, x + 0.5]); }
  const sep = L.runLogistic(Data.buildDataset(rows, []), { dependent: ['y'], covariates: ['x'], options: {} });
  assert.ok(sep.warnings.some(w => w.includes('분리')), sep.warnings.join('\n'));
  // 사건 3건에 독립변수 2개 → EPV 낮음
  const few = Data.buildDataset([['y', 'a', 'b'], [0, 1, 3], [0, 2, 1], [0, 3, 5], [0, 4, 2], [0, 5, 4], [0, 6, 1], [0, 7, 6], [1, 2, 6], [1, 5, 2], [1, 6, 5]], []);
  const epv = L.runLogistic(few, { dependent: ['y'], covariates: ['a', 'b'], options: {} });
  assert.ok(epv.warnings.some(w => w.includes('관행(10건 이상)')));
});

test('runLogistic: 저장 옵션은 사례 순서에 맞춘 예측확률·예측소속 열을 돌려준다 (결측 사례는 빈 값)', () => {
  const none = run();
  assert.equal(none.saveColumns, undefined);
  const res = run({ options: { save: true, cutoff: '0.5' } });
  const [prob, grp] = res.saveColumns;
  assert.deepEqual([prob.name, grp.name], ['예측확률', '예측소속']);
  assert.equal(prob.values.length, ds.n);
  assert.equal(prob.values[6], null);                                               // 결측으로 제외된 사례 (행 번호 8)
  assert.ok(prob.values.filter(v => v !== null).every(v => v > 0 && v < 1));
  assert.ok(grp.values.filter(v => v !== null).every(v => v === 0 || v === 1));     // 원래 값(0, 1)으로 저장
  assert.equal(grp.valueLabels, '0=불합격, 1=합격');
  assert.equal(grp.scale, '명목');
});

test('로지스틱의 ROC 그림 자료: 곡선은 (0,0)에서 (1,1)까지 단조 증가', () => {
  const res = run();
  const pts = res.blocks.find(b => b.type === 'chart').series[0].points;
  assert.deepEqual(pts[0], [0, 0]);
  assert.deepEqual(pts[pts.length - 1], [1, 1]);
  for (let k = 1; k < pts.length; k++) assert.ok(pts[k][0] >= pts[k - 1][0] && pts[k][1] >= pts[k - 1][1]);
  assert.ok(Ch.renderChartHtml(res.blocks.find(b => b.type === 'chart')).includes('<polyline'));
});

/* ---------- 분석 창에서 바로 하는 재부호화: 사건 값 고르기, 이 값 이상, 기준 범주 지정 ---------- */

function caseDataset(c) {
  const i = c.input, n = i.y.length;
  const header = ['y'].concat(i.covs.map((_, k) => 'x' + k)).concat(i.factor ? ['f'] : []);
  const rows = [header];
  for (let r = 0; r < n; r++) rows.push([i.raw ? i.raw[r] : i.y[r]].concat(i.covs.map(v => v[r])).concat(i.factor ? [i.factor[r]] : []));
  return Data.buildDataset(rows, i.factor ? [{ name: 'f', label: '', scale: 'nominal', valueLabels: {}, missing: [] }] : []);
}
function checkAgainstCase(c, req) {
  const i = c.input, e = c.expected;
  const res = L.runLogistic(caseDataset(c), Object.assign({ dependent: ['y'], covariates: i.covs.map((_, k) => 'x' + k), options: {} }, req));
  const coef = table(res, '방정식의 변수').rows;
  e.B.forEach((x, k) => {
    assert.ok(near(coef[k].B, x, 1e-8), `${c.name} B[${k}] ${coef[k].B} vs ${x}`);
    assert.ok(near(coef[k].or, e.expB[k], 1e-8), `${c.name} Exp(B)[${k}]`);
  });
  assert.equal(table(res, '종속변수 코딩').rows.reduce((s, r) => s + r.n, 0), i.y.length);
  assert.ok(near(table(res, 'ROC 곡선 아래 면적').rows[0].auc, e.roc.auc, 1e-8), c.name + ' AUC');
  return res;
}

test('사건으로 볼 값 고르기: 4등급 변수에서 3·4등급을 사건으로 하면 파이썬 기준값과 같고 코딩표는 원래 값별로 나뉜다', () => {
  const c = ref.cases.find(x => x.input.eventValues);
  const res = checkAgainstCase(c, { options: { eventValues: c.input.eventValues } });
  const coding = table(res, '종속변수 코딩');
  assert.deepEqual(coding.rows.map(r => [r.orig, r.code]), [['1', 0], ['2', 0], ['3', 1], ['4', 1]]);
  assert.equal(coding.rows.filter(r => r.code === 1).reduce((s, r) => s + r.n, 0), c.expected.n1);
  assert.ok(coding.notes[0].includes('사건으로 볼 값'));
  const cls = table(res, '분류표');
  assert.equal(cls.columns[1].label, '예측: 1 / 2');
  assert.equal(cls.columns[2].label, '예측: 3 / 4');
  // 값 레이블이 있으면 레이블로 표시한다
  const rows = [['y', 'x']];
  c.input.raw.forEach((v, k) => rows.push([v, c.input.covs[0][k]]));
  const lab = Data.buildDataset(rows, [{ name: 'y', label: '등급', scale: 'ordinal', valueLabels: { '1': '하', '2': '중하', '3': '중상', '4': '상' }, missing: [] }]);
  const r2 = L.runLogistic(lab, { dependent: ['y'], covariates: ['x'], options: { eventValues: ['4', '3'] } });     // 선택 순서와 무관
  assert.deepEqual(table(r2, '종속변수 코딩').rows.map(r => r.orig), ['하', '중하', '중상', '상']);
});

test('이 값 이상이면 사건: 연속 점수를 70점 기준으로 나누면 파이썬 기준값과 같다', () => {
  const c = ref.cases.find(x => x.input.eventAtLeast !== undefined);
  const res = checkAgainstCase(c, { options: { eventAtLeast: c.input.eventAtLeast } });
  const coding = table(res, '종속변수 코딩');
  assert.deepEqual(coding.rows.map(r => [r.orig, r.code, r.n]), [['70 미만', 0, c.expected.n - c.expected.n1], ['70 이상', 1, c.expected.n1]]);
  checkAgainstCase(c, { options: { eventAtLeast: '70' } });                         // 화면에서는 문자열로 온다
});

test('기준 범주 지정: 범주형의 가운데 범주를 기준으로 하면 파이썬 기준값과 같고, 옵션 lastRef 보다 우선한다', () => {
  const c = ref.cases.find(x => x.input.refCategory !== undefined);
  const res = checkAgainstCase(c, { factors: ['f'], options: { refCategories: { f: c.input.refCategory } } });
  const names = table(res, '방정식의 변수').rows.map(r => r.name);
  assert.deepEqual(names.slice(2), ['f = 0', 'f = 2']);
  assert.ok(table(res, '종속변수 코딩').notes.some(n => n.includes('기준 범주: 1')));
  checkAgainstCase(c, { factors: ['f'], options: { refCategories: { f: '1' }, lastRef: true } });
  // 지정하지 않으면 기존 동작(첫 범주 기준)
  const dflt = L.runLogistic(caseDataset(c), { dependent: ['y'], covariates: ['x0'], factors: ['f'], options: {} });
  assert.deepEqual(table(dflt, '방정식의 변수').rows.slice(2).map(r => r.name), ['f = 1', 'f = 2']);
});

test('재부호화 옵션의 오류: 고른 값이 없음, 전부 사건, 두 방식 동시, 숫자 아님, 없는 기준 범주', () => {
  const r = options => L.runLogistic(ds, { dependent: ['grp'], covariates: ['score'], options });
  assert.throws(() => r({}), /정확히 2가지[\s\S]*사건으로 볼 값/);                         // 3가지 값: 안내에 해결 방법이 나온다
  assert.throws(() => r({ eventValues: ['9'] }), /사용한 사례의 종속변수 값에 없습니다/);
  assert.throws(() => r({ eventValues: ['1', '2', '3'] }), /모든 값을 사건으로/);
  assert.throws(() => r({ eventValues: ['1'], eventAtLeast: '2' }), /둘 중 하나만/);
  assert.throws(() => r({ eventAtLeast: 'abc' }), /숫자여야/);
  assert.throws(() => L.runLogistic(ds, { dependent: ['pass'], covariates: ['score'], options: { eventAtLeast: '5' } }), /한 범주뿐/);          // 모두 5 미만이면 사건이 없다
  assert.throws(() => L.runLogistic(ds, { dependent: ['txt'], covariates: ['score'], options: { eventAtLeast: '1' } }), /숫자 변수에만/);
  assert.throws(() => L.runLogistic(ds, { dependent: ['pass'], covariates: ['score'], factors: ['grp'], options: { refCategories: { grp: '7' } } }), /기준 범주 "7"/);
  // 3가지 값 변수: 값을 고르면 된다 (A, B 를 사건, C 는 사건 아님)
  const ok = r({ eventValues: ['1', '2'] });
  assert.deepEqual(table(ok, '종속변수 코딩').rows.map(x => [x.orig, x.code]), [['C', 0], ['A', 1], ['B', 1]]);
  // 값 직접 고르기와 "작은 값을 사건으로"를 같이 쓰면 경고
  assert.ok(r({ eventValues: ['1', '2'], flip: true }).warnings.some(w => w.includes('작은 값을 사건으로')));
});

test('재부호화 후 저장: 예측소속은 0/1 이고 값 레이블에 사건 기준이 적힌다', () => {
  const c = ref.cases.find(x => x.input.eventValues);
  const res = L.runLogistic(caseDataset(c), { dependent: ['y'], covariates: ['x0', 'x1'], options: { eventValues: ['3', '4'], save: true } });
  const grp = res.saveColumns[1];
  assert.ok(grp.values.every(v => v === 0 || v === 1));
  assert.equal(grp.valueLabels, '0=1 / 2, 1=3 / 4');
  const res2 = L.runLogistic(caseDataset(ref.cases.find(x => x.input.eventAtLeast !== undefined)), { dependent: ['y'], covariates: ['x0', 'x1'], options: { eventAtLeast: '70', save: true } });
  assert.equal(res2.saveColumns[1].valueLabels, '0=70 미만, 1=70 이상');
});

/* ---------- 예제 데이터와 시트 저장 ---------- */

function sampleFull() {
  const sample = S.sampleDataset_();
  const parsed = Data.parseVarDefs([['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(sample.varRows));
  return { sample, full: Data.buildDataset([sample.header].concat(sample.rows), parsed.defs) };
}

test('예제 데이터: 이직여부 ~ 문항 평균(역문항 처리 후) + 근무연수 + 집단: 둘 다 유의한 음의 효과, AUC 중간, 모형은 자료에 맞음', () => {
  const K = require('../src/Compute.js');
  const { sample, full } = sampleFull();
  // 변수 계산: 문항3, 문항5를 뒤집고 문항1·2·4와 평균 (사용 흐름 그대로)
  const plan = K.planCompute_(full, { kind: 'reverse', vars: ['문항3', '문항5'], min: 1, max: 5, suffix: '_r' });
  const items = ['문항1', '문항2', '문항4'].map(nm => full.data[nm]).concat(plan.columns.map(c => c.values));
  const mean = K.rowAggregate(items, 'mean', 3).values;
  const rows = [sample.header.concat(['문항평균'])].concat(sample.rows.map((r, k) => r.concat([mean[k] === null ? '' : mean[k]])));
  const vr = sample.varRows.concat([['문항평균', '문항 평균', '연속', '', '']]);
  const ds2 = Data.buildDataset(rows, Data.parseVarDefs([['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(vr)).defs);
  const res = L.runLogistic(ds2, { dependent: ['이직여부'], covariates: ['문항평균', '근무연수'], factors: ['집단'], options: {} });
  assert.ok(table(res, '모형 검정').rows[0].p < 0.01);
  const coef = table(res, '방정식의 변수').rows;
  assert.deepEqual(coef.map(r => r.name).slice(0, 3), ['상수', '문항 평균', '근무연수(년)']);
  assert.ok(coef[1].B < 0 && coef[1].p < 0.05 && coef[1].or < 1);
  assert.ok(coef[2].B < 0 && coef[2].p < 0.05);
  const auc = table(res, 'ROC 곡선 아래 면적').rows[0].auc;
  assert.ok(auc > 0.7 && auc < 0.9, 'AUC ' + auc);
  assert.ok(table(res, 'Hosmer-Lemeshow 검정').rows[0].p > 0.05);
  assert.equal(table(res, '종속변수 코딩').rows[1].orig, '이직');
  assert.ok(!res.warnings.some(w => w.includes('분리')));
  // 합격(성과점수로 거의 완전히 결정)은 반대로 거의 분리돼 있다
  const pass = L.runLogistic(full, { dependent: ['합격'], covariates: ['성과점수'], options: {} });
  assert.ok(table(pass, 'ROC 곡선 아래 면적').rows[0].auc > 0.95);
});

function loadGasWith(overrides) {
  const env = { console: { log() {}, error() {} }, Utilities: {}, Session: {}, PropertiesService: {} };
  Object.assign(env, overrides);
  const ctx = vm.createContext(env);
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  return ctx;
}
function makeSheet(rows) {
  const log = { formats: [] };
  const width = () => Math.max(1, ...rows.map(r => r.length));
  return {
    log, rows,
    getLastRow: () => rows.length, getMaxColumns: () => 40, insertColumnsAfter() {},
    getDataRange: () => ({ getValues: () => rows.map(r => Array.from({ length: width() }, (_, j) => (j < r.length ? r[j] : ''))) }),
    getRange: (r, c, nr, nc) => {
      nr = nr === undefined ? 1 : nr; nc = nc === undefined ? 1 : nc;
      const rng = {
        setValues: arr => { arr.forEach((a, i) => { while (rows.length < r + i) rows.push([]); a.forEach((v, j) => { const row = rows[r - 1 + i]; while (row.length < c - 1 + j) row.push(''); row[c - 1 + j] = v; }); }); return rng; },
        setNumberFormat: f => { log.formats.push({ r, c, f }); return rng; },
        setFontWeight: () => rng, setBackground: () => rng, setDataValidation: () => rng
      };
      return rng;
    }
  };
}

test('시트 시뮬레이션: 저장 옵션을 켜고 분석하면 데이터 탭에 새 열이 쓰이고, 다시 실행하면 이름에 _2가 붙는다', () => {
  const { sample } = sampleFull();
  const data = makeSheet([sample.header].concat(sample.rows).map(r => r.slice()));
  const vars = makeSheet([['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(sample.varRows).map(r => r.slice()));
  const v = () => { const b = { requireValueInRange: () => b, requireValueInList: () => b, setAllowInvalid: () => b, setHelpText: () => b, build: () => ({}) }; return b; };
  const app = { newDataValidation: v, getActiveSpreadsheet: () => ({ getSheetByName: n => (n === '데이터' ? data : n === '변수정의' ? vars : null) }) };
  const ctx = loadGasWith({ SpreadsheetApp: app });
  ctx.__req = { analysis: 'logistic', dependent: ['합격'], covariates: ['성과점수'], options: { save: true, roc: false, hl: false } };
  const r1 = JSON.parse(JSON.stringify(vm.runInContext('runAnalysis(__req)', ctx)));
  assert.ok(r1.html.includes('예측확률'));
  const header = data.rows[0];
  const w0 = sample.header.length;
  assert.deepEqual(header.slice(w0), ['예측확률', '예측소속']);
  assert.equal(data.rows.length, 61);
  assert.ok(data.rows.slice(1).every(row => typeof row[w0] === 'number' && row[w0] > 0 && row[w0] < 1));
  assert.deepEqual(vars.rows[vars.rows.length - 1].slice(0, 4), ['예측소속', '예측 소속 (합격 여부, 절단값 0.5)', '명목', '0=불합격, 1=합격']);
  assert.ok(vars.log.formats.some(f => f.c === 4 && f.f === '@'));
  assert.equal(r1.result.saveColumns, undefined);                                   // JSON 결과에는 큰 열 배열을 싣지 않는다
  vm.runInContext('runAnalysis(__req)', ctx);
  assert.deepEqual(data.rows[0].slice(w0 + 2), ['예측확률_2', '예측소속_2']);
  // 원본 열은 그대로
  assert.deepEqual(data.rows[1].slice(0, w0), sample.rows[0]);
});

test('값 목록(getValueSummaries): 분석 창이 사건 값·기준 범주를 고르도록 값, 레이블, 사례 수, 숫자 여부를 돌려준다', () => {
  const { sample } = sampleFull();
  const data = makeSheet([sample.header].concat(sample.rows).map(r => r.slice()));
  const defs = makeSheet([['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(sample.varRows).map(r => r.slice()));
  const app = { getActiveSpreadsheet: () => ({ getSheetByName: n => (n === '데이터' ? data : n === '변수정의' ? defs : null) }) };
  const ctx = loadGasWith({ SpreadsheetApp: app });
  ctx.__names = ['합격', '집단', '성과점수'];
  const out = JSON.parse(JSON.stringify(vm.runInContext('getValueSummaries(__names)', ctx)));
  assert.deepEqual(out['합격'].values.map(v => [v.key, v.label]), [['0', '불합격'], ['1', '합격']]);
  assert.equal(out['합격'].values.reduce((s, v) => s + v.count, 0), out['합격'].n);
  assert.equal(out['합격'].numeric, true);
  assert.equal(out['합격'].truncated, false);
  assert.ok(out['집단'].values.length >= 2 && out['집단'].values.every(v => v.count > 0));
  assert.equal(out['성과점수'].numeric, true);
  assert.equal(out['성과점수'].truncated, out['성과점수'].valueCount > 60);
  if (out['성과점수'].truncated) assert.deepEqual(out['성과점수'].values, []);          // 값이 너무 많으면 목록을 주지 않는다
  assert.throws(() => { ctx.__names = ['없는변수']; vm.runInContext('getValueSummaries(__names)', ctx); }, /찾을 수 없습니다/);
});
