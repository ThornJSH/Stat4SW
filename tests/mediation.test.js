const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const S = require('../src/Sample.js');
const M = require('../src/Stats_Mediation.js');
const Ch = require('../src/Charts.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'mediation.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-12, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));

/* ---------- 난수와 분위수 ---------- */

test('mulberry32_: 시드가 같으면 같은 수열, 0~1 범위, 알려진 첫 값', () => {
  const a = M.mulberry32_(12345), b = M.mulberry32_(12345), c = M.mulberry32_(12346);
  const xa = Array.from({ length: 5 }, a), xb = Array.from({ length: 5 }, b), xc = Array.from({ length: 5 }, c);
  assert.deepEqual(xa, xb);
  assert.notDeepEqual(xa, xc);
  const many = M.mulberry32_(1);
  let sum = 0;
  for (let i = 0; i < 20000; i++) { const v = many(); assert.ok(v >= 0 && v < 1); sum += v; }
  assert.ok(Math.abs(sum / 20000 - 0.5) < 0.01);
  // 시드 0의 첫 값: 기준값 생성기(python 이식)와 같은 값이다 (부트스트랩 대조의 전제)
  assert.ok(Math.abs(M.mulberry32_(0)() - 0.26642920868471265) < 1e-15, String(M.mulberry32_(0)()));
});

test('quantileSorted_: numpy 선형 보간과 같다', () => {
  const s = [1, 2, 4, 8, 16];
  assert.equal(M.quantileSorted_(s, 0), 1);
  assert.equal(M.quantileSorted_(s, 1), 16);
  assert.equal(M.quantileSorted_(s, 0.5), 4);
  assert.ok(Math.abs(M.quantileSorted_(s, 0.025) - 1.1) < 1e-12);       // pos = 0.1 → 1 + 0.1·(2−1)
  assert.ok(Math.abs(M.quantileSorted_(s, 0.9) - 12.8) < 1e-12);        // pos = 3.6 → 8 + 0.6·(16−8)
});

test('choleskySolve_: 대칭 양정치 연립방정식, 특이행렬은 null', () => {
  const A = [[4, 2], [2, 3]];
  const [x] = M.choleskySolve_(A, [[2, 1]]);
  assert.ok(Math.abs(4 * x[0] + 2 * x[1] - 2) < 1e-12 && Math.abs(2 * x[0] + 3 * x[1] - 1) < 1e-12);
  assert.equal(M.choleskySolve_([[1, 2], [2, 4]], [[1, 1]]), null);
});

/* ---------- numpy 독립 계산과 대조 ---------- */

test('mediationFit: numpy 기준값 대조 (경로, 간접효과, Sobel, 부트스트랩 신뢰구간 — 같은 난수로 정확히 재현)', () => {
  assert.ok(ref.cases.length >= 5);
  for (const c of ref.cases) {
    const i = c.input, e = c.expected;
    const r = M.mediationFit({
      x: i.x, xName: 'X', y: i.y, yName: 'Y', conf: i.conf, boot: i.boot, seed: i.seed,
      ms: i.ms.map((v, k) => ({ name: 'M' + k, values: v })), covs: i.covs.map((v, k) => ({ name: 'C' + k, values: v }))
    });
    assert.ok(!r.error, c.name);
    const where = c.name;
    const keys = ['B', 'se', 't', 'p', 'lo', 'hi'];
    e.a.forEach((x, k) => keys.forEach(key => assert.ok(near(r.a[k][key], x[key], 1e-9), `${where} a${k}.${key}`)));
    e.b.forEach((x, k) => keys.forEach(key => assert.ok(near(r.b[k][key], x[key], 1e-9), `${where} b${k}.${key}`)));
    keys.forEach(key => {
      assert.ok(near(r.cDirect[key], e.cDirect[key], 1e-9), `${where} c' ${key}`);
      assert.ok(near(r.cTotal[key], e.cTotal[key], 1e-9), `${where} c ${key}`);
    });
    e.effects.forEach((x, k) => assert.ok(near(r.effects[k].ab, x.ab, 1e-10) && near(r.effects[k].std, x.std, 1e-10), `${where} 효과 ${k}`));
    e.sobel.forEach((x, k) => ['ab', 'se', 'z', 'p'].forEach(key => assert.ok(near(r.sobel[k][key], x[key], 1e-9), `${where} sobel ${k}.${key}`)));
    assert.ok(near(r.totalIndirect, e.totalIndirect, 1e-10) && near(r.totalStd, e.totalStd, 1e-10), where);
    assert.equal(r.boot.invalid, 0, where);
    e.bootInd.forEach((x, k) => ['lo', 'hi', 'se'].forEach(key => assert.ok(near(r.boot.ind[k][key], x[key], 1e-9), `${where} 부트스트랩 ${k}.${key} ${r.boot.ind[k][key]} vs ${x[key]}`)));
    ['lo', 'hi', 'se'].forEach(key => assert.ok(near(r.boot.total[key], e.bootTotal[key], 1e-9), `${where} 부트스트랩 합 ${key}`));
    // 항등식: 총효과 = 직접효과 + 간접효과의 합
    assert.ok(Math.abs(r.cTotal.B - (r.cDirect.B + r.totalIndirect)) < 1e-9, where);
    // 모형 R², F
    e.models.a.forEach((x, k) => assert.ok(near(r.aModels[k].r2, x.r2, 1e-9) && near(r.aModels[k].F, x.F, 1e-8), `${where} a 모형 ${k}`));
    assert.ok(near(r.yModel.r2, e.models.y.r2, 1e-9) && near(r.tModel.r2, e.models.t.r2, 1e-9), where);
  }
});

test('mediationBootstrap: 같은 시드는 같은 결과, 다른 시드는 다른 결과', () => {
  const c = ref.cases[0].input;
  const run = seed => M.mediationBootstrap(c.x, c.ms, c.covs, c.y, 300, seed);
  assert.deepEqual(run(5).total, run(5).total);
  assert.notDeepEqual(run(5).total, run(6).total);
  assert.equal(run(5).total.length, 300);
});

/* ---------- 결과 객체 ---------- */

const defs = [
  { name: 'x', label: '교육', scale: 'scale', valueLabels: {}, missing: [] },
  { name: 'm1', label: '효능', scale: 'scale', valueLabels: {}, missing: [] },
  { name: 'm2', label: '동기', scale: 'scale', valueLabels: {}, missing: [] },
  { name: 'y', label: '성과', scale: 'scale', valueLabels: {}, missing: [999] },
  { name: 'c', label: '공변량', scale: 'scale', valueLabels: {}, missing: [] },
  { name: 't', label: '글자', scale: 'nominal', valueLabels: {}, missing: [] }
];
function makeData(n) {
  let seed = 11;
  const u = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return (seed + 0.5) / 4294967296; };
  const norm = () => Math.sqrt(-2 * Math.log(u())) * Math.cos(2 * Math.PI * u());
  const rows = [['x', 'm1', 'm2', 'y', 'c', 't']];
  for (let i = 0; i < n; i++) {
    const x = 20 + 8 * norm(), c = norm();
    const m1 = 30 + 0.5 * x + 6 * norm(), m2 = 40 + 0.05 * x + 7 * norm();
    const y = 20 + 0.2 * x + 0.5 * m1 + 0.02 * m2 + 5 * norm();
    rows.push([x, m1, m2, i === 5 ? 999 : y, c, 'a']);
  }
  return Data.buildDataset(rows, defs);
}
const ds = makeData(80);
const req = (extra) => Object.assign({ dependent: ['y'], independent: ['x'], mediators: ['m1', 'm2'], options: { boot: '1000', seed: '3' } }, extra);

test('runMediation: 표 구성(모형 요약, 경로 계수, 효과 분해, Sobel)과 경로도, 결측 경고', () => {
  const res = M.runMediation(ds, req());
  assert.equal(res.title, '매개효과 분석');
  assert.equal(table(res, '모형 요약').rows.length, 4);                              // M1, M2, Y, 총효과
  assert.deepEqual(table(res, '경로 계수').rows.map(r => r.path.split(':')[0]), ['a1', 'a2', 'b1', 'b2', "c' (직접효과)", 'c (총효과)']);
  const eff = table(res, '효과 분해');
  assert.deepEqual(eff.rows.map(r => r.effect.split(':')[0].split(' (')[0]), ['총효과', '직접효과', '총 간접효과', '간접효과', '간접효과']);
  assert.deepEqual(eff.rows.map(r => !!r._bold), [false, false, true, true, true]);
  const sig = eff.rows.map(r => r.sig);
  assert.equal(sig[3], '*');                      // 효능을 통한 간접효과는 유의
  assert.equal(sig[4], '');                       // 동기를 통한 간접효과는 유의하지 않음
  assert.equal(table(res, 'Sobel 검정').rows.length, 2);
  assert.ok(res.blocks.some(b => b.type === 'chart' && b.kind === 'path' && b.spec.mediators.length === 2));
  assert.ok(res.warnings.some(w => w.includes('목록별 결측 제거로 1건')));
  assert.ok(R.renderResultHtml(res).includes('marker-end'));
  assert.ok(Ex.resultToGrid(res).rows.length > 15);
});

test('runMediation: 같은 시드는 같은 결과, 옵션(Sobel·경로도 끄기, 신뢰수준), 매개변수 1개는 "총 간접효과" 행이 없다', () => {
  const a = M.runMediation(ds, req()), b = M.runMediation(ds, req());
  assert.deepEqual(table(a, '효과 분해').rows, table(b, '효과 분해').rows);
  const c = M.runMediation(ds, req({ options: { boot: '1000', seed: '4', sobel: false, diagram: false, conf: '90' } }));
  assert.ok(!table(c, 'Sobel') && !c.blocks.some(x => x.type === 'chart'));
  assert.ok(table(c, '효과 분해').columns[3].label.includes('90%'));
  assert.notEqual(table(a, '효과 분해').rows[3].lo, table(c, '효과 분해').rows[3].lo);
  const one = M.runMediation(ds, req({ mediators: ['m1'] }));
  assert.deepEqual(table(one, '효과 분해').rows.map(r => r.effect.split(' (')[0]), ['총효과', '직접효과', '간접효과']);
  assert.equal(table(one, '경로 계수').rows[0].path.split(':')[0], 'a');            // 매개변수가 하나면 번호를 붙이지 않는다
});

test('runMediation: 공변량, 큰 표본에서 부트스트랩 횟수 자동 감축 경고', () => {
  const withCov = M.runMediation(ds, req({ covariates: ['c'] }));
  assert.ok(table(withCov, '모형 요약').rows[0].model.includes('공변량'));
  // 사례가 많으면 횟수를 줄인다 (예산 2.5e7 / 사례 수, 최소 1000)
  const big = makeData(30000);
  const res = M.runMediation(big, req({ options: { boot: '5000', seed: '1' } }));
  assert.ok(res.warnings.some(w => w.includes('부트스트랩 횟수를 5000회에서 1000회로')), res.warnings.join('\n'));
  assert.ok(table(res, '효과 분해').notes.some(n => n.includes('복원추출 1000회')));
});

test('runMediation: 입력 검증', () => {
  const run = r => M.runMediation(ds, Object.assign({ options: { boot: '500' } }, r));
  assert.throws(() => run({ dependent: [], independent: ['x'], mediators: ['m1'] }), /변수를 선택/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], mediators: [] }), /변수를 선택|매개변수를 하나 이상/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x', 'm2'], mediators: ['m1'] }), /독립변수는 하나만/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], mediators: ['x'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], mediators: ['m1'], covariates: ['m1'] }), /동시에 쓸 수 없습니다/);
  assert.throws(() => run({ dependent: ['t'], independent: ['x'], mediators: ['m1'] }), /숫자가 아닌 값/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], mediators: ['m1'], options: { boot: '100' } }), /부트스트랩 횟수/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], mediators: ['m1'], options: { boot: '600.5' } }), /부트스트랩 횟수/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], mediators: ['m1'], options: { boot: '500', seed: 'a' } }), /시드/);
  assert.throws(() => run({ dependent: ['y'], independent: ['x'], mediators: ['m1'], options: { boot: '500', conf: '0' } }), /신뢰수준/);
  const tiny = Data.buildDataset([['y', 'x', 'm'], [1, 2, 3], [2, 3, 5], [3, 5, 4]], []);
  assert.throws(() => M.runMediation(tiny, { dependent: ['y'], independent: ['x'], mediators: ['m'], options: { boot: '500' } }), /사례 수/);
  const same = Data.buildDataset([['y', 'x', 'm'], [1, 2, 3], [2, 2, 5], [3, 2, 4], [4, 2, 9], [5, 2, 1], [6, 2, 2]], []);
  assert.throws(() => M.runMediation(same, { dependent: ['y'], independent: ['x'], mediators: ['m'], options: { boot: '500' } }), /모두 같아/);
});

test('경로도: SVG 에 상자, 화살표, 계수 문구가 있다', () => {
  const block = Ch.pathDiagramBlock('경로도', { x: '교육', y: '성과', mediators: [{ label: '효능', a: '0.352***', b: '0.447***' }, { label: '동기', a: '0.334**', b: '0.108' }], direct: '0.333**', total: '0.526***' });
  const html = Ch.renderChartHtml(block);
  ['교육', '성과', '효능', '동기', 'a1 = 0.352***', 'b2 = 0.108', "c&#39; = 0.333**", 'c = 0.526***'].forEach(t => assert.ok(html.includes(t), t));
  assert.equal((html.match(/marker-end/g) || []).length, 5);                    // 화살표: a×2, b×2, 직접 1
  assert.ok(Ch.renderChartHtml(Ch.pathDiagramBlock('빈', { x: 'X', y: 'Y', mediators: [] })).includes('그릴 데이터가 없습니다'));
});

/* ---------- 예제 데이터 ---------- */

test('예제 데이터: 교육시간 → (자기효능감, 동기수준) → 직무성과: 효능은 유의한 간접효과, 동기는 아님, 부분 매개', () => {
  const sample = S.sampleDataset_();
  const parsed = Data.parseVarDefs([['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(sample.varRows));
  const full = Data.buildDataset([sample.header].concat(sample.rows), parsed.defs);
  const res = M.runMediation(full, { dependent: ['직무성과'], independent: ['교육시간'], mediators: ['자기효능감', '동기수준'], options: { boot: '5000', seed: '12345' } });
  const eff = table(res, '효과 분해').rows;
  const by = name => eff.find(r => r.effect.includes(name));
  assert.equal(by('자기효능감').sig, '*');
  assert.equal(by('동기 수준').sig, '');
  assert.equal(by('총 간접효과').sig, '*');
  assert.ok(by("직접효과").sig.length >= 1);                                     // 직접효과도 유의 → 부분 매개
  assert.ok(by('총효과').B > by('직접효과').B);
  assert.deepEqual(sample.header.slice(22, 26), ['교육시간', '자기효능감', '동기수준', '직무성과']);
});
