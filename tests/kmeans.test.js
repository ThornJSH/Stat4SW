const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const K = require('../src/Stats_Kmeans.js');
const Dd = require('../src/Stats_Distance.js');
const Med = require('../src/Stats_Mediation.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'kmeans.json'), 'utf8'));
const near = (a, e, tol) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1e-10, Math.abs(e)));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const names = ['연령', '소득', '점수'];
const ds = (rows) => Data.buildDataset([names].concat(rows), names.map(n => ({ name: n, label: '', scale: 'scale', valueLabels: {}, missing: [] })));
const dsFull = ds(ref.data);
const workSpace = (std) => {
  const cols = names.map((_, j) => ref.data.map(r => r[j]));
  const sc = Dd.standardizeColumns(cols, std, names);
  return ref.data.map((_, i) => sc.map(c => c[i]));
};

test('kmeansRun: 파이썬 직접 구현(+scipy kmeans2)과 같은 초기 중심, 군집 배정, 중심, 제곱합, 반복 수 — 3가지 설정', () => {
  for (const c of ref.cases) {
    const X = workSpace(c.standardize);
    const r = K.kmeansRun(X, c.K, c.seed, c.starts, 100);
    assert.deepEqual(r.init, c.init, c.name + ' 초기 중심');
    assert.equal(r.bestStart, c.best_start, c.name);
    assert.deepEqual(r.labels, c.labels, c.name + ' 배정');
    assert.deepEqual(r.sizes, c.sizes, c.name);
    r.centers.forEach((row, g) => row.forEach((v, j) => assert.ok(near(v, c.centers[g][j], 1e-9), `${c.name} 중심 ${g},${j}`)));
    assert.ok(near(r.wss, c.wss, 1e-10), c.name + ' WSS');
    assert.equal(r.iterations, c.iters, c.name + ' 반복 수');
    assert.equal(r.converged, c.converged);
    assert.ok(near(K.kmeansSilhouette(X, r.labels, c.K), c.silhouette, 1e-9), c.name + ' 실루엣');
  }
});

test('mulberry32 와 k-means++: 같은 시드는 같은 수열(파이썬 구현과 같은 첫 값들)', () => {
  const a = Med.mulberry32_(12345), b = Med.mulberry32_(12345);
  for (let i = 0; i < 5; i++) assert.equal(a(), b());
  const X = workSpace('z');
  assert.deepEqual(K.kmeansPlusPlus(X, 3, Med.mulberry32_(99)), K.kmeansPlusPlus(X, 3, Med.mulberry32_(99)));
  assert.notDeepEqual(K.kmeansPlusPlus(X, 3, Med.mulberry32_(99)), K.kmeansPlusPlus(X, 3, Med.mulberry32_(100)));
});

test('runKmeans: 표 구성과 값이 기준값과 같다 (중심, 사례 수, 제곱합, 군집 사이 거리, 변수별 F, 엘보)', () => {
  const c = ref.cases[0];
  const res = K.runKmeans(dsFull, { vars: names, options: { k: '3', seed: '12345', starts: '10', standardize: 'z', elbow: true } });
  assert.equal(res.title, 'k-평균 군집분석');
  const sum = table(res, '군집 결과 요약').rows;
  assert.equal(sum[0].value, 3);
  assert.ok(near(sum[4].value, c.tss, 1e-9) && near(sum[5].value, c.wss, 1e-9) && near(sum[6].value, 1 - c.wss / c.tss, 1e-9));
  const cen = table(res, '최종 군집 중심 (표준화한 값)').rows;
  c.centers.forEach((row, g) => row.forEach((v, j) => assert.ok(near(cen[g]['v' + j], v, 1e-9))));
  const orig = table(res, '최종 군집 중심 (원래 단위)').rows;
  c.centers_orig.forEach((row, g) => row.forEach((v, j) => assert.ok(near(orig[g]['v' + j], v, 1e-9))));
  const cnt = table(res, '군집별 사례 수').rows;
  assert.deepEqual(cnt.slice(0, 3).map(r => r.n), c.sizes);
  c.wss_by_cluster.forEach((v, g) => assert.ok(near(cnt[g].wss, v, 1e-9)));
  const cd = table(res, '군집 중심 사이의 거리').rows;
  c.center_dist.forEach((row, g) => row.forEach((v, h) => assert.ok(near(cd[g]['c' + h], v, 1e-9) || Math.abs(cd[g]['c' + h] - v) < 1e-12)));
  const f = table(res, '변수별 군집 간 차이').rows;
  c.f.forEach((e, j) => assert.ok(near(f[j].F, e.F, 1e-8) && near(f[j].p, e.p, 1e-6) && f[j].df1 === e.df1 && f[j].df2 === e.df2, '변수 ' + j));
  const el = table(res, '군집 수 판단 보조').rows;
  assert.equal(el.length, 9);
  c.elbow.forEach((e, i) => assert.ok(near(el[i].wss, e.wss, 1e-9) && near(el[i].sil, e.silhouette, 1e-9) && el[i].k === e.k, 'elbow k=' + e.k));
  assert.ok(el[1]._bold && !el[0]._bold);                                                       // k=3 줄이 굵게
  assert.equal(res.blocks.find(b => b.type === 'chart').kind, 'line');
  assert.ok(R.renderResultHtml(res).includes('군집 중심'));
  assert.ok(Ex.resultToGrid(res).rows.length > 30);
});

test('runKmeans: 저장 열(군집번호, 중심거리)은 사례 순서에 맞고 결측 사례는 빈 값, 같은 시드는 같은 결과', () => {
  const rows = ref.data.map(r => r.slice());
  rows[4][1] = '';
  const d2 = ds(rows);
  const opts = { k: '3', seed: '5', starts: '3', save: true };
  const a = K.runKmeans(d2, { vars: names, options: opts }), b = K.runKmeans(d2, { vars: names, options: opts });
  assert.deepEqual(a.saveColumns[0].values, b.saveColumns[0].values);
  const [lab, dis] = a.saveColumns;
  assert.deepEqual([lab.name, dis.name], ['군집번호', '군집중심거리']);
  assert.equal(lab.values.length, 90);
  assert.equal(lab.values[4], null);
  assert.ok(lab.values.filter(v => v !== null).every(v => v >= 1 && v <= 3) && dis.values.filter(v => v !== null).every(v => v >= 0));
  assert.equal(lab.valueLabels, '1=군집 1, 2=군집 2, 3=군집 3');
  assert.ok(a.warnings.some(w => w.includes('1건을 제외')));
  assert.equal(K.runKmeans(dsFull, { vars: names, options: { k: '3' } }).saveColumns, undefined);
});

test('runKmeans: 입력 검증과 경고 (변수 수, k, 시작 횟수, 표준화 불가, 표준화 안 한 단위 차이, 작은 군집, 비수렴)', () => {
  const run = (options, d) => K.runKmeans(d || dsFull, { vars: names, options: options });
  assert.throws(() => K.runKmeans(dsFull, { vars: ['연령'], options: {} }), /2개 이상/);
  assert.throws(() => run({ k: '1' }), /2~10/);
  assert.throws(() => run({ k: '2.5' }), /정수/);
  assert.throws(() => run({ k: '11' }), /2~10/);
  assert.throws(() => run({ starts: '0' }), /시작 횟수/);
  assert.throws(() => run({ maxIter: '0' }), /최대 반복/);
  assert.throws(() => run({ seed: 'x' }), /시드/);
  assert.throws(() => run({ k: '3' }, ds([[1, 2, 3], [2, 3, 4], [3, 4, 5]])), /충분히 많아야/);
  assert.throws(() => run({}, ds([[1, 5, 3], [2, 5, 4], [3, 5, 1], [4, 5, 2], [5, 5, 9]])), /"소득"의 값이 모두 같아/);
  assert.ok(run({ standardize: 'none' }).warnings.some(w => w.includes('표준편차가 크게 다릅니다')));
  assert.ok(!run({ standardize: 'z' }).warnings.some(w => w.includes('표준편차가 크게 다릅니다')));
  assert.ok(run({ maxIter: '1', k: '4', starts: '1', seed: '3' }).warnings.some(w => w.includes('수렴하지 않았습니다')));
  // 이상치 하나가 군집이 되면 작은 군집 경고
  const out = ref.data.map(r => r.slice());
  out[0] = [200, 5000, 400];
  assert.ok(run({ k: '2' }, ds(out)).warnings.some(w => w.includes('가장 작은 군집')));
  // 사례가 많으면 엘보를 생략하고 안내
  const big = [];
  for (let i = 0; i < 5200; i++) big.push([i % 17, (i * 7) % 29, (i * 3) % 13 + (i % 2)]);
  const res = run({ k: '2', starts: '1', elbow: true }, ds(big));
  assert.ok(res.warnings.some(w => w.includes('엘보')) && !table(res, '군집 수 판단 보조'));
});

test('kmeansLloyd: 빈 군집은 가장 먼 사례로 채운다 (같은 점이 겹친 초기 중심)', () => {
  const X = [[0, 0], [0, 1], [10, 10], [10, 11], [20, 0]];
  const r = K.kmeansLloyd(X, [[0, 0], [0, 0], [10, 10]], 50);               // 두 중심이 같아 한 군집이 비게 된다
  const counts = [0, 0, 0];
  r.labels.forEach(l => counts[l]++);
  assert.ok(counts.every(c => c > 0));
  assert.ok(r.converged);
});
