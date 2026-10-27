const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const S = require('../src/Sample.js');
const K = require('../src/Compute.js');
const R = require('../src/Render.js');
const Ex = require('../src/Export.js');
const L = require('../src/Stats_Reliability.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'reliability.json'), 'utf8'));
const table = (res, caption) => res.blocks.find(b => b.type === 'table' && b.caption.startsWith(caption));
const near = (a, e, tol = 1e-12) => (a === null || e === null ? a === e : Math.abs(a - e) <= tol * Math.max(1, Math.abs(e)));

/* ---------- 정확 산술 기준값 대조 ---------- */

test('cronbachAlpha: 정확 산술 기준값 대조 (α, 표준화 α, 문항 통계, 문항 제외 시 α, 수정된 상관, 신뢰구간)', () => {
  assert.ok(ref.cases.length >= 10);
  for (const c of ref.cases) {
    const e = c.expected, r = L.cronbachAlpha(c.columns, 95);
    const where = c.name;
    assert.equal(r.n, e.n, where);
    assert.equal(r.nAll, e.nAll, where);
    assert.ok(near(r.alpha, e.alpha, 1e-11), `${where} α ${r.alpha} vs ${e.alpha}`);
    assert.ok(near(r.alphaStd, e.alphaStd, 1e-11), `${where} 표준화 α`);
    assert.ok(near(r.scaleMean, e.scaleMean, 1e-11) && near(r.scaleVar, e.scaleVar, 1e-11), `${where} 척도 통계`);
    e.items.forEach((x, i) => {
      const y = r.items[i];
      for (const key of ['mean', 'sd', 'scaleMeanIfDeleted', 'scaleVarIfDeleted', 'corrected', 'alphaIfDeleted']) {
        assert.ok(near(y[key], x[key], 1e-10), `${where} 문항${i + 1} ${key} ${y[key]} vs ${x[key]}`);
      }
    });
    if (e.interitem) {
      assert.ok(near(r.interitem.mean, e.interitem.mean, 1e-11) && near(r.interitem.min, e.interitem.min, 1e-11) && near(r.interitem.max, e.interitem.max, 1e-11), where);
      e.r.forEach((row, i) => row.forEach((v, j) => assert.ok(near(r.r[i][j], v, 1e-11), `${where} r[${i}][${j}]`)));
    } else {
      assert.equal(r.interitem, null, where);
    }
    if (e.ci) assert.ok(near(r.ci.lo, e.ci.lo, 1e-9) && near(r.ci.hi, e.ci.hi, 1e-9), `${where} 신뢰구간 ${r.ci && r.ci.lo} vs ${e.ci.lo}`);
    else assert.equal(r.ci, null, where);
  }
});

test('cronbachAlpha: 알려진 값과 성질', () => {
  // 완전히 같은 문항 → α = 1, 표준화 α = 1
  const same = [[1, 2, 3, 4, 5], [1, 2, 3, 4, 5], [1, 2, 3, 4, 5]];
  const r = L.cronbachAlpha(same);
  assert.ok(near(r.alpha, 1) && near(r.alphaStd, 1));
  assert.ok(near(r.ci.lo, 1) && near(r.ci.hi, 1));
  // 문항을 상수만큼 옮기거나 같은 배수를 곱해도 α 는 변하지 않는다
  const a = [[1, 3, 2, 5, 4, 2], [2, 3, 1, 5, 5, 3], [1, 4, 2, 4, 3, 1]];
  const base = L.cronbachAlpha(a).alpha;
  assert.ok(near(L.cronbachAlpha(a.map(c => c.map(v => v * 3 + 7))).alpha, base, 1e-12));
  // 문항 순서와 무관
  assert.ok(near(L.cronbachAlpha([a[2], a[0], a[1]]).alpha, base, 1e-12));
  // 계산 불가: 사례 1건, 문항 합계 분산 0
  assert.equal(L.cronbachAlpha([[1], [2]]).alpha, null);
  assert.equal(L.cronbachAlpha([[2, 2, 2], [3, 3, 3]]).alpha, null);
  // 결측: 목록별 제거 (한 문항이라도 비면 그 사례 제외)
  const m = L.cronbachAlpha([[1, 2, null, 4, 5], [2, 1, 3, null, 5], [1, 3, 2, 4, 4]]);
  assert.deepEqual([m.n, m.nAll], [3, 5]);
  // 문항이 2개면 문항 제외 시 α 는 null
  assert.equal(L.cronbachAlpha([[1, 2, 3, 4], [2, 1, 4, 3]]).items[0].alphaIfDeleted, null);
});

/* ---------- 결과 객체 ---------- */

const defs = [
  { name: 'a', label: '문항A', scale: 'ordinal', valueLabels: {}, missing: [99] },
  { name: 'b', label: '', scale: 'ordinal', valueLabels: {}, missing: [99] },
  { name: 'c', label: '', scale: 'ordinal', valueLabels: {}, missing: [99] },
  { name: 'd', label: '', scale: 'ordinal', valueLabels: {}, missing: [99] },
  { name: 't', label: '', scale: 'nominal', valueLabels: {}, missing: [] }
];
const rows = [['a', 'b', 'c', 'd', 't'],
  [5, 4, 5, 2, 'x'], [4, 4, 4, 3, 'x'], [2, 1, 2, 4, 'y'], [1, 2, 1, 5, 'y'], [3, 3, 4, 3, 'x'], [5, 5, 4, 1, 'y'], [2, 3, 2, 4, 'x'], [4, 99, 5, 2, 'y'], [3, 2, 3, 3, 'x'], [1, 1, 2, 5, 'y']];
const ds = () => Data.buildDataset(rows, defs);

test('runReliability: 표 구성(사례 처리 요약, 신뢰도 통계량, 문항 통계량, 문항-전체), 목록별 제외 경고', () => {
  const res = L.runReliability(ds(), { vars: ['a', 'b', 'c'], options: { itemTotal: true, itemDesc: true } });
  const proc = table(res, '사례 처리 요약');
  assert.deepEqual(proc.rows.map(r => r.n), [9, 1, 10]);
  assert.ok(res.warnings.some(w => w.includes('목록별 결측 제거로 1건')));
  assert.deepEqual(table(res, '신뢰도 통계량').columns.map(c => c.key), ['alpha', 'alphaStd', 'k']);
  assert.equal(table(res, '문항 통계량').rows.length, 3);
  const it = table(res, '문항-전체 통계량');
  assert.deepEqual(it.columns.map(c => c.key), ['name', 'meanDel', 'varDel', 'corrected', 'alphaDel']);
  assert.equal(it.rows[0].name, '문항A');                    // 라벨 사용
  assert.ok(!table(res, '문항 간 상관'));
  assert.ok(R.renderResultHtml(res).includes('크론바흐 α'));
  assert.ok(Ex.resultToGrid(res).rows.length > 10);
});

test('runReliability: 옵션 — 신뢰구간, 문항 간 상관, 척도 통계량, 끄기', () => {
  const res = L.runReliability(ds(), { vars: ['a', 'b', 'c'], options: { ci: true, interitem: true, scale: true, conf: '90' } });
  const st = table(res, '신뢰도 통계량');
  assert.deepEqual(st.columns.map(c => c.key), ['alpha', 'alphaStd', 'lo', 'hi', 'k']);
  assert.ok(st.rows[0].lo < st.rows[0].alpha && st.rows[0].alpha < st.rows[0].hi);
  assert.ok(st.columns[2].label.includes('90%'));
  assert.equal(table(res, '문항 간 상관행렬').rows.length, 3);
  assert.ok(table(res, '문항 간 상관 요약'));
  assert.ok(table(res, '척도 통계량'));
  assert.ok(!table(res, '문항 통계량') && !table(res, '문항-전체'));            // 옵션을 끄면 없음
  assert.throws(() => L.runReliability(ds(), { vars: ['a', 'b'], options: { conf: '120' } }), /신뢰수준/);
});

test('runReliability: 문항 제외 시 α 가 올라가는 문항은 굵게', () => {
  const res = L.runReliability(ds(), { vars: ['a', 'b', 'c', 'd'], options: { itemTotal: true } });
  const it = table(res, '문항-전체 통계량');
  const alpha = table(res, '신뢰도 통계량').rows[0].alpha;
  it.rows.forEach(r => assert.equal(!!r._bold, r.alphaDel > alpha, r.name));
  assert.ok(it.rows.some(r => r._bold), 'd 는 역방향 문항이라 빼면 α 가 올라가야 한다');
  assert.ok(it.notes.some(n => n.includes('굵게 표시한 문항')));
});

test('runReliability: 역문항을 뒤집지 않은 경우 음의 α 경고와 역문항 처리 안내', () => {
  const res = L.runReliability(ds(), { vars: ['a', 'b', 'c', 'd'], options: {} });
  const alpha = table(res, '신뢰도 통계량').rows[0].alpha;
  // d 가 역방향이지만 나머지 셋이 서로 잘 맞아 α 는 양수일 수 있다: 이 경우 음의 문항-전체 상관 경고
  const neg = res.warnings.find(w => w.includes('음수'));
  assert.ok(neg && neg.includes('역문항 처리'), res.warnings.join('\n'));
  // 문항을 모두 역방향으로 섞으면 α 자체가 음수
  const mixed = Data.buildDataset([['p', 'q', 'r', 's'], [1, 5, 1, 5], [2, 4, 2, 4], [3, 2, 3, 2], [4, 2, 4, 2], [5, 1, 5, 1], [2, 4, 3, 4]], []);
  const m = L.runReliability(mixed, { vars: ['p', 'q'], options: {} });
  assert.ok(table(m, '신뢰도 통계량').rows[0].alpha < 0);
  assert.ok(m.warnings.some(w => w.includes('α가 음수') && w.includes('역문항 처리')));
  assert.ok(alpha !== undefined);
});

test('runReliability: 입력 검증과 계산 불가 안내', () => {
  assert.throws(() => L.runReliability(ds(), { vars: ['a'], options: {} }), /2개 이상/);
  assert.throws(() => L.runReliability(ds(), { vars: ['a', 't'], options: {} }), /숫자가 아닌 값/);
  assert.throws(() => L.runReliability(ds(), { vars: [], options: {} }), /변수를 선택/);
  const few = Data.buildDataset([['x', 'y'], [1, ''], ['', 2], [3, '']], []);
  assert.throws(() => L.runReliability(few, { vars: ['x', 'y'], options: {} }), /응답한 사례가 0건/);
  // 상수 문항: 경고, 문항 간 상관행렬 불가
  const flat = Data.buildDataset([['x', 'y', 'z'], [1, 3, 2], [2, 3, 4], [3, 3, 3], [4, 3, 5]], []);
  const r = L.runReliability(flat, { vars: ['x', 'y', 'z'], options: { interitem: true, itemTotal: true } });
  assert.ok(r.warnings.some(w => w.includes('"y"') && w.includes('모든 값이 같아')));
  assert.ok(r.warnings.some(w => w.includes('문항 간 상관행렬을 만들 수 없습니다')));
  assert.ok(!table(r, '문항 간 상관행렬'));
});

/* ---------- 예제 데이터의 문항으로 ---------- */

test('예제 데이터: 문항1~5는 그대로는 α가 낮고 역문항을 뒤집으면 높다 (변수 계산 → 신뢰도 분석)', () => {
  const sample = S.sampleDataset_();
  const parsed = Data.parseVarDefs([['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(sample.varRows));
  const full = Data.buildDataset([sample.header].concat(sample.rows), parsed.defs);
  const raw = L.runReliability(full, { vars: ['문항1', '문항2', '문항3', '문항4', '문항5'], options: { itemTotal: true } });
  const rawAlpha = table(raw, '신뢰도 통계량').rows[0].alpha;
  assert.ok(raw.warnings.some(w => w.includes('음수') && w.includes('역문항 처리')), raw.warnings.join('\n'));

  const plan = K.planCompute_(full, { kind: 'reverse', vars: ['문항3', '문항5'], min: 1, max: 5, suffix: '_r' });
  const cols = ['문항1', '문항2', '문항4'].map(n => full.data[n]).concat(plan.columns.map(c => c.values));
  const fixed = L.cronbachAlpha(cols, 95);
  assert.ok(fixed.alpha > 0.7, `뒤집은 뒤 α ${fixed.alpha}`);
  assert.ok(fixed.alpha > rawAlpha + 0.3, `${rawAlpha} → ${fixed.alpha}`);
  assert.ok(fixed.items.every(it => it.corrected > 0.3), '뒤집은 뒤에는 모든 문항의 수정된 문항-전체 상관이 양수');
  assert.ok(fixed.ci.lo > 0.5 && fixed.ci.hi < 1);
});
