const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const S = require('../src/Stats_Descriptive.js');
const R = require('../src/Render.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'descriptive.json'), 'utf8'));

function relErr(a, e) {
  if (a === e) return 0;
  return Math.abs(a - e) / Math.max(Math.abs(e), 1e-300);
}

/** 변수정의 없이 헤더+행으로 데이터셋을 만든다 */
function ds(values, defs) {
  return Data.buildDataset(values, defs || []);
}

function findTable(result, caption) {
  return result.blocks.find(b => b.type === 'table' && b.caption === caption);
}

/* ---------- describeColumn: 정확 산술 기준값 대조 ---------- */

test('describeColumn: 정확 산술 기준값(scipy와 교차 확인됨) 대조', () => {
  const keys = ['n', 'min', 'max', 'range', 'sum', 'mean', 'variance', 'sd', 'seMean', 'skewness', 'kurtosis'];
  let worst = 0;
  for (const c of ref.cases) {
    const r = S.describeColumn(c.values);
    for (const k of keys) {
      if (c[k] === null) {
        assert.equal(r[k], null, `${c.name} ${k}`);
        continue;
      }
      const e = relErr(r[k], c[k]);
      worst = Math.max(worst, e);
      assert.ok(e < 1e-12, `${c.name} ${k}: ${r[k]} vs ${c[k]} (상대오차 ${e.toExponential(2)})`);
    }
  }
  assert.ok(worst < 1e-12);
});

test('describeColumn: 왜도·첨도 표준오차는 알려진 표 값과 일치 (n=10, 100)', () => {
  const fill = n => Array.from({ length: n }, (_, i) => i + 1);
  const r10 = S.describeColumn(fill(10));
  assert.equal(r10.seSkewness.toFixed(3), '0.687');
  assert.equal(r10.seKurtosis.toFixed(3), '1.334');
  const r100 = S.describeColumn(fill(100));
  assert.equal(r100.seSkewness.toFixed(3), '0.241');
  assert.equal(r100.seKurtosis.toFixed(3), '0.478');
});

test('describeColumn: 표본이 작거나 분산이 0이면 계산 불가 값은 null', () => {
  assert.equal(S.describeColumn([]).mean, null);
  const one = S.describeColumn([5]);
  assert.equal(one.mean, 5);
  assert.equal(one.sd, null);
  const flat = S.describeColumn([3, 3, 3, 3, 3]);
  assert.equal(flat.sd, 0);
  assert.equal(flat.skewness, null);
  assert.equal(flat.kurtosis, null);
  assert.equal(S.describeColumn([1, 2, 3]).kurtosis, null); // n ≤ 3
});

/* ---------- getVars ---------- */

test('getVars: 빈 선택, 없는 변수, 숫자 아닌 값, 중복', () => {
  const d = ds([['a', 'b'], [1, 'x'], [2, 'y']]);
  assert.throws(() => Data.getVars(d, []), /변수를 선택/);
  assert.throws(() => Data.getVars(d, ['zzz']), /찾을 수 없습니다/);
  assert.throws(() => Data.getVars(d, ['b'], { numeric: true }), /숫자가 아닌 값/);
  assert.deepEqual(Data.getVars(d, ['a', 'a', 'b']).map(v => v.name), ['a', 'b']);
});

/* ---------- 빈도분석 ---------- */

test('빈도분석: 값 레이블, 퍼센트, 누적, 결측 행', () => {
  const defs = Data.parseVarDefs([
    ['h'],
    ['성별', '성별', '명목', '1=남, 2=여', ''],
    ['점수', '', '연속', '', 99]
  ]).defs;
  const d = ds([
    ['성별', '점수'],
    [1, 10], [2, 20], [1, 99], [1, 10], [2, ''], [2, 30], [1, 20], [2, 10]
  ], defs);
  const res = S.runFrequencies(d, { vars: ['성별', '점수'] });

  const stat = findTable(res, '통계량');
  assert.deepEqual(stat.rows.map(r => [r.name, r.valid, r.missing]), [['성별', 8, 0], ['점수', 6, 2]]);

  const sex = findTable(res, '성별');
  assert.deepEqual(sex.rows.map(r => r.value), ['남', '여', '']);
  assert.deepEqual(sex.rows.map(r => r.freq), [4, 4, 8]);
  assert.equal(sex.rows[0].pct, 50);
  assert.equal(sex.rows[1].cumPct, 100);
  assert.equal(sex.rows[2]._kind, 'total');
  assert.equal(sex.rows.length, 3); // 결측이 없으면 소계·결측 행이 없다

  const score = findTable(res, '점수');
  // 유효 값: 10(3건), 20(2건), 30(1건), 결측 2건, 전체 8건
  assert.deepEqual(score.rows.map(r => r.freq), [3, 2, 1, 6, 2, 8]);
  assert.deepEqual(score.rows.map(r => r.value), ['10', '20', '30', '합계', '결측', '']);
  assert.equal(score.rows[0].pct, 3 / 8 * 100);
  assert.equal(score.rows[0].validPct, 50);
  assert.equal(score.rows[2].cumPct, 100);
  assert.equal(score.rows[3].validPct, 100);
  assert.equal(score.rows[4].pct, 25);
  assert.equal(score.rows[5].pct, 100);
});

test('빈도분석: 문자 값은 사전 순, 숫자가 앞', () => {
  const d = ds([['v'], ['b'], ['a'], [2], ['b'], [1]], Data.parseVarDefs([['h'], ['v', '', '명목']]).defs);
  const t = findTable(S.runFrequencies(d, { vars: ['v'] }), 'v');
  assert.deepEqual(t.rows.slice(0, 4).map(r => r.value), ['1', '2', 'a', 'b']);
});

test('빈도분석: 고유값이 100개를 넘으면 경고', () => {
  const rows = [['x']];
  for (let i = 0; i < 150; i++) rows.push([i]);
  const res = S.runFrequencies(ds(rows), { vars: ['x'] });
  assert.equal(res.warnings.length, 1);
});

test('빈도분석: 데이터가 없으면 오류', () => {
  assert.throws(() => S.runFrequencies(ds([['x']]), { vars: ['x'] }), /데이터가 없습니다/);
});

/* ---------- 기술통계 ---------- */

test('기술통계: 변수별 N, 마지막 행은 목록별 N', () => {
  const d = ds([['a', 'b'], [1, 10], [2, ''], [3, 30], ['', 40], [5, 50]]);
  const res = S.runDescriptives(d, { vars: ['a', 'b'] });
  const t = findTable(res, '기술통계량');
  assert.deepEqual(t.rows.map(r => r.n), [4, 4, 3]); // a: 1,2,3,5  b: 10,30,40,50  목록별: 3건
  assert.equal(t.rows[0].mean, 2.75);
  assert.equal(t.rows[2].name, '유효 N(목록별)');
  assert.equal(t.rows[2]._kind, 'total');
});

test('기술통계: 목록별 결측 제거 옵션은 모든 변수에 같은 사례만 쓰고 경고', () => {
  const d = ds([['a', 'b'], [1, 10], [2, ''], [3, 30], ['', 40], [5, 50]]);
  const res = S.runDescriptives(d, { vars: ['a', 'b'], options: { listwise: true } });
  const t = findTable(res, '기술통계량');
  assert.deepEqual(t.rows.map(r => r.n), [3, 3, 3]);
  assert.equal(t.rows[0].mean, 3); // (1+3+5)/3
  assert.equal(res.warnings.length, 1);
});

test('기술통계: 옵션에 따라 열이 추가된다', () => {
  const d = ds([['a'], [1], [2], [3], [4], [10]]);
  const base = findTable(S.runDescriptives(d, { vars: ['a'] }), '기술통계량').columns.map(c => c.key);
  assert.deepEqual(base, ['name', 'n', 'min', 'max', 'mean', 'sd']);
  const all = findTable(S.runDescriptives(d, { vars: ['a'], options: { range: true, sum: true, seMean: true, variance: true, skewness: true, kurtosis: true } }), '기술통계량').columns.map(c => c.key);
  assert.deepEqual(all, ['name', 'n', 'range', 'min', 'max', 'sum', 'mean', 'seMean', 'sd', 'variance', 'skewness', 'seSkewness', 'kurtosis', 'seKurtosis']);
});

test('기술통계: 숫자가 아닌 변수는 거부, 유효값이 없으면 경고', () => {
  assert.throws(() => S.runDescriptives(ds([['t'], ['x'], ['y']]), { vars: ['t'] }), /숫자가 아닌 값/);
  const empty = ds([['a', 'b'], [1, ''], [2, '']]);
  const res = S.runDescriptives(empty, { vars: ['a', 'b'] });
  assert.ok(res.warnings.some(w => w.includes('b')));
});

/* ---------- 렌더러 ---------- */

test('렌더러: 숫자 형식', () => {
  const f = R.formatCell_;
  assert.equal(f(3.14159, { type: 'fixed', digits: 3 }), '3.142');
  assert.equal(f(-0.0001, { type: 'fixed', digits: 3 }), '0.000');
  assert.equal(f(2.5, { type: 'trim', digits: 3 }), '2.5');
  assert.equal(f(20, { type: 'trim', digits: 3 }), '20');
  assert.equal(f(7.4, { type: 'int' }), '7');
  assert.equal(f(0.0004, { type: 'pvalue' }), '<.001');
  assert.equal(f(0.0432, { type: 'pvalue' }), '.043');
  assert.equal(f(1, { type: 'pvalue' }), '1.000');
  assert.equal(f(null, { type: 'fixed' }), '');
  assert.equal(f(NaN, { type: 'fixed' }), '.');
  assert.equal(f('텍스트', { type: 'fixed' }), '텍스트');
});

test('렌더러: HTML 이스케이프와 구조', () => {
  const html = R.renderResultHtml({
    analysis: 'x',
    title: '<b>제목</b>',
    blocks: [
      { type: 'table', caption: 'A&B', columns: [{ key: 'k', label: '열<1>' }, { key: 'n', label: 'N', format: { type: 'int' } }], rows: [{ k: '<script>', n: 3 }, { _kind: 'total', k: '합계', n: 3 }], notes: ['각주 "1"'] },
      { type: 'text', text: '<i>문단</i>' }
    ],
    warnings: ['경고 & 주의']
  });
  assert.ok(!html.includes('<script>') && !html.includes('<b>') && !html.includes('<i>'));
  assert.ok(html.includes('&lt;script&gt;') && html.includes('A&amp;B') && html.includes('경고 &amp; 주의'));
  assert.ok(html.includes('<tr class="total">'));
  assert.ok(html.includes('<td class="right">3</td>'));
  assert.ok(html.includes('<td class="left">&lt;script&gt;</td>'));
  assert.ok(html.includes('class="note"') && html.includes('class="warn"'));
});

test('렌더러: 실제 분석 결과를 끝까지 렌더링', () => {
  const d = ds([['a'], [1], [2], [3], [4], [10]]);
  const html = R.renderResultHtml(S.runDescriptives(d, { vars: ['a'] }));
  assert.ok(html.includes('기술통계량') && html.includes('유효 N(목록별)') && html.includes('4.000'));
});

test('렌더러: 행 단위 서식(_fmt)은 그 행의 숫자 셀에 열 형식 대신 적용된다', () => {
  const html = R.renderResultHtml({
    analysis: 'x', title: 't',
    blocks: [{
      type: 'table', caption: 'c',
      columns: [{ key: 'name', label: '' }, { key: 'v', label: 'v', format: { type: 'fixed', digits: 3 } }],
      rows: [
        { name: '상관', v: 0.5, _fmt: { type: 'fixed', digits: 3 } },
        { name: '유의확률', v: 0.0004, _fmt: { type: 'pvalue' } },
        { name: 'N', v: 12, _fmt: { type: 'int' }, _kind: 'group' },
        { name: '기본', v: 2 }
      ]
    }],
    warnings: []
  });
  assert.ok(html.includes('<td class="right">0.500</td>'));
  assert.ok(html.includes('<td class="right">&lt;.001</td>'));
  assert.ok(html.includes('<td class="right">12</td>'));
  assert.ok(html.includes('<td class="right">2.000</td>'));
  assert.ok(html.includes('<tr class="group">'));
});

test('빈도분석: 고유한 값이 500개를 넘으면 빈도표를 만들지 않고 안내, 다른 변수는 계속', () => {
  const rows = [['x', 'g']];
  for (let i = 0; i < 600; i++) rows.push([i, i % 3]);
  const res = S.runFrequencies(ds(rows), { vars: ['x', 'g'], options: { barChart: true } });
  assert.ok(res.warnings.some(w => w.includes('"x"') && w.includes('600개') && w.includes('빈도표를 만들지 않았습니다')));
  const captions = res.blocks.map(b => b.caption);
  assert.ok(!captions.includes('x') && !captions.includes('막대그래프: x'));
  assert.ok(captions.includes('g') && captions.includes('막대그래프: g'));
  assert.equal(res.blocks.find(b => b.caption === '통계량').rows.length, 2);       // 통계량 표에는 두 변수 모두
});
