const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const Data = require('../src/Data.js');
const Rel = require('../src/Stats_Relation.js');
const R = require('../src/Render.js');

const ref = JSON.parse(fs.readFileSync(path.join(__dirname, 'reference', 'relation.json'), 'utf8'));

function relErr(a, e) {
  if (a === e) return 0;
  return Math.abs(a - e) / Math.max(Math.abs(e), 1e-300);
}
function close(a, e, tol, msg, absTol = 1e-10) {
  assert.ok(relErr(a, e) < tol || Math.abs(a - e) < absTol, `${msg}: ${a} vs ${e}`);
}
const ds = (values, defs) => Data.buildDataset(values, defs || []);
const tbl = (res, caption) => res.blocks.find(b => b.caption && b.caption.startsWith(caption));

/* ========== 상관분석 ========== */

test('pearsonR, spearmanR: scipy 대조 (r, p)', () => {
  for (const c of ref.pairs) {
    const p = Rel.pearsonR(c.x, c.y), s = Rel.spearmanR(c.x, c.y);
    close(p.r, c.pearson.r, 1e-9, `${c.name} pearson r`);
    close(p.p, c.pearson.p, 1e-8, `${c.name} pearson p`);
    close(s.r, c.spearman.r, 1e-9, `${c.name} spearman r`);
    close(s.p, c.spearman.p, 1e-8, `${c.name} spearman p`);
    assert.equal(p.df, c.x.length - 2);
  }
});

test('pearsonCI: scipy Fisher z 신뢰구간 대조', () => {
  for (const c of ref.pairs) {
    for (const conf of Object.keys(c.ci)) {
      const r = Rel.pearsonR(c.x, c.y);
      const ci = Rel.pearsonCI(r.r, r.n, Number(conf));
      close(ci.lo, c.ci[conf].lo, 1e-8, `${c.name} ${conf} lo`);
      close(ci.hi, c.ci[conf].hi, 1e-8, `${c.name} ${conf} hi`);
    }
  }
});

test('상관계수: 완전 상관, 상수, 표본이 작은 경우, 순위(동순위 평균)', () => {
  assert.equal(Rel.pearsonR([1, 2, 3, 4], [2, 4, 6, 8]).r, 1);
  assert.equal(Rel.pearsonR([1, 2, 3, 4], [2, 4, 6, 8]).p, 0);
  assert.equal(Rel.pearsonR([1, 2, 3, 4], [8, 6, 4, 2]).r, -1);
  assert.equal(Rel.pearsonR([1, 1, 1], [1, 2, 3]).r, null);
  assert.equal(Rel.pearsonR([1, 2], [1, 2]).r, null);
  assert.equal(Rel.pearsonCI(0.5, 3, 95).lo, null);   // n ≤ 3
  assert.equal(Rel.pearsonCI(1, 20, 95).lo, null);    // |r| = 1
  assert.deepEqual(Rel.rankAverage_([10, 30, 20, 20, 5]), [2, 5, 3.5, 3.5, 1]);
  assert.equal(Rel.spearmanR([1, 2, 3, 4, 5], [1, 4, 9, 16, 25]).r, 1);   // 단조 증가면 스피어만 1
});

function matrixData() {
  const m = ref.matrix;
  const n = m.data[m.names[0]].length;
  const rows = [m.names.slice()];
  for (let i = 0; i < n; i++) rows.push(m.names.map(k => (m.data[k][i] === null ? '' : m.data[k][i])));
  return ds(rows);
}

test('runCorrelation: 쌍별 결측 행렬이 쌍마다 scipy와 일치 (r, p, N)', () => {
  const m = ref.matrix;
  for (const method of ['pearson', 'spearman']) {
    const res = Rel.runCorrelation(matrixData(), { vars: m.names, options: { pearson: method === 'pearson', spearman: method === 'spearman' } });
    const t = res.blocks[0];
    assert.equal(res.blocks.length, 1);
    assert.equal(t.rows.length, m.names.length * 3);
    m.names.forEach((a, i) => {
      m.names.forEach((b, j) => {
        const rRow = t.rows[i * 3], pRow = t.rows[i * 3 + 1], nRow = t.rows[i * 3 + 2];
        if (i === j) {
          assert.equal(rRow['c' + j], 1);
          assert.equal(pRow['c' + j], null);
          return;
        }
        const key = i < j ? `${a},${b}` : `${b},${a}`;
        const e = m.pairwise[method][key];
        close(rRow['c' + j], e.r, 1e-9, `${method} ${key} r`);
        close(pRow['c' + j], e.p, 1e-8, `${method} ${key} p`);
        assert.equal(nRow['c' + j], e.n, `${method} ${key} N`);
      });
    });
  }
});

test('runCorrelation: 목록별 결측 제거는 모든 쌍에 같은 N', () => {
  const m = ref.matrix;
  const res = Rel.runCorrelation(matrixData(), { vars: m.names, options: { pearson: true, listwise: true } });
  const t = res.blocks[0];
  m.names.forEach((a, i) => {
    m.names.forEach((b, j) => {
      if (i === j) return;
      const key = i < j ? `${a},${b}` : `${b},${a}`;
      const e = m.listwise.pearson[key];
      close(t.rows[i * 3]['c' + j], e.r, 1e-9, `listwise ${key} r`);
      close(t.rows[i * 3 + 1]['c' + j], e.p, 1e-8, `listwise ${key} p`);
      assert.equal(t.rows[i * 3 + 2]['c' + j], e.n);
    });
  });
  assert.ok(res.warnings.some(w => w.includes('목록별 결측 제거')));
});

test('runCorrelation: 표 구조(행 서식, 열), 두 종류, 신뢰구간 표', () => {
  const m = ref.matrix;
  const res = Rel.runCorrelation(matrixData(), { vars: ['a', 'b', 'd'], options: { pearson: true, spearman: true, ci: true, conf: 99 } });
  assert.deepEqual(res.blocks.map(b => b.caption), ['피어슨 상관 행렬', '스피어만 순위상관 행렬', '피어슨 상관계수의 99% 신뢰구간 (Fisher z 변환)']);
  const t = res.blocks[0];
  assert.deepEqual(t.columns.map(c => c.key), ['name', 'stat', 'c0', 'c1', 'c2']);
  assert.deepEqual(t.rows.slice(0, 3).map(r => r._fmt.type), ['fixed', 'pvalue', 'int']);
  assert.equal(t.rows[0]._kind, 'group');
  assert.deepEqual(t.rows.slice(0, 3).map(r => r.stat), ['Pearson 상관계수', '유의확률(양측)', 'N']);
  const ci = res.blocks[2];
  assert.equal(ci.rows.length, 3);                       // 쌍 3개
  assert.ok(ci.rows.every(r => r.lo < r.r && r.r < r.hi));
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('Pearson 상관계수') && html.includes('유의확률(양측)'));
});

test('runCorrelation: 입력 오류와 경고', () => {
  const d = ds([['a', 'b', 'c'], [1, 2, 'x'], [2, 3, 'y'], [3, 5, 'z'], [4, 4, 'w']]);
  assert.throws(() => Rel.runCorrelation(d, { vars: ['a'], options: { pearson: true } }), /2개 이상/);
  assert.throws(() => Rel.runCorrelation(d, { vars: ['a', 'b'], options: {} }), /하나 이상/);
  assert.throws(() => Rel.runCorrelation(d, { vars: ['a', 'c'], options: { pearson: true } }), /숫자가 아닌 값/);
  assert.throws(() => Rel.runCorrelation(d, { vars: ['a', 'b'], options: { pearson: true, conf: 100 } }), /신뢰수준/);
  const flat = ds([['a', 'b'], [1, 5], [2, 5], [3, 5]]);
  const res = Rel.runCorrelation(flat, { vars: ['a', 'b'], options: { pearson: true } });
  assert.equal(res.warnings.length, 1);
  assert.equal(Number.isNaN(res.blocks[0].rows[0].c1), true);          // 계산 불가는 NaN → '.'
  assert.ok(R.renderResultHtml(res).includes('>.<'));
});

/* ========== 교차분석 ========== */

/** 빈도 행렬을 (행 코드, 열 코드) 사례로 펼친 데이터셋. 코드는 1부터. */
function tableDataset(t, defs) {
  const rows = [['r', 'c']];
  t.forEach((row, i) => row.forEach((n, j) => { for (let k = 0; k < n; k++) rows.push([i + 1, j + 1]); }));
  return ds(rows, defs);
}
const ALL_CROSS = { expected: true, rowPct: true, colPct: true, totalPct: true, residual: true, effect: true };
const testRow = (res, label) => tbl(res, '카이제곱 검정').rows.find(r => r.test === label);

test('runCrosstab: scipy 대조 (Pearson, 우도비, Yates, Fisher, 선형 대 선형, 효과크기)', () => {
  for (const c of ref.tables) {
    const res = Rel.runCrosstab(tableDataset(c.table), { rowvar: ['r'], colvar: ['c'], options: ALL_CROSS });
    const m = c.name;
    const pe = testRow(res, 'Pearson 카이제곱');
    close(pe.value, c.chi2, 1e-10, `${m} chi2`);
    assert.equal(pe.df, c.dof);
    close(pe.p, c.p, 1e-8, `${m} p`);
    const lr = testRow(res, '우도비');
    close(lr.value, c.lr, 1e-10, `${m} LR`);
    close(lr.p, c.lr_p, 1e-8, `${m} LR p`);
    const lin = testRow(res, '선형 대 선형 연관');
    close(lin.value, c.lbl.value, 1e-9, `${m} 선형 대 선형`);
    close(lin.p, c.lbl.p, 1e-8, `${m} 선형 대 선형 p`);
    assert.equal(testRow(res, '유효 케이스 수').value, c.N);
    if (c.yates) {
      close(testRow(res, '연속성 수정(Yates)').value, c.yates.value, 1e-10, `${m} Yates`);
      close(testRow(res, '연속성 수정(Yates)').p, c.yates.p, 1e-8, `${m} Yates p`);
      const fe = testRow(res, 'Fisher의 정확검정');
      close(fe.p2, c.fisher.two, 1e-9, `${m} Fisher 양측`);
      close(fe.p1, c.fisher.one, 1e-9, `${m} Fisher 단측`);
    } else {
      assert.equal(testRow(res, '연속성 수정(Yates)'), undefined);
      assert.equal(testRow(res, 'Fisher의 정확검정'), undefined);
    }
    const eff = tbl(res, '대칭적 측도').rows;
    close(eff.find(r => r.name === 'Cramér의 V').value, c.cramer, 1e-10, `${m} V`);
    close(eff.find(r => r.name === '분할계수').value, c.contingency, 1e-10, `${m} C`);
    if (c.phi !== undefined) close(eff.find(r => r.name === '파이').value, c.phi, 1e-10, `${m} phi`);
    else assert.equal(eff.find(r => r.name === '파이'), undefined);
  }
});

test('runCrosstab: 교차표의 빈도, 기대빈도, 퍼센트, 조정 잔차', () => {
  for (const c of ref.tables) {
    const res = Rel.runCrosstab(tableDataset(c.table), { rowvar: ['r'], colvar: ['c'], options: ALL_CROSS });
    const t = tbl(res, '교차표');
    const nr = c.table.length, nc = c.table[0].length, N = c.N;
    const R = c.table.map(r => r.reduce((a, b) => a + b, 0));
    const C = c.table[0].map((_, j) => c.table.reduce((a, r) => a + r[j], 0));
    const per = 6; // 빈도, 기대빈도, 행 %, 열 %, 전체 %, 잔차
    assert.equal(t.rows.length, (nr + 1) * per);
    for (let i = 0; i < nr; i++) {
      const g = t.rows.slice(i * per, (i + 1) * per);
      assert.deepEqual(g.map(r => r.stat), ['빈도', '기대빈도', '행 %', '열 %', '전체 %', '조정된 표준화 잔차']);
      assert.equal(g[0]._kind, 'group');
      for (let j = 0; j < nc; j++) {
        assert.equal(g[0]['c' + j], c.table[i][j]);
        close(g[1]['c' + j], c.expected[i][j], 1e-12, 'expected');
        close(g[2]['c' + j], c.table[i][j] / R[i] * 100, 1e-12, 'row%');
        close(g[3]['c' + j], c.table[i][j] / C[j] * 100, 1e-12, 'col%');
        close(g[4]['c' + j], c.table[i][j] / N * 100, 1e-12, 'tot%');
        close(g[5]['c' + j], c.adj_resid[i][j], 1e-10, 'adj resid');
      }
      assert.equal(g[0].total, R[i]);
      assert.equal(g[2].total, 100);
      close(g[3].total, R[i] / N * 100, 1e-12, 'row total col%');
    }
    const tot = t.rows.slice(nr * per);
    assert.equal(tot[0]._kind, 'total');
    assert.equal(tot[0].total, N);
    assert.deepEqual(C, C.map((_, j) => tot[0]['c' + j]));
    assert.equal(tot[5].c0, null);
  }
});

test('runCrosstab: 옵션을 끄면 빈도만, 값 레이블과 열 제목', () => {
  const defs = Data.parseVarDefs([['h'], ['r', '성별', '명목', '1=남, 2=여', ''], ['c', '합격', '명목', '1=불합격, 2=합격', '']]).defs;
  const res = Rel.runCrosstab(tableDataset([[12, 5], [7, 16]], defs), { rowvar: ['r'], colvar: ['c'], options: {} });
  const t = tbl(res, '교차표');
  assert.deepEqual(t.rows.map(r => r.stat), ['빈도', '빈도', '빈도']);
  assert.deepEqual(t.rows.map(r => r.name), ['남', '여', '합계']);
  assert.deepEqual(t.columns.map(c => c.label), ['성별', '', '불합격', '합격', '합계']);
  assert.ok(t.caption.includes('성별 × 합격'));
  assert.equal(tbl(res, '대칭적 측도'), undefined);
  const html = R.renderResultHtml(res);
  assert.ok(html.includes('불합격') && html.includes('Fisher의 정확검정'));
});

test('fisherExact2x2: 알려진 값과 대칭성', () => {
  // 차 시음 실험: [[3, 1], [1, 3]] → 양측 p = 0.4857142857..., 단측 p = 0.2428571428...
  const f = Rel.fisherExact2x2(3, 1, 1, 3);
  close(f.two, 0.4857142857142857, 1e-12, 'two');
  close(f.one, 0.24285714285714285, 1e-12, 'one');
  const g = Rel.fisherExact2x2(1, 3, 3, 1);
  close(g.two, f.two, 1e-12, 'swap two');
  const z = Rel.fisherExact2x2(5, 0, 0, 5);
  close(z.two, 2 / 252, 1e-12, 'perfect association');
});

test('runCrosstab: 문자 범주는 선형 대 선형을 만들지 않고, 기대빈도 경고', () => {
  const d = ds([['r', 'c'], ['a', 'x'], ['a', 'x'], ['b', 'y'], ['b', 'x'], ['a', 'y'], ['b', 'y']]);
  const res = Rel.runCrosstab(d, { rowvar: ['r'], colvar: ['c'], options: {} });
  assert.equal(testRow(res, '선형 대 선형 연관'), undefined);
  assert.ok(res.warnings.some(w => w.includes('기대빈도가 5 미만')));
  assert.ok(tbl(res, '카이제곱 검정').notes[0].includes('최소 기대빈도'));
});

test('runCrosstab: 결측 사례 제외 안내와 입력 오류', () => {
  const d = ds([['r', 'c'], [1, 1], [1, 2], [2, 1], [2, 2], ['', 1], [1, ''], [2, 2]]);
  const res = Rel.runCrosstab(d, { rowvar: ['r'], colvar: ['c'], options: {} });
  assert.ok(tbl(res, '교차표').notes[0].includes('2건 제외'));
  assert.equal(testRow(res, '유효 케이스 수').value, 5);
  assert.throws(() => Rel.runCrosstab(d, { rowvar: ['r'], colvar: ['r'] }), /같은 변수/);
  assert.throws(() => Rel.runCrosstab(d, { rowvar: [], colvar: ['c'] }), /변수를 선택/);
  assert.throws(() => Rel.runCrosstab(d, { rowvar: ['r', 'c'], colvar: ['c'] }), /하나씩/);
  assert.throws(() => Rel.runCrosstab(ds([['r', 'c'], [1, 1], [1, 2], [1, 1]]), { rowvar: ['r'], colvar: ['c'] }), /2개 이상/);
  assert.throws(() => Rel.runCrosstab(ds([['r', 'c'], [1, ''], ['', 1]]), { rowvar: ['r'], colvar: ['c'] }), /유효한 사례가 없습니다/);
  const many = [['r', 'c']];
  for (let i = 0; i < 35; i++) many.push([i, i % 2], [i, (i + 1) % 2]);
  assert.throws(() => Rel.runCrosstab(ds(many), { rowvar: ['r'], colvar: ['c'] }), /너무 많습니다/);
});
