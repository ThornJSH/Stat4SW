const test = require('node:test');
const assert = require('node:assert/strict');
const D = require('../src/Data.js');

test('값 레이블 파싱', () => {
  assert.deepEqual(D.parseValueLabels_('1=남, 2=여'), { '1': '남', '2': '여' });
  assert.deepEqual(D.parseValueLabels_('1=매우 낮음;5=매우 높음'), { '1': '매우 낮음', '5': '매우 높음' });
  assert.deepEqual(D.parseValueLabels_('1=a\n2=b=c'), { '1': 'a', '2': 'b=c' });
  assert.deepEqual(D.parseValueLabels_(''), {});
  assert.deepEqual(D.parseValueLabels_(null), {});
});

test('결측값 코드 파싱', () => {
  assert.deepEqual(D.parseMissingCodes_('99, 999'), [99, 999]);
  assert.deepEqual(D.parseMissingCodes_(99), [99]);
  assert.deepEqual(D.parseMissingCodes_('NA'), ['NA']);
  assert.deepEqual(D.parseMissingCodes_(''), []);
});

test('변수정의 파싱: 척도 매핑, 중복·알 수 없는 척도 경고', () => {
  const rows = [
    ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'],
    ['성별', '성별', '명목', '1=남, 2=여', ''],
    ['만족도', '만족도', '서열', '', 99],
    ['', '빈 이름', '연속', '', ''],
    ['성별', '중복', '연속', '', ''],
    ['키', '신장', '연속?', '', '']
  ];
  const { defs, warnings } = D.parseVarDefs(rows);
  assert.deepEqual(defs.map(d => d.name), ['성별', '만족도', '키']);
  assert.equal(defs[0].scale, 'nominal');
  assert.equal(defs[1].scale, 'ordinal');
  assert.deepEqual(defs[1].missing, [99]);
  assert.equal(defs[2].scale, null);
  assert.equal(warnings.length, 2);
});

test('데이터셋: 결측 코드·빈 셀 → null, 완전 빈 행 제외', () => {
  const defs = D.parseVarDefs([
    ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'],
    ['성별', '', '명목', '1=남, 2=여', ''],
    ['점수', '', '연속', '', 999]
  ]).defs;
  const values = [
    ['성별', '점수', '메모'],
    [1, 80, 'a'],
    [2, 999, ''],
    ['', '', ''],
    [1, '', 'b'],
    [2, '75.5', 'c']
  ];
  const ds = D.buildDataset(values, defs);
  assert.equal(ds.n, 4);
  assert.deepEqual(ds.data['성별'], [1, 2, 1, 2]);
  assert.deepEqual(ds.data['점수'], [80, null, null, 75.5]);
  assert.equal(ds.vars.find(v => v.name === '메모').inferred, true);
  assert.equal(ds.vars.find(v => v.name === '메모').scale, 'nominal');
  assert.equal(ds.vars.find(v => v.name === '점수').inferred, false);
});

test('데이터셋: 정의 없는 변수의 척도 자동 판정', () => {
  const ds = D.buildDataset([['a', 'b'], [1, 'x'], [2, 'y']], []);
  assert.equal(ds.vars[0].scale, 'scale');
  assert.equal(ds.vars[1].scale, 'nominal');
});

test('데이터셋: 연속 변수의 문자는 결측 처리하고 경고', () => {
  const defs = D.parseVarDefs([['h'], ['x', '', '연속']]).defs;
  const ds = D.buildDataset([['x'], [1], ['abc'], [3]], defs);
  assert.deepEqual(ds.data['x'], [1, null, 3]);
  assert.equal(ds.warnings.length, 1);
});

test('데이터셋: 데이터에 없는 정의와 중복 헤더 경고', () => {
  const defs = D.parseVarDefs([['h'], ['없는변수', '', '연속']]).defs;
  const ds = D.buildDataset([['a', 'a'], [1, 2]], defs);
  assert.equal(ds.vars.length, 1);
  assert.equal(ds.warnings.length, 2);
});

test('목록별 결측 제거', () => {
  const ds = D.buildDataset([['x', 'y', 'z'], [1, 1, 1], [2, '', 2], [3, 3, 3], ['', 4, 4]], []);
  assert.deepEqual(D.listwiseIndices(ds, ['x', 'y']), [0, 2]);
  const sel = D.selectComplete(ds, ['x', 'y', 'z']);
  assert.equal(sel.n, 2);
  assert.equal(sel.dropped, 2);
  assert.deepEqual(sel.columns['z'], [1, 3]);
});
