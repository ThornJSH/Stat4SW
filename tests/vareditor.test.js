const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const Data = require('../src/Data.js');
const V = require('../src/VarEditor.js');

const SRC = path.join(__dirname, '..', 'src');
const json = (x) => JSON.parse(JSON.stringify(x));

/* ========== 순수 함수 ========== */

test('formatValueLabels_: "값=레이블, ..." 형식, 빈 레이블은 건너뜀, 순서 유지', () => {
  const r = V.formatValueLabels_([{ key: '1', label: '남' }, { key: '2', label: ' 여 ' }, { key: '3', label: '' }, { key: '9', label: '   ' }]);
  assert.equal(r.text, '1=남, 2=여');
  assert.deepEqual(r.warnings, []);
  assert.equal(V.formatValueLabels_([]).text, '');
  assert.equal(V.formatValueLabels_([{ key: 'a', label: '가' }, { key: 'b', label: '나' }]).text, 'a=가, b=나');
});

test('formatValueLabels_: 쉼표·세미콜론·줄바꿈은 공백으로, 값에 구분 문자가 있으면 건너뛰고 경고', () => {
  const r = V.formatValueLabels_([
    { key: '1', label: '매우, 낮음' }, { key: '2', label: '보통;중간\n정도' }, { key: 'a,b', label: '못 저장' }, { key: 'x=y', label: '못 저장' }, { key: '5', label: '등호=포함 가능' }
  ]);
  assert.equal(r.text, '1=매우 낮음, 2=보통 중간 정도, 5=등호=포함 가능');
  assert.equal(r.warnings.length, 4);
  assert.ok(r.warnings.some(w => w.includes('"a,b"')) && r.warnings.some(w => w.includes('"x=y"')));
});

test('formatValueLabels_ → parseValueLabels_ 왕복: 저장한 문자열을 분석이 같은 값 레이블로 읽는다', () => {
  const pairs = [{ key: '1', label: '매우 낮음' }, { key: '2.5', label: '중간 =값' }, { key: '남', label: 'male' }, { key: '-1', label: '해당 없음' }];
  const parsed = Data.parseValueLabels_(V.formatValueLabels_(pairs).text);
  assert.deepEqual(parsed, { '1': '매우 낮음', '2.5': '중간 =값', '남': 'male', '-1': '해당 없음' });
});

test('formatMissingCodes_: 구분자 정리, 중복 제거, 문자 코드도 허용', () => {
  assert.equal(V.formatMissingCodes_('99 ;999,99'), '99, 999');
  assert.equal(V.formatMissingCodes_(''), '');
  assert.equal(V.formatMissingCodes_(null), '');
  assert.equal(V.formatMissingCodes_(99), '99');
  assert.equal(V.formatMissingCodes_('NA, -1'), 'NA, -1');
  assert.deepEqual(Data.parseMissingCodes_(V.formatMissingCodes_('99, 999')), [99, 999]);       // 저장한 문자열을 파서가 같은 코드로 읽는다
});

function dsOf(rows, defRows) {
  const defs = defRows ? Data.parseVarDefs([['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(defRows)).defs : [];
  return Data.buildDataset(rows, defs);
}

test('buildEditorData_: 데이터의 값과 사례 수, 기존 레이블, 데이터에 없는 레이블, 척도·결측값', () => {
  const ds = dsOf([['s', 'a'], [2, 1], [1, 1], [2, 2], [1, 99], [3, ''], [2, 1]],
    [['s', '성별', '명목', '1=남, 2=여, 7=기타', ''], ['a', '', '', '', 99]]);
  const d = V.buildEditorData_(ds, 's');
  assert.deepEqual(d.values.map(v => [v.key, v.count, v.label, v.inData]), [
    ['1', 2, '남', true], ['2', 3, '여', true], ['3', 1, '', true], ['7', 0, '기타', false]    // 데이터에 없는 레이블 7 도 보존
  ]);
  assert.equal(d.label, '성별');
  assert.equal(d.scale, '명목');
  assert.equal(d.inferred, false);
  assert.equal(d.n, 6);
  assert.equal(d.truncated, false);
  const a = V.buildEditorData_(ds, 'a');
  assert.equal(a.scale, '');                                   // 변수정의에 척도가 비어 있으면 자동 판정
  assert.equal(a.autoScale, '연속');                            // 자동 판정 결과(모두 숫자)
  assert.equal(a.inferred, true);
  assert.equal(a.missingCodes, '99');
  assert.equal(a.missing, 2);                                   // 99(결측 코드) + 빈 칸
  assert.deepEqual(a.values.map(v => [v.key, v.count]), [['1', 3], ['2', 1]]);    // 결측 코드 99는 값으로 나오지 않는다
  assert.throws(() => V.buildEditorData_(ds, 'zz'), /찾을 수 없습니다/);
});

test('buildEditorData_: 값이 50가지를 넘으면 표를 만들지 않는다 (값 개수는 알려 준다)', () => {
  const rows = [['x']];
  for (let i = 0; i < 80; i++) rows.push([i]);
  const d = V.buildEditorData_(dsOf(rows), 'x');
  assert.equal(d.truncated, true);
  assert.equal(d.valueCount, 80);
  assert.deepEqual(d.values, []);
  assert.equal(V.buildEditorData_(dsOf([['x'], ...Array.from({ length: 50 }, (_, i) => [i])]), 'x').truncated, false);   // 정확히 50가지는 표를 만든다
});

test('findDefRow_: 있는 변수는 그 행, 없으면 마지막 정의 다음 행 (시트 행 번호)', () => {
  const header = ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'];
  const vals = [header, ['a', '', '', '', ''], ['b', '', '', '', ''], ['', '', '', '', ''], ['c', '', '', '', '']];
  assert.deepEqual(V.findDefRow_(vals, 'b'), { row: 3, isNew: false });
  assert.deepEqual(V.findDefRow_(vals, ' c '.trim()), { row: 5, isNew: false });
  assert.deepEqual(V.findDefRow_(vals, 'zz'), { row: 6, isNew: true });                // c(5행) 다음
  assert.deepEqual(V.findDefRow_([header], 'a'), { row: 2, isNew: true });              // 정의가 하나도 없으면 2행
  assert.deepEqual(V.findDefRow_([header, ['', '', '', '', '']], 'a'), { row: 2, isNew: true });
  assert.deepEqual(V.findDefRow_([header, [' a ', 1, 2, 3, 4]], 'a'), { row: 2, isNew: false });   // 앞뒤 공백은 무시
});

test('syncPlan_: 추가할 변수(중복 제거, 데이터 순서)와 데이터에 없는 정의', () => {
  assert.deepEqual(V.syncPlan_(['a', 'b', 'c', 'b'], ['b', 'x']), { toAdd: ['a', 'c'], orphans: ['x'] });
  assert.deepEqual(V.syncPlan_([], ['x']), { toAdd: [], orphans: ['x'] });
  assert.deepEqual(V.syncPlan_(['a'], ['a']), { toAdd: [], orphans: [] });
});

test('buildDefRow_: 저장할 다섯 칸과 입력 검증', () => {
  const headers = ['s', 'a'];
  const ok = V.buildDefRow_({ name: 's', label: ' 성별 ', scale: '명목', missing: '99 ;99', values: [{ key: '1', label: '남' }, { key: '2', label: '여' }] }, headers);
  assert.deepEqual(ok.cells, ['s', '성별', '명목', '1=남, 2=여', '99']);
  assert.deepEqual(V.buildDefRow_({ name: 'a', scale: '', values: [] }, headers).cells, ['a', '', '', '', '']);
  assert.throws(() => V.buildDefRow_({ name: 'q', values: [] }, headers), /데이터 탭 1행에 없는 변수/);
  assert.throws(() => V.buildDefRow_(null, headers), /데이터 탭 1행에 없는 변수/);
  assert.throws(() => V.buildDefRow_({ name: 's', scale: '연속형', values: [] }, headers), /척도는 명목, 서열, 연속/);
  assert.throws(() => V.buildDefRow_({ name: 's', label: 'x'.repeat(201), values: [] }, headers), /200자/);
  assert.throws(() => V.buildDefRow_({ name: 's', values: new Array(501).fill({ key: '1', label: 'a' }) }, headers), /너무 많습니다/);
  // keepValueLabels: 값이 너무 많아 표를 만들지 않은 변수는 기존 값 레이블을 그대로 유지한다
  const keep = V.buildDefRow_({ name: 's', keepValueLabels: true, values: [] }, headers, '1=남, 2=여');
  assert.equal(keep.cells[3], '1=남, 2=여');
  assert.equal(V.buildDefRow_({ name: 's', keepValueLabels: true, values: [] }, headers).cells[3], '');
});

/* ========== Apps Script 전역 범위 시뮬레이션: 시트 읽기·쓰기 ========== */

function loadGasWith(overrides) {
  const env = { console: { log() {}, error() {} }, Utilities: {}, Session: {}, PropertiesService: {} };
  Object.assign(env, overrides);
  const ctx = vm.createContext(env);
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  return ctx;
}

/** 값을 실제로 기억하는 변수정의 시트 mock (getDataRange 는 실제처럼 직사각형으로 돌려준다) */
function makeVarSheet(rows) {
  const log = { validations: [], formats: [] };
  const width = () => Math.max(1, ...rows.map(r => r.length));
  const sheet = {
    log, rows,
    getLastRow: () => rows.length,
    getDataRange: () => ({ getValues: () => rows.map(r => Array.from({ length: width() }, (_, j) => (j < r.length ? r[j] : ''))) }),
    getRange: (r, c, nr, nc) => {
      nr = nr === undefined ? 1 : nr; nc = nc === undefined ? 1 : nc;
      const rng = {
        setValues: (arr) => {
          if (arr.length !== nr || arr.some(a => a.length !== nc)) throw new Error('setValues 크기 불일치');
          arr.forEach((a, i) => {
            while (rows.length < r + i) rows.push([]);
            a.forEach((v, j) => { const row = rows[r - 1 + i]; while (row.length < c - 1 + j) row.push(''); row[c - 1 + j] = v; });
          });
          return rng;
        },
        setNumberFormat: (f) => { log.formats.push({ r, c, f }); return rng; },
        setDataValidation: (v) => { log.validations.push({ r, c, nr, nc, v }); return rng; }
      };
      return rng;
    }
  };
  return sheet;
}

function validationBuilder() {
  const b = { _: {}, requireValueInRange(range, show) { b._.range = range; return b; }, requireValueInList(list) { b._.list = list; return b; },
    setAllowInvalid(a) { b._.allowInvalid = a; return b; }, setHelpText(t) { b._.help = t; return b; }, build() { return b._; } };
  return b;
}

function editorApp(dataRows, varRows) {
  const varSheet = makeVarSheet(varRows);
  const dataSheet = {
    getLastRow: () => dataRows.length,
    getMaxColumns: () => 26,
    getDataRange: () => ({ getValues: () => JSON.parse(JSON.stringify(dataRows)) }),
    getRange: (r, c, nr, nc) => ({ tag: 'data-range', r, c, nr, nc })
  };
  const app = {
    newDataValidation: validationBuilder,
    getActiveSpreadsheet: () => ({ getSheetByName: (n) => (n === '데이터' ? dataSheet : n === '변수정의' ? varSheet : null) })
  };
  return { app, varSheet };
}

const HEADER = ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'];
const DATA = [['성별', '집단', '점수', '메모'], [1, 'A', 80, 'x'], [2, 'B', 75, ''], [1, 'A', 99, 'y'], [2, 'C', 90, ''], [1, 'B', 70, 'z']];

test('편집기 시뮬레이션: 변수 목록(정의 여부, 데이터에 없는 정의)과 한 변수의 값 목록', () => {
  const { app } = editorApp(DATA, [HEADER, ['성별', '성별', '명목', '1=남, 2=여', ''], ['예전변수', '', '', '', '']]);
  const ctx = loadGasWith({ SpreadsheetApp: app });
  const state = json(vm.runInContext('varEditorState()', ctx));
  assert.equal(state.n, 5);
  assert.deepEqual(state.vars.map(v => [v.name, v.defined]), [['성별', true], ['집단', false], ['점수', false], ['메모', false]]);
  assert.deepEqual(state.orphans, ['예전변수']);
  const d = json(vm.runInContext("varEditorLoad('성별')", ctx));
  assert.deepEqual(d.values.map(v => [v.key, v.count, v.label]), [['1', 3, '남'], ['2', 2, '여']]);
  assert.throws(() => vm.runInContext("varEditorLoad('없는변수')", ctx), /찾을 수 없습니다/);
});

test('편집기 시뮬레이션: 저장은 그 변수의 A~E열만 쓰고 다른 행·열은 건드리지 않는다 (새 변수는 마지막 정의 다음 행)', () => {
  const { app, varSheet } = editorApp(DATA, [HEADER.concat(['메모열']), ['성별', '', '', '', '', '내 메모'], ['집단', '', '', '', '', '집단 메모']]);
  const ctx = loadGasWith({ SpreadsheetApp: app });
  ctx.__p1 = { name: '성별', label: '성별', scale: '명목', missing: '', values: [{ key: '1', label: '남' }, { key: '2', label: '여' }] };
  const r1 = json(vm.runInContext('varEditorSave(__p1)', ctx));
  assert.deepEqual([r1.row, r1.isNew, r1.valueLabels, r1.warnings], [2, false, '1=남, 2=여', []]);
  assert.deepEqual(varSheet.rows[1].slice(0, 6), ['성별', '성별', '명목', '1=남, 2=여', '', '내 메모']);      // F열 메모 유지
  assert.deepEqual(varSheet.rows[2].slice(0, 6), ['집단', '', '', '', '', '집단 메모']);                      // 다른 행 그대로
  ctx.__p2 = { name: '점수', label: '시험 점수', scale: '연속', missing: '99', values: [] };
  const r2 = json(vm.runInContext('varEditorSave(__p2)', ctx));
  assert.deepEqual([r2.row, r2.isNew], [4, true]);
  assert.deepEqual(varSheet.rows[3].slice(0, 5), ['점수', '시험 점수', '연속', '', '99']);
  // 값 레이블 칸은 텍스트 서식으로 써서 "1=남" 같은 문자열이 수식·숫자로 해석되지 않게 한다
  assert.ok(varSheet.log.formats.some(f => f.c === 4 && f.f === '@'));
  // 입력 보조(변수명 드롭다운, 척도 드롭다운)가 걸린다
  const names = varSheet.log.validations.find(v => v.c === 1), scales = varSheet.log.validations.find(v => v.c === 3);
  assert.equal(names.v.range.tag, 'data-range');                     // 데이터 탭 1행을 목록으로
  assert.deepEqual([names.v.range.r, names.v.range.c, names.v.range.nr, names.v.range.nc], [1, 1, 1, 26]);
  assert.equal(names.v.allowInvalid, true);
  assert.deepEqual(json(scales.v.list), ['명목', '서열', '연속']);
});

test('편집기 시뮬레이션: 저장한 값 레이블과 결측값 코드가 분석에 그대로 반영된다 (끝에서 끝으로)', () => {
  const { app } = editorApp(DATA, [HEADER]);
  const ctx = loadGasWith({ SpreadsheetApp: app });
  ctx.__p = { name: '집단', label: '실험 집단', scale: '명목', missing: '', values: [{ key: 'A', label: '대조군' }, { key: 'B', label: '처치1' }, { key: 'C', label: '처치2' }] };
  vm.runInContext('varEditorSave(__p)', ctx);
  ctx.__q = { name: '점수', label: '', scale: '연속', missing: '99', values: [] };
  vm.runInContext('varEditorSave(__q)', ctx);
  ctx.__r = { analysis: 'frequencies', vars: ['집단'], options: {} };
  const res = json(vm.runInContext('runAnalysis(__r)', ctx).result);
  const table = res.blocks.find(b => b.caption === '실험 집단');
  assert.deepEqual(table.rows.slice(0, 3).map(r => [r.value, r.freq]), [['대조군', 2], ['처치1', 2], ['처치2', 1]]);
  ctx.__d = { analysis: 'descriptives', vars: ['점수'], options: {} };
  const desc = json(vm.runInContext('runAnalysis(__d)', ctx).result);
  assert.equal(desc.blocks[0].rows[0].n, 4);                          // 99 는 결측: 5건 중 4건
});

test('편집기 시뮬레이션: 값이 너무 많은 변수는 저장해도 기존 값 레이블을 유지한다', () => {
  const rows = [['x']];
  for (let i = 0; i < 80; i++) rows.push([i]);
  const { app, varSheet } = editorApp(rows, [HEADER, ['x', '', '', '0=영, 1=일', '']]);
  const ctx = loadGasWith({ SpreadsheetApp: app });
  ctx.__p = { name: 'x', label: '값 많음', scale: '연속', missing: '', values: [], keepValueLabels: true };
  vm.runInContext('varEditorSave(__p)', ctx);
  assert.deepEqual(varSheet.rows[1].slice(0, 5), ['x', '값 많음', '연속', '0=영, 1=일', '']);
});

test('편집기 시뮬레이션: 정의 없는 변수 추가는 이름만 채우고, 두 번째 실행은 추가할 것이 없다', () => {
  const { app, varSheet } = editorApp(DATA, [HEADER, ['성별', '성별', '명목', '1=남, 2=여', ''], ['예전변수', '', '', '', '']]);
  const ctx = loadGasWith({ SpreadsheetApp: app });
  const r = json(vm.runInContext('varEditorSync()', ctx));
  assert.deepEqual(r.added, ['집단', '점수', '메모']);
  assert.deepEqual(r.orphans, ['예전변수']);
  assert.deepEqual(varSheet.rows.slice(1).map(x => x[0]), ['성별', '예전변수', '집단', '점수', '메모']);
  assert.deepEqual(varSheet.rows[1].slice(0, 5), ['성별', '성별', '명목', '1=남, 2=여', '']);          // 기존 정의는 그대로
  assert.equal(varSheet.rows[3][0], '집단');                                                           // 추가한 행은 이름만 (나머지 칸은 비어 있다)
  assert.ok(varSheet.rows[3].slice(1, 5).every(x => x === '' || x === undefined));
  const again = json(vm.runInContext('varEditorSync()', ctx));
  assert.deepEqual(again.added, []);
  assert.ok(varSheet.log.validations.length >= 2);
});

test('편집기 시뮬레이션: 데이터에 없는 변수 저장은 거부, 데이터가 비어 있으면 안내', () => {
  const { app } = editorApp(DATA, [HEADER]);
  const ctx = loadGasWith({ SpreadsheetApp: app });
  ctx.__bad = { name: '없는변수', label: '', scale: '', missing: '', values: [] };
  assert.throws(() => vm.runInContext('varEditorSave(__bad)', ctx), /데이터 탭 1행에 없는 변수/);
  const empty = editorApp([['']], [HEADER]);
  const ctx2 = loadGasWith({ SpreadsheetApp: empty.app });
  assert.throws(() => vm.runInContext('varEditorState()', ctx2), /1행에 변수명을 입력하고/);
});

test('템플릿 만들기: 변수정의 탭은 머리글과 변수명·척도 드롭다운을 갖고, 데이터 탭 A1은 값 없이 메모만 (예제 채우기가 만든 탭도 같다)', () => {
  const made = [];
  function newSheet(name) {
    const log = { notes: [], values: [], validations: [], bold: [] };
    const sh = {
      name, log, getMaxColumns: () => 26, setFrozenRows() {}, setColumnWidth() {},
      getRange: (r, c, nr, nc) => {
        const rng = {
          setNote: (n) => { log.notes.push({ r, c, n }); return rng; },
          setValues: (v) => { log.values.push({ r, c, v }); return rng; },
          setFontWeight: () => rng, setBackground: () => rng,
          setDataValidation: (v) => { log.validations.push({ c, v }); return rng; },
          setValue: (v) => { log.values.push({ r, c, v }); return rng; }
        };
        return rng;
      }
    };
    made.push(sh);
    return sh;
  }
  const existing = {};
  const app = {
    newDataValidation: validationBuilder,
    getActiveSpreadsheet: () => ({
      getSheetByName: (n) => existing[n] || null,
      insertSheet: (n) => (existing[n] = newSheet(n))
    }),
    getUi: () => ({ alert() {}, ButtonSet: { OK: 'OK', YES_NO: 'YES_NO' }, Button: { YES: 'YES' } })
  };
  const ctx = loadGasWith({ SpreadsheetApp: app });
  vm.runInContext('createTemplate()', ctx);
  const data = existing['데이터'], vars = existing['변수정의'];
  assert.equal(data.log.values.length, 0);                                    // A1에 값을 쓰지 않는다 (변수 이름으로 읽히므로)
  assert.equal(data.log.notes.length, 1);
  assert.ok(data.log.notes[0].n.includes('1행에 변수명'));
  assert.deepEqual(json(vars.log.values[0].v), [HEADER]);                     // 표준 머리글
  const cols = vars.log.validations.map(v => v.c).sort();
  assert.deepEqual(cols, [1, 3]);                                             // 변수명(A)과 척도(C) 드롭다운
  assert.equal(vars.log.validations.find(v => v.c === 1).v.range.tag, undefined);   // 이 mock 의 getRange 는 tag 가 없다: 호출만 확인
});
