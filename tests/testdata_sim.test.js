const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const TD = require('../src/TestData.js');

const SRC = path.join(__dirname, '..', 'src');
const HEADER = ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'];
const plain = x => JSON.parse(JSON.stringify(x));

/* ========== 시트 시뮬레이션: 시험 데이터 넣기, 백업, 되돌리기, 그 데이터로 분석 ========== */

function loadGasWith(overrides) {
  const env = { console: { log() {}, error() {} } };
  Object.assign(env, overrides);
  const ctx = vm.createContext(env);
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  return ctx;
}

/** 시트(탭)와 문서 속성을 기억하는 작은 통합 문서 mock */
function makeWorkbook(initial) {
  const sheets = [];
  const props = {};
  function mkSheet(name, rows) {
    const sh = {
      rows: rows.map(r => r.slice()), name: name, maxCols: 26, calls: [],
      getName: () => sh.name,
      setName(n) { sh.name = n; return sh; },
      setTabColor() { return sh; },
      copyTo() { return mkSheet(sh.name + '의 사본', sh.rows); },
      activate() { sh.calls.push('activate'); return sh; },
      clearContents() { sh.rows = []; return sh; },
      setFrozenRows() {}, setColumnWidth() {},
      getMaxColumns: () => sh.maxCols, getMaxRows: () => 1000,
      insertColumnsAfter(pos, n) { sh.maxCols += n; sh.calls.push('insertColumns'); },
      getLastRow: () => sh.rows.length,
      getDataRange: () => ({
        getValues: () => {
          const w = Math.max(1, ...sh.rows.map(r => r.length));
          return sh.rows.map(r => Array.from({ length: w }, (_, j) => (j < r.length ? r[j] : '')));
        }
      }),
      getRange(r, c, nr, nc) {
        nr = nr === undefined ? 1 : nr; nc = nc === undefined ? 1 : nc;
        const rng = {
          setValues(arr) {
            arr.forEach((a, i) => {
              while (sh.rows.length < r + i) sh.rows.push([]);
              a.forEach((v, j) => { const row = sh.rows[r - 1 + i]; while (row.length < c - 1 + j) row.push(''); row[c - 1 + j] = v; });
            });
            return rng;
          },
          clearContent() {
            for (let i = r - 1; i < Math.min(sh.rows.length, r - 1 + nr); i++) {
              for (let j = c - 1; j < c - 1 + nc; j++) if (j < sh.rows[i].length) sh.rows[i][j] = '';
            }
            return rng;
          }
        };
        ['setNumberFormat', 'setFontWeight', 'setBackground', 'setFontSize', 'setWrap', 'setVerticalAlignment', 'setNote', 'setDataValidation'].forEach(m => { rng[m] = () => rng; });
        return rng;
      }
    };
    sheets.push(sh);
    return sh;
  }
  Object.entries(initial).forEach(([n, rows]) => mkSheet(n, rows));
  const ss = {
    getSheetByName: n => sheets.find(s => s.name === n) || null,
    insertSheet: n => mkSheet(n, []),
    deleteSheet: s => { sheets.splice(sheets.indexOf(s), 1); },
    getSheets: () => sheets
  };
  const v = () => { const b = { requireValueInRange: () => b, requireValueInList: () => b, setAllowInvalid: () => b, setHelpText: () => b, build: () => ({}) }; return b; };
  return {
    sheets, ss, props,
    env: {
      SpreadsheetApp: { newDataValidation: v, getActiveSpreadsheet: () => ss },
      Utilities: { formatDate: () => '20261023_1530' },
      Session: { getScriptTimeZone: () => 'Asia/Seoul' },
      PropertiesService: {
        getDocumentProperties: () => ({ getProperty: k => (k in props ? props[k] : null), setProperty: (k, val) => { props[k] = val; }, deleteProperty: k => { delete props[k]; } }),
        getUserProperties: () => ({ getProperty: () => null, setProperty() {} })
      }
    }
  };
}
const ORIGINAL_DATA = [['a', 'b'], [1, 2], [3, 4]];
const ORIGINAL_VARS = [HEADER, ['a', '', '연속', '', ''], ['b', '', '연속', '', '']];

test('시험 데이터 도구: 목록에는 기본 예제가 앞에 오고 기대값이 문자열로 정리된다', () => {
  const wb = makeWorkbook({ '데이터': ORIGINAL_DATA, '변수정의': ORIGINAL_VARS });
  const ctx = loadGasWith(wb.env);
  const r = plain(vm.runInContext('listTestDatasets()', ctx));
  assert.equal(r.items[0].id, 'sample');
  assert.equal(r.items.length, TD.TEST_DATASETS_.length + 1);
  assert.equal(r.backup, null);
  const lg = r.items.find(i => i.id === 'logistic-1');
  assert.equal(lg.analysisLabel, '로지스틱 회귀');
  assert.ok(lg.expected.length >= 5 && lg.expected.every(e => typeof e.value === 'string'));
  assert.ok(lg.n === 100 && lg.vars.join() === 'y,x0,x1');
});

test('시험 데이터 도구: 넣으면 백업 탭을 만들고 데이터·변수정의·시험안내가 바뀌며, 다시 넣어도 백업은 처음 것이 유지된다', () => {
  const wb = makeWorkbook({ '데이터': ORIGINAL_DATA, '변수정의': ORIGINAL_VARS, '내 메모': [['중요']] });
  const ctx = loadGasWith(wb.env);
  const entry = TD.TEST_DATASETS_.find(d => d.id === 'mw-1');
  const r1 = plain(vm.runInContext("insertTestDataset('mw-1')", ctx));
  assert.equal(r1.rows, 12);
  const names = () => wb.sheets.map(s => s.name);
  assert.deepEqual(names().sort(), ['내 메모', '데이터', '데이터_백업_20261023_1530', '변수정의', '변수정의_백업_20261023_1530', '시험안내'].sort());
  const data = wb.ss.getSheetByName('데이터'), vars = wb.ss.getSheetByName('변수정의');
  assert.deepEqual(plain(data.rows[0]), entry.header);
  assert.equal(data.rows.length, 13);
  assert.deepEqual(plain(vars.rows.slice(1)), entry.varRows);
  assert.deepEqual(plain(wb.ss.getSheetByName('데이터_백업_20261023_1530').rows), ORIGINAL_DATA);
  const guide = wb.ss.getSheetByName('시험안내');
  assert.ok(guide.rows[0][0].includes('Mann-Whitney U'));
  assert.ok(guide.rows.some(r => String(r[0]).includes('Mann-Whitney U') && r[1] !== ''));          // 기대값 행
  assert.ok(guide.calls.includes('activate'));
  assert.ok(wb.props.STAT4SW_TEST_BACKUP);
  // 다른 시험 데이터로 바꿔 넣어도 백업은 늘지 않고 원래 내용 그대로
  vm.runInContext("insertTestDataset('logistic-1')", ctx);
  assert.equal(wb.sheets.length, 6);
  assert.deepEqual(plain(wb.ss.getSheetByName('데이터_백업_20261023_1530').rows), ORIGINAL_DATA);
  assert.equal(wb.ss.getSheetByName('데이터').rows[0].join(), 'y,x0,x1');
  assert.deepEqual(plain(wb.ss.getSheetByName('내 메모').rows), [['중요']]);                       // 다른 탭은 그대로
});

test('시험 데이터 도구: 되돌리기는 백업을 원래 이름으로 복원하고 시험 탭을 지우며, 두 번 하면 안내 오류', () => {
  const wb = makeWorkbook({ '데이터': ORIGINAL_DATA, '변수정의': ORIGINAL_VARS });
  const ctx = loadGasWith(wb.env);
  vm.runInContext("insertTestDataset('kw-1')", ctx);
  const r = plain(vm.runInContext('restoreTestBackup()', ctx));
  assert.ok(r.message.includes('되돌렸습니다'));
  assert.deepEqual(wb.sheets.map(s => s.name).sort(), ['데이터', '변수정의']);
  assert.deepEqual(plain(wb.ss.getSheetByName('데이터').rows), ORIGINAL_DATA);
  assert.deepEqual(plain(wb.ss.getSheetByName('변수정의').rows), ORIGINAL_VARS);
  assert.equal(wb.props.STAT4SW_TEST_BACKUP, undefined);
  assert.throws(() => vm.runInContext('restoreTestBackup()', ctx), /되돌릴 백업이 없습니다/);
});

test('시험 데이터 도구: 원래 데이터·변수정의 탭이 없던 시트에서는 새로 만들었다가 되돌릴 때 지운다', () => {
  const wb = makeWorkbook({ 'Sheet1': [['x']] });
  const ctx = loadGasWith(wb.env);
  vm.runInContext("insertTestDataset('reliability-1')", ctx);
  assert.deepEqual(wb.sheets.map(s => s.name).sort(), ['Sheet1', '데이터', '변수정의', '시험안내'].sort());
  assert.equal(wb.ss.getSheetByName('데이터').rows.length, 61);
  vm.runInContext('restoreTestBackup()', ctx);
  assert.deepEqual(wb.sheets.map(s => s.name), ['Sheet1']);
});

test('시험 데이터 도구: 기본 예제는 변수가 기본 26열을 넘으면 열을 늘리고, 모르는 id 는 거부한다', () => {
  const wb = makeWorkbook({ '데이터': ORIGINAL_DATA, '변수정의': ORIGINAL_VARS });
  const ctx = loadGasWith(wb.env);
  const r = plain(vm.runInContext("insertTestDataset('sample')", ctx));
  assert.equal(r.columns, 27);
  assert.ok(wb.ss.getSheetByName('데이터').calls.includes('insertColumns'));
  assert.equal(wb.ss.getSheetByName('데이터').rows[0].length, 27);
  assert.throws(() => vm.runInContext("insertTestDataset('없는id')", ctx), /알 수 없는 시험 데이터/);
});

/** 시험 데이터를 시트에 넣은 뒤 같은 요청으로 분석을 돌려 기대값과 비교한다 (시트 읽기와 변수정의 해석까지 포함한 끝에서 끝) */
for (const id of ['logistic-2', 'logistic-4', 'logistic-5', 'logistic-6', 'reliability-2', 'moderation-2', 'mw-2', 'distance-2', 'psm-1', 'kmeans-1', 'cart-1', 'cart-3', 'mcnemar-1', 'effect-1', 'friedman-1', 'rmanova-1', 'hierarchical-2', 'efa-1', 'efa-4', 'anova2-2', 'mixed-1', 'kappa-3', 'icc-2', 'ologit-1', 'partcorr-1', 'km-1', 'km-2', 'anova-1', 'anova-2', 'binom-1', 'binom-2', 'cochran-1', 'cochran-2', 'kendall-1', 'kendall-2', 'ipa-1', 'ipa-2']) {
  test(`시험 데이터 ${id}: 시트에 넣은 뒤 runAnalysis 로 돌려도 기준값과 일치한다 (끝에서 끝)`, () => {
    const wb = makeWorkbook({ '데이터': ORIGINAL_DATA, '변수정의': ORIGINAL_VARS });
    const ctx = loadGasWith(wb.env);
    const entry = TD.TEST_DATASETS_.find(d => d.id === id);
    vm.runInContext("insertTestDataset('" + id + "')", ctx);
    ctx.__req = plain(entry.request);
    const out = plain(vm.runInContext('runAnalysis(__req)', ctx));
    let checked = 0;
    entry.expected.forEach(spec => {
      if (spec.info) return;
      const t = out.result.blocks.find(b => b.type === 'table' && b.caption.startsWith(spec.table));
      let got = t.rows[spec.row][spec.key], want = spec.value;
      if (spec.abs) { got = Math.abs(got); want = Math.abs(want); }
      assert.ok(Math.abs(got - want) <= spec.tol * Math.max(1e-9, Math.abs(want)), id + ' ' + spec.label + ': ' + got + ' vs ' + want);
      checked++;
    });
    assert.ok(checked >= 2);
  });
}
