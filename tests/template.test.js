const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const SRC = path.join(__dirname, '..', 'src');
function loadGasWith(overrides) {
  const ctx = vm.createContext(Object.assign({ console: { log() {}, error() {} }, Utilities: {}, Session: {}, PropertiesService: {} }, overrides));
  for (const f of fs.readdirSync(SRC).filter(x => x.endsWith('.js')).sort()) vm.runInContext(fs.readFileSync(path.join(SRC, f), 'utf8'), ctx, { filename: f });
  return ctx;
}

/** 내용이 있는 데이터·변수정의 탭을 흉내 낸 mock. 지우기(clear)와 다시 만들기 호출을 기록한다 */
function makeApp(state) {
  const mkSheet = (name, rows, cols) => {
    const sh = {
      name, rows, cols, log: { clears: 0, notes: 0, headers: [], validations: 0, frozen: 0 },
      getName: () => name, getLastRow: () => sh.rows, getLastColumn: () => sh.cols, getMaxColumns: () => 26,
      clear() { sh.log.clears++; sh.rows = 0; sh.cols = 0; return sh; },
      setFrozenRows() { sh.log.frozen++; }, setColumnWidth() {},
      getRange: (r, c, nr, nc) => {
        const rng = {
          setNote: () => { sh.log.notes++; return rng; }, setValues: (v) => { sh.log.headers.push(v[0]); return rng; },
          setFontWeight: () => rng, setBackground: () => rng, setDataValidation: () => { sh.log.validations++; return rng; }
        };
        return rng;
      }
    };
    return sh;
  };
  const sheets = {};
  if (state.data) sheets['데이터'] = mkSheet('데이터', state.data[0], state.data[1]);
  if (state.vars) sheets['변수정의'] = mkSheet('변수정의', state.vars, 5);
  if (state.result) sheets['결과'] = mkSheet('결과', state.result, 3);
  (state.matrix || []).forEach(n => { sheets[n] = mkSheet(n, 5, 5); });
  state.deleted = [];
  state.sheets = sheets;
  const v = () => { const b = { requireValueInRange: () => b, requireValueInList: () => b, setAllowInvalid: () => b, setHelpText: () => b, build: () => ({}) }; return b; };
  return {
    newDataValidation: v,
    getUi: () => ({
      ButtonSet: { OK: 'OK', YES_NO: 'YES_NO', OK_CANCEL: 'OK_CANCEL' }, Button: { YES: 'YES', NO: 'NO', OK: 'OK', CANCEL: 'CANCEL' },
      alert: (title, msg, buttons) => { state.alerts.push({ title, msg, buttons }); return buttons === 'YES_NO' ? state.answer : 'OK'; },
      prompt: (title, msg, buttons) => {
        (state.prompts = state.prompts || []).push({ title, msg, buttons });
        return { getSelectedButton: () => state.promptButton || 'OK', getResponseText: () => (state.promptText === undefined ? '지우기' : state.promptText) };
      }
    }),
    getActiveSpreadsheet: () => ({
      getSheets: () => Object.keys(sheets).map(k => sheets[k]),
      deleteSheet: (sh) => { state.deleted.push(sh.getName()); delete sheets[sh.getName()]; },
      getSheetByName: (n) => sheets[n] || null,
      insertSheet: (n) => { state.inserted.push(n); return (sheets[n] = mkSheet(n, 0, 0)); }
    })
  };
}
const base = (extra) => Object.assign({ data: [61, 27], vars: [28], alerts: [], inserted: [], answer: 'YES' }, extra);

test('템플릿 시트 만들기: 내용이 있으면 지워질 것을 정확히 안내하고, 예를 누르면 두 탭을 비우고 템플릿으로 다시 만든다', () => {
  const s = base();
  const ctx = loadGasWith({ SpreadsheetApp: makeApp(s) });
  vm.runInContext('createTemplate()', ctx);
  assert.equal(s.alerts[0].title, '템플릿 시트 만들기');
  assert.equal(s.alerts[0].buttons, 'YES_NO');
  assert.equal(s.alerts[0].msg, [
    '아래 탭의 내용(값, 서식, 메모)을 모두 지우고 빈 템플릿으로 다시 만듭니다.',
    '· "데이터" 탭: 변수 27개, 사례 60행',
    '· "변수정의" 탭: 정의 27개',
    '',
    '· 이 도구에는 되돌리기가 없습니다. 필요하면 먼저 시트를 복사해 두세요 (탭 이름을 우클릭 → 복사).',
    '· 실수로 지웠다면 구글 시트의 "파일 > 버전 기록 > 버전 기록 보기"에서 이전 시점으로 되돌릴 수 있습니다 (시트 전체가 그 시점으로 돌아갑니다).',
    '· 분석 결과로 데이터 탭에 저장한 열(예측값, 군집 번호 등)도 함께 지워집니다.',
    '· "결과" 탭과 다른 탭은 그대로입니다.',
    '',
    '계속할까요?'
  ].join('\n'));
  const { 데이터: data, 변수정의: vars } = s.sheets;
  assert.equal(data.log.clears, 1);
  assert.equal(vars.log.clears, 1);
  assert.equal(data.log.notes, 1);                                                   // 템플릿 모양으로 다시 만든다 (A1 메모)
  assert.deepEqual(JSON.parse(JSON.stringify(vars.log.headers[0])), ['변수명', '라벨', '척도', '값 레이블', '결측값 코드']);
  assert.ok(vars.log.validations >= 2);
  assert.deepEqual(s.inserted, []);                                                  // 탭을 새로 만들지는 않는다
  assert.ok(s.alerts[1].msg.includes('"데이터" 탭을 비우고 템플릿으로 다시 만듦') && s.alerts[1].msg.includes('"변수정의" 탭을 비우고'));
});

test('템플릿 시트 만들기: 아니오를 누르면 아무것도 바꾸지 않고, 비어 있거나 없는 탭은 확인 없이 만든다', () => {
  let s = base({ answer: 'NO' });
  vm_run(s);
  assert.equal(s.alerts.length, 1);
  assert.equal(s.sheets['데이터'].log.clears + s.sheets['변수정의'].log.clears, 0);

  s = base({ data: [0, 0], vars: [1] });                                              // 머리글만 있는 빈 템플릿: 확인 없이 다시 만든다
  vm_run(s);
  assert.deepEqual(s.alerts.map(a => a.buttons), ['OK']);
  assert.equal(s.sheets['데이터'].log.clears, 1);

  s = base({ data: null, vars: null });                                               // 탭이 없으면 만든다
  vm_run(s);
  assert.deepEqual(s.inserted.slice().sort(), ['데이터', '변수정의']);
  assert.deepEqual(s.alerts.map(a => a.buttons), ['OK']);

  s = base({ data: [5, 3], vars: null });                                             // 데이터만 있으면 그것만 안내하고 변수정의는 새로 만든다
  vm_run(s);
  assert.ok(s.alerts[0].msg.includes('"데이터" 탭: 변수 3개, 사례 4행') && !s.alerts[0].msg.includes('"변수정의" 탭: 정의'));
  assert.deepEqual(s.inserted, ['변수정의']);
  function vm_run(st) { vm.runInContext('createTemplate()', loadGasWith({ SpreadsheetApp: makeApp(st) })); }
});

/* ---------- 새 작업 시작 (데이터 모두 지우기) ---------- */

const runNew = (s) => vm.runInContext('startNewWork()', loadGasWith({ SpreadsheetApp: makeApp(s) }));
const cleared = (s) => ['데이터', '변수정의'].reduce((n, k) => n + (s.sheets[k] ? s.sheets[k].log.clears : 0), 0);

test('새 작업 시작: 내용이 있으면 지워질 것을 안내하고 "지우기"를 입력해야 두 탭을 비운다 (예/아니오 확인창이 아니라 입력 확인)', () => {
  const s = base();
  runNew(s);
  assert.equal(s.prompts.length, 1);
  assert.equal(s.prompts[0].title, '새 작업 시작 (데이터 모두 지우기)');
  assert.equal(s.prompts[0].buttons, 'OK_CANCEL');
  assert.ok(s.prompts[0].msg.includes('· "데이터" 탭: 변수 27개, 사례 60행') && s.prompts[0].msg.includes('· "변수정의" 탭: 정의 27개'));
  assert.ok(s.prompts[0].msg.includes('아래 칸에 지우기 라고 입력하고') && s.prompts[0].msg.includes('버전 기록') && s.prompts[0].msg.includes('먼저 시트를 복사'));
  assert.equal(cleared(s), 2);                                                         // 두 탭이 비워지고 템플릿으로 다시 만들어진다
  assert.deepEqual(s.alerts.map(a => [a.title, a.buttons]), [['새 작업 시작', 'OK']]);  // 결과 탭이 없으면 결과를 묻지 않는다
});

test('새 작업 시작: 취소하거나 문구를 다르게 입력하면 아무것도 바꾸지 않는다', () => {
  let s = base({ promptButton: 'CANCEL' });
  runNew(s);
  assert.equal(cleared(s), 0);
  assert.equal(s.alerts.length, 0);
  s = base({ promptText: '지우' });
  runNew(s);
  assert.equal(cleared(s), 0);
  assert.equal(s.alerts.length, 1);
  assert.ok(s.alerts[0].msg.includes('정확히 입력하지 않아 취소했습니다') && s.alerts[0].msg.includes('아무것도 바뀌지 않았습니다'));
  s = base({ promptText: ' 지우기 ' });                                                // 앞뒤 공백은 허용
  runNew(s);
  assert.equal(cleared(s), 2);
});

test('새 작업 시작: 탭이 없거나(새 시트, 실수로 삭제) 비어 있으면 묻지 않고 만들고, 한 탭만 지워졌으면 그것만 새로 만든다', () => {
  let s = base({ data: null, vars: null });
  runNew(s);
  assert.deepEqual(s.inserted.slice().sort(), ['데이터', '변수정의']);
  assert.ok(!s.prompts || s.prompts.length === 0);
  assert.deepEqual(s.alerts.map(a => a.buttons), ['OK']);
  s = base({ data: [0, 0], vars: [1] });                                               // 머리글만 있는 빈 템플릿
  runNew(s);
  assert.ok(!s.prompts || s.prompts.length === 0);
  assert.equal(cleared(s), 2);
  s = base({ data: null, vars: [28] });                                                // 데이터 탭을 실수로 지움, 변수정의에는 내용이 있음
  runNew(s);
  assert.equal(s.prompts.length, 1);
  assert.ok(s.prompts[0].msg.includes('"변수정의" 탭: 정의 27개') && !s.prompts[0].msg.includes('"데이터" 탭:'));
  assert.deepEqual(s.inserted, ['데이터']);
  assert.equal(s.sheets['변수정의'].log.clears, 1);
});

test('새 작업 시작: 끝나고 이전 분석 결과(결과 탭, 행렬 탭)도 지울지 묻고, 아니오면 그대로 둔다', () => {
  let s = base({ result: 120, matrix: ['거리행렬', '유사도행렬_2'], answer: 'NO' });
  runNew(s);
  const q = s.alerts.find(a => a.buttons === 'YES_NO');
  assert.ok(q && q.msg.includes('이전 분석 결과도 함께 지울까요?') && q.msg.includes('"결과" 탭의 내용 전체 (120행') && q.msg.includes('"거리행렬" 탭') && q.msg.includes('"유사도행렬_2" 탭'));
  assert.equal(s.sheets['결과'].log.clears, 0);
  assert.deepEqual(s.deleted, []);
  s = base({ result: 120, matrix: ['거리행렬'], answer: 'YES' });
  runNew(s);
  assert.equal(s.sheets['결과'].log.clears, 1);
  assert.deepEqual(s.deleted, ['거리행렬']);
  assert.ok(s.alerts[s.alerts.length - 1].msg.includes('이전 분석 결과를 지웠습니다'));
  s = base({ result: 0 });                                                             // 결과 탭이 비어 있고 행렬 탭도 없으면 묻지 않는다
  runNew(s);
  assert.ok(!s.alerts.some(a => a.buttons === 'YES_NO'));
});
