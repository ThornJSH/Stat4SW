/**
 * 분석 창 UI를 브라우저에서 로컬로 확인하기 위한 미리보기 페이지를 만든다.
 * 실행: node tests/ui/build_preview.js  →  tests/ui/preview.html (git 제외)
 *
 * 실제 서버 코드(Data, Stats_*, Render, Api)를 <script>로 그대로 불러오고,
 * SpreadsheetApp 과 google.script.run 만 흉내 낸다. GAS처럼 전역 범위에서 돌아가므로
 * 로드 순서나 전역 함수 참조 문제도 함께 걸러진다. (Apps Script 환경 자체의 검증은 아니다)
 */
const fs = require('node:fs');
const path = require('node:path');

const root = path.join(__dirname, '..', '..', 'src');
const read = f => fs.readFileSync(path.join(root, f), 'utf8');

// GAS처럼 src/*.js 전체를 한 전역 범위로 불러온다. --reverse 는 역순으로 불러 로드 순서 의존을 점검한다.
// (Apps Script가 실제로 어떤 순서로 합치는지는 보장되지 않으므로 어느 순서에서도 동작해야 한다)
const serverFiles = fs.readdirSync(root).filter(f => f.endsWith('.js')).sort();
if (process.argv.includes('--reverse')) serverFiles.reverse();

const sample = require('../../src/Sample.js').sampleDataset_();
const mock = `
// ---- 시트 흉내: Sample.js 의 예제 데이터 (분석 창과 같은 값) ----
var __sample = ${JSON.stringify(sample)};
var __data = [__sample.header].concat(__sample.rows);
var __vars = [['변수명', '라벨', '척도', '값 레이블', '결측값 코드']].concat(__sample.varRows);
// 결과 탭: 쓴 내용을 window.__exported 에 기록하는 mock (실제 서식은 시트에서만 확인할 수 있다)
window.__exported = [];
window.__images = [];   // insertImage 로 들어온 그림 (PNG 바이트 길이, 시그니처, 위치)
var __resultSheet = {
  getLastRow: function () { return 0; }, getMaxRows: function () { return 5000; }, getMaxColumns: function () { return 30; },
  insertRowsAfter: function () {}, insertColumnsAfter: function () {}, setColumnWidths: function () {}, activate: function () {},
  insertImage: function (blob, col, row) {
    var rec = { col: col, row: row, type: blob.contentType, name: blob.name, bytes: blob.data.length,
      sig: Array.prototype.slice.call(blob.data, 0, 4).join(','),
      b64: btoa(String.fromCharCode.apply(null, Array.prototype.slice.call(blob.data))) };
    window.__images.push(rec);
    var img = { setWidth: function (w) { rec.width = w; return img; }, setHeight: function (h) { rec.height = h; return img; } };
    return img;
  },
  getRange: function (r, c, nr, nc) {
    var rng = {};
    ['setNumberFormats', 'setFontWeights', 'setFontSizes', 'setFontColors', 'setFontStyles', 'setBackgrounds', 'setHorizontalAlignments', 'setBorder'].forEach(function (m) { rng[m] = function () { return rng; }; });
    rng.setValues = function (arr) { window.__exported.push({ row: r, values: arr }); return rng; };
    rng.activate = function () { return rng; };
    return rng;
  }
};
// 변수정의 탭: 편집기가 쓴 값을 기억한다 (window.__varRows 로 확인)
window.__varRows = __vars;
// 데이터 탭: 변수 계산이 쓴 새 열도 기억한다 (window.__dataRows 로 확인)
window.__dataRows = __data;
var __dataSheet = {
  getLastRow: function () { return __data.length; },
  getMaxColumns: function () { return 40; },
  insertColumnsAfter: function () {},
  getDataRange: function () {
    var w = Math.max.apply(null, __data.map(function (r) { return r.length; }));
    return { getValues: function () { return __data.map(function (r) { var o = r.slice(); while (o.length < w) o.push(''); return o; }); } };
  },
  getRange: function (r, c, nr, nc) {
    nr = nr === undefined ? 1 : nr; nc = nc === undefined ? 1 : nc;
    var rng = {
      setValues: function (arr) {
        arr.forEach(function (a, i) {
          while (__data.length < r + i) __data.push([]);
          a.forEach(function (v, j) { var row = __data[r - 1 + i]; while (row.length < c - 1 + j) row.push(''); row[c - 1 + j] = v; });
        });
        return rng;
      },
      setFontWeight: function () { return rng; }, setBackground: function () { return rng; }, setNumberFormat: function () { return rng; }
    };
    return rng;
  }
};
var __varSheet = {
  getLastRow: function () { return __vars.length; },
  getDataRange: function () {
    var w = Math.max.apply(null, __vars.map(function (r) { return r.length; }));
    return { getValues: function () { return __vars.map(function (r) { var o = r.slice(); while (o.length < w) o.push(''); return o; }); } };
  },
  getRange: function (r, c, nr, nc) {
    nr = nr === undefined ? 1 : nr; nc = nc === undefined ? 1 : nc;
    var rng = {
      setValues: function (arr) {
        arr.forEach(function (a, i) {
          while (__vars.length < r + i) __vars.push([]);
          a.forEach(function (v, j) { var row = __vars[r - 1 + i]; while (row.length < c - 1 + j) row.push(''); row[c - 1 + j] = v; });
        });
        return rng;
      },
      setNumberFormat: function () { return rng; },
      setDataValidation: function () { return rng; }
    };
    return rng;
  }
};
function __validation() {
  var b = {};
  ['requireValueInRange', 'requireValueInList', 'setAllowInvalid', 'setHelpText'].forEach(function (m) { b[m] = function () { return b; }; });
  b.build = function () { return {}; };
  return b;
}
var SpreadsheetApp = {
  newDataValidation: __validation,
  BorderStyle: { SOLID: 'SOLID', SOLID_MEDIUM: 'SOLID_MEDIUM' },
  getActiveSpreadsheet: function () { return {
    insertSheet: function () { return __resultSheet; },
    getSheetByName: function (n) {
      if (n === '변수정의') return __varSheet;
      return n === '데이터' ? __dataSheet : null;
    }
  }; }
};
var Utilities = {
  formatDate: function () { return '2026-10-21 12:00:00'; },
  base64Decode: function (s) { var bin = atob(s), a = new Uint8Array(bin.length); for (var i = 0; i < bin.length; i++) a[i] = bin.charCodeAt(i); return a; },
  newBlob: function (data, contentType, name) { return { data: data, contentType: contentType, name: name }; }
};
var Session = { getScriptTimeZone: function () { return 'Asia/Seoul'; } };
var __props = {};
var PropertiesService = { getUserProperties: function () { return {
  getProperty: function (k) { return k in __props ? __props[k] : null; },
  setProperty: function (k, v) { __props[k] = v; window.__savedProps = __props; }
}; } };
// ---- google.script.run 흉내: 서버 전역 함수를 비동기로 호출 ----
function __makeRunner(ok, fail) {
  return new Proxy({}, { get: function (_, name) {
    if (name === 'withSuccessHandler') return function (f) { return __makeRunner(f, fail); };
    if (name === 'withFailureHandler') return function (f) { return __makeRunner(ok, f); };
    return function () {
      var args = arguments;
      setTimeout(function () {
        try { var r = window[name].apply(null, args); if (ok) ok(r); }
        catch (e) { if (fail) fail(e); else throw e; }
      }, 30);
    };
  } });
}
var google = { script: { run: __makeRunner(null, null) } };
`;

const include = name => read(name + '.html');
let page = read('index.html')
  .replace("<?!= include('style') ?>", include('style'))
  .replace("<?!= include('client') ?>",
    '<script>' + mock + '</script>\n' +
    serverFiles.map(f => '<script>\n' + read(f) + '\n</script>').join('\n') + '\n' +
    include('client'));

const out = path.join(__dirname, 'preview.html');
fs.writeFileSync(out, page);
console.log('wrote', out);

// 변수 설정 편집기(varedit.html) 미리보기: 같은 mock 과 서버 코드를 앞에 넣는다
const editorPage = read('varedit.html').replace('<body>', '<body>\n<script>' + mock + '</script>\n' +
  serverFiles.map(f => '<script>\n' + read(f) + '\n</script>').join('\n'));
const out2 = path.join(__dirname, 'preview_varedit.html');
fs.writeFileSync(out2, editorPage);
console.log('wrote', out2);

// 변수 계산(compute.html) 미리보기
const computePage = read('compute.html').replace('<body>', '<body>\n<script>' + mock + '</script>\n' +
  serverFiles.map(f => '<script>\n' + read(f) + '\n</script>').join('\n'));
const out3 = path.join(__dirname, 'preview_compute.html');
fs.writeFileSync(out3, computePage);
console.log('wrote', out3);
