/**
 * 일반/고급 분석 모드. 초급 사용자에게는 자주 쓰는 기본 분석만 보이게 하고, 메뉴 "설정·관리 > 고급 분석 켜기/끄기"로 나머지(고급 분석)를 보이게 한다.
 * 비밀번호는 쓰지 않는다 (스크립트를 볼 수 있는 사람은 어차피 우회할 수 있어 보안이 아니라 화면을 단순하게 하는 장치이다).
 * 켜짐 여부는 사용자별(PropertiesService 사용자 속성)로 기억한다. 관리 기능(시험 데이터)은 켠 뒤에만 실행된다. 서버의 runAnalysis 는 모드와 상관없이 모든 분석을 실행한다 (화면에서만 숨긴다).
 *
 * 분석을 추가하면 아래 두 목록 중 한 곳에 키를 넣어야 한다 (tests/mode.test.js 가 client.html 의 모든 분석이 분류됐는지 확인한다).
 */

/** 기본(초급) 분석: 자주 쓰고 해석이 쉬운 것 */
var BASIC_ANALYSES_ = [
  'frequencies', 'descriptives', 'normality', 'levene', 'binom',
  'ttest1', 'ttestInd', 'ttestPaired', 'anova',
  'mannWhitney', 'wilcoxonPaired', 'kruskalWallis', 'mcnemar',
  'correlation', 'crosstab', 'regression', 'reliability', 'ipa'
];

/** 고급 분석: 여러 변수·가정·해석이 필요한 것 */
var ADVANCED_ANALYSES_ = [
  'wilcoxonOne', 'friedman', 'cochran', 'anova2', 'rmanova',
  'hierarchical', 'logistic', 'ologit', 'moderation', 'mediation',
  'cart', 'efa', 'kappa', 'icc', 'kmeans', 'distance', 'psm', 'km', 'power'
];

var ADVANCED_MODE_KEY_ = 'stat4sw.advanced';
var ADVANCED_MENU_HINT_ = 'Stat4SW 통계분석 > 설정·관리 > 고급 분석 켜기/끄기';

/** 고급 분석이 켜져 있는지 (사용자별). 읽을 수 없으면 꺼짐 */
function isAdvancedMode_() {
  try {
    return PropertiesService.getUserProperties().getProperty(ADVANCED_MODE_KEY_) === '1';
  } catch (e) {
    return false;
  }
}

function setAdvancedMode_(on) {
  PropertiesService.getUserProperties().setProperty(ADVANCED_MODE_KEY_, on ? '1' : '0');
}

/** 분석 창·도우미에 알려 줄 모드 정보 */
function getModeInfo_() {
  return { advanced: isAdvancedMode_(), advancedKeys: ADVANCED_ANALYSES_.slice() };
}

/** 분석 창·도우미가 데이터를 읽지 못해도(데이터 탭이 없거나 비어 있을 때) 모드를 알 수 있게 따로 부르는 함수 */
function getModeInfo() {
  return getModeInfo_();
}

/** 분석 창 HTML(index.html)이 열릴 때 끼워 넣는 모드 정보 (JSON 문자열). 첫 화면부터 목록이 걸러진다 */
function modeBootJson_() {
  return JSON.stringify(getModeInfo_()).replace(/</g, '\\u003c');
}

/** 메뉴: 고급 분석 켜기/끄기 */
function toggleAdvancedMode() {
  var ui = SpreadsheetApp.getUi();
  var next = !isAdvancedMode_();
  setAdvancedMode_(next);
  ui.alert('고급 분석 ' + (next ? '켜짐' : '꺼짐'),
    next
      ? '고급 분석 ' + ADVANCED_ANALYSES_.length + '개가 분석 목록에 나타납니다. 이미 열려 있는 분석 창은 닫았다가 다시 열거나 새로고침하세요.\n(시험 데이터도 이제 쓸 수 있습니다.)'
      : '분석 목록에 기본 분석 ' + BASIC_ANALYSES_.length + '개만 보입니다. 이미 열려 있는 분석 창은 닫았다가 다시 열거나 새로고침하세요.',
    ui.ButtonSet.OK);
}

/** 관리 기능은 고급 분석을 켠 뒤에만 쓸 수 있다. 꺼져 있으면 안내하고 false. why: 잠그는 이유(선택) */
function requireAdvanced_(what, why) {
  if (isAdvancedMode_()) return true;
  SpreadsheetApp.getUi().alert(what, (why ? why + '\n' : '') + '이 기능은 고급 분석을 켠 뒤에 쓸 수 있습니다.\n메뉴 ' + ADVANCED_MENU_HINT_ + '를 누르세요.', SpreadsheetApp.getUi().ButtonSet.OK);
  return false;
}

/** 메뉴 항목: 시험 데이터 (고급 분석을 켠 뒤에만) */
function menuOpenTestData() {
  if (requireAdvanced_('시험 데이터')) openTestData();
}

if (typeof module !== 'undefined') {
  module.exports = { BASIC_ANALYSES_: BASIC_ANALYSES_, ADVANCED_ANALYSES_: ADVANCED_ANALYSES_ };
}
