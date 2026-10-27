/**
 * 메뉴 "Stat4SW 통계분석": 자주 쓰는 것(분석 창, 변수 설정·계산, 결과 정리, 예제 데이터)은 위에 두고,
 * 고급 분석 켜기/끄기와 관리 기능(시험 데이터)은 하위 메뉴 "설정·관리"로 묶는다 (시험 데이터는 고급 분석을 켠 뒤에만 실행된다).
 * 새 작업을 시작하는 "새 작업 시작 (데이터 모두 지우기)"는 초급 사용자도 쓰므로 위에 두고, 지우기 전에 "지우기"를 직접 입력받는다.
 */
function onOpen() {
  var ui = SpreadsheetApp.getUi();
  ui.createMenu('Stat4SW 통계분석')
    .addItem('분석 창 열기', 'openAnalysis')
    .addSeparator()
    .addItem('변수 설정 편집', 'openVarEditor')
    .addItem('변수 계산 (역문항, 합계·평균, z점수)', 'openCompute')
    .addItem('결과 탭 지우기 (거리·유사도 행렬 탭 포함)', 'clearResultSheet')
    .addItem('새 작업 시작 (데이터 모두 지우기)', 'startNewWork')
    .addItem('예제 데이터 채우기', 'fillSampleData')
    .addItem('분석 선택 도우미 (넓게 보기)', 'openAdvisor')
    .addSeparator()
    .addSubMenu(ui.createMenu('설정·관리')
      .addItem('고급 분석 켜기/끄기', 'toggleAdvancedMode')
      .addItem('시험 데이터 (기법별 점검용)', 'menuOpenTestData'))
    .addSeparator()
    .addItem('사용 안내', 'showHelp')
    .addItem('정보·저작권', 'showAbout')
    .addToUi();
}

/** 개인 사본 안내: 사용자마다 처음 분석 창을 열 때 한 번만 보여 준다 (데이터는 이 시트에만 있고, 시트를 공유하면 데이터도 공유된다) */
var PRIVACY_NOTICE_KEY_ = 'stat4sw.privacyNoticed';
function privacyNoticeText_() {
  return '이 시트는 당신의 개인 사본입니다.\n\n' +
    '· 분석할 데이터는 이 시트에만 저장되고 외부로 보내지 않습니다.\n' +
    '· 시트 링크를 다른 사람에게 공유하면 데이터도 함께 보입니다. 공유가 필요하면 익명화한 자료만 쓰세요.\n' +
    '· 다른 사람의 시트에 편집자로 들어가 함께 작업하지 말고, 각자 사본을 만들어 쓰세요.\n\n' +
    '(이 안내는 한 번만 보입니다. 메뉴 "정보·저작권"과 "사용 안내"에서 다시 볼 수 있습니다.)';
}
function noticePrivacyOnce_() {
  try {
    var props = PropertiesService.getUserProperties();
    if (props.getProperty(PRIVACY_NOTICE_KEY_) === '1') return;
    SpreadsheetApp.getUi().alert('개인 사본 안내', privacyNoticeText_(), SpreadsheetApp.getUi().ButtonSet.OK);
    props.setProperty(PRIVACY_NOTICE_KEY_, '1');
  } catch (e) {
    // 안내를 못 보여 줘도 분석 창은 열려야 한다
  }
}

function openAnalysis() {
  noticePrivacyOnce_();
  var html = HtmlService.createTemplateFromFile('index')
    .evaluate()
    .setWidth(1100)
    .setHeight(700);
  SpreadsheetApp.getUi().showModelessDialog(html, 'Stat4SW');
}

function openVarEditor() {
  var html = HtmlService.createHtmlOutputFromFile('varedit').setWidth(680).setHeight(640);
  SpreadsheetApp.getUi().showModalDialog(html, '변수 설정 편집');
}

function openCompute() {
  var html = HtmlService.createHtmlOutputFromFile('compute').setWidth(700).setHeight(680);
  SpreadsheetApp.getUi().showModalDialog(html, '변수 계산');
}

/** 어떤 분석을 쓸지 넓은 화면에서 고르는 도우미. 고른 분석은 분석 창에서 (변수가 채워진 채) 열린다. */
function openAdvisor() {
  var html = HtmlService.createHtmlOutputFromFile('advisor').setWidth(1000).setHeight(720);
  SpreadsheetApp.getUi().showModalDialog(html, '분석 선택 도우미');
}

function openTestData() {
  var html = HtmlService.createHtmlOutputFromFile('testdata').setWidth(760).setHeight(680);
  SpreadsheetApp.getUi().showModalDialog(html, '시험 데이터');
}

function showHelp() {
  var html = HtmlService.createTemplateFromFile('help').evaluate().setWidth(560).setHeight(560);
  SpreadsheetApp.getUi().showModalDialog(html, 'Stat4SW 사용 안내');
}

function showAbout() {
  var html = HtmlService.createTemplateFromFile('about').evaluate().setWidth(640).setHeight(680);
  SpreadsheetApp.getUi().showModalDialog(html, 'Stat4SW 정보·저작권');
}

/** HTML 템플릿에서 <?!= include('파일명') ?> 로 다른 html 파일을 끼워 넣는다. */
function include(filename) {
  return HtmlService.createHtmlOutputFromFile(filename).getContent();
}
