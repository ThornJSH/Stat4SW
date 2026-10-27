var APP_VERSION = '1.1.1';

// 새 버전과 업데이트 방법을 알리는 페이지(GitHub 릴리스) 주소. 비어 있으면 정보 화면에 링크 대신 안내 문구가 나온다.
var RELEASES_URL_ = 'https://github.com/ThornJSH/Stat4SW/releases/latest';
var SHEET_DATA = '데이터';
var SHEET_VARS = '변수정의';
var SHEET_RESULT = '결과';
// 분석할 수 있는 최대 사례(행) 수. 근거: 수치 정확도는 20만 건 이상에서도 확인했고, 한계는 실행 시간(6분)과 시트 읽기 시간이다.
// Apps Script에서의 실제 속도는 편집기에서 benchmarkAnalyses 를 실행해 확인하고, 필요하면 이 값을 조정한다.
var MAX_ROWS = 50000;
var VARS_HEADER = ['변수명', '라벨', '척도', '값 레이블', '결측값 코드'];
