/**
 * 예제 데이터 (메뉴 "예제 데이터 채우기", 미리보기, 테스트가 함께 쓴다).
 * 지금까지 만든 모든 분석을 연습·검증할 수 있도록 변수마다 의도한 특성을 넣었다. 60명, 고정 난수라 항상 같은 값이 나온다.
 *
 * 앞의 다섯 열(성별, 집단, 연령, 점수, 만족도)은 처음 예제와 값이 똑같다 (확인사항의 기준 결과가 이 값에 의존하므로 바꾸지 않는다).
 * 뒤의 열은 별도 난수열로 만들어 앞 열의 값에 영향을 주지 않는다 (맨 뒤 일곱 열은 각각 또 다른 난수열이라 그 앞 열도 그대로다).
 *
 * 열과 용도:
 *  성별, 집단, 연령, 점수, 만족도  기존 예제 (점수 999, 만족도 99는 결측값 코드)
 *  근무연수     연령과 거의 같이 움직임 → 상관분석, 회귀의 공선성(VIF)
 *  학력         서열 1~4 → 서열 척도, 값 레이블, Kruskal-Wallis 등
 *  지역         글자 값의 명목 변수(서울/경기/부산/기타) → 교차분석, 집단변수
 *  사전점수, 사후점수  같은 사람의 두 측정 (집단마다 향상 폭이 다름, 일부 빈 칸) → 대응표본 t검정, Wilcoxon
 *  월소득       오른쪽으로 크게 치우침 (9999는 결측값 코드) → 정규성 기각, 비모수·스피어만 추천
 *  불안점수     성별마다 분산이 크게 다름 → 독립표본 t검정에서 "등분산을 가정하지 않음" 행, Levene 기각
 *  스트레스     집단마다 분산이 다름 → 분산분석의 Levene 기각, Welch, Games-Howell
 *  문항1~5      같은 특성을 재는 5점 문항 (문항3, 문항5는 역문항, 99는 결측값 코드) → 역문항 처리, 합계·평균
 *  성과점수     근무연수, 학력, 집단, 문항 특성에 영향받음 → 선형 회귀
 *  합격         성과점수가 높을수록 1 → 교차분석 (로지스틱 회귀 예정)
 *  직무자율성, 직무몰입  근무연수가 직무몰입에 미치는 효과가 자율성이 높을수록 커지고 집단마다 기울기가 다름 → 조절효과 (연속형·범주형 조절변수)
 *  교육시간, 자기효능감, 동기수준, 직무성과  교육시간 → (자기효능감, 동기수준) → 직무성과. 자기효능감은 유의한 매개, 동기수준은 효과가 약함 → 매개효과 (병렬 다중 매개)
 *  이직여부     0/1. 문항 점수(역문항을 뒤집은 평균)가 낮고 근무연수가 짧을수록 1이 되기 쉬움 (집단 차이는 약하고 무작위 잡음이 큼, 합격과 달리 완전 분리가 아님) → 로지스틱 회귀
 */

var SAMPLE_N_ = 60;
var SAMPLE_VAR_ROWS_ = [
  ['성별', '성별', '명목', '1=남, 2=여', ''],
  ['집단', '실험 집단', '명목', '1=A, 2=B, 3=C', ''],
  ['연령', '연령(세)', '연속', '', ''],
  ['점수', '시험 점수', '연속', '', 999],
  ['만족도', '만족도(5점)', '서열', '1=매우 낮음, 5=매우 높음', 99],
  ['근무연수', '근무연수(년)', '연속', '', ''],
  ['학력', '학력', '서열', '1=고졸, 2=전문대졸, 3=대졸, 4=대학원졸', ''],
  ['지역', '거주 지역', '명목', '', ''],
  ['사전점수', '사전 점수', '연속', '', ''],
  ['사후점수', '사후 점수', '연속', '', ''],
  ['월소득', '월소득(만원)', '연속', '', 9999],
  ['불안점수', '불안 점수', '연속', '', ''],
  ['스트레스', '스트레스 점수', '연속', '', ''],
  ['문항1', '업무에 만족한다', '서열', '1=전혀 그렇지 않다, 2=그렇지 않다, 3=보통이다, 4=그렇다, 5=매우 그렇다', 99],
  ['문항2', '동료와 잘 지낸다', '서열', '1=전혀 그렇지 않다, 2=그렇지 않다, 3=보통이다, 4=그렇다, 5=매우 그렇다', 99],
  ['문항3', '업무가 지루하다 (역문항)', '서열', '1=전혀 그렇지 않다, 2=그렇지 않다, 3=보통이다, 4=그렇다, 5=매우 그렇다', 99],
  ['문항4', '성장하고 있다고 느낀다', '서열', '1=전혀 그렇지 않다, 2=그렇지 않다, 3=보통이다, 4=그렇다, 5=매우 그렇다', 99],
  ['문항5', '이직하고 싶다 (역문항)', '서열', '1=전혀 그렇지 않다, 2=그렇지 않다, 3=보통이다, 4=그렇다, 5=매우 그렇다', 99],
  ['성과점수', '성과 점수', '연속', '', ''],
  ['합격', '합격 여부', '명목', '0=불합격, 1=합격', ''],
  ['직무자율성', '직무자율성 점수', '연속', '', ''],
  ['직무몰입', '직무몰입 점수', '연속', '', ''],
  ['교육시간', '교육 이수 시간', '연속', '', ''],
  ['자기효능감', '자기효능감 점수', '연속', '', ''],
  ['동기수준', '동기 수준 점수', '연속', '', ''],
  ['직무성과', '직무 성과 점수', '연속', '', ''],
  ['이직여부', '이직 여부', '명목', '0=잔류, 1=이직', '']
];

/** 예제 데이터와 변수정의: { header, rows, varRows } (rows: 머리글 제외 60행) */
function sampleDataset_() {
  // 앞 다섯 열: 처음 예제와 같은 난수열 (바꾸지 말 것)
  var seed = 20260101;
  function rnd() {
    seed = (seed * 1664525 + 1013904223) % 4294967296;
    return seed / 4294967296;
  }
  // 뒤의 열: 별도 난수열
  var seed2 = 20261025;
  function rnd2() {
    seed2 = (seed2 * 1664525 + 1013904223) % 4294967296;
    return (seed2 + 0.5) / 4294967296;
  }
  function norm2() { return Math.sqrt(-2 * Math.log(rnd2())) * Math.cos(2 * Math.PI * rnd2()); }
  // 맨 뒤 두 열(직무자율성, 직무몰입): 또 다른 난수열
  var seed3 = 20261028;
  function rnd3() {
    seed3 = (seed3 * 1664525 + 1013904223) % 4294967296;
    return (seed3 + 0.5) / 4294967296;
  }
  function norm3() { return Math.sqrt(-2 * Math.log(rnd3())) * Math.cos(2 * Math.PI * rnd3()); }
  // 맨 뒤 네 열(교육시간, 자기효능감, 동기수준, 직무성과): 또 다른 난수열
  var seed4 = 20261029;
  function rnd4() {
    seed4 = (seed4 * 1664525 + 1013904223) % 4294967296;
    return (seed4 + 0.5) / 4294967296;
  }
  function norm4() { return Math.sqrt(-2 * Math.log(rnd4())) * Math.cos(2 * Math.PI * rnd4()); }
  // 맨 뒤 한 열(이직여부): 또 다른 난수열
  var seed5 = 20261031;
  function rnd5() {
    seed5 = (seed5 * 1664525 + 1013904223) % 4294967296;
    return (seed5 + 0.5) / 4294967296;
  }
  function clamp(x, lo, hi) { return Math.max(lo, Math.min(hi, x)); }
  function round1(x) { return Math.round(x * 10) / 10; }

  var header = SAMPLE_VAR_ROWS_.map(function (r) { return r[0]; });
  var rows = [], perf = [], noise = [], tail = [];
  for (var i = 0; i < SAMPLE_N_; i++) {
    var group = 1 + Math.floor(rnd() * 3);
    var score = Math.round(60 + group * 4 + rnd() * 30);
    var row = [
      1 + Math.floor(rnd() * 2),
      group,
      20 + Math.floor(rnd() * 30),
      rnd() < 0.05 ? 999 : score,
      rnd() < 0.08 ? 99 : 1 + Math.floor(rnd() * 5)
    ];
    var sex = row[0], age = row[2];

    var tenure = clamp(Math.round((age - 22) * 0.8 + norm2() * 2), 0, age - 18);
    var u = rnd2();
    var edu = u < 0.15 ? 1 : u < 0.35 ? 2 : u < 0.8 ? 3 : 4;
    u = rnd2();
    var region = u < 0.35 ? '서울' : u < 0.65 ? '경기' : u < 0.85 ? '부산' : '기타';
    var pre = clamp(Math.round(65 + norm2() * 8), 30, 100);
    var post = clamp(Math.round(pre + 1 + (group - 1) * 3.5 + norm2() * 5), 30, 100);
    var preMissing = rnd2() < 0.03, postMissing = rnd2() < 0.05;
    var income = Math.round(Math.exp(5.4 + 0.02 * tenure + 0.55 * norm2()));
    var incomeMissing = rnd2() < 0.06;
    var anxiety = round1(50 + (group - 2) * 2 + norm2() * (sex === 1 ? 4 : 12));
    var stress = round1(40 + norm2() * [4, 8, 16][group - 1]);
    var trait = norm2();
    var items = [];
    for (var k = 1; k <= 5; k++) {
      var raw = clamp(Math.round(3 + 0.9 * trait + 0.7 * norm2()), 1, 5);
      if (k === 3 || k === 5) raw = 6 - raw;                // 역문항
      items.push(rnd2() < 0.05 ? 99 : raw);
    }
    var performance = round1(35 + 1.5 * tenure + 4 * edu + 2 * (group - 2) + 3 * trait + norm2() * 6);
    perf.push(performance);
    noise.push(norm2() * 4);

    // 조절효과 연습용: 근무연수 → 직무몰입의 기울기가 자율성이 높을수록 커지고 집단마다 다르다
    var autonomy = round1(50 + 10 * norm3());
    var tc = tenure - 10, ac = autonomy - 50;
    var engagement = round1(50 + 0.3 * tc + 0.1 * ac + 0.04 * tc * ac + [0, 0.3, 0.6][group - 1] * tc + norm3() * 4);
    // 매개효과 연습용: 교육시간 → (자기효능감, 동기수준) → 직무성과
    var training = clamp(round1(20 + 8 * norm4()), 0, 40);
    var efficacy = round1(30 + 0.5 * training + 6 * norm4());
    var motivation = round1(40 + 0.3 * training + 7 * norm4());
    var jobPerf = round1(20 + 0.2 * training + 0.5 * efficacy + 0.1 * motivation + 5 * norm4());
    // 로지스틱 회귀 연습용: 문항 특성(trait)이 낮고 근무연수가 짧을수록 이직(1)하기 쉽다 (집단 차이는 약하다)
    var quitEta = -0.3 - 1.0 * trait - 0.1 * (tenure - 10) + [0, 0.4, 0.9][group - 1];
    var quit = rnd5() < 1 / (1 + Math.exp(-quitEta)) ? 1 : 0;
    tail.push([autonomy, engagement, training, efficacy, motivation, jobPerf, quit]);

    rows.push(row.concat([
      tenure, edu, region,
      preMissing ? '' : pre, postMissing ? '' : post,
      incomeMissing ? 9999 : income,
      anxiety, stress
    ], items, [performance, 0]));
  }
  // 합격: 성과점수(+약간의 잡음)가 중앙값 이상이면 1 (두 범주가 비슷한 수가 되도록)
  var sorted = perf.slice().sort(function (a, b) { return a - b; });
  var median = (sorted[SAMPLE_N_ / 2 - 1] + sorted[SAMPLE_N_ / 2]) / 2;
  rows.forEach(function (r, idx) { r[r.length - 1] = perf[idx] + noise[idx] >= median ? 1 : 0; });
  rows.forEach(function (r, idx) { tail[idx].forEach(function (v) { r.push(v); }); });
  return { header: header, rows: rows, varRows: SAMPLE_VAR_ROWS_.map(function (r) { return r.slice(); }) };
}

if (typeof module !== 'undefined') {
  module.exports = { sampleDataset_: sampleDataset_, SAMPLE_N_: SAMPLE_N_, SAMPLE_VAR_ROWS_: SAMPLE_VAR_ROWS_ };
}
