/**
 * 분석 선택 도우미. 두 가지 방식으로 "어떤 분석을 쓸지"를 알려 준다.
 *  1) 질문 나무(ADVISOR_TREE_): 목적 → 집단 수 → 결과 변수 유형을 객관식으로 묻고 추천 분석과 비모수 대안을 보여 준다.
 *     화면(client.html)이 getAdvisor()로 나무를 한 번 받아 서버 왕복 없이 따라간다.
 *  2) 변수 기반(adviseFromVariables): 분석 창에서 체크한 변수의 척도·집단 수로 바로 추천하고 변수를 역할 칸에 넣을 정보를 돌려준다.
 *
 * 규칙은 Stat_Assistant 의 LogicDB 구조(목적, 집단 수, 변수 유형)를 참고하되 새로 만들었다 (표본 크기 N≥30 한 줄 규칙은 쓰지 않는다).
 * 모수·비모수 선택은 분석 결과의 가정 점검(Stats_Assumptions.js)이 정규성을 보고 알려 주므로, 여기서는 표본이 작을 때의 안내만 덧붙인다.
 * 이 도구에 없는 분석은 { name, unsupported: true, fallback } 으로 표시하고 쓸 수 있는 가까운 분석을 안내한다.
 */

var ADVISOR_SMALL_N_ = 30;
var ADVISOR_SMALL_NOTE_ = '사례가 적습니다(30건 미만). 정규성이 의심되면 비모수 대안을 함께 보세요. 분석 결과의 "가정 점검"이 정규성을 판단해 알려 줍니다.';
var ADVISOR_MAX_GROUP_VALUES_ = 12;       // 변수 기반 추천에서 집단변수로 보는 최대 범주 수

/** 추천 분석: key 는 client.html 의 분석 키. parametric 이면 표본이 작을 때 비모수 대안 안내가 붙는다. */
function advRec_(key, why, parametric) { return { key: key, why: why, parametric: !!parametric }; }
function advAlt_(key, why) { return { key: key, why: why }; }
/** 이 도구에 없는 분석 */
function advMissing_(name, why, fallback) { return { name: name, unsupported: true, fallback: fallback, why: why }; }
function advLeaf_(rec, alts, notes) { return { rec: rec, alt: alts || [], notes: notes || [] }; }
function advAsk_(q, choices) { return { q: q, choices: choices }; }
function advChoice_(label, next) { return { label: label, next: next }; }

function advBuildTree_() {
  var KIND_NOM = '범주 (성별, 합격 여부처럼 값이 몇 가지로 나뉨)';
  var KIND_ORD = '서열·리커트 (5점 척도 한 문항처럼 순서만 의미)';
  var KIND_NUM = '숫자 (점수, 나이처럼 크기가 의미, 여러 문항의 합계·평균 포함)';
  var ordNote = '여러 문항을 합산·평균한 점수는 숫자로 보고 평균 비교(t검정, ANOVA)를 써도 됩니다.';

  var describe = advAsk_('변수는 어떤 종류인가요?', [
    advChoice_(KIND_NOM, advLeaf_(advRec_('frequencies', '범주마다 몇 건인지와 비율을 센다.'), [], [])),
    advChoice_(KIND_NUM, advLeaf_(advRec_('descriptives', '평균, 표준편차, 최솟값·최댓값 등으로 요약한다.'),
      [advAlt_('normality', '값이 정규분포에 가까운지 확인한다 (평균 비교를 하기 전에).')], []))
  ]);

  var compare2 = advAsk_('비교할 결과(종속) 변수는 어떤 종류인가요?', [
    advChoice_(KIND_NUM, advLeaf_(advRec_('ttestInd', '두 집단의 평균이 다른지 본다.', true),
      [advAlt_('mannWhitney', '정규분포가 아니거나 표본이 작을 때 순위로 비교한다.'), advAlt_('levene', '두 집단의 흩어진 정도(분산)가 같은지 먼저 확인한다.')], [])),
    advChoice_(KIND_ORD, advLeaf_(advRec_('mannWhitney', '순서만 의미 있는 점수는 순위로 두 집단을 비교한다.'), [advAlt_('ttestInd', '여러 문항의 합계·평균 점수라면 평균 비교도 가능하다.')], [ordNote])),
    advChoice_(KIND_NOM, advLeaf_(advRec_('crosstab', '집단별 비율이 다른지 교차표와 카이제곱으로 본다.'), [], ['2×2 표에서 기대 빈도가 작으면 결과의 Fisher 정확 검정을 보세요.']))
  ]);
  var compare3 = advAsk_('비교할 결과(종속) 변수는 어떤 종류인가요?', [
    advChoice_(KIND_NUM, advLeaf_(advRec_('anova', '세 집단 이상의 평균이 다른지 보고, 다르면 사후검정으로 어느 집단 사이인지 본다.', true),
      [advAlt_('kruskalWallis', '정규분포가 아니거나 표본이 작을 때 순위로 비교한다.'), advAlt_('levene', '집단들의 분산이 같은지 먼저 확인한다.')], [])),
    advChoice_(KIND_ORD, advLeaf_(advRec_('kruskalWallis', '순서만 의미 있는 점수는 순위로 세 집단 이상을 비교한다.'), [advAlt_('anova', '여러 문항의 합계·평균 점수라면 평균 비교도 가능하다.')], [ordNote])),
    advChoice_(KIND_NOM, advLeaf_(advRec_('crosstab', '집단별 비율이 다른지 교차표와 카이제곱으로 본다.'), [], []))
  ]);
  var compareGroups = advAsk_('비교할 집단은 몇 개인가요?', [advChoice_('2개 (남·여, 참여·비참여)', compare2), advChoice_('3개 이상 (연령대, 지역)', compare3),
    advChoice_('집단을 나누는 기준이 2가지 (예: 성별 × 프로그램)', advLeaf_(advRec_('anova2', '두 요인의 주효과와 상호작용(조합 효과)을 한 번에 검정한다. 종속변수는 숫자(점수)여야 한다.', true), [advAlt_('anova', '기준을 하나만 볼 때는 일원 분산분석을 쓴다.')], ['요인 두 개는 각각 범주형이고 모든 조합(셀)에 사례가 있어야 합니다.', '상호작용이 유의하면 주효과보다 단순 주효과와 그림으로 해석하세요.']))]);

  var repeatedNote = '변수를 측정한 순서대로 보내세요 (사전, 중간, 사후). 한 사람이 모든 측정에 값이 있어야 합니다.';
  var paired3 = advAsk_('측정한 결과는 어떤 종류인가요?', [
    advChoice_(KIND_NUM, advLeaf_(advRec_('rmanova', '같은 사람의 3회 이상 측정 평균이 달라졌는지 보고, 구형성 검정과 교정, 쌍별 비교까지 한다.', true), [advAlt_('friedman', '정규분포가 아니거나 표본이 작을 때 순위로 비교한다.')], [repeatedNote, '실험군·대조군처럼 집단이 있어 집단에 따라 변화가 다른지 보려면 "집단변수" 칸에 집단을 넣으세요 (집단 × 측정 혼합 설계).'])),
    advChoice_(KIND_ORD, advLeaf_(advRec_('friedman', '순서만 의미 있는 점수의 3회 이상 반복 측정을 순위로 비교하고, 달라졌다면 어느 측정 사이인지 사후비교로 본다.'), [], [repeatedNote, ordNote])),
    advChoice_(KIND_NOM, advLeaf_(advRec_('cochran', '같은 사람의 3회 이상 예/아니오 반응의 성공 비율이 달라졌는지 Cochran Q 검정으로 보고, 쌍별 McNemar 검정으로 어느 측정 사이가 다른지 본다.'), [advAlt_('friedman', '0/1로 코딩한 자료에 Friedman 검정을 하면 Cochran의 Q와 같은 값이 나온다.')], [repeatedNote, '분석 창에서 성공으로 볼 값(예: 1)을 정하세요. 성공이 아닌 응답은 모두 실패로 처리합니다.']))
  ]);
  var paired2 = advAsk_('측정한 결과는 어떤 종류인가요?', [
    advChoice_(KIND_NUM, advLeaf_(advRec_('ttestPaired', '같은 사람의 사전·사후 평균이 달라졌는지 본다.', true), [advAlt_('wilcoxonPaired', '정규분포가 아니거나 표본이 작을 때 순위로 비교한다.'), advAlt_('rmanova', '실험군·대조군처럼 집단이 있어 집단에 따라 변화가 다른지 보려면 반복측정 분산분석의 "집단변수" 칸에 집단을 넣는다 (집단 × 측정 혼합 설계).')], [])),
    advChoice_(KIND_ORD, advLeaf_(advRec_('wilcoxonPaired', '순서만 의미 있는 점수의 사전·사후 변화를 순위로 본다.'), [advAlt_('ttestPaired', '여러 문항의 합계·평균 점수라면 평균 비교도 가능하다.')], [ordNote])),
    advChoice_(KIND_NOM, advLeaf_(advRec_('mcnemar', '같은 사람의 사전·사후 예/아니오 변화(예: 취업 여부)가 달라졌는지, 바뀐 쌍만으로 본다.'), [advAlt_('crosstab', '전·후 응답이 어떻게 겹치는지 교차표로만 볼 수도 있다 (변화의 유의성은 McNemar가 맞다).')], ['전·후 두 변수를 각각 한 열로 두고, 두 변수의 값이 합쳐서 2가지(0/1, 아니오/예)여야 합니다.']))
  ]);
  var paired = advAsk_('같은 사람을 몇 번 측정했나요?', [advChoice_('2번 (사전·사후)', paired2), advChoice_('3번 이상 (사전·중간·사후 등)', paired3)]);

  var onesample = advAsk_('측정한 변수는 어떤 종류인가요?', [
    advChoice_(KIND_NUM, advLeaf_(advRec_('ttest1', '한 집단의 평균이 기준값(목표치, 전국 평균)과 다른지 본다.', true), [advAlt_('wilcoxonOne', '정규분포가 아니거나 표본이 작을 때 중앙값으로 비교한다.')], [])),
    advChoice_(KIND_ORD, advLeaf_(advRec_('wilcoxonOne', '순서만 의미 있는 점수가 기준값과 다른지 순위로 본다.'), [], [ordNote])),
    advChoice_(KIND_NOM, advLeaf_(advRec_('binom', '두 가지로 나뉘는 범주(예: 합격/불합격)의 성공 비율이 기준 비율과 다른지 정확한 이항 분포로 검정하고 비율의 신뢰구간도 구한다.'), [advAlt_('frequencies', '검정 없이 범주마다 건수와 비율만 확인한다.')], ['분석 창에서 성공으로 볼 값(예: 1)과 검정 비율(예: 0.5)을 정하세요. 성공이 아닌 값은 모두 실패로 처리합니다.']))
  ]);

  var relation = advAsk_('두 변수는 어떤 종류인가요?', [
    advChoice_('둘 다 숫자', advLeaf_(advRec_('correlation', '두 숫자 변수가 함께 늘거나 줄어드는 정도(피어슨 r)를 본다.', true), [], ['정규분포가 아니거나 이상치가 있으면 옵션에서 스피어만 순위상관을 함께 켜세요.', '연령처럼 다른 변수의 영향을 빼고 보려면 "통제변수" 칸에 그 변수를 넣으세요 (편상관).'])),
    advChoice_('서열·리커트가 하나라도 있다', advLeaf_(advRec_('correlation', '순서 정보로 연관을 보는 스피어만 순위상관을 쓴다.'), [], ['옵션에서 "스피어만 순위상관"을 켜세요. 표본이 작거나 같은 값이 많으면(5점 척도 등) "켄달 타우-b"도 함께 켜서 확인하세요.'])),
    advChoice_('둘 다 범주', advLeaf_(advRec_('crosstab', '두 범주형 변수가 관련 있는지 교차표와 카이제곱, 효과크기로 본다.'), [], [])),
    advChoice_('한쪽은 숫자, 한쪽은 범주', compareGroups)
  ]);

  var predictNum = advAsk_('무엇을 알고 싶나요?', [
    advChoice_('여러 변수가 결과에 미치는 영향과 설명력', advLeaf_(advRec_('regression', '독립변수들이 숫자 결과에 미치는 영향과 설명력(R²)을 본다.'), [advAlt_('hierarchical', '통제변수를 먼저 넣고 핵심 변수가 설명력을 얼마나 더 높이는지 단계별로 본다.'), advAlt_('cart', '규칙(조건)으로 결과를 예측하는 나무를 만든다.')], [])),
    advChoice_('통제변수를 먼저 넣고, 핵심 변수의 추가 설명력을 단계별로 보고 싶다', advLeaf_(advRec_('hierarchical', '독립변수를 1~4단계로 나눠 차례로 넣고, 단계마다 R²가 유의하게 늘었는지 본다.'), [advAlt_('regression', '단계 구분 없이 한 번에 모두 넣고 영향을 본다.')], ['1단계에는 통제변수(성별, 연령 등)를, 2단계 이후에는 확인하려는 변수를 보내세요.'])),
    advChoice_('영향이 중간 변수를 거쳐 전달되는지 (매개)', advLeaf_(advRec_('mediation', 'X가 매개변수를 거쳐 Y에 미치는 간접효과를 본다.'), [], [])),
    advChoice_('영향이 조건(집단, 수준)에 따라 달라지는지 (조절)', advLeaf_(advRec_('moderation', 'X가 Y에 미치는 효과가 조절변수의 값에 따라 달라지는지 본다.'), [], []))
  ]);
  var predict = advAsk_('결과(예측하려는 변수)는 어떤 종류인가요?', [
    advChoice_('숫자 (점수, 금액)', predictNum),
    advChoice_('두 가지 범주 (합격/불합격, 취업 여부)', advLeaf_(advRec_('logistic', '어떤 변수가 사건이 일어날 가능성을 높이는지 오즈비로 본다.'), [advAlt_('cart', '규칙(조건)으로 분류하는 나무를 만든다.')], ['0/1 열을 미리 만들 필요 없이 분석 창에서 사건으로 볼 값을 고를 수 있습니다.'])),
    advChoice_('세 가지 이상 범주', advLeaf_(advRec_('cart', '규칙으로 여러 범주를 분류하는 분류 나무를 만든다.'), [], ['다항 로지스틱 회귀는 아직 지원하지 않습니다.'])),
    advChoice_('서열·리커트 한 문항', advLeaf_(advRec_('ologit', '순서가 있는 범주(예: 5점 척도 한 문항) 결과의 영향 요인을 서열 로지스틱(비례오즈)으로 본다. B가 양수이면 변수가 클수록 더 높은 범주가 될 가능성이 큽니다.'), [advAlt_('regression', '여러 문항의 합계·평균 점수처럼 거의 연속인 점수라면 선형 회귀도 쓸 수 있다.'), advAlt_('cart', '규칙으로 분류하는 나무를 만든다.')], ['결과의 값은 3~10가지여야 하고, 비례오즈 가정 점검(Brant) 결과를 꼭 확인하세요.']))
  ]);

  var agreement = advAsk_('평가한 결과는 어떤 종류인가요?', [
    advChoice_(KIND_NUM, advLeaf_(advRec_('icc', '평가자(또는 반복 측정) 사이에서 점수가 얼마나 일치하는지 급내상관계수(ICC)로 본다.'), [], ['평가자마다 변수 한 열로 두고 대상(사례)은 행입니다. 어느 ICC 값을 볼지는 결과 아래 안내를 보세요.'])),
    advChoice_(KIND_ORD, advLeaf_(advRec_('kappa', '서열 범주(예: 5점 척도)의 평가 일치를 가중 kappa로 본다.'), [advAlt_('icc', '점수를 숫자로 보아 급내상관계수(ICC)로 볼 수도 있다.')], ['평가자가 2명이면 옵션의 "이차 가중"을 고르세요. 평가자마다 변수 한 열로 두세요.'])),
    advChoice_(KIND_NOM, advLeaf_(advRec_('kappa', '범주 평가(예: 있음/없음)의 일치도를 Cohen kappa(평가자 3명 이상은 Fleiss kappa)로 본다.'), [], ['평가자마다 변수 한 열로 두세요. 한쪽 범주에 평가가 몰리면 일치율도 함께 보고하세요.']))
  ]);

  var cluster = advAsk_('무엇을 하고 싶나요?', [
    advChoice_('비슷한 사례끼리 몇 개의 집단으로 묶고 싶다', advLeaf_(advRec_('kmeans', '숫자 변수로 사례를 정해 둔 개수의 군집으로 묶는다.'), [], ['결과의 군집 번호를 저장해 집단 비교에 이어 쓸 수 있습니다.'])),
    advChoice_('사례(또는 변수) 사이가 얼마나 다른지 재고 싶다', advLeaf_(advRec_('distance', '유클리드, 코사인, 상관 등의 거리·유사도 행렬을 만든다.'), [], [])),
    advChoice_('이상한 사례(다변량 이상치)를 찾고 싶다', advLeaf_(advRec_('distance', '옵션 "마할라노비스 이상치 점검"으로 변수 조합이 이상한 사례를 찾는다.'), [], []))
  ]);

  var treatment = advAsk_('처치(프로그램 참여 등)를 받을 사람을 무작위로 정했나요?', [
    advChoice_('아니오 (참여자가 스스로 선택했거나 기존 집단)', advLeaf_(advRec_('psm', '처치군과 대조군이 처음부터 다른 점을 성향점수로 맞춘 뒤 효과를 비교한다.'), [], ['처치 변수는 0/1(값 2가지), 공변량은 처치 선택에 영향을 줄 변수들입니다. 관찰하지 못한 변수는 보정되지 않습니다.'])),
    advChoice_('예 (무작위 배정)', compareGroups)
  ]);

  return advAsk_('무엇을 알고 싶나요?', [
    advChoice_('값이 어떻게 분포하는지 보고 싶다 (빈도, 평균, 정규성)', describe),
    advChoice_('집단끼리 비교하고 싶다 (남·여, 참여·비참여, 연령대)', compareGroups),
    advChoice_('같은 사람의 사전·사후 변화를 보고 싶다', paired),
    advChoice_('한 집단을 기준값(목표치, 평균)과 비교하고 싶다', onesample),
    advChoice_('두 변수가 서로 관련 있는지 보고 싶다', relation),
    advChoice_('무엇이 결과에 영향을 주는지 알고 싶다 (영향 요인, 예측)', predict),
    advChoice_('설문 문항들이 같은 것을 일관되게 재는지 확인하고 싶다', advLeaf_(advRec_('reliability', '문항 사이 일관성(크론바흐 α)과 문항을 빼면 α가 어떻게 되는지 본다.'), [], ['문항은 같은 방향으로 점수가 매겨져야 합니다. 거꾸로 묻는 문항은 변수 계산의 역문항 처리를 먼저 하세요.'])),
    advChoice_('여러 문항이 몇 개의 묶음(요인)으로 나뉘는지 알고 싶다', advLeaf_(advRec_('efa', '문항 사이의 상관 구조에서 몇 개의 요인이 있는지, 각 문항이 어느 요인에 속하는지 찾는다.'), [advAlt_('reliability', '이미 묶음이 정해져 있다면 묶음별로 문항 일관성(크론바흐 α)을 확인한다.')], ['사례 수는 문항 수의 5배 이상(100건 이상이면 좋음)이어야 안정적입니다.', '요인을 찾은 뒤 묶음별로 신뢰도 분석을 이어서 하세요.'])),
    advChoice_('두 사람(또는 두 번)의 평가가 얼마나 일치하는지 알고 싶다 (평가자 간 일치도, 검사-재검사)', agreement),
    advChoice_('사건(사망·재발·이탈 등)이 일어날 때까지의 시간을 보고 싶다 (생존분석, 중도절단 있음)', advLeaf_(advRec_('km', '사건까지의 시간을 중도절단을 반영해 생존곡선으로 나타내고, 집단변수를 넣으면 log-rank 검정으로 집단 사이를 비교한다.'), [advAlt_('ttestInd', '중도절단이 없고 모두 사건을 겪었다면 시간의 평균을 비교할 수도 있지만 분포가 치우치기 쉽다.')], ['시간 변수와 사건 변수(사건 여부)가 필요합니다. 사건이 일어나지 않고 관찰이 끝난 사례는 사건 변수에서 사건값이 아닌 값으로 두세요(중도절단).', '위험이 시간에 따라 비례한다는 가정(비례위험)이 있어 곡선이 교차하면 log-rank 결과를 조심해서 보세요.'])),
    advChoice_('서비스·프로그램에서 무엇을 먼저 개선할지 알고 싶다 (문항별 중요도와 만족도 비교)', advLeaf_(advRec_('ipa', '문항마다 중요도 평균과 수행도(만족도) 평균을 기준선으로 4사분면에 나눠 집중 개선, 현상 유지, 낮은 우선순위, 과잉 투자 문항을 찾는다.'), [advAlt_('ttestPaired', '한 문항의 중요도와 만족도가 다른지만 보려면 대응표본 t검정을 쓴다.')], ['중요도 문항들과 수행도(만족도) 문항들을 같은 개수, 같은 순서로 짝지어 보내세요.', '기준선 방식과 이동량은 결과를 바꾸므로 보고서에 어떤 기준선을 썼는지 적으세요.'])),
    advChoice_('연구를 시작하기 전에 몇 명이 필요한지, 정해진 인원의 검정력을 알고 싶다 (표본 크기, 검정력)', advLeaf_(advRec_('power', '검정 종류와 예상 효과크기를 넣으면 필요한 표본 크기, 검정력, 탐지 가능한 효과크기를 구한다 (데이터 불필요).'), [], ['효과크기는 선행연구의 값이나 의미 있다고 볼 최소한의 차이로 정하세요. 모르면 비워 두면 중간 효과로 계산합니다.', '탈락이 예상되면 세부 옵션의 예상 탈락률을 넣으면 모집 인원이 나옵니다.'])),
    advChoice_('사례를 묶거나 이상한 사례·거리를 보고 싶다', cluster),
    advChoice_('프로그램(처치)의 효과를 비교하고 싶다', treatment)
  ]);
}

var ADVISOR_TREE_ = advBuildTree_();

/** 분석 키 → 이름 (client.html 의 분석 이름과 같아야 한다: 테스트가 확인). 넓게 보는 도우미 대화상자가 쓴다. */
var ADVISOR_LABELS_ = {
  frequencies: '빈도분석',
  descriptives: '기술통계',
  normality: '정규성 검정 (Shapiro-Wilk, Jarque-Bera)',
  levene: '등분산 검정 (Levene)',
  ttest1: '일표본 t검정',
  ttestInd: '독립표본 t검정',
  ttestPaired: '대응표본 t검정',
  anova: '일원 분산분석 (ANOVA)',
  anova2: '이원 분산분석 (요인 2개)',
  mannWhitney: 'Mann-Whitney U 검정 (비모수, 2집단)',
  wilcoxonPaired: 'Wilcoxon 부호순위 검정 (비모수, 대응표본)',
  wilcoxonOne: 'Wilcoxon 부호순위 검정 (비모수, 일표본)',
  mcnemar: 'McNemar 검정 (비모수, 전·후 예/아니오)',
  friedman: 'Friedman 검정 (비모수, 3회 이상 반복측정)',
  rmanova: '반복측정 분산분석 (집단변수는 선택)',
  kruskalWallis: 'Kruskal-Wallis 검정 (비모수, 3집단 이상)',
  correlation: '상관분석',
  crosstab: '교차분석 (카이제곱)',
  reliability: '신뢰도 분석 (크론바흐 α)',
  efa: '탐색적 요인분석 (EFA)',
  kappa: '평가자 일치도 (kappa)',
  icc: '급내상관 (ICC)',
  logistic: '로지스틱 회귀 (이분형)',
  moderation: '조절효과 분석 (상호작용)',
  mediation: '매개효과 분석 (간접효과)',
  distance: '거리·유사도 (마할라노비스 이상치)',
  cart: '의사결정나무 (CART)',
  kmeans: 'k-평균 군집분석',
  psm: '성향점수 매칭 (처치 효과 추정)',
  regression: '선형 회귀',
  ologit: '서열 로지스틱 회귀 (비례오즈)',
  power: '검정력·표본 크기 계산',
  km: 'Kaplan–Meier 생존분석 (log-rank)',
  binom: '이항 검정 (성공 비율 검정)',
  cochran: 'Cochran Q 검정 (예/아니오, 3회 이상 반복측정)',
  ipa: 'IPA 중요도–수행도 분석',
  hierarchical: '위계적 회귀 (단계별 변수 추가)'
};

/** 대화상자용: 질문 나무, 분석 이름, 표본 크기 안내 문구 */
function getAdvisor() {
  // Mode.js 의 전역 함수 (Node 단위 테스트에서는 없으므로 꺼짐으로 본다)
  var mode = typeof getModeInfo_ === 'function' ? getModeInfo_() : { advanced: false, advancedKeys: [] };
  var hint = typeof ADVANCED_MENU_HINT_ === 'string' ? ADVANCED_MENU_HINT_ : '';
  return { tree: ADVISOR_TREE_, labels: ADVISOR_LABELS_, smallN: ADVISOR_SMALL_N_, smallNote: ADVISOR_SMALL_NOTE_, advanced: mode.advanced, advancedKeys: mode.advancedKeys, advancedHint: hint };
}

var ADVISOR_PENDING_KEY_ = 'stat4sw.pendingAdvice';
var ADVISOR_PENDING_TTL_MS_ = 2 * 60 * 1000;

/**
 * 넓게 보는 도우미 대화상자에서 호출: 추천한 분석(과 변수 배치)을 저장하고 분석 창을 연다.
 * 분석 창이 열릴 때 takePendingAdvice 로 받아 그 분석을 고르고 변수를 역할 칸에 넣는다.
 */
function openAnalysisWith(key, roles) {
  try {
    if (typeof key !== 'string' || !/^[A-Za-z0-9]{1,30}$/.test(key)) throw new Error('알 수 없는 분석입니다.');
    var clean = {};
    if (roles && typeof roles === 'object') {
      Object.keys(roles).slice(0, 8).forEach(function (rk) {
        if (!/^[A-Za-z0-9]{1,30}$/.test(rk) || !Array.isArray(roles[rk])) return;
        clean[rk] = roles[rk].slice(0, 60).filter(function (n) { return typeof n === 'string' && n.length <= 100; });
      });
    }
    PropertiesService.getUserProperties().setProperty(ADVISOR_PENDING_KEY_, JSON.stringify({ key: key, roles: clean, t: Date.now() }));
    openAnalysis();
  } catch (e) {
    throw asUserError_(e, 'openAnalysisWith');
  }
}

/** 분석 창이 열릴 때 호출: 대기 중인 추천이 있으면 한 번만 돌려주고 지운다 (2분이 지난 것은 버린다). 없으면 null */
function takePendingAdvice() {
  try {
    var props = PropertiesService.getUserProperties();
    var raw = props.getProperty(ADVISOR_PENDING_KEY_);
    if (!raw) return null;
    props.deleteProperty(ADVISOR_PENDING_KEY_);
    var p = JSON.parse(raw);
    if (!p || typeof p.key !== 'string' || Date.now() - Number(p.t) > ADVISOR_PENDING_TTL_MS_) return null;
    return { key: p.key, roles: p.roles || {} };
  } catch (e) {
    return null;
  }
}

/**
 * 체크한 변수의 척도로 바로 추천 (순수 함수). items: [{ name, scale('scale'|'ordinal'|'nominal'), groups(서로 다른 값의 수) }]
 * 반환: { key, roles: { 역할키: [변수] }, why, alt: [{ key, why, roles }], notes } 또는 { none: true, message }
 */
function adviseFromItems_(items) {
  var none = function (m) { return { none: true, message: m }; };
  if (!items || !items.length) return none('변수를 체크하면 그 변수로 바로 추천해 드립니다. 아니면 아래 질문으로 고르세요.');
  var nm = function (list) { return list.map(function (x) { return x.name; }); };
  var isNum = function (x) { return x.scale === 'scale' || x.scale === 'ordinal'; };
  var small = function (x) { return x.groups >= 2 && x.groups <= ADVISOR_MAX_GROUP_VALUES_; };
  if (items.length === 1) {
    var one = items[0];
    return one.scale === 'nominal'
      ? { key: 'frequencies', roles: { vars: nm(items) }, why: '범주형 변수의 값별 사례 수와 비율을 센다.', alt: [], notes: [] }
      : { key: 'descriptives', roles: { vars: nm(items) }, why: '숫자 변수의 평균, 표준편차 등으로 요약한다.', alt: [{ key: 'normality', why: '정규분포에 가까운지 확인한다.', roles: { vars: nm(items) } }], notes: [] };
  }
  if (items.length === 2) {
    var a = items[0], b = items[1];
    var num = isNum(a) ? a : (isNum(b) ? b : null), grp = a.scale === 'nominal' ? a : (b.scale === 'nominal' ? b : null);
    if (num && grp && num !== grp) {
      if (!small(grp)) return none('집단변수 "' + grp.name + '"의 값이 ' + grp.groups + '가지라 집단 비교로 보기 어렵습니다. 아래 질문으로 고르세요.');
      var roles = { dependent: [num.name], factor: [grp.name] };
      if (grp.groups === 2) {
        return num.scale === 'scale'
          ? { key: 'ttestInd', roles: roles, why: '"' + grp.name + '" 두 집단의 "' + num.name + '" 평균을 비교한다.', alt: [{ key: 'mannWhitney', why: '정규분포가 아니거나 표본이 작으면 순위로 비교한다.', roles: roles }], notes: [] }
          : { key: 'mannWhitney', roles: roles, why: '서열 점수 "' + num.name + '"을(를) "' + grp.name + '" 두 집단에서 순위로 비교한다.', alt: [{ key: 'ttestInd', why: '여러 문항의 합계·평균 점수라면 평균 비교도 가능하다.', roles: roles }], notes: [] };
      }
      return num.scale === 'scale'
        ? { key: 'anova', roles: roles, why: '"' + grp.name + '" ' + grp.groups + '개 집단의 "' + num.name + '" 평균을 비교한다.', alt: [{ key: 'kruskalWallis', why: '정규분포가 아니거나 표본이 작으면 순위로 비교한다.', roles: roles }], notes: [] }
        : { key: 'kruskalWallis', roles: roles, why: '서열 점수 "' + num.name + '"을(를) "' + grp.name + '" ' + grp.groups + '개 집단에서 순위로 비교한다.', alt: [{ key: 'anova', why: '여러 문항의 합계·평균 점수라면 평균 비교도 가능하다.', roles: roles }], notes: [] };
    }
    if (a.scale === 'nominal' && b.scale === 'nominal') {
      var alts = [];
      if (a.groups === 2 && b.groups === 2) alts.push({ key: 'mcnemar', why: '같은 사람의 사전·사후 예/아니오 측정이라면 변화가 있었는지 McNemar 검정으로 본다.', roles: { first: [a.name], second: [b.name] } });
      return { key: 'crosstab', roles: { rowvar: [a.name], colvar: [b.name] }, why: '두 범주형 변수가 관련 있는지 교차표와 카이제곱으로 본다.', alt: alts, notes: [] };
    }
    if (isNum(a) && isNum(b)) {
      var any = a.scale === 'ordinal' || b.scale === 'ordinal';
      return {
        key: 'correlation', roles: { vars: nm(items) }, why: '두 숫자 변수의 연관 정도(상관계수)를 본다.' + (any ? ' 서열 변수가 있어 옵션에서 스피어만 순위상관을 함께 켜세요.' : ''),
        alt: [{ key: 'ttestPaired', why: '같은 사람의 사전·사후 측정이라면 평균이 달라졌는지 비교한다.', roles: { first: [a.name], second: [b.name] } }], notes: []
      };
    }
    return none('이 변수 조합은 바로 추천하기 어렵습니다. 아래 질문으로 고르세요.');
  }
  if (items.every(isNum)) {
    return {
      key: 'correlation', roles: { vars: nm(items) }, why: '숫자 변수들 사이의 상관행렬을 본다.',
      alt: [
        { key: 'reliability', why: '같은 것을 재는 설문 문항이라면 문항 사이 일관성(크론바흐 α)을 확인한다.', roles: { vars: nm(items) } },
        { key: 'rmanova', why: '같은 사람을 3회 이상 측정한 값(사전·중간·사후)이라면 평균이 달라졌는지 반복측정 분산분석으로 본다. 변수는 측정 순서대로 보내진다.', roles: { vars: nm(items) } },
        { key: 'friedman', why: '같은 반복 측정을 순위로 비교하는 비모수 검정이다 (정규분포가 아닐 때).', roles: { vars: nm(items) } }
      ], notes: []
    };
  }
  return none('변수가 여러 개이고 종류가 섞여 있어 바로 추천하기 어렵습니다. 변수를 두 개 이하로 체크하거나 아래 질문으로 고르세요.');
}

/** 대화상자용: 체크한 변수 이름으로 추천. 데이터에서 척도와 집단 수를 읽는다. */
function adviseFromVariables(names) {
  try {
    var ds = loadDataset_();
    var items = (names || []).map(function (name) {
      var v = null;
      for (var i = 0; i < ds.vars.length; i++) if (ds.vars[i].name === name) { v = ds.vars[i]; break; }
      if (!v) throw new Error('변수 "' + name + '"을(를) 찾을 수 없습니다. 분석 창의 새로고침을 누르세요.');
      var seen = {}, count = 0;
      ds.data[name].forEach(function (x) { if (x !== null && !seen[typeof x + ':' + x]) { seen[typeof x + ':' + x] = true; count++; } });
      return { name: name, scale: v.scale, groups: count };
    });
    var res = adviseFromItems_(items);
    res.n = ds.n;
    return res;
  } catch (e) {
    throw asUserError_(e, 'adviseFromVariables');
  }
}

if (typeof module !== 'undefined') {
  module.exports = { ADVISOR_TREE_: ADVISOR_TREE_, ADVISOR_LABELS_: ADVISOR_LABELS_, getAdvisor: getAdvisor, adviseFromItems_: adviseFromItems_, ADVISOR_SMALL_N_: ADVISOR_SMALL_N_ };
}
