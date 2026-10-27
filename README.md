# Stat4SW — Google Sheets + GAS 통계 패키지

구글 시트 안에서 통계 분석을 하는 오픈소스(MIT) 도구입니다. 데이터는 시트의 탭에 두고, 메뉴 **Stat4SW 통계분석 > 분석 창 열기**에서 변수를 고르면 결과가 표(와 차트)로 쌓입니다.
시트에 바인딩된 Apps Script(GAS)로 만들었고, 외부 라이브러리 없이 직접 구현했습니다.

## 사용하기
### 방법 1. 배포된 템플릿 쓰기 (가장 쉬움)
[릴리스 페이지](../../releases/latest)의 설명에 있는 **템플릿 링크**를 열고 "사본 만들기"를 누르면 내 구글 계정에 내 시트가 생깁니다. 스크립트도 함께 복사되고 데이터는 내 시트에만 있습니다. 처음 메뉴를 실행할 때 권한 승인이 필요합니다. 자세한 사용법은 [doc/사용안내.md](doc/사용안내.md)를 보세요.

### 방법 2. 소스로 직접 만들기
1. 새 구글 시트를 만들고 **확장 프로그램 > Apps Script**를 엽니다.
2. `src/`의 `.js` 파일은 같은 이름의 스크립트 파일로, `.html` 파일은 같은 이름의 HTML 파일로 추가하고, `appsscript.json` 내용을 프로젝트 설정의 매니페스트에 반영합니다. (또는 아래 "개발 환경 설정"의 clasp로 올립니다.)
3. 시트를 새로고침하면 메뉴 **Stat4SW 통계분석**이 생깁니다. `새 작업 시작 (데이터 모두 지우기)`으로 `데이터`·`변수정의` 탭을 만듭니다.

## 라이선스
MIT 라이선스 ([LICENSE](LICENSE)). 저작권자 표기와 계산 방법 출처는 [doc/저작권·라이선스.md](doc/저작권·라이선스.md)와 앱의 `정보·저작권` 화면에 있습니다.

## 문서
| 문서 | 내용 |
| --- | --- |
| [doc/사용안내.md](doc/사용안내.md) | 사용자용: 시작하기, 데이터 준비, 분석 창 사용법, 문제 해결 |
| [doc/인수인계.md](doc/인수인계.md) | 현재 상태, 남은 일, 확인할 사항, 작업 방법 (다음 세션 시작점) |
| [doc/릴리스절차.md](doc/릴리스절차.md) | 배포용 템플릿 만들기, 공유, 업데이트 안내 |
| [CHANGELOG.md](CHANGELOG.md) | 버전별 변경 이력 |
| [doc/개발계획.md](doc/개발계획.md) | 단계별 계획, 설계 결정, 진행 상황, 나중에 검토할 일 |
| [doc/설계-결과객체.md](doc/설계-결과객체.md) | 분석 결과 객체와 렌더러(표·차트) 형식 |
| [doc/확인사항.md](doc/확인사항.md) | **점검 절차서**: 기능별 수행 방법·기대 결과·확인 상태(요약표), 릴리스 전 점검 묶음, 알려진 한계, 작업 환경 메모 |

## 기능
- 기초: 빈도분석, 기술통계, 정규성(Shapiro-Wilk, Jarque-Bera), 등분산(Levene), 이항 검정(정확 유의확률, Clopper–Pearson·Wilson 신뢰구간)
- 평균 비교: 일표본·독립표본·대응표본 t검정, 일원 분산분석(Welch, 효과크기), 사후검정(Tukey, Games-Howell, Bonferroni, Scheffé, LSD, Šidák, Dunnett), 이원 분산분석(유형 III·II, 상호작용 그림, 단순 주효과)
- 비모수 검정: Mann-Whitney U, Wilcoxon 부호순위(대응표본·일표본), Kruskal-Wallis(+ Dunn 사후검정)
- 반복측정 분산분석(Mauchly 구형성, Greenhouse–Geisser·Huynh–Feldt 교정, 쌍별 비교; 집단변수를 넣으면 집단 × 측정 혼합 설계: 단순 효과, 상호작용 그림)
- Cochran Q 검정(예/아니오 3회 이상 반복측정, 정확 McNemar 쌍별 비교)
- Friedman 검정(3회 이상 반복측정, Dunn–Bonferroni 사후비교, Kendall W)
- t검정의 효과크기(Cohen's d, Hedges' g, 비중심 t 정확 신뢰구간)
- McNemar 검정(정확 유의확률, 오즈비 정확 구간), 분석 선택 도우미(변수로 추천, 질문, 넓게 보기)
- 의사결정나무(CART): 분류·회귀, 비용복잡도 가지치기와 교차검증, 나무 그림, 변수 중요도
- k-평균 군집분석: k-means++, 군집 수 판단 보조(엘보·실루엣), 군집 번호 저장
- 성향점수 매칭: 로지스틱 성향점수, 최근접·마할라노비스 매칭(캘리퍼, 1:k, 복원), 균형 진단(SMD), 효과(ATT)
- 거리·유사도(14가지 측도, 표준화, 근접도 행렬 저장)와 마할라노비스 이상치 점검
- 서열(비례오즈) 로지스틱 회귀: 문턱값·B·Exp(B), 비례오즈 가정 점검(Brant), 분류표, 예측 저장
- 로지스틱 회귀: 오즈비, 모형 검정, Hosmer-Lemeshow, 분류표, ROC·AUC, 예측값 저장 (사건 값·기준 범주를 분석 창에서 바로 고름)
- 매개효과: 병렬 다중 매개, 간접효과의 부트스트랩 신뢰구간과 Sobel, 경로도
- 조절효과: 상호작용 검정(ΔR²), 단순기울기, Johnson-Neyman, 상호작용 그림 (연속·범주형 조절변수)
- 신뢰도: 크론바흐 α (문항-전체 통계량, 문항 제외 시 α, Feldt 신뢰구간)
- Kaplan–Meier 생존분석(중도절단, 생존율·신뢰구간, 중앙값, log-rank·Breslow·Tarone–Ware, 쌍별 비교)
- 검정력·표본 크기 계산(t검정, 분산분석, 상관, 두 비율, 회귀; 데이터 불필요, 탈락 감안 모집 인원)
- IPA 중요도–수행도 분석(4사분면 그림, 기준선 방식·이동·직접 지정, 격차 검정)
- 평가자 일치도: Cohen·Fleiss kappa(가중 kappa, 쌍별 kappa), 급내상관 ICC 6유형(신뢰구간)
- 탐색적 요인분석(주축분해·주성분, 베리맥스·프로맥스, KMO·Bartlett, 스크리 그림, 병렬분석, 요인 점수 저장)
- 관계: 상관분석(피어슨, 스피어만, 켄달 타우-b, 통제변수를 넣으면 편상관·순위 편상관), 교차분석(카이제곱, Fisher 정확검정, 효과크기)
- 회귀: 선형 회귀(범주형 변수는 더미로 자동 변환, 표준화 계수, VIF, 공선성 진단, 더빈-왓슨), 위계적 회귀(1~4단계, R² 변화량과 F 변화 검정)
- 시험 데이터: 기법별 파이썬 기준값 데이터를 시트에 넣고 기대 결과와 비교 (메뉴 Stat4SW 통계분석 > 설정·관리 > 시험 데이터)
- 예제 데이터: 변수 27개(60명), 분석마다 연습할 변수가 있음 (`Sample.js`)
- 변수 계산: 역문항 처리, 문항 합계·평균, z점수 (데이터 탭의 새 열과 변수정의에 등록)
- 차트(옵션): 막대그래프, 히스토그램, 산점도, 상자그림
- 결과를 `결과` 탭에 표 서식 그대로 내보내기 (차트는 그림으로 삽입), 분석 옵션 기억

## 기여와 커밋
- 코드를 바꾸면 `node --test tests/*.test.js`가 모두 통과해야 하고, 새 분석은 독립된 기준값 테스트와 시험 데이터, 점검 절차서 항목이 함께 있어야 합니다 ([doc/인수인계.md](doc/인수인계.md)의 "작업 방법").
- `commit.bat`은 테스트를 돌린 뒤 변경 파일을 커밋하는 보조 도구입니다 (올리는 `git push`는 하지 않습니다).

## 배포 방식과 주의
- 이 시트를 **사본 만들기**로 복사해서 쓰는 템플릿 방식입니다. 스크립트도 함께 복사되며, 복사한 시트에서 처음 메뉴를 실행할 때 권한 승인이 필요합니다("확인되지 않은 앱" 경고 → 고급 → 이동). 릴리스 방법은 [doc/릴리스절차.md](doc/릴리스절차.md)를 보세요.
- **익명화된 데이터를 전제로 합니다.** 이름·연락처·주민번호 등 식별정보는 시트에 넣지 마세요. 코드는 별도의 개인정보 보호 처리를 하지 않습니다.
- 이미 복사한 시트에는 이후 코드 수정이 자동 반영되지 않습니다. 메뉴 `Stat4SW 통계분석 > 정보·저작권`에서 버전을 확인하고, 새 버전은 릴리스 페이지의 새 템플릿으로 옮겨 가세요 ([doc/사용안내.md](doc/사용안내.md)).

## 폴더 구조
```
Stat4SW/
├─ .clasp.json.example # clasp 설정 예시 (실제 .clasp.json 은 로컬 전용)
├─ CHANGELOG.md
├─ doc/               # 계획, 설계, 사용 안내, 릴리스 절차, 확인사항
├─ src/               # clasp가 push/pull 하는 Apps Script 소스
│  ├─ Menu.js, Template.js, Release.js, Config.js      메뉴, 템플릿·예제 데이터, 배포용 점검
│  ├─ Api.js, Safety.js, Export.js                       분석 창이 호출하는 서버 함수, 옵션 기억·오류 안내, 시트 내보내기
│  ├─ Data.js                                            변수정의·데이터셋·결측 처리
│  ├─ Distributions.js                                   정규, t, F, χ², Studentized range
│  ├─ Stats_*.js                                         분석별 계산 (Descriptive, Assumptions, Compare, Anova, Relation, Regression)
│  ├─ Render.js, Charts.js                               결과 객체 → HTML 표, SVG 차트
│  ├─ index.html, client.html, style.html, help.html     분석 창과 사용 안내
│  └─ TestDistributions.js                               시트 내장 함수와 비교하는 검증 함수 (편집기에서 실행)
└─ tests/             # Node 테스트, scipy 기준값 생성 스크립트와 기준값(reference/), UI 미리보기(ui/)
```

## 개발 환경 설정 (최초 1회)
1. Node.js 24 이상 확인: `node -v`
2. clasp 설치: `npm install -g @google/clasp`
3. Apps Script API 켜기: https://script.google.com/home/usersettings → "Google Apps Script API" ON
4. 로그인 (스크립트 소유 계정으로): `clasp login`
5. `.clasp.json.example`을 `.clasp.json`으로 복사하고 scriptId 를 채운 뒤, 이 폴더에서 `clasp push`로 올리거나 `clasp pull`로 받아 `src/`와 맞춥니다.
6. 기준값을 다시 만들 일이 있으면 Python 3.12와 scipy, numpy가 필요합니다 (Windows에서는 `py -3.12`).

## 평소 작업
| 작업 | 명령 |
|---|---|
| 웹 편집기 변경사항 가져오기 | `clasp pull` |
| 로컬 수정사항 올리기 | `clasp push` |
| 올리기 전 변경 파일 확인 | `clasp status` |
| 웹 편집기 열기 | `clasp open-script` (구버전: `clasp open`) |

- `clasp push`는 원격 파일을 **로컬 내용으로 통째로 덮어씁니다.** 웹 편집기에서 고친 게 있으면 먼저 `clasp pull`.
- 웹 편집기와 VS Code를 동시에 수정하지 말고, 작업 시작 시 pull → 끝나면 push 습관을 들이세요.
- 커밋 메시지는 파일로 만들어 `git commit -F`로 쓰는 것이 안전합니다 (PowerShell에서 큰따옴표가 깨지는 것을 피함).
- 템플릿 방식이라 `clasp create-version`이나 배포(deployments)는 쓰지 않습니다.

## 테스트
Node 24 이상, 외부 패키지 없음.
```
node --test tests/*.test.js
```
- `tests/reference/*.json`은 기준값입니다. scipy로 만들고, scipy가 부정확한 경우(평균이 큰 자료 등)에는 분수 정확 산술로 만들어 scipy와 교차 확인합니다. 다시 만들려면 `py -3.12 tests/gen_reference_<이름>.py` (`gen_reference.py`는 확률분포).
- `tests/gas_sim.test.js`는 Apps Script의 "모든 파일이 한 전역 범위"를 흉내 내 파일 사이의 전역 이름 충돌과 로드 순서 의존, 모든 분석의 실행·내보내기·옵션 저장을 확인합니다 (실제 Apps Script는 아닙니다).
- 시트 내장 함수와의 비교는 Apps Script 편집기에서 `compareDistributionsWithSheet`를 실행합니다.
- 배포용 템플릿 점검은 편집기에서 `checkTemplateReady`를 실행합니다.
- 분석 창 UI는 `node tests/ui/build_preview.js`(`--reverse`: 파일을 역순으로 불러 로드 순서 의존 점검)로 만든 `tests/ui/preview.html`을 브라우저에서 열어 로컬로 확인할 수 있습니다 (실제 서버 코드를 그대로 쓰고 시트만 흉내 냅니다).
