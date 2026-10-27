@echo off
setlocal enabledelayedexpansion
chcp 65001 >nul
cd /d "%~dp0"
title Stat4SW 로컬 커밋

rem ============================================================
rem  Stat4SW 로컬 커밋 (이 컴퓨터의 git 에만 저장합니다)
rem  - 이 파일은 push / pull / 원격 저장소 연결을 하지 않습니다 (GitHub 등으로 보내지 않음).
rem  - 사용법
rem      commit.bat                  : 테스트 실행 -> 변경 파일 저장(커밋)
rem      commit.bat "메시지"         : 메시지를 직접 지정 (따옴표 안에 느낌표 ! 는 쓰지 마세요)
rem      commit.bat /notest          : 테스트를 건너뜀
rem      commit.bat "메시지" /notest : 둘 다
rem  - 긴 메시지는 같은 폴더에 commit_message.txt (UTF-8로 저장)를 만들어 두면 그 내용을 쓰고, 커밋 뒤 지웁니다.
rem  - Stat_Assistant 폴더는 참고 자료라 기본으로 제외합니다. 포함하려면 아래 EXCLUDE 값을 지우세요.
rem ============================================================
set "EXCLUDE=Stat_Assistant"

echo.
echo === Stat4SW 로컬 커밋 (원격 저장소로 보내지 않습니다) ===
echo.

where git >nul 2>nul
if errorlevel 1 (
  echo [오류] git 을 찾을 수 없습니다. git 을 설치하거나 PATH 에 추가하세요.
  goto :end
)
git rev-parse --is-inside-work-tree >nul 2>nul
if errorlevel 1 (
  echo [오류] 이 폴더는 git 저장소가 아닙니다.
  goto :end
)

set "MSG="
set "NOTEST="
for %%A in (%*) do (
  if /i "%%~A"=="/notest" (set "NOTEST=1") else if not defined MSG set "MSG=%%~A"
)

rem --- 원격 저장소 알림 (이 파일은 보내지 않음) ---
set "HASREMOTE="
for /f "delims=" %%R in ('git remote') do set "HASREMOTE=1"
if defined HASREMOTE (
  echo [알림] 원격 저장소가 설정돼 있습니다. 이 파일은 push 를 하지 않으니 안심하세요.
  echo.
)

rem --- 테스트 ---
if defined NOTEST (
  echo [건너뜀] 테스트를 실행하지 않습니다.
) else (
  where node >nul 2>nul
  if errorlevel 1 (
    echo [건너뜀] node 가 없어 테스트를 실행하지 않습니다.
  ) else (
    echo 테스트 실행 중... 잠시 기다려 주세요.
    node --test "tests/*.test.js" >"%TEMP%\stat4sw_test.log" 2>&1
    if errorlevel 1 (
      echo [실패] 테스트가 실패했습니다. 커밋하지 않습니다. 자세한 내용: %TEMP%\stat4sw_test.log
      echo         테스트를 건너뛰고 커밋하려면 commit.bat /notest 로 실행하세요.
      goto :end
    )
    echo [통과] 테스트 통과
  )
)
echo.

rem --- 변경 파일 담기 (push 하지 않음) ---
if defined EXCLUDE (
  git add -A -- . ":(exclude)%EXCLUDE%"
) else (
  git add -A
)

git diff --cached --quiet
if not errorlevel 1 (
  echo 커밋할 변경이 없습니다.
  goto :show
)

echo 커밋할 파일:
git diff --cached --stat
echo.

rem --- 커밋 메시지: commit_message.txt > 인자 > 입력 > 날짜 순 ---
if exist commit_message.txt (
  git commit -F commit_message.txt
  if errorlevel 1 (
    echo [오류] 커밋에 실패했습니다.
  ) else (
    del commit_message.txt >nul 2>nul
    echo [완료] 로컬 커밋을 만들었습니다.
  )
  goto :show
)

if not defined MSG (
  set /p "MSG=커밋 메시지 (그냥 Enter 를 누르면 날짜로 자동 작성): "
)
if not defined MSG set "MSG=작업 저장 %date% %time:~0,8%"
set "MSG=!MSG:"=!"
git commit -m "!MSG!"
if errorlevel 1 (
  echo [오류] 커밋에 실패했습니다.
) else (
  echo [완료] 로컬 커밋을 만들었습니다.
)

:show
echo.
echo --- 최근 커밋 ---
git log --oneline -5
echo.
echo --- 남은 변경 (커밋 안 된 파일) ---
git status --short
echo.
echo (이 컴퓨터에만 저장되었고 온라인으로 보내지 않았습니다)

:end
if not defined NOPAUSE pause
endlocal
