@echo %off
setlocal enabledelayedexpansion

set "RootDir=%~dp0"

pushd "%RootDir%"
set "FullPath=%CD%"
popd

set "RPathFound=false"

:: First, try to find R in the portable directory
if exist "%RootDir%R\bin\R.exe" (
    set "R=%RootDir%R\bin\R.exe"
    set "RPathFound=true"
)

@REM :: If not found, search for R in default installation paths
@REM if not !RPathFound! == true (
@REM     for %%d in ("C:\Program Files\R" "C:\Program Files (x86)\R") do (
@REM         for /d %%i in (%%d\R-*) do (
@REM             if exist "%%i\bin\x64\R.exe" (
@REM                 set "R=%%i\bin\x64\R.exe"
@REM                 set "RPathFound=true"
@REM                 goto :findDir
@REM             )
@REM         )
@REM     )
@REM )

:: If R is still not found, prompt the user
:promptRPath
if not !RPathFound! == true (
    set /p "R=Cannot detect R installation location. Please specify the full path to your R installation (R.exe): "
    if not exist "!R!" (
        echo The specified path does not exist. Please try again.
        goto :promptRPath
    ) else (
        set "RPathFound=true"
        echo R found at !R!
    )
)

:: Define the app directory
set "AppDir=%RootDir%app"


:runApp
REM Set debug mode (true or false)
set "DEBUG=true"

"%R%" --no-save --slave -f "%RootDir%run.R" --args "%FullPath%\app" %DEBUG%

pause
exit /b 0

