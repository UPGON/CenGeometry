@echo off
REM ===========================================================
REM   THIS FILE IS FOR WINDOWS ONLY.
REM   On a Mac, double-click "START HERE (Mac).command" instead
REM   -- a .bat will just open in a text editor and do nothing.
REM ===========================================================
REM Double-click this file to open CenGeometry in your web browser.
REM First run takes a couple of minutes while it installs; after that it is quick.

cd /d "%~dp0"
cls
echo ===============================================
echo   CenGeometry
echo ===============================================
echo.

where python >nul 2>&1
if errorlevel 1 (
    echo Python was not found on this computer.
    echo.
    echo Install it from https://www.python.org/downloads/
    echo Be sure to tick "Add Python to PATH" during installation,
    echo then double-click this file again.
    echo.
    pause
    exit /b 1
)

if not exist ".venv" (
    echo First-time setup - this happens once and takes a few minutes.
    echo.
    python -m venv .venv
    if errorlevel 1 ( echo Could not create the workspace. & pause & exit /b 1 )
    call .venv\Scripts\python.exe -m pip install --quiet --upgrade pip
    call .venv\Scripts\python.exe -m pip install --quiet -r requirements.txt
    if errorlevel 1 ( echo Could not install the required packages. & pause & exit /b 1 )
    echo Setup complete.
    echo.
)

echo Starting... your browser will open automatically.
echo Leave this window open while you work.
echo To stop, close this window or press Control-C.
echo.
call .venv\Scripts\streamlit.exe run app.py
