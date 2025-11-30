@echo off
setlocal
for /f "usebackq tokens=*" %%i in (`wsl wslpath "%~dp0"`) do set WSL_DIR=%%i
start "qPCR analysis" wsl -e bash -lc "cd \"%WSL_DIR%\" && source .venv/bin/activate && streamlit run app.py --server.port 8501 --server.headless true"
timeout /t 6 >nul
start "" http://localhost:8501
endlocal
