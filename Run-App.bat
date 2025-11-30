@echo off
setlocal
for /f "usebackq tokens=* delims=" %%i in (`wsl wslpath "%~dp0"`) do set WSL_DIR=%%i
wsl -e bash -lc "netstat -tlnp | grep -E ':(8501) ' | awk '{print $7}' | cut -d/ -f1 | xargs -r kill -9"
wsl -e bash -lc "cd '%WSL_DIR%' && source .venv/bin/activate && streamlit run app.py --server.port 8501 --server.headless true"
timeout /t 6 >nul
start "" http://localhost:8501
endlocal
