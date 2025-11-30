@echo off
setlocal
set "WSL_PATH='/mnt/d/coding projects/qPCR-analysis-app'"
start "qPCR analysis" wsl -e bash -lc "cd %WSL_PATH% && source .venv/bin/activate && streamlit run app.py --server.port 8501 --server.headless true"
timeout /t 6 >nul
start "" http://localhost:8501
endlocal
