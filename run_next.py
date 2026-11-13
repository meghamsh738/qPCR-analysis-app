"""Start this qPCR worktree only, with separate local state and a strict port."""
import os
import socket
import subprocess
import sys
from pathlib import Path
from qpcr_storage import storage_paths, runtime_root


def main():
    root = runtime_root()
    os.environ['EASYLAB_RUNTIME_ROOT'] = str(root)
    paths = storage_paths()
    for key, variable in [('data','EASYLAB_DATA_PATH'),('settings','EASYLAB_SETTINGS_PATH'),('exports','EASYLAB_EXPORTS_PATH')]:
        os.environ[variable] = str(paths[key])
    port = 18501
    with socket.socket() as probe:
        try:
            probe.bind(('127.0.0.1', port))
        except OSError:
            raise SystemExit(f'Port {port} is already in use. No existing server was reused or stopped.')
    for variable in ['PYTHONHOME','PYTHONPATH','PYTHONUSERBASE','PYTHONPYCACHEPREFIX','PYTHON_EGG_CACHE','STREAMLIT_SERVER_PORT','STREAMLIT_SERVER_ADDRESS']:
        os.environ.pop(variable,None)
    os.environ['PYTHONNOUSERSITE'] = '1'
    os.environ['MPLCONFIGDIR'] = str(root / 'cache' / 'matplotlib')
    os.environ['STREAMLIT_BROWSER_GATHER_USAGE_STATS'] = 'false'
    command = [sys.executable, '-m', 'streamlit', 'run', str(Path(__file__).with_name('app.py')), '--server.address', '127.0.0.1', '--server.port', str(port), '--server.headless', 'true', '--browser.gatherUsageStats', 'false']
    os.chdir(str(Path(__file__).resolve().parent))
    os.execve(sys.executable, command, os.environ.copy())

if __name__ == '__main__':
    main()
