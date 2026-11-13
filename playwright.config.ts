import { defineConfig } from '@playwright/test'
import path from 'node:path'

const parentWorkspace = path.resolve(__dirname, '../..')
const testRuntime = path.basename(parentWorkspace) === 'easylab-next'
  ? path.join(parentWorkspace, 'runtime', 'browser-tests')
  : path.join(__dirname, '.next-runtime', 'easylab-next', 'browser-tests')
const python = process.platform === 'win32' ? '.venv/Scripts/python.exe' : '.venv/bin/python'
export default defineConfig({
  testDir: './tests',
  timeout: 120000,
  expect: { timeout: 20000 },
  use: { baseURL: 'http://127.0.0.1:18501', viewport: { width: 1400, height: 900 }, colorScheme: 'light', reducedMotion: 'reduce' },
  reporter: [['list']],
  webServer: {
    command: `${python} run_next.py`,
    env: { EASYLAB_RUNTIME_ROOT: testRuntime, STREAMLIT_BROWSER_GATHER_USAGE_STATS:'false' },
    url: 'http://127.0.0.1:18501/_stcore/health',
    reuseExistingServer: false,
    timeout: 120000,
  },
})
