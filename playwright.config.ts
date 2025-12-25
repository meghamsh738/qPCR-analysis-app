import { defineConfig } from '@playwright/test'

export default defineConfig({
  testDir: './tests',
  timeout: 120000,
  expect: {
    timeout: 15000,
  },
  use: {
    baseURL: 'http://127.0.0.1:8501',
    viewport: { width: 1440, height: 900 },
    colorScheme: 'light',
  },
  reporter: [['list']],
  webServer: {
    command: 'python3 -m streamlit run app.py --server.headless true --server.port 8501 --server.address 127.0.0.1',
    url: 'http://127.0.0.1:8501',
    reuseExistingServer: !process.env.CI,
    timeout: 120000,
  },
})
