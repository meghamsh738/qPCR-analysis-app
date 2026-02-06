const { app, BrowserWindow, ipcMain, dialog } = require('electron')
const { spawn } = require('child_process')
const fs = require('fs')
const net = require('net')
const path = require('path')

const isDev = !app.isPackaged
const rootDir = path.join(__dirname, '..', '..')
const iconPath = path.join(__dirname, '..', 'build', 'icon.png')

const SERVER_PORT = 8501
let serverProcess = null

const ensureDirectories = (paths) => {
  const targets = Object.values(paths || {}).filter((val) => typeof val === 'string' && val.trim())
  targets.forEach((target) => fs.mkdirSync(target, { recursive: true }))
}

const getDefaultPaths = () => {
  const base = path.join(app.getPath('documents'), 'Easylab', 'qPCR Analysis')
  return {
    dataPath: path.join(base, 'data'),
    attachmentsPath: path.join(base, 'attachments'),
    exportsPath: path.join(base, 'exports'),
    syncPath: path.join(base, 'sync'),
  }
}

const waitForPort = (port, timeoutMs = 10000) => new Promise((resolve, reject) => {
  const start = Date.now()
  const check = () => {
    const socket = net.createConnection({ port }, () => {
      socket.end()
      resolve(true)
    })
    socket.on('error', () => {
      socket.destroy()
      if (Date.now() - start > timeoutMs) reject(new Error('timeout'))
      else setTimeout(check, 400)
    })
  }
  check()
})

const resolvePythonCandidates = () => {
  const candidates = []
  if (process.env.APP_PYTHON_PATH) candidates.push(process.env.APP_PYTHON_PATH)
  candidates.push('python', 'python3', 'py')
  return Array.from(new Set(candidates))
}

const spawnServer = async () => {
  if (serverProcess) return
  const appPath = isDev ? path.join(rootDir, 'app.py') : path.join(process.resourcesPath, 'app.py')
  const cwd = isDev ? rootDir : process.resourcesPath

  const env = {
    ...process.env,
    PYTHONPATH: cwd,
  }

  const candidates = resolvePythonCandidates()
  for (const candidate of candidates) {
    try {
      const proc = spawn(candidate, [
        '-m',
        'streamlit',
        'run',
        appPath,
        '--server.headless',
        'true',
        '--server.port',
        String(SERVER_PORT),
        '--server.address',
        '127.0.0.1',
      ], {
        cwd,
        env,
        stdio: 'ignore',
        windowsHide: true,
      })

      const ready = await Promise.race([
        waitForPort(SERVER_PORT, 8000),
        new Promise((_, reject) => proc.once('error', reject)),
      ])

      if (ready) {
        serverProcess = proc
        proc.on('exit', () => {
          serverProcess = null
        })
        return
      }
    } catch (err) {
      continue
    }
  }

  dialog.showMessageBox({
    type: 'warning',
    title: 'Server not started',
    message: 'The qPCR analysis server could not start. Install Python 3.10+ or set APP_PYTHON_PATH, then restart.',
  })
}

const stopServer = () => {
  if (serverProcess) {
    serverProcess.kill()
    serverProcess = null
  }
}

const createWindow = () => {
  const win = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 1024,
    minHeight: 640,
    resizable: true,
    backgroundColor: '#0F172A',
    title: app.getName(),
    icon: iconPath,
    webPreferences: {
      preload: path.join(__dirname, 'preload.cjs'),
      contextIsolation: true,
      nodeIntegration: false,
    },
  })

  const url = `http://127.0.0.1:${SERVER_PORT}`
  win.loadURL(url)
}

app.whenReady().then(async () => {
  await spawnServer()
  createWindow()

  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) createWindow()
  })
})

app.on('before-quit', () => {
  stopServer()
})

app.on('window-all-closed', () => {
  stopServer()
  if (process.platform !== 'darwin') app.quit()
})

ipcMain.handle('select-directory', async (_event, options = {}) => {
  const { title, defaultPath } = options
  const result = await dialog.showOpenDialog({
    title: title || 'Select folder',
    defaultPath,
    properties: ['openDirectory', 'createDirectory'],
  })
  if (result.canceled || result.filePaths.length === 0) return null
  return result.filePaths[0]
})

ipcMain.handle('ensure-directories', async (_event, paths) => {
  try {
    ensureDirectories(paths)
    return { ok: true }
  } catch (err) {
    return { ok: false, message: err instanceof Error ? err.message : String(err) }
  }
})

ipcMain.handle('get-app-info', () => ({
  name: app.getName(),
  version: app.getVersion(),
  platform: process.platform,
}))

ipcMain.handle('get-default-paths', () => getDefaultPaths())
