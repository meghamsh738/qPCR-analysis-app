// This branch's desktop shell is isolated from the current standalone installation.
const { app, BrowserWindow, dialog } = require('electron')
const { spawn, execFileSync } = require('node:child_process')
const fs = require('node:fs')
const path = require('node:path')
const net = require('node:net')
const root = path.resolve(__dirname, '../..')
const developmentRoot = path.resolve(root, '../..')
const runtime = !app.isPackaged && path.basename(developmentRoot) === 'easylab-next'
  ? path.join(developmentRoot, 'runtime')
  : path.join(app.getPath('documents'), 'EasyLab Next', 'runtime')
const profile = path.join(runtime, 'profiles', 'qpcr-analysis')
// Fixed isolated paths: inherited stable profile/data overrides are deliberately ignored.
for (const p of [runtime, profile]) {
  let existing = p
  while (!fs.existsSync(existing)) existing = path.dirname(existing)
  if (fs.realpathSync(existing) !== existing) throw new Error('Next storage cannot use symbolic links')
}
fs.mkdirSync(profile, { recursive: true })
app.setName('EasyLab Next qPCR Analysis')
app.setPath('userData', profile)
app.setPath('sessionData', path.join(profile, 'session'))
app.setPath('logs', path.join(profile, 'logs'))
const port = 18501
let child
let window
const free = () => new Promise((resolve, reject) => {
  const probe = net.createServer()
  probe.once('error', reject)
  probe.listen(port, '127.0.0.1', () => probe.close(resolve))
})
const wait = () => new Promise((resolve, reject) => {
  const start = Date.now()
  const check = () => {
    if (!child || child.exitCode !== null) return reject(new Error('Next Python server stopped'))
    const socket = net.createConnection({ host: '127.0.0.1', port })
    socket.once('connect', () => {
      socket.destroy()
      try {
        const output = process.platform === 'win32'
          ? execFileSync('powershell.exe', ['-NoProfile','-Command',`(Get-NetTCPConnection -LocalAddress 127.0.0.1 -LocalPort ${port} -State Listen).OwningProcess`], {encoding:'utf8',timeout:5000})
          : execFileSync('/usr/sbin/lsof', ['-nP',`-tiTCP:${port}`,'-sTCP:LISTEN'], {encoding:'utf8',timeout:5000})
        const owners = output.trim().split(/\s+/).filter(Boolean).map(Number)
        if (owners.length !== 1 || owners[0] !== child.pid) throw new Error('The port belongs to another process; Next refused to open it')
        resolve()
      } catch(error) { reject(error) }
    })
    socket.once('error', () => { socket.destroy(); Date.now()-start > 30000 ? reject(new Error('Next server timed out')) : setTimeout(check, 200) })
  }
  check()
})
app.whenReady().then(async () => {
  try {
    await free()
    const source = app.isPackaged ? process.resourcesPath : root
    const python = path.join(source, '.venv', process.platform === 'win32' ? 'Scripts/python.exe' : 'bin/python')
    if (!fs.existsSync(python)) throw new Error('Create this Next checkout’s Python environment first, or launch qPCR through EasyLab Next Suite. No existing installation was used.')
    const env = { ...process.env, EASYLAB_RUNTIME_ROOT: runtime, EASYLAB_DATA_PATH:path.join(runtime,'data','qPCR Analysis','data'), EASYLAB_SETTINGS_PATH:path.join(runtime,'settings','qpcr-analysis'), EASYLAB_EXPORTS_PATH:path.join(runtime,'data','qPCR Analysis','exports'), MPLCONFIGDIR:path.join(runtime,'cache','matplotlib'), PYTHONPATH:source, PYTHONNOUSERSITE:'1' }
    child = spawn(python, [path.join(source, 'run_next.py')], { cwd:source, env, stdio:'ignore', windowsHide:true })
    child.on('error', () => {})
    await wait()
    window = new BrowserWindow({width:1400,height:900,minWidth:600,minHeight:500,title:'EasyLab Next · qPCR',webPreferences:{contextIsolation:true,nodeIntegration:false,sandbox:true}})
    window.webContents.setWindowOpenHandler(() => ({action:'deny'}))
    window.webContents.on('will-navigate', (event, url) => { if (!url.startsWith(`http://127.0.0.1:${port}/`)) event.preventDefault() })
    await window.loadURL(`http://127.0.0.1:${port}`)
  } catch (error) {
    dialog.showErrorBox('EasyLab Next could not start', error.message)
    app.quit()
  }
})
app.on('before-quit', () => { if (child && child.exitCode === null) child.kill() })
app.on('window-all-closed', () => app.quit())
