const { contextBridge, ipcRenderer } = require('electron')

contextBridge.exposeInMainWorld('electronAPI', {
  selectDirectory: (options) => ipcRenderer.invoke('select-directory', options),
  ensureDirectories: (paths) => ipcRenderer.invoke('ensure-directories', paths),
  getAppInfo: () => ipcRenderer.invoke('get-app-info'),
  getDefaultPaths: () => ipcRenderer.invoke('get-default-paths'),
})
