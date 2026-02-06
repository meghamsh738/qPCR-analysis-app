import { readFile, writeFile, access } from 'node:fs/promises'
import path from 'node:path'
import sharp from 'sharp'
import pngToIco from 'png-to-ico'

const cwd = path.resolve(process.cwd())
const directPath = path.join(cwd, 'build', 'icon.svg')
const desktopPath = path.join(cwd, 'desktop', 'build', 'icon.svg')

let baseDir = path.join(cwd, 'desktop')
try {
  await access(directPath)
  baseDir = cwd
} catch {
  await access(desktopPath)
}

const svgPath = path.join(baseDir, 'build', 'icon.svg')
const pngPath = path.join(baseDir, 'build', 'icon.png')
const icoPath = path.join(baseDir, 'build', 'icon.ico')

const svgBuffer = await readFile(svgPath)
const pngBuffer = await sharp(svgBuffer).resize(1024, 1024).png().toBuffer()
await writeFile(pngPath, pngBuffer)
const icoBuffer = await pngToIco(pngBuffer)
await writeFile(icoPath, icoBuffer)

console.log(`Generated ${pngPath} and ${icoPath}`)
