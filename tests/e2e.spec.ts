import { test, expect } from '@playwright/test'
import fs from 'fs/promises'
import path from 'path'
import { fileURLToPath } from 'url'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)
const screenshotDir = path.join(__dirname, '..', 'screenshots')

async function snapSection(page, headingRegex: RegExp, filename: string) {
  const heading = page.getByRole('heading', { name: headingRegex })
  await expect(heading).toBeVisible()
  await heading.scrollIntoViewIfNeeded()
  await page.waitForTimeout(600)
  await page.screenshot({ path: path.join(screenshotDir, filename) })
}

test('example workflow screenshots', async ({ page }) => {
  await fs.mkdir(screenshotDir, { recursive: true })

  await page.goto('/')
  const mainHeading = page.getByTestId('stMainBlockContainer').getByRole('heading', { name: 'qPCR Analysis' })
  await expect(mainHeading).toBeVisible()
  await expect(page.getByText('Wells loaded')).toBeVisible()

  await page.waitForTimeout(1000)

  await page.screenshot({ path: path.join(screenshotDir, 'example_run.png'), fullPage: true })

  await snapSection(page, /1\)\s*Review & clean wells/i, 'overview.png')
  await snapSection(page, /2\)\s*Replicate averages/i, 'replicates.png')
  await snapSection(page, /3\)\s*Standards map/i, 'standards.png')
  await snapSection(page, /4\)\s*Fit standard curves/i, 'curves.png')
  await snapSection(page, /6\)\s*Normalize to reference gene/i, 'quant_normalize.png')
  await snapSection(page, /7\)\s*Export/i, 'export.png')
})
