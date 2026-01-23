import { test, expect } from '@playwright/test'

async function snapSection(page, headingRegex: RegExp, filename: string, scrollOffset = 0) {
  const heading = page.getByRole('heading', { name: headingRegex })
  await expect(heading).toBeVisible()
  await heading.scrollIntoViewIfNeeded()
  if (scrollOffset) {
    await page.evaluate((offset) => window.scrollBy(0, offset), scrollOffset)
  }
  await page.waitForTimeout(600)
  await expect(page).toHaveScreenshot(filename)
}

async function scrollToLargeImage(page) {
  const images = page.locator('img')
  const count = await images.count()
  for (let i = 0; i < count; i += 1) {
    const img = images.nth(i)
    const box = await img.boundingBox()
    if (box && box.width > 320 && box.height > 200) {
      await img.scrollIntoViewIfNeeded()
      return true
    }
  }
  return false
}

test('example workflow screenshots', async ({ page }) => {
  await page.goto('/')
  await page.addStyleTag({ content: '* { transition: none !important; animation: none !important; }' })
  const mainHeading = page.getByTestId('stMainBlockContainer').getByRole('heading', { name: 'qPCR Analysis' })
  await expect(mainHeading).toBeVisible()
  await expect(page.getByText('Wells loaded')).toBeVisible()

  await page.waitForTimeout(1000)

  await expect(page).toHaveScreenshot('example_run.png', { fullPage: true })

  await snapSection(page, /1\)\s*Review & clean wells/i, 'overview.png')
  await snapSection(page, /2\)\s*Replicate averages/i, 'replicates.png')
  await snapSection(page, /3\)\s*Standards map/i, 'standards.png')
  const showToggle = page.getByLabel('Show curve plots')
  if (await showToggle.isVisible()) {
    if (!(await showToggle.isChecked())) {
      await showToggle.check()
      await page.waitForTimeout(400)
    }
  }
  await page.waitForTimeout(800)
  const foundPlot = await scrollToLargeImage(page)
  if (!foundPlot) {
    await snapSection(page, /4\)\s*Fit standard curves/i, 'curves.png', 720)
  } else {
    await page.waitForTimeout(400)
    await expect(page).toHaveScreenshot('curves.png')
  }
  await snapSection(page, /6\)\s*Normalize to reference gene/i, 'quant_normalize.png')
  await snapSection(page, /7\)\s*Export/i, 'export.png')
})
