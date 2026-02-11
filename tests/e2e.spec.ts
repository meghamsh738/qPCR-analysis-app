import { test, expect } from '@playwright/test'

async function snapSection(page, headingRegex: RegExp, filename: string, scrollOffset = 0) {
  const banner = page.getByRole('banner')
  const heading = page.getByRole('heading', { name: headingRegex })
  await expect(heading).toBeVisible()
  await heading.scrollIntoViewIfNeeded()
  if (scrollOffset) {
    await page.evaluate((offset) => window.scrollBy(0, offset), scrollOffset)
  }
  await page.waitForTimeout(600)
  // Streamlit's top toolbar includes a "Running..." indicator that can animate and cause
  // tiny diffs between consecutive screenshots; mask it for stable visual snapshots.
  await expect(page).toHaveScreenshot(filename, { mask: [banner] })
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
  const banner = page.getByRole('banner')
  const setupHeading = page.getByRole('heading', { name: /First-run setup/i })
  try {
    await setupHeading.waitFor({ timeout: 8000 })
    await page.getByRole('button', { name: 'Finish setup' }).click()
  } catch {
    // Setup already completed.
  }

  const mainHeading = page.getByTestId('stMainBlockContainer').getByRole('heading', { name: 'qPCR Analysis' })
  await expect(mainHeading.first()).toBeVisible({ timeout: 60000 })
  await expect(page.getByText('Wells loaded')).toBeVisible()

  // Give Streamlit time to settle (tables render async); keep this short to avoid flakiness.
  await page.waitForTimeout(1500)

  // Avoid fullPage screenshots here; off-screen components can still be rendering and cause
  // subtle diffs between consecutive captures. Section snapshots below cover the full flow.
  await expect(page).toHaveScreenshot('example_run.png', { mask: [banner] })

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
