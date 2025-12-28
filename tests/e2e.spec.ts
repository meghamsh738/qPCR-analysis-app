import { test, expect } from '@playwright/test'

async function snapSection(page, headingRegex: RegExp, filename: string) {
  const heading = page.getByRole('heading', { name: headingRegex })
  await expect(heading).toBeVisible()
  await heading.scrollIntoViewIfNeeded()
  await page.waitForTimeout(600)
  await expect(page).toHaveScreenshot(filename)
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
  await snapSection(page, /4\)\s*Fit standard curves/i, 'curves.png')
  await snapSection(page, /6\)\s*Normalize to reference gene/i, 'quant_normalize.png')
  await snapSection(page, /7\)\s*Export/i, 'export.png')
})
