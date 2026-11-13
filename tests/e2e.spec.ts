import { test, expect } from '@playwright/test'

test('practice review keeps its step and retains explicit sample identities', async ({ page }) => {
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  await page.goto('/')
  await expect(page).toHaveTitle(/EasyLab Next/)
  await page.getByRole('button', { name:'Open practice experiment', exact:true }).click()
  await page.getByRole('button', { name:'Confirm columns', exact:true }).click()
  await expect(page.getByRole('heading',{name:'Map samples',exact:true})).toBeVisible()
  await page.getByRole('button', { name:'Confirm sample map', exact:true }).click()
  await expect(page.getByRole('heading',{name:'Review wells',exact:true})).toBeVisible()
  await page.getByRole('button', { name:'Confirm well review', exact:true }).click()
  await expect(page.getByRole('heading',{name:'Methods',exact:true})).toBeVisible()
  // No references/calibrators selected: the method must explain missing requirements.
  await page.getByRole('button', { name:'Confirm methods and standards', exact:true }).click()
  await page.getByRole('button', { name:'Calculate reviewed results', exact:true }).click()
  await expect(page.getByRole('heading',{name:'Results',exact:true})).toBeVisible()
  await expect(page.getByText('0 eligible sample–assay endpoints of 12. Counts are method-specific.',{exact:true})).toBeVisible()
  await expect(page.getByTestId('stException')).toHaveCount(0)
  expect(errors).toEqual([])
})

test('narrow-window collapsed navigation keeps keyboard focus visible', async ({ page }) => {
  await page.setViewportSize({width:768,height:900})
  await page.goto('/')
  await expect(page.getByRole('button',{name:'Open practice experiment',exact:true})).toBeVisible()
  await page.keyboard.press('Tab')
  const bounds = await page.locator(':focus').boundingBox()
  expect(bounds).not.toBeNull()
  expect(bounds!.x).toBeGreaterThanOrEqual(0)
  expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(768)
})
