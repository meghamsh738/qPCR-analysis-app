# EasyLab Next qPCR phase 1 — development review

13 November 2026. Fictional practice data; no real lab observations or saved accounts.

The screenshots show the guided Next workflow in background Chromium at 1400 × 900 and 768 × 900. They belong to the phase-1 code committed with this folder. The tested server was launched from this Next checkout on strict port 18501 with separate data/settings under `easylab-next/runtime`. All task-owned preview servers and browsers were closed after checking.

Verified: native Results preview and invalid-input handling; contextual sample identity; well/method review; both endpoint methods; biological graph points and selected basic comparisons; undo/invalidation; local and portable save/reopen; labels, colours, endpoint and comparison settings restored; narrow-window keyboard navigation; CSV/TSV, figures and Excel download. No browser console errors in the final walkthrough. Two automated browser tests passed.

The 60 Python tests passed, including legacy scientific regression checks, new identity/replicate/range cases, reference/calibrator weighting, statistics against SciPy reference routines, package/storage isolation and spreadsheet recalculation checks. The accepted review XLSX was recalculated using a separate headless LibreOffice profile and its three printed pages inspected. Native Excel, packaged macOS/installed-app behavior and Windows remain unverified.

The desktop screenshot is the current visual handoff; narrow screenshot shows the same workflow with the collapsible step list. Historical screenshot tests elsewhere in this repository describe the earlier application.
