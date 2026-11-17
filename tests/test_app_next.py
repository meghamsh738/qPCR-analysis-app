"""Guided-review flow and invalidation without a browser or user profile."""
from pathlib import Path
from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parents[1] / 'app.py')

def button(at,label):
    return next(x for x in at.button if x.label == label)

def field(elements,label):
    return next(x for x in elements if x.label == label)

def clean(at):
    assert not at.exception, [e.message for e in at.exception]

def fresh_page(at):
    # AppTest 1.50 retains stale form widgets after a rerun with navigation.
    # Recreate the page with the application's explicit state, like a clean rerender.
    clean(at)
    page=AppTest.from_file(APP,default_timeout=30)
    for key in ['run','undo','generation','step']:
        page.session_state[key]=at.session_state[key]
    return page.run()

def test_guided_practice_results_save_and_undo(tmp_path,monkeypatch):
    monkeypatch.setenv('EASYLAB_RUNTIME_ROOT', str(tmp_path/'next'))
    monkeypatch.setenv('MPLCONFIGDIR',str(tmp_path/'mpl'))
    for key in ['EASYLAB_DATA_PATH','EASYLAB_SETTINGS_PATH','EASYLAB_EXPORTS_PATH']:
        monkeypatch.delenv(key,raising=False)
    at=AppTest.from_file(APP,default_timeout=30).run();clean(at)
    button(at,'Open practice experiment').click().run();at=fresh_page(at);clean(at)
    button(at,'Confirm columns').click().run();at=fresh_page(at);clean(at)
    button(at,'Confirm sample map').click().run();at=fresh_page(at);clean(at)
    button(at,'Confirm well review').click().run();at=fresh_page(at);clean(at)
    field(at.multiselect,'Reference genes').set_value(['Reference'])
    field(at.multiselect,'Calibrator biological samples').set_value(['Practice-1','Practice-2','Practice-3'])
    field(at.checkbox,'I accept the near-100%, similar-efficiency assumption for ΔΔCq in this experiment').check()
    button(at,'Confirm methods and standards').click().run();at=fresh_page(at);clean(at)
    button(at,'Calculate reviewed results').click().run();at=fresh_page(at);clean(at)
    result=at.session_state['run']['results']
    assert result['analysisKey']
    assert len(result['ddcq_results'])==6
    assert sum(row['Eligible'] for row in result['ddcq_results'])==6
    button(at,'Continue to graphs and statistics').click().run();at=fresh_page(at);clean(at)
    field(at.checkbox,'The biological samples and design meet this test’s assumptions').check()
    button(at,'Run selected comparisons').click().run();at=fresh_page(at);clean(at)
    assert at.session_state['run']['statistics']['table']
    field(at.radio,'Analysis steps').set_value('Save & export').run();clean(at)
    button(at,'Save experiment locally').click().run();at=fresh_page(at);clean(at)
    assert list((tmp_path/'next').rglob('*.easylab.json'))
    field(at.radio,'Analysis steps').set_value('Methods').run();clean(at)
    field(at.multiselect,'Reference genes').set_value([])
    button(at,'Confirm methods and standards').click().run();at=fresh_page(at);clean(at)
    assert not at.session_state['run']['results']
    button(at,'Undo').click().run();at=fresh_page(at);clean(at)
    assert at.session_state['run']['results']['analysisKey']


def test_two_reference_sets_select_the_graph_and_open_export(monkeypatch, tmp_path):
    from qpcr_example import practice_data
    from qpcr_io import records
    from qpcr_storage import new_run, decode_package, encode_package
    from qpcr_workflow import prepare_wells, make_standard_map
    import pandas as pd

    monkeypatch.setenv('EASYLAB_RUNTIME_ROOT', str(tmp_path/'next'))
    monkeypatch.setenv('MPLCONFIGDIR', str(tmp_path/'mpl'))
    raw = practice_data()
    alternate = raw.loc[raw.Gene.eq('Reference')].copy()
    alternate['Gene'] = 'Alternate'
    alternate['Well'] = 'Alt-' + alternate.Well
    alternate = alternate.loc[~alternate.SampleId.eq('Practice-6')]
    wells = prepare_wells(pd.concat([raw, alternate], ignore_index=True))
    run = new_run()
    run.update(wells=records(wells), standards=records(make_standard_map(wells)))
    run['approvals'] = {'0': 1, '1': 1, '2': 1}
    at = AppTest.from_file(APP, default_timeout=45)
    at.session_state.run = run
    at.session_state.undo = []
    at.session_state.generation = 0
    at.session_state.step = 'Methods'
    at.run()
    field(at.number_input, 'Number of reference sets').set_value(2).run()
    field(at.text_input, 'Set name').set_value('Main reference')
    field(at.text_input, 'Set 2 name').set_value('Alternate reference')
    field(at.multiselect, 'Target assays').set_value(['Target'])
    field(at.multiselect, 'Reference genes').set_value(['Reference'])
    field(at.multiselect, 'Reference genes for set 2').set_value(['Alternate'])
    field(at.multiselect, 'Calibrator biological samples').set_value(['Practice-1', 'Practice-2', 'Practice-3'])
    field(at.checkbox, 'I accept the near-100%, similar-efficiency assumption for ΔΔCq in this experiment').check()
    field(at.checkbox, 'I accept the ΔΔCq efficiency assumption for set 2').check()
    button(at, 'Confirm methods and standards').click().run()
    at = fresh_page(at)
    button(at, 'Calculate reviewed results').click().run()
    at = fresh_page(at)
    stored = at.session_state.run
    reopened = decode_package(encode_package(stored))
    assert reopened == stored
    second = stored['settings']['reference_sets'][1]['id']
    field(at.selectbox, 'Reference set to review').set_value(second).run()
    button(at, 'Continue to graphs and statistics').click().run()
    at = fresh_page(at)
    assert field(at.selectbox, 'Reference set').value == second
    field(at.selectbox, 'Analysis method').set_value('ΔΔCq').run()
    at = fresh_page(at)
    field(at.selectbox, 'Assay').set_value('Target').run()
    at = fresh_page(at)
    assert any('Alternate reference · 5 biological-sample dots' in item.value for item in at.caption)
    field(at.radio, 'Analysis steps').set_value('Save & export').run()
    clean(at)
    field(at.selectbox, 'Reference set for the review workbook').set_value(second).run()
    clean(at)
