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
