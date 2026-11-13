import io
import json
from pathlib import Path
import pandas as pd
import pytest
from openpyxl import load_workbook
from qpcr_io import read_table, suggested_mapping, map_columns, excel_results
from qpcr_storage import new_run, encode_package, decode_package, save_run, save_preset, load_presets, storage_paths
from qpcr_session import changed, fresh_run
from qpcr_graphs import result_figure, static_figure

@pytest.fixture(autouse=True)
def isolated(tmp_path,monkeypatch):
    monkeypatch.setenv('EASYLAB_RUNTIME_ROOT',str(tmp_path/'next'))
    for key in ['EASYLAB_DATA_PATH','EASYLAB_SETTINGS_PATH','EASYLAB_EXPORTS_PATH']:
        monkeypatch.delenv(key,raising=False)

def test_native_results_preserves_undetermined_and_invalid():
    raw,details=read_table(b'Instrument metadata\nWell\tTarget Name\tSample Name\tCt\nA1\tIL1b\tAnimal 1\tUndetermined\nA2\tIL1b\tAnimal 1\tbad\n','results.txt')
    mapping=suggested_mapping(raw)
    assert mapping['Gene']=='Target Name' and mapping['Label']=='Sample Name'
    mapped=map_columns(raw,mapping)
    assert list(mapped.Cq)==['Undetermined','bad'] and len(mapped)==2

def test_duplicate_headers_blocked():
    with pytest.raises(ValueError,match='Duplicate'):
        read_table(b'Well,Gene,Label,Cq,Cq\nA1,T,S,1,2','bad.csv')

def test_package_restores_decisions_and_graphs_and_detects_edits():
    run=new_run();run.update(graph={'scale':'Log','title':'Expression'},wells=[{'Keep':False,'Reason':'Incorrect identity','OriginalCq':'Undetermined'}])
    decoded=decode_package(encode_package(run))
    assert decoded==run
    corrupted=json.loads(encode_package(run));corrupted['run']['name']='changed'
    with pytest.raises(ValueError,match='checksum'):
        decode_package(json.dumps(corrupted).encode())
    first=save_run(run);run['revision']+=1;run['name']='Changed';second=save_run(run)
    assert first!=second and first.exists()

def test_storage_rejects_external_and_symlink_paths(tmp_path,monkeypatch):
    monkeypatch.setenv('EASYLAB_DATA_PATH',str(tmp_path/'stable'))
    with pytest.raises(ValueError,match='separate'):
        storage_paths()
    monkeypatch.delenv('EASYLAB_DATA_PATH')
    alias=tmp_path/'alias';alias.symlink_to(tmp_path/'next',target_is_directory=True)
    monkeypatch.setenv('EASYLAB_RUNTIME_ROOT',str(alias))
    with pytest.raises(ValueError,match='symbolic'):
        storage_paths()

def test_new_run_and_presets_do_not_approve_experimental_decisions():
    run=new_run();run['settings'].update(references=['R'],calibrators=['S'],efficiency_confirmed=True);run['approvals']={'3':4};run['results']={'x':[1]}
    save_preset('Lab',run['settings'],{'scale':'Log'})
    assert 'references' not in load_presets()['Lab']['settings']
    fresh=fresh_run(run)
    assert not fresh['approvals'] and not fresh['settings']['references'] and not fresh['settings']['efficiency_confirmed']
    result=changed(run,{'wells':[{'Cq':15}]},2)
    assert not result['results'] and '3' not in result['approvals']
    result=changed(result,{'results':{'samples':[1]}},4)
    assert result['results']['samples']==[1]

def test_excel_treats_imported_formula_as_text():
    book=load_workbook(io.BytesIO(excel_results({'Samples':pd.DataFrame({'Label':['=1+1']})},{})))
    assert book['Samples']['A2'].value=='=1+1' and book['Samples']['A2'].data_type=='s'

def test_plot_and_exports_use_same_eligible_sample_values():
    frame=pd.DataFrame({'SampleId':['a','b','c','d'],'Label':['A','B','C','D'],'Group':['C','C','T','T'],'Value':[0.,1.,2.,9.],'Eligible':[True,True,True,False],'Status':['ok']*4})
    fig,points,omitted=result_figure(frame,{'scale':'Log','palette':{'C':'#123456'}})
    assert len(points)==2 and omitted==1 and list(points.Value)==[1.,2.]
    assert fig.layout.yaxis.type=='log'
    assert fig.data[0].marker.color=='#123456'
    assert len(static_figure(frame,{'scale':'Linear'},'png'))>1000
    assert len(static_figure(frame,{'scale':'Linear'},'svg'))>1000
    assert len(static_figure(frame,{'scale':'Linear'},'pdf'))>1000

def test_launcher_refuses_existing_port(tmp_path,monkeypatch):
    import os
    import socket
    import subprocess
    import sys
    probe=socket.socket()
    try:
        probe.bind(('127.0.0.1',18501))
    except OSError:
        probe.close()
        pytest.skip('Port already used by another task; leave it untouched')
    try:
        probe.listen()
        root=Path(__file__).resolve().parents[1]
        result=subprocess.run([sys.executable,str(root/'run_next.py')],cwd=root,capture_output=True,text=True,timeout=15,env=os.environ.copy())
        assert result.returncode != 0
        assert 'already in use' in result.stderr
    finally:
        probe.close()

def test_oversized_package_is_never_saved(tmp_path,monkeypatch):
    import qpcr_storage
    monkeypatch.setattr(qpcr_storage,'MAX_PACKAGE_BYTES',1024)
    run=new_run();run['sources']=[{'base64':'x'*2048}]
    with pytest.raises(ValueError,match='package limit'):
        encode_package(run)
    with pytest.raises(ValueError,match='package limit'):
        save_run(run)
    assert not list(tmp_path.rglob('*.easylab.json'))

def test_png_dimensions_match_settings():
    from PIL import Image
    frame=pd.DataFrame({'SampleId':['a','b'],'Label':['A','B'],'Group':['C','C'],'Value':[1.,2.],'Eligible':[True,True],'Status':['Ready','Ready']})
    data=static_figure(frame,{'width':1000,'height':600,'scale':'Linear'},'png')
    assert Image.open(io.BytesIO(data)).size == (1000,600)

def test_graph_auto_labels_follow_endpoint_but_custom_labels_are_preserved():
    from qpcr_session import graph_defaults
    original=graph_defaults({},'Standard curve','T','Target/reference ratio',['C','T'])
    changed_plot=graph_defaults(original,'Standard curve','T','Calibrator-relative fold change',['C','T'])
    assert changed_plot['ylabel']=='Calibrator-relative fold change'
    custom=dict(original,ylabel='Target expression',title='Treatment effect')
    changed_plot=graph_defaults(custom,'ΔΔCq','Other','Calibrator-relative fold change',['C','T'])
    assert changed_plot['ylabel']=='Target expression' and changed_plot['title']=='Treatment effect'

def test_revisiting_sample_map_preserves_assay_specific_standard_decisions():
    from qpcr_session import preserve_standard_decisions
    old=[{'Plate':'P','Gene':'G','Label':'S2','Concentration':123.,'Include':False,'Reason':'Preparation error'}]
    generated=[{'Plate':'P','Gene':'G','Label':'S1','Concentration':16.,'Include':True,'Reason':''},{'Plate':'P','Gene':'G','Label':'S2','Concentration':4.,'Include':True,'Reason':''}]
    merged=preserve_standard_decisions(generated,old)
    assert merged[1]==old[0]
    assert merged[0]==generated[0]
    series_change=preserve_standard_decisions(generated,old,preserve_concentrations=False)
    assert series_change[1]['Concentration']==4. and series_change[1]['Include'] is False

def test_nested_storage_links_cannot_redirect_experiment_writes(tmp_path):
    paths=storage_paths()
    stable=tmp_path/'existing-data';stable.mkdir()
    (paths['data']/'experiments').symlink_to(stable,target_is_directory=True)
    with pytest.raises(ValueError,match='Linked'):
        save_run(new_run())
    assert not list(stable.iterdir())
