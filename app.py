"""EasyLab Next qPCR: explicit review, biological samples and local experiments."""
from __future__ import annotations
import copy
import itertools
import json
from pathlib import Path
import pandas as pd
import streamlit as st
from qpcr_io import ALIASES, read_table, suggested_mapping, map_columns, records, source_record, excel_results
from qpcr_storage import new_run, storage_paths, encode_package, decode_package, save_run, load_presets, save_preset
from qpcr_session import STEPS, changed, fresh_run, analysis_key, graph_defaults, preserve_standard_decisions
from qpcr_example import practice_data
from qpcr_workflow import prepare_wells, summarise_samples, make_standard_map, fit_curves, analyze
from qpcr_statistics import run_statistics
from qpcr_graphs import result_figure, static_figure, curve_figure, PALETTE, analysis_report

st.set_page_config(page_title='EasyLab Next · qPCR', page_icon='🧬', layout='wide')
st.markdown('<style>.block-container{max-width:1300px;padding-top:4rem} h1{letter-spacing:-.035em} [data-testid="stSidebar"]{background:#eef5f5} .stCaption{line-height:1.5} [data-testid="stSidebar"][aria-expanded="false"]{visibility:hidden} button[data-testid="stExpandSidebarButton"]{width:auto!important;gap:4px;padding-right:8px} button[data-testid="stExpandSidebarButton"]::after{content:"Steps";font-size:.85rem}</style>', unsafe_allow_html=True)
try:
    paths = storage_paths()
except ValueError as exc:
    st.error(str(exc))
    st.stop()
if 'run' not in st.session_state:
    st.session_state.run = new_run()
    st.session_state.undo = []
    st.session_state.generation = 0
if 'next_step' in st.session_state:
    st.session_state.step = st.session_state.pop('next_step')
run = st.session_state.run
key = f"{run['runId']}-{run['revision']}-{st.session_state.generation}"

def frame(name):
    return pd.DataFrame(run.get(name, []))

def apply(updates, stage, next_step=None, approve=True):
    candidate = changed(run, updates, stage)
    if not approve:
        candidate['approvals'].pop(str(stage),None)
    try:
        encode_package(candidate)
    except ValueError as exc:
        st.error(str(exc))
        st.stop()
    st.session_state.undo = (st.session_state.undo + [copy.deepcopy(run)])[-15:]
    st.session_state.run = candidate
    if next_step is not None:
        st.session_state.next_step = STEPS[next_step]
    st.rerun()

def replace_run(value):
    st.session_state.undo = (st.session_state.undo + [copy.deepcopy(run)])[-15:]
    st.session_state.run = value
    st.session_state.generation += 1
    st.session_state.next_step = STEPS[0]
    st.rerun()

def need_wells():
    if not run['wells']:
        st.info('Import and confirm a well table first.')
        st.stop()

def result_tables():
    return {name: pd.DataFrame(value) for name, value in run['results'].items() if isinstance(value, list)}

def csv_download(label, data, filename, sep=','):
    st.download_button(label, data.to_csv(index=False,sep=sep), filename, 'text/csv' if sep==',' else 'text/tab-separated-values')

def import_data(files):
    tables, sources = [], []
    names = [f[0] for f in files]
    if len(names) != len(set(names)):
        raise ValueError('Files have repeated names. Give each plate file a distinct name before importing.')
    for filename, data in files:
        table, details = read_table(data, filename)
        table['OriginalFile'] = filename
        if not suggested_mapping(table)['Plate']:
            table['Plate'] = Path(filename).stem
        tables.append(table)
        sources.append(source_record(filename, data, details))
    raw = pd.concat(tables, ignore_index=True).fillna('')
    settings = copy.deepcopy(run['settings'])
    settings.update(references=[],calibrators=[],allow_extrapolation=False,efficiency_confirmed=False)
    apply({'sources': sources, 'raw': records(raw), 'wells': [], 'standards': [], 'settings':settings, 'practice':False}, 0, approve=False)

with st.sidebar:
    st.caption('EASYLAB NEXT · SEPARATE WORKSPACE')
    st.title('qPCR Analysis')
    st.caption(f"{run['name']} · revision {run['revision']}")
    st.radio('Analysis steps', STEPS, key='step', format_func=lambda s: f'{STEPS.index(s)+1}. {s}')
    a,b = st.columns(2)
    if a.button('Undo', disabled=not st.session_state.undo, width='stretch'):
        st.session_state.run = st.session_state.undo.pop()
        st.session_state.generation += 1
        st.rerun()
    with b.popover('Reset'):
        st.caption('Reset review status and dependent results; keep imported values.')
        if st.button('Reset step review'):
            stage = STEPS.index(st.session_state.step)
            result = changed(run, {}, stage)
            result['approvals'].pop(str(stage), None)
            if stage <= 4:
                result['results'] = {}
                result['statistics'] = {}
            if stage == 0:
                result.update(raw=[], sources=[], wells=[], standards=[])
            replace_run(result)
        confirm_reset = st.checkbox('Start a new experiment')
        if st.button('Confirm new experiment', disabled=not confirm_reset):
            replace_run(fresh_run(run))
    with st.expander('Open saved experiment'):
        package = st.file_uploader('Experiment package', type=['json'], key='open_package')
        if package and st.button('Open package'):
            try:
                replace_run(decode_package(package.getvalue()))
            except (ValueError, TypeError) as exc:
                st.error(str(exc))
        directory = paths['data'] / 'experiments'
        saved = sorted([p for p in directory.glob('*.easylab.json') if not p.is_symlink()], key=lambda p:p.stat().st_mtime, reverse=True) if directory.exists() else []
        if saved:
            selected = st.selectbox('Saved locally', [p.name for p in saved])
            if st.button('Open local experiment'):
                try:
                    replace_run(decode_package((directory / selected).read_bytes()))
                except (ValueError, TypeError) as exc:
                    st.error(str(exc))
    with st.expander('Local storage'):
        st.caption('Experiments and preferences are stored only in the separate Next folders below. Download a package for a portable backup.')
        st.code(str(paths['data']), language=None)
        st.code(str(paths['settings']), language=None)
    st.caption('Practice data are fictional. Originals are kept inside each saved experiment package.')

stage = STEPS.index(st.session_state.step)
st.caption(' · '.join(f"{'✓' if str(i) in run['approvals'] else '○'} {name}" for i,name in enumerate(STEPS)))
st.caption(f'STEP {stage+1} OF 7')
st.title(STEPS[stage])
if run.get('practice'):
    st.info('Practice experiment · fictional observations for learning the workflow.')

if stage == 0:
    st.write('Bring in the instrument export, then check which columns describe each well. Every row stays available for review.')
    left,right = st.columns([3,1])
    with left:
        files = st.file_uploader('Instrument Results or well tables', type=['csv','tsv','txt','xlsx'], accept_multiple_files=True)
        if st.button('Preview files', disabled=not files):
            try:
                import_data([(f.name, f.getvalue()) for f in files])
            except (ValueError, UnicodeError, pd.errors.ParserError) as exc:
                st.error(f'Could not read the table: {exc}')
    with right:
        st.write('Try the full workflow')
        if st.button('Open practice experiment'):
            example = practice_data()
            value = fresh_run(run)
            value.update(name='Practice · target expression', practice=True, raw=records(example), sources=[source_record('practice.csv',example.to_csv(index=False).encode(),'Fictional example')])
            replace_run(value)
    with st.expander('Paste a CSV or tab-separated table'):
        pasted = st.text_area('Table with header row')
        if st.button('Preview pasted table', disabled=not pasted.strip()):
            try:
                import_data([('Pasted plate.tsv',pasted.encode())])
            except (ValueError, pd.errors.ParserError) as exc:
                st.error(str(exc))
    raw = frame('raw')
    if not raw.empty:
        st.subheader(f'{len(raw)} imported rows')
        st.dataframe(raw, hide_index=True, width='stretch')
        mapping = suggested_mapping(raw)
        with st.form(f'columns-{key}'):
            selections = {}
            columns = st.columns(3)
            for i,(target,default) in enumerate(mapping.items()):
                options = [''] + list(raw.columns)
                selections[target] = columns[i%3].selectbox(target, options, index=options.index(default), help='Blank optional fields remain unknown. Biological SampleId is different from the instrument label.')
            st.caption('Required: Well, Gene, Label and Cq. If Type is blank, labels suggesting standards or controls are proposed automatically; confirm them at the next step.')
            if st.form_submit_button('Confirm columns', type='primary'):
                try:
                    wells = prepare_wells(map_columns(raw, selections))
                    apply({'wells': records(wells), 'standards': []}, 0, 1)
                except ValueError as exc:
                    st.error(str(exc))

elif stage == 1:
    need_wells()
    st.write('Give each biological sample a stable identity. Reused labels on different plates stay separate unless you deliberately assign the same SampleId. Use a new identity for a different tissue or batch. For paired measurements, give each condition its own SampleId and share a PairId across the matched samples.')
    wells = frame('wells')
    identity_columns = ['Plate','Label','Tissue','Batch','Type','SampleId','Group','PairId']
    for context in ['Tissue','Batch']:
        if context not in wells:
            wells[context] = ''
    identity = wells[identity_columns].fillna('').drop_duplicates().reset_index(drop=True)
    with st.expander('Add a separate sample-information table'):
        metadata = st.file_uploader('Sample information (Plate, Label, SampleId, Group, PairId, Tissue, Batch)', type=['csv','tsv','txt','xlsx'])
        if metadata:
            try:
                info,_ = read_table(metadata.getvalue(),metadata.name)
                st.dataframe(info,hide_index=True)
                mode = st.radio('Match records by', ['Plate and Label','Plate and Well','Plate, Label, Tissue and Batch','Label across plates (explicit)'])
                on = {'Plate and Label':['Plate','Label'],'Plate and Well':['Plate','Well'],'Plate, Label, Tissue and Batch':['Plate','Label','Tissue','Batch'],'Label across plates (explicit)':['Label']}[mode]
                if not all(c in info for c in on):
                    st.warning('The matching columns must have the exact headings shown above.')
                elif info.duplicated(on).any():
                    st.error('More than one sample-information row matches an identity. Resolve these duplicates first.')
                else:
                    joined = wells[on].drop_duplicates().merge(info,on=on,how='left',indicator=True)
                    unmatched = int((joined['_merge'] == 'left_only').sum())
                    st.caption(f'{len(joined)-unmatched} mapped identities; {unmatched} without a matching record.')
                    if st.button('Apply sample information'):
                        lookup = info.set_index(on)
                        for index,row in wells.iterrows():
                            token = tuple(row[c] for c in on) if len(on)>1 else row[on[0]]
                            if token in lookup.index:
                                context_changed = any(c in info and str(lookup.loc[token,c]) != str(row.get(c,'')) for c in ['Tissue','Batch'])
                                for c in ['SampleId','Group','PairId','Type','Tissue','Batch']:
                                    if c in info:
                                        wells.at[index,c] = str(lookup.loc[token,c])
                                if context_changed and 'SampleId' not in info and str(row.SampleId).startswith('auto:'):
                                    wells.at[index,'SampleId'] = ''
                        sources = run['sources']+[source_record(metadata.name,metadata.getvalue(),'Sample information')]
                        apply({'wells':records(prepare_wells(wells)), 'sources':sources},0)
            except (ValueError, KeyError, pd.errors.ParserError) as exc:
                st.error(f'Sample information needs correction: {exc}')
    with st.form(f'identities-{key}'):
        edited = st.data_editor(identity, disabled=['Plate','Label'], num_rows='fixed', hide_index=True, width='stretch',
            column_config={'Type':st.column_config.SelectboxColumn(options=['Sample','Standard','Negative','Positive','Control','Blank']), 'SampleId':st.column_config.TextColumn(help='The same biological sample across target and reference assays must have the same identity.'), 'PairId':st.column_config.TextColumn(help='Only fill this for explicitly matched subjects in a paired comparison.')})
        if st.form_submit_button('Confirm sample map', type='primary'):
            original_wells = wells[identity_columns].fillna('')
            for index, row in edited.iterrows():
                prior = identity.loc[index]
                mask = (original_wells == prior).all(axis=1)
                for c in ['Type','SampleId','Group','PairId','Tissue','Batch']:
                    wells.loc[mask,c] = row[c]
                if any(row[c] != prior[c] for c in ['Tissue','Batch']) and row.SampleId == prior.SampleId and str(prior.SampleId).startswith('auto:'):
                    wells.loc[mask,'SampleId'] = ''
            prepared = prepare_wells(wells)
            generated = records(make_standard_map(prepared,top=run['settings']['top'],factor=run['settings']['factor']))
            apply({'wells':records(prepared),'standards':preserve_standard_decisions(generated,run['standards'])},1,2)
    st.caption('Standards, controls and samples are counted separately. Unresolved rows remain visible in well review.')
    st.dataframe(wells.groupby(['Plate','Type'],dropna=False).size().reset_index(name='Wells'),hide_index=True)

elif stage == 2:
    need_wells()
    st.write('Review technical wells before combining biological samples. A warning asks for your attention; it does not remove a measurement.')
    wells = frame('wells')
    st.caption('Exclusions below apply only to the selected well. Give a short reason for each exclusion; excluded and undetermined observations remain in the package.')
    columns = [c for c in ['RowId','Plate','Well','Gene','Type','Label','SampleId','OriginalCq','CqStatus','Issue','Keep','Reason'] if c in wells]
    with st.form(f'wells-{key}'):
        review_display = wells[columns].copy()
        review_display['OriginalCq'] = review_display['OriginalCq'].fillna('').astype(str)
        edited = st.data_editor(review_display, disabled=[c for c in columns if c not in ['Keep','Reason']],num_rows='fixed',hide_index=True,width='stretch')
        if st.form_submit_button('Confirm well review',type='primary'):
            wells[['Keep','Reason']] = edited[['Keep','Reason']]
            if ((~wells.Keep.astype(bool)) & wells.Reason.fillna('').str.strip().eq('')).any():
                st.error('Give a reason for every excluded well.')
            else:
                apply({'wells':records(prepare_wells(wells))},2,3)
    with st.expander('Original instrument fields'):
        st.dataframe(frame('raw').astype(str),hide_index=True,width='stretch')
    st.subheader('Biological sample summaries')
    st.dataframe(summarise_samples(wells,threshold=run['settings']['threshold']), hide_index=True, width='stretch')
    st.caption('Technical wells are averaged first. One retained numeric well may be used with a note; all-undetermined samples have no endpoint.')

elif stage == 3:
    need_wells()
    st.write('Choose the standards, reference genes and calibrators for this experiment. Standard-curve ratios and calibrator-relative fold changes are separate results.')
    with st.expander('How the two methods calculate expression'):
        st.write('Standard curves: average retained technical Cqs, invert the plate-specific curve, then divide target quantity by the geometric mean of the selected reference quantities. A separate fold-change endpoint divides this ratio by the arithmetic mean of the selected biological calibrators.')
        st.write('ΔΔCq: subtract the mean reference Cq from each target mean Cq, subtract the equally weighted calibrator mean ΔCq, then calculate 2 to the power of minus ΔΔCq. This assumes near-100%, similar amplification efficiencies. The calibrator fold changes have a geometric mean of one; their arithmetic mean need not be one.')
        st.write('A missing or withheld reference prevents that sample’s normalized endpoint. Range and efficiency warnings remain visible and must be reviewed before conclusions are drawn.')
    wells = frame('wells')
    settings = copy.deepcopy(run['settings'])
    with st.expander('Preset and editable standard series'):
        try:
            presets = load_presets()
        except (ValueError, OSError):
            presets = {}
            st.warning('Saved presets could not be read; the current experiment is preserved.')
        names = ['Meghamsh — reviewed qPCR','Legacy — 256, fourfold']+list(presets)
        name = st.selectbox('Preset',names)
        if st.button('Apply preset settings'):
            selected = presets.get(name,{}).get('settings',{'top':256. if name.startswith('Legacy') else 16384.,'factor':4.,'threshold':.75})
            settings.update(selected,preset=name)
            apply({'settings':settings,'standards':preserve_standard_decisions(records(make_standard_map(wells,top=settings['top'],factor=settings['factor'])),run['standards'],preserve_concentrations=False)},2,approve=False)
        with st.form(f'series-{key}'):
            top = st.number_input('Top standard (arbitrary units)',min_value=.000001,value=float(settings['top']),format='%.6f')
            factor = st.number_input('Serial dilution factor',min_value=1.000001,value=float(settings['factor']))
            threshold = st.number_input('Technical-well Cq difference warning',min_value=0.,value=float(settings['threshold']),step=.05)
            efficiency_min = st.number_input('Lower efficiency warning (%)',min_value=0.,value=float(settings['efficiency_min']))
            efficiency_max = st.number_input('Upper efficiency warning (%)',min_value=0.,value=float(settings['efficiency_max']))
            efficiency_difference = st.number_input('Between-assay efficiency difference warning (percentage points)',min_value=0.,value=float(settings['efficiency_difference']))
            if st.form_submit_button('Apply series and warning threshold'):
                if efficiency_max <= efficiency_min:
                    st.error('The upper efficiency warning must be greater than the lower warning.')
                    st.stop()
                settings.update(top=top,factor=factor,threshold=threshold,efficiency_min=efficiency_min,efficiency_max=efficiency_max,efficiency_difference=efficiency_difference)
                apply({'settings':settings,'standards':preserve_standard_decisions(records(make_standard_map(wells,top=top,factor=factor)),run['standards'],preserve_concentrations=False)},2,approve=False)
    standards = frame('standards')
    if standards.empty:
        standards = make_standard_map(wells,top=settings['top'],factor=settings['factor'])
    curves, points = fit_curves(wells,standards)
    if not curves.empty:
        efficiencies = pd.to_numeric(curves.EfficiencyPercent,errors='coerce').dropna()
        if ((efficiencies < settings['efficiency_min']) | (efficiencies > settings['efficiency_max'])).any():
            st.warning('One or more fitted efficiencies fall outside your warning limits. Check dilution preparation and assay behavior before accepting the ΔΔCq assumptions.')
        if len(efficiencies)>1 and efficiencies.max()-efficiencies.min()>settings['efficiency_difference']:
            st.warning('Assay efficiencies differ by more than your warning threshold. The equal-efficiency assumption needs explicit scientific review; accepting it does not correct that difference.')
    with st.form(f'methods-{key}'):
        st.subheader('Assay-specific standard selection')
        st.caption('Edit quantities and include any suitable subset, including non-contiguous levels. Excluded standards need a reason.')
        edited = st.data_editor(standards,disabled=[c for c in standards if c not in ['Concentration','Include','Reason']],num_rows='fixed',hide_index=True,width='stretch')
        genes = sorted(wells.Gene.dropna().unique().tolist())
        refs = st.multiselect('Reference genes',genes,default=[r for r in settings['references'] if r in genes])
        sample_ids = sorted(wells.loc[wells.Type.eq('Sample'),'SampleId'].dropna().unique().tolist())
        calibrators = st.multiselect('Calibrator biological samples',sample_ids,default=[r for r in settings['calibrators'] if r in sample_ids],help='Select every biological sample in the intended baseline group. Each sample receives equal weight.')
        st.caption('Several reference genes are combined by mean Cq for ΔΔCq and geometric mean quantity for standard curves. Choose stable references for the tissues and conditions being compared.')
        efficiency = st.checkbox('I accept the near-100%, similar-efficiency assumption for ΔΔCq in this experiment',value=settings['efficiency_confirmed'])
        st.caption('This records your decision; it does not validate efficiency. Review assay efficiencies and calibrated ranges below. A mismatch needs scientific review.')
        extrapolation = st.checkbox('Allow explicitly labelled provisional extrapolation outside the standard curve',value=settings['allow_extrapolation'])
        if st.form_submit_button('Confirm methods and standards',type='primary'):
            if not edited.empty and ((~edited.Include.astype(bool))&edited.Reason.fillna('').str.strip().eq('')).any():
                st.error('Give a reason for every excluded standard level.')
            else:
                settings.update(references=refs,calibrators=calibrators,efficiency_confirmed=efficiency,allow_extrapolation=extrapolation)
                apply({'settings':settings,'standards':records(edited)},3,4)
    if not standards.empty:
        curves, points = fit_curves(wells,standards)
        st.subheader('Curve diagnostics')
        st.dataframe(curves,hide_index=True,width='stretch')
        st.dataframe(points,hide_index=True,width='stretch')
        for index,curve in curves.iterrows():
            selected = points.copy()
            for c in ['Plate','Gene']:
                if c in selected and c in curve:
                    selected = selected.loc[selected[c]==curve[c]]
            try:
                st.plotly_chart(curve_figure(selected,curve),use_container_width=True,key=f'curve-{index}')
            except (KeyError,ValueError,StopIteration) as exc:
                st.caption(f'Curve graph unavailable: {exc}')
    with st.expander('Save these defaults as a named preset'):
        preset_name = st.text_input('Preset name')
        if st.button('Save preset'):
            try:
                save_preset(preset_name,settings,run['graph'])
                st.success('Preset saved. Sample maps, references, calibrators and review approvals are not inherited by new experiments.')
            except ValueError as exc:
                st.error(str(exc))

elif stage == 4:
    need_wells()
    ready = all(str(s) in run['approvals'] for s in [0,1,2,3])
    st.write('Compare both methods at the biological-sample level. Missing requirements and range flags stay beside the results.')
    if not ready:
        st.warning('Confirm the column map, sample identities, well review and methods before calculating results.')
    if st.button('Calculate reviewed results',type='primary',disabled=not ready):
        try:
            s = run['settings']
            calculated = analyze(frame('wells'),frame('standards'),s['references'],s['calibrators'],allow_extrapolation=s['allow_extrapolation'],efficiency_confirmed=s['efficiency_confirmed'],threshold=s['threshold'])
            output = {k:records(v) if isinstance(v,pd.DataFrame) else v for k,v in calculated.items()}
            output['analysisKey'] = analysis_key(run)
            apply({'results':output},4)
        except ValueError as exc:
            st.error(str(exc))
    tables = result_tables()
    if tables:
        for title,name in [('Standard-curve results','standard_results'),('ΔΔCq results','ddcq_results')]:
            st.subheader(title)
            table = tables.get(name,pd.DataFrame())
            st.dataframe(table,hide_index=True,width='stretch')
            if not table.empty and 'Eligible' in table:
                st.caption(f'{int(table.Eligible.sum())} eligible sample–assay endpoints of {len(table)}. Counts are method-specific.')
        for name in ['issues','samples','quantities']:
            if name in tables:
                with st.expander(name.replace('_',' ').title()):
                    st.dataframe(tables[name],hide_index=True,width='stretch')
        if st.button('Continue to graphs and statistics'):
            st.session_state.next_step = STEPS[5]
            st.rerun()

elif stage == 5:
    tables = result_tables()
    if not tables:
        st.info('Calculate reviewed results first.')
        st.stop()
    st.write('Each dot represents one biological sample. Changing the axes changes the view, not the values used for statistical testing.')
    spec = copy.deepcopy(run['graph'])
    methods = {'Standard curve':'standard_results','ΔΔCq':'ddcq_results'}
    method = st.selectbox('Analysis method',list(methods),index=list(methods).index(spec.get('method','Standard curve')))
    values = tables.get(methods[method],pd.DataFrame())
    if values.empty or 'Gene' not in values:
        st.info('This method has no results. Revisit its requirements and reference choices.')
        st.stop()
    genes = sorted(values.Gene.unique())
    gene = st.selectbox('Assay',genes,index=genes.index(spec['gene']) if spec.get('gene') in genes else 0)
    values = values[values.Gene.eq(gene)].copy()
    # Some methods expose both a normalized ratio and a calibrator-relative endpoint.
    if 'Endpoint' in values:
        endpoints = list(values.Endpoint.unique())
        endpoint = st.selectbox('Endpoint',endpoints,index=endpoints.index(spec['endpoint']) if spec.get('endpoint') in endpoints else 0)
        values = values[values.Endpoint.eq(endpoint)]
    else:
        endpoint = 'Expression'
    groups = sorted(values.Group.fillna('').unique().tolist())
    defaults = graph_defaults(spec,method,gene,endpoint,groups)
    if defaults != run['graph']:
        apply({'graph':defaults},5)
    spec = defaults
    if spec['ylabel'] != endpoint or spec['title'] != gene:
        st.caption(f'Custom labels are in use. Plotted endpoint: {method} · {gene} · {endpoint}.')
    with st.form(f'plot-{key}'):
        cols = st.columns(3)
        title = cols[0].text_input('Graph title',value=spec.get('title',gene))
        ylabel = cols[1].text_input('Y-axis label',value=spec.get('ylabel',endpoint))
        scale = cols[2].selectbox('Y-axis', ['Linear','Log'],index=['Linear','Log'].index(spec.get('scale','Linear')))
        error = cols[0].selectbox('Uncertainty', ['None','SD','SEM','95% CI'],index=['None','SD','SEM','95% CI'].index(spec.get('error','None')))
        point_size = cols[1].slider('Point size',4,16,int(spec.get('point_size',8)))
        paired_lines = cols[2].checkbox('Connect explicitly matched pairs',value=spec.get('paired',False))
        order = st.multiselect('Group order (select all groups in the desired order)',groups,default=spec.get('group_order',groups) if set(spec.get('group_order',groups))==set(groups) else groups)
        width = cols[0].number_input('Export width (px)',min_value=400,max_value=2400,value=int(spec.get('width',1000)),step=100)
        height = cols[1].number_input('Export height (px)',min_value=300,max_value=1800,value=int(spec.get('height',600)),step=100)
        st.caption('PNG uses the specified pixel dimensions; PDF print size uses 100 pixels per inch. SVG remains scalable.')
        palette = {}
        for i,g in enumerate(groups):
            palette[g] = cols[i%3].color_picker(f'{g or "Ungrouped"} colour',value=spec.get('palette',{}).get(g,PALETTE[i%len(PALETTE)]))
        if st.form_submit_button('Apply graph settings'):
            if set(order) != set(groups):
                st.error('Include every group in the display order.')
            else:
                apply({'graph':dict(method=method,gene=gene,endpoint=endpoint,title=title,ylabel=ylabel,scale=scale,error=error,point_size=point_size,paired=paired_lines,group_order=order,width=width,height=height,palette=palette)},5)
    view = dict(spec,method=method,gene=gene,endpoint=endpoint)
    try:
        fig, points, omitted = result_figure(values,view)
        st.plotly_chart(fig,use_container_width=True)
        st.caption(f'{len(points)} biological-sample dots. Excluded, missing or withheld endpoints do not contribute.')
        if omitted:
            st.warning(f'{omitted} non-positive values cannot appear on a logarithmic axis. They remain in the data and statistical tests.')
        caption = f"{run['name']} · revision {run['revision']} · {method} · {gene} · {endpoint}"
        cols = st.columns(4)
        for col,fmt in zip(cols,['png','svg','pdf']):
            col.download_button(f'Download {fmt.upper()}',static_figure(values,view,fmt,caption),f'qpcr-{gene}.{fmt}')
        with cols[3]:
            csv_download('Graph data',points,'qpcr-graph-data.csv')
    except ValueError as exc:
        st.error(str(exc))
    st.subheader('Basic statistical comparisons')
    st.caption('Choose one endpoint and the complete comparison family before testing. Log axes do not log-transform the analysis. Multigroup repeated-measures and two-factor designs should be exported to a suitable statistics application.')
    with st.form(f'stats-{key}'):
        saved_stats = run['statistics'] if all(run['statistics'].get(k)==v for k,v in [('method',method),('gene',gene),('endpoint',endpoint)]) else {}
        designs = ['Welch · two independent groups','Paired · two explicitly matched groups','One-way ANOVA · independent groups, one factor']
        test = st.selectbox('Study design', designs,index={'welch':0,'paired':1,'anova':2}.get(saved_stats.get('test'),0))
        choices = list(itertools.combinations([g for g in groups if g],2))
        saved_comparisons = [tuple(pair) for pair in saved_stats.get('comparisons',[]) if tuple(pair) in choices]
        comparisons = st.multiselect('Selected comparisons (one family)',choices,default=saved_comparisons or choices[:1],format_func=lambda x:f'{x[0]} vs {x[1]}')
        assumptions = st.checkbox('The biological samples and design meet this test’s assumptions',value=bool(saved_stats))
        st.caption('Welch: independent samples and approximately normal values. Paired: explicit matching and approximately normal within-pair differences. ANOVA: independent observations, approximately normal residuals and similar group variances. Review small samples carefully.')
        if st.form_submit_button('Run selected comparisons',type='primary'):
            if not assumptions or not comparisons:
                st.error('Confirm the design and choose at least one comparison.')
            else:
                try:
                    code = 'welch' if test.startswith('Welch') else 'paired' if test.startswith('Paired') else 'anova'
                    if code == 'paired' and len({group for pair in comparisons for group in pair}) > 2:
                        raise ValueError('Repeated-measures comparisons across more than two groups need a suitable external statistics application. Select one matched pair of groups here.')
                    omnibus, comparison_table = run_statistics(values,method=code,comparisons=comparisons,group_order=groups)
                    apply({'statistics':{'method':method,'gene':gene,'endpoint':endpoint,'test':code,'comparisons':[list(p) for p in comparisons],'omnibus':omnibus,'table':records(comparison_table),'analysisKey':analysis_key(run)}},5)
                except ValueError as exc:
                    st.error(str(exc))
    statistics = run['statistics']
    if statistics and all(statistics.get(k)==v for k,v in [('method',method),('gene',gene),('endpoint',endpoint)]):
        st.caption(f"{statistics['test'].title()} · {method} · {gene} · {endpoint}")
        if statistics['test'] == 'anova':
            omnibus = statistics['omnibus']
            st.write(f"Overall ANOVA: F = {omnibus.get('Statistic')}, global P = {omnibus.get('PValue')}. Selected comparisons below use the pooled residual variance.")
        table = pd.DataFrame(statistics['table'])
        raw = st.checkbox('Show clearly labelled unadjusted P values')
        display = table if raw else table.drop(columns=[c for c in table if c.lower() in ['p','pvalue','p_value','p_unadjusted','rawp']],errors='ignore')
        st.dataframe(display,hide_index=True,width='stretch',column_config={'PAdjustedHolm':st.column_config.NumberColumn('P (Holm adjusted)',format='%.4g'),'PValue':st.column_config.NumberColumn('P (unadjusted)',format='%.4g'),'EffectMeanDifference':st.column_config.NumberColumn('Mean difference',format='%.5g'),'CI95Lower':st.column_config.NumberColumn('95% CI lower',format='%.5g'),'CI95Upper':st.column_config.NumberColumn('95% CI upper',format='%.5g')})
        st.caption('Holm-adjusted P values account for all selected comparisons. Confidence intervals describe each individual effect and are not simultaneous intervals.')
    with st.expander('Reference-gene review'):
        samples = tables.get('samples',pd.DataFrame())
        if not samples.empty:
            references = samples[samples.Gene.isin(run['settings']['references'])]
            st.dataframe(references,hide_index=True,width='stretch')
            st.caption('Compare retained mean Cq and disagreement across groups. Stable reference expression needs experimental support; these summaries alone do not establish suitability.')
            if not references.empty:
                reference = st.selectbox('Reference assay to plot', list(references.Gene.unique()))
                ref_data = references[references.Gene.eq(reference)].copy()
                ref_data['Value'] = ref_data['MeanCq']
                ref_fig, _, _ = result_figure(ref_data, {'title':f'{reference} · reference Cq by group','ylabel':'Mean Cq','scale':'linear'})
                st.plotly_chart(ref_fig,use_container_width=True)

else:
    st.write('Save a reopenable experiment with source copies, decisions, results and graph settings. Earlier saved revisions remain available.')
    with st.form(f'name-{key}'):
        name = st.text_input('Experiment name',value=run['name'])
        if st.form_submit_button('Apply name'):
            apply({'name':name.strip() or 'Untitled experiment'},6)
    left,right = st.columns(2)
    if left.button('Save experiment locally',type='primary'):
        try:
            path = save_run(run)
            st.success(f'Saved revision {run["revision"]}')
            st.code(str(path),language=None)
        except (ValueError,OSError) as exc:
            st.error(f'Could not save: {exc}')
    right.download_button('Download experiment package',encode_package(run),'qpcr-experiment.easylab.json','application/json')
    tables = result_tables()
    if tables:
        st.subheader('Reviewed results')
        if run['graph'].get('method'):
            st.download_button('Download selected-results PDF report',analysis_report(run),'qpcr-report.pdf','application/pdf')
        table_format = st.radio('Sample table format',['CSV','TSV'],horizontal=True)
        detailed = st.checkbox('Include detailed wells, curves and review material',value=True)
        selected = {k:v for k,v in tables.items() if detailed or k in ['standard_results','ddcq_results']}
        if detailed:
            selected['Reviewed wells'] = frame('wells')
            selected['Standards'] = frame('standards')
        selected['Statistics'] = pd.DataFrame(run['statistics'].get('table',[]))
        workbook_style = st.radio('Excel layout',['Assay-by-assay manual review','Detailed result tables'],horizontal=True)
        if workbook_style == 'Assay-by-assay manual review':
            from qpcr_review_export import cyc1_review_workbook
            workbook = cyc1_review_workbook(frame('wells'),frame('standards'),tables,run['settings'])
            st.caption('Editable assay sheets and a compact normalised sheet. Original source values remain below each review table. Excel edits do not change the saved app experiment. Open the experiment package to continue here. To bring spreadsheet edits back, export a reviewed flat well table and confirm its mappings again.')
        else:
            workbook = excel_results(selected,dict(run['settings'],statistics=run['statistics'],graph=run['graph'],revision=run['revision'],runId=run['runId']))
        st.download_button('Download Excel review workbook',workbook,'qpcr-review.xlsx','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        for name,table in selected.items():
            if not table.empty:
                csv_download(f'Download {name.replace("_"," ")}',table,f'{name}.{table_format.lower()}',sep=',' if table_format=='CSV' else '\t')
    else:
        st.info('You can save work in progress. Calculate reviewed results to enable result exports.')
    st.caption('The package checksum checks file integrity. Original imported files remain unchanged; their copies are included in the package. Portable packages may contain sample identities.')
