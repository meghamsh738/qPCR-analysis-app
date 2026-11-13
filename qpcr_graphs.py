"""Figures and static exports use the same reviewed rows and display specification."""
from __future__ import annotations
import hashlib
import io
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy import stats

PALETTE = ['#087f8c', '#c66b3d', '#5c62a5', '#73934c', '#a95883', '#58778d']


def graph_rows(frame, spec):
    sub = frame.copy()
    if 'Eligible' in sub:
        sub = sub[sub['Eligible'].fillna(False).astype(bool)]
    sub = sub[pd.to_numeric(sub['Value'], errors='coerce').notna()].copy()
    if sub['SampleId'].duplicated().any():
        raise ValueError('A graph must contain one value per biological sample. Select one assay and method.')
    omitted = 0
    if spec.get('scale', 'linear').lower() == 'log':
        omitted = int((sub['Value'] <= 0).sum())
        sub = sub[sub['Value'] > 0].copy()
    return sub, omitted


def _error(values, mode):
    if mode == 'None' or len(values) < 2:
        return 0.0
    sd = float(np.std(values, ddof=1))
    if mode == 'SD':
        return sd
    sem = sd / math.sqrt(len(values))
    return sem if mode == 'SEM' else sem * float(stats.t.ppf(.975, len(values)-1))


def _groups(frame, spec):
    present = frame['Group'].fillna('').replace('', 'Ungrouped').astype(str).unique().tolist()
    return list(dict.fromkeys([g for g in spec.get('group_order', []) if g in present] + present))


def _jitter(identifier):
    return (int(hashlib.sha256(str(identifier).encode()).hexdigest()[:8], 16) / 0xffffffff - .5) * .28


def result_figure(frame, spec):
    data, omitted = graph_rows(frame, spec)
    data['Group'] = data['Group'].fillna('').replace('', 'Ungrouped')
    groups = _groups(data, spec)
    fig = go.Figure()
    for index, group in enumerate(groups):
        sub = data[data['Group'] == group]
        color = spec.get('palette', {}).get(group, PALETTE[index % len(PALETTE)])
        fig.add_trace(go.Scatter(x=[index + _jitter(s) for s in sub.SampleId], y=sub.Value, mode='markers', name=group,
            marker={'color': color, 'size': spec.get('point_size', 9), 'line': {'width': .7, 'color': 'white'}},
            customdata=sub[['SampleId', 'Label', 'Status']].fillna('').to_numpy(),
            hovertemplate='Sample %{customdata[1]}<br>ID %{customdata[0]}<br>Value %{y:.5g}<br>%{customdata[2]}<extra></extra>'))
        mean = float(sub.Value.mean())
        error = _error(sub.Value.to_numpy(float), spec.get('error', 'None'))
        fig.add_trace(go.Scatter(x=[index-.21, index+.21], y=[mean, mean], mode='lines', line={'color':'#263536','width':2}, showlegend=False, hovertemplate='Mean %{y:.5g}<extra></extra>'))
        if error:
            fig.add_trace(go.Scatter(x=[index], y=[mean], mode='markers', marker={'size':1,'color':'#263536'}, error_y={'type':'data','array':[error],'visible':True,'color':'#263536'}, showlegend=False, hoverinfo='skip'))
    if spec.get('paired') and 'PairId' in data:
        for pair, sub in data[data.PairId.fillna('') != ''].groupby('PairId', sort=False):
            if not sub.Group.duplicated().any() and len(sub) > 1:
                ordered = sub.assign(_x=sub.Group.map({g:i for i,g in enumerate(groups)})).sort_values('_x')
                fig.add_trace(go.Scatter(x=ordered._x, y=ordered.Value, mode='lines', line={'color':'#a7b5b5','width':1}, opacity=.7, showlegend=False, hoverinfo='skip'))
    counts = data.groupby('Group').size().to_dict()
    fig.update_layout(template='plotly_white', title=spec.get('title', ''), height=int(spec.get('height', 480)),
        margin={'l':65,'r':25,'t':65,'b':65}, font={'family':'Arial, sans-serif','color':'#23383a'},
        xaxis={'tickmode':'array','tickvals':list(range(len(groups))),'ticktext':[f'{g}<br>n = {counts[g]}' for g in groups], 'title':spec.get('xlabel','Group')},
        yaxis={'title':spec.get('ylabel','Value'),'type':spec.get('scale','linear').lower()}, showlegend=False)
    return fig, data, omitted


def static_figure(frame, spec, fmt='png', caption='', _return_figure=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data, omitted = graph_rows(frame, spec)
    data['Group'] = data['Group'].fillna('').replace('', 'Ungrouped')
    groups = _groups(data, spec)
    fig, ax = plt.subplots(figsize=(int(spec.get('width', 900))/100, int(spec.get('height',480))/100))
    for index, group in enumerate(groups):
        sub = data[data.Group == group]
        ax.scatter([index+_jitter(s) for s in sub.SampleId], sub.Value, s=spec.get('point_size',9)**2/1.7, color=spec.get('palette',{}).get(group,PALETTE[index%len(PALETTE)]), edgecolors='white', linewidths=.5)
        mean = sub.Value.mean(); error = _error(sub.Value.to_numpy(float),spec.get('error','None'))
        ax.hlines(mean,index-.21,index+.21,color='#263536',linewidth=1.5)
        if error: ax.errorbar(index,mean,yerr=error,color='#263536',capsize=4,fmt='none')
    if spec.get('paired') and 'PairId' in data:
        for _, sub in data[data.PairId.fillna('')!=''].groupby('PairId',sort=False):
            if len(sub)>1 and not sub.Group.duplicated().any():
                ordered=sub.assign(_x=sub.Group.map({g:i for i,g in enumerate(groups)})).sort_values('_x')
                ax.plot(ordered._x,ordered.Value,color='#a7b5b5',linewidth=.7,zorder=0)
    ax.set_xticks(range(len(groups)),[f'{g}\nn = {int((data.Group==g).sum())}' for g in groups])
    ax.set_xlabel(spec.get('xlabel','Group')); ax.set_ylabel(spec.get('ylabel','Value')); ax.set_title(spec.get('title',''))
    ax.set_yscale(spec.get('scale','linear').lower()); ax.spines[['top','right']].set_visible(False)
    if caption: fig.text(.02,.01,caption[:240],fontsize=7,wrap=True)
    fig.tight_layout(rect=(0,.065,1,1))
    if _return_figure:
        return fig
    out=io.BytesIO(); fig.savefig(out,format=fmt,dpi=100); plt.close(fig)
    return out.getvalue()


def curve_figure(points, curve):
    usable = points[pd.to_numeric(points['Concentration'], errors='coerce').gt(0) & pd.to_numeric(points['MeanCq'], errors='coerce').notna()].copy()
    usable['_log'] = np.log10(usable['Concentration'].astype(float))
    selected = usable[usable['Include'].astype(bool)]
    excluded = usable[~usable['Include'].astype(bool)]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=selected['_log'],y=selected['MeanCq'],mode='markers',name='Selected standards',marker={'color':'#087f8c','size':9}))
    if not excluded.empty:
        fig.add_trace(go.Scatter(x=excluded['_log'],y=excluded['MeanCq'],mode='markers',name='Excluded standards',marker={'color':'#999999','symbol':'x','size':9}))
    slope, intercept = curve.get('Slope'), curve.get('Intercept')
    if len(selected) and pd.notna(slope):
        xp = np.linspace(selected['_log'].min(),selected['_log'].max(),80)
        fig.add_trace(go.Scatter(x=xp,y=slope*xp+intercept,mode='lines',name='Fitted curve',line={'color':'#253c40'}))
    fig.update_layout(template='plotly_white',height=360,xaxis_title='log₁₀(standard quantity)',yaxis_title='Mean Cq',margin={'t':20,'l':50,'r':20,'b':50})
    return fig

def analysis_report(run):
    """Printable selected-endpoint report, with all contributing/withheld samples."""
    import textwrap
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from qpcr_session import analysis_key
    spec = run['graph']
    key = 'standard_results' if spec['method']=='Standard curve' else 'ddcq_results'
    values = pd.DataFrame(run['results'][key])
    values = values[values.Gene.eq(spec['gene']) & values.Endpoint.eq(spec['endpoint'])]
    output = io.BytesIO()
    with PdfPages(output) as pdf:
        page = plt.figure(figsize=(11.7,8.3))
        page.text(.06,.92,run['name'],fontsize=20,weight='bold')
        lines = [f"EasyLab Next qPCR · revision {run['revision']}",f"Method: {spec['method']} · Assay: {spec['gene']} · Endpoint: {spec['endpoint']}",f"References: {', '.join(run['settings']['references']) or 'None selected'}",f"Calibrators: {', '.join(run['settings']['calibrators']) or 'None selected'}",f"Provisional extrapolation: {'permitted' if run['settings']['allow_extrapolation'] else 'withheld'}",f"ΔΔCq efficiency assumption accepted: {run['settings']['efficiency_confirmed']}",f"Eligible biological samples: {int(values.Eligible.sum())} of {len(values)}",f"Experiment ID: {run['runId']}",f"Analysis identifier: {analysis_key(run)[:20]}","Technical wells are averaged before biological summaries. The next page shows the selected graph. The sample table retains withheld endpoints and flags. Axis scaling does not change statistical values.","The reopenable experiment package contains original source copies, reviewed wells, standards and settings. Keep it with this report."]
        wrapped='\n\n'.join(textwrap.fill(line,110) for line in lines)
        page.text(.06,.84,wrapped,fontsize=10,va='top',linespacing=1.4)
        pdf.savefig(page);plt.close(page)
        page=static_figure(values,spec,caption=f"{run['name']} · revision {run['revision']} · {spec['endpoint']}",_return_figure=True)
        pdf.savefig(page);plt.close(page)
        report=values[['SampleId','Group','Value','Eligible','Status']].copy()
        report['Value']=report.Value.map(lambda v:f'{v:.6g}' if pd.notna(v) else '')
        report['Status']=report.Status.fillna('').map(lambda s:'\n'.join(textwrap.wrap(s,40)))
        for column,width in [('SampleId',25),('Group',18)]:
            report[column] = report[column].fillna('').map(lambda text:'\n'.join(textwrap.wrap(str(text),width)))
        start = 0
        while start < len(report):
            stop, used = start, .04
            while stop < len(report):
                height = .043+.017*max(str(v).count('\n') for v in report.iloc[stop])
                if stop > start and used+height > .86:
                    break
                used += height
                stop += 1
            part=report.iloc[start:stop]
            page,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off')
            ax.set_title(f"{spec['gene']} · {spec['endpoint']} · samples {start+1}–{start+len(part)}",loc='left',pad=20)
            table=ax.table(cellText=part.values,colLabels=['Biological sample','Group','Value','Eligible','Status'],cellLoc='left',loc='center',colWidths=[.22,.16,.12,.09,.41])
            table.auto_set_font_size(False);table.set_fontsize(8)
            for (row,col),cell in table.get_celld().items():
                cell.set_edgecolor('#d5e0e2')
                cell.set_height(.037 if row==0 else .043+.017*max(str(v).count('\n') for v in part.iloc[row-1]))
                if row==0:cell.set_facecolor('#17495c');cell.set_text_props(color='white',weight='bold')
            page.tight_layout();pdf.savefig(page);plt.close(page)
            start = stop
        saved=run.get('statistics',{})
        if saved and all(saved.get(k)==spec.get(k) for k in ['method','gene','endpoint']):
            data=pd.DataFrame(saved['table'])
            columns=[c for c in ['Group1','Group2','n1','n2','N','EffectMeanDifference','CI95Lower','CI95Upper','PAdjustedHolm'] if c in data]
            data=data[columns].map(lambda v:f'{v:.5g}' if isinstance(v,(float,np.floating)) else str(v))
            names={'Group1':'Group 1','Group2':'Group 2','n1':'n (group 1)','n2':'n (group 2)','N':'Total n','EffectMeanDifference':'Mean\ndifference','CI95Lower':'95% CI\nlower','CI95Upper':'95% CI\nupper','PAdjustedHolm':'P (Holm\nadjusted)'}
            for start in range(0,len(data),15):
                part=data.iloc[start:start+15].copy()
                for group in ['Group1','Group2']:
                    part[group]=part[group].map(lambda text:'\n'.join(textwrap.wrap(text,14)))
                page,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off');ax.set_title(f"{saved['test'].upper()} · selected comparison family",loc='left')
                table=ax.table(cellText=part.values,colLabels=[names[c] for c in columns],loc='center',cellLoc='center');table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,2.5)
                page.text(.06,.12,'P values use Holm adjustment. Effects are Group 1 minus Group 2. 95% intervals describe individual effects; they are not simultaneous intervals.',fontsize=9,wrap=True)
                page.text(.06,.08,'The design assumptions were explicitly accepted before testing. This report does not establish that the assumptions hold.',fontsize=9)
                pdf.savefig(page);plt.close(page)
    return output.getvalue()
