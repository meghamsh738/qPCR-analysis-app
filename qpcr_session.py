"""Explicit experiment edits and dependent review invalidation."""
from copy import deepcopy
from qpcr_storage import fingerprint, new_run

STEPS = ['Import', 'Map samples', 'Review wells', 'Methods', 'Results', 'Graphs & statistics', 'Save & export']

def analysis_key(run):
    return fingerprint({k: run[k] for k in ['wells', 'standards', 'settings']})

def changed(run, updates, stage):
    result = deepcopy(run)
    result.update(deepcopy(updates))
    result['revision'] += 1
    for key in list(result['approvals']):
        if int(key) >= stage:
            del result['approvals'][key]
    if stage < 4:
        result['results'] = {}
        result['statistics'] = {}
    result['approvals'][str(stage)] = result['revision']
    return result

def fresh_run(previous):
    result = new_run()
    for key in ['preset', 'top', 'factor', 'threshold', 'efficiency_min', 'efficiency_max', 'efficiency_difference']:
        result['settings'][key] = previous['settings'][key]
    result['graph'] = {k: v for k, v in previous['graph'].items() if k in ['scale', 'error', 'point_size', 'width', 'height', 'palette']}
    return result

def graph_defaults(saved, method, gene, endpoint, groups):
    spec = deepcopy(saved)
    if not spec.get('title') or spec.get('title') == spec.get('gene'):
        spec['title'] = gene
    if not spec.get('ylabel') or spec.get('ylabel') == spec.get('endpoint'):
        spec['ylabel'] = endpoint
    spec.update(method=method,gene=gene,endpoint=endpoint)
    for key,value in {'scale':'Linear','error':'None','point_size':8,'width':1000,'height':600,'paired':False,'palette':{}}.items():
        spec.setdefault(key,value)
    if set(spec.get('group_order',[])) != set(groups):
        spec['group_order'] = list(groups)
    return spec

def preserve_standard_decisions(generated, existing, preserve_concentrations=True):
    """Retain explicit selections when a sample map is revisited or regrouped."""
    by_key = {}
    for row in existing:
        key = tuple(row.get(k) for k in ['Plate','Gene','Label'])
        by_key.setdefault(key,[]).append(row)
    output=[]
    for row in generated:
        key=tuple(row.get(k) for k in ['Plate','Gene','Label'])
        output.extend([dict(row,**{k:v for k,v in saved.items() if preserve_concentrations or k != 'Concentration'}) for saved in by_key.get(key,[])] or [row])
    return output
