"""Small fictional practice experiment, never a user's research data."""
import pandas as pd

def practice_data():
    rows = []
    n = 0
    def add(gene, label, kind, cq, group='', sample='', pair=''):
        nonlocal n
        n += 1
        rows.append({'Plate':'Practice plate','Well':f'{chr(65+(n-1)//24)}{(n-1)%24+1}', 'Gene':gene,'Type':kind,'Label':label,'Cq':str(cq),'Group':group,'SampleId':sample,'PairId':pair})
    for gene, offset in [('Target', 0), ('Reference', 1)]:
        for level in range(8):
            for delta in [-.05, .05]:
                add(gene, f'Std {level+1}', 'Standard', 14+2*level+offset+delta)
        for i in range(6):
            group = 'Control' if i < 3 else 'Treatment'
            value = 21 + (i%3)*.35 + offset - ((1.2 if i >= 3 else 0) + (i%3)*.12 if gene == 'Target' else 0)
            for delta in [-.07,.07]:
                add(gene,f'Sample {i+1}','Sample',value+delta,group,f'Practice-{i+1}')
        add(gene,'NTC','NTC','Undetermined')
    return pd.DataFrame(rows)
