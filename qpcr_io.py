"""Preview-first import: no inferred numeric values and no discarded data rows."""
from __future__ import annotations
import base64
import csv
import io
import json
from pathlib import Path
import pandas as pd

ALIASES = {
    'Plate': ['plate', 'plate id', 'plate name'], 'Well': ['well', 'well position', 'well id'],
    'Gene': ['gene', 'assay', 'target name', 'detector name', 'target'],
    'Label': ['label', 'sample name', 'sample', 'sample id'],
    'Cq': ['cq', 'ct', 'c[t]', 'ct value', 'cq value'],
    'Type': ['type', 'task', 'sample type'], 'Group': ['group', 'treatment', 'condition'],
    'SampleId': ['biological sample id', 'biologicalsampleid', 'sampleid'],
    'PairId': ['pairid', 'pair id', 'subject id'],
    'Tissue': ['tissue','tissue type','region'], 'Batch': ['batch','batch id','cdna batch'],
}


def records(frame):
    return json.loads(frame.to_json(orient='records', date_format='iso'))


def suggested_mapping(frame):
    names = {str(c).strip().casefold(): str(c) for c in frame.columns}
    return {target: next((names[a] for a in aliases if a in names), '') for target, aliases in ALIASES.items()}


def _header_index(rows):
    for index, row in enumerate(rows[:100]):
        tokens = {str(x).strip().casefold() for x in row if x is not None}
        if tokens.intersection(ALIASES['Cq']) and (tokens.intersection(ALIASES['Gene']) or tokens.intersection(ALIASES['Label'])):
            return index
    return 0


def read_table(data, filename):
    if len(data) > 50 * 1024 * 1024:
        raise ValueError('Please import a file smaller than 50 MB.')
    if filename.lower().endswith('.xlsx'):
        book = pd.ExcelFile(io.BytesIO(data), engine='openpyxl')
        props = getattr(book.book,'custom_doc_props',[])
        if any(prop.name == 'EasyLabWorkbookKind' and prop.value == 'qpcr-review' for prop in props):
            raise ValueError('This is an editable review workbook, not a well-table export. Open its experiment package to continue in the app. To bring spreadsheet edits back, export a reviewed flat well table with Well, Gene, Label and Cq columns and confirm its mappings again.')
        sheet = next((s for s in book.sheet_names if s.strip().casefold() == 'results'), book.sheet_names[0])
        preview = pd.read_excel(book, sheet_name=sheet, header=None, dtype=str, keep_default_na=False)
        offset = _header_index(preview.values.tolist())
        header = [str(x).strip() for x in preview.iloc[offset].tolist()]
        if len(header) != len(set(header)):
            raise ValueError('Duplicate column headings need correction.')
        frame = pd.read_excel(book, sheet_name=sheet, header=offset, dtype=str, keep_default_na=False)
        details = f'{sheet} · header row {offset + 1}'
    else:
        text = data.decode('utf-8-sig')
        lines = text.splitlines()
        candidate = next((line for line in lines if any(s in line.casefold() for s in ['\tct\t', '\tcq\t', ',cq', ',ct'])), '\n'.join(lines[:20]))
        counts = {delimiter: candidate.count(delimiter) for delimiter in ['\t', ',', ';', '|']}
        separator = max(counts, key=counts.get)
        if counts[separator] == 0:
            raise ValueError('Use a CSV or tab-separated table with a header row.')
        rows = list(csv.reader(io.StringIO(text), delimiter=separator))
        offset = _header_index(rows)
        header = [x.strip() for x in rows[offset]]
        if len(header) != len(set(header)):
            raise ValueError('Duplicate column headings need correction.')
        frame = pd.read_csv(io.StringIO(text), sep=separator, header=offset, dtype=str, keep_default_na=False, skip_blank_lines=False)
        details = f'header row {offset + 1}'
    frame.columns = [str(c).strip() for c in frame.columns]
    if frame.empty or frame.columns.duplicated().any():
        raise ValueError('The table is empty or has duplicate column headings.')
    frame = frame.fillna('')
    # Completely blank lines are retained and flagged during mapping, never guessed as samples.
    return frame, details


def map_columns(frame, mapping, plate_name='Plate 1'):
    for required in ['Well', 'Gene', 'Label', 'Cq']:
        if not mapping.get(required):
            raise ValueError(f'Select the {required} column before continuing.')
    selected = [v for v in mapping.values() if v]
    if len(selected) != len(set(selected)):
        raise ValueError('Each input column must have one role; review the repeated selections.')
    out = frame.copy()
    for target, column in mapping.items():
        out[target] = frame[column] if column else (plate_name if target == 'Plate' else '')
    return out


def source_record(filename, data, details=''):
    return {'name': Path(filename).name, 'base64': base64.b64encode(data).decode('ascii'), 'description': details}


def excel_results(tables, settings):
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        for name, frame in tables.items():
            frame.to_excel(writer, sheet_name=name[:31], index=False)
        pd.DataFrame([{'Setting': k, 'Value': json.dumps(v, ensure_ascii=False)} for k, v in settings.items()]).to_excel(writer, sheet_name='Analysis settings', index=False)
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = 'A2'
            sheet.auto_filter.ref = sheet.dimensions
            for column in sheet.columns:
                sheet.column_dimensions[column[0].column_letter].width = min(42, max(14, max(len(str(cell.value or '')) for cell in column[:100]) + 2))
            # IDs/text supplied by input are literal spreadsheet text, never executable formulas.
            for row in sheet.iter_rows():
                for cell in row:
                    if cell.data_type == 'f':
                        cell.data_type = 's'
    return output.getvalue()
