"""Local Next storage and versioned, non-executing experiment packages."""
from __future__ import annotations
import hashlib
import json
import os
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path

FORMAT = 'easylab.qpcr'
SCHEMA_VERSION = 1
MAX_PACKAGE_BYTES = 50 * 1024 * 1024


def runtime_root():
    configured = os.environ.get('EASYLAB_RUNTIME_ROOT')
    source_root = Path(__file__).resolve().parent.parents[1]
    default = source_root / 'runtime' if source_root.name == 'easylab-next' else Path.home() / 'Documents' / 'EasyLab Next' / 'runtime'
    root = Path(configured or default).expanduser().absolute()
    stable = [Path.home() / '.easylab', Path.home() / 'Documents' / 'Easylab', Path.home() / 'Library' / 'Application Support' / 'easylab-suite']
    resolved = root.resolve()
    names = {part.casefold() for part in resolved.parts}
    if not names.intersection({'next','easylab-next','easylab next','easylab suite next dev'}) or names.intersection({'easylab-suite','easylab-modules'}):
        raise ValueError('Use a separate folder named EasyLab Next for the runtime; existing Suite workspaces are protected.')
    if any(resolved == p.resolve() or p.resolve() in resolved.parents for p in stable):
        raise ValueError('EasyLab Next cannot use the current app’s storage folder.')
    if root != resolved:
        raise ValueError('Use a real Next runtime folder, not a symbolic link.')
    return resolved


def storage_paths(create=True):
    root = runtime_root()
    defaults = {'data': root / 'data' / 'qPCR Analysis' / 'data', 'settings': root / 'settings' / 'qpcr-analysis', 'exports': root / 'data' / 'qPCR Analysis' / 'exports'}
    envs = {'data': 'EASYLAB_DATA_PATH', 'settings': 'EASYLAB_SETTINGS_PATH', 'exports': 'EASYLAB_EXPORTS_PATH'}
    paths = {}
    for name, fallback in defaults.items():
        value = Path(os.environ.get(envs[name], str(fallback))).expanduser().absolute()
        resolved = value.resolve()
        if root not in resolved.parents or resolved != value:
            raise ValueError(f'{name.title()} must be inside the separate Next runtime folder.')
        paths[name] = resolved
        if create:
            resolved.mkdir(parents=True, exist_ok=True)
    for nested in [paths['data'] / 'experiments', paths['settings'] / 'presets.json']:
        if nested.is_symlink() or nested.absolute() != nested.resolve():
            raise ValueError('Linked experiment or preset storage is not allowed in EasyLab Next.')
    return paths


def json_safe(value):
    # Dataframes are converted by callers with JSON orient=records, preserving nulls.
    return json.loads(json.dumps(value, ensure_ascii=False, allow_nan=False))


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def new_run():
    return {'engineVersion': 'next-qpcr-1', 'runId': str(uuid.uuid4()), 'name': 'Untitled experiment', 'revision': 0, 'createdAt': datetime.now(timezone.utc).isoformat(), 'sources': [], 'raw': [], 'wells': [], 'standards': [], 'settings': {'preset': 'Meghamsh — reviewed qPCR', 'top': 16384.0, 'factor': 4.0, 'threshold': 0.75, 'efficiency_min':90., 'efficiency_max':110., 'efficiency_difference':5., 'references': [], 'calibrators': [], 'allow_extrapolation': False, 'efficiency_confirmed': False}, 'approvals': {}, 'results': {}, 'graph': {}, 'statistics': {}}


def encode_package(run):
    clean = json_safe(run)
    payload = {'format': FORMAT, 'schemaVersion': SCHEMA_VERSION, 'run': clean, 'sha256': fingerprint(clean)}
    data = json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False).encode('utf-8')
    if len(data) > MAX_PACKAGE_BYTES:
        raise ValueError('This experiment exceeds the 50 MB package limit. Split the input into smaller experiments; the previous experiment is preserved.')
    return data


def decode_package(data):
    if len(data) > MAX_PACKAGE_BYTES:
        raise ValueError('The experiment package exceeds 50 MB.')
    try:
        payload = json.loads(data)
    except (ValueError, UnicodeError) as exc:
        raise ValueError('This file is not a readable experiment package.') from exc
    if not isinstance(payload, dict) or payload.get('format') != FORMAT or payload.get('schemaVersion') != SCHEMA_VERSION:
        raise ValueError('Unsupported experiment format or version; the existing experiment was preserved.')
    run = payload.get('run')
    if not isinstance(run, dict) or not all(k in run for k in new_run()):
        raise ValueError('The experiment package is incomplete.')
    try:
        uuid.UUID(run['runId'])
        valid = fingerprint(run) == payload.get('sha256')
    except (ValueError, TypeError, KeyError) as exc:
        raise ValueError('The experiment package has invalid identifiers or values.') from exc
    if not valid:
        raise ValueError('The package contents do not match their saved checksum.')
    for key in ['sources', 'raw', 'wells', 'standards']:
        if not isinstance(run[key], list) or not all(isinstance(row, dict) for row in run[key]):
            raise ValueError(f'Invalid {key} table.')
    for key in ['settings', 'approvals', 'results', 'graph', 'statistics']:
        if not isinstance(run[key], dict):
            raise ValueError(f'Invalid {key} settings.')
    if type(run['revision']) is not int or run['revision'] < 0:
        raise ValueError('Invalid saved revision.')
    defaults = new_run()['settings']
    if any(k not in run['settings'] for k in defaults):
        raise ValueError('Required method settings are missing.')
    if not isinstance(run['name'], str):
        raise ValueError('Invalid experiment name.')
    for field in ['references', 'calibrators']:
        if not isinstance(run['settings'][field], list) or not all(isinstance(x,str) for x in run['settings'][field]):
            raise ValueError('Invalid method selection.')
    for field in ['top', 'factor', 'threshold']:
        if not isinstance(run['settings'][field], (int,float)) or not __import__('math').isfinite(run['settings'][field]):
            raise ValueError('Invalid numerical method setting.')
    if run['settings']['top'] <= 0 or run['settings']['factor'] <= 1 or run['settings']['threshold'] < 0:
        raise ValueError('Invalid standard series or warning threshold.')
    if any(k not in [str(i) for i in range(7)] for k in run['approvals']):
        raise ValueError('Invalid review step.')
    return run


def save_run(run):
    data = encode_package(run)
    directory = storage_paths()['data'] / 'experiments'
    directory.mkdir(exist_ok=True)
    # UUID and numeric revision avoid filenames supplied by imported data.
    identifier = str(uuid.UUID(run['runId']))
    if not isinstance(run['revision'], int) or run['revision'] < 0:
        raise ValueError('Invalid revision.')
    stem = f"{identifier}-r{run['revision']}"
    target = directory / f'{stem}.easylab.json'
    count = 1
    while target.exists():
        if target.is_symlink():
            raise ValueError('Linked experiment files cannot be reused.')
        if target.read_bytes() == data:
            return target
        target = directory / f'{stem}-{count}.easylab.json'
        count += 1
    # Exclusive creation preserves all previous saved revisions.
    with target.open('xb') as f:
        f.write(data)
    return target


def load_presets():
    path = storage_paths()['settings'] / 'presets.json'
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    return data if isinstance(data, dict) else {}


def save_preset(name, settings, graph):
    if not name.strip():
        raise ValueError('Name the preset first.')
    presets = load_presets()
    # Never carry experimental identities, review decisions or references into a new run.
    presets[name.strip()] = {'settings': {k: settings[k] for k in ['top', 'factor', 'threshold', 'efficiency_min', 'efficiency_max', 'efficiency_difference']}, 'graph': {k: v for k, v in graph.items() if k in ['scale', 'error', 'point_size', 'width', 'height', 'palette']}}
    target = storage_paths()['settings'] / 'presets.json'
    fd, tmp = tempfile.mkstemp(dir=str(target.parent), suffix='.tmp')
    with os.fdopen(fd, 'w') as f:
        json.dump(presets, f, indent=2)
    os.replace(tmp, target)
