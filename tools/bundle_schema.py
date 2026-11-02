"""Fail-closed validation for the public input map and release manifest."""
from pathlib import Path, PurePosixPath
import re

HASH=re.compile(r'[0-9a-f]{64}\Z')

def relative_path(value):
    if not isinstance(value,str) or not value or '\\' in value or ':' in value or any(ord(c)<32 or ord(c)==127 for c in value):
        raise ValueError('Invalid relative path')
    q=PurePosixPath(value)
    if q.is_absolute() or value in {'.','..'} or any(x in {'.','..',''} for x in value.split('/')):
        raise ValueError('Unsafe relative path: '+value)
    return q

def safe(base, relative):
    relative_path(relative)
    base=Path(base).resolve();target=(base/relative).resolve()
    if not target.is_relative_to(base):raise ValueError('Path escapes its root')
    return target

def identity(row):
    if not isinstance(row.get('sha256'),str) or not HASH.fullmatch(row['sha256']):
        raise ValueError('Invalid SHA-256')
    if type(row.get('bytes')) is not int or row['bytes']<0:raise ValueError('Invalid byte count')

def validate(table, manifest):
    if not isinstance(table.get('files'),list) or not isinstance(manifest.get('assets'),list):
        raise ValueError('Invalid manifest lists')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+',manifest.get('repository','')):
        raise ValueError('Invalid repository identifier')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+',manifest.get('tag','')):
        raise ValueError('Invalid release tag')
    assets={};folded=set()
    for row in manifest['assets']:
        identity(row);name=row.get('name');q=relative_path(name)
        if len(q.parts)!=1 or name.casefold() in folded:raise ValueError('Duplicate or nonlocal asset name')
        folded.add(name.casefold());assets[name]=row
    paths=set();lengths={};used=set()
    for row in table['files']:
        identity(row);name=row.get('path');relative_path(name)
        if name.casefold() in paths:raise ValueError('Duplicate logical path: '+name)
        paths.add(name.casefold());h=row['sha256']
        if h in lengths and lengths[h]!=row['bytes']:raise ValueError('Alias length disagreement')
        lengths[h]=row['bytes']
        stored=row.get('repository_path')
        if stored is not None:
            relative_path(stored)
            if not stored.startswith(('source/','data/objects/')):
                raise ValueError('Unexpected repository object location')
        else:
            asset=row.get('release_asset')
            if asset not in assets:raise ValueError('Undeclared release asset')
            if row.get('member')!='objects/'+h:raise ValueError('Unexpected ZIP member')
            used.add(asset)
    if used!=set(assets):raise ValueError('Unused or missing release asset declaration')
    return assets
