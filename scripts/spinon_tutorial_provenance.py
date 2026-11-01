"""Explicit reviewed source compatibility for assembled reference data.

This validates the record and exact source identities, not the scientific
sufficiency of its cited evidence. It never fetches or executes that evidence.
Without a compatibility record, combining different sources remains an error.
"""
import re

SCHEMA = 'bethe-spinon-source-compatibility-v1'
SOURCE_KEYS = ('Bethe revision', 'Uni20 revision')


def source_identity(metadata):
    if not isinstance(metadata, dict):
        raise ValueError('Expected source metadata')
    identity = tuple(metadata.get(key) for key in SOURCE_KEYS)
    if any(not isinstance(value, str) or not re.fullmatch('[0-9a-f]{40}', value)
           for value in identity):
        raise ValueError('Exact clean Bethe and Uni20 revisions are required')
    return identity


def compatibility_sources(record):
    if not isinstance(record, dict) or record.get('schema') != SCHEMA:
        raise ValueError('Unknown source-compatibility record')
    if not isinstance(record.get('reason'), str) or not record['reason'].strip():
        raise ValueError('Source compatibility needs an explicit reason')
    sources = record.get('sources')
    if not isinstance(sources, list) or len(sources) < 2:
        raise ValueError('Source compatibility needs at least two exact source identities')
    identities = [source_identity(source) for source in sources]
    if len(set(identities)) != len(identities):
        raise ValueError('Duplicate compatible source identity')
    sectors = record.get('sector_sources')
    if not isinstance(sectors, dict) or not sectors or set(sectors) - {'2,4', '6'}:
        raise ValueError('Compatibility needs explicit sector-scoped sources')
    for allowed in sectors.values():
        if not isinstance(allowed, list) or not allowed:
            raise ValueError('Each sector needs its allowed source identities')
        pairs = [source_identity(source) for source in allowed]
        if len(set(pairs)) != len(pairs) or not set(pairs) <= set(identities):
            raise ValueError('Invalid sector-scoped source identities')
    evidence = record.get('evidence')
    if not isinstance(evidence, list) or not evidence:
        raise ValueError('Source compatibility needs reviewed evidence')
    hashes = set()
    for item in evidence:
        if (not isinstance(item, dict) or not isinstance(item.get('description'), str)
                or not item['description'].strip()
                or not isinstance(item.get('sha256'), str)
                or not re.fullmatch('[0-9a-f]{64}', item['sha256'])):
            raise ValueError('Compatibility evidence needs a description and SHA256')
        if item['sha256'] in hashes:
            raise ValueError('Duplicate compatibility evidence')
        hashes.add(item['sha256'])
    return frozenset(identities)


def validate_sources(metadata, compatibility=None):
    """Check exact source pairs; callers must also check physical conventions.

An explicit record may cover a declared cohort with no available points yet.
The collection reader separately rejects unused identities in that record.
"""
    metadata = list(metadata)
    observed = frozenset(source_identity(item) for item in metadata)
    if not observed:
        raise ValueError('No source identities to validate')
    if compatibility is None:
        if len(observed) != 1:
            raise ValueError('Cannot combine different source revisions without reviewed compatibility')
    else:
        if not observed <= compatibility_sources(compatibility):
            raise ValueError('Source revision is outside the explicit compatibility record')
        for item in metadata:
            allowed = compatibility['sector_sources'].get(item.get('Selected spinon sectors'), [])
            if source_identity(item) not in {source_identity(source) for source in allowed}:
                raise ValueError('Source revision is not approved for these sectors')
    return observed
