"""Explicit reviewed source compatibility for assembled reference data.

This validates the record and exact source identities, not the scientific
sufficiency of its cited evidence. It never fetches or executes that evidence.
Without a compatibility record, combining different sources remains an error.
"""
from decimal import Decimal, InvalidOperation
import re

SCHEMA = 'bethe-spinon-source-compatibility-v1'
SOURCE_KEYS = ('Bethe revision', 'Uni20 revision')
PHYSICAL_KEYS = ('Program', 'Hamiltonian', 'Family', 'Channel', 'Exchange J', 'Density convention')


def physical_identity(metadata):
    """Require the physical identity, independently of source/numerical controls.

    Exchange is compared numerically; descriptive fields must agree exactly.
    This compares declared physics, not the correctness of a calculation.
    """
    if not isinstance(metadata, dict):
        raise ValueError('Expected physical metadata')
    identity = []
    for key in PHYSICAL_KEYS:
        value = metadata.get(key)
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f'Missing or invalid physical identity field: {key}')
        if key == 'Exchange J':
            try:
                value = Decimal(value)
            except InvalidOperation as error:
                raise ValueError('Invalid physical identity field: Exchange J') from error
            if not value.is_finite() or value <= 0:
                raise ValueError('Invalid physical identity field: Exchange J')
        identity.append(value)
    return tuple(identity)


def validate_physics(metadata):
    """Require one physical identity, independently of the source-compatibility record."""
    reference = None
    for item in metadata:
        identity = physical_identity(item)
        if reference is None:
            reference = identity
        elif identity != reference:
            fields = ', '.join(key for key, actual, expected in zip(PHYSICAL_KEYS, identity, reference)
                               if actual != expected)
            raise ValueError(f'Cannot combine different physical identities: {fields}')
    if reference is None:
        raise ValueError('No physical identities to validate')
    return reference


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
    if set(record) != {'schema', 'sources', 'sector_sources', 'reason', 'evidence'}:
        raise ValueError('Unexpected or missing cross-sector compatibility fields')
    if not isinstance(record.get('reason'), str) or not record['reason'].strip():
        raise ValueError('Source compatibility needs an explicit reason')
    sources = record.get('sources')
    if not isinstance(sources, list) or len(sources) != 2:
        raise ValueError('Source compatibility needs exactly two source identities')
    identities = [source_identity(source) for source in sources]
    if len(set(identities)) != len(identities):
        raise ValueError('Duplicate compatible source identity')
    sectors = record.get('sector_sources')
    if not isinstance(sectors, dict) or set(sectors) != {'2,4', '6'}:
        raise ValueError('Compatibility needs explicit sector-scoped sources')
    for allowed in sectors.values():
        if not isinstance(allowed, list) or len(allowed) != 1:
            raise ValueError('Each sector needs one exact source identity')
        pairs = [source_identity(source) for source in allowed]
        if len(set(pairs)) != len(pairs) or not set(pairs) <= set(identities):
            raise ValueError('Invalid sector-scoped source identities')
    declared = {source_identity(source) for allowed in sectors.values() for source in allowed}
    if declared != set(identities):
        raise ValueError('Compatibility sources must match the two sector declarations')
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

A record permits one two/four source and one six source. Every declared source
must occur in the actual data; it cannot approve an uncomputed future cohort.
"""
    metadata = list(metadata)
    observed = frozenset(source_identity(item) for item in metadata)
    if not observed:
        raise ValueError('No source identities to validate')
    if compatibility is None:
        if len(observed) != 1:
            raise ValueError('Cannot combine different source revisions without reviewed compatibility')
    else:
        if observed != compatibility_sources(compatibility):
            raise ValueError('Actual source revisions differ from the explicit compatibility record')
        for item in metadata:
            allowed = compatibility['sector_sources'].get(item.get('Selected spinon sectors'), [])
            if source_identity(item) not in {source_identity(source) for source in allowed}:
                raise ValueError('Source revision is not approved for these sectors')
    return observed
