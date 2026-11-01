# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Read typed Bethe JSON tables without converting spectral values to floats.

This validates transport/schema, not physical convergence. Failed rows, absent
values, precision tags, root identities and all metadata remain available.
No campaign paths, preferred recovery records or channel conversions live here.
"""
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
from types import MappingProxyType


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _object(pairs):
    result = {}
    for key, value in pairs:
        _require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError(f"nonstandard JSON number: {value}")


def _freeze(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(map(_freeze, value))
    return value


def _integer(value, minimum, message):
    _require(type(value) is int and value >= minimum, message)
    return value


@dataclass(frozen=True)
class RealValue:
    """Original decimal spelling and binary working bits (or the 'exact' tag).

    Nonfinite diagnostic values are retained, but finite_decimal()/float64()
    reject them. Working precision is not an accuracy estimate.
    """
    text: str
    precision: int | str

    def decimal(self):
        _require('/' not in self.text, "exact rational: use fraction() and choose rounding explicitly")
        return Decimal(self.text)

    def fraction(self):
        """Exact rational value of the exported text; rejects nonfinite values."""
        return Fraction(self.text)

    def finite_decimal(self):
        value = self.decimal()
        _require(value.is_finite(), "nonfinite numerical value")
        return value

    def float64(self):
        value = self.fraction() if self.precision == 'exact' else self.finite_decimal()
        try:
            result = float(value)
        except OverflowError as error:
            raise ValueError("value is outside the binary64 range") from error
        _require(math.isfinite(result) and (value == 0 or result != 0),
                 "value is outside the nonzero binary64 range")
        return result


@dataclass(frozen=True)
class Table:
    title: str
    columns: tuple
    rows: tuple
    metadata: object
    summary: object
    first_row: int
    descriptor: object

    def require_metadata(self, expected):
        """Check caller-selected normalization, family and provenance fields."""
        for key, value in expected.items():
            _require(key in self.metadata and self.metadata[key] == value,
                     f"missing or different metadata: {key}")
        return self


@dataclass(frozen=True)
class TableExport:
    tables: object
    sha256: str
    source: str
    envelope: object


def _column(column):
    _require(isinstance(column, dict) and isinstance(column.get('id'), str)
             and bool(column['id']), "invalid column identifier")
    _require(type(column.get('nullable')) is bool, "missing column nullability")
    kind, encoding = column.get('type'), column.get('encoding')
    if kind == 'real':
        _require(column.get('radix') == 2, "only binary real columns are supported")
        if encoding == 'decimal_string':
            _integer(column.get('precision_bits'), 1, "invalid fixed precision")
            _require('precision' not in column, "ambiguous real precision")
        else:
            _require(encoding == 'value_precision' and column.get('precision') == 'per_value'
                     and 'precision_bits' not in column, "unsupported real encoding")
    elif kind == 'half_int':
        _require(encoding == 'number', "unsupported half-integer encoding")
        _integer(column.get('storage_bits'), 2, "invalid half-integer storage")
    elif kind in ('string', 'bool'):
        _require(encoding is None, "unsupported scalar encoding")
    else:
        _require(isinstance(kind, str) and re.fullmatch(r'u?int(8|16|32|64|128)', kind),
                 "unsupported column type")
        bits = int(re.search(r'\d+', kind)[0])
        digits = bits - (not kind.startswith('u'))
        _require(encoding == ('decimal_string' if digits > 53 else 'number'),
                 "unsupported integer encoding")


def _cell(value, column):
    if value is None:
        _require(column['nullable'], "null in nonnullable column")
        return None
    kind = column['type']
    if kind == 'real':
        if column['encoding'] == 'value_precision':
            _require(isinstance(value, dict) and set(value) == {'value', 'precision'},
                     "invalid runtime real cell")
            precision = value['precision']
            if precision != 'exact':
                _integer(precision, 1, "invalid runtime precision")
            text = value['value']
        else:
            text, precision = value, column['precision_bits']
        if precision == 'exact':
            _require(isinstance(text, str) and re.fullmatch(r'[+-]?\d+(?:/[1-9]\d*)?', text),
                     "invalid exact rational spelling")
            return RealValue(text, precision)
        _require(isinstance(text, str) and re.fullmatch(
            r'[+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|inf(?:inity)?|nan)',
            text, re.IGNORECASE), "invalid decimal real spelling")
        try:
            Decimal(text)
        except InvalidOperation as error:
            raise ValueError("invalid real value") from error
        return RealValue(text, precision)
    if kind == 'string':
        _require(isinstance(value, str), "nonstring in string column")
    elif kind == 'bool':
        _require(type(value) is bool, "nonboolean in bool column")
    elif kind == 'half_int':
        _require(type(value) in (int, Decimal), "invalid half-integer value")
        _require(not isinstance(value, Decimal) or value.is_finite(), "nonfinite half-integer")
        twice = 2 * Fraction(value)
        bits = column['storage_bits']
        _require(twice.denominator == 1 and -(1 << (bits-1)) <= twice < (1 << (bits-1)),
                 "half-integer outside storage range")
        return Decimal(value)
    else:
        if column['encoding'] == 'decimal_string':
            _require(isinstance(value, str) and re.fullmatch(r'-?(?:0|[1-9]\d*)', value),
                     "invalid decimal integer")
            value = int(value)
        _require(type(value) is int, "invalid integer value")
        bits = int(re.search(r'\d+', kind)[0])
        lo, hi = (0, 1 << bits) if kind.startswith('u') else (-(1 << (bits-1)), 1 << (bits-1))
        _require(lo <= value < hi, "integer outside storage range")
    return value


def _text_map(value, name):
    _require(isinstance(value, dict) and all(isinstance(v, str) for v in value.values()),
             f"invalid {name}")
    return _freeze(value)


def _table(value):
    _require(isinstance(value, dict) and isinstance(value.get('title'), str), "invalid table")
    metadata = _text_map(value.get('metadata'), 'metadata')
    summary = _text_map(value['summary'], 'summary') if 'summary' in value else None
    first_row = _integer(value.get('first_row'), 0, "invalid first row")
    columns, rows = value.get('columns'), value.get('rows')
    _require(isinstance(columns, list) and isinstance(rows, list), "missing schema or rows")
    for column in columns:
        _column(column)
    names = [column['id'] for column in columns]
    _require(len(set(names)) == len(names), "duplicate table column")
    decoded = []
    for row in rows:
        _require(isinstance(row, list) and len(row) == len(columns), "ragged table row")
        decoded.append(MappingProxyType({name: _cell(cell, column)
                       for name, cell, column in zip(names, row, columns)}))
    return Table(value['title'], _freeze(columns), tuple(decoded), metadata, summary,
                 first_row, _freeze({k: v for k, v in value.items() if k != 'rows'}))


def read_tables(path, *, sha256=None, table_name=None, expected_metadata=None):
    """Read a combined --format json export or a single --json-table file.

    A single-table file needs an explicit table_name. expected_metadata maps
    table names to required metadata fields; missing named tables are errors.
    The file hash binds exactly the bytes read, including incomplete results.
    The envelope 'status' describes output transport, not solver success.
    """
    path = Path(path)
    payload = path.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    _require(sha256 is None or sha256 == digest, "source checksum mismatch")
    document = json.loads(payload.decode('utf-8'), parse_float=Decimal,
                          parse_constant=_reject_constant, object_pairs_hook=_object)
    _require(isinstance(document, dict), "invalid table document")
    if 'tables' in document:
        _require(table_name is None and isinstance(document['tables'], dict), "invalid table envelope")
        raw = document['tables']
        envelope = {key: value for key, value in document.items() if key != 'tables'}
    else:
        _require(isinstance(table_name, str) and bool(table_name), "single table requires table_name")
        raw, envelope = {table_name: document}, {}
    _require(all(isinstance(name, str) and name for name in raw), "invalid table name")
    tables = {name: _table(value) for name, value in raw.items()}
    for name, expected in (expected_metadata or {}).items():
        _require(name in tables, f"missing required table: {name}")
        tables[name].require_metadata(expected)
    return TableExport(MappingProxyType(tables), digest, str(path.resolve()), _freeze(envelope))


def group_rows(table, keys):
    """Group all rows by explicit identity fields; nothing is silently replaced.

    Use (candidate_id, root_index) to index roots, or add run/precision identity
    when combining repeated calculations. No recovery precedence is assumed.
    """
    _require(not isinstance(keys, str), "identity keys must be a sequence")
    keys = tuple(keys)
    names = {column['id'] for column in table.columns}
    _require(keys and len(set(keys)) == len(keys) and set(keys) <= names, "invalid identity fields")
    groups = {}
    for row in table.rows:
        key = tuple(row[name] for name in keys)
        _require(all(value is not None and not isinstance(value, RealValue) for value in key),
                 "identity fields must be nonnull strings, integers or half-integers")
        groups.setdefault(key, []).append(row)
    return MappingProxyType({key: tuple(rows) for key, rows in groups.items()})


def index_rows(table, keys):
    """Require a unique compound identity instead of overwriting duplicate rows."""
    groups = group_rows(table, keys)
    _require(all(len(rows) == 1 for rows in groups.values()), "duplicate row identity")
    return MappingProxyType({key: rows[0] for key, rows in groups.items()})
