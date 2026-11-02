#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Fresh analytic neighbour/cutoff calculations; optionally inspect an export."""
import argparse
import json
import numpy as np
from numpy.polynomial.legendre import leggauss

from spectral_spatial import cumulative_weights, spatial_weights
from spectral_tables import read_tables
from spectral_time_transform import statistics


def example():
    x, w = leggauss(24)
    q, measure = np.pi*(x+1)/2, w/2
    # m0(q)=1/4-0.10*cos(q)+0.04*cos(2q), a positive test spectrum.
    static = .25-.10*np.cos(q)+.04*np.cos(2*q)
    operator = spatial_weights(q, [0, 1, 2], measure, even=True)
    neighbours = operator@static
    np.testing.assert_allclose(neighbours, [.25, -.05, .02], atol=4e-16)
    # Four independent replicate vectors. Each vector contains correlated bins.
    replicates = np.array([[1., 2., 3., 4.], [2., 1., 5., 0.],
                           [0., 3., 4., 2.], [3., 0., 2., 3.]])
    mean = replicates.mean(axis=0)
    covariance = np.cov(replicates, rowvar=False, ddof=1)/len(replicates)
    cuts = [.1, .5, 1.]
    operator = cumulative_weights([.1, .4, .7, 1.], cuts, [1/4]*4)
    cumulative = statistics(operator, mean, covariance)
    shell = statistics(operator[[2]]-operator[[1]], mean, covariance)
    return {'analytic_static_neighbours': dict(zip(('C(0)', 'C(1)', 'C(2)'), neighbours.tolist())),
            'cutoffs': cuts, 'cumulative': cumulative, 'outer_minus_inner': shell,
            'uncertainty': 'sampling only; four independent synthetic replicates; no grid or tail bound'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', help='Optional existing JSON table export to validate and inspect')
    parser.add_argument('--table-name', help='Required only for a single-table JSON file')
    parser.add_argument('--sha256', help='Optional expected hash of the exact file bytes')
    args = parser.parse_args()
    if (args.table_name or args.sha256) and not args.tables:
        parser.error('--table-name/--sha256 require --tables')
    result = example()
    if args.tables:
        export = read_tables(args.tables, table_name=args.table_name, sha256=args.sha256)
        result['export'] = {'source': export.source, 'sha256': export.sha256,
                            'tables': {name: {'rows': len(table.rows), 'first_row': table.first_row,
                                             'metadata': dict(table.metadata),
                                             'summary': None if table.summary is None else dict(table.summary)}
                                       for name, table in export.tables.items()}}
    print(json.dumps(result, indent=2, allow_nan=False))
