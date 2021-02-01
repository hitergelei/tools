#!/usr/bin/env python
import numpy as np

def run_KDE(
    x, y,
    nbin=100,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
    ):
    """
    Imported from Phonopy by Togo.
    """
    #
    if not xmin:
        xmin = np.min(x) *0.9
    if not xmax:
        xmax = np.max(x) *1.1
    if not ymin:
        ymin = np.min(y) *0.9
    if not ymax:
        ymax = np.max(y) *1.1

    #
    values = np.vstack([x.ravel(), y.ravel()])
    from scipy import stats
    kernel = stats.gaussian_kde(values)

    #
    xi, yi = np.mgrid[xmin:xmax:nbin*1j, ymin:ymax:nbin*1j]
    positions = np.vstack([xi.ravel(), yi.ravel()])
    zi = np.reshape(kernel(positions).T, xi.shape)
    return (xi, yi, zi), nbin, (xmin, xmax, ymin, ymax)
