#!/usr/bin/env python
import numpy as np

def QR_decomp(A):
    Q, R = np.linalg.qr(A)
    return Q, R

def QL_decomp(A):
    Q, L = np.linalg.qr(A[::-1, ::-1])
    return Q[::-1, ::-1], L[::-1, ::-1]
