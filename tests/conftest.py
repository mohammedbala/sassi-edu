"""pytest configuration shared by all SASSI-EDU tests.

Limit the threads of the BLAS/LAPACK libraries before numpy is imported: the verification problems are
small enough that multithreaded BLAS brings little, and several pytest processes running at once
(parallel developers or CI jobs) otherwise oversubscribe the CPU.  Values already set in the
environment are respected.
"""
import os

for _var in ("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "2")
