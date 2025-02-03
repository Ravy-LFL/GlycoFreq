# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def compute_contacts(np.ndarray[double, ndim=2] prot not None,
                     np.ndarray[double, ndim=2] carb not None,
                     double thr):
    """
    For each protein atom, returns the index of the first carbohydrate atom
    whose squared distance is less than thr**2, or -1 if there is no contact.
    """
    cdef int n_prot = prot.shape[0]
    cdef int n_carb = carb.shape[0]
    cdef double thr2 = thr * thr
    cdef int i, j
    cdef double dx, dy, dz, d2
    # Use np.intp_t, which corresponds to np.intp's C type
    cdef np.ndarray[np.intp_t, ndim=1] result = np.empty(n_prot, dtype=np.intp)
    cdef double[:, :] p = prot
    cdef double[:, :] c = carb

    # Initialize the result array
    for i in range(n_prot):
        result[i] = -1

    # Compute contacts while releasing the GIL for performance optimization
    with nogil:
        for i in range(n_prot):
            for j in range(n_carb):
                dx = p[i, 0] - c[j, 0]
                dy = p[i, 1] - c[j, 1]
                dz = p[i, 2] - c[j, 2]
                d2 = dx*dx + dy*dy + dz*dz
                if d2 < thr2:
                    result[i] = j
                    break
    return result

