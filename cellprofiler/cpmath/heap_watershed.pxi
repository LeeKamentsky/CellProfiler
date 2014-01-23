"""
CellProfiler is distributed under the GNU General Public License,
but this file is licensed under the more permissive BSD license.
See the accompanying file LICENSE for details.

Copyright (c) 2003-2009 Massachusetts Institute of Technology
Copyright (c) 2009-2012 Broad Institute
All rights reserved.

Please see the AUTHORS file for credits.

Website: http://www.cellprofiler.org
"""
import numpy as np
cimport numpy as np
cimport cython

cdef struct Heapitem:
    np.int32_t value
    np.int32_t age
    np.int32_t index[0]
    
cdef struct Heapitem1:
    np.int32_t value
    np.int32_t age
    np.int32_t index[1]
    
cdef struct Heapitem4:
    np.int32_t value
    np.int32_t age
    np.int32_t index[4]

cdef inline int smaller(Heapitem *a, Heapitem *b) nogil:
    if a.value <> b.value:
      return a.value < b.value
    return a.age < b.age

include "heap_general.pxi"
