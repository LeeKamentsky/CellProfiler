"""watershed.pyx - scithon implementation of guts of watershed

CellProfiler is distributed under the GNU General Public License,
but this file is licensed under the more permissive BSD license.
See the accompanying file LICENSE for details.

Copyright (c) 2003-2009 Massachusetts Institute of Technology
Copyright (c) 2009-2013 Broad Institute
All rights reserved.

Please see the AUTHORS file for credits.

Website: http://www.cellprofiler.org
"""

import numpy as np
cimport numpy as np
cimport cython

cdef extern from "numpy/arrayobject.h":
    ctypedef class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
    cdef void import_array()
    cdef int  PyArray_ITEMSIZE(np.ndarray) nogil
    cdef void * PyArray_DATA(np.ndarray) nogil
    cdef np.npy_intp * PyArray_STRIDES(np.ndarray) nogil
    cdef np.npy_intp PyArray_STRIDE(np.ndarray, int) nogil
    cdef void * PyArray_BYTES(np.ndarray) nogil
    cdef np.npy_intp PyArray_NDIM(np.ndarray) nogil
    cdef np.npy_intp *PyArray_SHAPE(np.ndarray) nogil
import_array()

DTYPE_INT32 = np.int32
ctypedef np.int32_t DTYPE_INT32_t
DTYPE_BOOL = np.bool
ctypedef np.int8_t DTYPE_BOOL_t

include "heap_watershed.pxi"

cdef inline DTYPE_INT32_t GETELEMENT2(np.ndarray obj, int i, int j) nogil:
    return (<DTYPE_INT32_t *>(PyArray_BYTES(obj) + PyArray_STRIDE(obj, 0)*i + \
            PyArray_STRIDE(obj, 1)*j))[0]

cdef inline DTYPE_INT32_t GETELEMENT1(np.ndarray obj, int i) nogil:
    return (<DTYPE_INT32_t *>(PyArray_BYTES(obj) + PyArray_STRIDE(obj, 0)*i))[0]
    
cdef inline DTYPE_INT32_t GETELEMENT(np.ndarray obj, int i) nogil:
    return (<DTYPE_INT32_t *>PyArray_BYTES(obj))[i]
    
cdef inline DTYPE_BOOL_t GETBOOL1(np.ndarray obj, int i) nogil:
    return (<DTYPE_BOOL_t*>(PyArray_BYTES(obj) + PyArray_STRIDE(obj, 0)*i))[0]
    
cdef inline DTYPE_BOOL_t GETBOOL2(np.ndarray obj, int i, int j) nogil:
    return (<DTYPE_BOOL_t *>(PyArray_BYTES(obj) + PyArray_STRIDE(obj, 0)*i +
                            PyArray_STRIDE(obj, 1)*j))[0]
cdef inline DTYPE_BOOL_t GETBOOL(np.ndarray obj, int i) nogil:
    return (<DTYPE_BOOL_t *>PyArray_BYTES(obj))[i]
cdef inline void PUTELEMENT1(np.ndarray obj, int i, DTYPE_INT32_t value) nogil:
    (<DTYPE_INT32_t *>(PyArray_BYTES(obj) + PyArray_STRIDE(obj, 0)*i))[0] = value
cdef inline void PUTELEMENT(np.ndarray obj, int i, DTYPE_INT32_t value) nogil:
    (<DTYPE_INT32_t *>PyArray_BYTES(obj))[i] = value
    
cdef inline void ELEM2LD(
    np.ndarray lower_distance, Heapitem4 *elem, int ld_level_stride) nogil:
    #
    # Put the lower propagation distance into index 0 of the lower_distance
    # array and the level into index 1.
    #
    # lower_distance - the M+1 dimension lower distance
    #
    # elem - the heap item whose age and level is being updated
    #
    cdef:
        int idx_age = elem.index[TW_LD_IDX]
    PUTELEMENT(lower_distance, idx_age, elem.age)
    PUTELEMENT(lower_distance, idx_age + ld_level_stride, elem.value)
    
cdef inline void PUTAGE(np.ndarray lower_distance, Heapitem4 *elem) nogil:
    #
    # Only store the age in the lower distance, having previously figured
    # out that the target level is the same as in the LD
    #
    cdef:
        int idx_age = elem.index[TW_LD_IDX]
    PUTELEMENT(lower_distance, idx_age, elem.age)
    
cdef inline void LD2ELEM(
    np.ndarray lower_distance, Heapitem4 *elem, int ld_level_stride) nogil:
    cdef:
        int idx_age = elem[0].index[TW_LD_IDX]
    elem[0].age = GETELEMENT(lower_distance, idx_age)
    elem[0].value = GETELEMENT(lower_distance, idx_age + ld_level_stride)
    
@cython.boundscheck(False)
def watershed(np.ndarray[DTYPE_INT32_t,ndim=1,negative_indices=False, mode='c'] image,
              np.ndarray[DTYPE_INT32_t,ndim=2,negative_indices=False, mode='c'] pq,
              DTYPE_INT32_t age,
              np.ndarray[DTYPE_INT32_t,ndim=2,negative_indices=False, mode='c'] structure,
              DTYPE_INT32_t ndim,
              np.ndarray[DTYPE_BOOL_t,ndim=1,negative_indices=False, mode='c'] mask,
              np.ndarray[DTYPE_INT32_t,ndim=1,negative_indices=False, mode='c'] image_shape,
              np.ndarray[DTYPE_INT32_t,ndim=1,negative_indices=False, mode='c'] output):
    """Do heavy lifting of watershed algorithm
    
    image - the flattened image pixels, converted to rank-order
    pq    - the priority queue, starts with the marked pixels
            the first element in each row is the image intensity
            the second element is the age at entry into the queue
            the third element is the index into the flattened image or labels
            the remaining elements are the coordinates of the point
    age   - the next age to assign to a pixel
    structure - a numpy int32 array containing the structuring elements
                that define nearest neighbors. For each row, the first
                element is the stride from the point to its neighbor
                in a flattened array. The remaining elements are the
                offsets from the point to its neighbor in the various
                dimensions
    ndim  - # of dimensions in the image
    mask  - numpy boolean (char) array indicating which pixels to consider
            and which to ignore. Also flattened.
    image_shape - the dimensions of the image, for boundary checking,
                  a numpy array of np.int32
    output - put the image labels in here
    """
    cdef Heapitem elem
    cdef Heapitem new_elem
    cdef DTYPE_INT32_t nneighbors = structure.shape[0] 
    cdef DTYPE_INT32_t i = 0
    cdef DTYPE_INT32_t index = 0
    cdef DTYPE_INT32_t old_index = 0
    cdef DTYPE_INT32_t max_index = image.shape[0]

    cdef Heap *hp = <Heap *> heap_from_numpy2(1)

    for i in range(pq.shape[0]):
        elem.value = pq[i, 0]
        elem.age = pq[i, 1]
        elem.index[0] = pq[i, 2]
        heappush(hp, <Heapitem *>&elem)

    while hp.items > 0:
        #
        # Pop off an item to work on
        #
        heappop(hp, <Heapitem *>&elem)
        ####################################################
        # loop through each of the structuring elements
        #
        old_index = elem.index[0]
        for i in range(nneighbors):
            # get the flattened address of the neighbor
            index = structure[i,0]+old_index
            if index < 0 or index >= max_index or output[index] or not mask[index]:
                continue

            new_elem.value   = image[index]
            new_elem.age   = elem.age + 1
            new_elem.index[0]   = index
            age          += 1
            output[index] = output[old_index]
            #
            # Push the neighbor onto the heap to work on it later
            #
            heappush(hp, <Heapitem *>&new_elem)
    heap_done(hp)

cdef enum:
    TW_IMAGE_IDX = 0
    TW_MASK_IDX = 1
    TW_OUTPUT_IDX = 2
    TW_LD_IDX = 3
    TW_IDX_COUNT = 4
    
########################################################################
#
# Tiled watershed is based on a parallel watershed algorithm:
#    Moga et al, "Parallel watershed transformation algorithms for 
#    image segmentation", Parallel Computing, 1998, pp 1981-2001
#
# The same qualities that allow a parallel execution allow us to
# break the image into tiles / cubes / hypercubes and process the
# cubes locally.
#
# In order to work across tiles, the seeded watershed has to maintain the
# following for each pixel:
#
# * The lowest level of the neighboring labeling pixel. The seed pixels'
#   level is the image value. A labeling pixel's level propagates to the 
#   labeled pixel unless the labeled pixel's image value is higher in which
#   case the labeled pixel's level is its image value. A neigboring labeled
#   pixel can propagate into a pixel of higher lowest level.
# * The age of the wavefront when propagating from pixels of the same level.
#   A pixel's age is zero if its labeling pixel is of different level than
#   itself. A pixel's age is one more than its labeling pixel if the pixel
#   has the same level as itself.
#
########################################################################

@cython.boundscheck(False)
def tw_initial_segmentation(
    np.ndarray image,
    np.ndarray seeds,
    np.ndarray seed_labels,
    np.ndarray structure,
    np.ndarray mask,
    np.ndarray output,
    np.ndarray lower_distance):
    """Compute the initial segmentation on the image tile
    
    image - the image tile / cube / hypercube unraveled. The image should be
            padded with masked pixels in order to simplify the boundary case
            
    seeds - the initial seeds, an NxM array of coordinates of seeds
            (M is the image dimension)
    
    seed_labels - an N-element vector (same size as seeds) the labels for the seeds
            
    structure - a vector of offsets to adjacent pixels in the unraveled image
    
    mask - a mask whose value is zero if a pixel is not to be labeled
    
    output - an array to receive the labels.
    
    lower_distance - an M+1 dimensional array with the last dimension being 2.
                     This is a scratchpad to hold the level (last index = 0)
                     and the age (last index = 1) of each pixel.
    """
    cdef:
        Heap *hp = heap_from_numpy2(TW_IDX_COUNT)
        Heapitem4 elem, newelem
        int i
        DTYPE_INT32_t label
        np.ndarray fs_image, fs_mask, fs_output, fs_lower_distance
        int ld_level_stride = (
            PyArray_STRIDE(lower_distance, PyArray_NDIM(lower_distance)-1) /
            PyArray_ITEMSIZE(lower_distance))
        

    #
    # Prevent segfaults using assertions of input assumptions
    #
    assert(len(seeds) == len(seed_labels))
    assert(seeds.shape[1] == image.ndim)
    assert(structure.shape[1] == image.ndim)
    __assert_same_shape(image, mask)
    __assert_same_shape(image, output)
    assert(image.ndim+1 == lower_distance.ndim)
    for i from 0<= i <image.ndim:
        assert(image.shape[i] == lower_distance.shape[i])
    assert(lower_distance.shape[lower_distance.ndim-1] >= 2)
    assert(np.all(~ mask[0, :]))
    assert(np.all(~ mask[-1, :]))
    assert(np.all(~ mask[:, 0]))
    assert(np.all(~ mask[:, -1]))
    
    lower_distance.fill(np.iinfo(lower_distance.dtype).max)
    fs_image = np.zeros(len(structure), DTYPE_INT32)
    fs_mask = np.zeros(len(structure), DTYPE_INT32)
    fs_output = np.zeros(len(structure), DTYPE_INT32)
    fs_lower_distance = np.zeros(len(structure), DTYPE_INT32)
    compute_flat_structure(image, structure, fs_image)
    compute_flat_structure(mask, structure, fs_mask)
    compute_flat_structure(output, structure, fs_output)
    compute_flat_structure(lower_distance, structure, fs_lower_distance)
    #with nogil:
    if True:
        for i from 0 <= i < PyArray_SHAPE(seeds)[0]:
            elem.index[TW_IMAGE_IDX] = compute_flat(image, seeds, i)
            elem.index[TW_MASK_IDX] = compute_flat(mask, seeds, i)
            elem.index[TW_OUTPUT_IDX] = compute_flat(output, seeds, i)
            elem.index[TW_LD_IDX] = compute_flat(lower_distance, seeds, i)
            #
            # Mark a seed as an absolute minimum.
            #
            elem.value = 0
            elem.age = 0
            ELEM2LD(lower_distance, &elem, ld_level_stride)
            PUTELEMENT(output, elem.index[TW_OUTPUT_IDX], 
                       GETELEMENT1(seed_labels, i))
            elem.value = GETELEMENT(image, elem.index[TW_IMAGE_IDX])
            heappush(hp, <Heapitem *>&elem)
        
        while hp.items > 0:
            heappop(hp, <Heapitem *>&elem)
            label = GETELEMENT(output, elem.index[TW_OUTPUT_IDX])
            for i from 0 <= i < PyArray_SHAPE(structure)[0]:
                if conditional_newelem(
                    &elem, &newelem, i, label, ld_level_stride,
                    image, mask, output, lower_distance,
                    fs_image, fs_mask, fs_output, fs_lower_distance):
                    heappush(hp, <Heapitem *>&newelem)
    heap_done(hp)
        
@cython.boundscheck(False)
def tw_propagate_border(
    np.ndarray image,
    np.ndarray seeds,
    np.ndarray seed_mask_index,
    np.ndarray seed_structure_mask,
    np.ndarray structure,
    np.ndarray mask,
    np.ndarray output,
    np.ndarray lower_distance):
    """Propagate labels from the image border into the image
    
    image - the tile / cube / hypercube
    
    seeds - an N x M matrix of coordinates of seeds
    
    seed_mask_index - an index into the seed structure mask for every seed
    
    seed_structure_mask - an N x M boolean array that supplies the mask to use
                          when iterating through the structure. Out-of-bounds
                          structure offsets have a value of False. N is the
                          number of different masks, M is the dimensions of
                          the structure vector.
    
    structure - the connectivity structure, an NxM matrix of N rows of
                structure coordinate offsets.
    
    mask - mask of pixels not to process, including the edge
    
    output - the partially-completed labels, which may change
    
    lower_distance - the distance to the closest lower pixel
    
    returns True if an internal value changed.
    """
    cdef:
        Heap *hp = heap_from_numpy2(TW_IDX_COUNT)
        Heapitem4 elem, newelem
        int i, j
        int seed_structure_index
        DTYPE_INT32_t label, value
        np.ndarray fs_image, fs_mask, fs_output, fs_lower_distance
        int ld_level_stride = (
            PyArray_STRIDE(lower_distance, PyArray_NDIM(lower_distance)-1) /
            PyArray_ITEMSIZE(lower_distance))
        
    #
    # Prevent segfaults using assertions of input assumptions
    #
    assert(seeds.shape[1] == image.ndim)
    assert(structure.shape[1] == image.ndim)
    assert(seed_structure_mask.shape[1] == structure.shape[0])
    __assert_same_shape(image, mask)
    __assert_same_shape(image, output)
    assert(image.ndim+1 == lower_distance.ndim)
    for i from 0<=i<image.ndim:
        assert(image.shape[i] == lower_distance.shape[i])
    assert(np.all(~ mask[0, :]))
    assert(np.all(~ mask[-1, :]))
    assert(np.all(~ mask[:, 0]))
    assert(np.all(~ mask[:, -1]))
    
    fs_image = np.zeros(len(structure), DTYPE_INT32)
    fs_mask = np.zeros(len(structure), DTYPE_INT32)
    fs_output = np.zeros(len(structure), DTYPE_INT32)
    fs_lower_distance = np.zeros(len(structure), DTYPE_INT32)
    compute_flat_structure(image, structure, fs_image)
    compute_flat_structure(mask, structure, fs_mask)
    compute_flat_structure(output, structure, fs_output)
    compute_flat_structure(lower_distance, structure, fs_lower_distance)
    
    with nogil:
        for i from 0 <= i < PyArray_SHAPE(seeds)[0]:
            elem.index[TW_IMAGE_IDX] = compute_flat(image, seeds, i)
            elem.index[TW_MASK_IDX] = compute_flat(mask, seeds, i)
            elem.index[TW_OUTPUT_IDX] = compute_flat(output, seeds, i)
            elem.index[TW_LD_IDX] = compute_flat(lower_distance, seeds, i)
            seed_structure_index = GETELEMENT1(seed_mask_index, i)
            LD2ELEM(lower_distance, &elem, ld_level_stride)
            #
            # either propagate the current level or reset if we've gone up a level
            #
            value = GETELEMENT(image, elem.index[TW_IMAGE_IDX])
            if value > elem.value:
                elem.age = 0
                elem.value = value
            label = GETELEMENT(output, elem.index[TW_OUTPUT_IDX])
            for j from 0 <= j < PyArray_SHAPE(structure)[0]:
                if GETBOOL2(seed_structure_mask, seed_structure_index, j):
                    if conditional_newelem(
                        &elem, &newelem, j, label, ld_level_stride,
                        image, mask, output, lower_distance,
                        fs_image, fs_mask, fs_output, fs_lower_distance):
                        heappush(hp, <Heapitem *>&newelem)
    if hp.items == 0:
        heap_done(hp)
        return False
    
    with nogil:
        while hp.items > 0:
            heappop(hp, <Heapitem *>&elem)
            label = GETELEMENT(output, elem.index[TW_OUTPUT_IDX])
            for i from 0 <= i < PyArray_SHAPE(structure)[0]:
                if conditional_newelem(
                    &elem, &newelem, i, label, ld_level_stride,
                    image, mask, output, lower_distance,
                    fs_image, fs_mask, fs_output, fs_lower_distance):
                     heappush(hp, <Heapitem *>&newelem)
    heap_done(hp)
    return True

def __assert_same_shape(np.ndarray a, np.ndarray b):
    assert (PyArray_NDIM(a) == PyArray_NDIM(b))
    for i in range(PyArray_NDIM(a)):
        assert (PyArray_SHAPE(a)[i] == PyArray_SHAPE(b)[i])
        
@cython.cdivision(True)        
cdef inline void compute_flat_structure(
    np.ndarray a, np.ndarray structure, np.ndarray flat_structure) nogil:
    #
    # Compute the flat structure for the given array data - the offset
    # to the next element.
    #
    # a - pointer to the array to be addressed
    # structure - an NxM array where N is the number of structuring elements
    #             and M is the number of dimensions of "a". Each element
    #             is the offset in the indexed dimension. The array type should
    #             be np.int32
    # flat_structure - an N element vector to receive the structure indexes
    #                  for a typed (e.g. DTYPE_INT32_t) pointer into the
    #                  array data.
    #
    cdef:
        DTYPE_INT32_t i
        DTYPE_INT32_t index
        int d
        int n_dim = PyArray_SHAPE(structure)[1]
    for i from 0 <= i < PyArray_SHAPE(structure)[0]:
        index=0
        for d from 0 <= d < n_dim:
            index += PyArray_STRIDES(a)[d] * GETELEMENT2(structure, i, d)
        index = index / PyArray_ITEMSIZE(a)
        PUTELEMENT1(flat_structure, i, index)

@cython.cdivision(True)        
cdef inline DTYPE_INT32_t compute_flat(np.ndarray a, np.ndarray coords, int i) nogil:
    #
    # Compute a flat index using the strides of an array
    #
    cdef:    
        int d
        int n_dim = PyArray_SHAPE(coords)[1]
        DTYPE_INT32_t acc = 0
    for d from 0 <= d < n_dim:
        acc += PyArray_STRIDES(a)[d] * GETELEMENT2(coords, i, d)
    return acc / PyArray_ITEMSIZE(a)
    
cdef inline int conditional_newelem(
    Heapitem4 *elem, Heapitem4 *newelem, 
    int i,
    DTYPE_INT32_t label,
    int ld_level_stride,
    np.ndarray image,
    np.ndarray mask,
    np.ndarray output,
    np.ndarray lower_distance,
    np.ndarray fs_image,
    np.ndarray fs_mask,
    np.ndarray fs_output,
    np.ndarray fs_lower_distance) nogil:
    cdef:
        DTYPE_INT32_t value, newvalue
        
    newelem.index[TW_MASK_IDX] = \
        elem.index[TW_MASK_IDX] + GETELEMENT1(fs_mask, i)
    if GETBOOL(mask, newelem.index[TW_MASK_IDX]):
        newelem.index[TW_IMAGE_IDX] = \
            elem.index[TW_IMAGE_IDX] + GETELEMENT1(fs_image, i)
        newelem.index[TW_OUTPUT_IDX] = \
            elem.index[TW_OUTPUT_IDX] + GETELEMENT1(fs_output, i)
        newelem.index[TW_LD_IDX] = \
            elem.index[TW_LD_IDX] + GETELEMENT1(fs_lower_distance, i)
        LD2ELEM(lower_distance, newelem, ld_level_stride)
        #
        # If the current best for the new element is better thant or same as
        # what would be replaced by the propagating element, quit.
        #
        if newelem.value < elem.value or (
           newelem.value == elem.value and newelem.age <= elem.age+1):
            return 0
        #
        # Save the new best value and age
        #
        newelem.value = elem.value
        newelem.age = elem.age + 1
        ELEM2LD(lower_distance, newelem, ld_level_stride)
        #
        # Decide whether to propagate the propagating element's value or
        # the new element's value. Use the new element's value if higher.
        # When entering onto a higher level, set the age to 0.
        # 
        value = GETELEMENT(image, newelem.index[TW_IMAGE_IDX])
        if value > elem.value:
            newelem.value = value
            newelem.age = 0
        else:
            newelem.value = elem.value
        PUTELEMENT(output, newelem.index[TW_OUTPUT_IDX], label)
        return 1
    else:
        return 0

@cython.cdivision(True)        
cdef inline int index_limit(np.ndarray a):
    cdef:
        int i
        int acc = 0
    for i from 0 <= i < a.ndim:
        acc += (PyArray_SHAPE(a)[i] -1) * PyArray_STRIDE(a, i) / PyArray_ITEMSIZE(a)
    return acc + 1
