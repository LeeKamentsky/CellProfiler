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
cdef extern from "stdlib.h":
   ctypedef unsigned long size_t
   void free(void *ptr) nogil
   void *malloc(size_t size) nogil
   void *realloc(void *ptr, size_t size) nogil
   void *memcpy(void *dest, void *src, size_t size) nogil

cdef struct Heap:
    unsigned int items
    unsigned int space
    unsigned int n_indexes
    void *data
    Heapitem **ptrs

cdef inline Heap *heap_from_numpy2(int n_indexes) nogil:
    cdef unsigned int k
    cdef Heap *heap 
    cdef size_t elem_size = sizeof(Heapitem) + sizeof(np.int32_t) * n_indexes
    heap = <Heap *> malloc(sizeof (Heap))
    heap.items = 0
    heap.space = 1000
    heap.n_indexes = n_indexes
    heap.data = malloc(heap.space * elem_size)
    heap.ptrs = <Heapitem **> malloc(heap.space * sizeof(Heapitem *))
    for k from 0 <= k < heap.space:
        heap.ptrs[k] = <Heapitem *>(<char *>heap.data + k * elem_size)
    return heap

cdef inline void heap_done(Heap *heap) nogil:
   free(heap.data)
   free(heap.ptrs)
   free(heap)

cdef inline void swap(unsigned int a, unsigned int b, Heap *h) nogil:
    h.ptrs[a], h.ptrs[b] = h.ptrs[b], h.ptrs[a]


######################################################
# heappop - inlined
#
# pop an element off the heap, maintaining heap invariant
# 
# Note: heap ordering is the same as python heapq, i.e., smallest first.
######################################################
cdef inline void heappop(Heap *heap,
                  Heapitem *dest) nogil:
    cdef unsigned int i, smallest, l, r # heap indices
    
    #
    # Start by copying the first element to the destination
    #
    memcpy(dest, heap.ptrs[0], 
           sizeof(Heapitem) + heap.n_indexes * sizeof(np.int32_t))
    heap.items -= 1

    # if the heap is now empty, we can return, no need to fix heap.
    if heap.items == 0:
        return

    #
    # Move the last element in the heap to the first.
    #
    swap(0, heap.items, heap)

    #
    # Restore the heap invariant.
    #
    i = 0
    smallest = i
    while True:
        # loop invariant here: smallest == i
        
        # find smallest of (i, l, r), and swap it to i's position if necessary
        l = i*2+1 #__left(i)
        r = i*2+2 #__right(i)
        if l < heap.items:
            if smaller(heap.ptrs[l], heap.ptrs[i]):
                smallest = l
            if r < heap.items and smaller(heap.ptrs[r], heap.ptrs[smallest]):
                smallest = r
        else:
            # this is unnecessary, but trims 0.04 out of 0.85 seconds...
            break
        # the element at i is smaller than either of its children, heap invariant restored.
        if smallest == i:
                break
        # swap
        swap(i, smallest, heap)
        i = smallest
        
##################################################
# heappush - inlined
#
# push the element onto the heap, maintaining the heap invariant
#
# Note: heap ordering is the same as python heapq, i.e., smallest first.
##################################################
cdef inline void heappush(Heap *heap,
                          Heapitem *new_elem) nogil:
  cdef unsigned int child         = heap.items
  cdef unsigned int parent
  cdef unsigned int k
  cdef void *new_data
  cdef size_t elem_size = sizeof(Heapitem) + sizeof(np.int32_t) * heap.n_indexes
  cdef size_t index

  # grow if necessary
  if heap.items == heap.space:
      heap.space = heap.space * 2
      new_data = <Heapitem *> realloc(<void *> heap.data, <size_t> (heap.space * elem_size))
      heap.ptrs = <Heapitem **> realloc(<void *> heap.ptrs, <size_t> (heap.space * sizeof(Heapitem *)))
      for k from 0 <= k < heap.items:
          index = <char *>heap.ptrs[k] - <char*>heap.data
          heap.ptrs[k] = <Heapitem *>(<char *>new_data+index)
      for k from heap.items <= k < heap.space:
          heap.ptrs[k] = <Heapitem *>(<char *>new_data + k * elem_size)
      heap.data = new_data

  # insert new data at child
  memcpy(heap.ptrs[child], new_elem, elem_size)
  heap.items += 1

  # restore heap invariant, all parents <= children
  while child>0:
      parent = (child + 1) / 2 - 1 # __parent(i)
      
      if smaller(heap.ptrs[child], heap.ptrs[parent]):
          swap(parent, child, heap)
          child = parent
      else:
          break
