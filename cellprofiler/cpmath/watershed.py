"""watershed.py - watershed algorithm

This module implements a watershed algorithm that apportions pixels into
marked basins. The algorithm uses a priority queue to hold the pixels
with the metric for the priority queue being pixel value, then the time
of entry into the queue - this settles ties in favor of the closest marker.

Some ideas taken from
Soille, "Automated Basin Delineation from Digital Elevation Models Using
Mathematical Morphology", Signal Processing 20 (1990) 171-182

The most important insight in the paper is that entry time onto the queue
solves two problems: a pixel should be assigned to the neighbor with the
largest gradient or, if there is no gradient, pixels on a plateau should
be split between markers on opposite sides.

CellProfiler is distributed under the GNU General Public License,
but this file is licensed under the more permissive BSD license.
See the accompanying file LICENSE for details.

Copyright (c) 2003-2009 Massachusetts Institute of Technology
Copyright (c) 2009-2013 Broad Institute
All rights reserved.

Please see the AUTHORS file for credits.

Website: http://www.cellprofiler.org
"""


from _heapq import heapify, heappush, heappop
import numpy as np
import scipy.ndimage
from cellprofiler.cpmath.rankorder import rank_order

import _watershed

def __get_strides_for_shape(shape):
    """Get the amount to multiply at each coord when converting to flat"""
    lshape = list(shape)
    lshape.reverse()
    stride = [1]
    for i in range(len(lshape)-1):
        stride.append(lshape[i]*stride[i])
    stride.reverse()
    return np.array(stride)

def __heapify_markers(markers,image):
    """Create a priority queue heap with the markers on it"""
    stride = __get_strides_for_shape(image.shape)
    coords = np.argwhere(markers != 0)
    ncoords= coords.shape[0]
    if ncoords > 0:
        pixels = image[markers != 0]
        age    = np.array(range(ncoords))
        offset = np.zeros(coords.shape[0],int)
        for i in range(image.ndim):
            offset = offset + stride[i]*coords[:,i]
        pq = np.column_stack((pixels, age, offset, coords))
        ordering = np.lexsort((age,pixels)) # pixels = top priority, age=second
        pq = pq[ordering,:]
    else:
        pq = np.zeros((0,markers.ndim+3),int)
    return (pq,ncoords)
    

def fast_watershed(image, markers, connectivity=None, offset=None, mask=None):
    """Return a matrix labeled using the watershed algorithm
    
    image - an array where the lowest value points are
            labeled first.
    markers - an array marking the basins with the values
              to be assigned in the label matrix. Zero means not a marker.
              This array should be of an integer type.
    connectivity - an array whose non-zero elements indicate neighbors
                   for connection.
                   Following the scipy convention, default is a one-connected
                   array of the dimension of the image.
    offset  - offset of the connectivity (one offset per dimension)
    mask    - don't label points in the mask
    
    Returns a labeled matrix of the same type and shape as markers
    
    This implementation converts all arguments to specific, lowest common
    denominator types, then passes these to a C algorithm that operates
    as above.
    """
    
    if connectivity == None:
        c_connectivity = scipy.ndimage.generate_binary_structure(image.ndim, 1)
    else:
        c_connectivity = np.array(connectivity,bool)
        if c_connectivity.ndim != image.ndim:
            raise ValueError,"Connectivity dimension must be same as image"
    if offset == None:
        if any([x%2==0 for x in c_connectivity.shape]):
            raise ValueError,"Connectivity array must have an unambiguous center"
        #
        # offset to center of connectivity array
        #
        offset = np.array(c_connectivity.shape)/2

    # pad the image, markers, and mask so that we can use the mask to keep from running off the edges
    pads = offset

    def pad(im):
        new_im = np.zeros([i + 2*p for i,p in zip(im.shape, pads)], im.dtype)
        new_im[[slice(p, -p,None) for p in pads]] = im
        return new_im

    if mask is not None:
        mask = pad(mask)
    else:
        mask = pad(np.ones(image.shape, bool))
    image = pad(image)
    markers = pad(markers)

    c_image = rank_order(image)[0].astype(np.int32)
    c_markers = np.ascontiguousarray(markers,dtype=np.int32)
    if c_markers.ndim!=c_image.ndim:
        raise ValueError,\
            "markers (ndim=%d) must have same # of dimensions "\
            "as image (ndim=%d)"%(c_markers.ndim, c_image.ndim)
    if not all([x==y for x,y in zip(c_markers.shape, c_image.shape)]):
        raise ValueError("image and markers must have the same shape")
    if mask!=None:
        c_mask = np.ascontiguousarray(mask,dtype=bool)
        if c_mask.ndim!=c_markers.ndim:
            raise ValueError, "mask must have same # of dimensions as image"
        if not all([x==y for x,y in zip(c_markers.shape, c_mask.shape)]):
            raise ValueError, "mask must have same shape as image"
        c_markers[np.logical_not(mask)]=0
    else:
        c_mask = None
    c_output = c_markers.copy()

    #
    # We pass a connectivity array that pre-calculates the stride for each
    # neighbor.
    #
    # The result of this bit of code is an array with one row per
    # point to be considered. The first column is the pre-computed stride
    # and the second through last are the x,y...whatever offsets
    # (to do bounds checking).
    c = []
    image_stride = __get_strides_for_shape(image.shape)
    for i in range(np.product(c_connectivity.shape)):
        multiplier = 1
        offs = []
        indexes = []
        ignore = True
        for j in range(len(c_connectivity.shape)):
            elems = c_image.shape[j]
            idx   = (i / multiplier) % c_connectivity.shape[j]
            off   = idx - offset[j]
            if off:
                ignore = False
            offs.append(off)
            indexes.append(idx)
            multiplier *= c_connectivity.shape[j]
        if (not ignore) and c_connectivity.__getitem__(tuple(indexes)):
            stride = np.dot(image_stride, np.array(offs))
            offs.insert(0,stride)
            c.append(offs)
    c = np.array(c,np.int32)

    pq,age = __heapify_markers(c_markers, c_image)
    pq = np.ascontiguousarray(pq,dtype=np.int32)
    if np.product(pq.shape) > 0:
        # If nothing is labeled, the output is empty and we don't have to
        # do anything
        c_output = c_output.flatten()
        if c_mask == None:
            c_mask = np.ones(c_image.shape,np.int8).flatten()
        else:
            c_mask = c_mask.astype(np.int8).flatten()
        _watershed.watershed(c_image.flatten(),
                             pq, age, c, 
                             c_image.ndim, 
                             c_mask,
                             np.array(c_image.shape,np.int32),
                             c_output)
    c_output = c_output.reshape(c_image.shape)[[slice(1,-1,None)] * image.ndim]
    try:
        return c_output.astype(markers.dtype)
    except:
        return c_output

watershed = fast_watershed

class Tile(object):
    def __init__(self, slices):
        '''Initialize a tile out of a bigger image
        
        slices - a sequence of two-tuples giving the begin and end of a slice
        '''
        self.slices = np.array(slices)
        self.slice_tuple = tuple([slice(s, e) for s, e in self.slices])
        self.big_slices = np.array(slices)
        self.big_slices[:, 0] -= 1
        self.big_slices[:, 1] += 1
        self.big_slice_tuple = tuple([slice(s, e) for s, e in self.slices])
        self.big_dim = self.big_slices[:, 1] - self.big_slices[:, 0]
        self.active = True
        
    def initial_segmentation(self, image, seeds, labels, structure, mask, 
                             output, lower_distance):
        seeds = (seeds - self.big_slices[:, 0]).astype(np.int32)
        labels = labels.astype(np.int32)
        structure = structure.astype(np.int32)
        if self.__needs_padding(image, mask, output, lower_distance):
            self.__initial_segmentation_padded(
                image, seeds, labels, structure, mask, output, lower_distance)
        else:
            self.__initial_segmentation_mapped(
                image, seeds, labels, structure, mask, output, lower_distance)

    def __needs_padding(self, image, mask, output, lower_distance):
        if np.any(self.slices[:, 0] == 0):
            return True
        if np.any(self.slices[:, 1] == np.array(image.shape)):
            return True
        if any([~ isinstance(x, np.ndarray) 
                for x in image, mask, output, lower_distance]):
            return True
        if any([x.dtype != np.int32 for x in image, output, lower_distance]):
            return True
        return mask is not None and mask.dtype == np.bool8

    def __initial_segmentation_padded(
        self, image, seeds, labels, structure, mask, output, lower_distance):
        pimage, poutput = [
            self.__get_padded(x) for x in image, output]
        plower_distance = self.__get_padded(
            lower_distance, ndim = image.ndim)
        pmask = self.__get_mask(mask)
        _watershed.tw_initial_segmentation(
            pimage, seeds, labels, structure, pmask, poutput, plower_distance)
        self.__put(output, poutput)
        self.__put(lower_distance, plower_distance, lower_distance.ndim-1)

    def __get_padded(self, a, dtype = np.int32, ndim = None):
        if ndim is None:
            ndim = a.ndim
        dim = list(self.big_dim) + list(a.shape[ndim:])
        out = np.zeros(dim, dtype)
        src_slices = []
        dest_slices = []
        for i in range(ndim):
            if self.slices[i, 0] == 0:
                if self.slices[i, 1] == a.shape[i]:
                    src_slices.append(slice(0, self.slices[i, 1]))
                    dest_slices.append(slice(1, -1))
                else:
                    src_slices.append(slice(0, self.slices[i, 1]+1))
                    dest_slices.append(slice(1, out.shape[i]))
            elif self.slices[i, 1] == a.shape[i]:
                src_slices.append(slice(self.slices[i, 0] - 1,
                                        self.slices[i, 1]))
                dest_slices.append(slice(0, out.shape[i]-1))
            else:
                src_slices.append(slice(self.slices[i, 0] -1, 
                                        self.slices[i, 1] + 1))
                dest_slices.append(slice(0, out.shape[i]))
        out[tuple(dest_slices)] = a[tuple(src_slices)]
        return out
    
    def __get_mask(self, mask):
        '''Get a padded slice of the mask and mask the edges'''
        if mask is None:
            pmask = np.ones(self.big_dim, np.bool8)
        else:
            pmask = self.__get_padded(mask, np.bool8)
        for i in range(pmask.ndim):
            slices = [slice(None)] * pmask.ndim
            slices[i] = 0
            pmask[tuple(slices)] = False
            slices[i] = pmask.shape[i] - 1
            pmask[tuple(slices)] = False
        return pmask

    def __put(self, dest, src, ndim=None):
        if ndim is None:
            ndim = dest.ndim
        dest[self.slice_tuple] = src[tuple([slice(1,-1)] * ndim)]
        
    def __initial_segmentation_mapped(
        self, image, seeds, labels, structure, mask, output, lower_distance):
        pmask = self.__get_mask(mask)
        _watershed.tw_initial_segmentation(
            image[self.big_slice_tuple], seeds, labels, structure,
            pmask, output[self.big_slice_tuple],
            lower_distance[self.big_slice_tuple])
        
    def propagate_border(self, image, structure, mask, output, lower_distance):
        if self.__needs_padding(image, mask, output, lower_distance):
            self.__propagate_border_padded(
                image, structure, mask, output, lower_distance)
        else:
            self.__propagate_border_mapped(
                image, structure, mask, output, lower_distance)
        return self.active
        
    def __propagate_border_padded(
        self, image, structure, mask, output, lower_distance):
        seeds, seed_mask_indexes, structure_masks = self.__find_edges(
            image, structure, mask, lower_distance)
        if seeds is None:
            # There were no seed coordinates
            self.active = False
            return
        pimage, poutput = [
            self.__get_padded(x) for x in image, output]
        plower_distance = self.__get_padded(
            lower_distance, ndim=image.ndim)
        pmask = self.__get_mask(mask)
        self.active = _watershed.tw_propagate_border(
            pimage, seeds, seed_mask_indexes, structure_masks, structure,
            pmask, poutput, plower_distance)
        self.__put(output, poutput)
        self.__put(lower_distance, plower_distance, lower_distance.ndim-1)
            
    def __propagagate_border_mapped(
        self, image, structure, mask, output, lower_distance):
        seeds, seed_mask_indexes, structure_masks = self.__find_edges(
            image, structure, mask, lower_distance)
        pmask = self.__get_mask(mask)
        self.active = _watershed.tw_propagate_border(
            image, seeds, seed_mask_indexes, structure_masks, structure,
            pmask, output, lower_distance)
        
    def __find_edges(self, image, structure, mask, lower_distance):
        structure_masks = []
        edges = []
        #
        # 3 states for each dimension: internal = None, 0 or -1. Create
        # all possible combinations of 3 states across all dimensions.
        #
        slice_none = slice(None)
        state = [None] * image.ndim
        slices = [slice_none] * image.ndim
        while True:
            skip = False
            for i, symbol in enumerate(state):
                if symbol is None:
                    state[i] = 0
                    slices[i] = 0
                    if self.slices[i, 0] == 0:
                        skip = True
                    break
                elif symbol is 0:
                    state[i] = -1
                    slices[i] = -1
                    if self.slices[i, 1] == image.shape[i]:
                        skip = True
                    break
                else:
                    state[i] = None
                    slices[i] = slice_none
            else:
                break
            if skip:
                continue
            #
            # Compute the mask of allowable structuring elements for these
            # edges
            #
            structure_mask = np.ones(structure.shape[0], bool)
            for i in range(structure.shape[0]):
                for j, symbol in enumerate(state):
                    if symbol is 0 and structure[i, j] <= 0:
                        structure_mask[i] = False
                        break
                    elif symbol is -1 and structure[i, j] >= 0:
                        structure_mask[i] = False
            #
            # Create slices for each of the dimensions in the edge.
            # 0 and -1 yield single-coordinate slices at the extrema
            # None yields a slice running the length of the edge.
            # "edges" can have zero (= corner) or more (2 = face) dimensions
            #
            coords = []
            for i, symbol in enumerate(state):
                if symbol == 0:
                    coords.append(slice(self.big_slices[i, 0], 
                                        self.big_slices[i, 0]+1))
                elif symbol == -1:
                    coords.append(slice(self.big_slices[i, 1]-1, 
                                        self.big_slices[i, 1]))
                else:
                    coords.append(slice(self.big_slices[i, 0]+1,
                                        self.big_slices[i, 1]-1))
            coords = tuple(coords)
            #
            # Use the coordinates to extract the mask values. The result
            # can be used to exclude masked pixels from the edge
            #
            if mask is not None:
                emask = (mask[coords] & (lower_distance[coords] != -1)).flatten()
            else:
                emask = lower_distance[coords] != -1
            if np.any(emask):
                #
                # flatten the coordinate mgrid and adjust coordinates to
                # the tile's padded reference frame
                #
                coord_columns = tuple([
                    x.flatten() - self.big_slices[i, 0]
                    for i, x in enumerate(np.mgrid[coords])])
                edge = np.column_stack(coord_columns)
                if mask is not None:
                    edge = edge[emask]
                    structure_mask = structure_mask[emask]
                structure_masks.append(structure_mask)
                edges.append(edge)
        if len(edges) == 0:
            return None, None, None
        edge_array = np.vstack(edges)
        edge_indexes = np.hstack([np.ones(len(e), np.int32) * i
                                  for i, e in enumerate(edges)])
        structure_mask_array = np.vstack(structure_masks)
        return edge_array, edge_indexes, structure_mask_array
            
def tiled_watershed(image, seed_locations, labels,
                    structure, output, lower_distance = None, 
                    tile_shape = None, mask=None):
    """Label an image using the watershed algorithm and tiles
    
    image - a matrix where the lowest value points are
            labeled first.
    seed_locations - an Nx(image.ndim) array of locations of the watershed
                     seeds, one seed per row
    seed_labels - the labels to be assigned to those seeds
    structure - an Nx(image.ndim) array of offsets for the structuring element.
                For instance, a 2-d four-connected structure is
                [[-1, 0], [0, 1], [1, 0], [0, -1]]
    output - the array that receives the labeling
    lower_distance - a scratchpad for storing one intermediate np.int32 per pixel
                     if None, allocate one
    tile_shape - the shape of a tile. If None, we use megapixel tiles.
    mask    - don't label points in the mask
    """
    if tile_shape is None:
        # Break everything into 16 mb blocks by default.
        # 2d = 4K x 4k
        # 3d = 256 x 256 x 256
        # 4d = 64 x 64 x 64 x 64
        # 5d cry
        assert image.ndim < 5, "Please supply a tile shape"
        tile_shape = tuple([2**(24/image.ndim)]*image.ndim)
    tile_shape = np.array(tile_shape)
    if lower_distance == None:
        lower_distance = np.zeros(list(image.shape)+[2], np.int32)
    n_tiles = ((np.array(image.shape) + tile_shape - 1) / tile_shape).astype(int)
    tile_grid = np.mgrid[
        tuple([slice(0, end, ts) for end, ts in zip(image.shape, tile_shape)])]
    tile_coords = np.column_stack([x.flatten() for x in tile_grid])
    tile_ends = np.minimum(tile_coords + tile_shape, np.array(image.shape))
    tiles = [
        Tile(np.column_stack((c, e))) for c, e in zip(tile_coords, tile_ends)]
    tile_array = np.array(tiles).reshape(n_tiles)
    
    stride = np.hstack(([1], np.cumprod(n_tiles[-2::-1])))[::-1]
    seed_address = np.sum(stride * (seed_locations / tile_shape).astype(int), 1)
    for i, tile in enumerate(tiles):
        tseed_locs = np.where(seed_address == i)
        tseeds = seed_locations[tseed_locs]
        tlabels = labels[tseed_locs]
        tile.initial_segmentation(image, tseeds, tlabels, structure, mask,
                                  output, lower_distance)
    
    active_tiles = np.ones(n_tiles, bool)
    while np.any(active_tiles):
        at_coords = np.argwhere(active_tiles)
        active_tiles.fill(False)
        for coords in at_coords:
            tile = tile_array[tuple(coords)]
            if tile.propagate_border(image, structure, mask, output,
                                     lower_distance):
                for scoords in np.sign(structure).astype(int):
                    target_coords = scoords + coords
                    if np.any((target_coords < 0) | 
                              (target_coords >= np.array(active_tiles.shape))):
                        continue
                    active_tiles[tuple(target_coords)] = True
