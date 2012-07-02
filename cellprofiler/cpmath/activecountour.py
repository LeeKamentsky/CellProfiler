'''activecontour.py Evolve a segmentation using active contours

This module uses the method described in:

Chan, Vese, "Active Contours Without Edges", 
IEEE Transactions on Image Processing, vol. 10, # 2, Febuary 2001, p 266

The method takes a binary segmentation and an image as input. It assumes
that the mean intensity of the image in the background differs from
that in the foreground and it weights the balance between keeping the
original segmentation and maximizing the foreground / background difference
with a single term, "alpha", which is is higher to pay more attention to
image intensity and lower to honor the original segmentation.

The code in the module was derived from a Matlab implementation by
Shawn Lankton (www.shawnlankton.com) which was contributed to Matlab Central
without copyright.

CellProfiler is distributed under the GNU General Public License,
but this file is licensed under the more permissive BSD license.
See the accompanying file LICENSE for details.

Copyright (c) 2003-2009 Massachusetts Institute of Technology
Copyright (c) 2009-2012 Broad Institute
All rights reserved.

Please see the AUTHORS file for credits.

Website: http://www.cellprofiler.org
'''

import numpy as np
from scipy.ndimage import distance_transform_edt

def active_contour(image, init_mask, max_its, alpha, fixed_threshold = False):
    '''Adjust a segmentation using active contours
    
    image - intensity image
    
    init_mask - binary segmentation of foreground from background
    
    max_its - maximum # of iterations of gradient descent
    
    alpha - weight of smoothing term. Higher attempts to maximize difference more.
    
    fixed_threshold - if true, "u" and "v" = the mean values of the foreground
                      and background are calculated on the initial segmentation
                      in order to stabilize the threshold. If false,
                      "u" and "v" are updated with each iteration, so the
                      threshold can drift. (an addition to the original code)
    
    returns the evolved mask.
    
    Please see the module header for citations.
    '''
    
    phi = mask2phi(init_mask)
    
    for iteration in range(max_its):
        # Get the curve's narrow band.
        idx = (phi <= 1.2) & (phi >= -1.2)

        # The current segmentation - interior == True
        current_segmentation = (phi <= 0)
        if np.all(current_segmentation) or not np.any(current_segmentation):
            return current_segmentation # all foreground or background
        
        if iteration > 0:
            print "Pass %d: bg -> fg = %d, fg -> bg = %d " % (
                iteration, np.sum(current_segmentation & ~ last_mask),
                np.sum(last_mask & ~ current_segmentation))
            
        
        if iteration == 0 or not fixed_threshold:
            u = np.mean(image[current_segmentation])
            v = np.mean(image[~current_segmentation])
        
        # The force from image information
        F = (image[idx] - u) ** 2 - (image[idx] - v) ** 2
        
        curvature = get_curvature(phi, idx)
        
        # Gradient descent to minimize energy
        dphidt = F / np.max(np.abs(F)) + alpha * curvature

        # Maintain the CFL condition, exiting in rare case that we are
        # at the exact local minimum
        
        max_dphidt = np.max(dphidt)
        if max_dphidt == 0:
            return current_segmentation
        
        dt = .45 / max_dphidt
        
        phi[idx] += dphidt * dt
        
        # Keep SDF smooth
        phi = sussman(phi, .5)
        
        last_mask = current_segmentation
        
    return current_segmentation

def mask2phi(init_mask):
    '''Convert the mask to a distance from the foreground / background boundary
    
    init_mask - segmentation mask
    
    returns a distance image where negative distances are distances in
    the foreground from the background and positive are vice-versa.
    '''
    # Note that polarity of Scipy and Matlab distance transforms are reversed
    # The last term puts the edge halfway between the foreground and background
    # edge pixels.
    return (distance_transform_edt(~ init_mask) - 
            distance_transform_edt(init_mask) + 
            init_mask - .5)

def get_curvature(phi, idx):
    '''Compute the curvature in the evolving distance function
    
    phi - SDF
    idx - a mask of the points in the boundary
    '''
    #
    # This is a longhand computation of Kappa, the curvature.
    # The curvature can't be computed on the boundaries - there's always
    # one unknown derivative term. So we exclude and return a value of zero for
    # those.
    #
    i, j = np.mgrid[0:phi.shape[0], 0:phi.shape[1]]
    good = ((i > 0) & (i < phi.shape[0] - 1) &
            (j > 0) & (j < phi.shape[1] - 1))
    good_idx = good[idx]
    if np.all(~ good_idx):
        return np.zeros(len(idx))
    
    y = i[idx & good]
    x = j[idx & good]
    
    ctr = phi[y, x]
    up = phi[y+1, x]
    dn = phi[y-1, x]
    lt = phi[y, x-1]
    rt = phi[y, x+1]
    ul = phi[y+1, x-1]
    ur = phi[y+1, x+1]
    dl = phi[y-1, x-1]
    dr = phi[y-1, x+1]
    
    # Get central derivatives of SDF at x,y
    phi_x  = -lt + rt
    phi_y  = -dn + up
    phi_xx = lt - 2 * ctr + rt
    phi_yy = dn - 2 * ctr + up
    phi_xy = (-dl - ur + dr + ul) / 4.
    phi_x2 = phi_x * phi_x
    phi_y2 = phi_y * phi_y
    
    # Compute curvature - I'm differing from Lankton here - he divides
    # by the gradient squared (phi_x2 + phi_y2) ** .5 / (phi_x2 + phi_y2) ** 1.5
    # and there's no reason I can think of for computing it that way.
    # I've simplified the equation.
    #
    curvature = (
        (phi_x2 * phi_yy + phi_y2 * phi_xx - 2 * phi_x * phi_y * phi_xy) /
        (phi_x2 + phi_y2 + np.finfo(phi.dtype).eps))
    #
    # Account for items that fell on the boundary
    #
    result = np.zeros(len(good_idx), curvature.dtype)
    result[good_idx] = curvature
    return result

def sussman(D, dt):
    """Level set re-initialization by the Sussman method
    
    D - level set
    dt - evolution gradient
    
    I think this is the reference:
    Sussman, "An Adaptive Level Set Approach for Incompressible Two-phase Flows",
    https://seesar.lbl.gov/anag/publications/colella/LBNL-40327-02.pdf
    
    Lankton's code seems to be a direct copy of a version of Sussman
    contributed by Romeil Sandhu under a BSD license. The code seems to differ
    from Sussman who computes a = x[:-1,:] - 2*x + x[1:, :] and
    b = x[:-2, :] - 2*x[:-1, :] + x, taking the absolute maximum, but this
    approach seems to be simpler, to try and make the value at x be 1+ the
    value nearby.
    
    In section 3.5: Redistance Operation, Sussman describes how to reinitialize
    the level set to maintain it as a distance function. A positive value
    should be one more than the smallest positive value surrounding it
    and likewise for negative. The equation below maintains that to a first
    order.
    
    The original code doesn't filter negative values which will happen
    for pixels right on the foreground / backdround border. That's done
    by giving dD a floor of zero.
    """
    #
    # The padded difference in the X direction. Ends are zero.
    dx = np.zeros((D.shape[0], D.shape[1]+1), D.dtype)
    dx[:, 1:-1] = D[:, 1:] - D[:, :-1]
    a = dx[:, :-1]
    b = dx[:, 1:]

    dy = np.zeros((D.shape[0]+1, D.shape[1]), D.dtype)
    dy[1:-1, :] = D[1:, :] - D[:-1, :]
    c = dy[:-1, :]
    d = dy[1:, :]
    
    dD = (np.sqrt(np.maximum(a * a * ((a > 0) == (D > 0)),
                             b * b * ((b < 0) == (D > 0))) +
                  np.maximum(c * c * ((c > 0) == (D > 0)),
                             d * d * ((d < 0) == (D > 0)))) - 1)
    dD = np.maximum(dD, 0)
    return D - dt * dD * D / np.sqrt(D * D + 1)


if __name__=="__main__":
    import sys
    from matplotlib.image import pil_to_array
    import PIL.Image
    
    image = np.flipud(pil_to_array(PIL.Image.open(sys.argv[1])))
    if image.ndim == 3:
        image = image[:, :, 0]
    thresh = int(sys.argv[2])
    seg = image > thresh
    result = active_contour(image.astype(float), seg, int(sys.argv[3]), float(sys.argv[4]))
    pass