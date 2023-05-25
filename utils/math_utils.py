import numpy as np

def linterpsolve(xvals, yvals, xpoint):
    '''
    'solves' a monotonic function described by xvals and yvals by 
    linearly interpolating between the points above and below xpoint
    xvals, yvals: 1D arrays
    xpoint: float
    '''
    if np.all(np.diff(xvals) >= 0.):
        incr = True
    elif np.all(np.diff(xvals) <= 0.):
        incr = False
    else:
        raise ValueError('linterpsolve only works for monotonic functions')
    ind1 = np.where(xvals <= xpoint)[0]
    ind2 = np.where(xvals >= xpoint)[0]
    if len(ind2) == 0 or len(ind1) == 0:
        raise ValueError('xpoint is outside the bounds of xvals')
    if incr:
        ind1 = np.max(ind1)
        ind2 = np.min(ind2)
    else:
        ind1 = np.min(ind1)
        ind2 = np.max(ind2)
    if ind1 == ind2:
        ypoint = yvals[ind1]
    else:
        w = (xpoint - xvals[ind1]) / (xvals[ind2] - xvals[ind1])
        ypoint = yvals[ind2] * w + yvals[ind1] * (1. - w)
    return ypoint

def find_intercepts(yvals, xvals, ypoint, xydct=None):
    '''
    'solves' a monotonic function described by xvals and yvals by
    linearly interpolating between the points above and below ypoint 
    xvals, yvals: 1D arrays
    ypoint: float
    Does not distinguish between intersections separated by less than 2
    xvals points
    '''
    if xvals is None:
        xvals = xydct['x']
    if yvals is None:
        yvals = xydct['y']

    if not (np.all(np.diff(xvals) <= 0.) or np.all(np.diff(xvals) >= 0.)):
        print('linterpsolve only works for monotonic x values')
        return None
    zerodiffs = yvals - ypoint
    leqzero = np.where(zerodiffs <= 0.)[0]
    if len(leqzero) == 0:
        return np.array([])
    elif len(leqzero) == 1:
        edges = [[leqzero[0], leqzero[0]]]
    else:
        segmentedges = np.where(np.diff(leqzero) > 1)[0] + 1
        if len(segmentedges) == 0: # one dip below zero -> edges are intercepts
            edges = [[leqzero[0], leqzero[-1]]]
        else:
            parts = [leqzero[: segmentedges[0]] if si == 0 else \
                     leqzero[segmentedges[si - 1] : segmentedges[si]] 
                     if si < len(segmentedges) else\
                     leqzero[segmentedges[si - 1] :] \
                     for si in range(len(segmentedges) + 1)]
            edges = [[part[0], part[-1]] for part in parts]
    intercepts = [[linterpsolve(zerodiffs[ed[0]-1: ed[0] + 1], xvals[ed[0]-1: ed[0] + 1], 0.),\
                   linterpsolve(zerodiffs[ed[1]: ed[1] + 2],   xvals[ed[1]: ed[1] + 2], 0.)]  \
                  if ed[0] != 0 and ed[1] != len(yvals) - 1 else \
                  [None,\
                   linterpsolve(zerodiffs[ed[1]: ed[1] + 2],   xvals[ed[1]: ed[1] + 2], 0.)] \
                  if ed[1] != len(yvals) - 1 else \
                  [linterpsolve(zerodiffs[ed[0]-1: ed[0] + 1], xvals[ed[0]-1: ed[0] + 1], 0.),\
                   None]  \
                  if ed[0] != 0 else \
                  [None, None]
                 for ed in edges]
    intercepts = [i for i2 in intercepts for i in i2]
    if intercepts[0] is None:
        intercepts = intercepts[1:]
    if intercepts[-1] is None:
        intercepts = intercepts[:-1]
    return np.array(intercepts)

def percentiles_from_histogram(histogram, edgesaxis, axis=-1, 
        percentiles=np.array([0.1, 0.25, 0.5, 0.75, 0.9])):
    '''
    get percentiles from the histogram along axis
    edgesaxis are the bin edges along that same axis
    histograms can be weighted by something: 
    this function just solves
    cumulative distribution == percentiles
    '''
    percentiles = np.array(percentiles)
    if not np.all(percentiles >= 0.) and np.all(percentiles <= 1.):
        msg = ('Input percentiles shoudl be fractions in the range [0, 1].'
               f'They were {percentiles}')
        raise ValueError(msg)
    cdists = np.cumsum(histogram, axis=axis, dtype=np.float) 
    sel = list((slice(None, None, None),) * len(histogram.shape))
    sel2 = np.copy(sel)
    sel[axis] = -1
    sel2[axis] = np.newaxis
    # normalised cumulative dist: divide by total along axis
    cdists /= (cdists[tuple(sel)])[tuple(sel2)] 
    # bin-edge corrspondence: at edge 0, cumulative value is zero
    # histogram values are counts in cells 
    # -> hist bin 0 is what is accumulated between edges 0 and 1
    # cumulative sum: counts in cells up to and including the current one: 
    # if percentile matches cumsum in cell, 
    # the percentile value is it's right edge -> edge[cell index + 1]
    # effectively, if the cumsum is prepended by zeros, 
    # we get a hist bin matches edge bin matching

    oldshape1 = list(histogram.shape)[:axis] 
    oldshape2 = list(histogram.shape)[axis + 1:]
    newlen1 = int(np.prod(oldshape1))
    newlen2 = int(np.prod(oldshape2))
    axlen = histogram.shape[axis]
    cdists = cdists.reshape((newlen1, axlen, newlen2))
    cdists = np.append(np.zeros((newlen1, 1, newlen2)), cdists, axis=1)
    # should already be true, but avoids fp error issues
    cdists[:, -1, :] = 1.

    leftarr  = cdists[np.newaxis, :, :, :] <= \
        percentiles[:, np.newaxis, np.newaxis, np.newaxis]
    rightarr = cdists[np.newaxis, :, :, :] >= \
        percentiles[:, np.newaxis, np.newaxis, np.newaxis]
    
    leftbininds = np.array([[[np.max(
                                   np.where(leftarr[pind, ind1, :, ind2])[0])
                               for ind2 in range(newlen2)] 
                               for ind1 in range(newlen1)] 
                               for pind in range(len(percentiles))])
    # print leftarr.shape
    # print rightarr.shape
    rightbininds = np.array([[[np.min(
                                   np.where(rightarr[pind, ind1, :, ind2])[0])
                               for ind2 in range(newlen2)] 
                               for ind1 in range(newlen1)] 
                               for pind in range(len(percentiles))])
    # if left and right bins are the same, effictively just choose one
    # if left and right bins are separated by more than one (plateau 
    # edge), this will give the middle of the plateau
    lweights = np.array([[[(cdists[ind1, rightbininds[pind, ind1, ind2],
                                   ind2] \
                             - percentiles[pind]) \
                            / (cdists[ind1, rightbininds[pind, ind1, ind2],
                                      ind2] \
                               - cdists[ind1, leftbininds[pind, ind1, ind2],
                                        ind2]) \
                            if rightbininds[pind, ind1, ind2] \
                                != leftbininds[pind, ind1, ind2] \
                            else 1.
                           for ind2 in range(newlen2)] 
                           for ind1 in range(newlen1)] 
                           for pind in range(len(percentiles))])
                
    outperc = lweights * edgesaxis[leftbininds] \
              + (1. - lweights) * edgesaxis[rightbininds]
    outshape = (len(percentiles),) + tuple(oldshape1 + oldshape2)
    outperc = outperc.reshape(outshape)
    return outperc

def getminmax2d(bins, edges, axis=None, log=True, pixdens=False): 
    # axis = axis to sum over; None -> don't sum over any axes 
    # now works for histgrams of general dimensions
    if axis is None:
        imgtoplot = bins
    else:
        imgtoplot = np.sum(bins, axis=axis)
    if pixdens:
        if axis is None:
            naxis = range(len(edges))
        else:
            if not hasattr(axis, '__len__'):
                saxis = [axis]
            else:
                saxis = axis
             # axes not to sum over
            naxis = list(set(range(len(edges))) - set(saxis))
        naxis.sort() 
        numdims = len(naxis)
         # if bins are log, the log sizes are used 
         # and the enclosed log density is minimised
        binsizes = [np.diff(edges[axisi]) for axisi in naxis]
        baseinds = list((np.newaxis,)*numdims)
        normmatrix = np.prod([(binsizes[ind])[tuple(baseinds[:ind]
                                                    + [slice(None,None,None)]
                                                    + baseinds[ind+1:])
                                              ] 
                              for ind in range(numdims)])
        imgtoplot /= normmatrix
        del normmatrix
    finite = np.isfinite(imgtoplot)
    if log:
        imin = np.min(imgtoplot[np.logical_and(finite, imgtoplot > 0)])
        imax = np.max(imgtoplot[np.logical_and(finite, imgtoplot > 0)])
        imin = np.log10(imin)
        imax = np.log10(imax)
    else:
        imin = np.min(imgtoplot[np.isfinite(imgtoplot)])
        imax = np.max(imgtoplot[np.isfinite(imgtoplot)])
    return imin, imax

def combine_hists(h1, h2, e1, e2, rtol=1e-5, atol=1e-8, add=True):
    '''
    add histograms h1, h2 with the same dimension, after aligning edges
    e1, e2
    add = True -> add histograms, return sum
    add = False -> align histograms, return padded histograms and bins
    
    e1, e2 are sequences of arrays, h1, h2 are arrays
    edgetol specifies what relative/absolute (absolute if one is zero)
    differences between edge values are acceptable to call bin edges 
    equal
    
    if edges are not equal along some axis, they must be on a common,
    equally spaced grid.
    (this is meant for combining histograms run with the same float or
    fixed array axbins options)
    '''
    if len(h1.shape) != len(h2.shape):
        raise ValueError('Can only add histograms of the same shape')
    if not (np.all(np.array(h1.shape) == np.array([len(e) - 1 for e in e1]))
            and 
            np.all(np.array(h2.shape) == np.array([len(e) - 1 for e in e2]))
           ):
        raise ValueError('Histogram shape does not match edges')
       
    # iterate over edges, determine overlaps
    p1 = []
    p2 = []
    es = []

    for ei in range(len(e1)):
        e1t = np.array(e1[ei])
        e2t = np.array(e2[ei])
        p1t = [None, None]
        p2t = [None, None]
        
        # if the arrays happen to be equal, it's easy
        if len(e1t) == len(e2t):
            if np.allclose(e1t, e2t, rtol=rtol, atol=atol):
                p1t = [0, 0]
                p2t = [0, 0]
                es.append(0.5 * (e1t + e2t))
                p1.append(p1t)
                p2.append(p2t)
                continue
        
        # if not, things get messy fast. Assume equal spacing (check) 
        s1t = np.diff(e1t)
        s2t = np.diff(e2t)
        if not (np.allclose(s1t[0][np.newaxis], s1t) 
                and np.allclose(s2t[0][np.newaxis], s2t)):
            msg = ('Cannot deal with unequally spaced arrays'
                   f' that do not match (axis {ei})')
            raise RuntimeError(msg)
        if not np.isclose(np.average(s1t), np.average(s2t), 
                          atol=atol, rtol=rtol):
            msg = f'Cannot deal with differently spaced arrays (axis {ei})'
            raise RuntimeError(msg) 
        st = 0.5 * (np.average(s1t) + np.average(s2t))
        if st <= 0.:
            msg = f'Cannot deal with decreasing array values (axis {ei})'
            raise RuntimeError(msg)
        # check if the arrays share a zero point for their scales
        if not np.isclose(((e1t[0] - e2t[0]) / st + 0.5) % 1 - 0.5, 0., 
                          atol=atol, rtol=rtol):
            msg = f'Cannot deal with arrays not on a common grid (axis {ei})'
            raise RuntimeError(msg)

        g0 = 0.5 * ((e1t[0] / st + 0.5) % 1. - 0.5
                    + (e2t[0] / st + 0.5) % 1. - 0.5)        
        # calulate indices of the array endpoints on the common grid (zero point is g0)
        e1i0 = int(np.floor((e1t[0] - g0) / st + 0.5))
        e1i1 = int(np.floor((e1t[-1] - g0) / st + 0.5))
        e2i0 = int(np.floor((e2t[0] - g0) / st + 0.5))
        e2i1 = int(np.floor((e2t[-1] - g0) / st + 0.5))
        
        # set histogram padding based on grid indices
        p1t = [None, None]
        p2t = [None, None]
        if e1i0 > e2i0:
            p1t[0] = e1i0 - e2i0
            p2t[0] = 0
        else:
            p1t[0] = 0
            p2t[0] = e2i0 - e1i0
        if e1i1 > e2i1:
            p1t[1] = 0
            p2t[1] = e1i1 - e2i1
        else:
            p1t[1] = e2i1 - e1i1
            p2t[1] = 0
        # set up new edges based on the grid, initially
        esi0 = min(e1i0, e2i0)
        esi1 = max(e1i1, e2i1)
        est = np.arange(g0 + esi0 * st, g0 + (esi1 + 0.5) * st, st)
        # overwrite with old edges 
        # (2, then 1, to give preference to the histogram 1 edges)
        # meant to avoid accumulating round-off errors through st, g0
        est[e2i0 - esi0: e2i1 + 1 - esi0] = e2t
        est[e1i0 - esi0: e1i1 + 1 - esi0] = e1t
        
        p1.append(p1t)
        p2.append(p2t)
        es.append(est)
        
    h1 = np.pad(h1, mode='constant', constant_values=0, pad_width=p1)
    h2 = np.pad(h2, mode='constant', constant_values=0, pad_width=p2)
    if add:
        hs = h1 + h2
        return hs, es
    else:
        return h1, h2, es

def periodic_sel(array, edges, period):
    '''
    select the elements of array that are between the edges, with both
    in units periodic in period
    
    input:
    ------
    array: float array 
        contains the values to select from
    edges: size 2 indexable of floats
        the lower and upper bounds to select
    period: float
        the period of the range
    
    output:
    -------
    boolean array indicating which elements fall within the range: 
    (edges[0] <= value < edges[1])
    '''
    array = np.array(array)
    period = float(period)
    array %= period
    edges = np.array(edges)
    edges %= period
    
    if edges[0] <= edges[1]:
        out = array >= edges[0]
        out &= array < edges[1]
    else:
        out = array >= edges[0]
        out |= array < edges[1]
    return out

def pad_periodic(pos, margin, period, additional=None):
    '''
    add perioidc repetitions to the pos arrays for values within margin
    of the edges, as well as to additional
    
    input:
    ------
    pos: float array, shape (number of dimesions, number of points)
        positions
    margin:  float, values 0 -- period 
        distance from the edges to include in repetitions
    period:  float
        the periodic of the volume
    additional: other arrays for which to duplicate values (without 
             adding/subtracting the box size); second dimension sould
             match pos
    '''
    
    pos = np.array(pos)
    pos %= period
    ndims = pos.shape[0]
    if additional is not None:
        additional = np.array(additional)
    
    for i in range(ndims):
        e0 = pos[i] < margin
        e1 = pos[i] > period - margin
        diff = np.zeros((ndims, 1))
        diff[i, :] = period
        a2 = pos[:, e1] - diff
        pos = np.append(pos, pos[:, e0] + diff, axis=1)
        pos = np.append(pos, a2, axis=1)
        if additional is not None:
            _a2 = additional[:, e1]
            additional = np.append(additional, additional[:, e0], axis=1)
            additional = np.append(additional, _a2, axis=1)
    
    if additional is None:
        out = pos
    else:
        out = (pos, additional)
    return out

# John Helly's routine, via Peter
def match(arr1, arr2, arr2_sorted=False, arr2_index=None):
    """
    For each element in arr1 return the index of the element with the
    same value in arr2, or -1 if there is no element with the same 
    value. Setting arr2_sorted=True will save some time if arr2 is 
    already sorted into ascending order.

    A precomputed sorting index for arr2 can be supplied using the
    arr2_index parameter. This can save time if the routine is called
    repeatedly with the same arr2 but arr2 is not already sorted.

    It is assumed that each element in arr1 only occurs once in arr2.
    """

    # Workaround for a numpy bug (<=1.4): ensure arrays are native endian
    # because searchsorted ignores endian flag
    if not(arr1.dtype.isnative):
        arr1_n = np.asarray(arr1, dtype=arr1.dtype.newbyteorder("="))
    else:
        arr1_n = arr1
    if not(arr2.dtype.isnative):
        arr2_n = np.asarray(arr2, dtype=arr2.dtype.newbyteorder("="))
    else:
        arr2_n = arr2

    # Sort arr2 into ascending order if necessary
    tmp1 = arr1_n
    if arr2_sorted:
        tmp2 = arr2_n
        idx = slice(0,len(arr2_n))
    else:
        if arr2_index is None:
            idx = np.argsort(arr2_n)
            tmp2 = arr2_n[idx]
        else:
            # Use supplied sorting index
            idx = arr2_index
            tmp2 = arr2_n[arr2_index]

    # Find where elements of arr1 are in arr2
    ptr  = np.searchsorted(tmp2, tmp1)

    # Make sure all elements in ptr are valid indexes into tmp2
    # (any out of range entries won't match so they'll get set to -1
    # in the next bit)
    ptr[ptr>=len(tmp2)] = 0
    ptr[ptr<0]          = 0

    # Return -1 where no match is found
    ind  = tmp2[ptr] != tmp1
    ptr[ind] = -1

    # Put ptr back into original order
    ind = np.arange(len(arr2_n))[idx]
    ptr = np.where(ptr>= 0, ind[ptr], -1)
    
    return ptr

def peakfinder1d(array3d, axisedges, axis=2,
                 minsepdecr=0.5, mintotval=10**12.5,
                 maxcomponents=4, 
                 componentedge_meanfrac=0.05):
    '''
    find 1d peaks in array3d along axis
    (intended for absorption peaks in p-p-v histograms)

    returns 3 3d arrays, with peak properties
    replacing weights in the selected axis. The inds array
    has and extra dimension (last) for first/last inds.

    Returns:
    --------
    cens, tots, inds
    cens: component centers (weighted mean, axisedges units)
    tots: total weight in component
    inds: indices along axis marking each component. 
          start and stop attributes for 'slice', i.e.
          first is the first included bin and second is
          one more than the last included bin.
    NaN values mean the numbered component was not found.
    '''
    dims = list(range(len(array3d.shape)))
    dims.remove(axis)
    itershape = [array3d.shape[i] for i in dims]
    numdims = len(dims)
    lend1 = itershape[0]
    lend2 = itershape[1]
    pos1 = dims[0]
    pos2 = dims[1]
    outshape = list(array3d.shape)
    outshape[axis] = maxcomponents
    outcen = np.zeros(tuple(outshape), dtype=np.float32) / 0.
    outtot = np.zeros(tuple(outshape), dtype=np.float32) / 0.
    outind = np.zeros(tuple(outshape) + (2,), dtype=np.float32) / 0.
    for i1 in range(lend1):
        for i2 in range(lend2):
            subsel = [None for _ in range(numdims)]
            subsel[pos1] = i1
            subsel[pos2] = i2
            searcharr = array3d[tuple(subsel)]
            if np.sum(searcharr) < mintotval:
                continue
            ## find starting peak candidates (noisy)
            #peak: 1st derivative changes sign pos -> neg
            # e.g.: 
            # searcharr: [5,  5,  3, 2, 4,  6, 5,  6,  1, 0, 6]
            # delta1:    [0, -2, -1, 2, 2, -1, 1, -5, -1, 6]
            # delta1 > 0:[T,  F,  F, T, T,  F, T,  F,  F, T]
            # wherepos:  [0, 3, 4, 6, 9]
            # wherejump: [0, 2, 3]
            # ci_s:      [1, 5, 7]
            delta1 = np.diff(searcharr)
            wherepos = np.where(delta1 >= 0.)[0]
            wherejump = np.where(np.diff(wherepos) > 1.1)[0]
            ci_s = wherepos[wherejump] + 1
            # if last index is a peak, cannot find a jump after it.
            if delta1[-1] > 0:
                ci_s = np.append(ci_s, len(delta1))
            ## merge initial peaks if insufficient decrement between
            ## them


