#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:33:56 2019

@author: wijers

Contains a number of general utility functions for making plots
"""
import numpy as np
import matplotlib as mpl
import mpl_toolkits.axes_grid1 as axgrid
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.legend_handler as mlh
import matplotlib.patheffects as mppe 

import fire_an.makeplots.tol_colors as tc

# default
fontsize = 12

def handleinfedges(hist, setmin=-100., setmax=100.):
    for ei in range(len(hist['edges'])):
        if hist['edges'][ei][0] == -np.inf:
            hist['edges'][ei][0] = setmin
        if hist['edges'][ei][-1] == np.inf:
            hist['edges'][ei][-1] = setmax

def handleinfedges_dct(edges, setmin=-100., setmax=100.):
    for ei in edges.keys():
        if edges[ei][0] == -np.inf:
            edges[ei][0] = setmin
        if edges[ei][-1] == np.inf:
            edges[ei][-1] = setmax

### small plot helpers
def setticks(ax, fontsize, color='black', labelbottom=True, top=True,
             labelleft=True, labelright=False, right=True, labeltop=False,
             left=True, bottom=True):
    ax.minorticks_on()
    ax.tick_params(labelsize=fontsize, direction='in', right=right, top=top,
                   left=left, bottom=bottom,
                   axis='both', which='both', color=color,
                   labelleft=labelleft, labeltop=labeltop,
                   labelbottom=labelbottom, labelright=labelright)


### functions to make actual plots
def add_colorbar(ax, img=None, vmin=None, vmax=None, cmap=None, clabel=None,
                 newax=False, extend='neither', fontsize=fontsize,
                 orientation='vertical'):
    if img is None:
        cmap = mpl.cm.get_cmap(cmap)
        norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
        cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap,norm=norm,extend=extend,orientation=orientation)
    elif newax:
        div = axgrid.make_axes_locatable(ax)
        cax = div.append_axes("right", size="5%", pad=0.2)
        cbar = mpl.colorbar.Colorbar(cax,img,extend=extend)
    else:
        cbar = mpl.colorbar.Colorbar(ax,img,extend=extend,orientation=orientation)
    ax.tick_params(labelsize=fontsize-2)
    if clabel is not None:
        cbar.set_label(clabel,fontsize=fontsize)
    return cbar

def add_2dplot(ax, bins, edges, toplotaxes, log=True, usepcolor=False, pixdens=False, shiftx=0., shifty=0., **kwargs):
    # hist3d can be a histogram of any number >=2 of dimensions
    # like in plot1d, get the number of axes from the length of the edges array
    # usepcolor: if edges arrays are not equally spaced, imshow will get the ticks wrong
    summedaxes = tuple(list( set(range(len(edges)))-set(toplotaxes) )) # get axes to sum over
    toplotaxes= list(toplotaxes)
    #toplotaxes.sort()
    axis1, axis2 = tuple(toplotaxes)
    # sum over non-plotted axes
    if len(summedaxes) == 0:
        imgtoplot = bins
    else:
        imgtoplot = np.sum(bins, axis=summedaxes)
    
    
    if pixdens:
        numdims = 2 # 2 axes not already summed over 
        binsizes = [np.diff(edges[toplotaxes[0]]), np.diff(edges[toplotaxes[1]]) ] # if bins are log, the log sizes are used and the enclosed log density is minimised
        baseinds = list((np.newaxis,)*numdims)
        normmatrix = np.prod([(binsizes[ind])[tuple(baseinds[:ind] +\
                              [slice(None,None,None)] +\
                              baseinds[ind+1:])] for ind in range(numdims)])
        if axis1 > axis2:
            imgtoplot = imgtoplot.T
        imgtoplot /= normmatrix
        if axis1 > axis2:
            imgtoplot = imgtoplot.T
        del normmatrix
        
    if log:
        imgtoplot = np.log10(imgtoplot)
    # transpose plot if axes not in standard order; normally, need to use transposed array in image
    if axis1 < axis2:
        imgtoplot = imgtoplot.T
    if usepcolor:
        _kwargs = kwargs.copy()
        if 'rasterized' not in _kwargs:
            _kwargs['rasterized'] = True           
        img = ax.pcolormesh(edges[axis1] + shiftx, edges[axis2] + shifty, imgtoplot, **_kwargs)
    else:
        img = ax.imshow(imgtoplot, origin='lower', interpolation='nearest',\
                        extent=(edges[axis1][0] + shiftx, edges[axis1][-1] + shiftx,\
                                edges[axis2][0] + shifty, edges[axis2][-1] + shifty),\
                        **kwargs)
    if 'vmin' in kwargs.keys():
        vmin = kwargs['vmin']
    else:
        vmin = np.min(imgtoplot[np.isfinite(imgtoplot)])
    if 'vmax' in kwargs.keys():
        vmax = kwargs['vmax']
    else:
        vmax = np.max(imgtoplot[np.isfinite(imgtoplot)])
    return img, vmin, vmax

def add_2dhist_contours(ax, bins, edges, toplotaxes,
                        mins=None, maxs=None, histlegend=True, 
                        fraclevels=True, levels=None, legend=True, 
                        dimlabels=None, legendlabel=None,
                        legendlabel_pre=None, shiftx=0., shifty=0., 
                        dimshifts=None, **kwargs):
    '''
    colors, linestyles: through kwargs
    othersmin and othersmax should be indices along the corresponding histogram axes
    assumes xlim, ylim are already set as desired
    dimlabels can be used to override (long) dimension labels from hist
    '''
    # get axes to sum over; preserve order of other axes to match limits
    
        
    summedaxes = list(range(len(edges)))
    summedaxes.remove(toplotaxes[0])
    summedaxes.remove(toplotaxes[1])
    
    #print('min/max per edge: %s'%str([(mins[i], maxs[i], len(edges[i])) for i in [0,1,2]]))
    #print('mins: %s, maxs: %s'%(mins, maxs))
    
    if dimlabels is None:
        dimlabels = [''] * len(edges)    
    if mins is None:
        mins= (None,)*len(edges)
    if maxs is None:
        maxs = (None,)*len(edges)
    if dimshifts is None:
        dimshifts = (0.,) * len(edges)
	
    # get the selection of min/maxs and apply the selection, put axes in the desired order
    sels = [slice(mins[i], maxs[i], None) for i in range(len(edges))]
    sels = tuple(sels)
    
    if len(summedaxes) > 0:
        binsum = np.sum(bins[sels], axis=tuple(summedaxes))
    else:
        binsum = bins[sels]
    if toplotaxes[0] > toplotaxes[1]:
        binsum = binsum.transpose()
    #print('min/max binsum: %.4e, %.4e'%(np.min(binsum),np.max(binsum)))
    
    binfrac = np.sum(binsum) / np.sum(bins) # fraction of total bins selected
    # add min < dimension_quantity < max in legend label
    if legendlabel is None:
        labelparts = [r'%.1f $<$ %s $<$ %.1f, '%(edges[i][mins[i]] + dimshifts[i], dimlabels[i], edges[i][maxs[i]] + dimshifts[i]) if (mins[i] is not None and maxs[i] is not None) else\
                      r'%.1f $<$ %s, '%(edges[i][mins[i]] + dimshifts[i], dimlabels[i])                  if (mins[i] is not None and maxs[i] is None)     else\
		              r'%s $<$ %.1f, '%(dimlabels[i], edges[i][maxs[i]] + dimshifts[i])                  if (mins[i] is None and maxs[i] is not None)     else\
		              '' for i in range(len(edges))] #no label if there is no selection on that dimension
        legendlabel = ''.join(labelparts)
        # add percentage of total histogram selected
        if legendlabel[-2:] == ', ':
            legendlabel = legendlabel[:-2] + ': '
        legendlabel += '%.1f%%'%(100.*binfrac) 

    if legendlabel_pre is not None:
        legendlabel = legendlabel_pre + legendlabel
    
    #xlim = ax.get_xlim()
    #ylim = ax.get_ylim()
    if levels is None:
        if fraclevels:
            levels = [1., 0.9, 0.5] # enclosed fractions for each level (approximate)
        else:
            levels = [1e-3,3e-2,0.1,0.5]

    if fraclevels: # assumes all levels are between 0 and 1
        binsum = binsum/np.sum(binsum) # redo normalisation for the smaller dataset
        #print('min/max binsum: %.4e, %.4e'%(np.min(binsum),np.max(binsum)))
        
        # for sorting, normialise bins by bin size: peak finding depends on density, should not favour larger bins
        numdims = 2 # 2 axes not already summed over 
        binsizes = [np.diff(edges[toplotaxes[0]]), np.diff(edges[toplotaxes[1]]) ] # if bins are log, the log sizes are used and the enclosed log density is minimised
        baseinds = list((np.newaxis,)*numdims)
        normmatrix = np.prod([(binsizes[ind])[tuple(baseinds[:ind] + [slice(None,None,None)] + baseinds[ind+1:])] for ind in range(numdims)])

        binsumcopy = binsum.copy() # copy to rework
        bindens    = binsumcopy/normmatrix
        bindensflat= bindens.copy().reshape(np.prod(bindens.shape)) # reshape creates views; argsorting later will mess up the array we need for the plot
        binsumcopy = binsumcopy.reshape(np.prod(binsumcopy.shape))
        binsumcopy = binsumcopy[np.argsort(bindensflat)] # get all histogram values in order of histogram density (low to high)
        
        binsumcopy = np.flipud(binsumcopy) # flip to high-to-low
        cumul = np.cumsum(binsumcopy) # add values high-to-low 
        wherelist = [[(np.where(cumul<=level))[0],(np.where(cumul>=level))[0]] for level in levels] # list of max-lower and min-higher indices

        ### made for using bin counts -> binsumcopy is ordered y its own values
	    # sort out list: where arrays may be empty -> levels outside 0,1 range, probabaly
	    # set value level 0 for level == 1. -> just want everything (may have no cumulative values that large due to fp errors)
	    # if all cumulative values are too high (maxmimum bin has too high a fraction), set to first cumulative value (=max bin value)
	    # otherwise: interpolate values, or use overlap
	    #print binsumcopy, binsumcopy.shape
	    #print cumul
	    #print wherelist
	    #return levels, cumul, binsumcopy, wherelist
        if np.all(normmatrix == normmatrix[0,0]): # all bins are the same size
            valslist = [cumul[0]  if  wherelist[i][0].shape == (0,) else\
	                    0.        if (wherelist[i][1].shape == (0,) or levels[i] == 1) else\
		                np.interp([levels[i]], np.array([      cumul[wherelist[i][0][-1]],      cumul[wherelist[i][1][0]] ]),\
                                               np.array([ binsumcopy[wherelist[i][0][-1]], binsumcopy[wherelist[i][1][0]] ]) )[0]\
		                for i in range(len(levels))]
            pltarr = binsum
        else: # find a reasonable interpolation of bindens in stead; need to plot the contours in binsdens as well, in this case
            bindensflat.sort() # to match cumul array indices: sort, then make high to low
            bindensflat = bindensflat[::-1]
            valslist = [bindensflat[0]  if  wherelist[i][0].shape == (0,) else\
	                    0.        if (wherelist[i][1].shape == (0,) or levels[i] == 1) else\
		                np.interp([levels[i]], np.array([      cumul[wherelist[i][0][-1]],      cumul[wherelist[i][1][0]] ]),\
		                                                 np.array([ bindensflat[wherelist[i][0][-1]], bindensflat[wherelist[i][1][0]] ]))[0]\
		                for i in range(len(levels))]
            pltarr = bindens
        #print('min/max bindens: %.4e, %f, min/max flat: %.4e, %f'%(np.min(bindens),np.max(bindens),np.min(bindensflat),np.max(bindensflat)))
        #print('binsum shape: %s, bindens shape: %s, normmatrix shape: %s,  x: %i, y: %i'%(str(binsum.shape),str(bindens.shape), str(normmatrix.shape), len(hist['edges'][toplotaxes[0]]), len(hist['edges'][toplotaxes[1]])))
        #print('wherelist: %s'%wherelist)
        #plt.subplot(2,1,1)
        #plt.pcolor(hist['edges'][toplotaxes[0]], hist['edges'][toplotaxes[1]], np.log10(binsum.T), vmin = np.min( np.log10(binsum.T)[np.isfinite(np.log10(binsum.T))]))
        #plt.colorbar()
        #plt.subplot(2,1,2)
        #plt.pcolor(hist['edges'][toplotaxes[0]], hist['edges'][toplotaxes[1]], np.log10(bindens.T), vmin = np.min( np.log10(bindens.T)[np.isfinite(np.log10(bindens.T))]))
        #plt.colorbar()
        #plt.show()
        
        del normmatrix
        del binsumcopy
        del binsum
        del bindens
        del bindensflat
        #for i in range(len(levels)):
        #    if not (wherelist[i][0].shape == (0,) or wherelist[i][1].shape == (0,)):
	    #        print('interpolating (%f, %f) <- index %i and (%f, %f)  <- index %i to %f'\
	    #	 %(cumul[wherelist[i][0][-1]],binsumcopy[wherelist[i][0][-1]],wherelist[i][0][-1],\
	    #          cumul[wherelist[i][1][0]], binsumcopy[wherelist[i][1][0]], wherelist[i][1][0],\
    	#	   levels[i]) )
        #print(np.all(np.diff(binsumcopy)>=0.))
        uselevels = np.copy(valslist)
        # check for double values; fudge slightly if levels are the same
        anyequal = np.array([np.array(valslist) == val for val in valslist])
        if np.sum(anyequal) > len(valslist): # any levels equal to a *different* level
            eqvals = [np.where(anyequal[ind])[0] for ind in range(len(valslist))] # what levels any value is equal to
            eqgroups = set([tuple(list(eq)) for eq in eqvals]) # get the sets of unique values
            eqgroups = list(eqgroups)
            fudgeby = 1.e-8
            grouplist = [(np.where(np.array([ind in group for group in eqgroups]))[0])[0] for ind in range(len(valslist))] # which group is each uselevel index in
            groupindlist = [(np.where(ind == np.array(eqgroups[grouplist[ind]]))[0])[0] for ind in range(len(valslist))] # which group index corresponds to a goven uselevel index
            addto = [[valslist[group[0]]*fudgeby*ind for ind in range(len(group))] for group in eqgroups] #add nothing for single-element groups
                            
            valslist = [uselevels[ind] + addto[grouplist[ind]][groupindlist[ind]] for ind in range(len(valslist))]
            print('Desired cumulative fraction levels were %s; using value levels %s fudged from %s'%(levels, valslist, uselevels))
            uselevels = valslist
        else:
            print('Desired cumulative fraction levels were %s; using value levels %s'%(levels,uselevels))
    else:
        uselevels=levels
    
    removezerolevelprops = False
    if len(uselevels) > 1:
        if uselevels[0] == uselevels[1]:
            uselevels = uselevels[1:]
            removezerolevelprops = True
            
    #print binsum, binsum.shape
    if 'linestyles' in kwargs:        
        linestyles = kwargs['linestyles']
    else:
        linestyles = [] # to not break the legend search
    
    if removezerolevelprops: # a duplicate level was kicked out -> remove properties for that level
        if 'linestyles' in kwargs.keys():
            kwargs['linestyles'] = kwargs['linestyles'][1:]
        if 'colors' in kwargs.keys():
            kwargs['colors'] = kwargs['colors'][1:]
            
    # get pixel centres from edges
    centres0 = edges[toplotaxes[0]][:-1] + shiftx + 0.5 * np.diff(edges[toplotaxes[0]]) 
    centres1 = edges[toplotaxes[1]][:-1] + shifty + 0.5 * np.diff(edges[toplotaxes[1]])
    contours = ax.contour(centres0, centres1, pltarr.T, uselevels, **kwargs)
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    # make a legend to avoid crowding plot region
    #for i in range(len(levels)):
    #    contours.collections[i].set_label('%.0e'%levels[i])
    # color only legend; get a solid line in the legend
    
    #ax.tick_params(labelsize=fontsize,axis='both')
    if 'solid' in linestyles:
        i = np.where(np.array(linestyles)=='solid')[0][0]
        contours.collections[i].set_label(legendlabel)
    else: # just do the first one
        contours.collections[0].set_label(legendlabel)
    if histlegend:
        ax.legend(loc='lower right',title=r'$f_{\mathrm{O VII}}, f_H=0.752$')
    
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(np.linspace(minval, maxval, n)))
    return new_cmap    

def paste_cmaps(cmaplist, edges, trunclist=None, transwidths=None):
    if trunclist is None:
        trunclist = [(0., 1.)] * len(cmaplist)
    if transwidths is None:
        transwidths = [0.] * max((len(cmaplist) - 1), 1)
    # cobble together a color map
    nsample = 256
    cmaps = [mpl.cm.get_cmap(cmap) if isinstance(cmap, type('')) else cmap \
             for cmap in cmaplist]
    # the parts of each color bar to use
    cmaps = [truncate_colormap(cmaps[i], minval=trunclist[i][0], maxval=trunclist[i][1]) \
                               for i in range(len(cmaplist))]
    # the parts of the 0., 1. range to map each color bar to
    #edges = [float(ed) for ed in edges]
    vmin = edges[0]
    vmax = edges[-1]
    ivran = 1. / (vmax - vmin)
    edges_split = [[edges[i], 
                    edges[i + 1] - 0.5 * transwidths[i]] if i == 0 else \
                   [edges[i] + 0.5 * transwidths[i - 1], 
                    edges[i + 1]] if i == len(cmaplist) - 1 else \
                   [edges[i] + 0.5 * transwidths[i - 1], 
                    edges[i + 1] - 0.5 * transwidths[i]] \
                   for i in range(len(cmaplist))]
    #print(edges_split)
    
    ranges_mapto = [np.linspace((edges_split[i][0] - vmin) * ivran,\
                                (edges_split[i][1] - vmin) * ivran,\
                                nsample) for i in range(len(cmaplist))] 
    # occasionally off by some tiny fp error, causing range check errors
    ranges_mapto[0][0] = 0.
    ranges_mapto[-1][-1] = 1.
    #print(ranges_mapto)
    range_mapfrom = np.linspace(0., 1., nsample)
    maplist = [(ranges_mapto[ci][i], cmaps[ci](range_mapfrom[i])) \
               for ci in range(len(cmaplist)) for i in range(nsample)]
    name = '_'.join(['{name}-{vmin}-{vmax}'.format(name=cmap.name,
                                                   vmin=edges_split[i][0],
                                                   vmax=edges_split[i][1])\
                     for i, cmap in enumerate(cmaps)])
    #print(name)
    #return maplist
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
         name, maplist)
    cmap.set_under(cmap(0.))
    cmap.set_over(cmap(1.))
    return cmap


class HandlerDashedLines(mlh.HandlerLineCollection):
    """
    Custom Handler for LineCollection instances.
    """
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        # figure out how many lines there are
        numlines = len(orig_handle.get_segments())
        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)
        xdata = np.array(xdata)
        leglines = []
        # divide the vertical space where the lines will go
        # into equal parts based on the number of lines
        ydata = ((height) / (numlines + 1)) * np.ones(xdata.shape, float)
        # for each line, create the line at the proper location
        # and set the dash pattern
        for i in range(numlines):
            legline = mlines.Line2D(xdata, ydata * (numlines - i) - ydescent)
            self.update_prop(legline, orig_handle, legend)
            # set color, dash pattern, and linewidth to that
            # of the lines in linecollection
            try:
                color = orig_handle.get_colors()[i]
            except IndexError:
                color = orig_handle.get_colors()[0]
            try:
                dashes = orig_handle.get_dashes()[i]
            except IndexError:
                dashes = orig_handle.get_dashes()[0]
            try:
                lw = orig_handle.get_linewidths()[i]
            except IndexError:
                lw = orig_handle.get_linewidths()[0]
            if dashes[0] is not None and dashes[1] is not None:
                # seem to come out twice the input size when using dashes[1] -> fix
                legline.set_dashes([_d *0.5 for _d in dashes[1]])
            legline.set_color(color)
            legline.set_transform(trans)
            legline.set_linewidth(lw)
            leglines.append(legline)
        return leglines
    
def get_perc_and_points(xdata, ydata, xbins,\
                        percentiles=(5., 25., 50., 75., 95.),\
                        mincount_x=10,\
                        getoutliers_y=True, getmincounts_x=True,\
                        x_extremes_only=True):
    '''
    for xdata, ydata points, and bins xbins in the x direction:
    get percentiles in each bin, and the points outside the selections:
    - getoutliers_y: in each x bin, return the points with y more extreme than 
      the requested percentiles
    - getmincounts_x: return all the points in x bins with at least as many 
      elements as mincount_x
    - x_extremes_only: if getmincounts_x , only return points from x bins more 
      extreme than the first and last where the threshold is not met
    
    returns:
    --------
    the percentiles in each bin (list of arrays, sorted by increasing
    percentile)
    the outlier points (x, y arrays)
    the indices of points meeting the counts threshold (numpy array of indices)
    the x binning includes points outside the min/max range at the ends of the 
    arrays
    '''
    xdata = np.array(xdata)
    ydata = np.array(ydata)
    
    percentiles = np.sort(percentiles)
    bininds = np.digitize(xdata, xbins)
    nbins = len(xbins) + 1
    binlists = [[xdata[bininds == i], ydata[bininds == i]] \
                for i in range(nbins)]
    percs = np.array([np.percentile(binlists[i][1], percentiles) \
                      if len(binlists[i][1]) > 0 else \
                      np.NaN * np.ones(len(percentiles)) \
                      for i in range(nbins)])
    xmininds = np.where(np.array([len(_l[0]) \
                                  for _l in binlists]) >= mincount_x)[0]
    outliers = [np.array([]), np.array([])]
    if getmincounts_x:
        xsel_out = np.append(np.arange(xmininds[0]),\
                             np.arange(xmininds[-1] + 1, nbins)) \
                   if x_extremes_only else \
                   np.where(np.array([len(_l[0]) \
                                  for _l in binlists]) < mincount_x)[0]
        xout = np.array([x for i in xsel_out for x in binlists[i][0]])
        yout = np.array([y for i in xsel_out for y in binlists[i][1]])
        outliers[0] = np.append(outliers[0], xout)
        outliers[1] = np.append(outliers[1], yout)
    if getoutliers_y:
        if getmincounts_x:
            if x_extremes_only:
                xsel = np.arange(xmininds[0], xmininds[-1] + 1)
            else:
                xsel = xmininds
        else:    
            xsel = np.arange(0, nbins)
        sel_binlists = [np.logical_or(binlists[i][1] < percs[i][0],\
                                      binlists[i][1] > percs[i][-1]) \
                        for i in range(nbins)]
        
        xout = np.array([x for i in xsel for x in binlists[i][0][sel_binlists[i]]])
        yout = np.array([y for i in xsel for y in binlists[i][1][sel_binlists[i]]])
        outliers[0] = np.append(outliers[0], xout)
        outliers[1] = np.append(outliers[1], yout)
        
    return percs, outliers, xmininds
    
def getoutline(linewidth):
    patheff = [mppe.Stroke(linewidth=linewidth + 0.5, foreground="black"),
               #mppe.Stroke(linewidth=linewidth + 0.5, foreground="white"),
               mppe.Normal()]
    return patheff

# specific to redshift stuff, but I use this sort of thing a lot
def getzcolorfunc(zvals, ztol=1e-2):
    nz = len(zvals)
    zvals.sort()
    _colors = tc.tol_cmap('rainbow_discrete', nz)
    zticks = np.linspace(0.5 / nz, 1. - 0.5 / nz, nz)
    colors = _colors(zticks)
    colors = [mcolors.to_rgb(col) for col in colors]
    def getzcolor(zval):
        zmatch = np.where(np.isclose(zval, zvals, rtol=ztol, atol=ztol))[0]
        if len(zmatch) != 1:
            msg = (f'redshift {zval} is too far from any of {zvals}, '
                   f'or too close (within {ztol}) to more than one.')
            raise ValueError(msg)
        color = colors[zmatch[0]]
        return color
    return getzcolor
