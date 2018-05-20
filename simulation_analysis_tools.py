#!/usr/bin/env python
"""Analysis tools for astrophysical simulations"""
import os
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import yt
from yt.units import cm

MH = 1.6737236e-24     # hydrogen atom mass [g]
MU_H2 = 2.8            # molecular weight of H2, mu_H2


def func_powerlaw(x,m,c,a):
    """Powerlaw function."""
    return a * x**m + c


def func_powerlaw_log(x, m, c):
    """Powerlaw function in logarithmic scale, ie. a straight line."""
    return m*x + c


def save_fits(data, filename):
    """Save tabular data to fits file. """
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(filename,clobber=True)  # overwrite file
    hdulist.close()

    
def read_fits_data(filename):
    """Read existing fits file and return the data."""
    fitsfile = fits.open(filename)
    data = fitsfile[0].data
    return data


def smallest_cell_size(simulation):
    """Return the smallest cell size of the given simulation (for levels 6-)."""
    ds = yt.load(simulation)
    alldata = ds.all_data()   # read all data
    pixcm = np.min(alldata['dx'].in_units('cm').v)  # .v gives the value without units
    return pixcm
    

def make_zoomin_unigrid_densitycube(simulation,pixcm,flashlevel,xmin,xmax,ymin,ymax,zmin,zmax):
    """Make a uniform density grid of the chosen zoomin region. 

    Convert Flash HDF5 hierarchical simulation data into a uniform grid using yt. Select the zoomin region
    and save the density cube as a fits file.
    simulation: name of the simulation data cube (including path)
    pixcm: smallest cell size [cm]
    flashlevel: resolution level for the unigrid (starting from L1, growing towards higher resolutions), yt level is flashlevel-1
    xmin,xmax,ymin,ymax,zmin,zmax: edges of the zoomin region [cm]
    """
    # calculate the dimensions of the chosen zoomin region
    xp = int(np.ceil((xmax-xmin)/pixcm))
    yp = int(np.ceil((ymax-ymin)/pixcm))
    zp = int(np.ceil((zmax-zmin)/pixcm))

    # make a unigrid of the zoomin region
    lev = flashlevel-1  # in yt the levels are counted up from 0 -> yt level is flashlevel-1 (eg. L10 -> yt level 9)
    ds = yt.load(simulation)
    unigrid = ds.covering_grid(level=lev,left_edge=[xmin,ymin,zmin],dims=[xp, yp, zp])
    dens = unigrid['density']
    return dens


def make_columndensity(density,pixcm,ax,change_to_radmc_direction=1):
    """Calculate 2D column density along the given axis of the density cube.

    Assume density is given in g/cm**3 units. Convert column density to number of H2-molecules per cm**2 units. By default the resulting maps are rotated to the same direction as radmc3d gives for synthetic observation maps.
    """
    print('calculate column density for axis %s' % ax)
    cdens = np.sum(density, ax) * pixcm/(MU_H2*MH)
    if(change_to_radmc_direction):
        #same direction as the synthetic observations results given by radmc3d
        cdens = np.rot90(cdens.copy())
        cdens = np.flipud(cdens.copy())
    return cdens


def plot_cdens_map(cdensfile,titletxt,imagename):
    """Plot column density fits file map in logarithmic scale."""
    plt.clf()
    cdens = read_fits_data(cdensfile)
    plt.imshow(np.log10(cdens),origin='lower')
    plt.title('%s' % titletxt)
    cbar = plt.colorbar()
    cbar.set_label('Log N(H$_2$) [1/cm$^2$]')
    plt.savefig(imagename)


def plot_map(mapfile,titletxt,cbartxt,imagename,log10scale=False,contours=[],ccolors=[]):
    """Plot fits file map with given title and colorbar label, and option for logarithmic or linear scale."""
    plt.clf()
    data = read_fits_data(mapfile)
    if(log10scale):
        plt.imshow(np.log10(data),origin='lower')
    else:
        plt.imshow(data,origin='lower')
    plt.title('%s' % titletxt)
    cbar = plt.colorbar()
    cbar.set_label(cbartxt)
    if((len(contours) > 0) and (len(contours) == len(ccolors))):
        if(log10scale):
            CS = plt.contour(np.log10(data),contours,origin='lower',colors=ccolors,linewidths=0.9)
        else:
            CS = plt.contour(data,contours,origin='lower',colors=ccolors,linewidths=0.9)
    plt.savefig(imagename)


def plot_cdens_pdf_fit(cdensfile,titletxt,imagename,cdens_bins,norm=False,verticals=[],vcolors=[],fit=False,function=None,xmin_fit=None,xmax_fit=None):
    """Plot column density distribution as a normalized probability density function or a histogram in log10 scale.

    Parameters:
    cdensfile: column density fits file
    titletxt: title text
    imagename: output image name
    cdens_bins: histogram bin size
    norm: normalized histogram [True], default: False
    verticals: array of x-axis coordinates for marking vertical lines, default: []
    vcolors: array of colors for marking vertical lines, required to be the same length as verticals, default: []
    fit: True: fit function 'function' to the histogram in the range [xmin_fit,xmax_fit], assumes all necessary parameters are given, default: False
    function: name of the function used in fitting, default: None
    xmin_fit: lower limit for fitting, given as an exponent of 10 (ie. the limit value is 10**xmin_fit), default: None
    xmax_fit: upper limit for fitting, given as an exponent of 10 (ie. the limit value is 10**xmax_fit), default: None
    """
    plt.clf()
    cdens = read_fits_data(cdensfile)
    # hist returns: n: values of the histogram bins, bins: edges of the bins, patches: silent list of individual patches used to create the histogram
    if(norm):
        n, bins, patches = plt.hist(np.log10(np.ravel(cdens)),bins=cdens_bins,normed=True,histtype='step',log=True)
        plt.ylabel('Normalized pixel count', size=18)
    else:
        n, bins, patches = plt.hist(np.log10(np.ravel(cdens)),bins=cdens_bins,normed=False,histtype='step',log=True)
        plt.ylabel('Pixel count', size=18)
    plt.xlabel('Log N(H$_2$) [1/cm$^2$]', size=18)
    plt.title('%s' % titletxt)
    if((len(verticals) > 0) and (len(verticals) == len(vcolors))):
        for vv in range(len(verticals)):
            plt.plot([verticals[vv],verticals[vv]],[10**(-5),10**(-0.1)],lw=1.8,color=vcolors[vv])
            plt.text(verticals[vv], 10**(-0.9), '%s' % verticals[vv], size=14)
    if(fit):
        xx = bins[0:-1]+0.5*np.diff(bins)
        m1 = np.nonzero((xx>xmin_fit)&(xx<xmax_fit) & (n > 0))
        X = xx[m1]
        Y = np.log10(n[m1])
        popt, pcov = curve_fit(func_powerlaw_log, X, Y)
        print(popt)
        Yexp = func_powerlaw_log(X, *popt)
        Yfit = np.power(10,Yexp)     # use same scaling as the histogram
        plt.plot(X, Yfit, 'k--', lw=2)
        if(norm):
            plt.text(xmin_fit+0.5, 10**(-1), 'k = %.2f' % popt[0], size=14)
        else:
            plt.text(xmin_fit+0.5, 10**(4), 'k = %.2f' % popt[0], size=14)
    plt.savefig(imagename)

