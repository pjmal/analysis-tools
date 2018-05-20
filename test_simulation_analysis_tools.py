#!/usr/bin/env python
"""Test analysis tools for astrophysical simulations"""
import os
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import yt
from yt.units import cm

from simulation_analysis_tools import *


##############################
# example of the library usage
#TS = [200,300,400]       # time steps
TS = [400]       # time steps
#LS = [6,7,8,9,10]     # resolution level (from 1 forward towards higher resolutions)
LS = [6,7,10]    # resolution level (from 1 forward towards higher resolutions)
#AXES = [0,1,2]  # orientation axes for column density
AXES = [1]       # orientation axes for column density
cden_bins = 100  # bins for column density histograms

basedir = './'    # assume simulation data is in subdirectories below the working directory
database = 'SILCC_hdf5_plt_cnt_'   # base name for the simulations

mci = 'MC2'          # nametag for the cloud
xmin = -0.2e20 * cm  # zoomin region edge coordinates [cm]
xmax =  3.0e20 * cm
ymin =  4.7e20 * cm
ymax =  7.4e20 * cm
zmin = -1.2e20 * cm
zmax =  1.0e20 * cm

create_densityfits = 1       # 1: create density fits file, 0: read existing file
create_cdensfits = 1         # 1: create column density file, 0: read existing file
change_to_radmc_direction = 1  # 1: change the column density maps into same orientation as radmc3d maps

cdens_bins = 100

Ncontours = [21,22,23]  # log10 contours for column density, eg. 21 means the contour is drawn at 10**21
Nccolors = ['r','k','m']   # contour colors

for L in LS:  #simulation resolution level
    print(L)
    for T in TS:  #time step
        print(T)
        # read cell size from the simulation data
        simulation = '%sL%d/%s%.4d' % (basedir,L,database,T)
        print(simulation)
        pixcm = smallest_cell_size(simulation)
        print('pixcm: %s' % (pixcm))

        densityfile = 'density_L%d_%s_TS%d.fits' % (L,mci,T)
        if(create_densityfits):
            # create density cube
            dens = make_zoomin_unigrid_densitycube(simulation,pixcm,L,xmin,xmax,ymin,ymax,zmin,zmax)
            # write the data cube into a fits file
            save_fits(dens,densityfile)
        else:
            dens = read_fits_data(densityfile)

        cdensfile_base = 'colden_L%d_%s_TS%d' % (L,mci,T)
        for ax in AXES:
            cdensfile = '%s_ax%d.fits' % (cdensfile_base,ax)
            if(create_cdensfits):
                # create column density fits file
                cdens = make_columndensity(dens,pixcm,ax,change_to_radmc_direction)
                save_fits(cdens,cdensfile)
            else:
                cdens = read_fits_data(cdensfile)  # not needed here, shown as an example

            # make some example plots
            # normalized and not normalized histograms
            titletxt = 'NPDF_L%d_%s_TS%d_ax%d' % (L,mci,T,ax)
            imagename = 'NPDF_L%d_%s_TS%d_ax%d.png' % (L,mci,T,ax)
            plot_cdens_pdf_fit(cdensfile,titletxt,imagename,cdens_bins,norm=False)
            imagename = 'NPDF_L%d_%s_TS%d_ax%d_norm.png' % (L,mci,T,ax)
            plot_cdens_pdf_fit(cdensfile,titletxt,imagename,cdens_bins,norm=True)

            # add verticals lines to mark interesting changing points on the PDF
            imagename = 'NPDF_L%d_%s_TS%d_ax%d_norm_contours.png' % (L,mci,T,ax)
            plot_cdens_pdf_fit(cdensfile,titletxt,imagename,cdens_bins,norm=True,verticals=Ncontours,vcolors=Nccolors)

            # simple map
            titletxt = 'Nmap_L%d_%s_TS%d_ax%d' % (L,mci,T,ax)
            imagename = 'Nmap_L%d_%s_TS%d_ax%d.png' % (L,mci,T,ax)
            plot_cdens_map(cdensfile,titletxt,imagename)

            # add contours on map to show the regions where changes occur on the PDFs
            imagename = 'Nmap_L%d_%s_TS%d_ax%d_contours.png' % (L,mci,T,ax)
            cbartxt = 'Log N(H$_2$) [1/cm$^2$]'
            plot_map(cdensfile,titletxt,cbartxt,imagename,log10scale=True,contours=Ncontours,ccolors=Nccolors)

            # fit a line to the powerlaw part of the PDF and calculate slope
            titletxt = 'NPDF_L%d_%s_TS%d_ax%d' % (L,mci,T,ax)
            imagename = 'NPDF_L%d_%s_TS%d_ax%d_norm_fit.png' % (L,mci,T,ax)
            plot_cdens_pdf_fit(cdensfile,titletxt,imagename,cdens_bins,norm=True,fit=True,function=func_powerlaw_log,xmin_fit=21,xmax_fit=22.5)

