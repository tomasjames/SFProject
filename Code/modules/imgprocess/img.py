################################################################################
############################ 4th Year Project Code #############################
##############################  Functions to run  ##############################
################################################################################

############################# Import statements ################################
# Import RADMC-3D
import radmc3dPy

# Import OS manipulation tools
import os

# Import time tools
import time

# Import units
from radmc3dPy.natconst import *

# Import numpy
import numpy as np

# Import csv
import csv

# Import matplotlib
import matplotlib as mpl
mpl.use('Qt4Agg')
from matplotlib.pyplot import *
from matplotlib import gridspec

# Import astropy
from astropy.io import fits

# Import scipy
import scipy.interpolate
import scipy.ndimage

# Import glob for file scraping
import glob

# Import astrodendro
from astrodendro import Dendrogram

import random

######################### Define function to create maps of N and T #################

def mapMaker(data_type,region,B,dataBins,chiBins,imgwidth):

    '''
    '''

    ############################### Read in image data ##############################

    # Define a string that contains the correct folder name for us in the below statement
    beta_str = str('B=')+str(B)

    d2g = 100

    if data_type == 'radmc':
        # Read in both data types
        inp_data = np.loadtxt(('../../../curvefitting/cloud/')+str(beta_str)+str('/datafeed.txt'))
        chi_data = np.loadtxt(('../../../curvefitting/cloud/')+str(beta_str)+str('/chi_results.txt'),skiprows=1)

    elif data_type == 'arepo':
        # Read in both data types
        inp_data = np.loadtxt(('../../../curvefitting/sphdata/')+str(beta_str)+str('/datafeed.txt'))
        chi_data = np.loadtxt(('../../../curvefitting/sphdata/')+str(beta_str)+str('/chi_results.txt'),skiprows=1)

    elif data_type == 'herschel_snaps':

    # Define the first folder
    first_folder = 'herschel_snaps'
    second_folder = region

    # Read in both data types
    inp_data = np.loadtxt(('/export/home/c1158976/Code/simulations/curvefitting/herschel_snaps/{}/{}/datafeed.txt').format(region,beta_str))
    chi_data = np.loadtxt(('/export/home/c1158976/Code/simulations/curvefitting/herschel_snaps/{}/{}/chi_results.txt').format(region,beta_str),skiprows=1)

    # Also read in dust density information
    dust_density = np.loadtxt(('/export/home/c1158976/Code/simulations/herschel_snaps_psf/{}/B=2.0/blue/dust_density.inp').format(region), skiprows=3)
    dust_temperature = np.loadtxt(('/export/home/c1158976/Code/simulations/herschel_snaps_psf/{}/B=2.0/blue/dust_temperature.dat').format(region), skiprows=3)
    imag = fits.open(('/export/home/c1158976/Code/simulations/herschel_snaps_psf/{}/B=2.0/blue/blue_common_convolved.fits').format(region))

    # Determine data dimensions
    nx = len(imag[0].data[0])
    ny = len(imag[0].data[:,0])

    xpix, ypix, zpix = np.arange(0,nx), np.arange(0,ny), np.arange(0,(len(dust_density)/(nx*ny)))

    # Split data types into plottable quantities
    inp_N = inp_data[:,1]*d2g
    inp_T = inp_data[:,2]

    chi_N, chi_N_error = chi_data[:,1]*d2g, chi_data[:,4]*d2g
    chi_T, chi_T_error = chi_data[:,2], chi_data[:,5]

    # Reshape the data such that x and y pixels correspond
    N_chi_inp, N_chi_inp_error = np.reshape(chi_N, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data)))), np.reshape(chi_N_error, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data))))
    T_chi_inp, T_chi_inp_error = np.reshape(chi_T, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data)))), np.reshape(chi_T_error, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data))))

    N_data_inp = np.reshape(inp_N, (np.sqrt(len(inp_data)),np.sqrt(len(inp_data))))
    T_data_inp = np.reshape(inp_T, (np.sqrt(len(inp_data)),np.sqrt(len(inp_data))))

    # Determine data dimensions
    x_for_ticks = np.linspace(-nx/2,nx/2,6)
    y_for_ticks = np.linspace(-ny/2,ny/2,6)

    ############################ Determine basic image stats #########################

    # Find the lowest and largest N and T for normalising colorbars
    min_N_data, max_N_data = np.amin(np.log10(N_data_inp)), np.amax(np.log10(N_data_inp))
    min_T_data, max_T_data = np.amin(T_data_inp), np.amax(T_data_inp)

    min_N_chi, max_N_chi = np.amin(np.log10(N_chi_inp)), np.amax(np.log10(N_chi_inp))
    min_T_chi, max_T_chi = np.amin(T_chi_inp), np.amax(T_chi_inp)

    # Determine which of the min and max quantities is greatest
    if max(max_N_data, max_N_chi) == max_N_data:
        max_N = max_N_data
    elif max(max_N_data, max_N_chi) == max_N_chi:
        max_N = max_N_chi

    if min(min_N_data, min_N_chi) == min_N_data:
        min_N = min_N_data
    elif min(min_N_data, min_N_chi) == min_N_chi:
        min_N = min_N_chi

    if max(max_T_data, max_T_chi) == max_T_data:
        max_T = max_T_data
    elif max(max_T_data, max_T_chi) == max_T_chi:
        max_T = max_T_chi

    if min(min_T_data, min_T_chi) == min_T_data:
        min_T = min_T_data
    elif min(min_T_data, min_T_chi) == min_T_chi:
        min_T = min_T_chi

    #################################### Plot N PDFs ###################################

    # Determine number of bins
    bins_chi_N = np.linspace(np.log10(np.amin(chi_N)), np.log10(np.amax(chi_N)), chiBins)
    h_chi_N,b_chi_N,p_chi_N = hist(np.log10(chi_N),bins=bins_chi_N,histtype='step')
    close()

    # Determine a 'start' and 'end' point as histogram does not plot continuous curve to 0
    start_y_chi_N = [min(np.log10(h_chi_N)/np.log10(np.sum(h_chi_N))), np.log10(h_chi_N[0])/np.log10(np.sum(h_chi_N))]
    start_x_chi_N = [[b_chi_N[:-1][0]], [b_chi_N[:-1][0]]]

    end_y_chi_N = [min(np.log10(h_chi_N)/np.log10(np.sum(h_chi_N))), np.log10(h_chi_N[-1])/np.log10(np.sum(h_chi_N))]
    end_x_chi_N = [[b_chi_N[:-1][-1]], [b_chi_N[:-1][-1]]]

    # Determine number of bins
    bins_data_N = np.linspace(np.log10(np.amin(inp_N)), np.log10(np.amax(inp_N)), dataBins)
    h_inp_N,b_inp_N,p_inp_N = hist(np.log10(inp_N),bins=bins_data_N,histtype='step')
    close()

    # Determine a 'start' and 'end' point as histogram does not plot continuous curve to 0
    start_y_data_N = [min(np.log10(h_inp_N)/np.log10(np.sum(h_inp_N))), np.log10(h_inp_N[0])/np.log10(np.sum(h_inp_N))]
    start_x_data_N = [[b_inp_N[:-1][0]], [b_inp_N[:-1][0]]]

    end_y_data_N = [min(np.log10(h_inp_N)/np.log10(np.sum(h_inp_N))), np.log10(h_inp_N[-1])/np.log10(np.sum(h_inp_N))]
    end_x_data_N = [[b_inp_N[:-1][-1]], [b_inp_N[:-1][-1]]]

    # Plot
    fig, ax1 = subplots(figsize=(12,5))
    # Plot the chi squared N PDF
    ax1.plot(b_chi_N[:-1],np.log10(h_chi_N)/np.log10(np.sum(h_chi_N)),ls='steps',linewidth=1.0,color='blue',linestyle='--',label=str('$\chi^{2}$'))
    ax1.plot(start_x_chi_N,start_y_chi_N,linewidth=1.0,color='blue',linestyle='--')
    ax1.plot(end_x_chi_N,end_y_chi_N,linewidth=1.0,color='blue',linestyle='--')
    legend(loc='upper right')
    ax1.set_xlabel('$log_{10}(N)$',fontsize=12)
    # Change the colour of the y tick labels
    ax1.set_ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$',fontsize=12,color='blue')
    for tl in ax1.get_yticklabels():
        tl.set_color('blue')

    ax2 = ax1.twinx()
    # Plot the data N PDF
    ax2.plot(b_inp_N[:-1],np.log10(h_inp_N)/np.log10(np.sum(h_inp_N)),ls='steps',linewidth=1.0,color='red',linestyle='-',label=str('Data Derived'))
    ax2.plot(start_x_data_N,start_y_data_N,linewidth=1.0,color='red',linestyle='-')
    ax2.plot(end_x_data_N,end_y_data_N,linewidth=1.0,color='red',linestyle='-')
    legend(loc='lower right')
    # Change the colour of the y tick labels
    ax2.set_ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$\n\n',fontsize=12,color='red',rotation=270,labelpad=25)
    for tl in ax2.get_yticklabels():
        tl.set_color('red')

    #xlim(min_N,max_N)
    title('PDF of Data Derived N and $\chi^{2}$ Recovered $N$')
    #xlabel('$log_{10}(N)$',fontsize=10)
    #ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$',fontsize=10)
    #tick_params(axis='both',labelsize=8)

    savefig('N_PDF.png', dpi=300)
    close()

    #################################### Plot N PDFs ###################################

    # Determine number of bins
    bins_chi_T = np.linspace(np.amin(chi_T), np.amax(chi_T), chiBins)
    h_chi_T,b_chi_T,p_chi_T = hist(chi_T,bins=bins_chi_T,histtype='step')
    close()

    # Determine a 'start' and 'end' point as histogram does not plot continuous curve to 0
    start_y_chi_T = [min(h_chi_T), h_chi_T[0]]
    start_x_chi_T = [[b_chi_T[:-1][0]], [b_chi_T[:-1][0]]]

    end_y_chi_T = [min(h_chi_T), h_chi_T[-1]]
    end_x_chi_T = [[b_chi_T[:-1][-1]], [b_chi_T[:-1][-1]]]

    # Determine number of bins
    bins_data_T = np.linspace(np.amin(inp_T), np.amax(inp_T), dataBins)
    h_inp_T,b_inp_T,p_inp_T = hist(inp_T,bins=bins_data_T,histtype='step')
    close()

    # Determine a 'start' and 'end' point as histogram does not plot continuous curve to 0
    start_y_data_T = [min(h_inp_T), h_inp_T[0]]
    start_x_data_T = [[b_inp_T[:-1][0]], [b_inp_T[:-1][0]]]

    end_y_data_T = [min(h_inp_T), h_inp_T[-1]]
    end_x_data_T = [[b_inp_T[:-1][-1]], [b_inp_T[:-1][-1]]]

    # Plot
    fig, ax1 = subplots(figsize=(12,5))
    # Plot the chi squared N PDF
    ax1.plot(b_chi_T[:-1],h_chi_T,ls='steps',linewidth=1.0,color='blue',linestyle='--',label=str('$\chi^{2}$'))
    ax1.plot(start_x_chi_T,start_y_chi_T,linewidth=1.0,color='blue',linestyle='--')
    ax1.plot(end_x_chi_T,end_y_chi_T,linewidth=1.0,color='blue',linestyle='--')
    ax1.set_xlabel('$T\/(K)$',fontsize=12)
    legend(loc='upper left')
    # Change the colour of the y tick labels
    ax1.set_ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$',fontsize=12,color='blue')
    for tl in ax1.get_yticklabels():
        tl.set_color('blue')

    ax2 = ax1.twinx()
    # Plot the data N PDF
    ax2.plot(b_inp_T[:-1],h_inp_T,ls='steps',linewidth=1.0,color='red',linestyle='-',label=str('Data Derived'))
    ax2.plot(start_x_data_T,start_y_data_T,linewidth=1.0,color='red',linestyle='-')
    ax2.plot(end_x_data_T,end_y_data_T,linewidth=1.0,color='red',linestyle='-')
    legend(loc='lower left')
    # Change the colour of the y tick labels
    ax2.set_ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$\n\n',fontsize=12,color='red',rotation=270,labelpad=25)
    for tl in ax2.get_yticklabels():
        tl.set_color('red')

    #xlim(min_N,max_N)
    title('PDF of Data Derived T and $\chi^{2}$ Recovered $T$')
    #xlabel('$log_{10}(N)$',fontsize=10)
    #ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$',fontsize=10)
    #tick_params(axis='both',labelsize=8)

    savefig('T_PDF.png', dpi=300)
    close()

    ################################# Plot image data ################################

    # Format the tick marks
    xMark,yMark = [], []
    for i in range(0,len(x_for_ticks)):
        xMark.append(str(np.round(x_for_ticks[i]*((imgwidth*au/pc)/nx),2)))
        yMark.append(str(np.round(y_for_ticks[i]*((imgwidth*au/pc)/nx),2)))

    # Plot the data along with a PDF
    figure()
    imshow(np.log10(N_chi_inp),origin='lower',aspect='auto',vmin=min_N,vmax=max_N)
    colorbar(label='$log_{10}(N)$')
    xlabel('Width (pc)')
    ylabel('Height (pc)')
    xticks(np.linspace(0,255,6),xMark) # Change the tick labels
    yticks(np.linspace(0,255,6),yMark) # Change the tick labels
    title('A Map of the $\chi^{2}$ Recovered $N$\n')

    savefig('map_N_chi.png', dpi=300)
    close()

    '''
    # Plot error map
    figure()
    subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
    imshow(np.log10(N_chi_inp_error),origin='lower')
    colorbar(label='$log_{10}(\sigma_{N})$')
    xlabel('X (pixels)')
    ylabel('Y (pixels)')
    title('A Map of the $\chi^{2}$ Recovered $N$ Error Values\n')

    subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
    #hist(np.log10(N_chi_inp),bins=5,normed=True)
    h, b = np.histogram(np.log10(N_chi_inp_error), bins=np.linspace(np.log10(np.amin(N_chi_inp_error)), np.log10(np.amax(N_chi_inp_error)), len(np.unique(N_chi_inp_error))), density=False)
    #bar(b[:-1],h,width=(max(h)-min(h))/len(N_chi_inp))
    bar(b[:-1],h/np.sum(h),width=0.02)
    title('PDF of $\chi^{2}$ Recovered $N$')
    xlabel('$log_{10}N$')
    ylabel('$P(log_{10}N)$')
    tick_params(axis='both',labelsize=8)

    savefig('map_N_chi_error.png', dpi=300)
    close()
    '''

    # Plot the data along with a PDF

    figure()
    imshow(T_chi_inp,origin='lower',vmin=min_T,vmax=max_T)
    colorbar(label='$T\/(K)$')
    xlabel('Width (pc)')
    ylabel('Height (pc)')
    xticks(np.linspace(0,255,6),xMark) # Change the tick labels
    yticks(np.linspace(0,255,6),yMark) # Change the tick labels
    title('A Map of the $\chi^{2}$ Recovered $T$\n')

    savefig('map_T_chi.png', dpi=300)
    close()

    '''
    # Plot error map
    figure()
    subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
    imshow(T_chi_inp_error ,origin='lower')
    colorbar(label='$T\/(K)$')
    xlabel('X (pixels)')
    ylabel('Y (pixels)')
    title('A Map of the $\chi^{2}$ Recovered $T$\n')

    subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
    #hist(T_chi_inp,bins=5,normed=True)
    h, b = np.histogram(T_chi_inp_error, bins=np.linspace((np.amin(T_chi_inp_error)), (np.amax(T_chi_inp_error)), len(np.unique(T_chi_inp_error))), density=False)
    #bar(b[:-1],h,width=(max(h)-min(h))/len(T_chi_inp))
    bar(b[:-1],h/np.sum(h),width=0.06)
    xlim(min(inp_T),max(inp_T))
    title('PDF of $\chi^{2}$ Recovered $T$')
    xlabel(r'$T\/(K)$')
    ylabel(r'$P(T_{\chi^{2}})$')
    tick_params(axis='both',labelsize=8)

    savefig('map_T_chi_error.png', dpi=300)
    close()
    '''

    # Plot the data along with a PDF

    figure()
    N_data_inp[N_data_inp == 0] = np.nan
    imshow(np.log10(N_data_inp),origin='lower',vmin=min_N,vmax=max_N)
    colorbar(label='$log_{10}(N)$')
    xlabel('X (pixels)')
    ylabel('Y (pixels)')
    title('A Map of the Data Input $N$\n')

    savefig('map_N_data.png', dpi=300)
    close()

    # Plot the data along with a PDF

    figure()
    T_data_inp[T_data_inp == 0] = np.nan
    imshow(T_data_inp,origin='lower',vmin=min_T,vmax=max_T)
    colorbar(label='$T\/(K)$')
    xlabel('X (pixels)')
    ylabel('Y (pixels)')
    title('A Map of the Data Input $T$\n')

    savefig('map_T_data.png', dpi=300)
    close()

    ##################### Determine line of sight T variations #######################

    '''
    if data_type == 'radmc':
        dust_density = np.loadtxt(('../../../workingsims_psf/')+str(beta_str)+str('/blue/background_15K/dust_density.inp'), skiprows=3)
        dust_temperature = np.loadtxt(('../../../workingsims_psf/')+str(beta_str)+str('/blue/background_15K/dust_temperature.dat'), skiprows=3)
        imag = fits.open(('../../../workingsims_psf/{}/blue/background_15K/blue_common_convolved.fits').format(beta_str))

    elif data_type == 'arepo':
        dust_density = np.loadtxt(('../../../data_psf/')+str(beta_str)+str('/blue/dust_project/dust_density.inp'), skiprows=3)
        dust_temperature = np.loadtxt(('../../../data_psf/')+str(beta_str)+str('/blue/dust_project/dust_temperature.dat'), skiprows=3)
        imag = fits.open(('../../../data_psf/{}/blue/dust_project/blue_common_convolved.fits').format(beta_str))

    elif data_type == 'herschel_snaps':
        dust_density = np.loadtxt(('/export/home/c1158976/Code/simulations/herschel_snaps_psf/{}/B=2.0/blue/dust_density.inp').format(region), skiprows=3)
        dust_temperature = np.loadtxt(('/export/home/c1158976/Code/simulations/herschel_snaps_psf/{}/B=2.0/blue/dust_temperature.dat').format(region), skiprows=3)
        imag = fits.open(('/export/home/c1158976/Code/simulations/herschel_snaps_psf/{}/B=2.0/blue/blue_common_convolved.fits').format(region))

    # Determine data dimensions
    nx = len(imag[0].data[0])
    ny = len(imag[0].data[:,0])

    xpix, ypix, zpix = np.arange(0,nx), np.arange(0,ny), np.arange(0,(len(dust_density)/(nx*ny)))

    '''

    sigma_T_inp = []

    for i in range(0,len(xpix)*len(ypix)):

        # Finds the locations of every each pixel's Z values
        loc = np.arange(i,len(zpix)*(len(xpix)*len(ypix)),len(xpix)*len(ypix))

        # This is essentially the standard deviated weighted mean of the temperatures
        sigma_T_inp.append((np.sum(((dust_temperature[loc] - inp_T[i])**2)*inp_N[i]))/np.sum(inp_N))

    ################################# Save image data ################################
    '''
    # Collate all of the data into one array
    combined = [N_chi_inp, T_chi_inp, N_data_inp, T_data_inp]

    # Define the names of the individual data sets
    combined_names = ['N_chi_inp', 'T_chi_inp', 'N_data_inp', 'T_data_inp']

    for i in range(0,len(combined)):
        # Define a new header file
        hdul = fits.HDUList()

        # Append to a primary header
        hdul.append(fits.PrimaryHDU())

        # Append the data
        hdul.append(fits.ImageHDU(data=combined[i]))

        # Write the data to a fits file with name of the data array
        hdul.writeto(str(combined_names[i])+str('.fits'))
    '''
    ################################# Plot contours ################################

    figure(1)
    plot(chi_N, chi_T, 'g.', label=r'$\chi^{2}$ Recovered Values')
    xlabel(r'$N_{\chi^{2}}$')
    ylabel(r'$T_{\chi^{2}}$')
    title(r'A Plot Showing the $\chi^{2}$ Recovered Values of $T$ and $N$')
    legend(loc='best')
    savefig('contours.png',dpi=300)
    close()

    ######################## Plot further graphs for report #########################

    # Plot the recovered temperatures against the input temperatures
    figure(2)
    plot(inp_T, chi_T, 'b.', markersize=2)
    plot(np.linspace(np.min(inp_T),np.max(chi_T),10),np.linspace(np.min(inp_T),np.max(chi_T),10),'r--',label='Linear')
    xlabel(r'$T_{N}\/(K)$')
    ylabel(r'$T_{\chi^{2}}\/(K)$')
    xlim(min(inp_T),max(chi_T)+0.2)
    ylim(min(inp_T),max(chi_T)+0.2)
    title(r'A Graph Comparing the $\chi^{2}$ Recovered $T$ to the data input $T$')
    legend(loc='best')
    savefig('T.png', dpi=300)
    close()

    figure(2)
    plot(np.log10(inp_N), np.log10(chi_N), 'b.', markersize=2)
    plot(np.linspace(np.min(np.log10(inp_N)),np.max(np.log10(inp_N)),10),np.linspace(np.min(np.log10(inp_N)),np.max(np.log10(inp_N)),10),'r--',label='Linear')
    xlabel(r'$log_{10}(N_{data})$')
    ylabel(r'$log_{10}(N_{\chi^{2}})$')
    #xlim(min(inp_T),max(chi_T))
    #ylim(min(inp_T),max(chi_T))
    title(r'A Graph Comparing the $\chi^{2}$ Recovered $N$ to the data input $N$')
    legend(loc='best')
    savefig('N.png', dpi=300)
    close()

    figure(3)
    subplot(2,1,1)
    plot(np.log10(inp_N), sigma_T_inp, 'g.', markersize=2)
    xlabel(r'$log_{10}(N_{input})$')
    ylabel(r'$\sigma^{2}_{T_{N}}$')
    title('A Graph Comparing the Line of Sight Temperature Variations to the Input Column Density')

    subplot(2,1,2)
    plot(inp_T, sigma_T_inp, 'g.', markersize=2)
    xlabel(r'$T_{input}\/(K)$')
    ylabel(r'$\sigma^{2}_{T_{N}}$')
    title('A Graph Comparing the Line of Sight Temperature Variations to the Input Temperature')

    savefig('sigma_T_inp.png', dpi=300)
    close()

    figure(5)
    plot(np.log10(inp_N), (inp_T/chi_T), 'b.', markersize=2)
    xlabel(r'$log_{10}(N_{input})$')
    ylabel(r'$\frac{T_{input}}{T_{\chi^{2}}}$')
    title('A Graph Comparing the Input Column Density to the Temperature Ratio')
    savefig('T_ratio_inp.png', dpi=300)
    close()

    figure()
    plot(inp_T, (inp_N/chi_N), 'b.', markersize=2)
    xlabel(r'$T_{input}\/(K)$')
    ylabel(r'$\frac{N_{input}}{N_{\chi^{2}}}$')
    title('A Graph Comparing the Input Column Density to the Temperature Ratio')
    savefig('N_ratio_inp.png', dpi=300)
    close()

    figure(6)
    plot(np.log10(chi_N), (inp_T/chi_T), 'b.', markersize=2)
    xlabel(r'$N_{\chi^{2}}$')
    ylabel(r'$\frac{T_{input}}{T_{\chi^{2}}}$')
    title(r'A Graph Comparing the $\chi^{2}$ Column Density to the Temperature Ratio')
    savefig('T_ratio_chi.png', dpi=300)
    close()

    '''
    ############### Plot better comparisons for the RADMC-3D images #################

    # Read in all data
    blue_raw = np.loadtxt('../../workingsims/blue/background_15K/image_trans.out',skiprows=6)
    green_raw = np.loadtxt('../../workingsims/green/background_15K/image_trans.out',skiprows=6)
    plw_raw = np.loadtxt('../../workingsims/plw/background_15K/image_trans.out',skiprows=6)
    pmw_raw = np.loadtxt('../../workingsims/pmw/background_15K/image_trans.out',skiprows=6)
    psw_raw = np.loadtxt('../../workingsims/psw/background_15K/image_trans.out',skiprows=6)
    red_raw = np.loadtxt('../../workingsims/red/background_15K/image_trans.out',skiprows=6)

    eff = [68.92474, 97.90361, 496.10677, 346.71804, 247.12451, 153.94392]

    # Reformat the data
    blue = np.reshape(blue_raw, (np.sqrt(len(blue_raw)),np.sqrt(len(blue_raw))))
    green = np.reshape(green_raw, (np.sqrt(len(green_raw)),np.sqrt(len(green_raw))))
    plw = np.reshape(plw_raw, (np.sqrt(len(plw_raw)),np.sqrt(len(plw_raw))))
    pmw = np.reshape(pmw_raw, (np.sqrt(len(pmw_raw)),np.sqrt(len(pmw_raw))))
    psw = np.reshape(psw_raw, (np.sqrt(len(psw_raw)),np.sqrt(len(psw_raw))))
    red = np.reshape(red_raw, (np.sqrt(len(red_raw)),np.sqrt(len(red_raw))))

    store = [blue, green, plw, pmw, psw, red]
    names = ['blue', 'green', 'plw', 'pmw', 'psw', 'red']

    # Find max and min of each data and put into list
    max_store = [np.amax(blue), np.amax(green), np.amax(plw), np.amax(pmw), np.amax(psw), np.amax(red)]
    min_store = [np.amin(blue), np.amin(green), np.amin(plw), np.amin(pmw), np.amin(psw), np.amin(red)]

    # Plot each data type and normalise the max and min values for comparison
    for s in range(0,len(store)):
        imshow(store[s],origin='left',vmin=np.min(min_store),vmax=np.max(max_store))
        colorbar(label=r'Intensity $(erg/s/cm/cm/Hz/ster$')
        xlabel('X (pixels)')
        ylabel('Y (pixels)')
        title(str(names[s])+str(': $\lambda=$')+str(eff[s])+str('$\mu m$'))
        savefig(str(names[s])+str('.png'))
        close()
    '''

############################# Define function to create dendrograms ###################

def dendrogram():

    '''
    '''

    ############################# Read in the data ################################

    count, data, filenames = 0, [], []

    # Loop through all fits files
    for file in glob.glob('*.fits'):
        count += 1

        # Read in using astropy and then append the data to a list for use later
        contents = fits.open(str(file))[1]

        '''
        # Append the name of the file to the position (count+len(data)) followed by
        # the file contents
        data.insert(count+len(data), file)
        data.append(contents)
        filenames.append(file)
        '''

        # Define a dendrogram instance
        d = Dendrogram.compute(contents.data, verbose=True)

        # Let the dendrogram become an interactive plot
        p = d.plotter()

        # Add the data to a figure
        fig = figure()

        # Define subplots to plot with
        ax1 = subplot2grid((6,6), (0,0), colspan=4,rowspan=4)
        ax2 = subplot2grid((6,6), (5,0), colspan=4,rowspan=1)
        #ax13 = subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
        #ax14 = subplot2grid((3, 3), (1, 2), rowspan=2)

        leaves = []

        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Get the mask of the region highlighted
            mask = struct.get_mask()

            # Create FITS HDU objects
            mask_hdu = fits.PrimaryHDU(mask.astype('short'), contents.header)

            ax1.imshow(contents.data, origin='lower', interpolation='nearest')

            # Show contour for ``min_value``
            p.plot_contour(ax1, color='black')
            #p.colorbar()

            # Determine the type of structure
            if struct.is_leaf:
                leaves.append(struct)

                # Extract the indices of the core and its value
                indices = struct.indices(subtree=True)
                vals = struct.values(subtree=True)

                #print indices
                #print vals

                # Highlight two branches
                p.plot_contour(ax1, structure=struct, lw=3, colors="#%06x" % random.randint(0, 0xFFFFFF))

        # Plot the entire tree in black
        p.plot_tree(ax2, color='black')

        # Add labels to the plot
        if 'T_chi' in file:
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            ax1.set_title('A Map Showing the Structure\n of $T$ for the $\chi^{2}$ Recovered Values')

            ax2.set_xlabel('Structure')
            ax2.set_ylabel('$T\/(K)$')
            ax2.set_title('A Dendrogram Showing the Structure\n of $T$ for the $\chi^{2}$ Recovered Values')
        elif 'N_chi' in file:
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            ax1.set_title('A Map Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')

            ax2.set_xlabel('Structure')
            ax2.set_ylabel('$N\/(g\/cm^{-3})$')
            ax2.set_title('A Dendrogram Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')
        elif 'T_data' in file:
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            ax1.set_title('A Map Showing the Structure\n of $T$ for the Data Input Values')

            ax2.set_xlabel('Structure')
            ax2.set_ylabel('$T\/(K)$')
            ax2.set_title('A Dendrogram Showing the Structure\n of $T$ for the Data Input Values')
        elif 'N_data' in file:
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            ax1.set_title('A Map Showing the Structure\n of $N$ for the Data Input Values')

            ax2.set_xlabel('Structure')
            ax2.set_ylabel('$N\/(g\/cm^{-3})$')
            ax2.set_title('A Dendrogram Showing the Structure\n of $N$ for the Data Input Values')

        #ax1.imshow(contents, origin='lower', interpolation='nearest', cmap=cm.Blues, vmax1=4.0)

        # Plot black contours on to the data
        # Show contour for ``min_value``
        #p.plot_contour(ax1, color='black')

        # Save the image and then close to save memory and allow the next image to take its place
        fig.savefig(str(file)+str('.jpg'), dpi=300)
        close()
