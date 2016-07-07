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
from matplotlib.pyplot import *

# Import astropy
from astropy.io import fits

# Import scipy
import scipy.interpolate
import scipy.ndimage

# Import glob for file scraping
import glob

import random

######################### Define function to create maps of N and T #################

def mapMaker(data_type):

    '''
    '''

    ############################### Read in image data ##############################

    if data_type == 'radmc':
        # Read in both data types
        inp_data = np.loadtxt('../../curvefitting/cloud/datafeed.txt')
        chi_data = np.loadtxt('../../curvefitting/cloud/chi_coarse.txt',skiprows=1)

    elif data_type == 'arepo':
        # Read in both data types
        inp_data = np.loadtxt('../../curvefitting/sphdata/datafeed.txt')
        chi_data = np.loadtxt('../../curvefitting/sphdata/chi_coarse.txt',skiprows=1)

    # Split data types into plottable quantities
    inp_N = inp_data[:,1]
    inp_T = inp_data[:,2]

    chi_N, chi_N_error = chi_data[:,1], chi_data[:,4]
    chi_T, chi_T_error = chi_data[:,2], chi_data[:,5]

    ################################# Plot image data ################################

    # Reshape the data such that x and y pixels correspond
    N_chi_inp, N_chi_inp_error = np.reshape(chi_N, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data)))), np.reshape(chi_N_error, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data))))
    T_chi_inp, T_chi_inp_error = np.reshape(chi_T, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data)))), np.reshape(chi_T_error, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data))))

    N_data_inp = np.reshape(inp_N, (np.sqrt(len(inp_data)),np.sqrt(len(inp_data))))
    T_data_inp = np.reshape(inp_T, (np.sqrt(len(inp_data)),np.sqrt(len(inp_data))))

    # Plot the data along with a PDF
    figure()
    subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
    imshow(np.log10(N_chi_inp),origin='lower',vmin=np.log10(min(chi_N)),vmax=np.log10(max(chi_N)))
    colorbar(label='$log_{10}(N)$')
    xlabel('X (pixels)')
    ylabel('Y (pixels)')
    title('A Map of the $\chi^{2}$ Recovered $N$\n')

    subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
    #hist(np.log10(N_chi_inp),bins=5,normed=True)
    h, b = np.histogram(np.log10(N_chi_inp), bins=np.linspace(np.log10(np.amin(N_chi_inp)), np.log10(np.amax(N_chi_inp)), 100), density=False)
    #bar(b[:-1],h,width=(max(h)-min(h))/len(N_chi_inp))
    bar(b[:-1],np.log10(h),width=0.02)
    xlim(np.log10(min(chi_N)),np.log10(max(inp_N)))
    title('PDF of $\chi^{2}$ Recovered $N$')
    xlabel('$log_{10}(N)$',fontsize=10)
    ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$',fontsize=10)
    tick_params(axis='both',labelsize=8)

    savefig('map_N_chi.png', dpi=300)
    close()

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
    h, b = np.histogram(np.log10(N_chi_inp_error), bins=np.linspace(np.log10(np.amin(N_chi_inp_error)), np.log10(np.amax(N_chi_inp_error)), 100), density=False)
    #bar(b[:-1],h,width=(max(h)-min(h))/len(N_chi_inp))
    bar(b[:-1],h/np.sum(h),width=0.02)
    title('PDF of $\chi^{2}$ Recovered $N$')
    xlabel('$log_{10}N$')
    ylabel('$P(log_{10}N)$')
    tick_params(axis='both',labelsize=8)

    savefig('map_N_chi_error.png', dpi=300)
    close()

    # Plot the data along with a PDF
    figure()
    subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
    imshow(T_chi_inp,origin='lower',vmin=min(inp_T),vmax=max(chi_T))
    colorbar(label='$T\/(K)$')
    xlabel('X (pixels)')
    ylabel('Y (pixels)')
    title('A Map of the $\chi^{2}$ Recovered $T$\n')

    subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
    #hist(T_chi_inp,bins=5,normed=True)
    h, b = np.histogram(T_chi_inp, bins=np.linspace((np.amin(T_chi_inp)), (np.amax(T_chi_inp)), 100), density=False)
    #bar(b[:-1],h,width=(max(h)-min(h))/len(T_chi_inp))
    bar(b[:-1],h,width=0.1)
    xlim(min(inp_T),max(chi_T))
    title('PDF of $\chi^{2}$ Recovered $T$')
    xlabel('$T\/(K)$',fontsize=10)
    ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$',fontsize=10)
    tick_params(axis='both',labelsize=8)

    savefig('map_T_chi.png', dpi=300)
    close()

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
    h, b = np.histogram(T_chi_inp_error, bins=np.linspace((np.amin(T_chi_inp_error)), (np.amax(T_chi_inp_error)), 100), density=False)
    #bar(b[:-1],h,width=(max(h)-min(h))/len(T_chi_inp))
    bar(b[:-1],h/np.sum(h),width=0.06)
    xlim(min(inp_T),max(inp_T))
    title('PDF of $\chi^{2}$ Recovered $T$')
    xlabel(r'$T\/(K)$')
    ylabel(r'$P(T_{\chi^{2}})$')
    tick_params(axis='both',labelsize=8)

    savefig('map_T_chi_error.png', dpi=300)
    close()

    # Plot the data along with a PDF
    figure()
    N_data_inp[N_data_inp == 0] = np.nan
    subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
    imshow(np.log10(N_data_inp),origin='lower',vmin=np.log10(np.amin(N_chi_inp)),vmax=np.log10(np.amax(N_chi_inp)))
    colorbar(label='$log_{10}(N)$')
    xlabel('X (pixels)')
    ylabel('Y (pixels)')
    title('A Map of the Data Input $N$\n')

    subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
    #hist(np.log10(N_data_inp),bins=5,normed=True)
    h, b = np.histogram(np.log10(N_data_inp), bins=np.linspace(np.log10(np.amin(N_data_inp)), np.log10(np.amax(N_data_inp)), 100), density=False)
    #bar(b[:-1],h,width=(max(h)-min(h))/len(N_data_inp))
    bar(b[:-1],np.log10(h),width=0.02)
    xlim(np.log10(min(chi_N)),np.log10(max(inp_N)))
    title('PDF of Data Input $N$')
    xlabel('$log_{10}(N)$',fontsize=10)
    ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$',fontsize=10)
    tick_params(axis='both',labelsize=8)

    savefig('map_N_data.png', dpi=300)
    close()

    # Plot the data along with a PDF
    figure()
    T_data_inp[T_data_inp == 0] = np.nan
    subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
    imshow(T_data_inp,origin='lower',vmin=min(inp_T),vmax=max(chi_T))
    colorbar(label='$T\/(K)$')
    xlabel('X (pixels)')
    ylabel('Y (pixels)')
    title('A Map of the Data Input $T$\n')

    subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
    #hist(T_data_inp,bins=5,normed=True)
    h, b = np.histogram(T_data_inp, bins=np.linspace((np.amin(T_data_inp)), (np.amax(T_data_inp)), 100), density=False)
    #bar(b[:-1],h,width=(max(h)-min(h))/len(T_data_inp))
    bar(b[:-1],h,width=0.1)
    xlim(min(inp_T),max(inp_T))
    title('PDF of Data Input $T$')
    xlabel('$T\/(K)$',fontsize=10)
    ylabel('$log_{10}{(n_{pix})}\/(per\/bin)$',fontsize=10)
    tick_params(axis='both',labelsize=8)

    savefig('map_T_data.png', dpi=300)
    close()

    ##################### Determine line of sight T variations #######################

    dust_density = np.loadtxt('../../workingsims/blue/background_15K/dust_density.inp', skiprows=3)
    dust_temperature = np.loadtxt('../../workingsims/blue/background_15K/dust_temperature.dat', skiprows=3)

    imag = radmc3dPy.image.readImage('../../workingsims/blue/background_15K/image.out')
    xpix, ypix, zpix = np.arange(0,imag.nx), np.arange(0,imag.ny), np.arange(0,(len(dust_density)/(imag.nx*imag.ny)))

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
    title('A Graph Comparing the Line of Sight \n Temperature Variations to the Input Column Density')

    subplot(2,1,2)
    plot(inp_T, sigma_T_inp, 'g.', markersize=2)
    xlabel(r'$T_{input}\/(K)$')
    ylabel(r'$\sigma^{2}_{T_{N}}$')
    title('A Graph Comparing the Line of Sight \n Temperature Variations to the Input Temperature')
    tight_layout()
    savefig('sigma_T_inp.png', dpi=300)
    close()

    figure(5)
    plot(np.log10(inp_N), (inp_T/chi_T), 'b.', markersize=2)
    xlabel(r'$log_{10}(N_{input})$')
    ylabel(r'$\frac{T_{input}}{T_{\chi^{2}}}$')
    title('A Graph Comparing the Input Column Density to the Temperature Ratio')
    savefig('T_ratio_inp.png', dpi=300)
    close()

    figure(5)
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
            p.colorbar()

            # Determine the type of structure
            if struct.is_leaf:
                leaves.append(struct)

                # Extract the indices of the core and its value
                indices = struct.indices(subtree=True)
                vals = struct.values(subtree=True)

                print indices
                print vals

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
