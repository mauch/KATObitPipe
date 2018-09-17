"""Convert fits file to .jpeg file with useful contrast scaling options.
Optional input is contrast percentage.
"""

import pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import os, sys
from optparse import OptionParser
import numpy as np
import math

plt.switch_backend('Agg')

MAX_REJECT = 0.5 
MIN_NPIXELS = 5 
GOOD_PIXEL = 0 
BAD_PIXEL = 1 
KREJ = 2.5 
MAX_ITERATIONS = 5 

def zscale (image, nsamples=1000, contrast=0.5, bpmask=None, zmask=None): 
    """Implement IRAF zscale algorithm 
    nsamples=1000 and contrast=0.25 are the IRAF display task defaults 
    bpmask and zmask not implemented yet 
    image is a 2-d numpy array 
    returns (z1, z2) 
    """ 
    # Sample the image 
    samples = zsc_sample (image, nsamples, bpmask, zmask) 
    npix = len(samples) 
    samples.sort() 
    zmin = samples[0] 
    zmax = samples[-1] 
    # For a zero-indexed array 
    center_pixel = (npix - 1) // 2 
    if npix%2 == 1: 
        median = samples[center_pixel] 
    else: 
        median = 0.5 * (samples[center_pixel] + samples[center_pixel + 1]) 
 
    # 
    # Fit a line to the sorted array of samples 
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT)) 
    ngrow = max (1, int (npix * 0.01)) 
    ngoodpix, zstart, zslope = zsc_fit_line (samples, npix, KREJ, ngrow, 
                                             MAX_ITERATIONS) 
 
    if ngoodpix < minpix: 
        z1 = zmin 
        z2 = zmax 
    else: 
        if contrast > 0: zslope = zslope / contrast 
        z1 = max (zmin, median - (center_pixel - 1) * zslope) 
        z2 = min (zmax, median + (npix - center_pixel) * zslope) 
    return z1, z2 
 
def zsc_sample (image, maxpix, bpmask=None, zmask=None): 
     
    # Figure out which pixels to use for the zscale algorithm 
    # Returns the 1-d array samples 
    # Don't worry about the bad pixel mask or zmask for the moment 
    # Sample in a square grid, and return the first maxpix in the sample 
    nc = image.shape[0] 
    nl = image.shape[1] 
    stride = max (1.0, math.sqrt((nc - 1) * (nl - 1) / float(maxpix))) 
    stride = int (stride) 
    samples = image[::stride,::stride].flatten() 
    return samples[:maxpix] 
     
def zsc_fit_line (samples, npix, krej, ngrow, maxiter): 
 
    # 
    # First re-map indices from -1.0 to 1.0 
    xscale = 2.0 / (npix - 1) 
    xnorm = np.arange(npix) 
    xnorm = xnorm * xscale - 1.0 
 
    ngoodpix = npix 
    minpix = max (MIN_NPIXELS, int (npix*MAX_REJECT)) 
    last_ngoodpix = npix + 1 
 
    # This is the mask used in k-sigma clipping.  0 is good, 1 is bad 
    badpix = np.zeros(npix, dtype="int32") 
 
    # 
    #  Iterate 
 
    for niter in range(maxiter): 
 
        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix): 
            break 
         
        # Accumulate sums to calculate straight line fit 
        goodpixels = np.where(badpix == GOOD_PIXEL) 
        sumx = xnorm[goodpixels].sum() 
        sumxx = (xnorm[goodpixels]*xnorm[goodpixels]).sum() 
        sumxy = (xnorm[goodpixels]*samples[goodpixels]).sum() 
        sumy = samples[goodpixels].sum() 
        sum = len(goodpixels[0]) 
 
        delta = sum * sumxx - sumx * sumx 
        # Slope and intercept 
        intercept = (sumxx * sumy - sumx * sumxy) / delta 
        slope = (sum * sumxy - sumx * sumy) / delta 
         
        # Subtract fitted line from the data array 
        fitted = xnorm*slope + intercept 
        flat = samples - fitted 
 
        # Compute the k-sigma rejection threshold 
        ngoodpix, mean, sigma = zsc_compute_sigma (flat, badpix, npix) 
 
        threshold = sigma * krej 
 
        # Detect and reject pixels further than k*sigma from the fitted line 
        lcut = -threshold 
        hcut = threshold 
        below = np.where(flat < lcut) 
        above = np.where(flat > hcut) 
 
        badpix[below] = BAD_PIXEL 
        badpix[above] = BAD_PIXEL 
         
        # Convolve with a kernel of length ngrow 
        kernel = np.ones(ngrow,dtype="int32") 
        badpix = np.convolve(badpix, kernel, mode='same') 
 
        ngoodpix = len(np.where(badpix == GOOD_PIXEL)[0]) 
         
        niter += 1 
 
    # Transform the line coefficients back to the X range [0:npix-1] 
    zstart = intercept - slope 
    zslope = slope * xscale 
 
    return ngoodpix, zstart, zslope 
 
def zsc_compute_sigma (flat, badpix, npix): 
 
    # Compute the rms deviation from the mean of a flattened array. 
    # Ignore rejected pixels 
 
    # Accumulate sum and sum of squares 
    goodpixels = np.where(badpix == GOOD_PIXEL) 
    sumz = flat[goodpixels].sum() 
    sumsq = (flat[goodpixels]*flat[goodpixels]).sum() 
    ngoodpix = len(goodpixels[0]) 
    if ngoodpix == 0: 
        mean = None 
        sigma = None 
    elif ngoodpix == 1: 
        mean = sumz 
        sigma = None 
    else: 
        mean = sumz / ngoodpix 
        temp = sumsq / (ngoodpix - 1) - sumz*sumz / (ngoodpix * (ngoodpix - 1)) 
        if temp < 0: 
            sigma = 0.0 
        else: 
            sigma = math.sqrt (temp) 
 
    return ngoodpix, mean, sigma 

def get_background_variance(data,sigma_clip=5.0,tolerance=0.01):
    """Compute the variance by iteratively removing outliers greater than a given sigma
    until the mean changes by no more than tolerance.

    Inputs
    ------
    data - 1d numpy array of data to compute variance
    sigma_clip - the amount of sigma to clip the data before the next iteration
    tolerance - the fractional change in the mean to stop iterating

    Outputs
    -------
    variance - the final background variance in the sigma clipped image
    """

    #Initialise diff and data_clip and mean and std
    diff = 1
    mean = data.mean()
    data_clip = data
    while diff > tolerance:
        data_clip = data_clip[np.abs(data_clip)<mean+sigma_clip*data_clip.std()]
        newmean = data_clip.mean()
        diff = np.abs(mean-newmean)/(mean+newmean)
        mean = newmean
    return np.var(data_clip)


def writejpeg(data, filename, contrast=0.05, cmap='gray'):
    """Write a 2d array of data to a jpeg with contrast scaling.

    Inputs
    ------
    data - 2d array of data values
    filename - name of output jpeg file (.jpeg will be appended to the name)
    contrast - float between 0 and 100 that selects the fraction of the pixel distrubution to scale into the colormap
    cmap - the desired colourmap to be input to matplotlib
    """

    # Get grid size
    # Make a grid for the x,y coordinate values (just an integer grid)- again this should be reworked one day for a full WCS treatment
    y,x = np.mgrid[slice(1, data.shape[0]+1, 1), slice(1, data.shape[1]+1, 1)]
    # Get the values fot the contrast scaling        
    #datasort = np.sort(data[np.isfinite(data)])         #Remone NaNs often found in fits images.
    #percentile = contrast/100.0
    #arrayremove = int(len(datasort)*(1.0 - percentile)/2.0)
    #lowcut,highcut = datasort[arrayremove],datasort[-(arrayremove+1)]
    lowcut,highcut = zscale(data,nsamples=100000,contrast=contrast)
    # make a lookup table for colour scaling (todo- change the colour scaling from linear to other functions))
    levels = np.linspace(lowcut,highcut, num=150)
    # Colormap selection (todo- make the colormap a kwarg)
    cmapin = plt.get_cmap(cmap)
    norm = BoundaryNorm(levels, ncolors=cmapin.N, clip=True)
    im=plt.figure(figsize=(5,5))
    plt.axes([0,0,1,1])
    plt.pcolormesh(x,y,data, cmap=cmapin, norm=norm)
    #Axis size
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.axis('off')
    plt.savefig(filename+'.jpeg')
    plt.close(im)



def fits2jpeg(fitsfilename,contrast=99.9,cmap='jet',chans=None,imchans=False,forceaverage=False,weightaverage=False,area=1.0):
    """Convert FITS files to jpegs using matplotlib.

    Inputs
    ------
    fitsfilenames - String
                    The name of ths input fits file.
    contrast - float between 0 and 1 that selects the fraction of the pixel distrubution to scale into the colormap
    cmap - the desired colourmap to be input to matplotlib
    chans - the desired channel range to use
    imchans - produce an individual jpeg for each channel in the input file
    forceaverage - produce an average of the range of channels in chans
    weightaverage - weight the averages
    area - fraction of image (centered) to display
    """

    #Only work if the image exists
    if not os.path.exists(fitsfilename):
        print 'Specified fits file does not exist'
        sys.exit(-1)
    outname = os.path.splitext(fitsfilename)[0]
    #open file in pyfits. ANd get the image plane(s) and the header
    datahdu = pyfits.open(fitsfilename)
    imageheader = datahdu[0].header
    allimagedata = datahdu[0].data[0]
    #make a masked array to remove nans
    allimagedata = np.ma.masked_array(allimagedata, np.isnan(allimagedata))
    #Cut selected area
    xpixels = int(0.5*allimagedata.shape[1]*(1-area))
    ypixels = int(0.5*allimagedata.shape[2]*(1-area))
    allimagedata=allimagedata[:,xpixels:allimagedata.shape[1]-xpixels,ypixels:allimagedata.shape[2]-ypixels]
    #cut out desired area
    chan_range = chans
    if not chan_range: 
        imagedata=allimagedata
        chan_range='1,'+str(imagedata.shape[0])
    chan_range = str(chan_range).split(',')
    # Get the desired subset of the fits file to converty to jpeg
    if len(chan_range)==1: imagedata = allimagedata[int(chan_range[0])-1:int(chan_range[0])]
    else: imagedata = allimagedata[int(chan_range[0])-1:int(chan_range[1])-1]    
    # This will work on pipeline images- but needs to be reworked to work on any image you want
    # Loop over every channel and produce a jpeg for each one if the user asks
    if len(chan_range)>1:
    	if imchans==True or imageheader['CTYPE3'] != 'FREQ':
 	       	[writejpeg(imageplane,outname+'_'+str(num+int(chan_range[0])),float(contrast),cmap) for num,imageplane in enumerate(imagedata)]             
    if len(chan_range)==1:
    	if imchans==True or imageheader['CTYPE3'] != 'FREQ':
    		writejpeg(imagedata[0],outname,float(contrast),cmap)                 
    # Do a weighted or straight average image only if image cube or forced average 
    if forceaverage==True or imageheader['CTYPE3'] == 'FREQ':
        # Set a dummy weight array if not doing weights
        weightarray = np.ones(imagedata.shape[0])
        # Recompute weights if the user askes for variance weights
        if weightaverage==True : weightarray = [get_background_variance(imageplane.flatten()) for imageplane in imagedata]
        #Compute the average
        avdata = np.average(imagedata,axis=0,weights=weightarray)
        #Write out the averaged image
        writejpeg(avdata,outname,contrast,cmap)
    datahdu.close()

