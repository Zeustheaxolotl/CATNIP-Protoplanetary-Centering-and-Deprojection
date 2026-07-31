from astropy.io import fits
import numpy as np
from photutils.centroids import centroid_sources
from photutils.centroids import centroid_quadratic
import cv2
import matplotlib.pyplot as plt
from pathlib import Path


def contour_fits(data): 
    ''' This will take the initial fits file, clean it and return a list of contours found in the image.
    Data: is the variable holding your fits file in a 2d array.'''
    #Use the median as the guess for noise in the background
    median = np.nanmedian(data)
    #find the highest value in the data
    max_sig = np.nanmax(data)
    #Calculate signal to noise ratio
    sn_ratio = median/max_sig
    print(sn_ratio, 'Signal to Noise')

    #Set all nans in the array to 0
    data[np.isnan(data) == True] = 0

    #clean the background of your assumed noise, note this is not ideal for large disks. 
    data[data < median] = 0
    
    #These thresholds determine what amount of signal is needed to go into the contour and is when it is high is best for high s/n ratio disks. 
    threshold = np.percentile(data, 99.8)
    threshold_max = np.percentile(data, 100)
    #Calculate the levels of the contour so that there are 10 levels of brightness to choose from between the thresholds. 
    levels = np.linspace(threshold, threshold_max, 10)
    #Find the actual contours that go on these levels. 
    cs = plt.contour(data, levels)
    #print(cs)
    #show how the contours look over the image 
    plt.contour(data, levels= levels)
    plt.imshow(data, origin='lower', cmap='gray')
    plt.show()

    return cs

def fit_ellipse(data, contours): 
    ''' This function takes the data as well as the contours found in contour fits and matches elliipses to them.'''
    #Sort the contours so that only ones that could make ellipses are selected. 
    valid_contours = [
    cnt.astype(np.float32)
    for cnt in contours
    if len(cnt) >= 5
    ]
    #sort the contours by area in descending order. Largest area to smallest. 
    sorted_contours = sorted(
    valid_contours,
    key=cv2.contourArea,
    reverse=True
    )
    #make the empty array to keep ellipses. 
    ellipses = []
   #for each contour make an ellipse and add it to the ellipses. 
    for contour in sorted_contours: 
        ellipse = cv2.fitEllipse(contour)
        ellipses.append(ellipse)
    return ellipses

def photutils(data, ellipse_fit): 
    ''' Check if around the center of your ellipse fit there is a bright centroid source with another center. 
    Return the new center. '''

    xc = ellipse_fit[0][0]
    yc = ellipse_fit[0][1]
    w = ellipse_fit[1][0]
    h = ellipse_fit[1][1]
    #identify the area you are searching for this bright quadratic like bright source. 
    box = np.floor(max(w, h)/5)
    box = box.astype(int)

    if box%2 !=1: 
        #box has to be even for code to run. 
        box = box+1
    #calculate ideal center using centroid_quadratic. 
    x1, y1 = centroid_sources(data,xc, yc, box, centroid_func = centroid_quadratic)
    return x1, y1
        
    
def plot_idea(data, xc, yc, contours = True,  elpse_cnt=0, wmax=0, hmax=0, angle=0, idx =0, mask =None):
    ''' This is a function that handles all necessary plotting for your guesses. It can be plotted with estimations of the circular checks or without. 
    data: array of your fits file 
    xc, yc: proposed center candidates. 
    contours: the variable that chooses to plot contours or not 
    elpse_cnt: how many circles are checked 
    wmax, hmax: the estimated size of the disk
    angle: the estimated angle of the ellipse 
    mask: variable incase any of the pixels should not be included. 
    '''

    dim1 = data.shape[0]
    dim2 = data.shape[1]
    #the part of the plot that holds your contours
    grid = np.full((dim1, dim2), np.nan)
    index = np.indices(data.shape)
    #list of y indices for data
    y = index[0]
    #list of x indices for data
    x = index[1]
    #convert angle into radians
    
    rad_ang = np.radians(angle)
    #plot contours
    if contours == True: 
        if idx == 0: 
            for k in range(elpse_cnt):
                #go through all of the radii. 
                
                #say where you are in comparison to the end of the program
                status = f"{k+1}/{elpse_cnt}" 
                print(status)
                
                #find how the qualities of the ellipses plotted 
                ratio = (k+1)/elpse_cnt
                k_width = (wmax * ratio /2 )**2
                k_height = (hmax * ratio/2)**2
                sina = np.sin(rad_ang)
                cosa = np.cos(rad_ang)
                #calculate the distance for each coordinate from center.
                distx = ((x-xc)*np.cos(rad_ang)+(y-yc)*np.sin(rad_ang))
                disty = (-(x-xc)*np.sin(rad_ang)+(y-yc)*np.cos(rad_ang))
                if mask is not None: 
                    #gotta fix this cus right now if mask.shape != shape data never define idk
                    if mask.shape == data.shape:
                        #index where the dist is less than the ellipse we're checking+1 and greater than the ellipse func = 1
                        idx = np.where(((distx**2/k_width +disty**2/ k_height)<=1) & (distx**2/k_width +disty**2/ k_height > .9) & (mask > 0))
                else: 
                    #index where the radius is less than the radius we're checking+1 and greater than the radius
                    idx = np.where(((distx**2/k_width +disty**2/ k_height)<=1) & (distx**2/k_width +disty**2/ k_height > .9))
                    grid[idx] = 255
        else: 
            grid[idx]=255   
    #show plot for the data and the contours. 
    plt.imshow(data, origin='lower', cmap='gray')
    cmap = plt.get_cmap('viridis')
    cmap.set_bad(color = 'none')
    plt.imshow(grid, origin = 'lower', cmap = cmap)
    plt.plot(xc, yc, 'bo')
    plt.xlim(xc-100,xc+100)
    plt.ylim(yc-100,yc+100)
    plt.show()

def circfit_good(data, xc, yc, rmax, count, median = False, mask =None): 
    #INPUTS
    #data  :  the input image
    # xc: the beginning x center coordinate to branch out from calculated by above func to shorten work
    # yc: the beginning y center coordinate to branch out from calculated by above func to shorten work
    # rmax: the maximum estimated size of the disk. 
    # count: the number of circles evaluated 
    # median: if the standard deviation average should instead be median. 
    #mask: if any pixels need to be masked. 
    # angle: the angle at which the rings are offset as given by above func

    #OUTPUTS
    # grid : grid of results (xr vs yr vs stddev)
   
    #get dimensions of fits image
    dim1 = data.shape[0]
    dim2 = data.shape[1]

    index = np.indices(data.shape)
    #list of y indices for data
    y = index[0]
    #list of x indices for data
    x = index[1]
    grid = []
    r = np.sqrt((x - xc)**2 + (y - yc)**2)
    jump = rmax/count
    radii = np.arange(0, rmax, jump)
    for k in radii:
        #go through all of the radii. 

        #say where you are in comparison to the end of the program
        status = f"{k}/{rmax}" 
        #print(status)
        
        if mask is not None: 
            #gotta fix this cus right now if mask.shape != shape data never define idx
            if mask.shape == data.shape:
                #index where the radius is less than the radius we're checking+1 and greater than the radius
                idx = np.where((r >= k) & (r< k+1) & (mask >0))
        else: 
            #index where the radius is less than the radius we're checking+1 and greater than the radius
            idx = np.where((r >= k) & (r< k+1))
        if (idx[0].size == 0): 
            #if there are no pixels here than skip
            continue
        sd = np.std(data[idx]) #calculate stddev
        #print(sd)
        if not np.isfinite(sd): #make sure that stddev is not NaN or infinite
            sd = 0 # if so it's 0
            
        if abs(np.nanmedian(data[idx]))!=0:
            grid.append((sd/abs(np.nanmedian(data[idx])))**2)
    if median == False: 
        grid = np.nanmean(grid)
    else: 
        grid = np.nanmedian(grid)
    
    return grid


def deproject(data, xc, yc, PA, incl):
    ''' Deproject function that will take the estimated height of the disk and stretch it to match the width.
    PA: is the angle at which the disk is rotated on the 2d image 
    incl: is the inclination of the disk/the angle it is facing away from us. '''
    # PL added - error statment, if data shape is 0 report that the image is empty
    if data is None or data.size == 0:
        raise ValueError(f"Image is empty - cannot deproject. Check the Crop column in Image Data in your Google Sheet, you are likely cropping off the whole image. Data shape of image: {data.shape}.")
    # PL done
    data = np.nan_to_num(data)
    clean = data.astype(np.float64)
   
    ndimx,ndimy = data.shape[1],data.shape[0]
    M = cv2.getRotationMatrix2D((xc,yc), PA-90, 1.0)

    imrot = cv2.warpAffine(clean, M,
                            (data.shape[0], data.shape[1]), flags = cv2.INTER_CUBIC)
    #stretch image in y by cos(incl) to deproject and divide by that factor to preserve flux
    im_rebin = cv2.resize(imrot,
                            (ndimx,int(np.round(ndimy*(1/np.cos(incl))))),
                            interpolation = cv2.INTER_CUBIC)
    im_rebin = im_rebin/(1./np.cos(incl))
    ndimy2 = im_rebin.shape[0]
    ycen = int(round((ndimy2/ndimy)*yc))
    #rotate back to original orientation
    M = cv2.getRotationMatrix2D((xc, ycen), -1*(PA - 90), 1.0)
    im_rebin_rot = cv2.warpAffine(im_rebin, M,
                                    (ndimx, ndimy2),
                                    flags = cv2.INTER_CUBIC)
    data = im_rebin_rot
    return data, ycen




def find_center(data, count): 
    ''' The function that uses checks to find the true center of the disk''' 
    ellipse_fit = None 
    circle = False #variable to control whether it is circular or not. 
    circle_ex = False #variable to keep track of if circle has ever been true. 
    ellipses = []
    grid = 100000 # grid initial value should no grid be smaller than this it will return bad fit. 
    gridcirc = 1000000 # same as above but only for circular fits. 

    #find the contours  of the data
    contours = contour_fits(data)
    
    # for each contour make an ellipse that keeps them in order. 
    for i in range(len(contours.levels)):
        current_cnts = contours.allsegs[i]
        ellipse_lvls = fit_ellipse(data, current_cnts)
        ellipses.append(ellipse_lvls)
    # make one array of all ellipses. 
    all_ellipses = [e for level in ellipses for e in level]
    #find the ellipse with the maximum area 
    max_rad = max(
    all_ellipses,
    key=lambda x: x[1][0] * x[1][1]) 

    #record all information about this largest ellipse 
    max_x = max_rad[0][0]
    max_y = max_rad[0][1]
    maxw = max_rad[1][0]
    maxh = max_rad[1][1]
    radius = maxw #assume radius should be the width of the whole disk. 
    #assume that the initial data center is the center
    xc_prior = data.shape[0]/2
    yc_prior = data.shape[0]/2

# First pass through evaluate if the center of the image is the center of the disk. 
    for i in ellipses: 
        for j in i:
            #for each ellipse 
            circle = False 
            xc = j[0][0]
            yc = j[0][1]
            axis_ratio = min(j[1][0], j[1][1]) / max(j[1][0], j[1][1])
            if 1-axis_ratio < .2: 
                circle = True # is fit circular ish
            if (abs(xc-xc_prior) < 2) & (abs(yc-yc_prior) < 2): 
                #is it within our prior guess
                w = j[1][0]
                h = j[1][1]
                angle = j[2]
                if (np.isnan(w) == True) or (np.isnan(h) == True) or (np.isnan(angle) == True): 
                    break
                #find inclination
                inclin = np.arccos(w/h)
                #deproject
                deprojected1, yc_fit1 = deproject(data, xc, yc, angle, inclin)
                #find fit value
                grid1 = circfit_good(deprojected1, xc, yc_fit1, radius, count)

                if grid1< grid:
                    #if fit is better than previous fit record it  
                    ellipse_fit = j 
                    deproject_fit = deprojected1
                    yc_fit = yc_fit1
                    grid = grid1 
                if (circle == True) & (grid1<gridcirc): 
                    #if fit is better than previous and its circular record it. 
                    ellipse_fitcirc = j 
                    deproject_fitcirc = deprojected1
                    yc_fitcirc = yc_fit1
                    gridcirc = grid1 
                    circle_ex = True


    
    if (circle_ex == True)& (gridcirc<.2): # if circle was ever true and the fit is good: (this number has been found through trial and error.)
        return ellipse_fitcirc, ellipse_fitcirc[0][1], gridcirc, data, data

    if grid < 0.08:#if fit was good (this number has been found through trial and error.)
        return ellipse_fit, yc_fit, grid, deproject_fit, data

    xc_prior, yc_prior = photutils(data, max_rad)
    #assume that the center is at the center of the brightness of the image 
    # repeat previous system 
    for i in ellipses: 
        for j in i: 
            xc = j[0][0]
            yc = j[0][1]
            if (abs(xc-xc_prior) < 1) & (abs(yc-yc_prior) < 1): 
                w = j[1][0]
                h = j[1][1]
                angle = j[2]
                if (np.isnan(w) == True) or (np.isnan(h) == True) or (np.isnan(angle) == True): 
                    break
                inclin = np.arccos(w/h)
                deprojected1, yc_fit1 = deproject(data, xc, yc, angle, inclin)
                radius = maxw
                grid1 = circfit_good(deprojected1, xc, yc_fit1, radius, count)

                print(grid1)
                if grid1< grid: 
                    ellipse_fit = j 
                    deproject_fit = deprojected1
                    yc_fit = yc_fit1
                    grid = grid1 

    if grid < .01: #Is new fit good? (this number has been found through trial and error.)
        return ellipse_fit, yc_fit, grid, deproject_fit, data
    
    #Now cycle through every ellipse since we haven't found it yet. 
    #repeat previous pattern. 
    for i in ellipses: 
        for j in i: 

            xc = j[0][0]
            yc = j[0][1]
            if (xc<(max_x+maxw/2)) & (xc>(max_x-maxw/2)) & (yc>(max_y-maxh/2)) & (yc<(max_y+maxh/2)):
                w = j[1][0]
                h = j[1][1]
                angle = j[2]
                inclin = np.arccos(w/h)
                if (np.isnan(w) == True) or (np.isnan(h) == True) or (np.isnan(angle) == True): 
                    break
                deprojected1, yc_fit1 = deproject(data, xc, yc, angle, inclin)
                radius =maxw

                grid1 = circfit_good(deprojected1, xc, yc_fit1, radius, count)

                if grid1<grid: 
                    ellipse_fit = j 
                    deproject_fit = deprojected1
                    yc_fit = yc_fit1
                    grid = grid1 

    if grid > 30: #is the fit quite shockingly bad. Maybe there is too much skew and we should use the median. 
        #(this number has been found through trial and error.)
        #repeat previous process, but make sure that the grid returns the median of itself and not the mean. 
        for i in ellipses: 
            for j in i:
                xc = j[0][0]
                yc = j[0][1]
                #make sure that the ellipses are all within the largest ellipse. 
                if (xc<(max_x+maxw/2)) & (xc>(max_x-maxw/2)) & (yc>(max_y-maxh/2)) & (yc<(max_y+maxh/2)):
                    w = j[1][0]
                    h = j[1][1]
                    angle = j[2]
                    if (np.isnan(w) == True) or (np.isnan(h) == True) or (np.isnan(angle) == True): 
                        break
                    inclin = np.arccos(w/h)
                    deprojected1, yc_fit1 = deproject(data, xc, yc, angle, inclin)
                    radius =maxw

                    grid1 = circfit_good(deprojected1, xc, yc_fit1, radius, count, True)

                    if grid1<grid: 
                        ellipse_fit = j 
                        deproject_fit = deprojected1
                        yc_fit = yc_fit1
                        grid = grid1 


    if ellipse_fit == None: 
        #if somethings gone horribly wrong just return not a number values. 
        return 'nan', 'nan', 'nan', 'nan', 'nan'

    return ellipse_fit, yc_fit, grid, deproject_fit, data





images = []

folder = Path('path/images')
#for each image in images folder find center. 
for item in folder.iterdir(): 
    if item.name.startswith('.'):
        #make sure that the item is not a secret hidden folder. 
        continue
    images.append(item.absolute())

for im in images: 
    with fits.open(im, ignore_missing_simple=True ) as hdul: 
        data = hdul[0].data #assume the first extension is an image 
    if len(data.shape) == 4:
        instrument = 'ALMA'
        data = data[0][0]
    elif len(data.shape) == 3: 
        instrument = 'SPHERE'
        data = data[0]
    elif len(data.shape) == 2: 
        instrument = 'Subaru'
        data = data
    else: 
        raise ValueError('Data is not from ALMA, Sphere, or Subaru')
    count = 20 # the number of ellipses to fit. 
    #run the funtion
    ellipse, yc_dep, grid, deproject_fit, data = find_center(data, count)

    # if something hasn't gone terribly wrong plot the deprojection and the non-deprojection with the center overlayed. 
    if yc_dep != 'nan':
        if grid == 'nan': 
            ellipse = [[deproject_fit.shape[0]/2]]
        plot_idea(deproject_fit, ellipse[0][0], yc_dep, False)
        plot_idea(data, ellipse[0][0], ellipse[0][1], False)
   
