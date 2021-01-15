#!/usr/bin/env python
#Embedded file name: MLDfitscut.py
""" ===================================================================================
NAME:       MLDfitscut

PURPOSE:    Search through the Hubble Legacy Archive for a fitscut generated image of a
            lens that is in the MasterLensDatabase

USAGE:      BEFORE STARTING THE ROUTINE:
            You must download the latest XML file from the MasterLens Database
            (at http://admin.masterlens.org/xml.php?) and place it in the same
            directory as this routine.  It also requires a folder named 'fitscutimages'
            in its shared directory.

            You must also set your MLD login and username (for uploading images) under the
            first commmand below: uploadfitscutimages.
            (Look for the comment "USER NEEDS TO EDIT CODE HERE").
        This all needs to be modernized of course, but whatever for testing purposes!

            TO START THE ROUTINE:
            To automatically start looping through the entire database, run from terminal:
            python MLDfitscut COMMAND1 COMMAND2
            COMMAND1 indicates the HLA instrument: either 'ACS', 'WPC3', 'WFPC2'
            COMMAND2 indicates whether you want to loop over the *entire* database (COMMAND2='all')
                     or only over lenses that currently have no fitscut image (COMMAND2='new')

            To run the script on only a single specific lens, run in python:
            import MLDfitscut
            from MLDfitscut import getfitscut
            getfitscut([LENS-NAME], [INSTRUMENT])

            LENS-NAME should be the exact lens name from the MasterLensDatabase
            INSTRUMENT should again either be 'ACS', 'WPC3', or 'WFPC2'

            DURING THE ROUTINE:
            1) After finding a good fitscut image, the script will open an image and ask you to
            click on the center of the lens.  If the lens was very far off center you will be
            promted in terminal to overwrite the *MANUAL* RA/DEC coordinates in the database
            with the new values.

            2) You will be prompted for a 'zoom factor' based on how small/large the lens looked
            in the previous image.  (Most non-cluster lenses require a zoom of 2-4 for ACS).  You
            may also zoom out by providing a number between 0 and 1.

            3) The routine will then find the single image with the longest exposure, display a
            black and white image, and prompt (if the image looks good) to upload this image to the
            MasterLens Database.

            4) *IF* there are multiple filters of the same image, a color image will be construted,
            and you will again be prompted to whether or not the image should be uploaded to the database.
            
REVISION HISTORY:
02-Jan-2021 Modified by Michael S. Talbot, University of Utah. 
    -Updated "acsSIAP" to "hlaSIAP" in HLA query url to obtain target table.
    -Added "strip_bytes" function and stripped bytes in strings within the returned tables (this enabled the script to be ran on python3).
    -Added a minimal patch to the getlenscenter function to enable the plot to be displayed in python3, though at the cost of creating a useless extra panel in the plot creation. I did not edit this more since it appears to be deprecated. Thus I left the function command as commented out.
    -Added a temporary copy/paste of the file being processed by svg to the svg location for robustness on multi-systems.
    -For these edits to work, I had to edit the lenses.xml file via the following:
        -Remove all 'Description' sections from the lenses.xml file being used since minidom does not recognize HTML entry names and these sections were not always formatted correctly.
        -Add a space before Modified
        -Remove the several entry names unique to HTML and not minidom.
    -NOTE I REMOVED A BAKED IN USER AND PASSWORD for the masterlens database
08-Feb-2013  Written by Kyle R. Stewart, California Baptist University based on process and structure developed by Leonidas Moustakas and Tim Goodsall

REQUIRED PACKAGES:
(This is all old stuff, and we should be using conda or pip)

requests     Obtain it from http://docs.python-requests.org/en/latest/user/install/
             [or by macports via sudo port install py27-requests]
Tkinter      Macports: sudo port install py27-tkinter
astropy      Macports: sudo port install py27-astropy
PIL          Macports: sudo port install py27-pil

NOTE: This routine was designed for use on a Mac, and may not function on other operating systems.
(This shouldn't be a factor going forward)

=================================================================================== """
from numpy import array
import numpy as np
import xml.dom.minidom
import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import os
import sys
import astropy.io.votable as votable
import tkinter
from PIL import Image, ImageTk, ImageDraw
import warnings
import subprocess
import pandas as pd
from time import sleep

def uploadfitscutimages(uploadID, graphic_files, user = '', password = ''):
    """Upload the fitscut images to the masterlens database website"""
    from masterlens import masterlens
    m = masterlens()
    if m.ready:
        m.set_login(user, password)
        m.set_graphic_type('HLA-Derived Fitscut')
        m.set_lensID(uploadID)
        for graphic_file in graphic_files:
            m.set_graphic_file(graphic_file)
            m.upload_graphic()
            if m.graphicID != None:
                m.message('Uploaded %s as graphicID=%i' % (graphic_file, m.graphicID))
            else:
                m.error('Unable to load %s' % graphic_file)


def getlenscenter(image_file):
    """Opens an image, waits for mouse click, returns the x,y pixel of the click"""
    
    w = tkinter.Tk()
    
    lens_center = [256, 256]

    def callback(event):
        lens_center[0] = event.x
        lens_center[1] = event.y
        w.destroy()

    img = Image.open(image_file)
    width, height = img.size
    ca = tkinter.Canvas(w, width=width, height=height)
    ca.pack()
    photoimg = ImageTk.PhotoImage('RGB', img.size)
    photoimg.paste(img)
    ca.create_image(width // 2, height // 2, image=photoimg)
    ca.bind('<Button-1>', callback)
    ca.pack()
    tkinter.mainloop()
    return lens_center


def fitscutimagetest(f1, x, y, name):
    """tests the HLA fitscut script to see if there is a valid image, uses a 256x256 black and white, single filter image"""
    imagename = 'testimage_' + str(name) + '_' + str(f1) + '.jpg'
    urlstring = 'http://hla.stsci.edu/cgi-bin/fitscut.cgi?blue=' + str(f1) + '&size=256,256&x=' + str(x) + '&y=' + str(y) + '&wcs=1&zoom=0.5&compass=0&invert=1'
    try:
        print('Image test urlstring', urlstring)
        opener1 = urllib.request.build_opener()
        page1 = opener1.open(urlstring)
        my_picture = page1.read()
        filename = 'fitscutimages/testimages/' + imagename
        fout = open(filename, 'wb')
        fout.write(my_picture)
        fout.close()
        f = open('testimage_success_urls.txt', 'a')
        f.write(imagename + '     ' + urlstring + '\n')
        f.close()
        im1 = Image.open('fitscutimages/testimages/' + imagename)
        if im1.getbbox() == None:
            return False
        return True
    except (urllib.error.HTTPError, urllib.error.URLError, NameError) as e:
        return False


def fitscutimage(f1, f2, x, y, zoomlevel, imagename):
    """Generates the HLA fitcut image.
    If f2 is not an empty string, it will create a color image
    If f2 is an empty string, it will produce a black and white (inverted) single filter image"""
    pixsize = str(int(512 / float(zoomlevel)))
#    urlstring = 'http://hla.stsci.edu/cgi-bin/fitscut.cgi?red=' + str(f1) + '&green=' + str(f2) + '&size=' + pixsize + '&x=' + str(x) + '&y=' + str(y) + '&wcs=1&zoom=' + str(zoomlevel) + '&compass=1'

    urlstring = 'http://hla.stsci.edu/cgi-bin/fitscut.cgi?red={:s}&green={:s}&size={:s}&x={:s}&y={:s}&wcs=1&zoom={:s}&compass=1'.format(f1, f2, str(pixsize),
            str(x), str(y), str(zoomlevel))


    if f2 == '':
        urlstring = 'http://hla.stsci.edu/cgi-bin/fitscut.cgi?blue=' + str(f1) + '&size=' + pixsize + '&x=' + str(x) + '&y=' + str(y) + '&wcs=1&zoom=' + str(zoomlevel) + '&compass=1&invert=1'
    print('Getting FITSCUT image at: ' + urlstring)
    try:
        opener1 = urllib.request.build_opener()
        page1 = opener1.open(urlstring)
        my_picture = page1.read()
        filename = os.path.join('fitscutimages', 'raw', '{:s}_ZOOM{:s}.jpg'.format(imagename, str(zoomlevel)))
        fout = open(filename, 'wb')
        fout.write(my_picture)
        fout.close()
        if str(1) in str(zoomlevel):
            if f2 != '':
                fout_name = 'success_urls_color.txt'
            if f2 == '':
                fout_name = 'success_urls_grey.txt'
            f = open(fout_name, 'a')
            f.write(imagename + '     ' + urlstring + '\n')
            f.close()
    except (urllib.error.HTTPError, urllib.error.URLError, NameError) as e:
        print('Error Ocurred.  NO image will be created!')


def valuematch(x, y):
    """given two arrays, x and y, look for the first element of x
    that has an identical value in y"""
    for i in x:
        if i in y:
            return i


def makesvg(lensNAME, ins_name, f1, f2, pixscale, fname, zoomlevel):
    """add caption to image, including arcsec scale, lens name, filters used.  Assumes the jpg needing a caption is already in folder fitscutimages/ from current directory"""
    
    #Different functions might pass in pixscale different.
    if isinstance(pixscale, str): pass
    else:
        try: pixscale = float(pixscale[0])
        except: pass

    pixtext = str(420.0 + 0.5 / (float(pixscale) * 3600) * float(zoomlevel))
    pixscale = str(1.0 / (float(pixscale) * 3600) * float(zoomlevel))
    #pixtext = ''
    #pixscale = ''

    filter1 = f1
    filter2 = f2
    filename = fname
    instrument = ins_name
    cap = str(' ' + instrument + '-' + filter1 + '/' + filter2)
    textcolor = 'white'
    if filter2 == '':
        cap = str(' ' + instrument + '-' + filter1)
        textcolor = 'black'
    
    f = open(os.path.join('fitscutimages', 'svg', '{:s}.svg'.format(filename)), 'w')
    f.write('<?xml version="1.0"?>\n')
    f.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
    f.write('<svg height="512" width="512" version="1.1"\n')
    f.write('baseProfile="full"\n')
    f.write('xmlns="http://www.w3.org/2000/svg"\n')
    f.write('xmlns:xlink="http://www.w3.org/1999/xlink"\n')
    f.write('xmlns:ev="http://www.w3.org/2001/xml-events">\n')
    f.write('<defs>\n')
    f.write('<pattern id="img1" patternUnits="userSpaceOnUse" width="512" height="512">\n')
    
    #Michael Talbot. Python-3.7.3: Appears to only work when image is in same location as being saved. Swapped "/raw/" for "/svg/" but check!!! 
    f.write('    <image xlink:href="../svg/{:s}.jpg" x="0" y="0" width="512" height="512"/>\n'.format(filename))
    f.write('</pattern>\n')
    f.write('</defs>\n')
    f.write('<path d="M0,0\n')
    f.write('         L0,512 L512,512 L512,0\n')
    f.write('         z"\n')
    f.write('      fill="url(#img1)" />\n')
    f.write('<path d="M420,20 l0,3 l' + pixscale + ',0 l0,-3 z"\n')
    f.write('      fill="' + textcolor + '"  />\n')
    f.write('<text fill="' + textcolor + '" style="font-family:Sans;font-size:18px;text-anchor:middle;dominant-baseline:bottom" x="' + pixtext + '" y="40">1&#x201D;</text>\n')
    f.write('<text fill="' + textcolor + '" style="font-family:Sans;font-size:36px;text-anchor:middle;dominant-baseline:bottom" x="256" y="464">' + lensNAME + '</text>\n')
    f.write('<text fill="' + textcolor + '" style="font-family:Sans;font-size:36px;text-anchor:middle;dominant-baseline:bottom" x="256" y="500">' + cap + '</text>\n')
    f.write('</svg>')
    f.close()    

def gen_color_image(IDs, detectors, scales, lensRA, lensDEC, zoomlevel, lensNAME, instrument, filter1, filter2, f1_str, f2_str):

    if filter1.any() & filter2.any():

        green_ID = IDs[filter1][detectors[filter1] == valuematch(detectors[filter1], detectors[filter2])]
        red_ID = IDs[filter2][detectors[filter2] == valuematch(detectors[filter1], detectors[filter2])]
        scale = scales[filter1][detectors[filter1] == valuematch(detectors[filter1], detectors[filter2])]

        fname = '{:s}_{:s}_{:s}and{:s}'.format(lensNAME, instrument, f1_str, f2_str).replace(' ', '')
        zoom_fname = '{:s}_ZOOM{:s}'.format(fname, str(zoomlevel))

        fitscutimage(red_ID[0], green_ID[0], lensRA, lensDEC, str(zoomlevel), fname)
        makesvg(lensNAME, instrument, 'F125W', 'F160W', scale, zoom_fname, zoomlevel)
    else:
        zoom_fname = ''

    return zoom_fname
    
def strip_bytes(table):
    for i in range(len(table)):
        if isinstance(table.data[i], bytes): table.data[i] = table.data[i].decode()
    return table

def getfitscut(lensNAME, instrument, df, pargs):
    
    """Creates a captioned HLA fitscut image for a given lens and instrument.
    Lens Name must match the MasterLens Database.  Valid instruments include:
        ACS, WFPC2, WFC3"""

    # set up directory structures
    if not os.path.exists(os.path.join(pargs.od, 'fitscutimages', 'bylens', lensNAME)):
        os.makedirs(os.path.join(pargs.od, 'fitscutimages', 'bylens', lensNAME))


    # pull the values from the dataframe using the lens name
    lensRA = df.loc[lensNAME].ra
    lensDEC = df.loc[lensNAME].dec
    lensID = df.loc[lensNAME].id

    #Michael Talbot: Swapped 'acs' to 'hlaSIAP' in ###SIAP
    '''hla_url = 'http://hla.stsci.edu/cgi-bin/acsSIAP.cgi?pos=' + str(lensRA) + ',' + str(lensDEC) + '&size=0.02&format=FITS&inst=' + instrument'''
    
    hla_url = 'http://hla.stsci.edu/cgi-bin/hlaSIAP.cgi?pos=' + str(lensRA) + ',' + str(lensDEC) + '&size=0.02&format=FITS&inst=' + instrument
    filename = 'lens_votable.xml'
    print('Getting votable: ', hla_url)
    print('Getting HLA information...')

    urllib.request.urlretrieve(hla_url, filename)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tab = votable.parse(filename, pedantic=False, columns=['ExpTime',
         'Level',
         'Detector',
         'Aperture',
         'Spectral_Elt',
         'Dataset',
         'scale'])
        tab1 = tab.get_first_table()
    exposure = tab1.array['ExpTime']
    levels = tab1.array['Level']
    detectors = tab1.array['Detector']
    instruments = tab1.array['Aperture']
    filters = tab1.array['Spectral_Elt']
    IDs = tab1.array['Dataset']
    scales = tab1.array['scale']
    ind = np.argsort(exposure)[::-1]
    
    #Michael Talbot: Strip the byte format of the strings makes the program python3 compatible.
    exposure = strip_bytes(exposure[ind])
    levels = strip_bytes(levels[ind])
    detectors = strip_bytes(detectors[ind])
    instruments = strip_bytes(instruments[ind])
    filters = strip_bytes(filters[ind])
    IDs = strip_bytes(IDs[ind])
    scales = strip_bytes(scales[ind])
    for i in range(len(scales)):
        if len(scales[i]) > 1:
            scales[i] = scales[i][0]

    ind = ind[scales > 0]
    exposure = exposure[scales > 0]
    levels = levels[scales > 0]
    detectors = detectors[scales > 0]
    instruments = instruments[scales > 0]
    filters = filters[scales > 0]
    IDs = IDs[scales > 0]
    scales = scales[scales > 0]
    F125 = filters == 'F125W'
    F160 = filters == 'F160W'
    F475 = filters == 'F475X'
    F600 = filters == 'F600LP'
    F555 = filters == 'F555W'
    F606 = filters == 'F606W'
    F775 = filters == 'F775W'
    F814 = filters == 'F814W'
    F850 = filters == 'F850LP'
    print('Searching for all HLA images...')
    goodimage = [False] * len(IDs)
    if len(IDs) == 0:
        goodimage = False
        print('NO IMAGES FOUND')
    for j in range(len(IDs)):
        goodimage[j] = fitscutimagetest(IDs[j], lensRA, lensDEC, lensNAME)

    imageexists = array((F125 | F160 | F475 | F600 | F555 | F606 | F814 | F850) & goodimage, dtype=bool)

    if imageexists.any():
        repeatcentering = True

        blue_ID = IDs[imageexists][0]

        if pargs.center:

            fitscutimage(blue_ID, '', lensRA, lensDEC, str(1), (lensNAME + '_' + instrument + '_greyimage_' + filters[imageexists][0]).replace(' ', ''))
            print('\nPLEASE CLICK ON LENS CENTER.')
            lens_center = [256, 256]
            if pargs.center: lens_center = getlenscenter('fitscutimages/raw/' + (lensNAME + '_' + instrument + '_greyimage_' + filters[imageexists][0] + '_ZOOM1.jpg').replace(' ', ''))
            if abs(lens_center[0] - 256) > 20 or abs(lens_center[1] - 256) > 20:
                print('NEW RA/DEC: ', lensRA, lensDEC)
                print('(PLEASE UPDATE COORDINATES IN MASTER DATABASE!)\n')
                lensRA = lensRA + ((256 - lens_center[0]) * scales[imageexists][0])[0]
                lensDEC = lensDEC + ((256 - lens_center[1]) * scales[imageexists][0])[0]


        print('Creating images with captions...')
        outputnames = []
        if pargs.zoom is None:
            zoomlevel = eval(input('Zoom in on image by factor of: '))
        else:
            zoomlevel = pargs.zoom

        fitscutimage(blue_ID, '', lensRA, lensDEC, str(zoomlevel), (lensNAME + '_' + instrument + '_greyimage_' + filters[imageexists][0]).replace(' ', ''))

        makesvg(lensNAME, instrument, filters[imageexists][0], '', scales[imageexists][0], (lensNAME + '_' + instrument + '_greyimage_' + filters[imageexists][0] + '_ZOOM' + str(zoomlevel)).replace(' ', ''), zoomlevel)
        outputnames.append((lensNAME + '_' + instrument + '_greyimage_' + filters[imageexists][0] + '_ZOOM' + str(zoomlevel)).replace(' ', ''))
        
        F125 = (filters == 'F125W') & goodimage
        F160 = (filters == 'F160W') & goodimage
        F475 = (filters == 'F475X') & goodimage
        F600 = (filters == 'F600LP') & goodimage
        F555 = (filters == 'F555W') & goodimage
        F606 = (filters == 'F606W') & goodimage
        F814 = (filters == 'F814W') & goodimage
        F850 = (filters == 'F850LP') & goodimage


        zoom_fname = gen_color_image(IDs, detectors, scales, lensRA, lensDEC, zoomlevel, lensNAME, instrument, F125, F160, 'F125W', 'F160W')
        outputnames.append(zoom_fname)

        zoom_fname = gen_color_image(IDs, detectors, scales, lensRA, lensDEC, zoomlevel, lensNAME, instrument, F475, F600, 'F475X', 'F600LP')
        outputnames.append(zoom_fname)

        zoom_fname = gen_color_image(IDs, detectors, scales, lensRA, lensDEC, zoomlevel, lensNAME, instrument, F555, F814, 'F555W', 'F814W')
        outputnames.append(zoom_fname)

        zoom_fname = gen_color_image(IDs, detectors, scales, lensRA, lensDEC, zoomlevel, lensNAME, instrument, F555, F850, 'F555W', 'F850LP')
        outputnames.append(zoom_fname)

        zoom_fname = gen_color_image(IDs, detectors, scales, lensRA, lensDEC, zoomlevel, lensNAME, instrument, F606, F814, 'F606W', 'F814W')
        outputnames.append(zoom_fname)

        zoom_fname = gen_color_image(IDs, detectors, scales, lensRA, lensDEC, zoomlevel, lensNAME, instrument, F606, F850, 'F606W', 'F850LP')
        outputnames.append(zoom_fname)


#            if F125.any() & F160.any():
#                green_ID = IDs[F125][detectors[F125] == valuematch(detectors[F125], detectors[F160])]
#                red_ID = IDs[F160][detectors[F160] == valuematch(detectors[F125], detectors[F160])]
#                scale = scales[F125][detectors[F125] == valuematch(detectors[F125], detectors[F160])]
#                fitscutimage(red_ID[0], green_ID[0], lensRA, lensDEC, str(zoomlevel), (lensNAME + '_' + instrument + '_F125WandF160W').replace(' ', ''))
#                makesvg(lensNAME, instrument, 'F125W', 'F160W', scale, (lensNAME + '_' + instrument + '_F125WandF160W_ZOOM' + str(zoomlevel)).replace(' ', ''), zoomlevel)
#                outputnames.append((lensNAME + '_' + instrument + '_F125WandF160W_ZOOM' + str(zoomlevel)).replace(' ', ''))
#            if F475.any() & F600.any():
#                green_ID = IDs[F475][detectors[F475] == valuematch(detectors[F475], detectors[F814])]
#                red_ID = IDs[F600][detectors[F600] == valuematch(detectors[F475], detectors[F600])]
#                scale = scale[F600][detectors[F600] == valuematch(detectors[F475], detectors[F600])]
#                fitscutimage(red_ID[0], green_ID[0], lensRA, lensDEC, str(zoomlevel), (lensNAME + '_' + instrument + '_F475WandF600LP').replace(' ', ''))
#                makesvg(lensNAME, instrument, 'F475W', 'F600LP', scales, (lensNAME + '_' + instrument + '_F475WandF600LP_ZOOM' + str(zoomlevel)).replace(' ', ''), zoomlevel)
#                outputnames.append((lensNAME + '_' + instrument + '_F475WandF600LP_ZOOM' + str(zoomlevel)).replace(' ', ''))
#            if F555.any() & F814.any():
#                green_ID = IDs[F555][detectors[F555] == valuematch(detectors[F555], detectors[F814])]
#                red_ID = IDs[F814][detectors[F814] == valuematch(detectors[F555], detectors[F814])]
#                scale = scales[F814][detectors[F814] == valuematch(detectors[F555], detectors[F814])]
#                fitscutimage(red_ID[0], green_ID[0], lensRA, lensDEC, str(zoomlevel), (lensNAME + '_' + instrument + '_F555WandF814W').replace(' ', ''))
#                makesvg(lensNAME, instrument, 'F555W', 'F814W', scales, (lensNAME + '_' + instrument + '_F555WandF814W_ZOOM' + str(zoomlevel)).replace(' ', ''), zoomlevel)
#                outputnames.append((lensNAME + '_' + instrument + '_F555WandF814W_ZOOM' + str(zoomlevel)).replace(' ', ''))
#            if F555.any() & F850.any():
#                green_ID = IDs[F555][detectors[F555] == valuematch(detectors[F555], detectors[F850])]
#                red_ID = IDs[F850][detectors[F850] == valuematch(detectors[F555], detectors[F850])]
#                scale = scales[F850][detectors[F850] == valuematch(detectors[F555], detectors[F850])]
#                fitscutimage(red_ID[0], green_ID[0], lensRA, lensDEC, str(zoomlevel), (lensNAME + '_' + instrument + '_F555WandF850LP').replace(' ', ''))
#                makesvg(lensNAME, instrument, 'F555W', 'F850LP', scales, (lensNAME + '_' + instrument + '_F555WandF850LP_ZOOM' + str(zoomlevel)).replace(' ', ''), zoomlevel)
#                outputnames.append((lensNAME + '_' + instrument + '_F555WandF850LP_ZOOM' + str(zoomlevel)).replace(' ', ''))
#            if F606.any() & F814.any():
#                green_ID = IDs[F606][detectors[F606] == valuematch(detectors[F606], detectors[F814])]
#                red_ID = IDs[F814][detectors[F814] == valuematch(detectors[F606], detectors[F814])]
#                scale = scales[F814][detectors[F814] == valuematch(detectors[F606], detectors[F814])]
#                fitscutimage(red_ID[0], green_ID[0], lensRA, lensDEC, str(zoomlevel), (lensNAME + '_' + instrument + '_F606WandF814W').replace(' ', ''))
#                makesvg(lensNAME, instrument, 'F606W', 'F814W', scales, (lensNAME + '_' + instrument + '_F606WandF814W_ZOOM' + str(zoomlevel)).replace(' ', ''), zoomlevel)
#                outputnames.append((lensNAME + '_' + instrument + '_F606WandF814W_ZOOM' + str(zoomlevel)).replace(' ', ''))
#            if F606.any() & F850.any():
#                green_ID = IDs[F606][detectors[F606] == valuematch(detectors[F606], detectors[F850])]
#                red_ID = IDs[F850][detectors[F850] == valuematch(detectors[F606], detectors[F850])]
#                scale = scales[F850][detectors[F850] == valuematch(detectors[F606], detectors[F850])]
#                fitscutimage(red_ID[0], green_ID[0], lensRA, lensDEC, str(zoomlevel), (lensNAME + '_' + instrument + '_F606WandF850LP').replace(' ', ''))
#                makesvg(lensNAME, instrument, 'F606W', 'F850LP', scales, (lensNAME + '_' + instrument + '_F606WandF850LP_ZOOM' + str(zoomlevel)).replace(' ', ''), zoomlevel)
#                outputnames.append((lensNAME + '_' + instrument + '_F606WandF850LP_ZOOM' + str(zoomlevel)).replace(' ', ''))

        print('Making iPad/iPhone/Web/Main versions...')


        for basename in outputnames:

            if basename == '':
                continue
            convert_svg_png(os.path.join('fitscutimages', 'svg', basename), os.path.join('fitscutimages', 'annotated', basename))
            try:
                

                try:

                    im1 = Image.open('{:s}.png'.format(os.path.join('fitscutimages', 'annotated', basename)))

                    try:
                        gen_iphone_image(im1, basename)
                    except Exception as e:
                        print("Error: unable to generate iphone image")
                    try:
                        gen_ipad_image(im1, basename)
                    except Exception as e:
                        print("Error: unable to generate ipad image")
                    try:
                        gen_web_image(im1, basename)
                    except Exception as e:
                        print("Error: unable to generate web image")

                except Exception as e:
                    print(str(e))

    #                subprocess.call('open -g fitscutimages/' + basename + '.png', shell=True)

            except Exception as e:
                print("Unable to convert svg to png")
                print(str(e))
                print("Removing broken png...")
                try:
                    os.remove(os.path.join('fitscutimages', 'annotated', basename+'.png'))
                except:
                    pass


            if pargs.upload:

                continueYN = input('Upload ' + basename + ' to MASTERLENS DATABASE (y/n)?')

                if continueYN == 'y':
                    uploadfitscutimages(lensID, [ \
                     os.path.join(pargs.od, 'fitscutimages', 'iphone', '{:s}.iPhone.png'.format(basename)),
                     os.path.join(pargs.od, 'fitscutimages', 'ipad', '{:s}.iPad.png'.format(basename)),
                     os.path.join(pargs.od, 'fitscutimages', 'web', '{:s}.web.png'.format(basename)),
                     os.path.join(pargs.od, 'fitscutimages', 'web', '{:s}.weblarger.png'.format(basename)),
                     os.path.join(pargs.od, 'fitscutimages', 'annotated', '{:s}.png'.format(basename))])

            if pargs.center:

                if input('Repeat Centering and Zoom-in Process (y/n)? ') == 'n':
                    repeatcentering = False


def convert_svg_png(svg_fname, png_fname):

    #Michael Talbot: Robustness fix. It appears the authors have the correct patch to enable the svg process to link to images beyond its location. In the mean time, I modified the copy and remove commands (were commented out) to enable me to properly create the png image. Consider commenting out the cp and rm commands.
    image_file = os.path.join(svg_fname.replace('/svg/','/raw/'))
    
    subprocess.call('cp %s.jpg %s.jpg'%(image_file, svg_fname), shell=True)
    sleep(1)
    input('CHECK FILES EXISTS>%s.jpg'%svg_fname)
    subprocess.check_call('convert {:s}.svg {:s}.png'.format(svg_fname, png_fname), shell=True)
    sleep(1.5)
    subprocess.call('rm {:s}.jpg'.format(svg_fname), shell=True)

def gen_iphone_image(im, basename):
    """Generate image for iphone"""
    im_iphone = im.resize((480, 480), Image.BICUBIC).crop((0, 175, 480, 305))
    im_iphone.save(os.path.join('fitscutimages/', '{:s}.iPhone.png'.format(basename)))

def gen_ipad_image(im, basename):
    """Generate image for ipad"""
    im_ipad = im.resize((1024, 1024), Image.BICUBIC).crop((0, 373, 1024, 650))
    im_ipad.save(os.path.join('fitscutimages/', '{:s}.iPad.png'.format(basename)))

def gen_web_image(im, basename):
    """Generate image for web"""
    
    # small image
    im_web = im.resize((600, 600), Image.BICUBIC).crop((0, 200, 600, 400))
    im_web.save(os.path.join('fitscutimages/', '{:s}.web.png'.format(basename)))

    # large image
    im_web_large = im.resize((600, 600), Image.BICUBIC)
    im_web.save(os.path.join('fitscutimages/', '{:s}.weblarge.png'.format(basename)))


def main(pargs):
    """Loops over the Master Lens Database and tries to get HLA fitscut images for the lenses.
    First argument is the instrument used ("ACS", "WFPC2", "WFC3", "all").
    Second argument: "all" run through the entire database -- "new" only runs on lenses that don't currently have HLA images"""

    instrument = pargs.instrument
    lens = pargs.lens
    doc = xml.dom.minidom.parse(pargs.xmlfile)
    lenses = doc.getElementsByTagName('lens')
    names = []
    ra = []
    dec = []
    grade = []
    has_image = []
    counter = 0
    ids = []

    DEBUG = pargs.DEBUG
    print("DEBUG:", DEBUG)

    # parse the xml tree regardless of lens mode (i.e only wanting to look at a
    # single lens or 'all' of them

    for ix, lens_i in enumerate(lenses[:]):

        graphics = lens_i.getElementsByTagName('graphic')

        if len(graphics) == 0:
            has_image.append(False)
        else:
            # get the graphics types and if any of them are a 1 then a 'good' fitscut exists. 
            # Otherwise need to generate a fitscut image later on
            types = [g.getAttribute('type') for g in graphics]
            if '1' in types:
                has_image.append(True)
            else:
                has_image.append(False)

        names.append(str(lens_i.getElementsByTagName('system_name')[0].childNodes[0].nodeValue))
        ra.append(float(str(lens_i.getElementsByTagName('ra_coord')[0].childNodes[0].nodeValue)))
        dec.append(float(str(lens_i.getElementsByTagName('dec_coord')[0].childNodes[0].nodeValue)))
        ids.append(float(str(lens_i.getAttribute('lensID'))))


        if len(lens_i.getElementsByTagName('lensgrade')) != 0:
            grade.append(str(lens_i.getElementsByTagName('lensgrade')[0].childNodes[0].nodeValue))
        else:
            grade.append('-1')

    df = pd.DataFrame({'name':names, 'ra':ra, 'dec':dec, 'grade':grade, 'id':ids, 'image':has_image})
    df = df[['name', 'ra', 'dec', 'grade', 'id', 'image']]
    df.index = df['name']

    ra = np.array(ra)
    dec = np.array(dec)
    names = np.array(names, dtype=str)

    if instrument == 'all':
        for ix, row in df.iterrows():
            if  (pargs.lens not in ['all']) and row['name'] != pargs.lens:
                continue
            if row.image == False:
                print('Getting fitscut for lens {:s}: RA: {:.3f} DEC: {:.3f}'.format(str(row['name']), row.ra, row.dec))
                if DEBUG:
                    continue
                getfitscut(row.name, 'ACS', df, pargs)
                getfitscut(row.name, 'WFPC2', df, pargs)
                getfitscut(row.name, 'WPC3', df, pargs)
    else:
        for ix, row in df.iterrows():
            if  (pargs.lens not in ['all']) and row['name'] != pargs.lens:
                continue
            print('Getting fitscut for lens {:s}: RA: {:.3f} DEC: {:.3f}'.format(str(row['name']), row.ra, row.dec))
            if DEBUG:
                continue
            getfitscut(row['name'], instrument, df, pargs)

    return df


if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('instrument', help='Instrument to be used for generating fitscut images')
    parser.add_argument('lens', help='generate fitscut images on all lenses or only those without existing fitscut images')
    parser.add_argument('--DEBUG', action='store_true', default=False)
    parser.add_argument('--upload', action='store_true', default=False, help='Upload the generated images to the MLD site')
    parser.add_argument('--center', action='store_true', default=False, help='Interatively center the images')
    parser.add_argument('--zoom', default=None, help='use a default zoom value, otherwise specify interatively')
    parser.add_argument('--xmlfile', default='lenses.xml', help='full path to the lenses.xml file you want to use to run this code')
    parser.add_argument('--od', default='.', help='path to fitscutimages/ directory')

    pargs = parser.parse_args()

    fitsdir = os.path.join(pargs.od, 'fitscutimages')
    if not os.path.exists(fitsdir):
        os.mkdir(fitsdir)
    if not os.path.exists(os.path.join(fitsdir, 'raw')):
        os.mkdir(os.path.join(fitsdir, 'raw'))
    if not os.path.exists(os.path.join(fitsdir, 'annotated')):
        os.mkdir(os.path.join(fitsdir, 'annotated'))
    if not os.path.exists(os.path.join(fitsdir, 'iphone')):
        os.mkdir(os.path.join(fitsdir, 'iphone'))
    if not os.path.exists(os.path.join(fitsdir, 'ipad')):
        os.mkdir(os.path.join(fitsdir, 'ipad'))
    if not os.path.exists(os.path.join(fitsdir, 'web')):
        os.mkdir(os.path.join(fitsdir, 'web'))
    if not os.path.exists(os.path.join(fitsdir, 'svg')):
        os.mkdir(os.path.join(fitsdir, 'svg'))
    if not os.path.exists(os.path.join(fitsdir, 'testimages')):
        os.mkdir(os.path.join(fitsdir, 'testimages'))


    df = main(pargs)

