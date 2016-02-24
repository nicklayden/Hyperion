
################# Phoenix Spectral Lines #####################
#
#  This module imports a special .lines.7 phoenix spectrum file 
#  and finds the strongest absorbing species listed in the file
#  for each wavelength, and determines the element causing any 
#  spectral line at the line center.
#
#
#
import numpy as np
import os
import periodictable as pt
from scipy.signal import fftconvolve


def file_loading(directory, extension=None,teffs=None,contains=None,metallicity=None,notcontains=None):
    '''
        Loads PHX files into a list for all of my other modules to find spectra.
        Default extension is .7, for ALL phx files. 
        teffs is the list or array or temperature values to scan through, and the metallicity
        is default 0.5. Since these are all metal-poor spectra files, the negative is implied 
        on all metallicity values.
        
        contains    : if file contains this ____ string, include it.
        notcontains : if file contains this ____ string, do NOT include it.   
        teffs       : LIST object of effective temperature values to look for. Cannot input single float values.
        metallicity : singular metallicity value to look for. Not implementing metallicity range here.
        extension   : phx file extensions are long and specific, easy way to find exactly what you want with this.
                      default extension finds any .7 files.
        
        
        teffs:    Note, although you cannot input a single value like 'teffs=4250', you CAN input a single valued 
                  list object, i.e.   'teffs=[4250]'. This should work correctly.
        
    '''
    if extension == None:
        extension = '.7'
    if teffs == None:
        teffs = [4000,4125,4250,4375,4500]
    if metallicity ==None:
        metallicity = 0.5
    if contains == None:
        contains = '-2.0-'
    if notcontains == None:
        notcontains = 'GGG'
    
    file_list = list()    #Defines an empty list to put spectra filenames in.
    for file in os.listdir(directory):    #My spectra are in this directory.
       for i in range(0,len(teffs)):
           if file.endswith(extension) and ('{}-2.0-{}'.format(teffs[i],metallicity)) in file and contains in file and notcontains not in file:
               file_list.append(file)
    
    file_list.sort(key = lambda x: x.split('-')[0] )
    file_list.sort(key = lambda x: x.split('-')[1], reverse=True )   

    return file_list


def DtoE(spec):
    '''
        Changes the Fortran 'D' in phoenix .7 files into the 
        everyone else friendly 'e' for scientific notation.
    '''
    for j in range(0,len(spec)):
        spec[j] = spec[j].replace("D", "E")
        
#####  Creating a Gaussian curve for the convolution kernel #####
def convolve(fluxarray,Velocity=None,gridspace=None):
    ''' 
        Gaussian Velocity = 1 by default, unless specified.
        gridspace = 0.01 by default, this is the finest gridspacing in the
        non-uniform phoenix files.
    '''
    if gridspace==None:
        gridspace = 0.01
    if Velocity ==None:
        Velocity = 1.0
    
    
    FWHM = (5000.*Velocity)/3e5   # Velocity in kms-1.
    x = gridspace     #Uniform Grid spacing on both spectra
    sigma = float(FWHM)/np.sqrt(8.*np.log(2))/x       #Standard Deviation in pixels
    X = np.arange(int(-4.*sigma), int(4.*sigma)+1)    #Range that the gaussian covers
    gaussian = (1./(sigma*np.sqrt(2*np.pi)))/  ( np.exp(0.5*(X/sigma)**2) ) #Actual gaussian function
    final = fftconvolve(fluxarray, gaussian, mode='same')
    return final


################################################################################
#
#           (nlte-lte)/lte Relative Spectrum Differences
#               Smoothed by convolution
################################################################################
#
#   This function reads in the special lines.7 phoenix file for a 4250-2.0-0.5 
#   star to gather the strongest absorbing species at each wavelength in the 
#   spectrum, and outputs a list of the strongest spectral lines in the spectrum
#   with a specified optical depth difference.
#
#

def phoenix_line_list(opticaldepthdiff,elementcode):
    
    file_list2 = file_loading('PHOENIX',teffs=[4250],contains='lines.spec.7')
    for line in file_list2:
        allspec = ('PHOENIX/%(allspec)s' %{'allspec' : line}) 
        
        
        spectrum = np.genfromtxt(allspec, usecols= (0,1,3,4,5,6,8,9,2), dtype='string')
        DtoE(spectrum[:,0])
        DtoE(spectrum[:,1])
        DtoE(spectrum[:,8])
        spectrum = np.array(spectrum, dtype = float)
        
        for g in range(1,9):
            spectrum[:,g] = spectrum[:,g][np.argsort(spectrum[:,0])]
            
        spectrum[:,0] = np.sort(spectrum[:,0])
        spectrum[:,1] = 10.**spectrum[:,1]
        spectrum[:,8] = 10.**spectrum[:,8]
                
        lis = ['Mn','Fe','Mg','Sr','Ca','Ba']
        ionlist =['0','1']
        '''
        while True:
            try:
                
                element = raw_input('element to observe: ')
                assert element  in lis
                ionization = raw_input('ionization of element (0 or 1): ')
                assert ionization in ionlist
                        
                #opticaldepthdiff = input('Delta optical depth threshold for spectral lines: ')
                #print pt.elements.__getitem__(element)
                
                
            
                #print '0'.join(temp)
                #elemention = input('type phoenix code for element to observe: ')        
                #assert elemention in spectrum[:,5]                
                
                
                break                   
            except:
                print 'no'  
                
         '''       
        #element = pt.elements.__getattribute__(element).number
        #temp = [str(element), ionization]
        #
        #elemention = float('0'.join(temp) )   
        elemention = elementcode


        stronglines = spectrum[(spectrum[:,3] - spectrum[:,2] >= opticaldepthdiff) & (spectrum[:,5] == elemention)]
        
        sortlist = np.zeros(shape=(len(stronglines),3))
        for j in range(len(stronglines)):
            if stronglines[j,4] not in sortlist:
                sortlist[j,0] = stronglines[j,4]
                sortlist[j,1] = stronglines[j,5]
                sortlist[j,2] = stronglines[j,1]
          
          
        nozero = (sortlist == 0).sum(1)
        sortlist = sortlist[nozero == 0,:]

          
    return sortlist            
           
#### This function takes an element code from phoenix and turns it into the 
   # corresponding atomic symbol and ionization stage.  
           
def Atom_list(elementnum,ionstage):        
    element = pt.elements.__getitem__(elementnum)
    if ionstage == 0:
        ionize = 'I'
    elif ionstage == 1:
        ionize = 'II'
    else: 
        print 'This thing has an ionization stage other than I or II' 
    return 'PHX: {} - {}'.format(str(element), ionize)
    
def strong_line_list(opticdepth,elementcode):
    '''
        Generates a list containing the strong spectral lines that have
        an optical depth difference greater than or equal to opticdepth.
        elementcode is the 4 number string that phoenix gives to elements
        an their ionization stages.
    '''
    stronglines = phoenix_line_list(opticdepth,elementcode)
    atomicsymbol = stronglines[:,1].tolist()           
    Atoms = np.zeros(shape = (len(atomicsymbol),2))  
    
    for d in range(len(atomicsymbol)):
        Atoms[d,0] = str(atomicsymbol[d])[0:2]
        Atoms[d,1] = str(atomicsymbol[d])[2:4]    
            
    Atoms = Atoms.astype(int)          
    atomlist = []              
    for y in range(0,len(Atoms)):
        atomlist.append( Atom_list(Atoms[y,0],Atoms[y,1] ) )         
         
    return atomlist, stronglines





















