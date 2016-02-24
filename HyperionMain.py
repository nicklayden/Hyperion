##
#
#       PROGRAM:   HYPERION
#   
#       Hyperion: Titan god of light. Father of Helios (the sun)
#       -Greek Mythology
#
#
##

import os
import numpy as np
import matplotlib.pyplot as plt
from LineCompare import *
from PhoenixLines import convolve, file_loading
import time
from sklearn.metrics import mean_squared_error
from spec_control import *
import warnings


################################################################################
#
#   SETUP PARAMETERS FOR THIS MODULE
#       file paths, model temperatures to observe, Type of spectrum to run
#
#   Can sort through metallicity ranges or temperature ranges. Alpha enhanced
#   spectra or not.
#   Can handle observed spectra as well (Arcturus & others)
#
#       Directory: Path on the local machine to the Phoenix spectrum files.
#   contDirectory: Path on the local machine to the Phoenix continuum files.
#           teffs: The list of model temperatures you want to use in comparisons.
#            fehs: (optional) list of metallicities to compare.
#       file_list: Total list of all spectra to use.
#       cont_file: Total list of all continuum spectra for use in normalization.
#
################################################################################
'''
# The path to the phoenix spectra
Directory = 'PHOENIXnorm/'
# The path to the phoenix continuum spectra
contDirectory = 'PHOENIXcont/'        
        
teffs = [4000,4125,4250,4375,4500]

#teffs=[4000,4500]
fehs = [0.5,1.0]  
'''


'''
cont_file = list()
for file in os.listdir(contDirectory):
    for i in range(0,len(teffs)):
        if file.endswith('{}-2.0-0.5.sph.no_rad.ames.cont.7'.format(teffs[i])):
            cont_file.append(file)
            
file_list = list()    #Defines an empty list to put spectra filenames in.
for file in os.listdir(Directory):    #My spectra are in this directory.
   for i in range(0,len(teffs)):
       if file.endswith('norm.conv.7') and ('{}-2.0-{}'.format(teffs[i],0.5)) in file:
           file_list.append(file)

           
nonalpha_list = list()
for file in os.listdir(Directory):
    for i in range(0,len(teffs)):
        if file.endswith('norm.conv.nonalpha.7') and ('{}-2.0-{}'.format(teffs[i],0.5)) in file:
            nonalpha_list.append(file)
           
           
#This sorts by model, and then by temperature.     
cont_file.sort(key = lambda x: x.split('-')[1], reverse=True)           
file_list.sort(key = lambda x: x.split('-')[0] )
file_list.sort(key = lambda x: x.split('-')[1], reverse=True )   
nonalpha_list.sort(key = lambda x: x.split('-')[0] )
nonalpha_list.sort(key = lambda x: x.split('-')[1], reverse=True )   
'''


#######################################################################################################
#
#
#       FILE LOADER:
#                   Loads any phx files into a list for use in all functions here.
#
#
#
'''
def file_loading(directory, extension=None,teffs=None,contains=None,metallicity=None,notcontains=None):
    
        Loads PHX files into a list for all of my other modules to find spectra.
        Default extension is .7, for ALL phx files. 
        teffs is the list or array or temperature values to scan through, and the metallicity
        is default 0.5. Since these are all metal-poor spectra files, the negative is implied 
        on all metallicity values.
        
        contains    : if file contains this ____ string, include it.
        notcontains : if file contains this ____ string, do NOT include it.(Ambiguity with alpha/nonalpha)   
        teffs       : LIST object of effective temperature values to look for. Cannot input single float values.
        metallicity : singular metallicity value to look for. Not implementing metallicity range here.
        extension   : phx file extensions are long and specific, easy way to find exactly what you want with this.
                      default extension finds any .7 files.
        
        
        teffs:    Note, although you cannot input a single value like 'teffs=4250', you CAN input a single valued 
                  list object, i.e.   'teffs=[4250]'. This should work correctly.
        
    
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
'''    
    



################################################################################
# 
#   Root Mean Squared Error Calculation:
#   
#   Input:     modelin        :     Normalized model atmosphere 
#              observed       :     Observed spectrum
#              linecentre     :     Centre of desired spectral line. (Angstroms) 
#              dx             :     Uniform wavelength spacing for *both* spectra
#
def RMSE(modelin,observed,linecentre,dx):
 
    a = modelin.split('-')    
    specgrid = np.arange(3000,12000,0.01)
    model = load_norm_spec(modelin,Directory)
    ArcFlux   = np.interp(specgrid, observed[:,0], observed[:,1])
    TotalModel = model

    TotalArcturus = array((specgrid,ArcFlux)).T

    TotalModel = TotalModel[(TotalModel[:,0] <= linecentre + dx) & (TotalModel[:,0] >= linecentre - dx)]
    TotalArcturus = TotalArcturus[(TotalArcturus[:,0] <= linecentre + dx) & (TotalArcturus[:,0] >= linecentre - dx)]
    
    #print 'Calculating RMS error between Model={}, Teff={}K, and Arcturus on linecentre = {} +/- {} Angstroms'.format(a[0], a[1], linecentre, dx)
    
    error = np.sqrt(mean_squared_error(TotalArcturus[:,1],TotalModel[:,1]))
    
    print 'RMSE= {}, Model: {}-{}K'.format(error,a[0],a[1])

    
    return error, a[0], a[1], linecentre


################################################################################
#
#  Relative_spec:
#       Uses input spectra files to create a relative difference spectrum.
#       Formula = (F_nlte - F_lte)/F_lte
#
#
#       ARCTURUS SPECTRUM IS BROKEN BETWEEN 7500A AND 8000A APPROXIMATELY.
#
#
#
def Relative_spec(nltespec, ltespec):
        
    '''
            FUNCTION CURRENTLY DOES NOT WORK PROPERLY!!! 
            THIS IS AN ARTIFACT OF ME CREATING THE .norm.7 FILES
            TO INCREASE PROGRAM SPEED.


    '''        
        
        
        
    model = array(Spec_Normalize(nltespec,ltespec)).T
    

    a = nltespec.split('-')
    #b = ltespec.split('-')
    #Modelspecfile = ('PHOENIX/%(allspec)s' %{'allspec' : model}) 
    Arcturusspec = ltespec
    
    specgrid =np.arange(3000,12000,0.01)    
    
    ModelSpec = model
    Arcturus = Arcturusspec

    ModelSpecFlux = np.interp(specgrid, ModelSpec[:,0], ModelSpec[:,1])
    ArcturusFlux  = np.interp(specgrid, Arcturus[:,0]  , Arcturus[:,1])
    
    
        
    print 'Plotting {} model,  Teff= {}K '.format(a[0],a[1]), '[Fe/H]= -{}'.format(a[3][0:3])  ,'iter=', i
    
    
    
    ModelSpec2 = convolve(ModelSpecFlux, macroturb , 0.01) #Macroturbulence
    ModelSpec3 = convolve(ModelSpec2,    2 , 0.01)    #Spectrograph

    RelativeSpec = 100.*(ModelSpec3 - ArcturusFlux)/ArcturusFlux

    Unity = array((specgrid, RelativeSpec )).T
    
    Unity2 = Unity[(Unity[:,0] <= 7300 )]    
    
    
    
    SpecSigma = np.std(Unity2[:,1])       
    SpecMean = np.mean(Unity[:,1])
    print 'Spectrum Mean= ',SpecMean
    print 'Standard Deviation= ',SpecSigma      
    
    
    plt.plot(Unity2[:,0], Unity2[:,1], label = '[Teff] = {}K, model={}, cTeff= {}K'.format(a[1],a[0],b[1]))
    #plt.ylim(-100,1000)
    plt.xlim(linecentre - dx, linecentre + dx)
    plt.grid(True)
    plt.legend()
    plt.ylabel('$(F_{NLTE} - F_{OBS})/F_{OBS}$ (%)')
    plt.xlabel('Wavelength, ($\AA$)')
    plt.title('Relative Spectrum Differences \n Between Model and Observed Spectra for Arcturus \n ')
  

    return RelativeSpec

################################################################################
#
#   Spectrum Plotter:
#       Uses input spectra files to create a plot of the spectra.
#       
#
def Spec_Plot(specfile,linecentre,dx,ls=None):
    if ls ==None:
        ls = '-'
    element = return_element(linecentre)
    spectrum = load_norm_spec(specfile,Directory)

    a = specfile.split('-')
    print 'Plotting model {}, Teff= {}K, [Fe/H]= -{}'.format(a[0],a[1],a[3][0:3])
    '''
    plt.title('Observed & Synthetic Spectra \n'
              'Element:  {}     \n'  
              'Spectral line:  {} +/-{} \n'
              'With 6 km/s Macroturbulence Broadening'.format(element,linecentre,dx))
    '''
    plt.title(r'Atomic Species: {} : {}$\rm \AA$'.format(element,linecentre))
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Normalized Flux')
    plt.grid(True)
    plt.plot(spectrum[:,0], spectrum[:,1], label = '%(type)s, $T_{eff}$ = %(metal)sK'   %{ 'type'  :a[0], 'metal' :a[1]}, ls=ls)
    plt.xlim(linecentre - dx, linecentre + dx)
    #plt.legend(loc=2, bbox_to_anchor=(1.05,1))
    #plt.legend()
    
    








def Spec_Plot2(specfile,linecentre,dx):
 
    spectrum = load_norm_spec(specfile,Directory)

    a = specfile.split('-')
    print 'Plotting model {}, Teff= {}K, [Fe/H]= -{}'.format(a[0],a[1],a[3][0:3])

    plt.title('Observed & Synthetic Spectra \n'
              'Spectral line: {} +/- {}   \n'
              'With 6 km/s Macroturbulence Broadening'.format(linecentre,dx))
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Normalized Flux')
    plt.grid(True)
    plt.plot(spectrum[:,0], spectrum[:,1], label = '%(type)s, $T_{eff}$ = %(metal)sK'   %{ 'type'  :a[0], 'metal' :a[1]}, linestyle ='--')
    #plt.legend(loc=2, bbox_to_anchor=(1.05,1))
    plt.xlim(linecentre - dx, linecentre + dx)

    




################################################################################
#
#   Equivalent width calculator
#   Projects a line across a spectral line and integrates, to determine
#   the equivalent width of that spectral line.
#
#
def equiv_width(modelin,linecentre,dx,arcturus=None):   
    '''
        Default equivalent width selection is False.
        Just to make commands neater later on.
        
        Interpolates the spectra to a uniform grid, finds the specified line centre
        and line width, uses the end points of the spectral line to project a straight
        line across the range of values, then does a simple integration to calculate 
        the width.
        
        specgrid  : The uniform spacing grid across the spectrum range 3000A to 12000A        
        xmin,ymin : Point at the far left of the spectral line.
        xmax,ymax : Point at the far right "   "
        proj_x    : Range of x values from blue to red on any spectral line.
        proj_y    : Projection of the slope across the top of the spectral line.
        equivalent_width : Self-explanatory.
        
        
    '''
    if arcturus == True:
        model = modelin
        b = 'Arcturus'
        c = None

    elif arcturus == None:    
        model = load_norm_spec(modelin,Directory)
        a = modelin.split('-') 
        b = a[0]
        c = a[1]

    
    specgrid =np.arange(3000,12000,0.01)  

    modelflux = np.interp(specgrid, model[:,0], model[:,1])

    TotalModel = array((specgrid,modelflux)).T
    TotalModel = TotalModel[(TotalModel[:,0] <= linecentre + dx) & (TotalModel[:,0] >= linecentre - dx)]

    
    # (y - ymin) = slope(x - xmin) 
    xmin, xmax = TotalModel[0,0], TotalModel[-1,0]    
    ymin, ymax = TotalModel[0,1], TotalModel[-1,1]
    
    proj_x = np.arange(xmin,xmax, 0.01)
    slope  = (ymax-ymin)/(xmax-xmin)
     
    # slope intercept form to calculate projection line. 
    proj_y = slope*proj_x + (ymin - slope*xmin)
    
    
    #Uncomment this to see an example of the line projected across the spectral line.
    #plt.grid(True)
    #plt.plot(TotalModel[:,0], TotalModel[:,1], c='b')
    #plt.plot( proj_x   , proj_y     , c='r')
    #plt.legend()

    equivalent_width = ( (proj_y - TotalModel[:,1])*0.01 ).sum()  
    
    if arcturus == False:
        print 'W= {}, model: {} {}K'.format(equivalent_width, b,c)
    elif arcturus == True:
        print 'W= {},High-res Arcturus spectrum'.format(equivalent_width)

    return equivalent_width, b,c, linecentre







################################################################################
#
#   Transfigurator
#   Changes the generated typle from rmse and W calculations into arrays,
#   splitting nlte and lte model calculations into two arrays. Format is:
#   
#                           array = (tupletochange, Teffs)
#

def Transfigurator(sometuple,flag=None,width=None):
    '''
        Default flag is set to do RMSE arrangements.
        Changes the tuple generated from rmse and W calculations into arrays.
        Splitting the nlte and lte arrays apart for graphing later on.
    '''
    if flag==None:
        flag = 'RMSE'
    if width==None:
        width = 0        
        
    ltearray = list()
    nltearray = list()
    
    for i in range(0,len(sometuple)):
        if sometuple[i][1] == 'lte':
            ltearray.append((sometuple[i][0],sometuple[i][2]))
        elif sometuple[i][1] == 'nlte':
            nltearray.append((sometuple[i][0],sometuple[i][2]))


    ltearray = array(ltearray)
    nltearray = array(nltearray)
    
    if flag == ('W' or 'w'):
        print 'Equivalent width selection: Subtracting Arcturus Width from model widths.'
        ltearray[:,0] = ltearray[:,0] - arcwidth[0]
        nltearray[:,0] = nltearray[:,0] - arcwidth[0]
    elif flag == 'RMSE':
        print 'RMSE selection: Not doing anything fancy.'

    return ltearray,nltearray


def return_element(linecentre):
    '''
        Takes the linecentre wavelength for a given spectral line, and finds the offending 
        element in the special Linelist.
        
        
    
    '''    
    linecentreMatch = list()
    for i in range(len(Linelist)):
        selectLine = Linelist['wavelength'][i]
        a = 0.01
        if linecentre <= selectLine + a and linecentre >= selectLine - a:
            linecentreMatch.append((Linelist['element'][i],Linelist['ion'][i],Linelist['wavelength'][i] )) 
                   
    if len(linecentreMatch) > 1:
        warnings.warn('Warning: Input spectral line matched multiple tags in the list, this could be an ambiguity.')

            
    return str(linecentreMatch[0][0] + '-'  + linecentreMatch[0][1] )



def plot_rmse(rmselist, select=None):
        '''
        select commands = ['nlte','lte','both']
        
        '''
        line = str(linecentre).split('.')
        delta= str(dx).split('.')
        if select == None:
            select = 'both'
        if alphamode == True:
            alphaflag = 'alpha'
        else:
            alphaflag = 'nonalpha'
        ltermse, nltermse = Transfigurator(rmselist,'RMSE')
        plt.figure(2)
        
        if select == 'nlte':
            plt.plot(nltermse[:,1], nltermse[:,0], c='b', label='nlte RMSE')
        elif select == 'lte':
            plt.plot(ltermse[:,1], ltermse[:,0], c='g', label='lte RMSE')
        elif select == 'both':
            plt.plot(ltermse[:,1], ltermse[:,0], c='g', label='lte RMSE')
            plt.plot(nltermse[:,1], nltermse[:,0], c='b', label='nlte RMSE')
        
        plt.ylabel(r'RMSE')
        plt.xlabel(r' Effective Temperature (K)')
        plt.title(r'RMSE for {} - {}$\rm \AA$'.format(return_element(linecentre),linecentre) )
        #plt.title('Root Mean Squared Error Calculation, element: {} '
        #'\n at {} +/- {} $\AA$ \n For [Fe/H] = -0.5 models'.format(return_element(linecentre),linecentre,dx))
        plt.grid(True)
        plt.legend()
        plt.savefig('{}/{}-{}-{}-{}-{}-rmse.png'.format(outDir,line[0],line[1],delta[0],delta[1],alphaflag))

def rmse_output(rmselist,linecentre,dx):
    
    test1 = []
    test2 = []
    intermed = []
    intermed2 = []
    type1 = 'nlte'
    type2 = ' lte'
    for i in range(len(rmselist)):
        if type1 in rmselist[i]:
            intermed.append(rmselist[i][0])
            
    for i in range(len(rmselist)):
        if type1 not in rmselist[i]:
            intermed2.append(rmselist[i][0])
            
                    
    #header1 = list(str('RMSE Calculations for models in grid T = 4000->4500, dT = 125'))
    #header2 = list(str('element lambda dlambda alpha type 4000 4125 4250 4375 4500'))
    
    test1 = [return_element(linecentre),linecentre,dx,alphamode,type1, '%.4f' % intermed[0], 
                                                                 '%.4f' % intermed[1],
                                                                 '%.4f' % intermed[2],
                                                                 '%.4f' % intermed[3],
                                                                 '%.4f' % intermed[4]]
    
    test2 = [return_element(linecentre),linecentre,dx,alphamode,type2,'%.4f' % intermed2[0], 
                                                                 '%.4f' % intermed2[1],
                                                                 '%.4f' % intermed2[2],
                                                                 '%.4f' % intermed2[3],
                                                                 '%.4f' % intermed2[4]]
        
    values = []
    values.append(test1)
    values.append(test2)
    #Stack the lists for lte and nlte vertically, giving a list  that has each line
    #one of the outputs of the routine for each model type.
    valtest = np.vstack(values)
    
    #np.savetxt('rmsetest.dat',valtest, fmt='%s')
    
    
    
    opentest = np.genfromtxt('ThesisPlots/rmsevalues.dat', dtype=str)
    opentest1 = list(opentest)
    opentest1.append(valtest)
    opentest2 = np.vstack(opentest1)
    
    np.savetxt('ThesisPlots/rmsevalues.dat',opentest2, fmt='%s')

    return opentest2



################################################################################
#
#   MAIN BODY
#
################################################################################

def main():
    
    
    Timestart = time.time()
    
    global Arc,arcwidth,file_list,rmselist
    file_list = file_loading(Directory, contains='.alpha.')
    nonalpha_list = file_loading(Directory, contains='nonalpha')
    Arc = arcturus_link()


    if plot_arcturus == True:
        plt.figure(1)
        plt.plot(Arc[:,0], Arc[:,1], c='b', lw=3.5,label='Arcturus')
        #plt.legend()    
    
         
    if width_calc == True:    
        arcwidth = equiv_width(Arc,linecentre,dx,True )
        widthlist = list()
        if alphamode == True:
            for i in range(len(file_list)):
                widthlist.append(equiv_width(file_list[i],linecentre,dx))
        elif alphamode == False:
            for i in range(len(nonalpha_list)):
                widthlist.append(equiv_width(nonalpha_list[i],linecentre,dx))

    
    if rmse_calc == True:
        rmselist = list()
        if alphamode == True:
            for i in range(len(file_list)):
                rmselist.append(RMSE(file_list[i], Arc, linecentre,dx))

            plot_rmse(rmselist)
            rmse_output(rmselist,linecentre,dx)
            
        elif alphamode == False:
            for i in range(len(file_list)):
                rmselist.append(RMSE(nonalpha_list[i], Arc,linecentre,dx))

            plot_rmse(rmselist)
            rmse_output(rmselist,linecentre,dx)
    
    if alpha_plot == True:
        for i in range(len(file_list)):
            Spec_Plot(file_list[i],linecentre,dx+0.2)
    
    if nonalpha_plot == True:
        for i in range(len(file_list)):
            Spec_Plot(nonalpha_list[i],linecentre,dx+0.2,ls='--')
    
    if relative_plot == True:
        #for i in range(len(file_list)):
             Relative_spec(file_list[4],Arc)
    
   
    
    '''
    
    if rmse_calc == True:
        
        
        ltermse, nltermse = Transfigurator(rmselist,'RMSE')
        plt.figure(2)
        plt.plot(ltermse[:,1],ltermse[:,0], c='g', label='lte RMSE')
        plt.plot(nltermse[:,1], nltermse[:,0], c='b', label='nlte RMSE')
    
        plt.ylabel('RMSE (Flux units)')
        plt.xlabel('$T_{eff}$(K)')
        plt.title('Root Mean Squared Error Calculation, element: {} '
        '\n at {} +/- {} $\AA$ \n For [Fe/H] = -0.5 models'.format(return_element(linecentre),linecentre,dx))
        plt.grid(True)
        plt.legend()
    '''
    
    
    if width_calc == True:
            
        ltewidth, nltewidth = Transfigurator(widthlist,'W',width= arcwidth)
        plt.figure(3)    
        plt.plot(ltewidth[:,1],1000*ltewidth[:,0], c='g', label='lte W')
        plt.plot(nltewidth[:,1], 1000*nltewidth[:,0], c='b', label='nlte W')
        plt.ylabel('$W_{model} - W_{Arcturus}$, $m\AA$')
        plt.xlabel('$T_{eff}$(K)')
        plt.title('Difference in Equivalent Width \n From Model to Arcturus'
                            ' \n at {} +/- {} $\AA$ '.format(linecentre,dx))
        plt.grid(True)
        plt.legend()                    
    
    
    
    
    Timeend = time.time()
    Timelength = Timeend - Timestart
    print 'All calculations finished, time taken =',Timelength,' seconds.'


if __name__ == '__main__':
    main()




