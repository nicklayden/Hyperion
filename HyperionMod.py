################ Spectral Line Viewer ################
#
#
#  Various functions for analyzing the spectral lines in phoenix spectra
#  and the high resolution arcturus spectrum.
#
#
#
#
from PyAstronomy import pyasl
import numpy as np
import matplotlib.pyplot as plt
import os
from PhoenixLines import DtoE, convolve



def array(arr):
    return np.array(arr, dtype=float)

     

####                   LOAD_SPEC                                            ####    
#    Loads a spectrum into an array, turns D into E from Fortran files,        #
#    then sorts the spectra and turns the flux into cgs flux units.            #
#                                                                              #
####                                                                        ####
def load_spec(specfile,directory):
    spectrum = np.genfromtxt('{}/{}'.format(directory,specfile), usecols=(0,1), dtype= 'string')
    DtoE(spectrum[:,0])    
    DtoE(spectrum[:,1])
    spectrum = np.array(spectrum, dtype=float)
    spectrum[:,1] = spectrum[:,1][np.argsort(spectrum[:,0])]
    spectrum[:,0] = np.sort(spectrum[:,0])
    spectrum[:,1] = 10.0**(spectrum[:,1])
    
    return spectrum
    
def load_norm_spec(specfile,directory):
    '''
        Loads the normalized spectra created with the Normalize_Phoenix program.
    '''
    spectrum = np.genfromtxt('{}/{}'.format(directory,specfile), dtype='string')
    spectrum = np.array(spectrum,dtype=float)
    return spectrum


def regular_spectrum_graph(specfile,continuum):

    


    allspec = ('%(allspec)s' %{'allspec' : specfile})    
    contspec = ('%(cont)s' %{'cont' : continuum})

    scalefactor = [0.995,1.02,0.997,1.02,0.997,1.02]   
    a = specfile.split('-')   

    #Calls in the file from the directory
    spectrum = load_spec(allspec,'PHOENIX')
    conspec = load_spec(contspec, 'PHOENIXcont')    


    newgrid = np.arange(3000,12000,0.01)

    conspecflux = np.interp(newgrid, conspec[:,0], conspec[:,1])
    specflux = np.interp(newgrid, spectrum[:,0], spectrum[:,1])
    macroturb = 6
    spectrograph = 2
    conspecflux2 = convolve(conspecflux, macroturb,0.01)     #Macroturbulence
    conspecflux3 = convolve(conspecflux2, spectrograph,0.01)   #Spectrograph resolution
    
    specflux2 = convolve(specflux, FWHM,0.01)
    specflux3 = convolve(specflux2, FWHM2, 0.01)

    divflux = specflux3/conspecflux3

    plt.title('Observed Arcturus Spectrum \n With 4250K Alpha enhanced model comparison \n Smoothed by convolution with gaussian curve FWHM = 5/150 Angstroms')
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Normalized Flux (unitless)')
    plt.grid(True)
    plt.plot(newgrid, divflux/scalefactor[i], label = '%(type)s, $T_{eff}$ = %(metal)sK'   %{ 'type'  :a[0], 'metal' :a[1]})
    plt.xlim(7149,7150.8)
    plt.legend()
    
    

def Spec_Normalize(model,continuum):
    modspec = '{}'.format(model)    
    contspec = '{}'.format(continuum) 

    a = model.split('-')
    b = continuum.split('-')
    
    #print 'Normalizing {}-{}K model with {}K continuum.'.format(a[0],a[1],b[1])
    #Calls in the file from the directory
    spectrum = load_spec(modspec,'PHOENIX')
    conspec = load_spec(contspec,'PHOENIXcont')    


    newgrid = np.arange(3000,12000,0.01)

    conspecflux = np.interp(newgrid, conspec[:,0], conspec[:,1])
    specflux = np.interp(newgrid, spectrum[:,0], spectrum[:,1])
    
    
    divflux = specflux/conspecflux
    #print 'Normalization complete.'
    
    return newgrid, divflux



def arcturus_link():
    arc_list =list()
    for file in os.listdir('Arcturus'):    
        if 'ar' in file:
            arc_list.append(file)
    
    
    Arcturus = []
    for line in arc_list:
        arcspecfile = ('Arcturus/{}'.format(line) ) 
        ArcSpec = np.genfromtxt(arcspecfile,dtype=float)
        ArcSpec = np.array(ArcSpec,dtype=float)
        Arcturus.append(ArcSpec)
    
    Arc = np.vstack(Arcturus)
    Arc[:,1] = Arc[:,1][np.argsort(Arc[:,0])]
    Arc[:,2] = Arc[:,2][np.argsort(Arc[:,0])]
    Arc[:,3] = Arc[:,3][np.argsort(Arc[:,0])]
    Arc[:,0] = np.sort(Arc[:,0])
    
    Arc[:,0] = pyasl.airtovac2(Arc[:,0])
    np.savetxt('ArcturusSpectrum.txt',Arc)
    
    Arc = Arc[(Arc[:,1] > 0.) & (Arc[:,3] > 0.)]
    
    return Arc



####                  Relative Spectrum                                     ####
#   Creates a relative spectrum for various phoenix models with a high         # 
#   resolution arcturus spectrum.                                              #
#                                                                              #
####                                                                        ####
def Relative_spec(xmin, xmax, i):
    scales = [0.995,1.02,0.997,1.02,0.997,1.02]    
    ModelSpec = array(Spec_Normalize(file_list[i],cont_file[i//2])).T
    Arcturus = arcturus_link()
    a = file_list[i].split('-')
    #Modelspecfile = ('PHOENIX/%(allspec)s' %{'allspec' : model})     
    
    specgrid = np.arange(xmin, xmax, 0.01)    
    Arcturus = Arcturus[(Arcturus[:,0] >= xmin) & (Arcturus[:,0] <= xmax)]
    ModelSpec = ModelSpec[(ModelSpec[:,0] >= xmin) & 
                          (ModelSpec[:,0] <= xmax) ]    
    
    
    
    ModelSpecFlux = np.interp(specgrid, ModelSpec[:,0], ModelSpec[:,1])
    ArcturusFlux  = np.interp(specgrid, Arcturus[:,0]  , Arcturus[:,1])
    
    
    ModelSpec2 = convolve(ModelSpecFlux, 6 , 0.01)


    RelativeSpec = 100.*(ModelSpec2/scales[i] - ArcturusFlux)/ArcturusFlux
    
    Unity = array((specgrid, RelativeSpec )).T
    #Unity = Unity[(Unity[:,1] <= 1000 ) & (Unity[:,1] >= -1000)]
    plt.plot(Unity[:,0], Unity[:,1], label = 'Type: {}, T = {}K'.format(a[0],a[1]))
    plt.ylim(-100,1000)
    plt.grid(True)
    plt.legend()
    plt.ylabel('Relative Flux')
    plt.xlabel('Wavelength, ($\AA$)')
    #Arc = arcturus_link()    
       
        
    #Relative_spec(file_list[3],Arc)    
      
    

Linelist = np.genfromtxt('SpecConvolve/NLTEin.txt', dtype=([('element'   ,'a2'),
                                                           ('ion'        ,'a2'),
                                                           ('wavelength' ,'f8'),
                                                           ('potential'  ,'f8'),
                                                           ('oscillator' ,'f8')]) 
                                                               ,skip_header = 1)
Linelist['element']    = Linelist['element'][np.argsort(Linelist['wavelength'])]    
Linelist['ion']        = Linelist['ion'][np.argsort(Linelist['wavelength'])]    
Linelist['potential']  = Linelist['potential'][np.argsort(Linelist['wavelength'])]    
Linelist['oscillator'] = Linelist['oscillator'][np.argsort(Linelist['wavelength'])]    
Linelist['wavelength'] = np.sort(Linelist['wavelength'])
Linelist['wavelength'] = pyasl.airtovac2(Linelist['wavelength'])

'''


for r in range(len(Linelist)):

    plt.annotate('%(element)s - %(ion)s, %(wave)s'  
                  %{ 'element' : Linelist['element'][r], 
                         'ion' : Linelist['ion'][r], 
                        'wave' : Linelist['wavelength'][r]  }  , 
                 xy=(Linelist['wavelength'][r], 1.0), 
                 xytext=(Linelist['wavelength'][r],  1.05 ),
                 arrowprops=dict(facecolor='black',arrowstyle ='->'))
           
'''     

def main():
    
    temperature2 = [4250,4000,4500]
    
    file_list = list()    #Defines an empty list to put spectra filenames in.
    for file in os.listdir('PHOENIX'):    #My spectra are in this directory.
       for i in range(0,len(temperature2)):
           if file.endswith('.spec.7') and ('%(teff)s-2.0-%(feh)s' %{'teff' : temperature2[i] ,'feh' : 0.5}) in file and 'lines' not in file and 'lte-4250-2.0-0.5.sph.no_rad.ames.spec.7' not in file:
               file_list.append(file)
    file_list.sort(key = lambda x: x.split('-')[1] )
    global file_list
    cont_file = list()
    
    for file in os.listdir('PHOENIXcont'):
        if file.endswith('.cont.7') and '4000-2.0-0.5'in file or '4250-2.0-0.5' in file or '4000-2.0-0.5' in file:
            cont_file.append(file)
    
    global cont_file
    
    
    ltelist= list()
    for i in range(0,len(file_list)):
        if '.sph.no_rad.ames.spec.7' in file_list[i]:        
            ltelist.append(file_list[i])
    
    #nltelist= list()
    for i in range(0,len(file_list)):
        if '.sph.spec.7' in file_list[i] or 'sph.alpha.spec.7' in file_list[i]:        
            ltelist.append(file_list[i])
    
        
    Arc = arcturus_link()
    
    z = -0.00001666
    Arc[:,0] = Arc[:,0]/(z+1)
    
    plt.grid(True)
    plt.plot(Arc[:,0],Arc[:,1], label = 'Arcturus', linewidth = 2.0, c = 'b')
    plt.plot(Arc[:,0],Arc[:,1]/Arc[:,3], label = 'Telluric', linewidth = 2.0, c = 'k')
    plt.plot(Arc[:,0], Arc[:,2], label = 'Sol', linestyle = '--', c='k')
    plt.legend()
    plt.xlim(3000,9000)
    
    colourstring = ['k','c','y','m','pink','orange','brown' ]
    
    
    for i in range(0,len(file_list)):
        Relative_spec(3900,9300,i)
    
    for i in range(0,len(file_list)):
        
        regular_spectrum_graph(file_list[i],cont_file[i//2],colourstring[i],scales[i])
    regular_spectrum_graph(file_list[1],'r',1.016)
    

if __name__ == '__main__':
    main()


