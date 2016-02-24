###############################################################################
#
#   This module is for normalizing phoenix spectra with a corresponding
#   continuum, and scaling the normalized spectra up or down to the 
#   continuum (1.0).
#   
#   All scaling was done by eye, in the IR. any area with no spectral line
#   will work.
#
import time
import matplotlib.pyplot as plt
import os
import numpy as np
from LineCompare import load_spec, convolve,array,load_norm_spec,arcturus_link
from sklearn.metrics import mean_squared_error




# The path to the phoenix spectra
Directory = 'PHOENIX/'
# The path to the phoenix continuum spectra
contDirectory = 'PHOENIXcont/'        
        
        
        
        
teffs = [4000,4125,4250,4375,4500]

#  These are the values used to scale all spectra to the continuum,
#  in order of how the file_list is loaded.
alphascale = [0.997,1.01,0.997,1.009,0.998,1.012,1,1.0151,1.003,1.022]
noscale = [0.9988,1.0171,0.9990,1.0168,0.9991,1.0166,0.9991,1.0169,0.9990,1.0163]



#teffs=[4000,4500]
fehs = ['0.0','2.0','4.0']  




cont_file = list()
for file in os.listdir(contDirectory):
    for i in range(0,len(teffs)):
        if file.endswith('{}-2.0-0.5.sph.no_rad.ames.cont.7'.format(teffs[i])):
            cont_file.append(file)
            
file_list = list()    #Defines an empty list to put spectra filenames in.
for file in os.listdir(Directory):
   for i in range(0,len(teffs)):
       if file.endswith('spec.7') and ('{}-2.0-{}'.format(teffs[i],0.5)) in file and 'alpha' not in file:
           file_list.append(file)
           
           
#This sorts by model, and then by temperature.     
cont_file.sort(key = lambda x: x.split('-')[1], reverse=True)           
file_list.sort(key = lambda x: x.split('-')[0] )
file_list.sort(key = lambda x: x.split('-')[1], reverse=True )   


def normalizer(modelin,continuum,i):
    model = load_spec(modelin,'PHOENIX')
    cont  = load_spec(continuum,'PHOENIXcont')
    specgrid =np.arange(3000,12000,0.01)  
    
    
    a = modelin.split('-')    
    print 'Normalizing model:', a[0], a[1],'K', ',Continuum scaling:', noscale[i]
    
    modelflux = np.interp(specgrid,model[:,0], model[:,1])    
    contflux  = np.interp(specgrid,cont[:,0], cont[:,1])
    
    normflux = modelflux/contflux
    
    
    convnormflux = convolve(normflux, 6, 0.01)
    convnormflux2 = convolve(convnormflux,2,0.01)

    normalflux = convnormflux2/noscale[i]

    end = len(a[3]) - 2

    normalspec = array((specgrid,normalflux)).T
    np.savetxt('PHOENIXnorm/{}-{}-{}-{}.norm.conv.nonalpha.7'.format(a[0],a[1],a[2],a[3][0:end])
                ,normalspec, fmt='%.18e')





###############################################################################
#
#   MAIN
#
#

start = time.time()


for i in range(0,2*len(teffs)):
    normalizer(file_list[i],cont_file[i//2],i)


end = time.time()

totaltime = end-start

print 'Time taken: {}'.format(totaltime), 'seconds'



















