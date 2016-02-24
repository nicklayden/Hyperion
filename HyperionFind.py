#!usr/bin/env python
#Graphing the various spectra from nlte and lte models, with Feh = 0, -2, -4
#For the 4000k Stars.
#With ugriz filter overlay.

# This program is used to graph and analyze the PHOENIX model spectra, and compute
# the nlte and lte relative differences for models in the grid. Spectra are first 
# convolved with a gaussian of width equivalent to a rough observers spectral resolution
# with R ~ 30000, and therefore d(lambda) ~ 0.15 Angstroms. 
#
#
# Input files are the raw PHOENIX '.7' spectrum files, for lte and nlte stars,
# only the first two columns are necessary.
#
# Convolution is done by using a fast Fourier transform convolve function, found
# in the scipy.signal library.
#   
#
'''
NICK:  Need to add way to output only the strongest spectral lines in an array or file.








'''
#
#
#
#
#
#
import warnings
from PyAstronomy import pyasl
import periodictable as pt
import numpy as np
import matplotlib.pyplot as plt
import os
from PhoenixLines import strong_line_list, DtoE, convolve
from LineCompare import load_spec
from Spec_graphing import file_loading

#fehs = ['0.0','0.5','1.0','2.0','4.0','5.0','6.0' ]
#fehs2 = ['0.0','2.0','4.0']
#temps = [4000,4500,5000,5500,6000,6500]
colourstring = ['b','g','r','k','y','c','m' ]
temperature2 = [4125,4250,4375]
#temperature = input('Input model temperature: ')
#fehs3 = input('Input a list of metallicities to analyze, separated by commas: ')




def file_list():
    file_list = list()    #Defines an empty list to put spectra filenames in.
    for file in os.listdir('PHOENIX'):    #My spectra are in this directory.
       #for i in range(0,len(temperature2)):
        if file.endswith('.spec.7') and ('%(teff)s-2.0-%(feh)s' %{'teff' : 4250 ,'feh' : 0.5}) in file and 'alpha' in file and 'lines' not in file:
            file_list.append(file)
    file_list.sort(key = lambda x: x.split('-')[0] )
    file_list.sort(key = lambda x: x.split('-')[1], reverse=True )
    return file_list

file_list = file_list()
file_list2 = list()    #Defines an empty list to put spectra filenames in.
for file in os.listdir('PHOENIX'):    #My spectra are in this directory.
    if file.endswith('lines.spec.7') and 'lte-4250-2.0-0.5' in file:
        file_list2.append(file)

Directory = 'PHOENIX/'
#file_list = file_loading(Directory, contains='.alpha.', notcontains='lines', teffs=[4250])




#VOID ARRAY, USEFUL FOR STORING MULTIPLE TYPES INTO ONE OBJECT
#Putting the important spectral lines file into a multi-type array:
#Then sorting all elements by wavelength.
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

#print min(Linelist['wavelength']), max(Linelist['wavelength'])

#### Converting wavelengths in air to wavelengths in vacuum. ####
Linelist['wavelength'] = pyasl.airtovac2(Linelist['wavelength'])


allelementslist = []
for l in range(0,len(Linelist)):
    if Linelist['element'][l] not in allelementslist:
        allelementslist.append(Linelist['element'][l])

print 'Elements to choose from:'
print allelementslist

while True:
        
    try:
        element = raw_input('Type an element to analyze, or type All for all elements: ')
        if element != 'All':
            assert element in Linelist['element']
            Linelist = Linelist[(Linelist['element'] == element )]
        break  
    except:
        print 'The element entered is not in the important line list! Try again'
    
Linelist = Linelist[(Linelist['ion'] == 'I' )]


possiblestrength = np.arange(0,46,1)
while True:
    try:
        linestrength = input('Strength of lines to display (any # 0-45, increasing strength) :')
        assert linestrength in possiblestrength
        break
    except:
        print 'Value error: Input line strength is not possible.'







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



def phoenix_line_list(opticaldepthdiff):
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
        '''
        while True:
            try:
                elemention = input('type phoenix code for element to observe: ')        
                assert elemention in spectrum[:,5]                
                
                opticaldepthdiff = input('Delta optical depth threshold for spectral lines: ')
                break                   
            except:    
                print 'Element code entered was not found in the Phoenix line file. Try another'
        
        '''
        stronglines = spectrum[(spectrum[:,3] - spectrum[:,2] >= opticaldepthdiff)]
        
        sortlist = np.zeros(shape=(len(stronglines),3))
        for j in range(len(stronglines)):
            if stronglines[j,4] not in sortlist:
                sortlist[j,0] = stronglines[j,4]
                sortlist[j,1] = stronglines[j,5]
                sortlist[j,2] = stronglines[j,1]
          
          
        nozero = (sortlist == 0).sum(1)
        sortlist = sortlist[nozero == 0,:]

          
    return sortlist            

'''
#def strong_line_list():
stronglines = phoenix_line_list(35)
atomicsymbol = stronglines[:,1].tolist()           

#print min(stronglines[:,1]), max(stronglines[:,1])
        
        
Atoms = np.zeros(shape = (len(atomicsymbol),2))  

for d in range(len(atomicsymbol)):
    Atoms[d,0] = str(atomicsymbol[d])[0:2]
    Atoms[d,1] = str(atomicsymbol[d])[2:4]    
        
Atoms = Atoms.astype(int)          
'''       
    #return Atoms   
       
       
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
''' 
atomlist = []              
for y in range(0,len(Atoms)):
    atomlist.append( Atom_list(Atoms[y,0],Atoms[y,1] ) )         
'''


def array(arr):
    return np.array(arr, dtype=float)

     
     
def Spec_Graph(spectra):
    
    
    allspec1 = ('%(allspec)s' %{'allspec' : spectra[1]})   #nlte Spectrum
    allspec2 = ('%(allspec)s' %{'allspec' : spectra[0]})   #lte Spectrum
    
    g = allspec1.split('-')    
    gg = allspec2.split('-')
    
    print 'plotting F_{} - F_{} / F{} relative difference.'.format(g[0],gg[0],gg[0])    
    
    
    
    Spec1 = load_spec(allspec1,'PHOENIX')
    Spec2 = load_spec(allspec2,'PHOENIX')

    
    newgrid = np.arange(3000,12000, 0.01)   #Uniform grid space for each spectra
    
    NewSpec1 = np.interp(newgrid, Spec1[:,0], Spec1[:,1])
    NewSpec2 = np.interp(newgrid, Spec2[:,0], Spec2[:,1])
    
    

    ##### Fast Fourier Transform convolution function #####
    FWHM = 0.015
    smooth1 = convolve(NewSpec1,FWHM,0.01)
    smooth2 = convolve(NewSpec2,FWHM,0.01)

       

    RelativeSpectrum = ( (smooth1 - smooth2)/smooth2 ) * 100.  #Relative fluxes as a percent

    SpecAvg = np.mean(RelativeSpectrum)
    SpecSigma = np.std(RelativeSpectrum)
    #top lines
    
    #plt.hlines(SpecAvg + SpecSigma, 3000,9000, linestyle = '--', colors = colourstring[i], lw = 1.0)
    #plt.hlines(SpecAvg - SpecSigma, 3000,9000, linestyle = '--', colors = colourstring[i], lw = 1.0)
    #plt.hlines(SpecAvg, 3000,9000, linestyle = '--', colors = colourstring[i], lw = 2.0)    
    
    SmoothedSpec = np.array([newgrid, RelativeSpectrum], dtype= float)    
    SmoothedSpec = SmoothedSpec.T
        
    

    SpecFeatures = list()
    for q in range(0,len(SmoothedSpec)):
        if (SmoothedSpec[q,1] >= SpecAvg + SpecSigma) or (SmoothedSpec[q,1] <= SpecAvg - SpecSigma):
            SpecFeatures.append([SmoothedSpec[q,0],SmoothedSpec[q,1]])
        
    SpecFeatures = np.array(SpecFeatures, dtype=float)







    ##### This finds where the important lines are located on the convolved relative spectrum.   
    nothing = []   # Temporary list that is not useful after plotting
    for m in range(0,len(Linelist)):
        
        a = Linelist['wavelength'][m] - 0.005
        b = Linelist['wavelength'][m] + 0.005
        linefile = SmoothedSpec[:,0]
        nothing = np.append(nothing, np.where((linefile >= a) & (linefile <= b)))
    
        
    ##### Array made from y-value at line wavelength on each spectrum.   
    arrowlines = []
    for n in range(0,len(nothing)):
        #arrowlines.append([Linelist['wavelength'][n], SmoothedSpec[nothing[n],1], Linelist['element'][n], Linelist['ion'][n]])
        arrowlines.append([Linelist['wavelength'][n], SmoothedSpec[nothing[n],1]])

    #arrowlines = np.array(arrowlines, dtype= ([('wavelength','f8'),  ('percent','f8'),  ('element','a2'),  ('ion','a2') ]))
    arrowlines = np.array(arrowlines, dtype= float)    
    arrowlines[:,1] = arrowlines[:,1][np.argsort(arrowlines[:,0])]
    arrowlines[:,0] = np.sort(arrowlines[:,0])
    
    
    #plt.plot(arrowlines[:,0], arrowlines[:,1], c= colourstring[i], linewidth= 2.0, 
    #                           label = ' Fe/H = -%(fehs)s ' %{'fehs' : g[3][0:3]}  )
    
    #plt.scatter(arrowlines[:,0], arrowlines[:,1], c=colourstring[i])
    
    '''
    if (g[3][0:3] == '0.5'):
        for r in range(len(Linelist)):
            #if (arrowlines[r,1] >= 10.):
            plt.annotate('%(element)s - %(ion)s, %(wave)s'  
                          %{ 'element' : Linelist['element'][r], 
                                 'ion' : Linelist['ion'][r], 
                                'wave' : Linelist['wavelength'][r]  }  , 
                         xy=(Linelist['wavelength'][r], SmoothedSpec[nothing[r],1]), 
                         xytext=(Linelist['wavelength'][r], SmoothedSpec[nothing[r],1] + 60 ),
                         arrowprops=dict(facecolor='black',arrowstyle ='->'))
                   
    '''               
    
    '''
    for s in range(len(stronglines)):
        plt.annotate(atomlist[s] + ', {}'.format(stronglines[s,0]) , 
                     xy=(stronglines[s,0], 0), 
                     xytext=(stronglines[s,0],80),
                     arrowprops=dict(facecolor='black', arrowstyle= '->'))
    '''                      
   


    
    plt.title('Teff = {}K '' {}$\AA$ Convolution'.format(g[1],FWHM))    
    plt.grid(True)
    plt.ylabel('Relative Flux (%)')
    plt.xlabel('Wavelength ($\AA$)')
    
    
    #Overplotting the locations of important spectral lines in the spectral range.
    '''
    for j in range(len(Linelist)):
        plt.annotate('%(element)s - %(ion)s'  %{ 'element' : Linelist['element'][j], 'ion' : Linelist['ion'][j]  }  , 
                     xy=(Linelist['wavelength'][j], -10), 
                     xytext=(Linelist['wavelength'][j], 80),
                     arrowprops=dict(facecolor='black', headwidth = 0, linewidth = 0)
                )
    '''            
                
    #plt.fill_between(newgrid,y1=SpecAvg + SpecSigma, y2 = SpecAvg-SpecSigma,color = 'y')            
    plt.plot(newgrid, RelativeSpectrum, c= colourstring[i], label = ' Fe/H = -%(fehs)s ' %{'fehs' : g[3][0:3]}, linewidth = 2.0)
    
    #plt.scatter(SpecFeatures[:,0], SpecFeatures[:,1], c = 'r')
    #plt.plot(SpecFeatures[:,0], SpecFeatures[:,1], label = 'Significant lines')
    plt.legend()
    plt.xlim(3000, 9000)
    
    return array(SpecFeatures), array(arrowlines)
   
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
        

################################################################################
#
#       Line ID Files Reading and Plotting Spectral lines by element number
#        and ionization stage.
################################################################################

#def main():

ionization = '0'
element = pt.elements.__getattribute__(element).number
temp = [str(element), ionization]

elemention = float('0'.join(temp) )  


atomlist, stronglines = strong_line_list(linestrength,elemention)

aaron = list()
for i in range(0,len(file_list)/2):
    S = Spec_Graph(file_list)
    SpecFeatures = S[0]
    arrowlines = S[1]
   
 


for k in range(1,len(SpecFeatures)-1):
    Specup =  SpecFeatures[k + 1,1] - SpecFeatures[k,1]      
    Specdown = SpecFeatures[k,1] - SpecFeatures[k - 1,1]
    Specchange = SpecFeatures[k + 1,0] - SpecFeatures[k,0]
    Specchange2 = SpecFeatures[k,0] - SpecFeatures[k-1,0]
    
    if Specup * Specdown < 0.0 and Specchange == Specchange2:
        
        aaron.append((SpecFeatures[k,0],SpecFeatures[k,1]))  
        #plt.annotate('max here', xy=(SpecFeatures[k,0],SpecFeatures[k,1]),
         #                        xytext = (SpecFeatures[k,0], SpecFeatures[k,1] + 60),
          #                      arrowprops=dict(facecolor='red',arrowstyle='->'))
                             



aaron = array(aaron)

phoenix = []
nothing = []
for m in range(0,len(aaron)):
    
        a = aaron[m,0] - 0.01
        b = aaron[m,0] + 0.01
        linefile = stronglines[:,0]
        nothing = np.append(nothing, np.where((linefile >= a) & (linefile <= b)))

nothing = nothing.astype(int)
atomname = []

phx2 = []
#for g in range(0,len(nothing)):
#        phx2.append(stronglines[nothing[g],0])
        
nothing2 = []        
for j in range(0,len(stronglines)):
        
        a = stronglines[j,0] - 0.01
        b = stronglines[j,0] + 0.01
        linefile = aaron[:,0]
        nothing2 = np.append(nothing2, np.where((linefile >= a) & (linefile <= b)))

for q in range(0,len(nothing2)):
        phx2.append((stronglines[nothing[q],0], aaron[nothing2[q],1]))

phx2 = array(phx2)

'''
for s in range(len(phx2)):
    plt.annotate(atomlist[nothing[s]] + ', {}'.format(phx2[s,0]) , 
                 xy=(phx2[s,0], phx2[s,1] ), 
                 xytext=(phx2[s,0],  phx2[s,1] + 100),
                 arrowprops=dict(facecolor='black', arrowstyle= '->'))
      
'''




for q in range(0,len(nothing)):
        g = nothing[q]
        phoenix.append([aaron[g,0],aaron[g,1]])



connections = list()
for i in range(0,len(arrowlines)):
    for j in range(len(stronglines)):
        a = stronglines[j,0] - 0.05
        b = stronglines[j,0] + 0.05
        if arrowlines[i,0] >= a and arrowlines[i,0] <= b:
            connections.append((arrowlines[i,0],arrowlines[i,1]))
            
connections = array(connections)

for i in range(0,len(connections)):
    plt.annotate('{}'.format(connections[i,0]),
                 xy=(connections[i,0], connections[i,1]),
                 xytext=(connections[i,0], connections[i,1] - 10),
                 arrowprops = dict(facecolor='red', arrowstyle='->') )


#regular_spectrum_graph()

#if __name__ == '__main__':
#    main()





    
    
    
    
 
