#############################################################
#
#   CONTROL PARAMETERS FOR Spec_graphing.py
#
#
#
#
#
#from LineCompare import Linelist
'''
    Dependency list: (extra libraries needed)
            periodictable 
            PyAstronomy

    README UNDERNEATH!
    
    
'''

Directory = 'PHOENIXnorm/'
outDir = 'ThesisPlots'


#fullspec = [1910, 5640]
linecentre = 5396.17
dx = 0.30
macroturb = 6


alphamode = True
width_calc =   False
rmse_calc =  True
alpha_plot = False
nonalpha_plot = False
relative_plot = False
plot_arcturus = False



''' README:
    _______________________________________________________________________________________________


    Directory: Absolute location of the phoenix .norm.conv files that were created with 
               normalize module.

    outDir:   Absolute location to store saved graphs and calculation output files (RMSE, EW)

    linecentre: wavelength value for the desired spectral line centre.
    
    dx:      width into the continuum, from line centre. Ideally is an equivalent value on both sides.
             Currently have not implemented a dxleft/dxright dichotomy.

    macroturb: (defunct!) Macroturbulence velocity parameter for convolution, this is hard-coded
              into the normalize module. default = 6.
    _______________________________________________________________________________________________

    Various boolean flags for output from spec_graphing.
    
    alphamode: when True, this tells the program to use the alpha
               enhanced spectra when doing RMSE calculations. When False, it will 
               instead use the non-alpha spectra.
               
    width_calc: While True, the program will start doing equivalent width calculations,
                and plot them on a separate figure from the others.
                
    rmse_calc: While True, the program will calculate RMSE values for models in the grid,
               and is affected by the alphamode flag.
               
    alpha_plot: Plots the spectra for alpha-enhanced stars at the linecentre +/- dx.
    nonalpha_plot: Plots the spectra for 'regular' stars at the linecentre +/- dx.
    relative_plot: Plots the relative flux difference between for all models for their 
                   respective lte/nlte pairs.
                   
    plot_arcturus: Plots the observed hi-res arcturus spectrum, generally wanted True when 
                   using alpha/nonalpha_plot = True flags.                   
                   
                   
    _______________________________________________________________________________________________
               
'''