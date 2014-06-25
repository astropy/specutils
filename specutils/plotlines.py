import numpy as np
from scipy import signal
from scipy import interpolate
import pylab as pl


def oplotlines(bandname=None,linelist=None,angstrom=False,color='k',xlim=None,ylim=None,label=True,size=14):
    '''
    Overplots lines on top of a spectrum. Can select between different
    filters.  If there is a currently open plot, will try to detect
    the wavelength range appropriate for plotting.

    Optional Keywords:
    bandname -- name of the filter to plot (default: None)
    linelist -- an array of line positions to plot. This will overide
    the lines that this function uses (default: None).
    angstrom -- by default the lines are in microns. Set this keyword to
                switch to angstroms

    HISTORY:
    2014-02-19 -- T. Do
    '''

    if linelist is None:
        # hydrogen lines in microns
        hlines = np.array([1.00521, 1.09411, 1.28216, 1.52647, 1.58848, 1.61137,1.64117,1.68111,1.81791,1.87561,1.94509])
        hlinesNames = [r'HI (Pa $\delta$)', r'HI (Pa $\gamma$)',r'HI (Pa $\beta$)',
                       'HI', 'HI', 'HI', 'HI','HI','HI',r'HI (Pa $\alpha$)',
                       'HI (Br $\delta$)']
        helines = np.array([1.01264, 1.0833,1.16296, 1.16764,1.69230,1.70076])
        helinesNames = ['HeII', 'HeI','HeII','HeII','HeII','HeI','HeII']

    # should have a plot already so get the current axes
    ax = pl.gca()
    if xlim is None:
        xlim = ax.get_xlim()
    if ylim is None:
        ylim = ax.get_ylim()
    if angstrom:
        hlines = hlines*1e4
        helines = helines*1e4

    totalLines = np.append(hlines, helines)
    totalNames = np.append(hlinesNames, helinesNames)

    goodRange = np.where((totalLines >= xlim[0]) & (totalLines <= xlim[1]))[0]

    if len(goodRange) > 0:
        for i in goodRange:
            pl.plot([totalLines[i],totalLines[i]],ylim,color)
            if label:
                pl.text(totalLines[i],(ylim[1]-ylim[0])*0.05+ylim[0],totalNames[i],rotation='vertical',size=size,va='bottom')

def oplotskylines(band = 'H', linelist = None, xlim = None, ylim = None, color='k',angstrom=False):
    '''
    Plot OH skylines
    '''
    
    if band == 'Y':
        lines = np .array([
              9793.6294 , 9874.84889 , 9897.54143 , 9917.43821 , 10015.6207 ,
             10028.0978 , 10046.7027 , 10085.1622 , 10106.4478 , 10126.8684 ,
              10174.623 , 10192.4683 , 10213.6107 , 10289.3707 , 10298.7496 ,
             10312.3406 , 10350.3153 , 10375.6394 , 10399.0957 , 10421.1394 ,
             10453.2888 , 10471.829 , 10512.1022 , 10527.7948 , 10575.5123 ,
             10588.6942 , 10731.6768 , 10753.9758 , 10774.9474 , 10834.1592 ,
             10844.6328 , 10859.5264 , 10898.7224 , 10926.3765 , 10951.2749 ,
             10975.3784 , 11029.8517 , 11072.4773 , 11090.083 , 11140.9467 ,
             11156.0366 , ])


    if band == 'H':
        lines = np.array([
                14605.0225 , 14664.9975 , 14698.7767 , 14740.3346 , 14783.7537 ,
                14833.029 , 14864.3219 , 14887.5334 , 14931.8767 , 15055.3754 ,
                15088.2599 , 15187.1554 , 15240.922 , 15287.7652 ,
                15332.3843 , 15395.3014 , 15432.1242 , 15570.0593 , 15597.6252 ,
                15631.4697 , 15655.3049 , 15702.5101 , 15833.0432 , 15848.0556 ,
                15869.3672 , 15972.6151 , 16030.8077 , 16079.6529 , 16128.6053 ,
                16194.6497 , 16235.3623 , 16317.0572 , 16351.2684 , 16388.4977 ,
                16442.2868 , 16477.849 , 16502.395 , 16553.6288 , 16610.807 ,
                16692.2366 , 16708.8296 , 16732.6568 , 16840.538 , 16903.7002 ,
                16955.0726 , 17008.6989 , 17078.3519 , 17123.5694 , 17210.579 ,
                17248.5646 , 17282.8514 , 17330.8089 , 17386.0403 , 17427.0418 ,
                17449.9205 , 17505.7497 , 17653.0464 , 17671.843 , 17698.7879 ,
                17811.3826 , 17880.341 ])
  

    if band == 'J':
        lines = np.array([
            11538.7582 , 11591.7013 , 11627.8446 , 11650.7735 , 11696.3379 ,
            11716.2294 , 11788.0779 , 11866.4924 , 11988.5382 , 12007.0419 ,
            12030.7863 , 12122.4957 , 12135.8356 , 12154.9582 , 12196.3557 ,
            12229.2777 , 12257.7632 , 12286.964 , 12325.9549 , 12351.5321 ,
            12400.8893 , 12423.349 , 12482.8503 , 12502.43 , 12589.2998 ,
            12782.9052 , 12834.5202 , 12905.5773 , 12921.1364 , 12943.1311 ,
            12985.5595 , 13021.6447 , 13052.818 , 13085.2604 , 13127.8037 ,
            13156.9911 , 13210.6977 , 13236.5414 , 13301.9624 , 13324.3509 ,
             13421.579])

    if band == 'K':
        #drop: 19751.3895, 19736.4099
        print "the lines"
        lines = np.array([
        19518.4784 , 19593.2626 , 19618.5719 , 19642.4493 , 19678.046 ,
        19701.6455 , 19771.9063 , 19839.7764 ,
        20008.0235 , 20193.1799 , 20275.9409 , 20339.697 , 20412.7192 ,
         20499.237 , 20563.6072 , 20729.032 , 20860.2122 , 20909.5976 ,
        21176.5323 , 21249.5368 , 21279.1406 , 21507.1875 , 21537.4185 ,
        21580.5093 , 21711.1235 , 21802.2757 , 21873.507 , 21955.6857 ,
        22125.4484 , 22312.8204 , 22460.4183 , 22517.9267 , 22690.1765 ,
        22742.1907 , 22985.9156, 23914.55, 24041.62])

        ## for arc lines
        ##lines = np.array([20008.2,
        ##                  20275.9,
        ##                  20412.7,
        ##                  20563.6,
        ##                  20729.0,
        ##                  21802.2,
        ##                  21955.6,
        ##                  22125.5,
        ##                  22312.7])

    # default to microns unless angstrom is used. This is to make this consistent with oplotlines
    if angstrom is False:
        lines = lines/1e4

    if linelist is not None:
        lines = np.loadtxt(linelist)
        print 'Loaded Line file %s ' % linelist
    lines = np.array(lines)

    # should have a plot already so get the current axes
    ax = pl.gca()
    if xlim is None:
        xlim = ax.get_xlim()
    if ylim is None:
        ylim = ax.get_ylim()

    goodRange = np.where((lines >= xlim[0]) & (lines <= xlim[1]))[0]
    if len(goodRange) > 0:
        for i in goodRange:
            pl.plot([lines[i],lines[i]],ylim,color)
    
