import spytIO
import glob
from OpticalFlow import processProjectionSet
from utils import generateArrayOfRandomNumbers
import os
import sys
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2
from numpy.fft import ifft2 as ifft2
from numpy.fft import fftfreq as fftfreq

from scipy.ndimage.filters import gaussian_filter
from math import pi as pi
from math import floor as floor

import frankoChellappa as fc
import spytlabQT as qt
import corrections

import numpy as np


def saveresults(result,nbPoint):
    txtNbPoint='%4.4d'%nbPoint
    dx = result['dx']
    dy = result['dy']
    phi = result['phi']
    phi2 = result['phi2']
    phi3 = result['phi3']
    gradientNorm = result['gradientNorm']
    spytIO.saveEdf(dx, 'output/dx'+txtNbPoint+'.edf')
    spytIO.saveEdf(dy.real, 'output/dy'+txtNbPoint+'.edf')
    spytIO.saveEdf(phi.real, 'output/phi'+txtNbPoint+'.edf')
    spytIO.saveEdf(phi2.real, 'output/phiLarkinson'+txtNbPoint+'.edf')
    spytIO.saveEdf(phi3.real, 'output/phiKottler'+txtNbPoint+'.edf')
    spytIO.saveEdf(gradientNorm, 'output/gradientNorm'+txtNbPoint+'.edf')


if __name__ == "__main__":
    IrNames=glob.glob('/Users/embrun/Codes/specklematching/Experiments/MoucheSimapAout2017/ref/*.tif')
    IsNames= glob.glob('/Users/embrun/Codes/specklematching/Experiments/MoucheSimapAout2017/sample/*.tif')
    Ir=spytIO.openSeq(IrNames)
    Is= spytIO.openSeq(IsNames)

    nbPointsTotal=len(Is)
    for nbPoints in range(1,nbPointsTotal+1):
        pointsToTake=generateArrayOfRandomNumbers(nbPoints,0,nbPointsTotal-1)
        IrToTake=Ir[pointsToTake]
        IsToTake= Is[pointsToTake]
        result = processProjectionSet(IsToTake, IrToTake)
        saveresults(result,nbPoints)
