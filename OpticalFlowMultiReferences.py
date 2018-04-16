import spytIO
import glob
import os
import sys
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2
from numpy.fft import ifft2 as ifft2

from math import pi as pi
from math import floor as floor

import frankoChellappa as fc
import spytlabQT as qt

import numpy as np
import corrections as corr
import fastTomoExperiment as esrfTomo





def derivativesByOpticalflow(intensityImage,derivative,alpha=0,sig_scale=0):

    Nx, Ny = derivative.shape
    # fourier transfomm of the derivative and shift low frequencies to the centre
    ftdI = fftshift(fft2(derivative))
    # calculate frequencies
    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx)


    #building filters


    sigmaX = dqx / 1. * np.power(sig_scale,2)
    sigmaY = dqy / 1. * np.power(sig_scale,2)
    #sigmaX=sig_scale
    #sigmaY = sig_scale

    g = np.exp(-(((Qx)**2) / 2. / sigmaX + ((Qy)**2) / 2. / sigmaY))
    #g = np.exp(-(((np.power(Qx, 2)) / 2) / sigmaX + ((np.power(Qy, 2)) / 2) / sigmaY))
    beta = 1 - g;

    i = complex(0, 1)
    # fourier filters
    ftfiltX = (i * Qx / ((Qx**2 + Qy**2 + alpha))*beta)
    #ftfiltX = (1 * i * Qx / (np.power(Qx, 2) + np.power(Qy, 2) + alpha) * beta)
    ftfiltX[np.isnan(ftfiltX)] = 0

    ftfiltY = (i* Qy/ ((Qx**2 + Qy**2 + alpha))*beta)
    #ftfiltY = (1 * i * Qy / (np.power(Qx, 2) + np.power(Qy, 2) + alpha) * beta)
    ftfiltY[np.isnan(ftfiltY)] = 0

    # output calculation
    dImX = 1. / intensityImage * ifft2(ifftshift(ftfiltX * ftdI))
    dImY = 1. / intensityImage * ifft2(ifftshift(ftfiltY * ftdI))

    return dImX.real,dImY.real





def kottler(dX,dY):
    print('kottler')
    i = complex(0, 1)
    Nx, Ny = dX.shape
    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx)

    polarAngle = np.arctan2(Qy, Qx)
    ftphi = fftshift(fft2(dX + i * dY))*np.exp(i*polarAngle)
    ftphi[np.isnan(ftphi)] = 0
    phi3 = ifft2(fftshift(ftphi))
    return phi3



def LarkinAnissonSheppard(dx,dy,alpha =0 ,sigma=0):
    Nx, Ny = dx.shape
    i = complex(0, 1)
    G= dx + i*dy
    # fourier transfomm of the G function
    fourrierOfG = fftshift(fft2(G))


    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx)

    ftfilt = 1 / (i*Qx - Qy)
    ftfilt[np.isnan(ftfilt)] = 0
    phi=ifft2(ifftshift(ftfilt*fourrierOfG))
    phi=np.absolute(phi.real)
    return phi



def processOneProjection(Is,Ir):
    sigma = 0.95
    alpha = 0

    dI = (Is - Ir * (np.mean(Is) / np.mean(Ir)))
    dx, dy = derivativesByOpticalflow(Is, dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3}



def processProjectionSet(Is,Ir):
    sigma = 1
    alpha = 0

    subImage=Is-Ir
    subImage=np.mean(subImage,axis=0)

    dI = (subImage * (np.mean(Is) / np.mean(Ir)))
    dx, dy = derivativesByOpticalflow(np.mean(Is,axis=0), dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3}




if __name__ == "__main__":

    referenceFolder= '/Users/embrun/Codes/specklematching/Experiments/MoucheSimapAout2017/ref/'
    sampleFolder = '/Users/embrun/Codes/specklematching/Experiments/MoucheSimapAout2017/sample/'
    app = qt.QApplication(sys.argv)
    #app.connect(app, qt.SIGNAL("lastWindowClosed()"), app, qt.SLOT("quit()"))

    referenceFiles= (qt.QFileDialog.getOpenFileNames(None, 'Open a set of Images', '/','Image files (*.edf *.tif *.tiff)'))
    directory=os.path.dirname(str(referenceFiles[0]))
    print referenceFiles[0]
    sampleFiles = (qt.QFileDialog.getOpenFileNames(None, 'Open a set of Images', directory,'Image files (*.edf *.tif *.tiff)'))

    saveFolder = (qt.QFileDialog.getExistingDirectory(None,'Open directory to save the images',directory ))
    saveFolder=str(saveFolder)

    Ir=spytIO.openImage(str(referenceFiles[0]))
    Is=spytIO.openImage(str(sampleFiles[0]))
    result=processOneProjection(Is,Ir)

    dx = result['dx']
    dy = result['dy']
    phi = result['phi']
    phi2 = result['phi2']
    phi3 = result['phi3']
    spytIO.saveEdf(dx, saveFolder+'/dx.edf')
    spytIO.saveEdf(dy.real, saveFolder+'/dy.edf')
    spytIO.saveEdf(phi.real, saveFolder+'/phi.edf')
    spytIO.saveEdf(phi2.real, saveFolder+'/phi2.edf')
    spytIO.saveEdf(phi3.real, saveFolder+'/phi3.edf')


    Ir3D= spytIO.openSeq(referenceFiles)
    Is3D = spytIO.openSeq(sampleFiles)
    result = processProjectionSet(Is3D, Ir3D)

    dx = result['dx']
    dy = result['dy']
    phi = result['phi']
    phi2 = result['phi2']
    phi3 = result['phi3']
    spytIO.saveEdf(dx, saveFolder + '/dx.edf')
    spytIO.saveEdf(dy.real, saveFolder + '/dy.edf')
    spytIO.saveEdf(phi.real, saveFolder + '/phi.edf')
    spytIO.saveEdf(phi2.real, saveFolder + '/phi2.edf')
    spytIO.saveEdf(phi3.real, saveFolder + '/phi3.edf')

    #offsett=corr.registerRefAndSample(Ir,Is,1000)
    #print offsett





# open images



