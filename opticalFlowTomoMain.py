import spytIO
import glob
import os
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2
from numpy.fft import ifft2 as ifft2

from math import pi as pi
from math import floor as floor

import frankoChellappa as fc

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


def parseESRFTomoFolder(folderpath):
    print('ESRFTomoFolder')
    scanName=os.path.basename(folderpath)
    parametersScanFilename=folderpath+'/'+scanName+'.xml'
    print(parametersScanFilename)
    tomoExperiment=esrfTomo.FastTomoExperiment(parametersScanFilename)
    print('numberFlatField: ')
    print(tomoExperiment.numberFlatField)
    tomoExperiment.createAverageWfandDf()
    tomoExperiment.findCenterOfRotation()
    print('Cor Found at '+str(tomoExperiment.cor))
    projectionsFileNames=tomoExperiment.getProjectionsName()
    projectionsFileNames.sort()
    darkFieldFilename=tomoExperiment.darkOutputFile
    referenceFileNames= tomoExperiment.getReferencesFileNames()
    referenceFileNames.sort()
    print(referenceFileNames)
    return projectionsFileNames,referenceFileNames,darkFieldFilename



def processOneProjection(Is,Ir):
    sigma = 0.9
    alpha = 0

    dI = (Is - Ir * (np.mean(Is) / np.mean(Ir)))
    dx, dy = derivativesByOpticalflow(Ir, dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3}



def processTomoFolder(projectionsF,referencesF,darkFieldF,outputFolder):
    print('Process Tomo')
    Ir=spytIO.openImage(referencesF[0])
    darkField=spytIO.openImage(darkFieldF)
    Ir = corr.normalization2D(Ir, darkField)

    dxFolder=outputFolder+'/dx/'
    dyFolder = outputFolder + '/dy/'
    phiFolder = outputFolder + '/phi/'
    phi2Folder = outputFolder + '/phi2/'
    phi3Folder = outputFolder + '/phi3/'

    if not(os.path.exists(outputFolder)):
        os.mkdir(outputFolder)

    if not (os.path.exists(dxFolder)):
        os.mkdir(dxFolder)
    if not (os.path.exists(dyFolder)):
        os.mkdir(dyFolder)
    if not (os.path.exists(phiFolder)):
        os.mkdir(phiFolder)
    if not (os.path.exists(phi2Folder)):
        os.mkdir(phi2Folder)
    if not (os.path.exists(phi3Folder)):
        os.mkdir(phi3Folder)



    for proj in projectionsF:
        print(proj)
        Is=spytIO.openImage(proj)
        Is=corr.normalization2D(Is,darkField)

        #Ir=corr.registerRefAndSample(Ir,Is,1000)
        result=processOneProjection(Is,Ir)

        dx = result['dx']
        dx=np.asarray(dx.real,np.float32)
        dxFilename=dxFolder+'/dx_'+os.path.basename(proj)
        spytIO.saveEdf(dx, dxFilename)

        dy = result['dy']
        dy = np.asarray(dy.real, np.float32)
        dyFilename = dyFolder + '/dy_' + os.path.basename(proj)
        spytIO.saveEdf(dy, dyFilename)

        phi = result['phi']
        phi = np.asarray(phi.real, np.float32)
        phiFilename = phiFolder + '/phi_' + os.path.basename(proj)
        spytIO.saveEdf(phi, phiFilename)

        phi2 = result['phi2']
        phi2 = np.asarray(phi2.real, np.float32)
        phi2Filename = phi2Folder + '/phi2_' + os.path.basename(proj)
        spytIO.saveEdf(phi2, phi2Filename)

        phi3 = result['phi3']
        phi3 = np.asarray(phi3.real, np.float32)
        phi3Filename = phi3Folder + '/phi3' + os.path.basename(proj)
        spytIO.saveEdf(phi3, phi3Filename)
        print phi3Filename






if __name__ == "__main__":
    inputFolder='/Volumes/ID17/speckle/md1097/id17/SpeckleSacroIliaque/HA1100_sacroIlliaque_speckle_xss_6um_52kev__001__008_'
    outputFolder='/Volumes/ID17/speckle/md1097/id17/SpeckleSacroIliaque/OuputSacro_Floor0_'
    projectionsFileNames, referenceFileNames, darkFieldFilename= parseESRFTomoFolder(inputFolder)
    processTomoFolder(projectionsFileNames, referenceFileNames, darkFieldFilename,outputFolder)


    print(' Optical Flow Tomo ')
    print('Test One File')
    Ir=spytIO.openImage('ref1-1.edf')
    Is= spytIO.openImage('samp1-1.edf')


    result=processOneProjection(Is,Ir)
    dx = result['dx']
    dy = result['dy']
    phi = result['phi']
    phi2 = result['phi2']
    phi3 = result['phi3']
    spytIO.saveEdf(dx, 'output/dx.edf')
    spytIO.saveEdf(dy.real, 'output/dy.edf')
    spytIO.saveEdf(phi.real, 'output/phi.edf')
    spytIO.saveEdf(phi2.real, 'output/phi2.edf')
    spytIO.saveEdf(phi3.real, 'output/phi3.edf')




    #offsett=corr.registerRefAndSample(Ir,Is,1000)
    #print offsett





# open images



