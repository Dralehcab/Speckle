import spytIO
import glob
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

    # fourier filters
    ftfiltX = (1j * Qx / ((Qx**2 + Qy**2 + alpha))*beta)
    ftfiltX[np.isnan(ftfiltX)] = 0

    ftfiltY = (1j* Qy/ ((Qx**2 + Qy**2 + alpha))*beta)
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

    polarAngle = np.arctan(Qx, Qy)
    ftphi = fftshift(fft2(dX + i * dY))*np.exp(i*polarAngle)
    ftphi[np.isnan(ftphi)] = 0
    phi3 = ifft2(fftshift(ftphi))
    return phi3



def kottler2(dX,dY):
    print('kottler')
    NN, MM = dX.shape
    wx, wy = np.meshgrid(fftfreq(MM) * 2 * np.pi,
                         fftfreq(NN) * 2 * np.pi, indexing='xy')


    numerator = fft2(dX + 1j * dY)
    denominator=2.*np.pi*1j*wx - 2.*np.pi+1j*wy
    res = ifft2(numerator/denominator)
    #res -= np.mean(np.real(res))
    return res.real




def LarkinAnissonSheppard2(dX,dY):
    NN, MM = dX.shape
    wx, wy = np.meshgrid(fftfreq(MM) * 2 * np.pi,
                         fftfreq(NN) * 2 * np.pi, indexing='xy')

    # by using fftfreq there is no need to use fftshift
    numerator = fft2(dX + 1j*dY)
    denominator = 1j*wx - wy + np.finfo(dX.dtype).eps

    res = ifft2(numerator / denominator)
    res -= np.mean(np.real(res))
    return res.real


def LarkinAnissonSheppard(dx,dy,alpha =0 ,sigma=0):
    Nx, Ny = dx.shape
    i = complex(0, 1)
    G= dx + i*dy
    # fourier transfomm of the G function
    fourrierOfG = fftshift(fft2(G))
    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx)
    print(np.finfo(Qx.dtype).eps)
    ftfilt = 1 / (i*Qx - Qy + np.finfo(Qx.dtype).eps)
    ftfilt[np.isnan(ftfilt)] = 0
    phi=ifft2(ifftshift(ftfilt*fourrierOfG))
    phi=np.absolute(phi.real)
    return phi




def processOneProjection(Is,Ir):
    sigma = 0.95
    alpha = 0

    dI = (Is - Ir * (np.mean(gaussian_filter(Is,sigma=2.)) / np.mean(gaussian_filter(Ir,sigma=2.))))
    alpha=np.finfo(np.float32).eps
    dx, dy = derivativesByOpticalflow(Is, dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler2(dx, dy)
    phi2 = LarkinAnissonSheppard2(dx, dy)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3}

def processProjectionSetWithDarkFields(Is,Ir,dark):
    sigma = 1
    alpha = 0

    #Is=corrections.normalizationMultipleDarkField(Is,dark)
    #Ir=corrections.normalizationMultipleDarkField(Ir,dark)

    subImage=Is-Ir
    subImage=np.median(subImage,axis=0)

    dI = (subImage * (np.mean(Is) / np.mean(Ir)))
    dx, dy = derivativesByOpticalflow(np.mean(Ir,axis=0), dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler2(dx, dy)
    phi2 = LarkinAnissonSheppard2(dx, dy)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3}




def processProjectionSet(Is,Ir):
    sigma = 1
    alpha = 0

    subImage=Is-Ir
    subImage=np.median(subImage,axis=0)

    dI = (subImage * (np.mean(Is) / np.mean(Ir)))
    dx, dy = derivativesByOpticalflow(np.mean(Ir,axis=0), dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler2(dx, dy)
    phi2 = LarkinAnissonSheppard2(dx, dy)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3}


if __name__ == "__main__":
    Ir = spytIO.openImage('/Volumes/ID17/broncho/IHR_April2018/CigaleNuit/HA1000_Cigale_3um_gap90_75_Speckle__001_/ref0031_0000.edf')
    Is = spytIO.openImage('/Volumes/ID17/broncho/IHR_April2018/CigaleNuit/HA1000_Cigale_3um_gap90_75_Speckle__001_/HA1000_Cigale_3um_gap90_75_Speckle__001_0001.edf')

    result = processOneProjection(Is, Ir)
    dx = result['dx']
    dy = result['dy']
    phi = result['phi']
    phi2 = result['phi2']
    phi3 = result['phi3']
    spytIO.saveEdf(dx, '/Volumes/TMP_14_DAYS/embrun/dx.edf')
    spytIO.saveEdf(dy.real, '/Volumes/TMP_14_DAYS/embrun/dy.edf')
    spytIO.saveEdf(phi.real, '/Volumes/TMP_14_DAYS/embrun/phi.edf')
    spytIO.saveEdf(phi2.real, '/Volumes/TMP_14_DAYS/embrun/phiLarkinson.edf')
    spytIO.saveEdf(phi3.real, '/Volumes/TMP_14_DAYS/embrun/phiKottler.edf')
