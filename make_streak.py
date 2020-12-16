import astropy.io.fits as pyfits
import numpy as np
import math
import os
# To fix/add:
# width and height (x and y)?
# convolve lines with psf
# add point sources

class streakImage:
    def __init__(self, width, height):
        self.im_w = width
        self.im_h = height
        self.data = (np.zeros((height, width))).astype('int')

    def AddNoise(self, mean = 0., sigma = 1.):
            for i in range(self.im_h):
                for j in range(self.im_w):
                    self.data[i][j] += int(np.random.normal(mean, sigma))

    def GaussPSF(self, x, sigma): # to implement later
        return 1 / np.sqrt(2 * np.pi) / sigma * np.exp(- x * x / sigma / sigma / 2)
        # return 0

    def PointSrc(): # to implement later
        return 0

    def AddStreak(self, a, b, snr, width = 2):
        if np.absolute(a) >= 1:
            for i in range(self.im_w):
                for j in range(math.ceil(np.absolute(a))): # to avoid line breaks when a > 1
                    if(int(a * i + b) + j > 0 and int(a * i + b) + j < self.im_h):
                        self.data[int(a * i + b) + j][i] += 1. * snr
                        # b_n = int(a * i + b) + j - int(a_n * i)
                for k in range(math.ceil(np.absolute(a)), 30):
                    if(int(a * i + b) + k > 0 and int(a * i + b) + k < self.im_h):
                        self.data[int(a * i + b) + k][i] += 1. * snr * self.GaussPSF(k * np.cos(np.arctan(a)), width) / self.GaussPSF(0, width)

                for k in range(-30, 0):
                    if(int(a * i + b) + k  > 0 and int(a * i + b) + k < self.im_h):
                        self.data[int(a * i + b) + k][i] += 1. * snr * self.GaussPSF(k * np.cos(np.arctan(a)), width) / self.GaussPSF(0, width)

        if np.absolute(a) < 1:
            for i in range(self.im_w):
                # for j in range(math.ceil(np.absolute(1/a))): # to avoid line breaks when 1/a > 1
                if(int(a * i + b) > 0 and int(a * i + b) < self.im_h):
                    self.data[int(a * i + b)][i] += 1. * snr
                        # b_n = int(a * i + b) + j - int(a_n * i)

                for k in range(1, 30):
                    if(int(a * i + b) + k > 0 and int(a * i + b) + k < self.im_h):
                        self.data[int(a * i + b) + k][i] += 1. * snr * self.GaussPSF(k * np.cos(np.arctan(a)), width) / self.GaussPSF(0, width)

                for k in range(-30, 0):
                    if(int(a * i + b) + k  > 0 and int(a * i + b) + k < self.im_h):
                        self.data[int(a * i + b) + k][i] += 1. * snr * self.GaussPSF(k * np.cos(np.arctan(a)), width) / self.GaussPSF(0, width)

    def savefile(self, fname):
        if(os.path.isfile(fname)):
            os.remove(fname)
        hdu = pyfits.PrimaryHDU(self.data)
        hdu.writeto(fname)
        f1 = pyfits.open(fname, mode='update')
        hdr = f1[0].header
        hdr.set('NAXIS', 2)
        hdr.set('NAXIS1', self.im_w)
        hdr.set('NAXIS2', self.im_h)
        hdr.set('BINX', 1)
        hdr.set('BINY', 1)
        f1.flush()
        f1.close()
        f1 = pyfits.open(fname)
        hdr = f1[0].header
        print(repr(hdr))
        f1.close()

image_w = 300
image_h = 300
f_name = 'test.fits'
snr = 20
# streak: y = ax + b
a = 1
b = 0
width = 2
np.random.seed(0)

st = streakImage(image_w, image_h)
st.AddNoise()
st.AddStreak(a, b, snr, width)
st.savefile(f_name)
