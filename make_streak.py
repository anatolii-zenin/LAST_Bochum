import astropy.io.fits as pyfits
import numpy as np
# To fix/add:
# width and height (x and y)?
# convolve lines with psf
# add point sources

class streakImage:
    def __init__(self, width, height):
        self.im_w = width
        self.im_h = height
        self.noise = noise
        self.data = (np.zeros((height, width))).astype('int')

    def AddNoise(self, mean = 10., sigma = 2.):
            for i in range(self.im_h):
                for j in range(self.im_w):
                    self.data[i][j] += int(np.random.normal(mean, sigma))

    def GaussPSF(): # to implement later
        return 0

    def PointSrc(): # to implement later
        return 0

    def AddStreak(self, a, b, scale):
        for i in range(self.im_w):
            for j in range(int(a)): # to avoid line breaks when a > 1
                if(int(a * i + b) + j > 0 and int(a * i + b) + j < self.im_h):
                    self.data[int(a * i + b) + j][i] += 1. * scale

    def savefile(self, fname):
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
noise = 1
line_width = 5
line_scale = 10.
# streak: y = ax + b
a = 3
b = 10
np.random.seed(0)

st = streakImage(image_w, image_h)
if(noise > 0):
    st.AddNoise()
st.AddStreak(a, b, line_scale)
st.savefile(f_name)
