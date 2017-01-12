import numpy as np
import astropy.io.fits as pyfits
from six import string_types
from scipy.interpolate import InterpolatedUnivariateSpline

class ClusterImage(object):
    def __init__(self, source_img, model_img=None):
        self.source_img = self._get_image(source_img) 
        if model_img is None:
            self.model_img = self._make_model(self.source_img)
        else:
            self.model_img = self._get_image(model_img)
        self.resid_img = self.source_img/self.model_img-1.0

    def _get_image(self, img):
        if not isinstance(img, np.ndarray):
            if isinstance(img, string_types):
                f = pyfits.open(img)
                img = f[0].data
                f.close()
            elif isinstance(img, pyfits.HDUList):
                img = img[0].data
            elif isinstance(img, pyfits.ImageHDU) or \
                (isinstance(img, pyfits.PrimaryHDU) and img.is_image):
                img = img.data
        return img

    def _make_model(self, img):
        nx, ny = img.shape
        x, y = np.mgrid[0:nx, 0:ny]
        r = np.sqrt((x - x.mean())**2 + (y - y.mean())**2)
        rmax = 0.5*np.sqrt(nx*nx+ny*ny)
        nbins = max(nx, ny)//2
        rbin = np.linspace(0.0, rmax, nbins+1)
        p, _ = np.histogram(r, bins=rbin, weights=img)
        n, _ = np.histogram(r, bins=rbin)
        p /= n
        rmid = 0.5 * (rbin[1:] + rbin[:-1])
        pr = InterpolatedUnivariateSpline(rmid, p)
        return pr(r)

def add_images(img1, img2, cB1, cB2):
    X = cB1*img1.resid_img+cB2*img2.resid_img
    return pyfits.ImageHDU(X)
