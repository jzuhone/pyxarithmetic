import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline

from soxs import ApecGenerator, AuxiliaryResponseFile

class SpectrumMath(object):
    def __init__(self, band1, band2, abund, redshift, nbins=10000,
                 arf=None, apec_root=None, apec_vers="2.0.2"):
        emin = band1[0]
        emax = band2[1]
        self.band1 = band1
        self.band2 = band2
        self.abund = abund
        if arf is not None:
            self.arf = AuxiliaryResponseFile(arf)
        else:
            self.arf = arf
        agen = ApecGenerator(emin, emax, nbins, apec_root=apec_root,
                             apec_vers=apec_vers, broadening=False)
        self.redshift = redshift
        cspec, mspec = agen._get_table(list(range(agen.nT)), redshift, 0.0)
        range1 = np.logical_and(agen.emid >= band1[0], agen.emid <= band1[1])
        range2 = np.logical_and(agen.emid >= band2[0], agen.emid <= band2[1])
        spec1 = cspec[:,range1] + abund*mspec[:,range1]
        spec2 = cspec[:,range2] + abund*mspec[:,range2]
        if self.arf is not None:
            e1 = agen.emid[range1]
            e2 = agen.emid[range2]
            spec1 *= self.arf.interpolate_area(e1).value
            spec2 *= self.arf.interpolate_area(e2).value
        lambda1 = spec1.sum(axis=1)
        lambda2 = spec2.sum(axis=1)
        self.lT1 = InterpolatedUnivariateSpline(agen.Tvals, lambda1)
        self.lT2 = InterpolatedUnivariateSpline(agen.Tvals, lambda2)

    def get_emiss(self, kT):
        return self.lT1(kT), self.lT2(kT)

    def get_emiss_derivative(self, kT):
        dldt1 = self.lT1(kT, nu=1)
        dldt2 = self.lT2(kT, nu=1)
        l1 = self.lT1(kT)
        l2 = self.lT2(kT)
        return dldt1*(kT/l1), dldt2*(kT/l2)

    def get_flux_ratio(self, alpha, kT):
        if alpha == 0.0:
            return kT/kT
        else:
            dldt1, dldt2 = self.get_emiss_derivative(kT)
            return (2.0+alpha*dldt2)/(2.0+alpha*dldt1)

    def compute_coefficients(self, alpha_0, alpha_R, kT):
        R_0 = self.get_flux_ratio(alpha_0, kT)
        R_R = self.get_flux_ratio(alpha_R, kT)
        cB2 = 1.0/(R_R-R_0)
        cB1 = -cB2*R_0
        return cB1, cB2