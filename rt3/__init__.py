from _rt3 import run1, run2
from numpy import cos, deg2rad, loadtxt, zeros, vstack, rad2deg, arccos
import matplotlib.pyplot as plt

class RT3:

    def __init__(self, sza=40, nummu = 32, galbedo = 0.0):
        self.layf = 'atmos.lay'
        self.outf = 'rt3.out'
        self.midx = 1.5-0.0j
        self.r0 = 0.1
        self.r1 = 1.0
        self.gamma = -3.5
        self.npts = 101
        self.wl = 0.75
        self.taua = 0.1
        self.nmoms = 40
        self.direct_mu = cos(deg2rad(sza))
        self.quad_type = 'G'
        self.deltam = 'N'
        self.ground_albedo = galbedo
        self.nummu = nummu
        self.ssa_a = 1.0

    
    def __get_sza(self):
        return rad2deg(arccos(self.direct_mu))
    
    
    def __set_sza(self, sza):
        self.direct_mu = cos(deg2rad(sza))

    sza = property(__get_sza, __set_sza)
    
    def run(self):

        self.mu, self.Iv, self.Qv = run1(self.layf, self.outf, self.midx, 
            self.r0, self.r1, self.gamma, self.npts, self.wl, self.taua,
            self.nmoms, self.direct_mu, self.quad_type, self.deltam,
            self.ground_albedo, self.nummu)

    def run1(self):

        self.mu, self.Iv, self.Qv = run2(self.layf, self.outf, self.midx, 
            self.r0, self.r1, self.gamma, self.npts, self.wl, self.taua, self.ssa_a,
            self.nmoms, self.direct_mu, self.quad_type, self.deltam,
            self.ground_albedo, self.nummu)
                                                                                          
    def plotI(self):
        caption = "$"+f"I({self.sza})"+"$"
        plt.semilogy(self.mu, self.Iv, label=caption)
        

    def plotQ(self):
        caption = "$"+f"Q({self.sza})"+"$"
        plt.plot(self.mu, self.Iv)
        
    def show(self):
        plt.legend()
        plt.grid(True)
        plt.show()                                                           
                                
