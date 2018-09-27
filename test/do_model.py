#!/usr/bin/env python


import rt3
import numpy as np
import matplotlib.pyplot as plt


def main():
        SZA=0.1
        Rg = 0.3
        
        Rt = rt3.RT3()
        Rt.taum=0
        Rt.sza = SZA
        Rt.ground_albedo = Rg
        Rt.run1()
        mu, I, Q = Rt.mu, Rt.Iv, Rt.Qv
        Rt.Iv /= np.pi
        Rt.plotI()
        
        Rt.ground_albedo = 0	
        Rt.run1()
        mu, I0, Q0 = Rt.mu, Rt.Iv, Rt.Qv
        Rt.Iv /= np.pi
        Rt.plotI()
        

        plt.figure()
        plt.plot(mu, I-I0, label='$\Delta L$')
        plt.legend()
        Rt.show()
        

if __name__ == '__main__':
	main()
