#!/usr/bin/env python


import rt3
import numpy as np
import matplotlib.pyplot as plt


def main():
	SZA=65
	Rg = 0.1

	Rt = rt3.RT3()
	Rt.sza = SZA
	Rt.ground_albedo = Rg

	Rt.run()
	mu, I, Q = Rt.mu, Rt.Iv, Rt.Qv

	Rt.plotI()
	
	Rt.ground_albedo = 0	
	Rt.run()
	mu, I0, Q0 = Rt.mu, Rt.Iv, Rt.Qv
	
	Rt.plotI()
	

	plt.figure()
	plt.plot(mu, I-I0, label='$\Delta L$')
	plt.legend()
	Rt.show()


if __name__ == '__main__':
	main()