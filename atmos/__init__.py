from _stdatm import stdatm2, stdatm1
from numpy import array, linspace, arange, trapz, exp
from miev0 import make_scattering_file
import pylab as plt

def trend(z, t0=288.0, p0=101290.0):
	if not hasattr(z,'__iter__'):
		z = array(z)

	_, _, rho = stdatm2(z, t0, p0)
	_, _, rho0 = stdatm1(z, t0, p0) 

	return rho/rho0

def taum_wl(wl):
	return 0.008569/(wl**4)*(1.0+0.0113/(wl**2)+0.00013/(wl**4)) 


def make_alt(layfname='atmos.lay', r0=0.1, r1=1.0, npts=101, 
				gamma=-3.5, midx = 1.4-0.0j, nmoms=30, 
				nlays=10, wl=0.750, hpbl=3.0, taua=0.1):
	Hmol = 7.3354 # столб атмосферы до 100 км
	extm = taum_wl(wl) / Hmol
	alts = arange(0.0, nlays)
	ext=trend(alts*1000.0)*extm
	taum0=trapz(ext, alts)
	
	exta = taua/hpbl
	altitude = alts[::-1]
	extinction_aer = exp(-altitude/hpbl)*exta
	extinction_mol = ext[::-1]
	with open(layfname, 'wt') as fout:
		for i in range(extinction_mol.shape[0]-1):
			exta_i = extinction_aer[i]
			extm_i = extinction_mol[i]
			scat_name = f"scat_file{i}"
			make_scattering_file(fname=scat_name, midx=midx, r0=r0, r1=r1,
                         gamma=gamma, npts=npts,
                         wl=wl, taua=exta_i, taum=extm_i,
                         nmoms=nmoms)
#			print(f"{altitude[i]:7.2f}{0.0:7.2f}{0.0:7.3f}\t'{scat_name}'\n", end='')
			fout.write(f"{altitude[i]:7.2f}{0.0:7.2f}{0.0:7.3f}\t'{scat_name}'\n")
		i=extinction_mol.shape[0]-1
		fout.write(f"{altitude[i]:7.2f}{0.0:7.2f}{0.0:7.3f}\t'           '\n")
	print(f"Extinction = {extinction_aer}")
	print(f"taum={-trapz(extinction_mol, altitude):7.3f} of {taum_wl(wl):7.3f}, taua={-trapz(extinction_aer, altitude):7.2f}")

	plt.figure()
	plt.plot(altitude, extinction_mol, label='mol')
	plt.plot(altitude, extinction_aer, label='aer')
	plt.legend()
	plt.grid(True)
	plt.savefig('density.png')


