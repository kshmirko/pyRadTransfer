#!/usr/bin/env python
# -*- encoding: utf-8 -*-

#"""
#Проверка правильности работы радиационых кодов
#"""

import matplotlib
#matplotlib.use("PDF")
import matplotlib.pyplot as plt
from ssrt.singlescat import run_model, run_model1
from rt3 import RT3
from atmos import make_alt, prepare_work_files
import numpy as np


def test_1():
    """
    Test polradtran
    """
    taua = 0.1
    taum = None
    r0   = 0.1
    r1   = 1.0
    npts = 101
    gamma= -3.5
    layfile = 'atmoslay.lay'
    midx = 1.5-0.02j
    nmoms= 40
    nlays= 15
    wl   = 0.870
    hpbl = 3.0

    # prepare_work_files(r0=r0,
    #     r1=r1,
    #     npts=npts,
    #     gamma=gamma,
    #     layfname=layfile,
    #     midx = midx, 
    #     nmoms=nmoms, 
    #     nlays=nlays, 
    #     wl=wl,  
    #     taua=taua,
    #     taum=taum,
    #     H=[2, hpbl+1])


                
    make_alt(r0=r0,
        r1=r1,
        npts=npts,
        gamma=gamma,
        layfname=layfile,
        midx = midx, 
        nmoms=nmoms, 
        nlays=nlays, 
        wl=wl, 
        hpbl=hpbl, 
        taua=taua,
        taum=taum)
    

    nmu = 32
    outfile = 'rt3.out'
    inttype = 'G'
    deltam  = 'N'
    
    sza     = 65.0
    direct_flux = np.pi*np.cos(np.deg2rad(sza))
    galbedo = 0.0
    wavelen = wl
    numazi  = 2
    theta_a, I_a, Q_a = RT3(nmu, layfile, outfile, inttype, deltam, 
        direct_flux, sza, galbedo, wavelen, numazi)

    #th_a1, I_a1 = run_model(sza, layfile)
    th_a1, I_a1 = run_model1(sza, layfile, taua, 0.89, 0.015)

    taua=0.3
    make_alt(r0=r0,
        r1=r1,
        npts=npts,
        gamma=gamma,
        layfname=layfile,
        midx = midx, 
        nmoms=nmoms, 
        nlays=nlays, 
        wl=wl, 
        hpbl=hpbl, 
        taua=taua,
        taum=taum)
    

    nmu = 32
    outfile = 'rt3.out'
    inttype = 'G'
    deltam  = 'N'
    
    sza     = 65.0
    direct_flux = np.pi*np.cos(np.deg2rad(sza))
    galbedo = 0.0
    wavelen = wl
    numazi  = 2
    theta_a, I_a2, Q_a = RT3(nmu, layfile, outfile, inttype, deltam, 
        direct_flux, sza, galbedo, wavelen, numazi)

    th_a3, I_a3 = run_model1(sza, layfile, taua, 0.89, 0.015)
    #th_a3, I_a3 = run_model(sza, layfile)
    # taum = None
    # taua = 0.0
    # make_alt(r0=r0,
    #     r1=r1,
    #     npts=npts,
    #     gamma=gamma,
    #     layfname=layfile,
    #     midx = midx, 
    #     nmoms=nmoms, 
    #     nlays=nlays, 
    #     wl=wl, 
    #     hpbl=hpbl, 
    #     taua=taua,
    #     taum=taum)

    # theta_m, I_m, Q_m = RT3(nmu, layfile, outfile, inttype, deltam,
    #     direct_flux, sza, galbedo, wavelen, numazi)
    # th_m, I_m1 = run_model1(sza, layfile, 0.0, 0.86, 0.015)

    # taum = None
    # taua = 0.1
    # make_alt(r0=r0,
    #     r1=r1,
    #     npts=npts,
    #     gamma=gamma,
    #     layfname=layfile,
    #     midx = midx, 
    #     nmoms=nmoms, 
    #     nlays=nlays, 
    #     wl=wl, 
    #     hpbl=hpbl, 
    #     taua=taua,
    #     taum=taum)
    
    # theta_am, I_am, Q_am = RT3(nmu, layfile, outfile, inttype, deltam, 
    #     direct_flux, sza, galbedo, wavelen, numazi)
    # th_am, I_am1 = run_model1(sza, layfile, 0.1, 0.86, 0.015)

    # galbedo = 0.3
    # theta_amg, I_amg, Q_amg = RT3(nmu, layfile, outfile, inttype, deltam, 
    #     direct_flux, sza, galbedo, wavelen, numazi)

    plt.figure()
    plt.semilogy(theta_a, I_a, label='$L_{aer}(0)$')
    plt.semilogy(th_a1, I_a1, label='$L_{aer1}(0)$')
    plt.semilogy(theta_a, I_a2, label='$L_{aer}(0)$')
    plt.semilogy(th_a3, I_a3, label='$L_{aer1}(0)$')
    #plt.semilogy(theta_m, I_m, label='$L_{mol}(0)$')
    #plt.semilogy(th_m, I_m1, label='$L_{mol1}(0)$')
    #plt.semilogy(theta_am, I_am, label='$L_{tot}(0)$')
    #plt.semilogy(th_am, I_am1, label='$L_{tot1}(0)$')
    #plt.semilogy(theta_amg, I_amg, label='$L_{tot}'+f"({galbedo})"+'$')
    plt.legend()
    plt.grid(True)
    plt.figure()
    print(f"theta_a = {theta_a}")
    print(f"th_a = {th_a1}")
    plt.plot(theta_a, I_a3[::-1]/I_a2)
    plt.grid()
    
    
def test_2():
    """
    Test ssrt
    """
    taua = 0.1
    taum = 0.0
    r0   = 0.1
    r1   = 1.0
    npts = 101
    gamma= -3.5
    layfile = 'atmoslay.lay'
    midx = 1.5-0.02j
    nmoms= 40
    nlays= 15
    wl   = 0.870
    hpbl = 3.0

    make_alt(r0=r0,
        r1=r1,
        npts=npts,
        gamma=gamma,
        layfname=layfile,
        midx = midx, 
        nmoms=nmoms, 
        nlays=nlays, 
        wl=wl, 
        hpbl=hpbl, 
        taua=taua,
        taum=taum)
    

    
    sza     = 35.0
    theta_a, I_a = run_model1(sza, layfile, taua, 0.86, 0.0)

    theta_m, I_m = run_model1(sza, layfile, 0.0, 0.86, 0.015)

    theta_am, I_am = run_model1(sza, layfile, 0.1, 0.86, 0.015)

    plt.figure()
    plt.semilogy(theta_a, I_a, label='La(0)')
    plt.semilogy(theta_m, I_m, label='Lm(0)')
    plt.semilogy(theta_am, I_am, label='L(0)')
    plt.legend()
    plt.grid(True)
    


def test_3():
    taua = 0.05
    taum = None
    r0   = 0.1
    r1   = 1.0
    npts = 101
    gamma= -3.5
    layfile = 'atmoslay.lay'
    midx = 1.5-0.02j
    nmoms= 40
    nlays= 15
    wl   = 0.870
    hpbl = 2.0

    prepare_work_files(r0=r0,
        r1=r1,
        npts=npts,
        gamma=gamma,
        layfname=layfile,
        midx = midx, 
        nmoms=nmoms, 
        nlays=nlays, 
        wl=wl,  
        taua=taua,
        taum=taum,
        H=[5, hpbl+5])

    nmu = 32
    outfile = 'rt3.out'
    inttype = 'G'
    deltam  = 'N'
    direct_flux = 1.0
    sza     = 65.0
    galbedo = 0.0
    wavelen = wl
    numazi  = 2
    
    theta_a1, I_a1, Q_a1 = RT3(nmu, layfile, outfile, inttype, deltam, 
        direct_flux, sza, galbedo, wavelen, numazi)

    prepare_work_files(r0=r0,
        r1=r1,
        npts=npts,
        gamma=gamma,
        layfname=layfile,
        midx = midx, 
        nmoms=nmoms, 
        nlays=nlays, 
        wl=wl,  
        taua=taua,
        taum=taum,
        H=[1, hpbl+1])
    

    theta_a2, I_a2, Q_a2 = RT3(nmu, layfile, outfile, inttype, deltam, 
        direct_flux, sza, galbedo, wavelen, numazi)

    prepare_work_files(r0=r0,
        r1=r1,
        npts=npts,
        gamma=gamma,
        layfname=layfile,
        midx = midx, 
        nmoms=nmoms, 
        nlays=nlays, 
        wl=wl,  
        taua=taua,
        taum=taum,
        H=[2, hpbl+2])
    

    theta_a3, I_a3, Q_a3 = RT3(nmu, layfile, outfile, inttype, deltam, 
        direct_flux, sza, galbedo, wavelen, numazi)


    prepare_work_files(r0=r0,
        r1=r1,
        npts=npts,
        gamma=gamma,
        layfname=layfile,
        midx = midx, 
        nmoms=nmoms, 
        nlays=nlays, 
        wl=wl,  
        taua=taua,
        taum=taum,
        H=[3, hpbl+3])
    

    
    theta_a4, I_a4, Q_a4 = RT3(nmu, layfile, outfile, inttype, deltam, 
        direct_flux, sza, galbedo, wavelen, numazi)



    prepare_work_files(r0=r0,
        r1=r1,
        npts=npts,
        gamma=gamma,
        layfname=layfile,
        midx = midx, 
        nmoms=nmoms, 
        nlays=nlays, 
        wl=wl,  
        taua=taua,
        taum=taum,
        H=[4, hpbl+4])
    

    
    theta_a5, I_a5, Q_a5 = RT3(nmu, layfile, outfile, inttype, deltam, 
        direct_flux, sza, galbedo, wavelen, numazi)


    plt.figure()
    plt.title(f"sza={sza:4.1f}  wl={wavelen:5.3f}\n"+
                f"taua={taua}")
    plt.semilogy(theta_a1, I_a1, label='$L_{aer}(0), hpbl[5,'+f"{hpbl+5}]"+'$')
    plt.semilogy(theta_a2, I_a2, label='$L_{aer}(0), hpbl[1,'+f"{hpbl+1}]"+'$')
    plt.semilogy(theta_a3, I_a3, label='$L_{aer}(0), hpbl[2,'+f"{hpbl+2}]"+'$')
    plt.semilogy(theta_a4, I_a4, label='$L_{aer}(0), hpbl[3,'+f"{hpbl+3}]"+'$')
    plt.semilogy(theta_a5, I_a5, label='$L_{aer}(0), hpbl[4,'+f"{hpbl+4}]"+'$')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"test_3_{taua:5.3f}_I.pdf")
    
    plt.figure()
    plt.title(f"sza={sza:4.1f}  wl={wavelen:5.3f}\n"+
                f"taua={taua}")
    plt.plot(theta_a1, Q_a1, label='$L_{aer}(0), hpbl[5,'+f"{hpbl+5}]"+'$')
    plt.plot(theta_a2, Q_a2, label='$L_{aer}(0), hpbl[1,'+f"{hpbl+1}]"+'$')
    plt.plot(theta_a3, Q_a3, label='$L_{aer}(0), hpbl[2,'+f"{hpbl+2}]"+'$')
    plt.plot(theta_a4, Q_a4, label='$L_{aer}(0), hpbl[3,'+f"{hpbl+3}]"+'$')
    plt.plot(theta_a5, Q_a5, label='$L_{aer}(0), hpbl[4,'+f"{hpbl+4}]"+'$')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"test_3_{taua:5.3f}_Q.pdf")

    plt.figure()
    plt.title(f"sza={sza:4.1f}  wl={wavelen:5.3f}\n"+
                f"taua={taua}")
    plt.plot(theta_a1, -Q_a1/I_a1, label='$L_{aer}(0), hpbl[5,'+f"{hpbl+5}]"+'$')
    plt.plot(theta_a2, -Q_a2/I_a2, label='$L_{aer}(0), hpbl[1,'+f"{hpbl+1}]"+'$')
    plt.plot(theta_a3, -Q_a3/I_a3, label='$L_{aer}(0), hpbl[2,'+f"{hpbl+2}]"+'$')
    plt.plot(theta_a4, -Q_a4/I_a4, label='$L_{aer}(0), hpbl[3,'+f"{hpbl+3}]"+'$')
    plt.plot(theta_a5, -Q_a5/I_a5, label='$L_{aer}(0), hpbl[4,'+f"{hpbl+4}]"+'$')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"test_3_{taua:5.3f}_Pol.pdf")

def main():
    """
    главная функция, вызывает различные тестовые расчеты
    """
    test_1()
    #test_3()
    #test_2()
    plt.show()


if __name__ == '__main__':
	main()