#!/usr/bin/env python

import click
from ssrt.singlescat import run_model
from rt3 import RT3
from atmos import make_alt
import matplotlib.pyplot as plt
from numpy import savez_compressed, load, abs, argmin, pi

@click.group()
def main():
    pass
    
@main.command()
@click.option('--sza', default=30.0, type=float, 
                help='Зенитный угол солнца [градусы]')
@click.option('--layfile', default='layfile.lay',
                type=click.Path(exists=True), 
                help='Путь к файлу с описанием структуры атмосферы. '+
                    'Важно, чтоб все прилагающиеся scat файлы были там же.')
#@click.option('--plot', type=bool)
@click.option('--datafile',type=click.Path(exists=False),
    default='ssrt.npz',
    help="Имя файла с результатами")
def ssrt(sza, layfile, datafile):
    """
    Расчет радиационного переноса в приближении однократного рассеяния
    """
    theta, I = run_model(sza, layfile)
    #for i,_ in enumerate(theta):
    #    print(f"{theta[i]:7.2f}{I[i]:12.3e}")
        

    savez_compressed(datafile, sza=sza, layfile=layfile, theta=theta, I=I)
    #plt.figure()
    #plt.plot(theta, I)
    #plt.show()
     

        
@main.command()
@click.option('--nmu', type=int, default=32, 
                help='Количество отсчетов по зенитному углу')
@click.option('--layfile', default='layfile.lay',
                type=click.Path(exists=False), 
                help='Путь к файлу с описанием структуры атмосферы. '+
                    'Важно, чтоб все прилагающиеся scat файлы были там же.')
@click.option('--outfile', default='rt3.out',
                type=click.Path(exists=False), 
                help='Путь к файлу с результатами расчетов')
@click.option('--inttype', default='G',
                type=click.Choice(['G','D','L']), 
                help='Способ интегрирования')
@click.option('--deltam', type=click.Choice(['Y','N']), default='N',
                help='Delta-M scaling option')
@click.option('--direct_flux', type=float, default=1.0, 
                help='Поток на акрхней границе атмосферы')
@click.option('--sza', default=30.0, type=float, 
                help='Зенитный угол солнца [градусы]')    
@click.option('--galbedo', default=0.0, type=float, 
                help='Альбедо допстилающей поверхности') 
@click.option('--wavelen', default=0.750, type=float, 
                help='Длина волны падающего излучения [um]') 
@click.option('--numazi', type=int, default=32, 
                help='Количество азимутальных углов в выводе')
@click.option('--datafile',type=click.Path(exists=False),
    default='rt3.npz',
    help="Имя файла с результатами")
def rt3(nmu, layfile, outfile, inttype, deltam, direct_flux, sza, 
            galbedo, wavelen, numazi, datafile):
    theta, I, Q = RT3(nmu, layfile, outfile, inttype, deltam, direct_flux,
                     sza, galbedo, wavelen, numazi)
    
    savez_compressed(datafile,
        layfile=layfile,
        outfile=outfile,
        inttype=inttype,
        deltam=deltam,
        direct_flux=direct_flux,
        sza=sza,
        galbedo=galbedo,
        wavelen=wavelen,
        numazi=numazi,
        datafile=datafile,
        theta=theta,
        I=I,
        Q=Q)
    #plt.figure()
    #plt.plot(theta, I)
    #plt.show()

@main.command()
@click.option('--r0', type=float, default=0.1, 
                help='Левая граница размеров частиц')
@click.option('--r1', type=float, default=1.0, 
                help='Правая граница размеров частиц')
@click.option('--npts', type=int, default=101, 
                help='Количество отсчетов на интервал размеров частиц')
@click.option('--gamma', type=float, default=-3.5, 
                help='Показатель степени распределения')
@click.option('--wl', type=float, default=0.750, 
                help='Длина волны, [um]')
@click.option('--midx', type=complex, default=1.5-0.0j, 
                help='Показатель преломления частиц')
@click.option('--nmoms', type=int, default=40, 
                help='Число коэффициентов Лежандра')
@click.option('--nlays', type=int, default=10, 
                help='Количество слоев толщиной в 1 км')
@click.option('--hpbl', type=float, default=3.0, 
                help='Высота погранслоя')
@click.option('--taua', type=float, default=0.1, 
                help='Аэрозольная оптическая толща')
@click.option('--taum', type=float, default=None, 
                help='Молекулярная оптическая толща')
@click.option('--layfile', type=click.Path(exists=False),  
                help='Имя файла с описанием атмосферы', 
                required=True)
@click.option('--datafile',type=click.Path(exists=False),
    default='prepare.npz',
    help="Имя файла с параметрами")
def prepare_files(r0, r1, npts, gamma, wl, midx, nmoms, nlays, hpbl, taua, taum, layfile, datafile):
    """
    Подготавливает файлы с описанием атмосферы
    
    """
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

    savez_compressed(datafile,
        taua=taua,
        r0=r0,
        r1=r1,
        gamma=gamma,)

@main.command()
@click.option('--ssrt',type=click.Path(exists=True),
    default='ssrt.npz',
    help="Имя файла с результатами")
@click.option('--rt3',type=click.Path(exists=True),
    default='rt3.npz',
    help="Имя файла с результатами")
@click.option('--prepare',type=click.Path(exists=True),
    default='prepare.npz',
    help="Имя файла с результатами")
@click.option('--savef',type=click.Path(exists=False),
    default=None,
    help="Имя файла с графиком")
def vizualize(ssrt, rt3, prepare, savef):
    F1 = load(ssrt)
    F2 = load(rt3)
    F3 = load(prepare)
    plt.figure()
    #print(F1['theta'][:10], F2['theta'][:10])
    #print(F1['I'][:10], F2['I'][:10])
    
    #plt.semilogy(F1['theta'], F1['I'], label='ssrt')
    #plt.semilogy(F2['theta'][::-1], F2['I'][::-1]/2, label='rt3')
    plt.plot(F2['theta'][::-1], 2*F1['I']/F2['I'][::-1], label='ss/rt3')

    plt.legend()
    plt.grid(True)
    plt.title(f"sza  ={F1['sza']:4.1f} wl   ={F2['wavelen']:5.3f} tau_a={F3['taua']:5.3f}\n"+
              f"r0   ={F3['r0']:4.2f} r1   ={F3['r1']:4.2f} ɣ   ={F3['gamma']:4.2f}")
    plt.xlabel('Θ, deg.')
    plt.ylabel('Radiance, W/m^2/sr/um')
    if not savef is None:
        #save_fig = f"{savef}_g{F2['galbedo']:4.2f}_t{F2['taua']:5.3f}_ratio.png"
        save_fig=f"{savef}_L.png"
        plt.savefig(save_fig)
    
    plt.figure()
    idx1 = argmin(abs(F1['theta']-F1['sza']))
    idx2 = argmin(abs(F2['theta']-F2['sza']))
    print(idx1, idx2, F1['I'][idx1], F2['theta'][idx2] )
    plt.semilogy(F1['theta'], F1['I'], label='ssrt')
    plt.semilogy(F2['theta'], F2['I']/2, label='rt3')
    plt.legend()
    plt.grid(True)
    plt.title(f"sza  ={F1['sza']:4.1f} wl   ={F2['wavelen']:5.3f} tau_a={F3['taua']:5.3f}\n"+
              f"r0   ={F3['r0']:4.2f} r1   ={F3['r1']:4.2f} ɣ   ={F3['gamma']:4.2f}")
    plt.xlabel('Θ, deg.')
    plt.ylabel('Radiance, W/m^2/sr/um')
    if  savef is None:
        plt.show()
    else:
        save_fig = f"{savef}_LN.png"
        plt.savefig(save_fig)
        pass
    F1.close()
    F2.close()
    F3.close()


if __name__ == '__main__':
    main()