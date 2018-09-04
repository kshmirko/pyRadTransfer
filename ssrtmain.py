#!/usr/bin/env python

import click
from ssrt.singlescat import run_model
from rt3 import RT3
import matplotlib.pyplot as plt

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
@click.option('--plot', type=bool)
def ssrt(sza, layfile, plot):
    """
    Расчет радиационного переноса в приближении однократного рассеяния
    """
    theta, I = run_model(sza, layfile)
    for i,_ in enumerate(theta):
        print(f"{theta[i]:7.2f}{I[i]:12.3e}")
        
    plt.figure()
    plt.plot(theta, I)
    plt.show()
     

        
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
def rt3(nmu, layfile, outfile, inttype, deltam, direct_flux, sza, 
            galbedo, wavelen, numazi):
    theta, I, Q = RT3(nmu, layfile, outfile, inttype, deltam, direct_flux,
                     sza, galbedo, wavelen, numazi)
    
    plt.figure()
    plt.plot(theta, I)
    plt.show()


if __name__ == '__main__':
    main()