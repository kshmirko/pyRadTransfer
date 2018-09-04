from _rt3 import rt3
from numpy import cos, deg2rad, loadtxt, zeros, vstack, rad2deg, arccos


def RT3(nmu=32, layfile='layfile.lay', outfile='rt3.out', inttype='D', 
        deltam='N', direct_flux=1.0, sza=40, galbedo=0.0, wavelen=0.750,
        numazimuth=2):
    """
    Драйвер для вызова расчетной модели PolRadTran
        nmu         -   количество отсчетов косинуса зенитного угла
        layfile     -   имя файла с описанием структуры атмосферы
        outfile     -   имя выходного файла с результатами расчетов
        inttype     -   метод интегрирования
                        'G' -   Гаусса
                        'D' -   двойного Гаусса
                        'L' -   Лобатто
        deltam      -   включение/выключение режима масштабирования 
                        фазовой функции
                        'Y' -   включение
                        'N' -   выключение
        direct_flux -   поток на верхней границе атмосферы
        sza         -   зенитный угол солнца [градусы]
        galbedo     -   альбедо подстилающей поверхности
        wavelen     -   длина волны падающего излучения
        numazimuth  -   количество азимутальных углов для вывода
    """
    
    mu = zeros(2*nmu)
    src_code = 1 # solar radiation only
    rt3(2, 4, mu, src_code, layfile, outfile, inttype, deltam, direct_flux,
        cos(deg2rad(sza)), 0.0, 'L', galbedo, 0.0+0j, 0.0, wavelen,
         'W', 'IQ', numazimuth)
         
    r = loadtxt(outfile, skiprows=11)
    idx1=(r[:,0]==0)&(r[:,1]==0)&(r[:,2]>0)
    idx2=(r[:,0]==0)&(r[:,1]==180)&(r[:,2]>0)
    
    sl1 = r[idx1,:]
    sl2 = r[idx2,:]
    
    d = vstack((sl2, sl1[::-1,:]))
    phi=cos(deg2rad(d[:,1])) 
    umu = phi*rad2deg(arccos(d[:,2]))
    return umu, d[:,3], d[:,4]
    