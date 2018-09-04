from ctypes.util import find_library
from ctypes import CDLL, c_double, c_int, c_char_p

__all__=['do_calc']

try:
    
    
    #RTRANS_name=find_library("radtran")
    #print(RTRANS_name)
    RTRANS_lib = CDLL("lib/libradtran.dylib")

    do_calc1=RTRANS_lib.do_calc1
    do_calc1.argtypes=[c_double, c_double, c_int, c_double, 
                        c_double, c_double, c_double, c_double, 
                        c_double, c_double, c_int, c_double, 
                        c_double, c_int, c_int, c_int, c_int, c_char_p, c_char_p]
    
    

    def do_calc(r0=0.1, r1=1.0, npts=101, gamma=-3.5, wl=0.750,
                midx=1.6-0.01j, dens=300, hpbl=3.0, tau=0.2, numazim=2,
                galbedo=0.0, sza=40.0, nmu=32, nlays=10, aziorder=6,
                icalc=0,
                layer_file='atmoslay.lay',
                out_file='rt3.out'):
        """
        Вызов функции из DLL для расчета радиационного переноса
        Библиотека включает функцию `do_calc1`
        ей передаются следующие параметры
        ---------------------------------
    
        r0, r1, npts    -   начальный, конечный радиусы [um] и количество
                            интервалов разбиения 
        gamma           -   показатель степени степенного распределения
        wl              -   длина волны падающего излучения [um]
        midx            -   показатель преломления вещества частицы (complex)
        dens            -   концентрация частиц в 1 см3
        hpbl            -   высота погранслоя [km]
        tau             -   !аэрозольная! оптическая толща всей атмосферы
        numazim         -   число азимутальных углов
        aziorder        -   количество коэффициентов разложения в ряд по фурье
        icalc           -   если не 0, вычисляем радиационный перенос,
                            если нет - просто создаем необходимые файлы
        layer_file      -   имя файла со слоями
        out_file        -   имя файла
        """
                
        layerf=bytes(layer_file, 'utf-8')
        outf=bytes(out_file, 'utf-8')
        do_calc1(c_double(r0),
                c_double(r1),
                c_int(npts),
                c_double(wl),
                c_double(midx.real),
                c_double(midx.imag),
                c_double(gamma),
                c_double(dens),
                c_double(hpbl),
                c_double(tau),
                c_int(numazim),
                c_double(galbedo),
                c_double(sza),
                c_int(nmu),
                c_int(nlays),
                c_int(aziorder),
                c_int(icalc),
                c_char_p(layerf),
                c_char_p(outf))
except:
    print ("Не могу найти библиотеку и/или нужную функцию/")
