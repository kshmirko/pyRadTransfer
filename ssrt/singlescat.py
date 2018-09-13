import numpy as np
from scipy.integrate import  cumtrapz
from numpy.polynomial.legendre import Legendre
from . import sinsca, getmom


def get_ss_downward_intensity(phi, theta, theta0, cumtau,\
        pmoms, ssas):
    """
    Вычитяет интенсивность нисходящего излучения в приближении однократного
    рассеяния. Атмосфера задается параметрами `cumtau`, `pmoms` и `ssas`.
        
    cumtau(nlays)       -   вектор кумулятивной оптической толщи (сверху вниз)
    pmoms(nlays,nmom)   -   коэффициенты полиномов лежандра для 
                            фазовых функций для каждого слоя
    ssas(nlays)         -   альбедо однократного рассеяния для каждого слоя
    phi(npts)           -   азимутальный угол наблюдения
    theta(npts)         -   зенитный угол наблюдения (град)
    theta0              -   зенитный углы солнца (град)
        
    Возвращает:
        I(phi, theta)   -   интенсивность нисходящего излучения в приближении 
                            однократнго рассеяния
    """
    
    uphi = np.cos(np.deg2rad(phi))
    umu  = np.cos(np.deg2rad(theta))
    umu0 = np.cos(np.deg2rad(theta0))
    
    I = np.zeros_like(umu)
    
    fbeam = np.pi/umu0
   
    #print(len(pmoms))
    PF = [Legendre(pmoms[i]) for i in range(len(pmoms))]
    
    for i, e in enumerate(umu):
        th0 = e*umu0+\
                    np.sqrt(1.0-e*e)*np.sqrt(1.0-umu0*umu0)*uphi[i]
        
        phase = np.array([PF[j](th0) for j in range(cumtau.shape[0])])/2.0 # in polradtran phase function normalized to 2.0
        #print(phase)
        I[i] = sinsca(0.00000001, 1, phase, ssas, cumtau, -e, umu0,\
            cumtau[-1], fbeam, np.pi)
    return I/2
    

def load_atmos(layfile):
    """
    Загружаем свойства слоев из файлов нв диске.
    `layfile`   -   файл с описанием структуры атмосферы (идентичный для 
                    использования с моделью polradtran rt3)
    
    Срукткра файла следующая:
    Он состоит из N+1 зписей, определяющих N атмосферных слоев
    
    H       Temp        GasExtinct      ScatFile
    
    H(1)    -   высота верхней границы самого верхнего слоя
    H(N+1)  -   высота нижней границы самого нижнего слоя
    Temp(I), 
    GasExtinct(I), 
    ScatFile(I)     -   температура, газовая экстинкция и имя файла рассеяния
                        для I-го слоя
    
    Формат `layfile`
    ~~~~~~~~~~~~~~~~
    
    30      0       0       '.sca_001'
    20      0       0       '.sca_002'
    0       0       0       '        '
    
    Вышеприведенная запист поределяет атмосферу из двух слоев:
    1.  0-20 km
    2.  20-30 km
    
    Формат файлов рассеяния
    ~~~~~~~~~~~~~~~~~~~~~~~
    
    EXT
    SCA
    ALB
    NMOM
    0   P1  P2  P3  P4  P5  P5
    1   ..  ..  ..  ..  ..  ..
    2   ..  ..  ..  ..  ..  ..
    NMOM..  ..  ..  ..  ..  ..
    
    P1-P6 - элементы матрицы Эванса
                                            
    """
    H=[]
    Temp=[]
    RayExt=[]
    ScaFiles=[]
    LaySca=[]
    LayExt=[]
    LaySSA=[]
    LayP = []
    with open(layfile) as fin:
        while True:
            line = fin.readline()
            if len(line)<2: break # Last line in file
            parts = line.split()
            scafile = parts[-1].replace('\'','') 
            
            H.append(float(parts[0]))
            Temp.append(float(parts[1]))
            RayExt.append(float(parts[2]))
            ScaFiles.append(scafile)
            
            if scafile is '' :
                ext=0
                sca=0
                alb=0
                P1 = 0.0
            else:
                ext, sca, alb, P1 = load_scafile(scafile)
            LayExt.append(ext)
            LaySca.append(sca)
            LaySSA.append(alb)
            LayP.append(P1)
    
    # вычисляем кумулятивную оптическую толщу
    #cumTau=np.r_[0.0, -cumtrapz(LayExt, H)]
    cumTau=np.r_[np.cumsum(LayExt)]
    return cumTau, LayP, np.array(LaySSA)
        
    
def load_scafile(scafile):
    """
    Загрузка данных из scafile
    
    Формат .sca_001
    ~~~~~~~~~~~~~~~
    
    EXT
    SCA
    ALB
    NMOM
    0   P1  P2  P3  P4  P5  P5
    1   ..  ..  ..  ..  ..  ..
    2   ..  ..  ..  ..  ..  ..
    NMOM..  ..  ..  ..  ..  ..
    
    P1-P6 - элементы матрицы Эванса
    
    возвращает
    ext, sca, alb, P1
    """
    
    with open(scafile, 'rt') as fin:
        try:
            ext=float(fin.readline())
            sca=float(fin.readline())
            alb=float(fin.readline())
            nmoms=int(fin.readline())
            pmoms = np.loadtxt(fin)
            P1 = np.array(pmoms[:,1])
            #P1 = np.array([pmoms[i,1]/(2*i+1) for i in range(nmoms)])
            return ext, sca, alb, P1
        except:
            print("Ошибка чтения файла %s"%scafile)
        
mu32=np.array([0.024350292663424432509,
    0.0729931217877990394495,0.1214628192961205544704,
    0.1696444204239928180373,0.2174236437400070841497,
    0.264687162208767416374,0.311322871990210956158,
    0.3572201583376681159504,0.4022701579639916036958,
    0.446366017253464087985,0.489403145707052957479,
    0.531279464019894545658,0.5718956462026340342839,
    0.611155355172393250249,0.6489654712546573398578,
    0.6852363130542332425636,0.7198818501716108268489,
    0.7528199072605318966119,0.7839723589433414076102,
    0.8132653151227975597419,0.8406292962525803627517,
    0.8659993981540928197608,0.8893154459951141058534,
    0.9105221370785028057564,0.9295691721319395758215,
    0.9464113748584028160625,0.9610087996520537189186,
    0.9733268277899109637419,0.9833362538846259569313,
    0.9910133714767443207394,0.9963401167719552793469,
    0.9993050417357721394569])
theta32=np.rad2deg(np.arccos(mu32))
phi64=np.r_[np.zeros(32), np.ones(32)*180]
uphi64=np.cos(np.deg2rad(phi64))
theta64=np.r_[theta32,theta32[::-1]]


def run_model(sza0, layfile):
    cumtau, pmoms, ssas = load_atmos(layfile)
    
    I = get_ss_downward_intensity(phi64, theta64, sza0, cumtau,\
            pmoms, ssas)
    
    return theta64*uphi64, I#*(np.pi/2.0)


def run_model1(sza0, layfile, taua, omegaa, taum):
    cumtau, pmoms, ssas = load_atmos(layfile)
    
    I = LN1(theta64, phi64, sza0, taum, taua, omegaa, pmoms[-2])
    
    return theta64*uphi64, I


def LN1(theta, phi, theta0, taum, taua, oma,pmomsa):
    """
    """
    uphi = np.cos(np.deg2rad(phi))
    umu  = np.cos(np.deg2rad(theta))
    umu0 = np.cos(np.deg2rad(theta0))
    
    I = np.zeros_like(umu)

    PFa = Legendre(pmomsa)
    PFm = Legendre([1.0, 0.0, 0.5])
    mu = np.linspace(-1, 1, 1000)
    print(np.trapz(PFa(mu), mu))
    if uphi.shape != umu.shape:
        raise Exception("Arrays uv and uphi must have same shape")

    for i,e in enumerate(umu):
        uth0 = e*umu0+\
                    np.sqrt(1.0-e*e)*np.sqrt(1.0-umu0*umu0)*uphi[i]
        #print(uth0, e)
        I[i] = 1.0/(4.0*np.abs(umu[i]))*(taum*PFm(uth0)+taua*oma*PFa(uth0))

    return I/np.pi

    
def test():
    cumtau, pmoms, ssas = load_atmos('layfile.lay')
    
    # чтобы сделать расчет для плоскости солнечного вертикала,
    # создадим два вектора: phi и theta и разместим их так
    #phi = np.r_[np.zeros(32), np.ones(32)*180]
    #uphi = np.cos(np.deg2rad(phi))
    #theta = np.r_[np.linspace(89,1,32),np.linspace(1,89,32)]
    #theta = np.r_[theta,thetaA[::-1]]
    theta0 = 70
    I = get_ss_downward_intensity(phi64, theta64, theta0, cumtau,\
            pmoms, ssas)
    
    import pylab as plt
    
    #print(np.c_[uphi64*theta64, I])
    plt.plot(uphi64*theta64, I)
    
    plt.show()
    
        