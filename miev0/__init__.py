from _miev0 import distr, rayev, wiscombe2evans, makescatfile
from numpy import linspace, pi


def make_scattering_file(fname, midx=1.5-0.0j, r0=0.1, r1=1.0,
                         gamma=-3.5, npts=101,
                         wl=0.532, taua=0.1, taum=0.0, omega=None,
                         nmoms=30):
    """
    По заданным параметрам создает файл с матрицей эванса и оптическими
    характеристиками. 
    """
    tau_tot = taua+taum
    mol_ev = rayev(nmoms)

    r = linspace(r0, r1, npts)
    y = (r/r[0])**gamma
    #print(nmoms)
    aer_wis, ext_i, sca_i, asy_i, vol_i, ierr = distr(r, y, midx, wl,
                        nmoms)
    aer_ev = wiscombe2evans(aer_wis)
    if ierr != 0:
        raise Exception("Неверное значение nmoms="+nmoms)
    omega_a = sca_i/ext_i
    sca_a = taua*omega_a
    tau_tot = taua+taum
    
    sca_tot = taua*omega_a+taum
    
    Evans_tot = (aer_ev*taua*omega_a+mol_ev*taum)/(taua*omega_a+taum)
    omega_tot = sca_tot/tau_tot
    
    with open(fname,"wt") as fout:
        fout.write(f"{tau_tot:12.4e}\n")
        fout.write(f"{sca_tot:12.4e}\n")
        fout.write(f"{omega_tot:12.4e}\n")
        fout.write(f"{nmoms}\n")
        for i in range(nmoms+1):
            fout.write(f"{i:4}"+
            f"{Evans_tot[i,0]:12.4e}{Evans_tot[i,1]:12.4e}"+
            f"{Evans_tot[i,2]:12.4e}{Evans_tot[i,3]:12.4e}"+
            f"{Evans_tot[i,4]:12.4e}{Evans_tot[i,5]:12.4e}\n")
    return tau_tot, sca_tot, omega_tot, Evans_tot


