from numpy.distutils.core import Extension


ssrt = Extension(name='_ssrt',
                sources = ['src/ssrt/SINSCA.f',
                            'src/ssrt/GETMOM.f',
                            'src/ssrt/ErrPack.f'])

# disort=Extension(name='_disort',
#                 sources = ['src/disort/BDREF.f',
#                             'src/disort/ErrPack.f',
#                             'src/disort/GETMOM.f',
#                             'src/disort/LINPAK.f',
#                             'src/disort/PRTFIN.f',
#                             'src/disort/RDI1MACH.f',
#                             'src/disort/DISORT.f',
#                             'src/disort/Driver.f'])

disort=Extension(name='_disort',
                sources = ['src/disort_new/code.f',
                            'src/disort_new/Driver.f'])


rt3=Extension(name='_rt3',
            sources = ['src/rt3/radintg3.f',
                        'src/rt3/radmat.f',
                        'src/rt3/radscat3.f',
                        'src/rt3/radtran3.f',
                        'src/rt3/radutil3.f',
                        'src/rt3/rt3subs.f'])

miev0=Extension(name='_miev0',
            sources = ['src/miev0/ErrPack.f',
                        'src/miev0/RDI1MACH.f',
                        'src/miev0/MIEV0.f',
                        'src/miev0/MIEDRV.f',
                        'src/miev0/wis2ev.f',
                        'src/miev0/trapz.f',
                        'src/miev0/MIEDIST.f',
                        'src/miev0/RAYWIS.f',
                        'src/miev0/RAYEV.f'])

atmos=Extension(name='_stdatm',
            sources = ['src/atmos/STDATM.f'])


if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'ssrtmod',
          description   = "Single scattering radiative transfer code",
          author        = "Konstantin Shmirko",
          author_email  = "shmirko.konstantin(at)gmail.com",
          version       = "0.1.0",
          license       = "GPL",
          packages      = ['ssrt', 'disort','rt3', 'atmos'],
          ext_modules = [ssrt, disort, rt3, miev0, atmos],
          )
