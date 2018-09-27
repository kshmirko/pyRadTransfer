from numpy.distutils.core import Extension


ssrt = Extension(name='_ssrt',
                sources = ['src/ssrt/SINSCA.f',
                            'src/ssrt/GETMOM.f',
                            'src/ssrt/ErrPack.f'])

disort=Extension(name='_disort',
                sources = ['src/disort/BDREF.f',
                            'src/disort/ErrPack.f',
                            'src/disort/GETMOM.f',
                            'src/disort/LINPAK.f',
                            'src/disort/PRTFIN.f',
                            'src/disort/RDI1MACH.f',
                            'src/disort/DISORT.f',
                            'src/disort/Driver.f'])

# disort=Extension(name='_disort',
#                sources = ['src/disort_new/code.f',
#                            'src/disort_new/Driver_new.f'])


rt3=Extension(name='_rt3',
            sources = ['src/rt3/radintg3.f',
                        'src/rt3/radmat.f',
                        'src/rt3/radscat3.f',
                        'src/rt3/radtran3.f',
                        'src/rt3/radutil3.f',
                        'src/rt3/rt3subs.f',
                        'src/rt3/Driver.f',
                        'src/rt3/ErrPack.f',
                        'src/rt3/RDI1MACH.f',
                        'src/rt3/MIEV0.f',
                        'src/rt3/MIEDRV.f',
                        'src/rt3/wis2ev.f',
                        'src/rt3/TRAPZ.f',
                        'src/rt3/MIEDIST.f',
                        'src/rt3/RAYWIS.f',
                        'src/rt3/RAYEV.f',
                        'src/rt3/MKSCTFL.f',
                        'src/rt3/GETTAUM.f',
                        'src/rt3/MKLYFL.f',
                        'src/rt3/MKDISTR.f',
                        'src/rt3/MKFILES.f',
                        'src/rt3/SCTFLE.f',
                        'src/rt3/Run.f',
                        'src/rt3/RUNRT3.f',
                        'src/rt3/LINTERP.f',
                        'src/rt3/RDMEAS.f',
                        'src/rt3/LEGVAL.f',
                        'src/rt3/RunSub.f',
                        ])


# miev0=Extension(name='_miev0',
#             sources = ['src/miev0/ErrPack.f',
#                         'src/miev0/RDI1MACH.f',
#                         'src/miev0/MIEV0.f',
#                         'src/miev0/MIEDRV.f',
#                         'src/miev0/wis2ev.f',
#                         'src/miev0/TRAPZ.f',
#                         'src/miev0/MIEDIST.f',
#                         'src/miev0/RAYWIS.f',
#                         'src/miev0/RAYEV.f',
#                         'src/miev0/MKSCTFL.f',
#                         'src/miev0/GETTAUM.f',
#                         'src/miev0/MKLYFL.f'])

# atmos=Extension(name='_stdatm',
#             sources = ['src/atmos/STDATM.f'])


if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'ssrtmod',
          description   = "Single scattering radiative transfer code",
          author        = "Konstantin Shmirko",
          author_email  = "shmirko.konstantin(at)gmail.com",
          version       = "0.1.0",
          license       = "GPL",
          packages      = ['ssrt', 'disort','rt3', 
                            #'atmos', 'miev0',
                            ],
          ext_modules = [ssrt, disort, rt3, 
                          #miev0, atmos,
                          ],
          )
