#!/usr/bin/env python

# Joe McGovern, Marine Institute 17.08.2021
# Needed to replace forward slash / with r and backwards slash \ for windows platform

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pyroms_toolbox', parent_package, top_path)
    config.add_subpackage('BGrid_GFDL')
    config.add_subpackage('BGrid_POP')
    config.add_subpackage('BGrid_SODA')
    config.add_subpackage('CGrid_GLORYS')
    config.add_subpackage('seawater')
    config.add_subpackage('Grid_HYCOM')
    config.add_library('_average', sources=[r'src\average.f90']),
    config.add_library('_move_runoff', sources=[r'src\move_runoff.f90']),
    config.add_library('_move_river_t', sources=[r'src\move_river_t.f90']),
    config.add_library('creep', sources=[r'src\creeping_sea.f90']),
    config.add_extension('_average',
                         sources=[r'src\average.f90'],
                         libraries=['_average']
                         )
    config.add_extension('_move_runoff',
                         sources=[r'src\move_runoff.f90'],
                         libraries=['_move_runoff']
                         )
    config.add_extension('_move_river_t',
                         sources=[r'src\move_river_t.f90'],
                         libraries=['_move_river_t']
                         )
    config.add_extension('creep',
                         sources=[r'src\creeping_sea.f90'],
                         libraries=['creep']
                         )
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
