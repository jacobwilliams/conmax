# conmax

Unofficial mirror of CONMAX: http://www.netlib.org/opt/conmax.f


### Status

![Build Status](https://github.com/jacobwilliams/conmax/actions/workflows/CI.yml/badge.svg)

### Compiling

A [FoBiS](https://github.com/szaghi/FoBiS) configuration file (`conmax.fobis`) is provided that can also build the library and examples. Use the `mode` flag to indicate what to build. For example:

  * To build all the examples using gfortran: `FoBiS.py build -f conmax.fobis -mode tests-gnu`
  * To build all the examples using ifort:    `FoBiS.py build -f conmax.fobis -mode tests-intel`
  * To build a static library using gfortran: `FoBiS.py build -f conmax.fobis -mode static-gnu`
  * To build a static library using ifort:    `FoBiS.py build -f conmax.fobis -mode static-intel`

  The full set of modes are: `static-gnu`, `static-gnu-debug`, `static-intel`, `static-intel-debug`, `shared-gnu`, `shared-gnu-debug`, `shared-intel`, `shared-intel-debug`, `tests-gnu`, `tests-gnu-debug`, `tests-intel`, `tests-intel-debug`

  To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```FoBis.py rule --execute makedoc -f conmax.fobis```

### Documentation

The latest API documentation can be found [here](http://jacobwilliams.github.io/conmax/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### License

The conmax source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/conmax/blob/master/LICENSE) (BSD-style).
