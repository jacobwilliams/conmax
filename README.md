
### Status

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/conmax.svg)](https://github.com/jacobwilliams/conmax/releases/latest)
[![Build Status](https://github.com/jacobwilliams/conmax/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/conmax/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/conmax/branch/master/graph/badge.svg?token=43HK33CSMY)](https://codecov.io/gh/jacobwilliams/conmax)


This is a work in progress of a refactored version of [CONMAX](http://www.netlib.org/opt/conmax.f) in Modern Fortran.

### Compiling

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that the library and tests cases can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To use `conmax` within your FPM project, add the following to your `fpm.toml` file:
```toml
[dependencies]
conmax = { git="https://github.com/jacobwilliams/conmax.git" }
```

To generate the documentation using [FORD](https://github.com/Fortran-FOSS-Programmers/ford), run:

```
  ford conmax.md
```

### Documentation

The latest API documentation can be found [here](http://jacobwilliams.github.io/conmax/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### License

The conmax source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/conmax/blob/master/LICENSE) (BSD-style).

### References

 * E. H. Kaufman Jr., D. J. Leeming & G. D. Taylor, "An ODE-based approach to nonlinearly constrained minimax problems", Numerical Algorithms, Volume 9, pages 25-37 (1995)
 * Original CONMAX sourcecode at Netlib: http://www.netlib.org/opt/conmax.f