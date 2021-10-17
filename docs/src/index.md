# FSPS.jl

FSPS.jl provides access to the [Flexible Stellar Population Synthesis](https://github.com/cconroy20/fsps) (FSPS) Fortran routines (see [Conroy et al. 2009](https://doi.org/10.1088/0004-637X/699/1/486), [Conroy & Gunn 2010](https://doi.org/10.1088/0004-637X/712/2/833)) in Julia, and is roughly equivalent to [python-fsps](https://github.com/dfm/python-fsps) by [Foreman-Mackey et al](https://doi.org/10.5281/zenodo.591505).

## Installation

FSPS.jl requires that FSPS be installed on the system and that the `SPS_HOME` environment variable be set to the filepath of the installation. Then, FSPS.jl may be installed using

```julia-repl
] add FSPS
```

at the Julia REPL.
