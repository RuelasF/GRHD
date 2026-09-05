# GRHD

A finite-volume general-relativistic hydrodynamics simulator written in Fortran.
The code supports flat, Eddington-Finkelstein, and Kerr-Schild metrics through a
shared 3+1 interface.

## Current reproducible default

The default test is the equatorial Fishbone-Moncrief torus:

- geometry: spherical, logarithmic radial coordinate;
- metric: Kerr-Schild;
- black-hole mass: (M = 1);
- spin: (a = 0.9);
- grid: (400 x 1 x 200);
- radial domain: (r \in [1.2, 40]);
- final time: (t = 1000);
- reconstruction and Riemann solver: WENO5 and HLLE.

The remaining test cases preserve their own metric and output settings. Do not
change those settings to run the torus.

## Build and run

```bash
make clean
make -j
./grhd2
```

Run the available regression tests with:

```bash
make test-hllc
make test-gw
```

The numerical architecture and test case are currently set in
`initialization.f90`. Output is written to
`FM_KS_a0.9d0_data/weno5_hlle/` for the default torus configuration.

## Plot profiles

```bash
python3 plot_fm_profiles.py
```

Use an explicit pattern when plotting a different run:

```bash
python3 plot_fm_profiles.py 'FM_KS_a0.5d0_data/weno5_hlle/*.vtk' \
  --field Density --frames 100 --output figures/fm_a05_density.png --no-show
```

Run `python3 plot_fm_profiles.py --help` for all options.

## Repository hygiene

Compiled files, solver output, checkpoints, logs, editor history, and backups
are intentionally ignored by Git. Keep only source code, parameter definitions,
small reference data, and documentation under version control.

## Validation status

The smooth-wave campaign recovers the expected convergence order for Godunov,
TVD-MC, WENO3, MP5, and WENO5 with HLLC. The HLLC and Finn--Evans geometry and
stress regression tests pass on Zotz.

The principal PPI campaign completed 18 equatorial runs to (t=3000) with
Eddington--Finkelstein/Kerr--Schild backgrounds, six metric/spin configurations, and
WENO3/WENO5/MP5 with HLLE. These are research results in a one-polar-cell
2.5D model; they do not constitute a validated three-dimensional or
self-gravitating calculation. The current gravitational-wave output is a
Finn--Evans weak-field proxy, not a gauge-invariant physical strain.

See `CONTEXTO_GRHD_PARA_PRISM.md` for the complete technical and scientific
status, evidence, limitations, campaign locations, and proposed thesis scope.
