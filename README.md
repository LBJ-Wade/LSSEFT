# Overview

*LSSEFT* is a tool to perform calculations in the _effective field theory of large scale structure_. Currently it can perform the 1-loop calculations needed for the dark matter power spectrum in real and redshift-space. It can decompose the redshift-space power spectrum into Legendre modes, and it implements a resummation prescription that accounts for damping of the baryon acoustic oscillations due to displacements.

Future versions of *LSSEFT* may support more features. If you would like to contribute to this effort then any patches, corrections or enhancements are very welcome. The best way to start is to fork this repository and submit enhancements via [pull requests](https://help.github.com/articles/about-pull-requests).

This version of *LSSEFT* has been written by David Seery at the University of Sussex.

# Releases

The current release of *LSSEFT* is 2017.1, and can be identified via a DOI linking to a deposit at [zenodo.org](https://zenodo.org). The same .tar.gz archives for each release are available directly from GitHub, but for citations please use the [zenodo.org](https://zenodo.org) versions.

- 2017.1 (x March 2017) Source code DOI:xxx

# How to install and use *LSSEFT*

To obtain the source code, either clone this repository
or download a release as a tar archive. The build process is managed
by [CMake](https://cmake.org) and requires at least version 3.0.
It is easiest to build in a separate directory. The name doesn't
matter, and `build` is as good as any other.
Create this directory at the top level, change into it, and then run
`cmake` specifying that you wish to build with the `Release` configuration:
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
```
If you wish to build with a non-standard toolchain then you can specify
the C and C++ compilers here. For example, to use the Intel compiler
you should invoke `cmake` using
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc
```

## Dependencies

*LSSEFT* depends on a number of other libraries:

* An MPI library. Normally [OpenMPI](https://www.open-mpi.org) is a good
choice. Depending on your system it is usually easy to install,
either via a standard package management framework (on Linux) or
[MacPorts](https://www.macports.org)/[Homebrew](https://www.macports.org)
(on macOS).

* The [OpenSSH](https://www.openssh.com) library,
which is used to compute MD5 hashes. This is typically available under
the same package management systems.

* The [Boost](http://www.boost.org) library, with the MPI module
enabled and compiled against whatever MPI library you wish to use.

These dependencies should be installed on your system before you
invoke `cmake`. If `cmake` is unable to locate them, the build
configuration process will fail with an error.

*LSSEFT* also depends on two less common libraries. Since these are unlikely
to be installed on your system, the build process will download
and build them automatically:

* The [SPLINTER](https://github.com/bgrimstad/splinter) library
for calculating B-splines. This is distributed under the
[Mozilla Public License version 2](https://www.mozilla.org/en-US/MPL/2.0/).

* The [Cuba library](http://www.feynarts.de/cuba/) for multidimensional integration. Cuba is distributed under the [GNU Lesser Public General License](http://www.gnu.org/licenses/lgpl.html).

## Build

Once configuration is complete, build *LSSEFT* by invoking `make`.
There are many source files and it is faster to build in parallel
by specifying the `-j` switch.
For example, to launch `make` with 4 processes,
use:
```bash
make lsseft -j4
```
*LSSEFT* uses C++14 features and expects a C++14-aware compiler.
Any recent version of Clang, gcc or the Intel compiler should work
successfully.

## Performing calculations

For most calculations you will wish to modify the
`master_controller.cpp` file, which controls what is actually computed
by setting up C++ objects that model the calculation.
Most of these objects require a specification of the background
cosmological parameters. These are currently set up to match the
[MDR1](https://www.cosmosim.org/cms/simulations/mdr1/)
cosmology, because this was the choice used in our paper:
```C++
// fix the background cosmological model
// here, that's taken to have parameters matching the MDR1 simulation
FRW_model cosmology_model(MDR1::name, MDR1::omega_m, MDR1::omega_cc, MDR1::h, MDR1::T_CMB, MDR1::Neff,
                          MDR1::f_baryon, MDR1::z_star, MDR1::z_drag, MDR1::z_eq, MDR1::Acurv, MDR1::ns, MDR1::kpiv);
```
To do a real calculation requires setting up a grid of wavenumbers
at which the one-loop integrals (and hence power spectra) are to be
computed, and a grid of redshifts at which we wish to sample the results.
*LSSEFT* provides a `stepping_range<>` object to make it easy to
construct suitable ranges, spaced linearly or logarithmically:
```C++

    // set up a list of redshifts at which to sample the late-time growth functions
    stepping_range<double> lo_redshift_samples(0.0, 50.0, 250, 1.0, spacing_type::logarithmic_bottom);

    // set up a list of UV cutoffs, measured in h/Mpc, to be used with the loop integrals
    stepping_range<Mpc_units::energy> UV_cutoffs(1.4, 1.4, 0, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);

    // set up a list of IR cutoffs, measured in h/Mpc, to be used with the loop integrals
    stepping_range<Mpc_units::energy> IR_cutoffs(1E-4, 1E-4, 0, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);

    // set up a list of k at which to compute the loop integrals
    stepping_range<Mpc_units::energy> loop_k_samples(0.005, 1.0, 500, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);

    // set up a list of IR resummation scales, measured in h/Mpc
    stepping_range<Mpc_units::energy> IR_resummation(1.4, 1.4, 0, 1.0 / Mpc_units::Mpc, spacing_type::linear);
```
Ranges can be composed to build aggregate ranges. These don't have to
be contigious, or to use the same spacing.

Some parts of the computation are parametrized, such as the filtering
of an initial power spectrum into 'wiggle' and 'no-wiggle' components.
Mostly, it's fine to stick with default values for these:
```C++
// SET UP PARAMETERS

// set up parameters for filter
Pk_filter_params filter_params(0.25, 0.05 / Mpc_units::Mpc, 0.02);
std::unique_ptr<filter_params_token> filter_tok = dmgr.tokenize(filter_params);

// set up parameters for growth function
growth_params Df_params(this->arg_cache.use_EdS());
std::unique_ptr<growth_params_token> growth_tok = dmgr.tokenize(Df_params);

// set up parameters for Matsubara X&Y integral
MatsubaraXY_params XY_params;
std::unique_ptr<MatsubaraXY_params_token> XY_tok = dmgr.tokenize(XY_params);

// set up parameters for loop integral
loop_integral_params loop_params;
std::unique_ptr<loop_integral_params_token> loop_tok = dmgr.tokenize(loop_params);
```
The current version of `master_controller.cpp` then uses these
ranges and parameters to perform a computation of the one-loop
real-space and redshift-space power spectra. It decomposes the redshift-space power spectrum into multipoles and applies a resummation
prescription.

The current `master_controller.cpp` accepts some command-line switches:

* `--help`: Print brief usage description
* `--verbose` or `-v`: Emit verbose status messages
* `--no-colour` or `--no-color`: Suppress colourized terminal messages
* `--database` or `-d`: Specify database file [mandatory]
* `--initial-powerspectrum` or `-i`: Specify linear power spectrum to be used when computing the one-loop integrals. *LSSEFT* remembers which linear power spectrum was used for each result in the database, and computes the MD5 hash to be sure that you cannot inadvertently change the powerspectrum file without being aware that it may no longer be compatible with the database [mandatory for power-spectrum calculations]
* `--final-powerspectrum` or `-f`: Optionally specify a second linear power spectrum to be used to compute the tree-level part of each one-loop result. This is often slightly more accurate than using the growth factors computed by *LSSEFT* to rescale the initial power spectrum to the final redshift. (For example, *LSSEFT* ignores neutrinos and the influence of any residual radiation component at the initial time, which would be included in the power spectrum computed by a Boltzmann code.)

## Known issues

* Database performance is a bit slow, especially during inserts. This was
fixed to some extent by
[31d5630b](https://github.com/ds283/LSSEFT/commit/31d5630b49a652418276869c9029bb8ac8549d9b), which adjusted the database schema
so that SQLite was not required to perform so much checking during
inserts. Although this gave a significant
performance boost, there is still a residual bottleneck.
Most likely the problem is that *LSSEFT* dumps each MPI payload
into the database individually.
For the more complex calculations
(computation of the resummed real-space power spectrum,
computation of the multipole power spectrum)
each (k, z) pair is an individual MPI payload, so there are lots of them
and each is individually quite small. The result is that we
end up performing a lot of SQLite transactions, which is known to be
a performance bottleneck.
It would probably be better to buffer a number of these payloads in memory
before performing a single bulk insert wrapped in a single transaction.

# Licensing

*LSSEFT* is distributed under the
[GNU General Public License version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html),
or (at your option) any later version. This license is bundled with the source code as LICENSE.txt.

# How to cite *LSSEFT*

Further development of *LSSEFT* depends on demonstrating its usefulness to funding agencies. If you use *LSSEFT* to produce numerical results supporting your research then we would appreciate a citation to the paper

- _The matter power spectrum in redshift space using eﬀective ﬁeld theory_, Lucía Fonseca de la Bella, Donough Regan, David Seery and Shaun Hotchkiss. arXiv:xxxx.xxxxx DOI:xxx

# Acknowledgments

Development of the *LSSEFT* code has been supported by the grant _Precision tests of the inflationary scenario_, funded by the European Union's Seventh Framework Programme (FP/2007–2013) and ERC Grant Agreement No. 308082.

Some development was additionally funded by the UK Science and Technology Facilities Council via grants ST/I000976/1 and ST/L000652/1, which funded the science programme at the University of Sussex Astronomy Centre from April 2011–March 2014 and April 2014–March 2017, respectively.
