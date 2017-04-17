# Overview

*LSSEFT* is a tool to perform calculations in the _effective field theory of large scale structure_. For references to the original literature, see the [References](#references) section below. Currently, *LSSEFT* can perform the 1-loop calculations needed for the dark matter power spectrum in real and redshift-space. It can decompose the redshift-space power spectrum into Legendre modes, and it implements a resummation prescription that accounts for damping of the baryon acoustic oscillations due to displacements.

Future versions of *LSSEFT* may support more features. If you would like to contribute to this effort then any patches, corrections or enhancements are very welcome. The best way to start is to fork this repository and submit enhancements via [pull requests](https://help.github.com/articles/about-pull-requests).

This version of *LSSEFT* has been written by David Seery at the University of Sussex.

# Releases

The current release of *LSSEFT* is 2017.2.
It can be downloaded from the link below.

The same .tar.gz archives for each release are available from the GitHub
[Releases](https://github.com/ds283/LSSEFT/releases)
page, but for citations please use the
[zenodo.org](https://zenodo.org) DOI.

- 2017.2 (17 April 2017) Source code
[![DOI](missing)](missing)

- 2017.1 (8 April 2017) Source code [![DOI](https://zenodo.org/badge/85569146.svg)](https://zenodo.org/badge/latestdoi/85569146)

# How to install and use *LSSEFT*

*LSSEFT* will build and run on Linux or macOS. Windows is not
currently a supported platform.

To obtain the source code, either clone this repository
or download a tar archive for the most recent release.
The build process is managed by [CMake](https://cmake.org)
and requires at least version 3.0.
It is easiest to build in a separate directory; the name doesn't
matter, and `build` is as good as any other.
Create this directory at the top level, change into it, and then run
`cmake` specifying that you wish to build with the `Release` configuration:
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
```
By default, CMake will build *LSSEFT* with the default toolchain
for your platform.
If you wish to build with a non-standard toolchain then you can specify
the C and C++ compilers here. For example, to use the Intel compiler
you should invoke `cmake` using
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc
```
*LSSEFT* has been tested to build correctly with
[gcc](https://gcc.gnu.org),
[Clang](https://clang.llvm.org)
and the Intel compiler.

## Dependencies

*LSSEFT* depends on a number of other libraries:

* An MPI library. Normally [OpenMPI](https://www.open-mpi.org) is a good
choice. Depending on your system it may be installed already.
If not, it is usually easy to install,
either via a standard package management framework (on Linux) or
[MacPorts](https://www.macports.org)/[Homebrew](https://www.macports.org)
(on macOS).

* The [OpenSSH](https://www.openssh.com) library,
which is used to compute MD5 hashes. This is typically available under
the same package management system used for your MPI library.

* The [Boost](http://www.boost.org) library, with the MPI module
enabled and compiled against whatever MPI library you wish to use.

These dependencies are very common, or large, and are best installed
as system-wide shared resources. For this reason *LSSEFT* does not
bundle them: it expects them to be pre-installed.
They should be discoverable via one of CMake's standard
search strategies (eg. via `pkg-config`, or installation in a
standard location)
before you invoke CMake.
If it is unable to locate them, the build
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

The calculations performed by *LSSEFT* are controlled by
the `master_controller::execute()` method. To customize
*LSSEFT* for your own use, you would normally supply
a new implementation of this method; the other
methods in `master_controller` usually do not need to be changed.
There are two bundled implementations, in `controller/work_functions`.
One of these is for the Planck2015 cosmology, and the other is for
the [MDR1](https://www.cosmosim.org/cms/simulations/mdr1/)
model. These are the cosmologies that were used in the first
*LSSEFT* paper; see [below](#how-to-cite-lsseft).

In *LSSEFT* the desired computation is modelled
by setting up C++ objects that describe the calculation.
Most of these objects require a specification of the background
cosmological parameters.
For example, in the MDR `master_controller::execute()`
we use
```C++
// fix the background cosmological model
FRW_model cosmology_model(MDR1::name, MDR1::omega_m, MDR1::omega_cc, MDR1::h, MDR1::T_CMB, MDR1::Neff,
                          MDR1::f_baryon, MDR1::z_star, MDR1::z_drag, MDR1::z_eq, MDR1::Acurv, MDR1::ns, MDR1::kpiv);
```
Practical calculations require setting up a grid of wavenumbers
at which the one-loop integrals (and hence power spectra) are to be
computed, and also a
grid of redshifts at which we wish to sample the results.
*LSSEFT* provides a `stepping_range<>` object to make it easy to
construct suitable ranges, spaced linearly or logarithmically:
```C++
    // set up a list of redshifts at which to sample the late-time growth functions; we need the z=50 point
    // to define where the integration starts
    stepping_range<double> z50(50.0, 50.0, 0, 1.0, spacing_type::linear);
    stepping_range<double> z0(0.0, 0.0, 0, 1.0, spacing_type::linear);
    stepping_range<double> z025(0.25, 0.25, 0, 1.0, spacing_type::linear);
    stepping_range<double> z05(0.5, 0.5, 0, 1.0, spacing_type::linear);
    stepping_range<double> z075(0.75, 0.75, 0, 1.0, spacing_type::linear);
    stepping_range<double> z1(1.0, 1.0, 0, 1.0, spacing_type::linear);
    auto lo_redshift_samples = z0 + z025 + z05 + z075 + z1 + z50;

    // set up a list of UV cutoffs, measured in h/Mpc, to be used with the loop integrals
    stepping_range<Mpc_units::energy> UV_cutoffs(1.4, 1.4, 0, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);

    // set up a list of IR cutoffs, measured in h/Mpc, to be used with the loop integrals
    stepping_range<Mpc_units::energy> IR_cutoffs(1E-4, 1E-4, 0, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);

    // set up a list of k at which to compute the loop integrals
    stepping_range<Mpc_units::energy> loop_k_samples(0.005, 1.0, 500, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);

    // set up a list of IR resummation scales, measured in h/Mpc
    stepping_range<Mpc_units::energy> IR_resummation(1.4, 1.4, 0, 1.0 / Mpc_units::Mpc, spacing_type::linear);
```
Notice that
ranges can be composed to build aggregate ranges, as for
`lo_redshift_samples`, which is composed from six individual atomic ranges.
(Currently these have to be set up as shown, but in future it would be
nice to allow conversion from a `std::initializer_list`.)
These subranges don't have to be contiguous, or to use the same spacing.
The earliest redshift is taken as the starting point, and the linear
growth factor `D(z)` is normalized to unity at that time.

The range `loop_k_samples` will define the grid of k-modes at which
we sample the integrals.
The other ranges used here are `UV_cutoffs`, `IR_cutoffs`,
and `IR_resummation`, which defined parameters controlling how these
integrals are performed. If you specify multiple values then
*LSSEFT* will produce results for each possible combination.
This is a simple way to check for cutoff dependence of the loop integrals,
or even the final power spectrum, either in the ultraviolet or
infrared.

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
// defauts to qmin = 10 Mpc/h, qmax = 300 Mpc/h
MatsubaraXY_params XY_params;
std::unique_ptr<MatsubaraXY_params_token> XY_tok = dmgr.tokenize(XY_params);

// set up parameters for loop integral
loop_integral_params loop_params;
std::unique_ptr<loop_integral_params_token> loop_tok = dmgr.tokenize(loop_params);
```
The parameters of the filter may need to be changed, if the default
set do not give good results. You may also wish to change the
`qmin` and `qmax` set by `MatsubaraXY_params`.

The current version of `master_controller.cpp` uses these
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
* `--EdS-mode`: use the Einstein--de Sitter approximation for
growth functions and growth factors. In this approxmation the growth
functions D_i(z) can be written as powers of the
linear growth factor D(z),
and the growth factors f_i(z) can be written as multiples of the linear
growth factor f(z).

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

- _The matter power spectrum in redshift space using effective ﬁeld theory_, Lucía Fonseca de la Bella, Donough Regan, David Seery and Shaun Hotchkiss. arXiv:xxxx.xxxxx DOI:xxx

# References

The effective field theory of large-scale structure was proposed by

- Baumann, Nicolis, Senatore & Zaldarriaga, _Cosmological Non-Linearities as an Effective Fluid_, JCAP *1207* (2012) 051 [DOI:10.1088/1475-7516/2012/07/051](http://dx.doi.org/10.1088/1475-7516/2012/07/051), [arXiv:1004.2488](http://arxiv.org/abs/arXiv:1004.2488)

The theory was subsqeuently developed in

- Carrasco, Hertzberg & Senatore, _The Effective Field Theory of Cosmological Large Scale Structures_, JHEP *1209* (2012) 082 [DOI:10.1007/JHEP09(2012)082](http://dx.doi.org/10.1007/JHEP09(2012)082) [arXiv:1206.2926](http://arxiv.org/abs/arXiv:1206.2926)

- Carrasco, Foreman, Green & Senatore, _The Effective Field Theory of Large Scale Structures at Two Loops_, JCAP *1407* (2014) 057 [DOI:10.1088/1475-7516/2014/07/057](http://dx.doi.org/10.1088/1475-7516/2014/07/057) [arXiv:1310.0464](http://arxiv.org/abs/arXiv:1310.0464)

and other publications; for a fuller list, see the bibliography of
[the paper linked above](#how-to-cite-lsseft). Details of the velocity power spectrum and its renormalization, which is important in redshift space, were studied by Mercolli & Pajer

- Mercolli & Pajer, _On the velocity in the Effective Field Theory of Large Scale Structures_, JCAP *1403* (2014) 006 [DOI:10.1088/1475-7516/2014/03/006](http://dx.doi.org/10.1088/1475-7516/2014/03/006) [arXiv:1307.3220](http://arxiv.org/abs/arXiv:1307.3220)

The importance of a resummation scheme to damp the amplitude of baryon acoustic oscillations was discussed by Senatore & Zaldarriaga

- Senatore & Zaldarriaga, _The IR-resummed Effective Field Theory of Large Scale Structures_, JCAP *1502* (205) 013 [DOI:10.1088/1475-7516/2015/02/013](http://dx.doi.org/10.1088/1475-7516/2015/02/013) [arXiv:1404.5954](http://arxiv.org/abs/arXiv:1404.5954)

The particular resummation scheme implemented by *LSSEFT* was proposed by Vlah et al.

- Vlah, Seljak, Chu & Feng, _Perturbation theory, effective field theory, and oscillations in the power spectrum_, JCAP *1603* (2016) 057 [DOI:10.1088/1475-7516/2016/03/057](http://dx.doi.org/10.1088/1475-7516/2016/03/057) [arXiv:1509.02120](http://arxiv.org/abs/arXiv:1509.02120)

# Acknowledgments

Development of the *LSSEFT* code has been supported by the grant _Precision tests of the inflationary scenario_, funded by the European Union's Seventh Framework Programme (FP/2007–2013) and ERC Grant Agreement No. 308082.

Some development was additionally funded by the UK Science and Technology Facilities Council via grants ST/I000976/1 and ST/L000652/1, which funded the science programme at the University of Sussex Astronomy Centre from April 2011–March 2014 and April 2014–March 2017, respectively.
