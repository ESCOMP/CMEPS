# Unit tests for `cesm/flux_atmocn`

pFUnit unit tests for the CESM atmosphere/ocean bulk flux schemes. The build needs only
a Fortran compiler, CMake and pFUnit 4. It does not need CIME, ESMF, PIO or MPI.

## Running

### Derecho

```bash
cesm/flux_atmocn/test/run_tests_derecho.sh          # Intel ifx (default)
cesm/flux_atmocn/test/run_tests_derecho.sh gnu      # gfortran
```

The script loads the compiler, `cmake` and `pfunit/4.14.0` modules, then configures,
builds and runs the tests in `./build-flux-<compiler>`. An optional second argument sets a
different build directory.

### Any machine

You need gfortran 11 or newer, CMake 3.18 or newer, and a pFUnit 4 install. To build pFUnit
yourself (one time, about 5 minutes):

```bash
git clone --depth 1 -b v4.14.0 https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
cmake -S pFUnit -B pFUnit/build -DSKIP_MPI=YES -DSKIP_OPENMP=YES -DCMAKE_INSTALL_PREFIX=$HOME/pfunit
cmake --build pFUnit/build -j 8 --target install
```

Then, from the top of the CMEPS repository:

```bash
cmake -S cesm/flux_atmocn/test -B build-flux -DCMAKE_PREFIX_PATH=$HOME/pfunit
cmake --build build-flux -j 8
ctest --test-dir build-flux --output-on-failure
```

This runs every unit test. Add `-V` to the `ctest` command to see the full pFUnit output.

## CESM_share

Configuring the tests downloads the head of the [CESM_share](https://github.com/ESCOMP/CESM_share)
`main` branch. To pin a specific version instead, pass `-DCESM_SHARE_TAG=<tag/sha>` to
`cmake`. Only `shr_kind_mod`, `shr_const_mod` and
`shr_wv_sat_mod` are compiled from it. To use a local checkout instead (for example, offline
or to test in-progress share changes):

```bash
cmake ... -DFETCHCONTENT_SOURCE_DIR_CESM_SHARE=/path/to/CESM_share
```

`shr_sys_mod` is replaced by `stubs/shr_sys_mod_stub.F90`, because the real one needs ESMF
and MPI. In the stub, `shr_sys_abort` raises a pFUnit exception and then returns, so abort
paths can be tested with `@assertExceptionRaised`.

## Layout

| Path | Contents |
|---|---|
| `CMakeLists.txt` | Builds the code under test into one library and one test executable per `.pf` file |
| `stubs/` | Replacements for CESM_share modules that need ESMF or MPI |
| `helpers/flux_test_fixtures.F90` | Shared setup/teardown, canned input states, tolerances |
| `unit_tests/` | pFUnit test files (`.pf`) for the flux schemes (Large, COARE, ...) |
| `run_tests_derecho.sh` | Build-and-run wrapper for Derecho |

## Adding tests

- Drop a `.pf` file into a subdirectory of `test/` (for example `unit_tests/test_flux_ua.pf`, or
  a new directory such as `wtracers/`). CMake picks it up on the next build, and it runs
  with all the other tests. Test file names must be unique across subdirectories.
- Each test module should call `fixture_init()` in its `@before` routine and
  `fixture_reset()` in its `@after` routine. `fixture_init` sets up `shr_wv_sat_mod` the way
  the CESM driver does. `fixture_reset` restores the `shr_flux_mod` settings, so a test that
  changes them (for example with `shr_flux_adjust_constants`) does not affect other tests.
- If a scheme gains a new dependency, add the source file to `add_library(flux_under_test ...)`
  in `CMakeLists.txt`. If that dependency needs ESMF, MPI or PIO, write a stub in `stubs/`.
