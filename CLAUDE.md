# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this project does

Fills land (missing/NaN) points in a global sea-surface-temperature (SST) field so the
result is a smooth, gap-free global grid. The approach:

1. Load the SST grid (lat/lon/time/depth) from `SST_1day.nc`.
2. Pad the southern edge with 39 rows (−89.75° to −80.25° at 0.25° spacing) so the grid
   reaches the South Pole.
3. Build a land mask (`true` = NaN/land).
4. Compute a zonal (row) mean over ocean points, copy the first ocean row's mean down to
   the pole (`SP_cutoff`), then smooth it with a ±8-row running mean.
5. Seed land points with that smoothed zonal mean.
6. Relax the land points by explicit diffusion (`Filled_SST += CFL * dt * Laplace(...)`),
   updating only masked points, for `Nsteps = 3000` (CFL 0.01 for the first 100 steps,
   then 0.4). Abort if any value exceeds |100|.
7. Write `Filled_SST`, `lat`, `lon` to `filled_lands_file.nc` (NetCDF-4).

## Files

- `landfil.cpp` — `main()`: the whole pipeline above, plus `pad_top` and `concat_vectors` helpers.
- `laplace.hpp` / `laplace.cpp` — `Laplace(f, cosLat, dlat, dlon)`: spherical Laplacian on a
  lat/lon grid (radians), radius `R_E = 6371000 m`, periodic in longitude, interior rows only
  (first and last rows are left 0).
- `functions.hpp` — **not in the repo**. It must provide the `dataset` class
  (`load_time/load_depth/load_latitude/load_longitude`, members `Nlat`, `Nlon`, `latitude`,
  `longitude`) and bring `std::vector` into scope as `vector` (`landfil.cpp` uses
  unqualified `vector`).

## Build

No build system yet. Requires a C++17 compiler and the NetCDF C library:

```sh
g++ -std=c++17 -O2 landfil.cpp laplace.cpp -lnetcdf -o landfil
./landfil   # expects SST_1day.nc in the working directory
```

Install NetCDF on Debian/Ubuntu with `apt-get install libnetcdf-dev`.

## Current state: does not compile

`landfil.cpp` is a work-in-progress port (likely from Python) and has known errors.
Fix these before anything else:

- `functions.hpp` is missing (see above).
- `R_E` is defined in `laplace.hpp`, `laplace.cpp`, and `landfil.cpp`, which causes
  redefinition errors. Keep one definition (e.g. `inline constexpr double R_E` in the header)
  and have `laplace.cpp` include `laplace.hpp`.
- `Nlon` is declared twice in `main` (and `Nlat`/`Nlon` again as `int` in the output section).
- Mask loop reads `SST[i][j]` (unpadded, wrong size) instead of `padded_SST[i][j]`; `SST` is
  also never filled from the input file. No SST load call exists yet.
- `double dlat = Nlat_extended[1] - Nlat_extended[0];` should use `extended_lat`.
- `int SP_cutoff = 0` is missing `;`; `True` should be `true`; `sp_cutoff` vs `SP_cutoff`
  casing mismatch.
- `static_cast<int> i-8` needs parentheses: `static_cast<int>(i) - 8`.
- The time-step `for` loop has no braces, so only the `CFL_val` line is in the loop and
  `Istep` is out of scope afterwards. The step loop must wrap the RHS, update, stability check,
  and progress print; the NetCDF write goes after it.
- `return 0` at the end is missing `;`.
- `new_lat` is in radians but `source_data.latitude` units are unverified — make sure both
  are radians before concatenating (`Laplace` and `cos()` expect radians).
- NetCDF calls ignore return codes; check them (`nc_strerror`).

## Conventions

- Grids are `std::vector<std::vector<double>>` indexed `[lat][lon]`; latitude ascends from the
  South Pole.
- Angles in radians internally; distances in metres.
- Missing data is `NAN`; test with `std::isnan`.
- Keep the existing style: 4-space indentation, `std::`-qualified names in new code.
