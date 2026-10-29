# BinSwitch

A tempo2 plugin that switches the binary model of a par file between `DD`, `DDH`, `DDGR`, `BT`, `T2` and `ELL1`.

BinSwitch does **not** fit. It writes the same timing solution in the new parametrisation:

* Lines that the conversion does not involve are copied verbatim, with their fit flags and uncertainties.
* Converted values are written with full precision (21 significant digits).
* Nothing is invented: if a quantity cannot be converted, BinSwitch warns you or stops.

Refit the new par file with tempo2 afterwards.

## What it writes

| Target | Orbit | Shapiro delay / masses |
|---|---|---|
| `ELL1` | `TASC`, `EPS1`, `EPS2` (+ `EPS1DOT`, `EPS2DOT`) | `M2`, `SINI` |
| `DD`, `BT`, `T2` | `T0`, `ECC`, `OM` (+ `OMDOT`, `EDOT`) | `M2`, `SINI` |
| `DDH` | `T0`, `ECC`, `OM` (+ `OMDOT`, `EDOT`) | `H3`, `STIG` |
| `DDGR` | `T0`, `ECC`, `OM` (+ `OMDOT`, `EDOT`) | `M2`, `MTOT` |

* **Parametrisation unchanged** (DD → DDH, BT → DD, …): the orbital lines are kept as they are, and only the
  Shapiro/mass parameters are converted.
* **Parametrisation changed** (ELL1 ↔ the others): the orbital parameters are converted, and `PB` and `F0`
  are rescaled in place, keeping the rest of their line.
* A Shapiro parameter is written only if the input has it or it can be derived. A par file without `M2`/`SINI`
  gets none.
* Converted parameters are written at the end, under `# --- converted by binswitch from BINARY <model> ---`.

## Conversions

### ELL1 ↔ T0/ECC/OM

These are the small-eccentricity relations of the ELL1 model (Lange et al. 2001), including how tempo2's DD
advances ω: ω = ω₀ + k·(true anomaly), with k = ω̇/n.

```
PB_ELL1 = PB / (1 + k)                      w_asc = OM / (1 + k)
TASC    = T0 - w_asc / n                    e_asc = ECC + EDOT (TASC - T0)
EPS1    = e_asc sin(w_asc)                  EPS1DOT = EDOT sin(w_asc) + e_asc wdot cos(w_asc)
EPS2    = e_asc cos(w_asc)                  EPS2DOT = EDOT cos(w_asc) - e_asc wdot sin(w_asc)
F0_ELL1 = F0 (1 + 3/2 d(x EPS1)/dt)
```

The ELL1 → T0/ECC/OM direction is the exact inverse. The `F0` factor exists because the DD delay contains a
term −(3/2)·x·EPS1 that ELL1 leaves out. While EPS1 is constant this term is a constant offset, but with
`EDOT`, `OMDOT` or `A1DOT` it drifts in time, which is equivalent to a change of F0.

### Shapiro delay and masses

* **DDH** (Freire & Wex 2010): `STIG = sin i / (1 + cos i)` and `H3 = T☉ M2 STIG³`, and the inverse.
* **DDGR**: `MTOT` comes from the mass function with `M2` and `SINI`, or else from `OMDOT` assuming GR. If
  neither is possible, BinSwitch stops without writing anything.

## Accuracy

ELL1 is an approximation. Compared with DD, it neglects terms of order **x·e²** in the Roemer delay and
**x²·n·e** in tempo2's inverse-timing correction (x = `A1`, n = 2π/`PB`). BinSwitch prints a warning when
x·e² + x²·n·e exceeds the timing precision (median TOA error / √N). `--test` prints both terms.

Measured on simulated TOAs (`tests/`, x = 5.1 lt-s, PB = 1.53 d, e = 1.2 × 10⁻⁵, with `OMDOT`, `EDOT` and a
Shapiro delay):

| Conversion | Max residual difference |
|---|---|
| DD → ELL1 | 19.8 ns (x²·n·e = 15 ns) |
| DD → ELL1 → DD | 0.0000 ns; every value back to 10⁻¹⁵ |
| DD → DDH, DDH → DD | 0.0000 ns |
| circular DD → ELL1 | 0.1 ns |

The constant −(3/2)·x·EPS1 shifts the absolute pulse phase between DD and ELL1 (79 µs in the test). The
offset of any fit absorbs it. Refit before using the new par to predict absolute phase (folding, polycos).

## Build

tempo2 loads a plugin named `<name>_<arch>_plug.t2`. The arch is `linux-gnu` for the `pulsar_folder_cuda.sif`
build. Use `$LOGIN_ARCH` if your tempo2 needs it.

```
singularity exec -B /bEDD /fpra/timing/01/fazal/singularity_img/pulsar_folder_cuda.sif bash -c '
  g++ -O2 -fPIC -c binswitch_plug.C -o binswitch_plug.o -I$TEMPO2/include -I/software/tempo2
  g++ -shared binswitch_plug.o -o binswitch_linux-gnu_plug.t2'
```

## Install

Copy the plugin into a directory on tempo2's plugin search path. Running `tempo2` with no arguments lists the
search paths.

```
cp binswitch_linux-gnu_plug.t2 $TEMPO2/plugins/
```

Without installing it, `-gr2` loads the plugin from the current directory instead of `-gr`.

## Usage

```
tempo2 -gr binswitch -f in.par in.tim --binary <MODEL> --out out.par [--keep-fitflags] [--test]
```

| Option | Meaning |
|---|---|
| `-f in.par in.tim` | Input par and tim files (required; the TOAs are only used for the accuracy check) |
| `--binary`, `-b` | Target model: `DD`, `DDH`, `DDGR`, `BT`, `T2` or `ELL1`. Default: the model in the par file |
| `--out`, `-o` | Output par file. Default: `<PSRJ>.base.par` |
| `--keep-fitflags` | Converted parameters inherit the fit flags of the parameters they come from (otherwise 0) |
| `--test` | When the parametrisation changes: print the eccentricity, the two ELL1 error terms and the timing precision |

Examples:

```
tempo2 -gr binswitch -f myDD.par my.tim --binary ELL1 --out myELL1.par --keep-fitflags
tempo2 -gr binswitch -f my.par my.tim --binary DDH --out myDDH.par
tempo2 -gr binswitch -f my.par my.tim --binary DDGR --out myDDGR.par
```

## Tests

```
singularity exec -B /bEDD pulsar_folder_cuda.sif bash tests/run_tests.sh
```

The script simulates TOAs for a fictitious pulsar (`tests/J0000+0000.par`) with tempo2's `fake` plugin. It then
converts the par between models and compares the residuals TOA by TOA, the round-trip values, the untouched
lines, the fit flags, the warnings and the errors (19 checks).

## Not converted

* `H4` (ELL1H), `SHAPMAX` (DDS) and `KIN` are kept as they are, with a warning.
* `H3` without `STIG` cannot be converted to `M2`/`SINI`, so it is dropped with a warning.
* DDGR → another model: the post-Keplerian parameters that DDGR derives from GR (`OMDOT`, `GAMMA`, `PBDOT`,
  `SINI`) are not written. Fit them in the new model.
* `BT` keeps `M2`/`SINI`, but the BT model has no Shapiro delay.
* `PBDOT` is not rescaled by 1/(1+k); the difference is second order in k = ω̇/n.
