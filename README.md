# BinSwitch

A tempo2 plugin that switches the binary model of a par file. It refits the par file once, then writes a new par
file that keeps every line of the original except the binary parameters, followed by a rewritten binary block for
the model you ask for.

Target models: `DD`, `DDH`, `DDGR`, `BT`, `T2` and `ELL1`.

## What it writes

| Target | Orbital elements | Extra parameters |
|---|---|---|
| `ELL1` | `A1`, `PB`, `TASC`, `EPS1`, `EPS2` | `SINI`, `M2` |
| `DD`, `BT`, `T2` | `A1`, `PB`, `T0`, `ECC`, `OM` | `SINI`, `M2` |
| `DDH` | `A1`, `PB`, `T0`, `ECC`, `OM` | `H3`, `STIG` |
| `DDGR` | `A1`, `PB`, `T0`, `ECC`, `OM` | `M2`, `MTOT`, `SINI` |

For every target, `OMDOT`, `PBDOT`, `A1DOT`, `EDOT` and `GAMMA` are copied if they are set.

Conversions:

* **ELL1 ↔ T0/ECC/OM**: `EPS1 = ECC sin ω`, `EPS2 = ECC cos ω`, `TASC = T0 − PB ω/360°`, and the inverse.
  This is the usual small-eccentricity approximation. With `--test`, BinSwitch warns if ELL1 is requested
  for ECC > 0.01.
* **DDH**: `STIG = sin i / (1 + cos i)` and `H3 = T☉ M2 STIG³`, from `M2` and `SINI`.
* **DDGR**: `MTOT` comes from the mass function with `M2` and `SINI`, or else from `OMDOT`, assuming GR.

If `M2` or `SINI` cannot be found in the par file, BinSwitch uses `M2 = 1.0` and `SINI = 0.9` and prints a
warning.

## Build

tempo2 looks for a plugin named `<name>_<arch>_plug.t2`. The arch is `linux-gnu` for the
`pulsar_folder_cuda.sif` build. Use `$LOGIN_ARCH` if your tempo2 needs it.

Inside the container:

```
singularity exec -B /bEDD /fpra/timing/01/fazal/singularity_img/pulsar_folder_cuda.sif bash -c '
  g++ -fPIC -c binswitch_plug.C -o binswitch_plug.o -I$TEMPO2/include -I/software/tempo2
  g++ -shared binswitch_plug.o -o binswitch_linux-gnu_plug.t2'
```

## Install

Copy the plugin into a directory on tempo2's plugin search path (running `tempo2` with no arguments lists the
paths):

```
cp binswitch_linux-gnu_plug.t2 $TEMPO2/plugins/
```

Without installing it, `-gr2` loads the plugin from the current directory instead:

```
tempo2 -gr2 binswitch -f in.par in.tim --binary ELL1 --out out.par
```

## Usage

```
tempo2 -gr binswitch -f in.par in.tim --binary <MODEL> --out out.par [--keep-fitflags] [--test]
```

| Option | Meaning |
|---|---|
| `-f in.par in.tim` | Input par and tim files (required) |
| `--binary`, `-b` | Target model. Default: the model in the par file, or `DD` if there is none |
| `--out`, `-o` | Output par file. Default: `<PSRJ>.base.par` |
| `--keep-fitflags` | Keep the fit flags of `A1`, `PB` and the copied parameters (otherwise all are 0) |
| `--test` | Print a round-trip check of the ELL1 ↔ T0/ECC/OM conversion and the small-e warning |

Examples:

```
# DD -> ELL1
tempo2 -gr binswitch -f myDD.par my.tim --binary ELL1 --out myELL1.par

# current solution -> DDH, with the conversion check
tempo2 -gr binswitch -f my.par my.tim --binary DDH --out base_DDH.par --test

# -> DDGR (masses inferred if possible, otherwise defaults with a warning)
tempo2 -gr binswitch -f my.par my.tim --binary DDGR --out myDDGR.par
```

## Known issues

Check every output par with tempo2 before using it.

* **Precision.** Values are formatted with tempo2's `print_longdouble`, which uses `%Lg`, so only 6
  significant digits are written (for example `PB 0.290684`). For a real timing solution this breaks the
  phase connection.
* **Default masses.** If the input has no `M2`/`SINI`, the defaults `M2 = 1.0` and `SINI = 0.9` are written
  into the output, which adds a Shapiro delay the input did not have.
* `H3` and `STIG` in the input are not read, so a DDH input loses its Shapiro measurement.
* The new `T0`/`ECC`/`OM` or `TASC`/`EPS1`/`EPS2` always get fit flag 0, even with `--keep-fitflags`.
* The binary block uses post-fit values, but all other lines are copied from the input file unchanged
  (pre-fit), so the spin and astrometric parameters do not come from the same fit.
* `--test` converts with the same formulas both ways, so it shows 0 differences; it does not measure the error of the
  small-eccentricity approximation.
