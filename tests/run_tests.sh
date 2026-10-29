#!/bin/bash
# BinSwitch tests on simulated TOAs. Run from the repository directory, inside the container, after building:
#   singularity exec -B /bEDD pulsar_folder_cuda.sif bash tests/run_tests.sh
# A fictitious pulsar (tests/J0000+0000.par: DD, small e, OMDOT, EDOT, Shapiro delay) is simulated with tempo2's
# fake plugin; every converted par must predict the same residuals, TOA by TOA.
set -u
REPO=$(cd "$(dirname "$0")/.." && pwd)
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
cp "$REPO"/binswitch_*_plug.t2 "$REPO"/tests/J0000+0000.par "$W"/
cd "$W"
npass=0; nfail=0
ok()   { echo "PASS  $*"; npass=$((npass + 1)); }
bad()  { echo "FAIL  $*"; nfail=$((nfail + 1)); }
bs()   { tempo2 -gr2 binswitch "$@" 2>&1 | grep '^binswitch'; }
res()  { tempo2 -nofit -output general2 -f "$1" "$2" -s "RES {sat} {pre}\n" 2>/dev/null | awk '/^RES/{print $2, $3}'; }
dres() { paste <(res "$1" sim.tim) <(res "$2" sim.tim) | awk '{d=($2-$4)*1e9; if(d<0)d=-d; if(d>m)m=d; n++} END{if(n<100)print "nan"; else printf "%.4f", m}'; }
check() {  # name, par1, par2, tolerance (ns)
  local d; d=$(dres "$2" "$3")
  if awk -v d="$d" -v t="$4" 'BEGIN{exit !(d != "nan" && d <= t)}'; then ok "$1: max |dres| = $d ns (<= $4)"
  else bad "$1: max |dres| = $d ns (tolerance $4)"; fi
}

# simulated TOAs (1 us white noise), from the par without fit flags
awk '{print $1, $2}' J0000+0000.par > sim.par
tempo2 -gr fake -f sim.par -nobsd 1 -ndobs 5 -randha y -start 58000 -end 59000 -rms 0.001 > /dev/null 2>&1
mv sim.simulate sim.tim

# the ELL1 approximation error for this orbit, x^2 n e + x e^2 (see README): ~16 ns
ELL1_TOL=$(awk '$1=="A1"{x=$2} $1=="PB"{pb=$2} $1=="ECC"{e=$2} END{n=2*3.14159265358979/(pb*86400); printf "%.1f", 2*(x*x*n*e + x*e*e)*1e9}' J0000+0000.par)

bs -f J0000+0000.par sim.tim --binary ELL1 --out ell1.par --keep-fitflags > log_ell1.txt
bs -f ell1.par sim.tim --binary DD --out ell1_dd.par --keep-fitflags > /dev/null
bs -f J0000+0000.par sim.tim --binary DDH --out ddh.par > /dev/null
bs -f ddh.par sim.tim --binary DD --out ddh_dd.par > /dev/null
bs -f J0000+0000.par sim.tim --binary DD --out dd_dd.par > /dev/null

check "DD -> ELL1 (ELL1 model accuracy)" J0000+0000.par ell1.par "$ELL1_TOL"
check "DD -> ELL1 -> DD"                  J0000+0000.par ell1_dd.par 0.01
check "DD -> DDH"                         J0000+0000.par ddh.par 0.01
check "DD -> DDH -> DD"                   J0000+0000.par ddh_dd.par 0.01
check "DD -> DD"                          J0000+0000.par dd_dd.par 0.0

# values after DD -> ELL1 -> DD (full precision: relative difference < 1e-15)
for k in F0 PB T0 ECC OM OMDOT EDOT; do
  a=$(awk -v k=$k '$1==k{print $2}' J0000+0000.par); b=$(awk -v k=$k '$1==k{print $2}' ell1_dd.par)
  if awk -v a="$a" -v b="$b" 'BEGIN{d=(a-b)/a; if(d<0)d=-d; exit !(b != "" && d < 1e-15)}'; then ok "round trip $k: $a -> $b"
  else bad "round trip $k: $a -> $b"; fi
done

# every line that the conversion does not touch is copied verbatim
untouched='^(BINARY|F0|PB|T0|ECC|OM|OMDOT|EDOT|TASC|EPS1|EPS2|EPS1DOT|EPS2DOT|M2|SINI|H3|STIG|MTOT) '
if diff <(grep -vE "$untouched" J0000+0000.par) <(grep -vE "$untouched" ell1.par | grep -v '^# --- converted') > /dev/null
then ok "untouched lines copied verbatim"; else bad "untouched lines changed"; fi
if grep -qE '^PB +[0-9.]+ +1 +1.0e-11$' ell1.par; then ok "PB rewritten in place, fit flag and uncertainty kept"; else bad "PB line: $(grep '^PB' ell1.par)"; fi

# fit flags are inherited with --keep-fitflags
if awk '$1=="TASC" && $3==1{f=1} END{exit !f}' ell1.par; then ok "--keep-fitflags: TASC inherits the fit flag of T0"
else bad "--keep-fitflags: TASC fit flag"; fi

# circular orbit without Shapiro delay: ELL1 is exact, and no M2/SINI may appear
grep -vE '^(ECC|OMDOT|EDOT|M2|SINI) ' J0000+0000.par > circ.par
bs -f circ.par sim.tim --binary ELL1 --out circ_ell1.par > /dev/null
check "circular DD -> ELL1" circ.par circ_ell1.par 0.5
if ! grep -qE '^(M2|SINI) ' circ_ell1.par; then ok "no M2/SINI invented"; else bad "M2/SINI written for a par without them"; fi

# warnings and errors
sed 's/^ECC .*/ECC 0.1 1 1e-8/' J0000+0000.par > ecc.par
if bs -f ecc.par sim.tim --binary ELL1 --out ecc_ell1.par | grep -q "WARNING: the ELL1 approximation"; then ok "ELL1 warning for e = 0.1"
else bad "no ELL1 warning for e = 0.1"; fi
bs -f circ.par sim.tim --binary DDGR --out ddgr.par > /dev/null
if [ ! -e ddgr.par ]; then ok "DDGR without masses: refused, nothing written"; else bad "DDGR without masses was written"; fi

echo "passed $npass, failed $nfail"
[ "$nfail" -eq 0 ]
