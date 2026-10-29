// binswitch_plug.C
//
// TEMPO2 graphical plugin:
//   - switches the binary model of a par file between DD, DDH, DDGR, BT, T2 and ELL1
//   - keeps every line of the input par verbatim, except the ones the new model parametrises differently
//   - writes converted values with full longdouble precision
//   - does not fit: the output is the input solution, re-parametrised. Refit it with tempo2 afterwards.
//
// Build/install: see README.md.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <cctype>

#include "tempo2.h"

#ifdef LONGDOUBLE_IS_FLOAT128
#include <quadmath.h>
#endif

using std::string;

static inline string ltrim(const string &s) {
  size_t i = 0;
  while (i < s.size() && std::isspace((unsigned char)s[i])) i++;
  return s.substr(i);
}

// First word of a par line, upper-cased, with tempo2's aliases mapped to one name ("" for comments).
static inline string first_token(const string &line) {
  string t = ltrim(line);
  if (t.empty() || t[0] == '#') return "";
  size_t sp = t.find_first_of(" \t\r\n");
  string tok = (sp == string::npos) ? t : t.substr(0, sp);
  std::transform(tok.begin(), tok.end(), tok.begin(), [](unsigned char c){ return std::toupper(c); });
  if (tok == "E") tok = "ECC";
  if (tok == "XDOT") tok = "A1DOT";
  return tok;
}

static inline longdouble deg2rad_ld(longdouble d) {
  const longdouble pi = (longdouble)3.141592653589793238462643383279502884L;
  return d * (pi / (longdouble)180.0L);
}
static inline longdouble rad2deg_ld(longdouble r) {
  const longdouble pi = (longdouble)3.141592653589793238462643383279502884L;
  return r * ((longdouble)180.0L / pi);
}
static inline longdouble wrap_0_360(longdouble deg) {
  longdouble x = deg;
  while (x < 0)   x += (longdouble)360.0L;
  while (x >= 360) x -= (longdouble)360.0L;
  return x;
}

static inline bool param_is_set(const pulsar &p, int ip) {
  return (p.param[ip].paramSet != nullptr && p.param[ip].paramSet[0] == 1);
}
static inline longdouble param_val0(const pulsar &p, int ip) {
  return p.param[ip].val[0];
}
static inline int param_fitflag0(const pulsar &p, int ip) {
  return (param_is_set(p, ip) && p.param[ip].fitFlag) ? p.param[ip].fitFlag[0] : 0;
}

// Full precision. tempo2's print_longdouble() uses %Lg, i.e. only 6 significant digits.
static string fmt_ld(longdouble v) {
  char buf[128];
#if defined(LONGDOUBLE_IS_IEEE754)
  std::snprintf(buf, sizeof(buf), "%.21Lg", (long double)v);   // 21 digits round-trip an 80-bit long double
  return string(buf);
#elif defined(LONGDOUBLE_IS_FLOAT128)
  quadmath_snprintf(buf, sizeof(buf), "%.36Qg", v);
  return string(buf);
#else
  return print_longdouble(v);
#endif
}

// Tsun in seconds (GM_sun/c^3). Standard pulsar-timing constant ~4.925490947 us.
static inline longdouble Tsun_sec() { return (longdouble)4.925490947e-6L; }
static inline longdouble sec_per_day()  { return (longdouble)86400.0L; }
static inline longdouble sec_per_year() { return (longdouble)(365.25L * 86400.0L); }

static inline bool safe_positive(longdouble x) { return (x > (longdouble)0.0L); }

struct Kepler {
  longdouble pb_days = 0;
  longdouble a1_sec  = 0;
  longdouble ecc     = 0;
  longdouble om_deg  = 0;
  longdouble t0_mjd  = 0;

  longdouble omdot_degyr = 0;  // OMDOT (deg/yr)
  longdouble edot = 0;         // EDOT (1/s)

  longdouble tasc_mjd = 0;
  longdouble eps1 = 0;
  longdouble eps2 = 0;
  longdouble eps1dot = 0;      // EPS1DOT, EPS2DOT (1/s)
  longdouble eps2dot = 0;

  longdouble pb_classic = 0;   // PB as DD/BT/T2 define it (periastron to periastron)
  longdouble pb_ell1    = 0;   // PB as ELL1 defines it (ascending node to ascending node)

  bool has_classic = false; // T0 (+ ECC, OM)
  bool has_ell1    = false; // TASC (+ EPS1, EPS2)
  bool epsdot_lost = false; // EPS1DOT/EPS2DOT given with EPS1 = EPS2 = 0: no OMDOT/EDOT equivalent
};

// Both parametrisations, one of them read from the par and the other converted, with the small-e relations
// of the ELL1 model (Lange et al. 2001). tempo2's DD advances w as w0 + k*(true anomaly), k = OMDOT/n, so for
// e -> 0 the ascending node moves at (1+k) n:
//   PB_ell1 = PB / (1+k)          TASC = T0 - w_asc / n,  w_asc = w0 / (1+k)   (the epoch where w = w_asc)
//   EPS1 = e sin w_asc            EPS1DOT = EDOT sin w_asc + e wdot cos w_asc
//   EPS2 = e cos w_asc            EPS2DOT = EDOT cos w_asc - e wdot sin w_asc
// and the exact inverse for ELL1 -> T0/ECC/OM. Terms of order x e^2 are neglected, as in ELL1 itself.
static Kepler get_kepler_from_psr(const pulsar &psr) {
  Kepler k;
  const longdouble pi = (longdouble)3.141592653589793238462643383279502884L;
  const longdouble degyr = rad2deg_ld((longdouble)1.0L) * sec_per_year();   // rad/s -> deg/yr

  if (param_is_set(psr, param_pb)) k.pb_days = param_val0(psr, param_pb);
  if (param_is_set(psr, param_a1)) k.a1_sec  = param_val0(psr, param_a1);

  k.has_classic = param_is_set(psr, param_t0);
  k.has_ell1    = param_is_set(psr, param_tasc);
  if (!safe_positive(k.pb_days)) return k;

  if (k.has_classic) {
    k.ecc    = param_is_set(psr, param_ecc) ? param_val0(psr, param_ecc) : (longdouble)0.0L;
    k.om_deg = param_is_set(psr, param_om)  ? param_val0(psr, param_om)  : (longdouble)0.0L;
    k.t0_mjd = param_val0(psr, param_t0);
    k.omdot_degyr = param_is_set(psr, param_omdot) ? param_val0(psr, param_omdot) : (longdouble)0.0L;
    k.edot        = param_is_set(psr, param_edot)  ? param_val0(psr, param_edot)  : (longdouble)0.0L;

    longdouble n    = (longdouble)2.0L * pi / (k.pb_days * sec_per_day());
    longdouble wdot = k.omdot_degyr / degyr;                  // rad/s
    longdouble kk   = wdot / n;
    k.pb_classic = k.pb_days;
    k.pb_ell1    = k.pb_days / ((longdouble)1.0L + kk);

    longdouble w_asc = deg2rad_ld(k.om_deg) / ((longdouble)1.0L + kk);
    k.tasc_mjd = k.t0_mjd - w_asc / n / sec_per_day();
    longdouble e_asc = k.ecc + k.edot * (k.tasc_mjd - k.t0_mjd) * sec_per_day();
    longdouble s = (longdouble)sinl((long double)w_asc), c = (longdouble)cosl((long double)w_asc);
    k.eps1    = e_asc * s;
    k.eps2    = e_asc * c;
    k.eps1dot = k.edot * s + e_asc * wdot * c;
    k.eps2dot = k.edot * c - e_asc * wdot * s;
  } else if (k.has_ell1) {
    k.eps1     = param_is_set(psr, param_eps1) ? param_val0(psr, param_eps1) : (longdouble)0.0L;
    k.eps2     = param_is_set(psr, param_eps2) ? param_val0(psr, param_eps2) : (longdouble)0.0L;
    k.eps1dot  = param_is_set(psr, param_eps1dot) ? param_val0(psr, param_eps1dot) : (longdouble)0.0L;
    k.eps2dot  = param_is_set(psr, param_eps2dot) ? param_val0(psr, param_eps2dot) : (longdouble)0.0L;
    k.tasc_mjd = param_val0(psr, param_tasc);

    longdouble e_asc = (longdouble)hypotl((long double)k.eps1, (long double)k.eps2);
    longdouble w_asc = (e_asc > (longdouble)0.0L) ? (longdouble)atan2l((long double)k.eps1, (long double)k.eps2)
                                                  : (longdouble)0.0L;
    longdouble wdot = 0;
    if (e_asc > (longdouble)0.0L) {
      wdot   = (k.eps1dot * k.eps2 - k.eps2dot * k.eps1) / (e_asc * e_asc);
      k.edot = (k.eps1 * k.eps1dot + k.eps2 * k.eps2dot) / e_asc;
    } else if (k.eps1dot != (longdouble)0.0L || k.eps2dot != (longdouble)0.0L) {
      k.epsdot_lost = true;
    }
    longdouble n_ell = (longdouble)2.0L * pi / (k.pb_days * sec_per_day());
    longdouble q  = wdot / n_ell;
    longdouble kk = q / ((longdouble)1.0L - q);               // wdot / n with n the DD mean motion
    k.pb_ell1    = k.pb_days;
    k.pb_classic = k.pb_days * ((longdouble)1.0L + kk);

    longdouble n = (longdouble)2.0L * pi / (k.pb_classic * sec_per_day());
    k.t0_mjd = k.tasc_mjd + w_asc / n / sec_per_day();
    k.ecc    = e_asc + k.edot * (k.t0_mjd - k.tasc_mjd) * sec_per_day();
    k.om_deg = wrap_0_360(rad2deg_ld(w_asc * ((longdouble)1.0L + kk)));
    k.omdot_degyr = wdot * degyr;
  }
  return k;
}

// mass function in solar masses
static bool mass_function(longdouble a1_sec, longdouble pb_days, longdouble &f_out) {
  if (!safe_positive(a1_sec) || !safe_positive(pb_days)) return false;
  const longdouble pi = (longdouble)3.141592653589793238462643383279502884L;
  longdouble pb_sec = pb_days * sec_per_day();
  longdouble n = (longdouble)2.0L * pi / pb_sec; // rad/s
  longdouble f = (n*n) * (a1_sec*a1_sec*a1_sec) / Tsun_sec();
  f_out = f;
  return safe_positive(f);
}

// GR periastron advance -> MTOT
static bool mtot_from_omdot_gr(longdouble omdot_deg_per_yr, longdouble pb_days, longdouble ecc, longdouble &mtot_out) {
  if (!safe_positive(pb_days)) return false;
  if (ecc < (longdouble)0.0L) return false;

  const longdouble pi = (longdouble)3.141592653589793238462643383279502884L;
  longdouble pb_sec = pb_days * sec_per_day();
  longdouble omdot_rad_per_s = deg2rad_ld(omdot_deg_per_yr) / sec_per_year();

  longdouble one_minus_e2 = (longdouble)1.0L - ecc*ecc;
  if (one_minus_e2 <= (longdouble)0.0L) return false;

  longdouble A = (omdot_rad_per_s * one_minus_e2) / (longdouble)3.0L;
  if (A <= (longdouble)0.0L) return false;

  longdouble B = (pb_sec / ((longdouble)2.0L*pi));              // Pb/(2pi)
  longdouble B53 = (longdouble)powl((long double)B, (long double)(5.0L/3.0L));
  longdouble Ts23 = (longdouble)powl((long double)Tsun_sec(), (long double)(2.0L/3.0L));

  longdouble inside = A * B53 / Ts23;
  if (inside <= (longdouble)0.0L) return false;

  longdouble mtot = (longdouble)powl((long double)inside, (long double)(3.0L/2.0L));
  if (!safe_positive(mtot)) return false;
  mtot_out = mtot;
  return true;
}

// DDH orthometric conversions (assume cos i >= 0 branch; Shapiro is symmetric)
static bool sini_from_stig(longdouble stig, longdouble &sini_out) {
  if (stig <= (longdouble)0.0L) return false;
  longdouble s = (longdouble)2.0L * stig / ((longdouble)1.0L + stig*stig);
  if (s <= (longdouble)0.0L || s > (longdouble)1.0L) return false;
  sini_out = s;
  return true;
}
static bool stig_from_sini(longdouble sini, longdouble &stig_out) {
  if (sini <= (longdouble)0.0L || sini > (longdouble)1.0L) return false;
  longdouble cosi = (longdouble)sqrtl((long double)((longdouble)1.0L - sini*sini)); // + branch
  longdouble v = sini / ((longdouble)1.0L + cosi);
  if (v <= (longdouble)0.0L) return false;
  stig_out = v;
  return true;
}
static bool h3_from_m2_stig(longdouble m2, longdouble stig, longdouble &h3_out) {
  if (!safe_positive(m2) || stig <= (longdouble)0.0L) return false;
  longdouble h3 = Tsun_sec() * m2 * (stig*stig*stig);
  if (!safe_positive(h3)) return false;
  h3_out = h3;
  return true;
}
static bool m2_from_h3_stig(longdouble h3, longdouble stig, longdouble &m2_out) {
  if (!safe_positive(h3) || stig <= (longdouble)0.0L) return false;
  longdouble denom = Tsun_sec() * (stig*stig*stig);
  if (!safe_positive(denom)) return false;
  longdouble m2 = h3 / denom;
  if (!safe_positive(m2)) return false;
  m2_out = m2;
  return true;
}

// Derive MTOT from (f, m2, sini)
static bool mtot_from_f_m2_sini(longdouble f, longdouble m2, longdouble sini, longdouble &mtot_out) {
  if (!safe_positive(f) || !safe_positive(m2) || sini <= (longdouble)0.0L) return false;
  longdouble s32 = (longdouble)powl((long double)sini, (long double)(3.0L/2.0L));
  longdouble m232 = (longdouble)powl((long double)m2, (long double)(3.0L/2.0L));
  longdouble mtot = (m232 * s32) / (longdouble)sqrtl((long double)f);
  if (!safe_positive(mtot)) return false;
  mtot_out = mtot;
  return true;
}

// Shapiro / mass parameters: the ones in the par file, plus whatever follows from them.
struct ShapiroState {
  std::map<string, longdouble> val;   // M2, SINI, H3, STIG, MTOT
  std::map<string, int> fit;          // fit flag of each value (inherited for derived values)
  bool has(const string &k) const { return val.count(k) > 0; }
};

static ShapiroState get_shapiro_state(const pulsar &psr, const Kepler &k) {
  ShapiroState S;
  const std::pair<const char*, int> pars[] = {
    {"M2", param_m2}, {"SINI", param_sini}, {"H3", param_h3}, {"STIG", param_stig}, {"MTOT", param_mtot}};
  for (const auto &p : pars) {
    if (param_is_set(psr, p.second)) { S.val[p.first] = param_val0(psr, p.second); S.fit[p.first] = param_fitflag0(psr, p.second); }
  }
  longdouble v;
  // H3/STIG -> M2/SINI
  if (S.has("STIG") && !S.has("SINI") && sini_from_stig(S.val["STIG"], v)) { S.val["SINI"] = v; S.fit["SINI"] = S.fit["STIG"]; }
  if (S.has("H3") && S.has("STIG") && !S.has("M2") && m2_from_h3_stig(S.val["H3"], S.val["STIG"], v)) {
    S.val["M2"] = v; S.fit["M2"] = std::max(S.fit["H3"], S.fit["STIG"]);
  }
  // M2/SINI -> H3/STIG
  if (S.has("SINI") && !S.has("STIG") && stig_from_sini(S.val["SINI"], v)) { S.val["STIG"] = v; S.fit["STIG"] = S.fit["SINI"]; }
  if (S.has("M2") && S.has("STIG") && !S.has("H3") && h3_from_m2_stig(S.val["M2"], S.val["STIG"], v)) {
    S.val["H3"] = v; S.fit["H3"] = std::max(S.fit["M2"], S.fit["STIG"]);
  }
  // MTOT: from the mass function with M2 and SINI, else from OMDOT assuming GR
  longdouble f;
  if (!S.has("MTOT") && S.has("M2") && S.has("SINI") && mass_function(k.a1_sec, k.pb_days, f) &&
      mtot_from_f_m2_sini(f, S.val["M2"], S.val["SINI"], v)) {
    S.val["MTOT"] = v; S.fit["MTOT"] = 0;
  }
  if (!S.has("MTOT") && param_is_set(psr, param_omdot) &&
      mtot_from_omdot_gr(param_val0(psr, param_omdot), k.pb_days, k.ecc, v)) {
    S.val["MTOT"] = v; S.fit["MTOT"] = 0;
  }
  return S;
}

// The line with its value (second word) replaced, keeping the fit flag, the uncertainty and the spacing.
static string replace_value(const string &line, const longdouble &val) {
  size_t a = line.find_first_not_of(" \t");
  a = line.find_first_of(" \t", a);
  a = line.find_first_not_of(" \t", a);
  if (a == string::npos) return line;
  size_t b = line.find_first_of(" \t\r\n", a);
  return line.substr(0, a) + fmt_ld(val) + (b == string::npos ? "" : line.substr(b));
}

static void write_param_line(std::ostream &os, const string &key, const longdouble &val, int fitflag) {
  os << key << " " << fmt_ld(val) << " " << (fitflag ? 1 : 0) << "\n";
}

static void usage() {
  std::fprintf(stderr,
    "Usage: tempo2 -gr binswitch -f in.par in.tim --binary <DD|DDH|DDGR|BT|T2|ELL1> --out out.par "
    "[--keep-fitflags] [--test]\n");
}

extern "C" int graphicalInterface(int argc, char *argv[], pulsar *psr, int *npsr) {
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  std::memset(parFile, 0, sizeof(parFile));
  std::memset(timFile, 0, sizeof(timFile));

  *npsr = 1;

  string outPath = "";
  string targetBinary = "";
  bool keep_fitflags = false;
  bool do_test = false;

  for (int i = 1; i < argc; i++) {
    if (std::strcmp(argv[i], "-f") == 0 && i + 2 < argc) {
      std::strncpy(parFile[0], argv[i+1], MAX_FILELEN-1);
      std::strncpy(timFile[0], argv[i+2], MAX_FILELEN-1);
      i += 2;
    } else if ((std::strcmp(argv[i], "--out") == 0 || std::strcmp(argv[i], "-o") == 0) && i + 1 < argc) {
      outPath = argv[++i];
    } else if ((std::strcmp(argv[i], "--binary") == 0 || std::strcmp(argv[i], "-b") == 0) && i + 1 < argc) {
      targetBinary = argv[++i];
    } else if (std::strcmp(argv[i], "--keep-fitflags") == 0) {
      keep_fitflags = true;
    } else if (std::strcmp(argv[i], "--test") == 0) {
      do_test = true;
    }
  }

  if (parFile[0][0] == '\0' || timFile[0][0] == '\0') {
    std::fprintf(stderr, "binswitch: ERROR: need -f <par> <tim>\n");
    usage();
    return 1;
  }

  // Values exactly as in the par file: no preProcess, no fit.
  readParfile(psr, parFile, timFile, *npsr);
  readTimfile(psr, timFile, *npsr);

  pulsar &P = psr[0];

  string currentBinary = P.binaryModel;
  std::transform(targetBinary.begin(), targetBinary.end(), targetBinary.begin(), [](unsigned char c){ return std::toupper(c); });
  if (targetBinary.empty()) targetBinary = currentBinary;
  if (targetBinary.empty()) targetBinary = "DD";
  const std::set<string> supported = {"DD", "DDH", "DDGR", "BT", "T2", "ELL1"};
  if (!supported.count(targetBinary)) {
    std::fprintf(stderr, "binswitch: ERROR: unsupported target model '%s'\n", targetBinary.c_str());
    usage();
    return 1;
  }

  if (outPath.empty()) outPath = string(P.name) + ".base.par";

  std::ifstream inpar(parFile[0]);
  if (!inpar) {
    std::fprintf(stderr, "binswitch: ERROR: cannot open input par file: %s\n", parFile[0]);
    return 2;
  }
  std::vector<string> lines;
  { string line; while (std::getline(inpar, line)) lines.push_back(line); }

  Kepler k = get_kepler_from_psr(P);
  if (!(k.has_classic || k.has_ell1)) {
    std::fprintf(stderr, "binswitch: ERROR: the par file has neither T0 nor TASC: no binary orbit to convert\n");
    return 2;
  }
  if (!param_is_set(P, param_pb) || !param_is_set(P, param_a1)) {
    std::fprintf(stderr, "binswitch: ERROR: the par file needs PB and A1\n");
    return 2;
  }
  ShapiroState S = get_shapiro_state(P, k);

  const bool to_ell1   = (targetBinary == "ELL1");
  const bool from_ell1 = !k.has_classic;
  const bool switch_family = (to_ell1 != from_ell1);   // T0/ECC/OM <-> TASC/EPS1/EPS2

  // Parameters that depend on the model. The orbital ones are kept verbatim if the parametrisation does not
  // change and are all rewritten if it does. The Shapiro/mass ones are kept verbatim if the target uses them.
  // Whatever the target needs and the output does not have yet is written, converted, at the end.
  std::vector<string> orbital = {"T0", "ECC", "OM", "OMDOT", "EDOT", "TASC", "EPS1", "EPS2", "EPS1DOT", "EPS2DOT"};
  // Rescaled in place (the rest of their line is kept): PB, which ELL1 counts between ascending nodes, and F0.
  // The DD delay contains -(3/2) x EPS1, a constant that ELL1 leaves out; when EPS1 changes with time
  // (EPS1DOT, A1DOT) it is a drift, i.e. a change of F0 by the factor (1 - dC/dt), dC/dt = -(3/2) d(x EPS1)/dt.
  std::map<string, longdouble> inplace;
  if (switch_family && k.pb_ell1 != k.pb_classic) inplace["PB"] = to_ell1 ? k.pb_ell1 : k.pb_classic;
  {
    longdouble a1dot = param_is_set(P, param_a1dot) ? param_val0(P, param_a1dot) : (longdouble)0.0L;
    longdouble dcdt = (longdouble)-1.5L * (k.a1_sec * k.eps1dot + a1dot * k.eps1);
    longdouble f0 = param_val0(P, param_f);
    if (switch_family && dcdt != (longdouble)0.0L && param_is_set(P, param_f))
      inplace["F0"] = to_ell1 ? f0 * ((longdouble)1.0L - dcdt) : f0 / ((longdouble)1.0L - dcdt);
  }
  const std::vector<string> shapiro = {"M2", "SINI", "H3", "STIG", "MTOT"};
  std::vector<string> needed = to_ell1 ? std::vector<string>{"TASC", "EPS1", "EPS2"}
                                       : std::vector<string>{"T0", "ECC", "OM"};
  if (switch_family) {
    if (to_ell1 && (k.eps1dot != (longdouble)0.0L || k.eps2dot != (longdouble)0.0L)) {
      needed.push_back("EPS1DOT"); needed.push_back("EPS2DOT");
    }
    if (!to_ell1 && k.omdot_degyr != (longdouble)0.0L) needed.push_back("OMDOT");
    if (!to_ell1 && k.edot != (longdouble)0.0L) needed.push_back("EDOT");
  } else if (to_ell1) {
    if (param_is_set(P, param_eps1dot)) needed.push_back("EPS1DOT");
    if (param_is_set(P, param_eps2dot)) needed.push_back("EPS2DOT");
  } else {
    if (param_is_set(P, param_omdot)) needed.push_back("OMDOT");
    if (param_is_set(P, param_edot)) needed.push_back("EDOT");
  }
  std::vector<string> shap;
  if (targetBinary == "DDH")       shap = {"H3", "STIG"};
  else if (targetBinary == "DDGR") shap = {"M2", "MTOT"};
  else                             shap = {"M2", "SINI"};
  for (const auto &s : shap) if (S.has(s)) needed.push_back(s);

  // Errors and warnings about what cannot be converted
  if (targetBinary == "DDGR" && !(S.has("M2") && S.has("MTOT"))) {
    std::fprintf(stderr, "binswitch: ERROR: DDGR needs M2 and MTOT, which cannot be found or derived from this par "
                         "(give M2 and SINI, H3 and STIG, or OMDOT). Nothing written.\n");
    return 4;
  }
  if (targetBinary == "DDH" && !(S.has("H3") && S.has("STIG"))) {
    std::fprintf(stdout, "binswitch: WARNING: no H3/STIG (or M2/SINI) in the par: the DDH par has no Shapiro delay.\n");
  }
  if (param_is_set(P, param_h3) && !param_is_set(P, param_stig) && targetBinary != "DDH") {
    std::fprintf(stdout, "binswitch: WARNING: H3 without STIG cannot be converted to M2/SINI: H3 is dropped.\n");
  }
  if (param_is_set(P, param_h4)) {
    std::fprintf(stdout, "binswitch: WARNING: H4 is not converted (kept as it is).\n");
  }
  if (param_is_set(P, param_kin)) {
    std::fprintf(stdout, "binswitch: WARNING: KIN is not converted (kept as it is).\n");
  }
  if (param_is_set(P, param_shapmax)) {
    std::fprintf(stdout, "binswitch: WARNING: SHAPMAX is not converted (kept as it is).\n");
  }
  if (switch_family && k.epsdot_lost) {
    std::fprintf(stdout, "binswitch: WARNING: EPS1DOT/EPS2DOT with EPS1 = EPS2 = 0 have no OMDOT/EDOT equivalent: "
                         "they are dropped.\n");
  }
  if (currentBinary == "DDGR" && targetBinary != "DDGR") {
    std::fprintf(stdout, "binswitch: NOTE: DDGR derives the post-Keplerian parameters from GR; they are not written. "
                         "Fit OMDOT/GAMMA/PBDOT/SINI in the new model.\n");
  }

  // Accuracy of ELL1: it neglects terms of order x e^2 in the Roemer delay and x^2 n e in tempo2's inverse
  // timing correction (checked against DD on simulated TOAs: the difference scales as x^2 n e).
  if (switch_family) {
    const longdouble pi = (longdouble)3.141592653589793238462643383279502884L;
    longdouble n = (longdouble)2.0L * pi / (k.pb_classic * sec_per_day());
    longdouble xe2_us  = k.a1_sec * k.ecc * k.ecc * (longdouble)1e6L;
    longdouble x2ne_us = k.a1_sec * k.a1_sec * n * k.ecc * (longdouble)1e6L;
    std::vector<double> err;
    for (int i = 0; i < P.nobs; i++) if (!P.obsn[i].deleted) err.push_back((double)P.obsn[i].toaErr);
    double med = 0;
    if (!err.empty()) { std::nth_element(err.begin(), err.begin() + err.size()/2, err.end()); med = err[err.size()/2]; }
    double limit = err.empty() ? 0 : med / std::sqrt((double)err.size());
    double ell1_err = (double)(xe2_us + x2ne_us);
    if (do_test) {
      std::fprintf(stdout, "binswitch: test: ECC = %.3Le; ELL1 error scale x*e^2 = %.3Le us, x^2*n*e = %.3Le us; "
                           "median TOA error %.3g us over %zu TOAs (median/sqrt(N) = %.3g us)\n",
                   (long double)k.ecc, (long double)xe2_us, (long double)x2ne_us, med, err.size(), limit);
    }
    if (!err.empty() && ell1_err > limit) {
      std::fprintf(stdout, "binswitch: WARNING: the ELL1 approximation error (x*e^2 + x^2*n*e = %.3g us) is not small "
                           "compared with the timing precision (%.3g us): ELL1 is not accurate for this orbit.\n",
                   ell1_err, limit);
    }
  }

  // Values and inherited fit flags of everything that may need to be written
  auto ff = [&](int ip) { return keep_fitflags ? param_fitflag0(P, ip) : 0; };
  std::map<string, longdouble> value = {
    {"T0", k.t0_mjd}, {"ECC", k.ecc}, {"OM", k.om_deg}, {"OMDOT", k.omdot_degyr}, {"EDOT", k.edot},
    {"TASC", k.tasc_mjd}, {"EPS1", k.eps1}, {"EPS2", k.eps2}, {"EPS1DOT", k.eps1dot}, {"EPS2DOT", k.eps2dot}};
  const int ff_eps = std::max(ff(param_eps1), ff(param_eps2)), ff_eom = std::max(ff(param_ecc), ff(param_om));
  const int ff_epsdot = std::max(ff(param_eps1dot), ff(param_eps2dot));
  const int ff_omedot = std::max(ff(param_omdot), ff(param_edot));
  std::map<string, int> flag = {
    {"T0", ff(param_tasc)}, {"TASC", ff(param_t0)}, {"ECC", ff_eps}, {"OM", ff_eps}, {"EPS1", ff_eom}, {"EPS2", ff_eom},
    {"OMDOT", ff_epsdot}, {"EDOT", ff_epsdot}, {"EPS1DOT", ff_omedot}, {"EPS2DOT", ff_omedot}};
  for (const auto &kv : S.val) { value[kv.first] = kv.second; flag[kv.first] = keep_fitflags ? S.fit[kv.first] : 0; }

  // Output: input lines with the BINARY line replaced and unused model parameters removed, then the converted ones
  std::set<string> present;
  std::vector<string> kept;
  bool wrote_binary = false;
  for (const auto &ln : lines) {
    string tok = first_token(ln);
    if (tok == "BINARY") {
      if (!wrote_binary) kept.push_back("BINARY " + targetBinary);
      wrote_binary = true;
      continue;
    }
    if (inplace.count(tok) && !present.count(tok)) {
      kept.push_back(replace_value(ln, inplace[tok]));
      present.insert(tok);
      continue;
    }
    bool is_orbital = std::find(orbital.begin(), orbital.end(), tok) != orbital.end();
    bool is_shapiro = std::find(shapiro.begin(), shapiro.end(), tok) != shapiro.end();
    if (is_orbital || is_shapiro) {
      bool keep = std::find(needed.begin(), needed.end(), tok) != needed.end() && !present.count(tok) &&
                  !(is_orbital && switch_family);
      if (!keep) continue;
      present.insert(tok);
    }
    kept.push_back(ln);
  }
  if (!wrote_binary) kept.push_back("BINARY " + targetBinary);

  std::ofstream out(outPath.c_str());
  if (!out) {
    std::fprintf(stderr, "binswitch: ERROR: cannot write output par file: %s\n", outPath.c_str());
    return 3;
  }
  for (const auto &ln : kept) out << ln << "\n";

  std::vector<string> converted;
  for (const auto &key : needed) if (!present.count(key)) converted.push_back(key);
  if (!converted.empty()) {
    out << "# --- converted by binswitch from BINARY " << (currentBinary.empty() ? "?" : currentBinary) << " ---\n";
    for (const auto &key : converted) write_param_line(out, key, value[key], flag[key]);
  }

  std::fprintf(stdout, "binswitch: wrote %s (BINARY %s -> %s; converted:", outPath.c_str(),
               currentBinary.c_str(), targetBinary.c_str());
  for (const auto &key : converted) std::fprintf(stdout, " %s", key.c_str());
  for (const auto &kv : inplace) std::fprintf(stdout, " %s", kv.first.c_str());
  std::fprintf(stdout, "%s)\n", (converted.empty() && inplace.empty()) ? " nothing" : "");
  return 0;
}
