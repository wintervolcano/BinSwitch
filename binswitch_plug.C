// binswitch_plug.C
//
// TEMPO2 graphical plugin:
//   - runs a standard TEMPO2 fit (if doFitAll is available; resolved via dlsym to avoid ABI issues)
//   - rewrites a "base" par file preserving all original lines except binary params
//   - converts between DD, DDH, T2, BT, ELL1, DDGR
//   - uses tempo2 longdouble + print_longdouble() for precision-safe formatting
//
// Build/install: see README.md.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <cctype>

#include <dlfcn.h>

#include "tempo2.h"

using std::string;

static inline string ltrim(const string &s) {
  size_t i = 0;
  while (i < s.size() && std::isspace((unsigned char)s[i])) i++;
  return s.substr(i);
}
static inline string rtrim(const string &s) {
  if (s.empty()) return s;
  size_t i = s.size();
  while (i > 0 && std::isspace((unsigned char)s[i-1])) i--;
  return s.substr(0, i);
}
static inline string trim(const string &s) { return rtrim(ltrim(s)); }

static inline string first_token(const string &line) {
  string t = ltrim(line);
  if (t.empty() || t[0] == '#') return "";
  size_t sp = t.find_first_of(" \t\r\n");
  return (sp == string::npos) ? t : t.substr(0, sp);
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
  return (p.param[ip].fitFlag ? p.param[ip].fitFlag[0] : 0);
}

// ---- runtime resolution of doFitAll (avoid ABI/mangling differences) ----
using doFitAll_fn = void (*)(pulsar*, int, int);

static doFitAll_fn resolve_doFitAll() {
  // Try unmangled (your lib exports: T doFitAll)
  void *p = dlsym(RTLD_DEFAULT, "doFitAll");
  if (p) return reinterpret_cast<doFitAll_fn>(p);

  // Try the mangled symbol (some builds export this)
  p = dlsym(RTLD_DEFAULT, "_Z8doFitAllP6pulsarii");
  if (p) return reinterpret_cast<doFitAll_fn>(p);

  return nullptr;
}

// Tsun in seconds (GM_sun/c^3). Standard pulsar-timing constant ~4.925490947 us.
static inline longdouble Tsun_sec() { return (longdouble)4.925490947e-6L; }
static inline longdouble sec_per_day()  { return (longdouble)86400.0L; }
static inline longdouble sec_per_year() { return (longdouble)(365.25L * 86400.0L); }

struct Kepler {
  longdouble pb_days = 0;
  longdouble a1_sec  = 0;
  longdouble ecc     = 0;
  longdouble om_deg  = 0;
  longdouble t0_mjd  = 0;

  longdouble tasc_mjd = 0;
  longdouble eps1 = 0;
  longdouble eps2 = 0;

  bool has_classic = false; // ECC/OM/T0
  bool has_ell1    = false; // EPS1/EPS2/TASC
};

static Kepler get_kepler_from_psr(const pulsar &psr) {
  Kepler k;
  if (param_is_set(psr, param_pb)) k.pb_days = param_val0(psr, param_pb);
  if (param_is_set(psr, param_a1)) k.a1_sec  = param_val0(psr, param_a1);

  bool has_classic = param_is_set(psr, param_ecc) && param_is_set(psr, param_om) && param_is_set(psr, param_t0);
  bool has_ell1    = param_is_set(psr, param_eps1) && param_is_set(psr, param_eps2) && param_is_set(psr, param_tasc);

  k.has_classic = has_classic;
  k.has_ell1    = has_ell1;

  if (has_classic) {
    k.ecc    = param_val0(psr, param_ecc);
    k.om_deg = param_val0(psr, param_om);
    k.t0_mjd = param_val0(psr, param_t0);

    // Approximate ELL1 conversion (valid for small-e); keep but warn elsewhere if needed.
    longdouble om_rad = deg2rad_ld(k.om_deg);
    k.eps1 = k.ecc * (longdouble)sin((long double)om_rad);
    k.eps2 = k.ecc * (longdouble)cos((long double)om_rad);
    k.tasc_mjd = k.t0_mjd - k.pb_days * (k.om_deg / (longdouble)360.0L);
  } else if (has_ell1) {
    k.eps1    = param_val0(psr, param_eps1);
    k.eps2    = param_val0(psr, param_eps2);
    k.tasc_mjd = param_val0(psr, param_tasc);

    longdouble ecc = (longdouble)hypotl((long double)k.eps1, (long double)k.eps2);
    k.ecc = ecc;

    longdouble om_rad = (ecc > (longdouble)0.0L) ? (longdouble)atan2l((long double)k.eps1, (long double)k.eps2)
                                                 : (longdouble)0.0L;
    k.om_deg = wrap_0_360(rad2deg_ld(om_rad));
    k.t0_mjd = k.tasc_mjd + k.pb_days * (k.om_deg / (longdouble)360.0L);
  }
  return k;
}

static inline bool safe_positive(longdouble x) { return (x > (longdouble)0.0L); }

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

struct ShapiroState {
  bool have_m2 = false;
  bool have_sini = false;
  bool have_mtot = false;
  bool have_h3 = false;
  bool have_stig = false;

  longdouble m2 = 0;
  longdouble sini = 0;
  longdouble mtot = 0;
  longdouble h3 = 0;
  longdouble stig = 0;

  bool warned = false;
};

static ShapiroState get_shapiro_state(const pulsar &psr) {
  ShapiroState s;
  if (param_is_set(psr, param_m2))   { s.have_m2 = true;   s.m2   = param_val0(psr, param_m2); }
  if (param_is_set(psr, param_sini)) { s.have_sini = true; s.sini = param_val0(psr, param_sini); }
  if (param_is_set(psr, param_mtot)) { s.have_mtot = true; s.mtot = param_val0(psr, param_mtot); }
  // We do NOT rely on param_h3/param_stig enums (often missing). We'll derive H3/STIG from M2/SINI when needed.
  return s;
}

static void infer_masses_and_shapiro(const pulsar &psr, const Kepler &k, ShapiroState &S) {
  longdouble f = 0;
  bool have_f = mass_function(k.a1_sec, k.pb_days, f);

  // If we can derive STIG/H3 from M2/SINI, do it.
  if ((!S.have_stig || !S.have_h3) && S.have_m2 && S.have_sini) {
    longdouble stig;
    if (!S.have_stig && stig_from_sini(S.sini, stig)) { S.stig = stig; S.have_stig = true; }
    longdouble h3;
    if (S.have_stig && !S.have_h3 && h3_from_m2_stig(S.m2, S.stig, h3)) { S.h3 = h3; S.have_h3 = true; }
  }

  // If MTOT missing and we have M2+SINI+f -> infer MTOT
  if (!S.have_mtot && S.have_m2 && S.have_sini && have_f) {
    longdouble mtot;
    if (mtot_from_f_m2_sini(f, S.m2, S.sini, mtot)) { S.mtot = mtot; S.have_mtot = true; }
  }

  // If have OMDOT -> infer MTOT (GR)
  if (!S.have_mtot && param_is_set(psr, param_omdot) && k.has_classic) {
    longdouble mtot;
    if (mtot_from_omdot_gr(param_val0(psr, param_omdot), k.pb_days, k.ecc, mtot)) {
      S.mtot = mtot;
      S.have_mtot = true;
    }
  }

  // If still missing core Shapiro bits, set sensible defaults and warn
  if (!S.have_m2 || !S.have_sini) {
    std::fprintf(stdout,
      "binswitch: WARNING: cannot infer complete Shapiro/mass set from current fit; using defaults M2=1.0, SINI=0.9\n");
    S.warned = true;
    if (!S.have_m2)   { S.m2 = (longdouble)1.0L; S.have_m2 = true; }
    if (!S.have_sini) { S.sini = (longdouble)0.9L; S.have_sini = true; }
  }

  // Ensure STIG/H3 if targeting DDH output
  if ((!S.have_stig || !S.have_h3) && S.have_m2 && S.have_sini) {
    longdouble stig;
    if (!S.have_stig && stig_from_sini(S.sini, stig)) { S.stig = stig; S.have_stig = true; }
    longdouble h3;
    if (S.have_stig && !S.have_h3 && h3_from_m2_stig(S.m2, S.stig, h3)) { S.h3 = h3; S.have_h3 = true; }
  }
}

static void maybe_print_test_conversions(const Kepler &k) {
  if (!(k.has_classic || k.has_ell1)) return;

  longdouble ecc2 = (longdouble)hypotl((long double)k.eps1, (long double)k.eps2);
  longdouble om_rad2 = (ecc2 > (longdouble)0.0L) ? (longdouble)atan2l((long double)k.eps1, (long double)k.eps2)
                                                 : (longdouble)0.0L;
  longdouble om_deg2 = wrap_0_360(rad2deg_ld(om_rad2));
  longdouble t0_2 = k.tasc_mjd + k.pb_days * (om_deg2 / (longdouble)360.0L);

  longdouble om_rad3 = deg2rad_ld(om_deg2);
  longdouble eps1_3 = ecc2 * (longdouble)sin((long double)om_rad3);
  longdouble eps2_3 = ecc2 * (longdouble)cos((long double)om_rad3);
  longdouble tasc_3 = t0_2 - k.pb_days * (om_deg2 / (longdouble)360.0L);

  auto maxabs = [](longdouble a, longdouble b){ longdouble d=a-b; return (d<0?-d:d); };

  std::fprintf(stdout,
    "binswitch: conversion test (ELL1<->classic approx): |dECC|=%g |dOMdeg|=%g |dT0|=%g(d) |dTASC|=%g(d) |dEPS1|=%g |dEPS2|=%g\n",
    (double)maxabs(k.ecc, ecc2),
    (double)maxabs(k.om_deg, om_deg2),
    (double)maxabs(k.t0_mjd, t0_2),
    (double)maxabs(k.tasc_mjd, tasc_3),
    (double)maxabs(k.eps1, eps1_3),
    (double)maxabs(k.eps2, eps2_3)
  );
}

static void write_param_line(std::ostream &os, const string &key, const longdouble &val, int fitflag) {
  os << key << " " << print_longdouble(val);
  if (fitflag == 0 || fitflag == 1) os << " " << fitflag;
  os << "\n";
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
    std::fprintf(stderr,
      "binswitch: ERROR: need -f <par> <tim>\n"
      "Usage: tempo2 -gr binswitch -f in.par in.tim --binary DDH --out out.par [--keep-fitflags] [--test]\n");
    return 1;
  }

  readParfile(psr, parFile, timFile, *npsr);
  readTimfile(psr, timFile, *npsr);
  preProcess(psr, *npsr, argc, argv);

  // Two-pass evaluate; do one fit if available.
  doFitAll_fn fitter = resolve_doFitAll();
  if (!fitter) {
    std::fprintf(stdout,
      "binswitch: WARNING: could not resolve doFitAll via dlsym; will NOT refit, will write from current loaded values.\n");
  }

  for (int it = 0; it < 2; it++) {
    formBatsAll(psr, *npsr);
    formResiduals(psr, *npsr, 1);
    if (it == 0 && fitter) fitter(psr, *npsr, 0);
  }

  pulsar &P = psr[0];

  string currentBinary = (P.binaryModel ? string(P.binaryModel) : "");
  if (targetBinary.empty()) targetBinary = currentBinary;
  if (targetBinary.empty()) targetBinary = "DD";

  if (outPath.empty()) outPath = string(P.name) + ".base.par";

  std::ifstream inpar(parFile[0]);
  if (!inpar) {
    std::fprintf(stderr, "binswitch: ERROR: cannot open input par file: %s\n", parFile[0]);
    return 2;
  }

  std::vector<string> lines;
  { string line; while (std::getline(inpar, line)) lines.push_back(line); }

  Kepler k = get_kepler_from_psr(P);
  ShapiroState S = get_shapiro_state(P);

  infer_masses_and_shapiro(P, k, S);

  if (do_test) {
    if (targetBinary == "ELL1" && k.has_classic && k.ecc > (longdouble)0.01L) {
      std::fprintf(stdout,
        "binswitch: WARNING: requested ELL1 but ECC=%g is not 'small-e'; ELL1 conversion is approximate and may be inappropriate.\n",
        (double)k.ecc);
    }
    maybe_print_test_conversions(k);
  }

  std::set<string> drop = {
    "BINARY",
    "PB","A1","A1DOT","XDOT","ECC","E","EDOT","ECCDOT","OM","OMDOT","T0",
    "TASC","EPS1","EPS2","EPS1DOT","EPS2DOT",
    "PBDOT","XPBDOT","GAMMA","DR","DTH","DTHETA",
    "M2","MTOT","SINI","SHAPMAX",
    "H3","STIG"
  };

  std::vector<string> kept;
  kept.reserve(lines.size());
  for (const auto &ln : lines) {
    string tok = first_token(ln);
    if (!tok.empty() && drop.find(tok) != drop.end()) continue;
    kept.push_back(ln);
  }

  std::ofstream out(outPath.c_str());
  if (!out) {
    std::fprintf(stderr, "binswitch: ERROR: cannot write output par file: %s\n", outPath.c_str());
    return 3;
  }

  for (const auto &ln : kept) out << ln << "\n";

  out << "\n# --- Binary block rewritten by binswitch (from post-fit values where available) ---\n";
  out << "BINARY " << targetBinary << "\n";

  int ff_a1 = keep_fitflags ? param_fitflag0(P, param_a1) : 0;
  int ff_pb = keep_fitflags ? param_fitflag0(P, param_pb) : 0;
  write_param_line(out, "A1", k.a1_sec, ff_a1);
  write_param_line(out, "PB", k.pb_days, ff_pb);

  if (targetBinary == "ELL1") {
    write_param_line(out, "TASC", k.tasc_mjd, 0);
    write_param_line(out, "EPS1", k.eps1, 0);
    write_param_line(out, "EPS2", k.eps2, 0);
  } else {
    write_param_line(out, "T0",  k.t0_mjd, 0);
    write_param_line(out, "ECC", k.ecc,    0);
    write_param_line(out, "OM",  k.om_deg, 0);
  }

  auto maybe_write_existing = [&](const char *key, int ip) {
    if (param_is_set(P, ip)) {
      int ff = keep_fitflags ? param_fitflag0(P, ip) : 0;
      write_param_line(out, key, param_val0(P, ip), ff);
    }
  };

  maybe_write_existing("OMDOT", param_omdot);
  maybe_write_existing("PBDOT", param_pbdot);
  maybe_write_existing("A1DOT", param_a1dot);
  maybe_write_existing("EDOT",  param_edot);
  maybe_write_existing("GAMMA", param_gamma);

  if (targetBinary == "DDH") {
    if (S.have_h3)   write_param_line(out, "H3",   S.h3,   0);
    if (S.have_stig) write_param_line(out, "STIG", S.stig, 0);
    if (!(S.have_h3 && S.have_stig)) {
      std::fprintf(stdout, "binswitch: WARNING: could not produce both H3 and STIG; wrote what was available.\n");
    }
  } else if (targetBinary == "DDGR") {
    if (S.have_m2)   write_param_line(out, "M2",   S.m2,   0);
    if (S.have_mtot) write_param_line(out, "MTOT", S.mtot, 0);
    if (S.have_sini) write_param_line(out, "SINI", S.sini, 0);
    if (!(S.have_m2 && S.have_mtot)) {
      std::fprintf(stdout, "binswitch: WARNING: DDGR requested but could not infer both M2 and MTOT robustly; wrote defaults where needed.\n");
    }
  } else {
    if (S.have_sini) write_param_line(out, "SINI", S.sini, 0);
    if (S.have_m2)   write_param_line(out, "M2",   S.m2,   0);
    if (S.warned) {
      std::fprintf(stdout, "binswitch: NOTE: wrote M2/SINI defaults; you may prefer a DDH fit for H3/STIG first.\n");
    }
  }

  out << "# --- End rewritten binary block ---\n";
  std::fprintf(stdout, "binswitch: wrote %s (target BINARY=%s)\n", outPath.c_str(), targetBinary.c_str());
  return 0;
}
