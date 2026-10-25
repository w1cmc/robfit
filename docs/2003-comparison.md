# 2003 srobunix vs current robfit — Fortran comparison

*Compared 2026-10-25: `~/src/robfit.2003` (132 Fortran files, 239 distinct routines) against this repo at 5302ad0 (branch robfitdeb; 804 Fortran files, ~1,080 distinct routines).*

Method: each SUBROUTINE/FUNCTION/PROGRAM was extracted and normalized (case, whitespace, comments, continuation lines removed) before hashing, so "identical" means same code, not just same name.

## Summary

| Result | Count |
|---|---|
| Identical in both | 18 |
| Same name, different body | 34 |
| Only in 2003 copy | 187 |
| Only in current repo | ~840 |
| Renamed but identical | 0 |

**Identical:** KIJ, BINDEX, CCHAR, CHANGE, CLBS, MAOPEN, GETARG, NARGS, RSEED, TIME, TLOC, SLEEP, plus DOS/`*QQ` shims.

**Modified:** core minimizer/linear algebra (ROBMIN, SMINV, GAUSSJ, GSOLVE, BRACK, NRAF, SMVALS, BNRAPH), MENURD/CTON, BLI, DCF, EFFFIT, G3LINK, RUNKPS, RESL, SORT, NLFIT.

## Only in the 2003 copy

1. **Old "book" robfit programs** — VRMAIN, RAWDD, STGEN, STDIS, FSPFIT, FSPDIS, XCALIBER, plus link routines (BKLINK, G4LINK, GSLINK, SHLINK, VRLINK), SMSQ, and the MACBKG background fit. Superseded by the current programs.
2. **Unkfit/Unkdis isotope-identification engine** — UNKFIT (1,250 lines), UNKLIN, DCF2, NLMASS, UNKFI, STANBK, DATAC, FIDDLE, DISMAIN. No equivalent in the current tree.
3. **RobWin DLL display code** (`fspdidll/`: FSPDIG, STDISG, RESOUT).
4. **Graphics/platform code** — tekplot, watcplot, `plots/`, DOS shims. Obsolete.

## Notes relevant to open TODOs

- **NLERREST / KIJ crash:** KIJ is identical in both trees (`KIJ = J*(J-1)/2 + I`, packed symmetric index). The fault is in the caller (zero/uninitialized index or undersized array), not KIJ.
- **EffCal missing symbols:** `escapep`, `funfcap`, `funfcap3` are not in the 2003 copy but are defined in this repo: `robfitdocs/templfit/funwp.for`, `robfitdocs/templfit/Efficiency/Escape/ESCAPE.FOR`, `escape (3).for`, and `Efficiency/EgyToEff.for`.
- **Duplicate routines:** the current repo often has several versions of one routine (e.g. ROBMIN in `Robmin/`, `Robmin/Orig/`, and several `templfit/attenuation/nlfit/*/` folders). Choosing canonical versions matters more than the 2003 deltas.

## Reference data in the 2003 copy

- `srobunix/tcases/ztcase1–10.sp` and old output files
- `unktest/ztcase9.sp` + `ZTCASE9.OUT`, `U238KLM.sp` + `U238KLM.OUT`, `DETCBA~1.OUT`, isotope libraries `*.IIS`
