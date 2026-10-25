# Regression test: ztcase9 (2003 reference vs current build)

*Run 2026-10-25 on Linux x86_64, gfortran 14.2 (xPack), repo at 5302ad0 plus the fixes below.*

## Reference data (from `~/src/robfit.2003`)

- Spectrum: `unktest/ztcase9.sp` (16013 ch, 3-day background, 120% HPGe, St. Petersburg FL, June 1994)
- Energy cal (`ztcase9.cal`): E = 17.665 + 0.146277·ch keV; width cal 1.8481, 0.0105018
- Continuum: `srobunix/tcases/ztcase9.cn` (old robfit bkgfit; `exp(cubic + Σ cᵢ·(kᵢ−x)₊³)`, x = channel, see `macbkg.for` POLY)
- Peaks: `srobunix/tcases/ztcase9.pk` (291 peaks, old fspfit)
- `unktest/ZTCASE9.OUT` is Unkfit isotope-ID output; there is no current equivalent, so it is not comparable.

## Procedure

1. Inject the 2003 calibration into the .sp header (line 4: `Energy calibration: 17.665 0.146277 0 0 0 0`).
2. `robimport` (creates COM.UNF + dat.unf from `[spectrum] file`). Originally done with the legacy `robfit` and a restored `robfit.mnu`; rerun with `robimport` gives an identical continuum.
3. `bkgfit ZTCASE9.SP` with `robfit.ini` continuum range 25–2350 keV.
4. Compare `ZTCASE9CnmFit.unf` with the 2003 continuum channel by channel.

## Result: continuum PASSES

| chn | keV | data | 2003 | 2026 | ratio |
|---:|---:|---:|---:|---:|---:|
| 300 | 61.5 | 35766 | 34644 | 35096 | 1.013 |
| 1000 | 163.9 | 26134 | 26284 | 26011 | 0.990 |
| 2000 | 310.2 | 7132 | 7269 | 7148 | 0.983 |
| 4000 | 602.8 | 2080 | 1977 | 2155 | 1.090 |
| 6000 | 895.3 | 1083 | 1056 | 1060 | 1.004 |
| 10000 | 1480.4 | 312 | 333 | 305 | 0.916 |
| 15000 | 2211.8 | 102 | 109 | 106 | 0.976 |

Over ch 200–15900: median ratio 1.006, 5th–95th percentile 0.92–1.29. The high tail is in dense-peak regions (e.g. ch 1500), where the 2003 continuum was fitted alongside peaks and the 2026 one (fit mode 2, no peaks) rides higher.

## Result: peak fit runs end to end; does not yet reproduce 2003

*Updated on branch `peakfit-bridge`, after the bridge fixes there.*

`peakfit_fit_all` (the PFIT loop) over ch 300–15900, starting from no
peaks, fit type 2 (data − continuum), 2003 resolution calibration set
with `peakfit_set_calibration`, repeated (PFIT adds ≤ 50 peaks per call):

| call | peaks | χ²/channel |
|---:|---:|---:|
| 1 | 50 | 5.14 |
| 2 | 100 | 3.09 |
| 3 | 150 | 2.97 |
| 4 | 200 | 3.40 |
| 5 | 250 | 2.83 |

Against the 2003 list (277 peaks in that range), after 5 calls:

- 57 of the 58 strong 2003 peaks (area/err > 20) are found; location
  differences median +0.34 ch, 90% within 1.2 ch.
- Only 151 of 277 are matched overall, and 118 of the 2026 peaks have
  no 2003 counterpart.
- Strong single lines are split in two: 609 keV (ch 4044) is fitted
  as 4042.3 + 4049.4. The pair's area (370k) agrees with the net counts
  (370–400k) and with 2003 (415k), but each peak alone does not.
  Fitted widths are about 0.35× the 2003 widths.
- The finder re-finds the same residual spot on later passes and PFIT
  adds it again: stacks of zero-height peaks at one location
  (e.g. five at 4260.96 in a 3800–4300 test). χ² can rise as peaks are
  added (call 4).
- Skew is fixed at 1.4 by default; with `peak_vary_skew = yes` it runs
  to the 2.0 clamp and the 609 keV line is still split.

Likely causes, to decide on before changing the physics code:
1. Width definition: 2003 `W` ≈ FWHM (14–15 ch at 609 keV, as measured);
   the 2026 fit's `Wth` converges near FWHM/3, and new peaks are seeded
   from `ChnToWth`, so the resolution calibration and the model's width
   parameter may not mean the same thing.
2. PFIT has commented-out code that rejected a new peak within 1 channel
   of an existing one; RESL21A only separates new finds from each other.
3. There is no cleanup of zero-height peaks between passes.

## Bugs fixed while running this (uncommitted)

1. `main/FromSp.for` — dat.unf was written one value per record as REAL*8. Every reader expects one record of REAL*4 (the `R4TOUNF` format). Now writes one REAL*4 record.
2. `main/FromSp.for` — never set `IBEGC/IENDC` in `/INTVALS/`, so .sp spectra had an empty data range (0..0). Now sets 1..N.
3. `main/FromSp.for` — wrote `…DAT.UNF` but readers ask for `…dat.unf` (fine on macOS, breaks on case-sensitive filesystems).
4. `CnmFit/CnmFuns.for` `ChnToCnmS` — read CnmValues.unf FI values as REAL*8, but `CnmBli` writes them as REAL*4. Now reads REAL*4.

## Remaining issues found

- ~~.sp/.n42 bootstrap needed robfit.mnu + legacy roblink~~ — fixed: `robimport`, run automatically by `robfit_menu`.
- ~~FromN42 range / `DAT.UNF` name~~ — fixed (its data was already REAL*4).
- ~~Mixed `DAT.UNF` / `dat.unf` names~~ — fixed in the built programs (xmlplot, g3linkt_ini, pklink_ini, peakfit_bridge, spout); unbuilt variants still say `DAT.UNF`.
- `xmldis` stops in PTPL at a leftover `STOP 'TESTING'` (xmlplot.for:553) when `view_template = yes`.
- `xmldis` stops if `<name>Nuclides.unf` is missing (true for any .sp spectrum).
- `/FITTEDCONS/` layout differs between files (some have PCAP/LPCAP, `CnmBli.for` doesn't).
- Makefile: on Linux the `peakfit` link needs `-lm`; FC and the libgfortran path are hardcoded to Homebrew.
