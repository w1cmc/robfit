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

## Result: peak fit NOT YET TESTABLE

Driven through the `libpeakfit` C API in 400-ch windows:
- `peakfit_fit` returns immediately (0 iterations, χ² = 0, all errors 0).
- With raw data (bridge hardcodes `IFIT = 1`), `peakfit_find_peaks` places evenly spaced fake peaks on top of the background.
- With `IFIT = 2` (data − continuum) it finds some real peaks (308, 398, 511 vs 2003's 310, 389/405, 512) but misses most.
- Out-of-bounds in `resl.for:119` (`AIDATMF` index iend+1); segfaults without bounds checking.
- `peakfit_set_calibration` is declared in `peakfit.h` but not implemented.
- COMMON `/INTVALS/` is 20 ints in `peakfit_bridge.for` vs 22 in `commons.for` (linker warns).

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
