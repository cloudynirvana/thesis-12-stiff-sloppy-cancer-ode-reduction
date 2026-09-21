# Stiff-sloppy spectra and systematic reduction of high-dimensional cancer-state ODEs under gated observation maps

**Thesis #12** (series label NP-04). Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Which directions of a high-dimensional cancer-state parameter space are sloppy under realistic observation maps, and can a documented MBAM / Fisher-geometry reduction preserve the gated claims of the full model without silently renaming leftover sloppy combinations as biology?

The calculations use one four-state toy with ten positive rates. Map B records log viable burden at eight times. Map L records the proliferating, quiescent and apoptotic states at those times. Stress is unobserved on both maps. The gated claim is the late log-change of viable burden, g = log B(40) - log B(11).

On this generator, map B has practical Fisher rank 3 of 10. Apoptotic clearance is an exact null coordinate. Map L has practical rank 6 of 10, and the two trailing eigenvectors lie nearly on the level set of the stress product kappa = h gp / gs. A quasi-steady reduction that keeps kappa reproduces g with relative error 0.0070. An eigenvector flow on map B can divide the proliferating-cell apoptosis rate by about 26 while the noise-free chi-square stays near 10<sup>-7</sup> and g does not move. A rebound term that spends the null clearance coordinate shifts g and is rejected by map B.

The trajectories are synthetic. They are not a fit to a cell line, and they are not a reduction of an unpublished 15-parameter field. This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/sloppy_reduction.py` | Seeded Fisher spectra, flows, and reductions (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Spectra, participation, flows, and burden curves |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/sloppy_reduction.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the spectra. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Stiff-sloppy spectra and systematic reduction of high-dimensional cancer-state ODEs under gated observation maps [Internet]. Thesis #12 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-12-stiff-sloppy-cancer-ode-reduction

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
