# Interactive 3D structure demo (Figure 7)

A self-contained, browser-based proof-of-concept that turns the **static PyMOL
renders of Figure 7** into a **rotatable 3D viewer**.

## What it shows

The 12 predicted CAZymes in Section 7 of `soda_rumen_metagenomics_analysis.Rmd`
(6 GH families × Soda Lake / Rumen) were structurally matched against the PDB
with Foldseek. This demo loads each one's **Foldseek best-hit PDB structure** so
you can rotate, zoom, and re-style it interactively.

| GH family | Soda Lake | Rumen |
|-----------|-----------|-------|
| GH1  | 2J75 | 4BCE |
| GH3  | 3SQL | 5XXN |
| GH9  | 3X17 | 6DHT |
| GH10 | 4K68 | 7D88 |
| GH28 | —    | 3JUR |

GH5_11 (both) and GH28 Soda Lake are omitted: their top Foldseek hits were
non-CAZyme scaffolds (marked `*` in the paper) rather than glycoside hydrolases.

## How to use

Open `figure7_3d_viewer.html` in any modern browser (double-click it).

- **Drag** to rotate, **scroll** to zoom.
- Switch between **cartoon / surface** representations.
- Colour by **rainbow** (N→C) or **secondary structure**.
- Toggle **auto-spin**.

> Requires an internet connection when opened: the page pulls `3Dmol.js` from a
> CDN and the coordinates from RCSB. (This repo's container has no outbound
> network, so the viewer can't be previewed server-side — open it locally.)

## Important caveat

These are the **experimental best-hit templates**, not the project's own
**AlphaFold3 models**. The actual `.cif` predictions were generated on Dardel
and are not committed to this repo. To view the real models instead, drop their
`.cif`/`.pdb` files next to this HTML and point `$3Dmol.download` /
`viewer.addModel` at the local file.
