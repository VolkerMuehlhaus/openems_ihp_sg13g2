These scripts support the workflow when using openEMS. 

**build_pypi_readme.py** regenerates `README_pypi.md` at the repo root from `README.md`, rewriting relative `./doc/...` links to absolute GitHub URLs so images render on the PyPI package page. Run this before `python -m build`, from the repo root: `python scripts/build_pypi_readme.py`

**sync_local_modules_copy.py** copies `workflow/modules/` over `more_examples/local_modules_copy/modules/`, so that self-contained example (demonstrating the local "modules" folder distribution method, as an alternative to the `gds2openEMS` PyPI package) doesn't silently drift out of sync with the canonical source. It's a real copy, not a symlink, since GitHub's "Download ZIP" doesn't reliably preserve symlinks across platforms. Run this whenever `workflow/modules/` changes, from the repo root: `python scripts/sync_local_modules_copy.py`

**build_userguide_pdf.py** regenerates `doc/Using_OpenEMS_Python_with_IHP_SG13G2_v3.pdf` from the Markdown source at `doc/userguide_md_format/Using_OpenEMS_Python_with_IHP_SG13G2_v3.md` — the Markdown is the source of truth, edit that and then run this script to refresh the PDF. Requires `pandoc` (<https://pandoc.org/>) on `PATH` and the `typst` Python package (`pip install typst`) as the PDF-rendering engine — no separate LaTeX or GTK/Pango install needed. Relative links to other repo files (e.g. `../XML_stackup_format.md`) are rewritten to absolute GitHub URLs, since they'd otherwise point at the wrong local path once the built file lives in `doc/` instead of `doc/userguide_md_format/`. Run this whenever the guide changes, from the repo root: `python scripts/build_userguide_pdf.py`

**build_userguide_odt.py** regenerates `doc/Using_OpenEMS_Python_with_IHP_SG13G2_v3.odt` from the same Markdown source, for anyone who wants to finish editing or export to PDF themselves in LibreOffice Writer. Requires only `pandoc` on `PATH` — no PDF-rendering engine needed, since pandoc writes ODT natively. Same relative-link rewriting as `build_userguide_pdf.py`. Run from the repo root: `python scripts/build_userguide_odt.py`

**deembed_openEMS.py** is a script to deembed parasitic inductance of ports by cascading negative series L at each port. Output is written to a new file with suffix "_deembedded". 

Usage:

`
python3 deembed_openEMS.py inputfile.s2p
`

The value of parasitic port inductance is calculated from port geometry, using thin sheet approximation. To do so, this script requires geometry information data created by the latest openEMS workflow version, located in the same directory as the EM simulation result file. There is not limit on the number of ports. This is an experimental feature.
