These scripts support the workflow when using openEMS. 

**build_pypi_readme.py** regenerates `README_pypi.md` at the repo root from `README.md`, rewriting relative `./doc/...` links to absolute GitHub URLs so images render on the PyPI package page. Run this before `python -m build`, from the repo root: `python scripts/build_pypi_readme.py`

**deembed_openEMS.py** is a script to deembed parasitic inductance of ports by cascading negative series L at each port. Output is written to a new file with suffix "_deembedded". 

Usage:

`
python3 deembed_openEMS.py inputfile.s2p
`

The value of parasitic port inductance is calculated from port geometry, using thin sheet approximation. To do so, this script requires geometry information data created by the latest openEMS workflow version, located in the same directory as the EM simulation result file. There is not limit on the number of ports. This is an experimental feature.
