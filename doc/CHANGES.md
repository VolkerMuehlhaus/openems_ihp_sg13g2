# Change list

This is an (incomplete) list of changes and new features.

## 20-Aug-2026
Fixed `write_snp()` writing a hardcoded 50 Ohm reference impedance into the Touchstone (`*.snp`) file header regardless of the actual port impedance used in the simulation. All model scripts now pass the real port reference impedance from `simulation_ports`, so the header correctly reflects non-50 Ohm setups (e.g. differential/GSG ports).

## 13-18-Aug-2026
Added support for derived layers in the stackup XML format, with a new IHP resistor example under `more_examples`.

The stackup reader (`util_stackup_reader.py`) now supports a `<Variables>` block and `=`-prefixed expressions, ported from gds2palace_ihp_sg13g2's reader to keep the two independent copies in sync. This allows a stackup XML file to define named values that other attributes can reference, instead of repeating the same physical value in multiple places.

Keyhole/hole polygons are now split before CSXCAD extrusion, fixing incorrect geometry for layers with holes. This requires the `shapely` package, now a documented dependency.

## 12-Jul-2026
Added `calculate_Zij()` for multiport Z-parameter extraction (previously only available for 2-port).

## 19-Jun-2026
New `numThreads` workflow option to force the number of openEMS solver threads, instead of relying on automatic thread count detection. Two examples are provided in the `more_examples` folder.
