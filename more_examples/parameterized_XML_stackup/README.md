# Parameterized XML Stackup

Demonstrates overriding a `<Variable>` from the XML stackup file (`total_thickness`,
`air_thickness`) directly from a Python model script, without editing the XML — via the
`variable_overrides` parameter of `stackup_reader.read_substrate()`.

Both scripts model the same 2-port RF-CMIM capacitor (`rfcmim_30x15x10_full.gds`), read
the same stackup (`SG13G2_parameterized.xml`), and override the same two variables to
`100` µm each. They produce an equivalent openEMS model (identical mesh) - only the
Python coding style differs, to show both syntax options supported by this workflow:

| File | Style |
|---|---|
| [`run_rfcmim_2port_parameters.py`](run_rfcmim_2port_parameters.py) | Loose top-level variables (`total_thickness = 100`), passed to `read_substrate()` and the `simulation_setup` functions individually |
| [`run_rfcmim_2port_parameters_settings.py`](run_rfcmim_2port_parameters_settings.py) | A single `settings = {}` dictionary (`settings['variable_overrides']['total_thickness'] = 100`), passed as one object - see `workflow/run_generic_nport.py` for background on this style |

## See also

- [XML stackup file format reference](../../doc/XML_stackup_format.md) - full attribute-by-attribute
  description of the format read by `stackup_reader.read_substrate()`, including `<Variables>` and
  `"="`-expressions.
- [setupEM Stackup Editor](https://github.com/VolkerMuehlhaus/setupEM/blob/main/doc/setupEM_userguide.md#the-stackup-editor) -
  graphical tool for creating/editing these XML stackup files (materials, dielectric stack, layers,
  derived layers, Variables, thermal tables) instead of hand-editing the XML. setupEM otherwise
  targets the gds2palace/Elmer FEM workflows, not openEMS - it is mentioned here only for its
  stackup editor.
- [EMStudio](https://github.com/IHP-GmbH/EMStudio) - IHP's own GUI for the openEMS FDTD workflow
  itself, including overriding stackup `<Variable>`s like this example, plus its own GUI stackup
  editor. A good alternative to setupEM if you want a GUI that targets openEMS directly rather than
  the FEM solvers.
