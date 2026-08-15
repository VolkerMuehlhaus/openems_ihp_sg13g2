# Stackup XML Format

This document describes the XML file format read by `util_stackup_reader.read_substrate()`.
It defines the process stackup used to build a 3D simulation model from GDSII: materials,
the dielectric stack, the drawn metal/via layers and their GDSII layer numbers, optional
derived (computed) layers, and optional thermal-conductivity tables.

Reference files in this repo: [`SG13G2_200um.xml`](../SG13G2_200um.xml) (basic stackup),
[`SG13G2_resistors_200um.xml`](../SG13G2_resistors_200um.xml) (adds resistor materials and
derived layers), [`SG13G2_200um_test_derived.xml`](../SG13G2_200um_test_derived.xml) (derived
layer examples including `SIZE`).

## Top-level structure

```xml
<?xml version="1.0" encoding="UTF-8" standalone="no" ?>
<Stackup schemaVersion="2.0">
  <Materials>
    ...
  </Materials>
  <ELayers LengthUnit="um">
    <Dielectrics>
      ...
    </Dielectrics>
    <Layers>
      ...
    </Layers>
    <DerivedLayers>       <!-- optional -->
      ...
    </DerivedLayers>
  </ELayers>
  <Tables>                 <!-- optional -->
    ...
  </Tables>
</Stackup>
```

`schemaVersion` is informational only — the reader does not currently check or branch on it.

## `<Materials>`

One `<Material>` element per material, referenced by name from `<Dielectric>` and `<Layer>`
elements later in the file.

| Attribute                 | Required | Default | Description |
|----------------------------|----------|---------|--------------|
| `Name`                     | yes      | —       | Material name, referenced elsewhere as `Material="..."`. |
| `Type`                     | yes      | —       | `Conductor`, `Dielectric`, `Semiconductor`, or `Resistor`. |
| `Permittivity`              | no       | `1`     | Relative permittivity (εr). |
| `DielectricLossTangent`     | no       | `0`     | Loss tangent. |
| `Conductivity`               | no       | `0`     | Conductivity in S/m. Used for `Conductor`/`Semiconductor` materials. |
| `Rs`                        | no       | `0`     | Sheet resistance in Ω/square. Used for `Resistor` materials (paired with a zero-thickness `Type="sheet"` layer). |
| `Density`                   | no       | `1`     | Mass density, used by thermal simulation setup. |
| `ThermalConductivity`        | no       | `0`     | Constant thermal conductivity, used by thermal simulation setup. |
| `ThermalConductivityTable`   | no       | —       | Name of a `<Table>` (see [Tables](#tables-thermal-conductivity) below) to use instead of a constant value, for temperature-dependent conductivity. |
| `Color`                     | no       | —       | Hex RGB color (no `#`), used for 3D preview / GUI display. |

```xml
<Materials>
  <Material Name="Metal1" Type="Conductor" Permittivity="1" DielectricLossTangent="0"
            Conductivity="21640000.0" Color="39bfff"/>
  <Material Name="SiO2" Type="Dielectric" Permittivity="4.1" DielectricLossTangent="0.0"
            Conductivity="0" Color="fffcad"/>
  <Material Name="Substrate" Type="Semiconductor" Permittivity="11.9" DielectricLossTangent="0"
            Conductivity="2.0" Color="01e0ff"/>
  <Material Name="RSIL" Type="Resistor" Rs="7" Color="d0d0d0"/>
</Materials>
```

## `<ELayers>`

Container for the dielectric stack and drawn layers. `LengthUnit` sets the unit for all
`Thickness`/`Zmin`/`Zmax`/`Offset`/`Oversize` values in this section (typically `"um"`).

### `<Dielectrics>`

Defines the vertical dielectric stack, independent of GDSII — these layers exist everywhere
and are not drawn as polygons. Order in the file is **top to bottom**; z-positions are
computed automatically by stacking `Thickness` values from the top down (`calculate_zpositions()`
walks the list in reverse).

| Attribute   | Required | Description |
|-------------|----------|-------------|
| `Name`      | yes      | Dielectric name. |
| `Material`  | yes      | References a `<Material>` by name. |
| `Thickness` | yes*     | Layer thickness. Used to auto-stack this dielectric relative to its neighbors. |
| `Zmin`/`Zmax` | no     | Explicit absolute z-position instead of `Thickness`. If both are given, they take precedence over `Thickness` and this dielectric is excluded from the auto-stacking pass. |
| `Boundary`  | no       | GDSII layer number that defines this dielectric's finite (x,y) extent in the model. Optional — omit for a dielectric that simply fills the whole simulation domain. |

\* `Thickness` is required unless `Zmin`/`Zmax` are both given.

```xml
<Dielectrics>
  <Dielectric Name="AIR" Material="AIR" Thickness="200.0000" />
  <Dielectric Name="Passive" Material="Passive" Thickness="0.4000" />
  <Dielectric Name="SiO2" Material="SiO2" Thickness="15.7303" />
  <Dielectric Name="EPI" Material="EPI" Thickness="3.7500" />
  <Dielectric Name="Substrate" Material="Substrate" Thickness="180.0000" />
</Dielectrics>
```

With an explicit boundary layer instead of the full domain:

```xml
<Dielectric Name="Package" Material="MoldCompound" Thickness="500" Boundary="299"/>
```

### `<Layers>`

Drawn layers: metals, vias, sheet resistors, and dielectric "bricks" that are read from
specific GDSII layer numbers instead of being implied by the dielectric stack.

An optional `<Substrate Offset="..."/>` element shifts the z-position of every `<Layer>` by
a fixed amount (used e.g. to place the drawn stack on top of a backside/substrate region
that is itself modeled with negative z).

| Attribute | Required | Description |
|-----------|----------|--------------|
| `Name`    | yes      | Layer name, used for lookups (`getbylayername`) and in port definitions (`from_layername`/`to_layername` etc. elsewhere in the pipeline). |
| `Type`    | yes      | `conductor`, `via`, `dielectric`, or `sheet` (see below). |
| `Material`| yes      | References a `<Material>` by name. |
| `Zmin`    | yes      | Bottom z-position. |
| `Zmax`    | yes      | Top z-position. Equal to `Zmin` forces `Type` to `sheet` regardless of the stated `Type`. |
| `Layer`   | yes      | GDSII layer number. Also used as the target layer number for a [derived layer](#derivedlayers-boolean-operations-on-layers), in which case its geometry does not need to exist directly in the GDSII file. |

Layer `Type` meanings:
- **`conductor`** — a normal metal layer with thickness (`Zmax > Zmin`).
- **`via`** — a vertical connection between two conductor layers; eligible for via-array
  merging (`merge_polygon_size` in `read_gds()`).
- **`dielectric`** — a drawn (GDSII-sourced) dielectric region, as opposed to the implicit
  stack in `<Dielectrics>`.
- **`sheet`** — zero-thickness layer (`Zmin == Zmax`), typically paired with a `Resistor`
  material for sheet-resistance elements.

```xml
<Layers>
  <Substrate Offset="183.75"/>
  <Layer Name="Metal1" Type="conductor" Zmin="1.0400" Zmax="1.4600" Material="Metal1" Layer="8" />
  <Layer Name="Via1"   Type="via"       Zmin="1.4600" Zmax="2.0000" Material="Via1"   Layer="19" />
  <Layer Name="LBE"    Type="dielectric" Zmin="-183.75" Zmax="0"    Material="AIR"    Layer="157" />
  <Layer Name="RSIL"   Type="sheet"     Zmin="0.4"     Zmax="0.4"   Material="RSIL"   Layer="314"/>
</Layers>
```

### `<DerivedLayers>` (boolean operations on layers)

Defines new layer numbers whose geometry is *computed* from other layers (boolean
operations and/or resizing) instead of being read directly from GDSII — e.g. an overlap
between two drawn layers, used here to derive resistor geometry from poly/implant/contact
layers. Full reference, including chaining and self-touching polygon handling, is in
[`derived_layers.md`](derived_layers.md); summary:

- Give the derived layer's target `Layer` number a normal `<Layer>` entry too (as above),
  so it gets a Z-position/material and is included in `metals_list.getlayernumbers()`.
- `Operation` is one of `AND`, `OR`, `XOR`, `NOT` (≥2 `<Operand>` children, folded pairwise
  in order — order matters for `NOT`), or `SIZE` (exactly 1 operand, requires non-zero
  `Oversize`, pure resize with no boolean op).
- `Oversize` (optional, any operation) grows (positive) or shrinks (negative) the result by
  a fixed distance in layout units, applied after the boolean step.
- An `<Operand Layer="...">` can reference a native GDSII layer number or another derived
  layer's number; dependencies are resolved automatically regardless of definition order.

```xml
<DerivedLayers>
  <!-- layer 240 = layer 134 (TopMetal2) AND layer 126 (TopMetal1) -->
  <DerivedLayer Name="TM1_TM2_Overlap" Layer="240" Operation="AND">
    <Operand Layer="134"/>
    <Operand Layer="126"/>
  </DerivedLayer>

  <!-- layer 241 = layer 134 minus layer 126 -->
  <DerivedLayer Name="TM2_minus_TM1" Layer="241" Operation="NOT">
    <Operand Layer="134"/>
    <Operand Layer="126"/>
  </DerivedLayer>

  <!-- layer 242 = layer 126, grown by 2 units on every side, no boolean op -->
  <DerivedLayer Name="TM1_grown" Layer="242" Operation="SIZE" Oversize="2">
    <Operand Layer="126"/>
  </DerivedLayer>
</DerivedLayers>
```

Three or more operands are allowed for `AND`/`OR`/`XOR`, folded left to right — used here to
intersect three layers into one resistor recognition layer:

```xml
<DerivedLayer Name="Rhigh_or_Rppd" Layer="311" Operation="AND">
  <Operand Layer="28"/>
  <Operand Layer="14"/>
  <Operand Layer="310"/>  <!-- 310 is itself a derived layer, resolved first -->
</DerivedLayer>
```

## `<Tables>` (thermal conductivity)

Optional, sibling of `<Materials>` and `<ELayers>` directly under `<Stackup>`. Defines named
temperature/value lookup tables for materials whose thermal conductivity depends on
temperature (referenced from `<Material ThermalConductivityTable="...">` instead of a
constant `ThermalConductivity`).

```xml
<Tables>
  <Table Name="Si_k_vs_T">
    <Point Temperature="250" Value="170"/>
    <Point Temperature="300" Value="148"/>
    <Point Temperature="350" Value="130"/>
  </Table>
</Tables>
```

```xml
<Material Name="Substrate" Type="Semiconductor" Permittivity="11.9" Conductivity="2.0"
          ThermalConductivityTable="Si_k_vs_T" Color="01e0ff"/>
```

## What `read_substrate()` returns

```python
materials_list, dielectrics_list, metals_list = stackup_reader.read_substrate(XML_filename)
```

- `materials_list` — all `<Material>` entries (`stackup_materials_list`).
- `dielectrics_list` — all `<Dielectric>` entries, with z-positions resolved
  (`dielectric_layers_list`).
- `metals_list` — all `<Layer>` entries, sorted by z with above/below neighbors resolved
  (`metal_layers_list`). `metals_list.getlayernumbers()` is the layer list normally passed to
  `gds_reader.read_gds()`.
- `metals_list.derived_layers` — parsed `<DerivedLayers>` entries (empty if the XML has none).
  Not a separate return value, so this signature stays compatible with older 3-value call
  sites; `read_gds()` picks it up automatically from `metals_list` if not passed explicitly.
