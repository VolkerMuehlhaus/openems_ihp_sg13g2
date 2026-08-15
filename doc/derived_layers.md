# Derived Layers

A derived layer is a new layer number whose geometry is *computed* — via a boolean
operation and/or a resize (grow/shrink) — from other layers, instead of being read
directly from GDSII. Derived layers are defined in the stackup XML file and are
resolved automatically by `gds_reader.read_gds()` (see `util_stackup_reader.py` and
`util_gds_reader.py`).

## Enabling derived layers

Two things are needed in the stackup XML, inside `<ELayers>`:

1. A normal `<Layer>` entry for the derived layer number, like any other metal. This
   is what gives the derived layer its Z-height and material for the 3D model.
2. A `<DerivedLayers>` section (a sibling of `<Layers>`) with one `<DerivedLayer>`
   element per computed layer.

```xml
<ELayers LengthUnit="um">
  ...
  <Layers>
    ...
    <Layer Name="TM1_TM2_Overlap" Type="conductor" Zmin="4" Zmax="5" Material="TopMetal1" Layer="240" />
  </Layers>
  <DerivedLayers>
    <DerivedLayer Name="TM1_TM2_Overlap" Layer="240" Operation="AND">
      <Operand Layer="134"/>
      <Operand Layer="126"/>
    </DerivedLayer>
  </DerivedLayers>
</ELayers>
```

No changes are needed in the model `.py` script — once `metals_list` is read from the
XML via `stackup_reader.read_substrate()` and passed into `gds_reader.read_gds()`,
derived layers are resolved automatically. Just make sure the derived layer's number
is included in the `layernumbers` list passed to `read_gds()` (it will be, as long as
it has a `<Layer>` entry, since that list comes from `metals_list.getlayernumbers()`).

## `<DerivedLayer>` attributes

| Attribute   | Required | Description |
|-------------|----------|-------------|
| `Name`      | yes      | Name for this derived layer (used in error/debug messages). |
| `Layer`     | yes      | The new layer number that this derived layer is stored as. |
| `Operation` | yes      | One of `AND`, `OR`, `XOR`, `NOT`, `SIZE` (see below). |
| `Oversize`  | no       | Distance (in layout units) to grow (positive) or shrink (negative) the result. Default `0` (no resize). Required and must be non-zero for `Operation="SIZE"`. |

Each `<DerivedLayer>` has one or more `<Operand Layer="..."/>` child elements, in order.
An operand can be a native GDSII layer number or another derived layer's number —
dependencies between derived layers are resolved automatically, so operand order
across `<DerivedLayer>` entries does not matter.

## Operations

### `AND`, `OR`, `XOR`

Standard boolean set operations, folded pairwise across all operands in order:
`(op1 <op> op2) <op> op3 ...`. Needs at least 2 operands. Operand order does not
matter for these three.

```xml
<!-- layer 240 = layer 134 (TopMetal2) AND layer 126 (TopMetal1) -->
<DerivedLayer Name="TM1_TM2_Overlap" Layer="240" Operation="AND">
  <Operand Layer="134"/>
  <Operand Layer="126"/>
</DerivedLayer>
```

### `NOT`

Subtracts all operands after the first from the first: `op1 - op2 - op3 ...`. Needs
at least 2 operands. **Operand order matters** — the first operand is the base shape.

```xml
<!-- layer 241 = layer 134 (TopMetal2) minus layer 126 (TopMetal1) -->
<DerivedLayer Name="TM2_minus_TM1" Layer="241" Operation="NOT">
  <Operand Layer="134"/>
  <Operand Layer="126"/>
</DerivedLayer>
```

### `SIZE`

Grows or shrinks a single operand by `Oversize`, with no boolean operation. Needs
exactly 1 operand and a non-zero `Oversize`. This is the readable way to put a
resized copy of an existing layer onto a new layer number.

```xml
<!-- layer 242 = layer 126 (TopMetal1), grown by 2 units on every side -->
<DerivedLayer Name="TM1_grown" Layer="242" Operation="SIZE" Oversize="2">
  <Operand Layer="126"/>
</DerivedLayer>
```

## Combining a boolean operation with `Oversize`

`Oversize` can also be added to `AND`/`OR`/`XOR`/`NOT` — it is applied to the boolean
result as a final step:

```xml
<!-- layer 243 = (layer 134 AND layer 126), then grown by 2 units -->
<DerivedLayer Name="TM1_TM2_Overlap_grown" Layer="243" Operation="AND" Oversize="2">
  <Operand Layer="134"/>
  <Operand Layer="126"/>
</DerivedLayer>
```

## Chaining derived layers

An operand can reference another derived layer's number instead of a native GDSII
layer. Processing order is resolved automatically (topological sort), so it does not
matter whether the referenced derived layer is defined earlier or later in the XML.
A circular dependency between derived layers is reported as an error.

```xml
<DerivedLayer Name="A" Layer="240" Operation="AND">
  <Operand Layer="134"/>
  <Operand Layer="126"/>
</DerivedLayer>
<DerivedLayer Name="A_grown" Layer="241" Operation="SIZE" Oversize="2">
  <Operand Layer="240"/>  <!-- uses derived layer A as its operand -->
</DerivedLayer>
```

## Self-touching ("keyhole") polygon results

If a boolean result is not simply-connected (e.g. it has a real hole, or a `NOT`
result splits into disjoint islands), it gets encoded as a single point sequence
that revisits a vertex via a zero-width bridge. This is handled correctly in the
mesh-generation step (`create_surfaces_from_polygon` in `util_simulation_setup.py`),
which is downstream of both native GDSII polygons and derived-layer results: it
splits the point sequence back into simple loops, reconstructs holes with an OCC
boolean, and creates a separate surface per disjoint piece. This is independent of
the derived-layer computation itself and applies to any polygon reaching
`add_metal_volumes`, native or derived.
