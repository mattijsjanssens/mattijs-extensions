# Ensight format mesh reading

This release adds an Ensight Gold mesh importer is supported. It can handle all cell shapes supported by the `foamToEnsight` converter. It currently gets given the geometry file name (geo`extension`). It supports both ascii and binary formats.

Supported keywords:
- `extents`
- `node id 'given', 'ignore', 'assign'`
- `node_ids`
- `element id 'given', 'ignore', 'assign'`
- `element_ids`
- `part`
- `coordinates`
- `tetra4`
- `pyramid5`
- `penta6`
- `hexa8`
- `nfaced`
- `tria3`
- `quad4`
- `nsided`

It does not support
- 2D (finite-area) meshes
- `block` structured meshes
- quadratic elements (e.g. `twenty node hexahedron`)
- faceZones
- baffles (they probably get merged away unless they are in the first part - not tested)

It reads all parts, combines all the cells (`tetra3`, `hexa8` etc) and determines the outside faces. It merges all the points using a geometric test (see below) and uses all faces (`tria3` etc.) to patch any outside faces. Patch names are the original part names with any illegal word symbol replaced by '_'. Any remaining outside faces get added to a `defaultFaces` patch of type `empty`.

## Options
- `mergeTol` : supply optional merge tolerance to get the correspondence between points of different parts. Default is 1e-10 of the bounding box of all points. Specifying 0 disables any point merging (and hence patching).
- `scale` : specify optional scaling for the coordinates. Default is no scaling. Scaling can e.g. be used if the mesh is specified in [mm] instead of [m].
- `keepHandedness` : by default the mesh reader will flip (non-polyhedral) cells with negative volume. It will display warning messages of the form
```
zero or negative pyramid volume:
```
Use the flag to disable this check and use the normal vertex numbering.

Source code
- application/utilities/mesh/conversion/ensightToFoam
