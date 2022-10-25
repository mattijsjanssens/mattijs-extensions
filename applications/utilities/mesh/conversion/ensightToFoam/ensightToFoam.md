# Ensight format mesh reading

In this release an basic Ensight (mesh only) mesh importer is supported. It can handle all cell shapes supported by the `foamToEnsight` converter. It currently gets given the name of the `Ensight` directory and loads in there the `data/constant/geometry` file. E.g.

```
# Generate a mesh
blockMesh
# Create an Ensight file
foamToEnsight -ascii
# Read Ensight file and recreate the mesh
ensightToFoam Ensight
```
It currently
- does not run on a decomposed mesh
- expects above tree format to find the geometry file
- synthesises an empty patch to hold all unpatches faces
- geometrically matches patch parts to boundary faces
- can read both ascii and binary format


Source code
- application/utilities/mesh/conversion/ensightToFoam
