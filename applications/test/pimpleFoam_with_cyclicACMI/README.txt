2021-12-12:

- temperature implicit coupled through cyclicAMI
- but p,U e.g. normal wall
- hence needs to be on a lduInterface
1. make all fvPatches an lduInterface:
    - but then all get included in GAMGAgglomeration
2. or have mixed wall (ACMI)
    - and only have cyclicAMI on T field
