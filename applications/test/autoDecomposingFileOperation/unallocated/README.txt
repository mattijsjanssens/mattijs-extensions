unallocatedFvMesh : mesh class. Just holds objectRegistry reference
and mesh properties:
    - nCells, nFaces etc.
    - boundary with unallocatedGenericFvPatch each with
        - start, size
        - name

- all patches get read as a unallocatedGenericFvPatch
(with type 'unallocatedGeneric')
- including e.g. empty or processor patches!

- so now construct-from-patch will complain
