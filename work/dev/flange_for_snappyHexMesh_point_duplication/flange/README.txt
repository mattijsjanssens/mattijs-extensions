Run case with snappyHexMesh_debug. This has additional checking inside
processorPolyPatch::order to make sure number of points and edges stays same.
This fails after the final balancing.
Also:
- additional printing in fvMeshDistribute::mergeSharedPoints to make
sure there are no merged points
- meshRefinement::getCollocatedPoints reports zero found
- meshRefinement::dupNonManifoldPoints: disabled
- dumps mesh just before final distribution

Conclusion: is fvMeshDistribute. It adds the processors one-by-one
and doesn't remember the point connectivity, only the face connectivity.


Duplicate point is at

(-0.00374956564272283 0.00187499995809048 -0.0133268236124621)

Before balancing:
- processor0 : no duplicate points
- processor1 : far away
- processor2 : no duplicate points
- processor3 : far away

Before balancing:
- point is on boundary edge of processor on processors 0,2
  but also on multiple walls
