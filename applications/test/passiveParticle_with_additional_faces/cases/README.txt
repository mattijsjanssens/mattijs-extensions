20200405: two_tets:
    two overlapping tets. base is 0,0,0. Second triangle
    offset by z=0.5

    Inject particles say at (0.1 0.1 0.1) and track in +z direction.
    See if can make it jump into second triangle.


PROBLEM: even tets include the tet centre. Might as well use hexes.


E.g. on two_tets
- cell:
    (0 0 0)     // vertex0
    (1 0 0)
    (0 1 0)     // vertex2
    (0 0 1)
- centroid:
    (0.25 0.25 0.25)
  
- position (0.1 0.1 0.1) becomes
    0.4*centroid + 0.6*vertex0
transform:(0.25 0 1 0 0.25 0 0 0 0.25 0 0 1)
coord    :(0.4 0.6 0 2.68882e-17)


20200405: two_hex/
- problem: individual tets from hex decomposition should not have coincident
faces
- problem: cannot re-use tet-intersection since starting point might not
be inside the tet. Instead just use face/triangle intersection in absolute
coordinates?
