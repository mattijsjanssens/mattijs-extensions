Block moving out of cylinder

2020-09-04: problem with mesh motion. Points on cyclicACMI seem to move:

- plot cellDisplacement on e.g. central_top_blockage
  (e.g. with Paraview 'plot over line'). Varies linearly.
- plot pointDiplacement: not linear
- problem is volPointInterpolation:
    - points on central_top_blockage are on 6 boundary faces (=ok):
        - empty (front and back)
        - duplicate ACMI and blockage
        - and that all to left and right
    - this gets accounted for in e.g. makeBoundaryWeights using isPatchPoint_:
      only two patch faces actually get used (=ok):
        weights: 6(0 0 26.3637 26.3637 0 0)
    - this gives sumWeights: 52.7273 (=ok)

    - but now addSeparated calls initSwapAddSeparated/swapAddSeparated
      on the ACMI patches and this additionally adds additional contributions
      to give sumWeights:238.423 !!!
    - problem is that sumWeights does not know we've overridden the bcs on
      pointDisplacement.
    - hacked for now by restoring the sumWeights on cyclicACMI

2020-09-16: volPointInterpolation solved in develop. Problem is now cyclicACMI
    since top baffle is connecting to different region.

