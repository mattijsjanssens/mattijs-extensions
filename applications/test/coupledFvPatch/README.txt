ACMI:

value:
- Could be weighted weighted according to overlap.
  (but then  left and right of ACMI might have different values)
- But then gaussGrad is incorrect!!!
- So HAS to be unweighted value

matrixInterfaceUpdate: unweighted (coefficient is already scaled by area)

snGrad: unweighted

valueXXXCoeffs: unweighted

So currently there is no difference to the behaviour of AMI, apart from
the weights.

To do:
- make value the area weighted value since interpolation is now incorrect.
- move area contributions (gaussGrad, gaussLaplacianScheme)
  to coupledFvPatch so can be overridden.
