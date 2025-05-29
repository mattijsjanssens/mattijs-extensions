# Expression templates

- OpenFOAM is field-based syntax
- hides complex structures in formula-like syntax
- all fields are a set of values for internal fields and
special handling of the boundary ('fvPatchFields')
- *lots* of intermediate fields: show example
- remarkably little overhead on CPU (few %)


## Expression templates
- https:// expression templates. Have operator[]
- show class structure
- advantages
  - avoid intermediate storage (memory allocation)
  - avoid out-of-band access
  - easier to multi-thread
- disadvantages
  - hard to debug
  - does not allow boundary conditions on e.g. surface fields

- timing of Field algebra


## GeometricFields as expression templates
- already has operator[] ...
- no storage for GeometricField
- no registration
- is just expression of field operation
- and boundary conditions
- how to do I/O?


## (v2412) Expression template wrapper
- reference to GeometricField
- both expression templates and direct operations
    (dimensions, oriented)
- constructs expressions on this GeometricField
  - internal field : (operator on) straight list of values
  - uncoupled patch field : straight list of patch values
  - coupled patch field : blending of
    - indirect access of internal field
    - straight list of patch values

- show wrapper syntax:


## 'fused' discretisation
- exactly same purpose:
  - avoid intermediate surface fields
  - write single loop over complex expression
    ```
    const auto snGrad = [&]
    (
        const vector& Sf,

        const scalar weight,
        const scalar ownGamma,
        const scalar neiGamma,

        const scalar dc,
        const Type& ownVal,
        const Type& neiVal
    ) -> Type
    {
        const auto snGrad(dc*(neiVal-ownVal));
        const scalar faceGamma(weight*(ownGamma-neiGamma)+neiGamma);
        return mag(Sf)*faceGamma*snGrad;
    };

    fvc::surfaceSnSum
    (
        weights,
        gamma,

        deltaCoeffs,
        vf,

        snGrad,

        result,
        false       // avoid boundary evaluation until volume division
    );
    ```
- only explicit finiteVolume discretisation
- has already small CPU benefit


pitzDaily tutorial:

| Gauss| fusedGauss|
|------|-----------|
| 8.22 |  7.62     |

(identical residuals at 6 digits)


## Other wrappers:
- constant
- tmp
- fvMatrix (linear operations only)
- discretisation : cell-to-face linear interpolation


## Wrapping it up
- .expr() syntax
- assignment
- example of F1 function


## More timings
- volField operations
  - no threads
  - compile with threads
- compare F1


## Adapt for GPU
- std::execution::par_unseq
- requires random-access iterators (or legacy-forward-iterators?)


## Future work
- apply .expr() to code
- have _expr functions for finiteVolume discretisation?
