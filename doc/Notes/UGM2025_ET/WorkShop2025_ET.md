---
marp: true
theme: ESI
paginate: true
_header: ""
footer: Mattijs Janssens, OpenFOAM Workshop 2025
title:  Expression Templates
header: Expression Templates
auto-scaling: false


---

# Expression templates

- OpenFOAM is field-based syntax
- hides complex structures in formula-like syntax:
```
volScalarField c(...);
c = a + sqrt(b);
```
- all fields are a set of values for internal fields and
special handling of the boundary ('fvPatchFields')
- above expression does
  - allocate new temporary field
  - fill with sqrt(b)
  - in-place add a
  - copy values to c
- *lots* of intermediate fields: show example
- remarkably little overhead on CPU (`kcachegrind`, few %)

---

## Expression templates
- [Expression Templates](https://en.wikipedia.org/wiki/Expression_templates) : define `operator[]` that encodes the whole operation:
  ```
  scalar add<field, sum<sqrt<field>>>::operator[i]
  {
       return a[i] + sqrt[b[i]];
  }
  ```
  and use it to evaluate:
  ```
  forAll(c, i)
  {
      c[i] = add<field, sum<sqrt<field>>>::operator[](i);
  }
  ```
---

- advantages
  - avoid intermediate storage (memory allocation)
  - avoid out-of-band access
  - easier to multi-thread
  - at evaluation time decide if bc's are required
- disadvantages
  - hard to debug
  - does not allow boundary conditions on removed intermediate fields (e.g. surface fields)

---

- simple field algebra:

  - adding two fields (c = a + b)
  - Intel E5-2620

  | Size     | Intermediate field| Expression templates|
  |----------|-------------------|---------------------|
  | 100000   |  0.000097         | 0.000083 |
  | 1000000  |  0.00287519       | 0.00110006 |
  | 10000000 |  0.0268587        | 0.014261 |

---
```
c = a + b
```
![Simple algebra](./figures/simple_algebra.png "Simple algebra")

---

- more complex algebra
```
c = cos(a + 0.5*sqrt(b-sin(a)))
```
![Complex algebra](./figures/complex_algebra.png "Complex algebra")

---

## GeometricFields as expression templates
- already has operator[] ...
- no storage for GeometricField
- no registration
- is just expression of field operation
- and boundary conditions
- how to do I/O?

---

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

---

- wrapper syntax:
```
// Construct expression
Expression::GeometricFieldConstRefWrap<volScalarField> wa(a);
Expression::GeometricFieldConstRefWrap<volScalarField> wb(b);
const auto expression(wa + sqrt(wb));

// Evaluate expression into field c
Expression::GeometricFieldRefWrap<volScalarField> wc(c);
expression.evaluate(wc);
```

---

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

---

- only explicit finiteVolume discretisation
- has already small CPU benefit

  pitzDaily tutorial:

  | Gauss| fusedGauss|
  |------|-----------|
  | 8.22 |  7.62     |

  (identical residuals at 6 digits)

---

<!----
## Other wrappers:
- constant
- tmp
- fvMatrix (linear operations only)
- discretisation : cell-to-face linear interpolation

---

--->

## Wrapping it up
- `.expr()` to create the wrapper
- assignment to evaluate the wrapper
  ```
  c = a.expr() + sqrt(b.expr());
  ```

---

- wrapper generators

| Class| Expression | Assignment |
|------|-----------|-------------|
| List | .expr() | 	yes |
| Field | 	.expr() | 	yes |
| GeometricField (e.g. volScalarField) | 	.expr() | 	yes |
| tmp\<Field\> | 	.expr() | 	no |
| tmp\<GeometricField\> |	.expr() | 	no |
| DimensionedType (for constant Field) | 	.expr(\<size\>) | 	no |
| DimensonedType (for constant GeometricField) | 	.expr(\<GeoField\>) | 	no |
| fvMatrix | 	.expr() | 	yes |

---

## Detail : constants

- when used in `Field` expression : only size (and value) required:
  ```
  const Expression::UniformListWrap<scalar> two(mesh.nCells(), 2.0);
  ```
  or
  ```
  const auto two(dimensionedScalar(dimless, 2.0).expr(mesh.nCells()));
  ```
- when used in `GeometricField` expression : needs GeometricField reference:
  ```
  const auto two(dimensionedScalar(dimless, 2.0).expr(mesh.magSf()));
  ```
  (GeometricField needs to be 'live' at evaluation time)

---

## Discretisation

| Function | 	Input   | Output |
|----------|----------|--------|
| Interpolation	| expression template | expression template |
| uncorrected Gauss Laplacian	| fvMatrix, expression template | fvMatrix |

---

## Detail : Gauss Laplacian

```
// Get expression for difference weights
const auto deltaCoeffs = this->tsnGradScheme_().deltaCoeffs(vf).expr();

// Get expression for interpolation weights
const auto weights = this->tinterpGammaScheme_().weights(gamma).expr();

// Interpolate gamma field and multiply with face-area magnitude
const auto gammaMagSf =
    Expression::lerp(gamma.expr(), weights, mesh)
  * mesh.magSf().expr();

// Create upper
fvm.upper() = deltaCoeffs.internalField()*gammaMagSf.internalField()
```

---

## More timings
- volField operations
  - no threads
  - compile with threads

---

## Typical application: kOmegaSST F1 function
```
auto arg1 = min
(
    min
    (
        max
        (
            (one/betaStar)*sqrt(k_.expr())/(omega_.expr()*y_.expr()),
            fiveHundred
           *(this→mu().expr()/this→rho_.expr())
           /(sqr(y_.expr())*omega_.expr())
        ),
        fourAlphaOmega2*k_.expr()/(CDkOmegaPlus*sqr(y_.expr()))
    ),
    ten
);
return tanh(pow4(arg1));
```

---

## Adapt for offloading
- bottom level evaluation is
  ```
  for (label i = 0; i < lst.size(); ++i)
  {
      lst[i] = operator[](i);
  }
- replace with
  ```
  std::copy
  (
      std::execution::par_unseq,
      static_cast<E const&>(*this).cbegin(),
      static_cast<E const&>(*this).cend(),
      lst.begin()
  );
  ```
- requires random-access iterators

---

## Future work
- extend to all operators, functions
- fix fvMatrix source terms
- handle type-changing code (scalar*vector)
- have _expr functions for finiteVolume discretisation?
- apply!

---
