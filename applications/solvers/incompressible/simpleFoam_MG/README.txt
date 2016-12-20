MGsimpleFoam:
- level 1 setError:
    Uerr = mapper.interpolate(level0residual_with_error)-level1residual_wo_error

    ..

- level N setError

- solve coarsest level:
    - solve UEqn+error+grad(p)
    - solve pEqn-error

- other levels correct&solve:
    - update p,U,phi with solution from N
    - solve level (see above, with error)
    ..
- level 0 correct&solve
    - update p,U,phi with solution from 1
    - solve level



simpleFoam_MG:
- solve level0 (without error):
    - solve UEqn+grad(p)
    - solve pEqn
    - interpolate residual to coarse level
- solve other levels (with error)
    - solve UEqn+error+grad(p)
    - solve pEqn-error
    - interpolate residual to coarse level

