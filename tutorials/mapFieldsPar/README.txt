20160617: test mapFieldsPar

- use case with fvPatchFields with stored data and zero sizes of patch
  on at least one processor.

- consistent mapping with destination fields

    forwardStep_fine/
        mpirunDebug -np 5 mapFieldsPar ../forwardStep -consistent -parallel

- consistent mapping without destination fields

- inconsistent mapping

- map with -subtract (assumes destination fields)
