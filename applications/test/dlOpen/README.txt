2017-05-26: with all libraries from $FOAM_LIBBIN it crashes trying 
to load libengine.so:

Crash:


Loading "libbarotropicCompressibilityModel.so"
Unloading "libbarotropicCompressibilityModel.so"
Loading "libfiniteVolume.so"
Unloading "libfiniteVolume.so"
Loading "libfileFormats.so"
Unloading "libfileFormats.so"
Loading "liblagrangianSpray.so"
Unloading "liblagrangianSpray.so"
Loading "libSLGThermo.so"
Unloading "libSLGThermo.so"
Loading "libsixDoFRigidBodyMotion.so"
Unloading "libsixDoFRigidBodyMotion.so"
Loading "libhelpTypes.so"
Unloading "libhelpTypes.so"
Loading "libmolecularMeasurements.so"
Unloading "libmolecularMeasurements.so"
Loading "libsolidSpecie.so"
Unloading "libsolidSpecie.so"
Loading "libphaseChangeTwoPhaseMixtures.so"
Unloading "libphaseChangeTwoPhaseMixtures.so"
Loading "libsolidParticle.so"
Unloading "libsolidParticle.so"
Loading "libengine.so"
==26216== Invalid read of size 4
==26216==    at 0x10B6B8B7: Foam::HashTable<Foam::autoPtr<Foam::motionSolver> (*)(Foam::polyMesh const&, Foam::IOdictionary const&), Foam::word, Foam::string::hash>::set(Foam::word const&, Foam::autoPtr<Foam::motionSolver> (* const&)(Foam::polyMesh const&, Foam::IOdictionary const&), bool) (in /home/penfold2/mattijs/OpenFOAM/OpenFOAM-dev.feature-globalFile/platforms/linux64GccDPInt32Opt/lib/libdynamicMesh.so)
==26216==    by 0x11D0AD07: Foam::motionSolver::adddictionaryConstructorToTable<Foam::displacementSBRStressFvMotionSolver>::adddictionaryConstructorToTable(Foam::word const&) (in /home/penfold2/mattijs/OpenFOAM/OpenFOAM-dev.feature-globalFile/platforms/linux64GccDPInt32Opt/lib/libfvMotionSolvers.so)
==26216==    by 0x11CF5C9D: _GLOBAL__sub_I_displacementSBRStressFvMotionSolver.C (in /home/penfold2/mattijs/OpenFOAM/OpenFOAM-dev.feature-globalFile/platforms/linux64GccDPInt32Opt/lib/libfvMotionSolvers.so)

So:

- app loads libfoamToVTK.so which loads libdynamicMesh.so (has RTS table for motionSolvers)
- app unloads libfoamToVTK.so
- libengine loads libfvMotionSolvers.so
- fvMotionSolvers.so adds displacementSBRStressFvMotionSolver to motionSolvers (in libdynamicMesh) and this fails
=
