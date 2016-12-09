fvMesh& mesh = solidRegions[i];

volScalarField& T = Ts[i];
const dimensionedScalar& DT = DTs[i];
fv::options& fvOptions = solidHeatSources[i];
