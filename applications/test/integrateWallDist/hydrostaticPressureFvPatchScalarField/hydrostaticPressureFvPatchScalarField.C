/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "hydrostaticPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "PatchEdgeFaceWave.H"
#include "patchIntegrateInfo.H"
#include "EdgeMap.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::hydrostaticPressureFvPatchScalarField::linearInterpolate
(
    const scalarField& faceCentres,
    const scalarField& edgeCentres,
    const scalarField& faceFld
) const
{
    const polyPatch& pp = patch().patch();

    tmp<scalarField> tsumFld(new scalarField(pp.nEdges(), 0));
    scalarField& sumFld = tsumFld.ref();
    scalarField sumWeights(pp.nEdges(), 0);

    // Sum local face contributions
    {
        const labelListList& edgeFaces = pp.edgeFaces();
        forAll(edgeFaces, edgei)
        {
            const labelList& eFaces = edgeFaces[edgei];

            forAll(eFaces, eFacei)
            {
                label facei = eFaces[eFacei];

                scalar d = mag(faceCentres[facei]-edgeCentres[edgei]);
                scalar w = 1.0/max(SMALL, d);
                sumFld[edgei] += w*faceFld[facei];
                sumWeights[edgei] += w;
            }
        }
    }


    // Sum coupled contributions
    {
        const globalMeshData& gd = pp.boundaryMesh().mesh().globalData();
        const indirectPrimitivePatch& cpp = gd.coupledPatch();

        // Fill coupled edge data
        scalarField cppSum(cpp.nEdges(), 0.0);
        scalarField cppWeight(cpp.nEdges(), 0.0);

        forAllConstIter(Map<label>, patchToCoupledEdge_, iter)
        {
            cppSum[iter()] = sumFld[iter.key()];
            cppWeight[iter()] = sumWeights[iter.key()];
        }

        // Add contributions
        globalMeshData::syncData
        (
            cppSum,
            gd.globalEdgeSlaves(),
            gd.globalEdgeTransformedSlaves(),
            gd.globalEdgeSlavesMap(),
            gd.globalTransforms(),
            plusEqOp<scalar>(),
            mapDistribute::transform()
        );
        globalMeshData::syncData
        (
            cppWeight,
            gd.globalEdgeSlaves(),
            gd.globalEdgeTransformedSlaves(),
            gd.globalEdgeSlavesMap(),
            gd.globalTransforms(),
            plusEqOp<scalar>(),
            mapDistribute::transform()
        );

        // Extract
        forAllConstIter(Map<label>, patchToCoupledEdge_, iter)
        {
            sumFld[iter.key()] = cppSum[iter()];
            sumWeights[iter.key()] = cppWeight[iter()];
        }
    }


    // Normalise
    forAll(sumWeights, edgei)
    {
        if (sumWeights[edgei] > VSMALL)
        {
            sumFld[edgei] /= sumWeights[edgei];
        }
    }

    return tsumFld;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hydrostaticPressureFvPatchScalarField::
hydrostaticPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    rhoName_("rho"),
    pRefValue_(0.0),
    pRefPoint_(Zero)
{}


Foam::hydrostaticPressureFvPatchScalarField::
hydrostaticPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    pRefValue_(readScalar(dict.lookup("pRefValue"))),
    pRefPoint_(dict.lookup("pRefPoint"))
{
    // Calculate the patch-specific information
    {
        const polyPatch& pp = patch().patch();
        const edgeList& edges = pp.edges();
        const labelList& mp = pp.meshPoints();

        // Build map from pp edges (in mesh indices) to edge indices
        EdgeMap<label> edgeToIndex(2*edges.size());
        {
            forAll(edges, edgei)
            {
                const edge& e = edges[edgei];
                edgeToIndex.insert(edge(mp[e[0]], mp[e[1]]), edgei);
            }
        }

        const polyMesh& mesh = pp.boundaryMesh().mesh();
        const globalMeshData& gd = mesh.globalData();
        const indirectPrimitivePatch& cpp = gd.coupledPatch();
        const edgeList& cppEdges = cpp.edges();
        const labelList& cppMp = cpp.meshPoints();

        patchToCoupledEdge_.resize(edges.size());

        forAll(cppEdges, cppEdgei)
        {
            const edge& e = cppEdges[cppEdgei];
            const edge meshE(cppMp[e[0]], cppMp[e[1]]);
            EdgeMap<label>::const_iterator iter = edgeToIndex.find(meshE);
            if (iter != edgeToIndex.end())
            {
                patchToCoupledEdge_.insert(iter(), cppEdgei);
            }
        }

DebugVar(patchToCoupledEdge_);
    }


    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::hydrostaticPressureFvPatchScalarField::
hydrostaticPressureFvPatchScalarField
(
    const hydrostaticPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_),
    pRefValue_(ptf.pRefValue_),
    pRefPoint_(ptf.pRefPoint_)
{}


Foam::hydrostaticPressureFvPatchScalarField::
hydrostaticPressureFvPatchScalarField
(
    const hydrostaticPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    rhoName_(ptf.rhoName_),
    pRefValue_(ptf.pRefValue_),
    pRefPoint_(ptf.pRefPoint_)
{}


Foam::hydrostaticPressureFvPatchScalarField::
hydrostaticPressureFvPatchScalarField
(
    const hydrostaticPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    rhoName_(ptf.rhoName_),
    pRefValue_(ptf.pRefValue_),
    pRefPoint_(ptf.pRefPoint_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hydrostaticPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    const polyPatch& pp = patch().patch();

    // Get component to integrate along
    const scalarField faceCentres(patch().Cf() & g.value());
    scalarField edgeCentres(pp.nEdges());
    {
        forAll(pp.edges(), edgei)
        {
            const edge& e = pp.edges()[edgei];
            const point& p0 = pp.points()[pp.meshPoints()[e[0]]];
            const point& p1 = pp.points()[pp.meshPoints()[e[1]]];
            edgeCentres[edgei] = ((0.5*(p0+p1)) & g.value());
        }
    }


     const fvPatchField<scalar>& faceRho =
         patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    //const scalarField faceRho(pp.size(), 1.0);


    const uniformDimensionedScalarField& hRef =
        db().lookupObject<uniformDimensionedScalarField>("hRef");

    dimensionedScalar ghRef
    (
        mag(g.value()) > small
      ? g & (cmptMag(g.value())/mag(g.value()))*hRef
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );



    // Integrate gh*rho from gh0
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    scalar gh0;
    scalarField pf(allFaceInfo.size());
    {
        tmp<scalarField> edgeRho
        (
            linearInterpolate
            (
                faceCentres,
                edgeCentres,
                faceRho
            )
        );


        List<patchIntegrateInfo> allEdgeInfo(pp.nEdges());
        List<patchIntegrateInfo> allFaceInfo(pp.size());

        // Seed edges
        DynamicList<label> changedEdges;
        DynamicList<patchIntegrateInfo> changedInfo;
        {
            // Find the starting edge
            label startEdgei = findMin(edgeCentres);
            gh0 = edgeCentres[startEdgei];
            changedEdges.append(startEdgei);
            changedInfo.append(patchIntegrateInfo(0.0));
        }

        // Static information
        patchIntegrateInfo::trackData td
        (
            faceCentres,
            edgeCentres,
            faceRho,
            edgeRho
        );


        // Walk
        PatchEdgeFaceWave
        <
            primitivePatch,
            patchIntegrateInfo,
            patchIntegrateInfo::trackData
        > calc
        (
            pp.boundaryMesh().mesh(),
            pp,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(pp.nEdges(), sumOp<label>()),
            td
        );

        scalarField pf(allFaceInfo.size());
        forAll(pf, facei)
        {
            pf[facei] = allFaceInfo[facei].value();
        }
    }

    Info<< type() << " : integrated rho*(gh-gh0) with gh0=" << gh0 << endl;


    //operator==(pf);

    operator==(pf - faceRho*(faceCentres - gh0));   //ghRef.value()));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::hydrostaticPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    os.writeKeyword("pRefValue") << pRefValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("pRefPoint") << pRefPoint_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        hydrostaticPressureFvPatchScalarField
    );
}

// ************************************************************************* //
