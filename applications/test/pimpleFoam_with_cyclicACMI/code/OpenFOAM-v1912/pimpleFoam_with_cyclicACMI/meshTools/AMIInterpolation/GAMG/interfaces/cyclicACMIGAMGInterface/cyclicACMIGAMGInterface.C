/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "AMIInterpolation.H"
#include "cyclicACMIGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "Map.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicACMIGAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicACMIGAMGInterface::cyclicACMIGAMGInterface
(
    const label interfacei,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& allNeighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
:
    GAMGInterface
    (
        interfacei,
        coarseInterfaces
    ),
    fineCyclicACMIInterface_
    (
        refCast<const cyclicACMILduInterface>(fineInterface)
    ),
    neighbSize_(-1)
{
Pout<< "Patch:" << interfacei << " fineSize:"  << fineInterface.faceCells().size()
    << endl;
//for (const auto& coarseIf : coarseInterfaces)
//{
//    Pout<< "    coarse size:" << coarseIf.faceCells().size() << endl;
//}
Pout<< "localRestrictAddressing:" << localRestrictAddressing << endl;
Pout<< "allNeighbourRestrictAddressing:" << allNeighbourRestrictAddressing
    << endl;


    // Construct face agglomeration from cell agglomeration
    {
        // From coarse face to cell
        DynamicList<label> dynFaceCells(localRestrictAddressing.size());

        // From face to coarse face
        DynamicList<label> dynFaceRestrictAddressing
        (
            localRestrictAddressing.size()
        );

        Map<label> masterToCoarseFace(localRestrictAddressing.size());

        for (const label curMaster : localRestrictAddressing)
        {
            const auto iter = masterToCoarseFace.cfind(curMaster);

            if (iter.found())
            {
                // Already have coarse face
                dynFaceRestrictAddressing.append(iter.val());
            }
            else
            {
                // New coarse face
                const label coarseI = dynFaceCells.size();
                dynFaceRestrictAddressing.append(coarseI);
                dynFaceCells.append(curMaster);
                masterToCoarseFace.insert(curMaster, coarseI);
            }
        }

        faceCells_.transfer(dynFaceCells);
        faceRestrictAddressing_.transfer(dynFaceRestrictAddressing);
    }


    // On the owner side construct the AMI

    const labelList& nbrIds = fineCyclicACMIInterface_.neighbPatchIDs();

    amiPtrs_.setSize(nbrIds.size());

    label n = 0;
    for (const label index : fineCyclicACMIInterface_.AMIIndices())
    {
        const auto tami(fineCyclicACMIInterface_.AMI(index));

        if (tami.valid())
        {
DebugVar(tami().srcAddress());
DebugVar(tami().tgtAddress());
            labelList nbrFaceRestrictAddressing;
            {
                //const auto& nbr = fineCyclicACMIInterface_.neighbPatch(index);
                const label nbrSize = tami().tgtAddress().size();
Pout<< "    offset:" << n
    << " size:" << nbrSize
    << endl;
                const SubList<label> neighbourRestrictAddressing
                (
                    allNeighbourRestrictAddressing,
                    nbrSize,
                    n
                );
                n += nbrSize;

                // From face to coarse face
                DynamicList<label> dynNbrFaceRestrictAddressing
                (
                    neighbourRestrictAddressing.size()
                );

                Map<label> masterToCoarseFace
                (
                    neighbourRestrictAddressing.size()
                );

                for (const label curMaster : neighbourRestrictAddressing)
                {
                    const auto iter = masterToCoarseFace.cfind(curMaster);

                    if (iter.found())
                    {
                        // Already have coarse face
                        dynNbrFaceRestrictAddressing.append(iter.val());
                    }
                    else
                    {
                        // New coarse face
                        const label coarseI = masterToCoarseFace.size();
                        dynNbrFaceRestrictAddressing.append(coarseI);
                        masterToCoarseFace.insert(curMaster, coarseI);
                    }
                }

                nbrFaceRestrictAddressing.transfer
                (
                    dynNbrFaceRestrictAddressing
                );
            }
DebugVar(faceRestrictAddressing_);
DebugVar(nbrFaceRestrictAddressing);
            amiPtrs_.set
            (
                index,
                new AMIPatchToPatchInterpolation
                (
                    tami(),
                    faceRestrictAddressing_,
                    nbrFaceRestrictAddressing
                )
            );
        }
    }

    Pout<< "** done patch:" << interfacei
        << " finesize:" << fineInterface.faceCells().size()
        << " coarsesize:" << faceCells().size()
        << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicACMIGAMGInterface::~cyclicACMIGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cyclicACMIGAMGInterface::neighbSize() const
{
    if (neighbSize_ == -1)
    {
        neighbSize_ = 0;
        for (const label index : AMIIndices())
        {
            neighbSize_ += neighbPatch(index).size();
        }
    }
    return neighbSize_;
}


Foam::tmp<Foam::labelField>
Foam::cyclicACMIGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    // Return internal field (e.g. cell agglomeration) in nbr patch index

    tmp<labelField> tpnf(new labelField(neighbSize()));
    labelField& pnf = tpnf.ref();

    label n = 0;
    for (const label index : AMIIndices())
    {
        const labelUList& nbrFaceCells = neighbPatch(index).faceCells();
        for (const auto celli : nbrFaceCells)
        {
            pnf[n++] = iF[celli];
        }
    }

    return tpnf;
}


// ************************************************************************* //
