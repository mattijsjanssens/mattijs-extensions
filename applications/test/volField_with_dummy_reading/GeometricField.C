/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "GeometricField.H"
#include "Time.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "data.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define checkField(gf1, gf2, op)                                    \
if ((gf1).mesh() != (gf2).mesh())                                   \
{                                                                   \
    FatalErrorInFunction                                            \
        << "different mesh for fields "                             \
        << (gf1).name() << " and " << (gf2).name()                  \
        << " during operatrion " <<  op                             \
        << abort(FatalError);                                       \
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::readFields
(
    const dictionary& dict
)
{
    Internal::readField(dict, "internalField");

    boundaryField_.readField(*this, dict.subDict("boundaryField"));

    if (dict.found("referenceLevel"))
    {
        Type fieldAverage(pTraits<Type>(dict.lookup("referenceLevel")));

        Field<Type>::operator+=(fieldAverage);

        forAll(boundaryField_, patchi)
        {
            boundaryField_[patchi] == boundaryField_[patchi] + fieldAverage;
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::readFields()
{
    IOobject io
    (
        this->name(),
        this->time().timeName(),
        this->db(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );


    boolList headerOk(Pstream::nProcs(), false);
    headerOk[Pstream::myProcNo()] = io.headerOk();

    if (headerOk[Pstream::myProcNo()])
    {
        const IOdictionary dict(io, this->readStream(typeName));

        this->close();

        readFields(dict);
    }


    if (!returnReduce(headerOk[Pstream::myProcNo()], andOp<bool>()))
    {
        Pstream::gatherList(headerOk);
        //Pstream::scatterList(headerOk);


        // Not all have the file. Get a processor which does have it.
        label masterProci = 0;  //(headerOk ? Pstream::myProcNo() : labelMax);

DebugVar(masterProci);


        PstreamBuffers pBufs(Pstream::nonBlocking);
        if (Pstream::master())  // masterProci
        {
            // Construct zero-sized mesh for exchanging PatchFields
            autoPtr<typename GeoMesh::Mesh> dummyPtr
            (
                GeoMesh::New
                (
                    IOobject("dummy", this->db()),
                    this->mesh(),
                    false       // no comms
                )
            );
            const typename GeoMesh::Mesh& dummy = dummyPtr();


            // Bit of trickery to clone the field onto the dummy mesh.

            tmp<GeometricField<Type, PatchField, GeoMesh>> tresF;
            {
                // Create boundaryField with zero sized patches

                Boundary patchFields(dummy.boundary());

                forAll(patchFields, patchi)
                {
                    patchFields.set
                    (
                        patchi,
                        PatchField<Type>::New
                        (
                            PatchField<Type>::calculatedType(),
                            dummy.boundary()[patchi],
                            DimensionedField<Type, GeoMesh>::null()
                        )
                    );
                }

                tresF = tmp<GeometricField<Type, PatchField, GeoMesh>>
                (
                    new GeometricField<Type, PatchField, GeoMesh>
                    (
                        IOobject
                        (
                            "subset"+this->name(),
                            dummy.time().timeName(),
                            dummy,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        dummy,
                        this->dimensions(),
                        Field<Type>(this->primitiveField(), labelList(0)),
                        patchFields
                    )
                );
                GeometricField<Type, PatchField, GeoMesh>& resF = tresF.ref();

                // 2. Change the PatchFields to the correct type using
                //    a mapper
                //    constructor (with reference to the now correct
                //    internal field)

                typename GeometricField<Type, PatchField, GeoMesh>::
                    Boundary& bf = resF.boundaryFieldRef();


                //GeoMesh vm(dummy);

                forAll(bf, patchi)
                {
                    bf.set
                    (
                        patchi,
                        PatchField<Type>::New
                        (
                            this->boundaryField()[patchi],
                            dummy.boundary()[patchi],
                            resF(),
                            GeoMesh::NewMapper(dummy, patchi)    //vm.mapper()
                        )
                    );
                }
            }
            GeometricField<Type, PatchField, GeoMesh>& resF = tresF.ref();
DebugVar(resF);



            // Send resF to all processors
            for (label proci = 1; proci < Pstream::nProcs(); proci++)
            {
                if (!headerOk[proci])
                {
                    UOPstream str(proci, pBufs);
                    str << resF;
                }
            }
        }
        pBufs.finishedSends();

        if (!headerOk[Pstream::myProcNo()])
        {
            // Read field from master
            UIPstream str(Pstream::masterNo(), pBufs);
            dictionary dict(str);
            DebugVar(dict);
            readFields(dict);
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::readIfPresent()
{
    if
    (
        this->readOpt() == IOobject::MUST_READ
     || this->readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        WarningInFunction
            << "read option IOobject::MUST_READ or MUST_READ_IF_MODIFIED"
            << " suggests that a read constructor for field " << this->name()
            << " would be more appropriate." << endl;
    }
    else if (this->readOpt() == IOobject::READ_IF_PRESENT && this->headerOk())
    {
        readFields();

        // Check compatibility between field and mesh
        if (this->size() != GeoMesh::size(this->mesh()))
        {
            FatalIOErrorInFunction(this->readStream(typeName))
                << "   number of field elements = " << this->size()
                << " number of mesh elements = "
                << GeoMesh::size(this->mesh())
                << exit(FatalIOError);
        }

        readOldTimeIfPresent();

        return true;
    }

    return false;
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::readOldTimeIfPresent()
{
    // Read the old time field if present
    IOobject field0
    (
        this->name()  + "_0",
        this->time().timeName(),
        this->db(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE,
        this->registerObject()
    );

    if (field0.headerOk())
    {
        if (debug)
        {
            InfoInFunction << "Reading old time level for field"
                << endl << this->info() << endl;
        }

        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            field0,
            this->mesh()
        );

        field0Ptr_->timeIndex_ = timeIndex_ - 1;

        if (!field0Ptr_->readOldTimeIfPresent())
        {
            field0Ptr_->oldTime();
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const word& patchFieldType
)
:
    Internal(io, mesh, ds, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal(io, mesh, ds, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
:
    Internal(io, mesh, dt, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    boundaryField_ == dt.value();

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal(io, mesh, dt, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    boundaryField_ == dt.value();

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Internal& diField,
    const PtrList<PatchField<Type>>& ptfl
)
:
    Internal(io, diField),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, ptfl)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from components" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& iField,
    const PtrList<PatchField<Type>>& ptfl
)
:
    Internal(io, mesh, ds, iField),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, ptfl)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from components" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh
)
:
    Internal(io, mesh, dimless, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary())
{
    readFields();

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalIOErrorInFunction(this->readStream(typeName))
            << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    readOldTimeIfPresent();

    if (debug)
    {
        InfoInFunction
            << "Finishing read-construction of" << endl << this->info() << endl;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dictionary& dict
)
:
    Internal(io, mesh, dimless, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary())
{
    readFields(dict);

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalErrorInFunction
            << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    readOldTimeIfPresent();

    if (debug)
    {
        InfoInFunction
            << "Finishing dictionary-construct of "
            << endl << this->info() << endl;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    Internal(gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy" << endl << this->info() << endl;
    }

    if (gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            *gf.field0Ptr_
        );
    }

    this->writeOpt() = IOobject::NO_WRITE;
}


#ifndef NoConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
:
    Internal
    (
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp" << endl << this->info() << endl;
    }

    this->writeOpt() = IOobject::NO_WRITE;

    tgf.clear();
}
#endif


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


#ifndef NoConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
:
    Internal
    (
        io,
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp resetting IO params"
            << endl << this->info() << endl;
    }

    tgf.clear();

    readIfPresent();
}
#endif


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    Internal(newName, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting name"
            << endl << this->info() << endl;
    }

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            newName + "_0",
            *gf.field0Ptr_
        );
    }
}


#ifndef NoConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
:
    Internal
    (
        newName,
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp resetting name"
            << endl << this->info() << endl;
    }

    tgf.clear();
}
#endif


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const word& patchFieldType
)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, patchFieldType)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes

)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    )
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params and patch types"
            << endl << this->info() << endl;
    }

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


#ifndef NoConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal
    (
        io,
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    )
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp resetting IO params and patch types"
            << endl << this->info() << endl;
    }

    boundaryField_ == tgf().boundaryField_;

    tgf.clear();
}
#endif


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::~GeometricField()
{
    deleteDemandDrivenData(field0Ptr_);
    deleteDemandDrivenData(fieldPrevIterPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal&
Foam::GeometricField<Type, PatchField, GeoMesh>::ref()
{
    this->setUpToDate();
    storeOldTimes();
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal::FieldType&
Foam::GeometricField<Type, PatchField, GeoMesh>::primitiveFieldRef()
{
    this->setUpToDate();
    storeOldTimes();
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary&
Foam::GeometricField<Type, PatchField, GeoMesh>::boundaryFieldRef()
{
    this->setUpToDate();
    storeOldTimes();
    return boundaryField_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storeOldTimes() const
{
    if
    (
        field0Ptr_
     && timeIndex_ != this->time().timeIndex()
     && !(
            this->name().size() > 2
         && this->name()(this->name().size()-2, 2) == "_0"
         )
    )
    {
        storeOldTime();
    }

    // Correct time index
    timeIndex_ = this->time().timeIndex();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storeOldTime() const
{
    if (field0Ptr_)
    {
        field0Ptr_->storeOldTime();

        if (debug)
        {
            InfoInFunction
                << "Storing old time field for field" << endl
                << this->info() << endl;
        }

        *field0Ptr_ == *this;
        field0Ptr_->timeIndex_ = timeIndex_;

        if (field0Ptr_->field0Ptr_)
        {
            field0Ptr_->writeOpt() = this->writeOpt();
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::label Foam::GeometricField<Type, PatchField, GeoMesh>::nOldTimes() const
{
    if (field0Ptr_)
    {
        return field0Ptr_->nOldTimes() + 1;
    }
    else
    {
        return 0;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime() const
{
    if (!field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + "_0",
                this->time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                this->registerObject()
            ),
            *this
        );
    }
    else
    {
        storeOldTimes();
    }

    return *field0Ptr_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime()
{
    static_cast<const GeometricField<Type, PatchField, GeoMesh>&>(*this)
        .oldTime();

    return *field0Ptr_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storePrevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        if (debug)
        {
            InfoInFunction
                << "Allocating previous iteration field" << endl
                << this->info() << endl;
        }

        fieldPrevIterPtr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            this->name() + "PrevIter",
            *this
        );
    }
    else
    {
        *fieldPrevIterPtr_ == *this;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::prevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        FatalErrorInFunction
            << "previous iteration field" << endl << this->info() << endl
            << "  not stored."
            << "  Use field.storePrevIter() at start of iteration."
            << abort(FatalError);
    }

    return *fieldPrevIterPtr_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::
correctBoundaryConditions()
{
    this->setUpToDate();
    storeOldTimes();
    boundaryField_.evaluate();
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::needReference() const
{
    // Search all boundary conditions, if any are
    // fixed-value or mixed (Robin) do not set reference level for solution.

    bool needRef = true;

    forAll(boundaryField_, patchi)
    {
        if (boundaryField_[patchi].fixesValue())
        {
            needRef = false;
            break;
        }
    }

    reduce(needRef, andOp<bool>());

    return needRef;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::relax(const scalar alpha)
{
    if (debug)
    {
        InfoInFunction
           << "Relaxing" << endl << this->info() << " by " << alpha << endl;
    }

    operator==(prevIter() + alpha*(*this - prevIter()));
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::relax()
{
    word name = this->name();

    if
    (
        this->mesh().data::template lookupOrDefault<bool>
        (
            "finalIteration",
            false
        )
    )
    {
        name += "Final";
    }

    if (this->mesh().relaxField(name))
    {
        relax(this->mesh().fieldRelaxationFactor(name));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::word Foam::GeometricField<Type, PatchField, GeoMesh>::select
(
    bool final
) const
{
    if (final)
    {
        return this->name() + "Final";
    }
    else
    {
        return this->name();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::writeMinMax
(
    Ostream& os
) const
{
    os  << "min/max(" << this->name() << ") = "
        << Foam::min(*this).value() << ", "
        << Foam::max(*this).value()
        << endl;
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::
writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::T() const
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> result
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + ".T()",
                this->instance(),
                this->db()
            ),
            this->mesh(),
            this->dimensions()
        )
    );

    Foam::T(result.ref().primitiveFieldRef(), primitiveField());
    Foam::T(result.ref().boundaryFieldRef(), boundaryField());

    return result;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
Foam::GeometricField<Type, PatchField, GeoMesh>::component
(
    const direction d
) const
{
    tmp<GeometricField<cmptType, PatchField, GeoMesh>> Component
    (
        new GeometricField<cmptType, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + ".component(" + Foam::name(d) + ')',
                this->instance(),
                this->db()
            ),
            this->mesh(),
            this->dimensions()
        )
    );

    Foam::component(Component.ref().primitiveFieldRef(), primitiveField(), d);
    Foam::component(Component.ref().boundaryFieldRef(), boundaryField(), d);

    return Component;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::replace
(
    const direction d,
    const GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
     >& gcf
)
{
    primitiveFieldRef().replace(d, gcf.primitiveField());
    boundaryFieldRef().replace(d, gcf.boundaryField());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::replace
(
    const direction d,
    const dimensioned<cmptType>& ds
)
{
    primitiveFieldRef().replace(d, ds.value());
    boundaryFieldRef().replace(d, ds.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::max
(
    const dimensioned<Type>& dt
)
{
    Foam::max(primitiveFieldRef(), primitiveField(), dt.value());
    Foam::max(boundaryFieldRef(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::min
(
    const dimensioned<Type>& dt
)
{
    Foam::min(primitiveFieldRef(), primitiveField(), dt.value());
    Foam::min(boundaryFieldRef(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::negate()
{
    primitiveFieldRef().negate();
    boundaryFieldRef().negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    if (this == &gf)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, gf, "=");

    // Only assign field contents not ID

    ref() = gf();
    boundaryFieldRef() = gf.boundaryField();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    if (this == &(tgf()))
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    const GeometricField<Type, PatchField, GeoMesh>& gf = tgf();

    checkField(*this, gf, "=");

    // Only assign field contents not ID

    this->dimensions() = gf.dimensions();

    // Transfer the storage from the tmp
    primitiveFieldRef().transfer
    (
        const_cast<Field<Type>&>(gf.primitiveField())
    );

    boundaryFieldRef() = gf.boundaryField();

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const dimensioned<Type>& dt
)
{
    ref() = dt;
    boundaryFieldRef() = dt.value();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator==
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    const GeometricField<Type, PatchField, GeoMesh>& gf = tgf();

    checkField(*this, gf, "==");

    // Only assign field contents not ID

    ref() = gf();
    boundaryFieldRef() == gf.boundaryField();

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator==
(
    const dimensioned<Type>& dt
)
{
    ref() = dt;
    boundaryFieldRef() == dt.value();
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op              \
(                                                                              \
    const GeometricField<TYPE, PatchField, GeoMesh>& gf                        \
)                                                                              \
{                                                                              \
    checkField(*this, gf, #op);                                                \
                                                                               \
    ref() op gf();            \
    boundaryFieldRef() op gf.boundaryField();                                  \
}                                                                              \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op              \
(                                                                              \
    const tmp<GeometricField<TYPE, PatchField, GeoMesh>>& tgf                  \
)                                                                              \
{                                                                              \
    operator op(tgf());                                                        \
    tgf.clear();                                                               \
}                                                                              \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op              \
(                                                                              \
    const dimensioned<TYPE>& dt                                                \
)                                                                              \
{                                                                              \
    ref() op dt;                                       \
    boundaryFieldRef() op dt.value();                                          \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    gf().writeData(os, "internalField");
    os  << nl;
    gf.boundaryField().writeEntry("boundaryField", os);

    // Check state of IOstream
    os.check
    (
        "Ostream& operator<<(Ostream&, "
        "const GeometricField<Type, PatchField, GeoMesh>&)"
    );

    return (os);
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    os << tgf();
    tgf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef checkField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GeometricBoundaryField.C"
#include "GeometricFieldFunctions.C"

// ************************************************************************* //
