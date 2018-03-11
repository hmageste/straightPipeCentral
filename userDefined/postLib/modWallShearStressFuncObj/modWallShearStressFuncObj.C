/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "modWallShearStressFuncObj.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(modWallShearStressFuncObj, 0);
    addToRunTimeSelectionTable(functionObject, modWallShearStressFuncObj, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::modWallShearStressFuncObj::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall shear stress");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    file() << endl;
}


void Foam::functionObjects::modWallShearStressFuncObj::calcShearStress
(
    const volSymmTensorField& Reff,
    volVectorField& shearStress
)
{
    shearStress.dimensions().reset(Reff.dimensions());

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        vectorField& ssp = shearStress.boundaryFieldRef()[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        ssp = (-Sfp/magSfp) & Reffp;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::modWallShearStressFuncObj::modWallShearStressFuncObj
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_()
{
    volVectorField* modWallShearStressFuncObjPtr
    (
        new volVectorField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "0",
                sqr(dimLength)/sqr(dimTime),
                Zero
            )
        )
    );

    mesh_.objectRegistry::store(modWallShearStressFuncObjPtr);

    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::modWallShearStressFuncObj::~modWallShearStressFuncObj()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::modWallShearStressFuncObj::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall shear stress on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::modWallShearStressFuncObj::execute()
{
    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    volVectorField& modWallShearStressFuncObj =
        mesh_.lookupObjectRef<volVectorField>(type());

    tmp<volSymmTensorField> Reff;
    if (mesh_.foundObject<cmpModel>(turbulenceModel::propertiesName))
    {
        const cmpModel& model =
            mesh_.lookupObject<cmpModel>(turbulenceModel::propertiesName);

        Reff = model.devRhoReff();
    }
    else if (mesh_.foundObject<icoModel>(turbulenceModel::propertiesName))
    {
        const icoModel& model =
            mesh_.lookupObject<icoModel>(turbulenceModel::propertiesName);

        Reff = model.devReff();
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    calcShearStress(Reff(), modWallShearStressFuncObj);

    return true;
}


bool Foam::functionObjects::modWallShearStressFuncObj::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volVectorField& modWallShearStressFuncObj =
        obr_.lookupObject<volVectorField>(type());

    const fvPatchList& patches = mesh_.boundary();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const vectorField& ssp = modWallShearStressFuncObj.boundaryField()[patchi];

        vector minSsp = gMin(ssp);
        vector maxSsp = gMax(ssp);

        if (Pstream::master())
        {
            file() << mesh_.time().value()
                << token::TAB << pp.name()
                << token::TAB << minSsp
                << token::TAB << maxSsp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
