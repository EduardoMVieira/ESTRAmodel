/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025 UFES
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

Class
    Foam::constitutiveEqs::ESTRALog

SourceFiles
    ESTRALog.C

\*---------------------------------------------------------------------------*/

#ifndef ESTRALog_H
#define ESTRALog_H

#include "ESTRA.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{

/*---------------------------------------------------------------------------*\
                           Class ESTRALog Declaration
\*---------------------------------------------------------------------------*/

class ESTRALog
:
    public ESTRA
{
    // Private Data

        //- Logarithm of the conformation tensor
        volSymmTensorField theta_;

        //- Eigenvalues of theta diagonalization (values on the tensor diagonal)
        volTensorField eigVals_;

        //- Eigenvectors of theta diagonalization
        volTensorField eigVecs_;


public:

    //- Runtime type information
    TypeName("ESTRALog");


    // Constructors

        //- Construct from components
        ESTRALog
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construct
        ESTRALog(const ESTRALog&) = delete;


    //- Destructor

        virtual ~ESTRALog();


    // Member Functions

        //- Correct the viscoelastic stress: alpha is the color function
        //  in two phase-flows.
        virtual void correct
        (
            const volScalarField* alpha = nullptr,
            const volTensorField* gradU = nullptr
        ) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ESTRALog&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constitutiveEqs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
