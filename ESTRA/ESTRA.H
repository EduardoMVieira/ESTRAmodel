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
    Foam::constitutiveEqs::ESTRA

Description
    ESTRA model.

SourceFiles
    ESTRA.C

\*---------------------------------------------------------------------------*/

#ifndef ESTRA_H
#define ESTRA_H

#include "constitutiveEq.H"
#include "thermoFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{

/*---------------------------------------------------------------------------*\
                            Class ESTRA Declaration
\*---------------------------------------------------------------------------*/

class ESTRA
:
    public constitutiveEq
{
    // Private Data

    


protected:

    // Protected Data

        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity
            dimensionedScalar etaP_;

            //- Relaxation time
            dimensionedScalar lambda_;

            //- Initial molecular adimensional extensibility
            dimensionedScalar L2_0_;

        // Convected parameters

            //- Molecular adimensional extensibility field
            volScalarField L2_;

            //- Gamma
            volScalarField gamma_;

            //- Transition to softening regime stress
            volScalarField tau0_;

            //- First difference of normal stresses
            volScalarField N1_;

        //- Transported viscoelastic stress
        volSymmTensorField tau_;

        //- Conformation tensor
        volSymmTensorField A_;

        //- varF_ = 1/f for tau models, and varF_ = f otherwise.
        volScalarField varF_;

        //- Actual iteration
        label iteration_;

        //- Number of iterations after which test the breaking
        label breakingInterval_;

        //- Thermo function for temperature dependence of lambda
        autoPtr<thermoFunction> thermoLambdaPtr_;

        //- Thermo function for temperature dependence of viscosity
        autoPtr<thermoFunction> thermoEtaPtr_;


    // Protected Member Functions

        //- Return the solvent viscosity
        virtual const dimensionedScalar etaS() const override
        {
            return etaS_;
        }

        //- Return the polymeric viscosity
        virtual const dimensionedScalar etaP() const override
        {
            return etaP_;
        }

        //- Return etaS corrected for temperature
        virtual tmp<volScalarField> etaSThermo() const override
        {
            return thermoEtaPtr_->createField(etaS_);
        }

        //- Return etaP corrected for temperature
        virtual tmp<volScalarField> etaPThermo() const override
        {
            return thermoEtaPtr_->createField(etaP_);
        }

        //- Is the model prepared to work in non-isothermal conditions
        virtual bool hasThermo() const override
        {
            return true;
        }


public:

    //- Runtime type information
    TypeName("ESTRA");


    // Constructors

        //- Construct from components
        ESTRA
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construct
        ESTRA(const ESTRA&) = delete;


    //- Destructor

        virtual ~ESTRA();


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const override
        {
            return tau_;
        }

        //- Return the density
        virtual const dimensionedScalar rho() const override
        {
            return rho_;
        }

        //- Is model GNF (non-elastic)
        virtual bool isGNF() const override
        {
            return false;
        };

        //- Return gamma
        tmp<volScalarField> gamma() const
        {
            return gamma_;
        }

        //- Correct the viscoelastic stress: alpha is the color function
        //  in two phase-flows.
        virtual void correct
        (
            const volScalarField* alpha,
            const volTensorField* gradU
        ) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ESTRA&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constitutiveEqs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
