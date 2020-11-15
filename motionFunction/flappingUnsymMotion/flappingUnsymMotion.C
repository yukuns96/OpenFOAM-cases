/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "flappingUnsymMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(flappingUnsymMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        flappingUnsymMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::flappingUnsymMotion::
flappingUnsymMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::flappingUnsymMotion::
transformation() const
{
    scalar t = time_.value();

    const vector displacement = H_amplitude_*Foam::sin(omega_*t + Foam::constant::mathematical::pi/2) - H_amplitude_;
    vector eulerAngles = P_amplitude_*Foam::sin(omega_*t);
    //vector eulerAngles = amplitude_*sin(omega_*t);

    // Convert the rotational motion from deg to rad
    eulerAngles *= degToRad();

    quaternion R(quaternion::XYZ, eulerAngles);
    
    septernion TR(septernion(-origin_ + -displacement)*R*septernion(origin_));
    //septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::flappingUnsymMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    /*
    SBMFCoeffs_.readEntry("origin", origin_);
    SBMFCoeffs_.readEntry("amplitude", amplitude_);
    SBMFCoeffs_.readEntry("omega", omega_);
    */

    SBMFCoeffs_.lookup("H_amplitude")>>H_amplitude_;
    SBMFCoeffs_.lookup("P_amplitude")>>P_amplitude_;
    SBMFCoeffs_.lookup("origin")>>origin_;
    SBMFCoeffs_.lookup("omega")>>omega_;

    return true;
}


// ************************************************************************* //