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
    scalar omega = 2.0*Foam::constant::mathematical::pi*f_;

    // Convert beta to rad
    scalar beta = beta_*Foam::constant::mathematical::pi/180.0;
  
    // Translational displacement
    scalar x = h_/Foam::tan(beta) * Foam::cos(omega*t);
    scalar y = h_ * Foam::cos(omega*t) - h_;
    scalar z = 0.0;

    // Translaional time derivative
    scalar xd = -omega*h_ /Foam::tan(beta) * Foam::sin(omega*t);
    scalar yd = -omega*h_ * Foam::sin(omega*t);

    //AOA profile
    //scalar alpha = alpha_m_ * (0.5-0.5*Foam::cos(2.0*omega*t));
    scalar alpha = alpha_m_ * Foam::sin(omega*t);

    // Rotational displacement
    scalar theta_x = 0.0;
    scalar theta_y = 0.0;
    scalar theta_z = -(Foam::atan2(yd,(magU_ - xd)) + alpha*Foam::constant::mathematical::pi/180.0);
    theta_z *= 180.0/Foam::constant::mathematical::pi;

    //const vector displacement = H_amplitude_*Foam::sin(omega_*t + Foam::constant::mathematical::pi/2) - H_amplitude_;
    const vector displacement(x,y,z);

    //vector eulerAngles = P_amplitude_*Foam::sin(omega_*t);
    //vector eulerAngles = amplitude_*sin(omega_*t);
    vector eulerAngles(theta_x,theta_y,theta_z);

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

    SBMFCoeffs_.lookup("magU")>>magU_;
    SBMFCoeffs_.lookup("f")>>f_;
    SBMFCoeffs_.lookup("origin")>>origin_;
    SBMFCoeffs_.lookup("h")>>h_;
    SBMFCoeffs_.lookup("alpha_max")>>alpha_m_;
    SBMFCoeffs_.lookup("beta")>>beta_;

    return true;
}


// ************************************************************************* //
