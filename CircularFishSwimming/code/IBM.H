/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Description
    Immersed boundary method.

SourceFiles
    IBM.C

\*---------------------------------------------------------------------------*/

#ifndef IBM_H
#define IBM_H

#include "fvCFD.H"

#include <fstream>
#include <iomanip>
#include <sstream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IBM Declaration
\*---------------------------------------------------------------------------*/

class IBM
{
public:

    // Public data types

        //- Runtime type information
        TypeName("IBM");

        //- Constructors

        IBM   
        (
        const fvMesh &mesh,
        const scalar &dx, 
        const scalar &dy,
        const scalar &dz,
        const label &nx, 
        const label &ny,
        const label &nz,
        const scalar &minX, 
        const scalar &minY,
        const scalar &minZ,
        const scalar &maxX, 
        const scalar &maxY,
        const scalar &maxZ,
        const label &objNum
        );

        //- Destructor
        virtual ~IBM();
        
        void initFish(scalar &cTime, label &index, scalar &l, label &n, 
                      scalar &lambda, scalar &waveT, scalar &rad, scalar &cycleT);
        
        //void updateIbpCoordinates(List<vector>& coorList, scalar time);
        void updateIbpCoordinate(scalar &time);
        void updateIbpVelocity(scalar &dt);
        void interpolateVelocity(vectorField &tU);
        void distributeForce(vectorField &f, scalar &dt);

        void writeVTK(label num, label timeStepIndex);
        void writeForce(label index, scalar &time, volScalarField &p, volTensorField &stressTensor);

        void calcDelta();

private:

    // Private Data

        
        //- The fluid cell size in x-direction
        const scalar dx_;
        
        //- The fluid cell size in y-direction
        const scalar dy_;
        
        //- The fluid cell size in z-direction
        const scalar dz_;

        //- The fluid cell size in x-direction
        const scalar nx_;
        
        //- The fluid cell size in y-direction
        const scalar ny_;
        
        //- The fluid cell size in z-direction
        const scalar nz_;

        const scalar minX_;
        
        const scalar minY_;
        
        const scalar minZ_;

        const scalar maxX_;
        
        const scalar maxY_;
        
        const scalar maxZ_;

        const fvMesh& mesh_;

        const label objNum_;

        scalar cTime_; // current simulation time

        //- parameter for fish
        label index_; // index of each fish
        scalar length_; // lengrh of one fish
        label nSection_; // number of divisions of one fish
        scalar lambda_; // wave length
        scalar waveT_; // wave period
        scalar rad_; // radius of cycling
        scalar cycleT_;

        List<vector> lastIbpCoor_;
        List<vector> currentIbpCoor_;
        List<vector> interpolatedIbpVel_;
        List<vector> desiredIbpVel_;
        List<vector> IbpForce_;
        List<scalar> IbpDs_;
        List<scalar> IbpPressure_;
        List<tensor> IbpstressTensor_;
        List<scalar> IbpVolume_;
        List<vector> IbpNorm_; 
        
        // list used for recording the stencil and delta for each ib point to avoid repeated mesh search
        labelListList IbpStencil_;
        scalarListList IbpDelta_;
        List<scalar> IbpDeltaSum_;


        List<label> offSet_;
        //- Calculate the delta function coefficient in each direction
        scalar deltaFunction_(const scalar &distance);

        label findCell_(vector &IbpCoordinate);


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
