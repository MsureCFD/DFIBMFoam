/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "IBM.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBM, 0);

    IBM::IBM
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
    )
    :
        
        dx_(dx),
        dy_(dy), 
        dz_(dz),
        nx_(nx),
        ny_(ny), 
        nz_(nz),  
        minX_(minX),
        minY_(minY), 
        minZ_(minZ),
        maxX_(maxX),
        maxY_(maxY), 
        maxZ_(maxZ),
        mesh_(mesh),
        objNum_(objNum)
    {

        Info << "The IBM is applied" << endl;
        

    }

    IBM::~IBM() {}

    void IBM::initFish(scalar &cTime, label &index, scalar &l, label &n, 
                      scalar &lambda, scalar &waveT, scalar &rad, scalar &cycleT)
    {
        cTime_ = cTime;
        index_ = index;
        length_ = l;
        nSection_ = n;
        lambda_ = lambda;
        waveT_ = waveT;
        rad_ = rad;
        cycleT_ = cycleT;
        Info << "current time is "<< cTime_ << endl;
        Info << "the index of current fish is "<< index_ << endl;
        Info << "the length of the fish is "<< length_ << endl;
        Info << "the sections number of the fish is "<< nSection_ << endl;
        Info << "the lambda of the fish is "<< lambda_ << endl;
        Info << "the wave period of the fish is "<< waveT_ << endl;
        Info << "the cycle radius of the fish is "<< rad_ << endl;
        Info << "the cycle period of the fish is "<< cycleT_ << endl;
        Info << "      "  << endl; 

        lastIbpCoor_.setSize(nSection_*2);
        currentIbpCoor_.setSize(nSection_*2);  
        interpolatedIbpVel_.setSize(nSection_*2);
        desiredIbpVel_.setSize(nSection_*2);
        IbpForce_.setSize(nSection_*2);
        IbpVolume_.setSize(nSection_*2);
        IbpNorm_.setSize(nSection_*2);
        IbpDs_.setSize(nSection_*2);
        IbpPressure_.setSize(nSection_*2);
        IbpstressTensor_.setSize(nSection_*2);

        /*offSet_.setSize(7);
        offSet_[0] = -3;
        offSet_[1] = -2;
        offSet_[2] = -1;
        offSet_[3] = 0;
        offSet_[4] = 1;
        offSet_[5] = 2;
        offSet_[6] = 3;*/
        
        offSet_.setSize(5);
        offSet_[0] = -2;
        offSet_[1] = -1;
        offSet_[2] = -0;
        offSet_[3] = 1;
        offSet_[4] = 2;

        /*offSet_.setSize(3);
        offSet_[0] = -1;
        offSet_[1] = 0;
        offSet_[2] = 1;*/




        updateIbpCoordinate(cTime);
        
    }

    void IBM::updateIbpCoordinate(scalar &time)
    {
        lastIbpCoor_ = currentIbpCoor_;
        
        scalar dL = length_/nSection_;
        scalar initAngle = index_*2.0*3.14/objNum_;
        scalar angle = initAngle + 2.0*3.14*time/cycleT_;
        scalar center_x = rad_ * Foam::cos(angle);
        scalar center_y = rad_ * Foam::sin(angle);
        scalar rotation_angle = angle - 3.14/2.0;
        for (label i = 0; i < nSection_; i++)
        {  
            
            scalar cx1 = i*dL;
            scalar cx2 = (i+1)*dL;

            scalar h1 = length_*(0.351*Foam::sin(cx1/length_ - 1.796) + 0.359)*Foam::sin(2.0*3.14*(cx1/lambda_ - time/waveT_));
            scalar h2 = length_*(0.351*Foam::sin(cx2/length_ - 1.796) + 0.359)*Foam::sin(2.0*3.14*(cx2/lambda_ - time/waveT_));

            scalar d1 = 0.125*length_ / 0.2 * (0.2969 * Foam::sqrt(cx1 / length_) - 0.1260 * (cx1 / length_) 
                      - 0.3516 * Foam::pow(cx1 / length_, 2.0) + 0.2843 * Foam::pow(cx1 / length_, 3.0) - 0.1015 * Foam::pow(cx1 / length_, 4.0));
            scalar d2 = 0.125*length_ / 0.2 * (0.2969 * Foam::sqrt(cx2 / length_) - 0.1260 * (cx2 / length_) 
                      - 0.3516 * Foam::pow(cx2 / length_, 2.0) + 0.2843 * Foam::pow(cx2 / length_, 3.0) - 0.1015 * Foam::pow(cx2 / length_, 4.0));

            scalar fish_x1_unrotated_upper = cx1;
            scalar fish_y1_unrotated_upper = h1 + d1;
            scalar fish_x1_unrotated_lower = cx1;
            scalar fish_y1_unrotated_lower = h1 - d1;

            scalar fish_x2_unrotated_upper = cx2;
            scalar fish_y2_unrotated_upper = h2 + d2;
            scalar fish_x2_unrotated_lower = cx2;
            scalar fish_y2_unrotated_lower = h2 - d2;

            scalar fish_x1_upper = center_x + Foam::cos(rotation_angle) * fish_x1_unrotated_upper 
                                            - Foam::sin(rotation_angle) * fish_y1_unrotated_upper;
            scalar fish_y1_upper = center_y + Foam::sin(rotation_angle) * fish_x1_unrotated_upper 
                                            + Foam::cos(rotation_angle) * fish_y1_unrotated_upper;
            scalar fish_x1_lower = center_x + Foam::cos(rotation_angle) * fish_x1_unrotated_lower 
                                            - Foam::sin(rotation_angle) * fish_y1_unrotated_lower;
            scalar fish_y1_lower = center_y + Foam::sin(rotation_angle) * fish_x1_unrotated_lower 
                                            + Foam::cos(rotation_angle) * fish_y1_unrotated_lower;

            scalar fish_x2_upper = center_x + Foam::cos(rotation_angle) * fish_x2_unrotated_upper 
                                            - Foam::sin(rotation_angle) * fish_y2_unrotated_upper;
            scalar fish_y2_upper = center_y + Foam::sin(rotation_angle) * fish_x2_unrotated_upper 
                                            + Foam::cos(rotation_angle) * fish_y2_unrotated_upper;
            scalar fish_x2_lower = center_x + Foam::cos(rotation_angle) * fish_x2_unrotated_lower 
                                            - Foam::sin(rotation_angle) * fish_y2_unrotated_lower;
            scalar fish_y2_lower = center_y + Foam::sin(rotation_angle) * fish_x2_unrotated_lower 
                                            + Foam::cos(rotation_angle) * fish_y2_unrotated_lower;

            currentIbpCoor_[i*2].x() = (fish_x1_upper + fish_x2_upper) / 2.0;
            currentIbpCoor_[i*2].y() = (fish_y1_upper + fish_y2_upper) / 2.0;
            currentIbpCoor_[i*2].z() = 0.0;
            currentIbpCoor_[i*2 + 1].x() = (fish_x1_lower + fish_x2_lower) / 2.0;
            currentIbpCoor_[i*2 + 1].y() = (fish_y1_lower + fish_y2_lower) / 2.0;
            currentIbpCoor_[i*2 + 1].z() = 0.0;

            IbpDs_[i*2] = Foam::sqrt((fish_x2_upper - fish_x1_upper) * (fish_x2_upper - fish_x1_upper) 
                        + (fish_y2_upper - fish_y1_upper) * (fish_y2_upper - fish_y1_upper));
            IbpVolume_[i*2] = IbpDs_[i*2] * Foam::sqrt(dx_*dy_);
            IbpDs_[i*2 + 1] = Foam::sqrt((fish_x2_lower - fish_x1_lower) * (fish_x2_lower - fish_x1_lower) 
                        + (fish_y2_lower - fish_y1_lower) * (fish_y2_lower - fish_y1_lower));
            IbpVolume_[i*2+1] = IbpDs_[i*2 + 1] * Foam::sqrt(dx_*dy_);

            vector t = vector::zero;
            t.x() = (fish_x2_upper - fish_x1_upper) / IbpDs_[i*2];
            t.y() = (fish_y2_upper - fish_y1_upper) / IbpDs_[i*2];
            IbpNorm_[i*2].x() = - t.y();
            IbpNorm_[i*2].y() = t.x();
            IbpNorm_[i*2].z() = 0.0;

            t.x() = (fish_x2_lower - fish_x1_lower) / IbpDs_[i*2 + 1];
            t.y() = (fish_y2_lower - fish_y1_lower) / IbpDs_[i*2 + 1];
            IbpNorm_[i*2 + 1].x() = t.y();
            IbpNorm_[i*2 + 1].y() = -t.x();
            IbpNorm_[i*2 + 1].z() = 0.0;


        }

        

    }

    void IBM::updateIbpVelocity(scalar &dt)
    {
        //- calculate the desired velocity for each ib point
        forAll(currentIbpCoor_, I)
        {
            
            desiredIbpVel_[I] = (currentIbpCoor_[I] - lastIbpCoor_[I])/dt;
            
        }

        //Info << desiredIbpVel_ << endl;
    }

    void IBM::writeVTK(label num, label timeStepIndex)
    {
        std::string folder = "output";
        std::ostringstream filename;
        filename << folder << "/fish" << num << "-" << std::setw(4) << std::setfill('0') << timeStepIndex << ".vtk";
        std::ofstream vtkFile(filename.str());


        vtkFile << "# vtk DataFile Version 3.0" << std::endl;
        vtkFile << "Fish shape at time step " << timeStepIndex << std::endl;
        vtkFile << "ASCII" << std::endl;
        vtkFile << "DATASET POLYDATA" << std::endl;


        vtkFile << "POINTS " << currentIbpCoor_.size() << " float" << std::endl;
        for (const auto& pt : currentIbpCoor_) {
            vtkFile << pt.x() << " " << pt.y() << " " << pt.z() << std::endl;
        }

        int numQuads = nSection_ - 1; 
        int numTriangles = numQuads * 2; 
        vtkFile << "POLYGONS " << numTriangles << " " << numTriangles * 4 << std::endl;
        for (label i = 0; i < numQuads; i++) {

            int idx1 = 2 * i;       
            int idx2 = 2 * i + 1;   
            int idx3 = 2 * i + 2;   
            int idx4 = 2 * i + 3;   

            vtkFile << "3 " << idx1 << " " << idx2 << " " << idx3 << std::endl;

            vtkFile << "3 " << idx2 << " " << idx4 << " " << idx3 << std::endl;
        }

        vtkFile.close();
    }

    void IBM::writeForce(label index, scalar &time, volScalarField &p, volTensorField &stressTensor)
    {   
        

        forAll(currentIbpCoor_, I)
        {
            scalar d = Foam::sqrt(dx_*dy_);
            vector pos1 = currentIbpCoor_[I] + 1*d*IbpNorm_[I];
            vector pos2 = currentIbpCoor_[I] + 2*d*IbpNorm_[I];
            label cellId1 = findCell_(pos1);
            label cellId2 = findCell_(pos2);
            IbpPressure_[I] = (2.0*d*p[cellId1] + 1*d*p[cellId2])/(1*d+2.0*d);
            IbpstressTensor_[I] = (2.0*d*stressTensor[cellId1] + 1*d*stressTensor[cellId2])/(1*d+2.0*d);
        }
        
        
        vector pressureForce = vector::zero;
        vector shearForce = vector::zero;
        vector totalForce = vector::zero;
        forAll(currentIbpCoor_, I)
        {
            pressureForce += IbpPressure_[I] * IbpNorm_[I] * IbpDs_[I];
            shearForce +=  (IbpstressTensor_[I] & IbpNorm_[I]) * IbpDs_[I];
        }
        pressureForce = pressureForce * 2.0/(5.0/120.0)/(5.0/120.0);
        shearForce = shearForce * 2.0/(5.0/120.0)/(5.0/120.0);
        totalForce = pressureForce + shearForce;
        
        std::string folder = "force";

         
        std::string filename = folder + "/fish" + std::to_string(index) + ".txt";
        std::ofstream ofs(filename, std::ios::app);
        if (!ofs.is_open())
        {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }
        ofs << std::fixed << std::setprecision(6);

        ofs << time << "  " << pressureForce.x() << "  " << pressureForce.y() << "  " << pressureForce.z() 
                   << "  " << shearForce.x() << "  " << shearForce.y() << "  " << shearForce.z() 
                   << "  " << totalForce.x() << "  " << totalForce.y() << "  " << totalForce.z() << std::endl;
        ofs.close();
    }

    label IBM::findCell_(vector &IbpCoordinate)
    {
        
        double epsilon = 1e-9;
        double tx = IbpCoordinate.x() - minX_ - epsilon;
        double ty = IbpCoordinate.y() - minY_ - epsilon;
        double tz = IbpCoordinate.z() - minZ_ - epsilon;

        if (tx < 0.0) tx = 0.0;
        if (ty < 0.0) ty = 0.0;
        if (tz < 0.0) tz = 0.0;

        int i = static_cast<int>(tx / dx_);
        int j = static_cast<int>(ty / dy_);
        int k = static_cast<int>(tz / dz_);

        int customIndex = i + nx_ * (j + ny_ * k);

        return customIndex;

    }

    void IBM::calcDelta
    (
    )
    {
        IbpStencil_.clear();
        IbpDelta_.clear();
        IbpDeltaSum_.clear();

        scalar rx;
        scalar ry;
        //scalar rz;
        scalar delta;
        scalar totDelta;
        

        label cellId;

        vector tmpPoint;
        vector nearestCellCoordinate;
        vector r;

        List<scalar> tmpDelta;
        List<label> tmpStencil;
        
        forAll(currentIbpCoor_, i)
        {
            /*xx = currentIbpCoor_[i].x();
            yy = currentIbpCoor_[i].y();
            zz = currentIbpCoor_[i].z();*/

            
            cellId = findCell_(currentIbpCoor_[i]);
            nearestCellCoordinate = mesh_.C()[cellId];

            totDelta = 0.0;
            tmpStencil.clear();
            tmpDelta.clear();

            //forAll(offSet_, zi)
            {
                //tmpPoint.z() = nearestCellCoordinate.z() + offSet_[zi]*dz_;
               // if ((tmpPoint.z() > minZ_) && (tmpPoint.z() < maxZ_))
                tmpPoint.z() = 0.0;
                {
                    forAll(offSet_, yi)
                    {
                        tmpPoint.y() = nearestCellCoordinate.y() + offSet_[yi]*dy_;
                        if ((tmpPoint.y() > minY_) && (tmpPoint.y() < maxY_))
                        {
                            forAll(offSet_, xi)
                            {
                                tmpPoint.x() = nearestCellCoordinate.x() + offSet_[xi]*dx_;
                                if ((tmpPoint.x() > minX_) && (tmpPoint.x() < maxX_))
                                {

                                    cellId = findCell_(tmpPoint);
                                    r = cmptMag(currentIbpCoor_[i] - mesh_.C()[cellId]);
                                    rx = r.x()/dx_;
                                    ry = r.y()/dy_;
                                    //rz = r.z()/dz_;
                                    //delta = deltaFunction_(rx)*deltaFunction_(ry)*deltaFunction_(rz);
                                    delta = deltaFunction_(rx)*deltaFunction_(ry);

                                    if (delta > 0.0)
                                    {
                                        tmpStencil.append(cellId);
                                        tmpDelta.append(delta);
                                        totDelta += delta;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            IbpStencil_.append(tmpStencil);
            IbpDelta_.append(tmpDelta);
            IbpDeltaSum_.append(totDelta);

            //Info << "I = " << i << "    "
            //     << "tmpStencil = " << tmpStencil 
            //     << "tmpDelta = " << tmpDelta<< endl;
        }

        //Info << IBpointsDeltaSum_ << endl;
        //Info << IbpDeltaSum_ << endl;

    }

    void IBM::interpolateVelocity
    (
        vectorField &tU
    )
    {
        forAll(currentIbpCoor_, I)
        {

            interpolatedIbpVel_[I] = vector::zero;

        }
        
        forAll(currentIbpCoor_, I)
        {
            List<scalar> tmpDelta;
            List<label> tmpStencil;
            tmpStencil = IbpStencil_[I];
            tmpDelta = IbpDelta_[I];
            forAll(tmpStencil, J)
            {
                interpolatedIbpVel_[I] += tmpDelta[J]*tU[tmpStencil[J]];
            }    
        }
    }

    void IBM::distributeForce
    (
        vectorField &f,
        scalar &dt
    )
    {
        forAll(currentIbpCoor_, I)
        {
            vector velocityDiff = 
                desiredIbpVel_[I] - interpolatedIbpVel_[I];
                
            List<scalar> tmpDelta;
            List<label> tmpStencil;
            tmpStencil = IbpStencil_[I];
            tmpDelta = IbpDelta_[I];
            forAll(tmpStencil, J)
            {
                IbpForce_[I] -= tmpDelta[J]*velocityDiff/dt*IbpVolume_[I];
                //f[tmpStencil[J]] += tmpDelta[J]*velocityDiff/dt/(dx_*dy_)*IBpointsVolume_[I];
                f[tmpStencil[J]] += tmpDelta[J]*velocityDiff/dt/(dx_*dy_)*IbpVolume_[I];
            }

            
        }
    }

    scalar IBM::deltaFunction_(const scalar &distance)
    {
        if (distance < 1.0)
        {
            return 0.125 * (3.0 - 2.0 * distance + 
                   sqrt(1.0 + 4.0 * distance - 4 * distance * distance));
        }
        else if (distance < 2.0)
        {
            return 0.125 * (5.0 - 2.0 * distance - sqrt(-7.0 + 12.0 * distance - 
                   4 * distance * distance));
        }

        return 0.0;    
    }




}
