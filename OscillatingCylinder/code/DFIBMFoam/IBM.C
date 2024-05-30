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

    void IBM::initCircle(scalar &cTime, label &index, vector &centerCoordinate, label &nSection, scalar &diameter)
    {
        cTime_ = cTime;
        index_ = index;
        nSection_ = nSection;
        diameter_ = diameter;
        centerCoordinate_ = centerCoordinate;

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
        IbpDeltaSum_.setSize(nSection_*2);
        
        offSet_.setSize(5);
        offSet_[0] = -2;
        offSet_[1] = -1;
        offSet_[2] = -0;
        offSet_[3] = 1;
        offSet_[4] = 2;

        updateIbpCoordinate(cTime);
        
    }

    void IBM::updateIbpCoordinate(scalar &time)
    {
        lastIbpCoor_ = currentIbpCoor_;
        scalar dL = diameter_/nSection_;
        
        for (label i = 0; i < nSection_; i++)
        {  
            
            scalar cx1 = i*dL - diameter_/2.0;
            scalar cx2 = (i+1)*dL - diameter_/2.0;

            scalar cy1 = Foam::sqrt((0.5 * diameter_) * (0.5 * diameter_) - cx1 * cx1);
            scalar cy2 = Foam::sqrt((0.5 * diameter_) * (0.5 * diameter_) - cx2 * cx2);

            currentIbpCoor_[i*2].x() = 0.5 * (cx1 + cx2) + centerCoordinate_.x() - 5.0/2.0/3.14*Foam::sin(2.0*3.14*1.0/120.0*time);
            currentIbpCoor_[i*2].y() = 0.5 * (cy1 + cy2) + centerCoordinate_.y();
            currentIbpCoor_[i*2].z() = 0.0;
            currentIbpCoor_[i*2 + 1].x() = 0.5 * (cx1 + cx2) + centerCoordinate_.x() - 5.0/2.0/3.14*Foam::sin(2.0*3.14*1.0/120.0*time);
            currentIbpCoor_[i*2 + 1].y() = - 0.5 * (cy1 + cy2) + centerCoordinate_.y();
            currentIbpCoor_[i*2 + 1].z() = 0.0;

            IbpDs_[i*2] = Foam::sqrt((cx2 - cx1) * (cx2 - cx1) + (cy2 - cy1) * (cy2 - cy1));
            IbpVolume_[i*2] = IbpDs_[i*2] * Foam::sqrt(dx_*dy_);
            IbpDs_[i*2 + 1] = Foam::sqrt((cx2 - cx1) * (cx2 - cx1) + (cy2 - cy1) * (cy2 - cy1));
            IbpVolume_[i*2 + 1] = IbpVolume_[i*2];

            vector t = vector::zero;
            scalar magt = Foam::sqrt((cx2 - cx1) * (cx2 - cx1) + (cy2 - cy1) * (cy2 - cy1));
            t.x() = (cx2 - cx1) / magt;
            t.y() = (cy2 - cy1) / magt;
            IbpNorm_[i*2].x() = - t.y();
            IbpNorm_[i*2].y() = t.x();
            IbpNorm_[i*2].z() = 0.0;

            t.x() = (cx2 - cx1) / magt;
            t.y() = (-cy2 + cy1) / magt;
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
            //desiredIbpVel_[I] = vector::zero;
            //desiredIbpVel_[I].x() = -5.0/12.0*Foam::cos(2.0*3.14*1.0/12.0*time);
            
        }


        //Info << desiredIbpVel_ << endl;
    }

    void IBM::writeVTK(label num, label timeStepIndex)
    {
        std::ostringstream filename;
        filename << "output/fish" << num << "-" << std::setw(4) << std::setfill('0') << timeStepIndex << ".vtk";
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
        
        /*forAll(currentIbpCoor_, I)
        {
            tensor zeroTensor(0,0,0,0,0,0,0,0,0);
            IbpstressTensor_[I] = zeroTensor;
            IbpPressure_[I] = 0.0;

            List<scalar> tmpDelta;
            List<label> tmpStencil;
            List<label> tmpOutside;
            tmpStencil = IbpStencil_[I];
            tmpDelta = IbpDelta_[I];
            tmpOutside = IbpOutside_[I];
            forAll(tmpStencil, J)
            {
                if (tmpOutside[J] == 1)
                {
                    IbpPressure_[I] += tmpDelta[J]*p[tmpStencil[J]]/IbpDeltaSum_[I];
                    IbpstressTensor_[I] += tmpDelta[J]*stressTensor[tmpStencil[J]]/IbpDeltaSum_[I];
                    //IbpPressure_[I] += tmpDelta[J]*p[tmpStencil[J]];
                    //IbpstressTensor_[I] += tmpDelta[J]*stressTensor[tmpStencil[J]];
                }
                
                //interpolatedIbpVel_[I] += tmpDelta[J]*tU[tmpStencil[J]];
            }

        }*/

        /*forAll(currentIbpCoor_, I)
        {
            scalar d1 = Foam::sqrt(dx_*dy_);
            vector pos1 = currentIbpCoor_[I] + d1*IbpNorm_[I];
            label cellId1 = findCell_(pos1);

            label cellId2 = findCell_(currentIbpCoor_[I]);
            scalar d2 = Foam::mag(mesh_.C()[cellId2] - currentIbpCoor_[I]);
            
            IbpPressure_[I] = (d2*p[cellId1] + d1*p[cellId2])/(d1+d2);
            IbpstressTensor_[I] = (d2*stressTensor[cellId1] + d1*stressTensor[cellId2])/(d1+d2);
        }*/
       /*forAll(currentIbpCoor_, I)
        {
            scalar d = Foam::sqrt(dx_*dy_);
            vector pos1 = currentIbpCoor_[I] + 1*d*IbpNorm_[I];
            label cellId1 = findCell_(pos1);
            IbpPressure_[I] = p[cellId1];
            IbpstressTensor_[I] = stressTensor[cellId1];
        }*/
        /*forAll(currentIbpCoor_, I)
        {
            scalar d = Foam::sqrt(dx_*dy_);
            vector pos1 = currentIbpCoor_[I] + 1*d*IbpNorm_[I];
            vector pos2 = currentIbpCoor_[I] + 2*d*IbpNorm_[I];
            vector pos3 = currentIbpCoor_[I] + 3*d*IbpNorm_[I];
            label cellId1 = findCell_(pos1);
            label cellId2 = findCell_(pos2);
            label cellId3 = findCell_(pos3);
            IbpPressure_[I] = (3.0*d*p[cellId1] + 2*d*p[cellId2] +1*d*p[cellId3])/(1*d + 2.0*d + 3*d);
            IbpstressTensor_[I] = (3.0*d*stressTensor[cellId1] + 2*d*stressTensor[cellId2] +1*d*stressTensor[cellId3])/(1*d + 2.0*d + 3*d);
        }*/

        /*forAll(currentIbpCoor_, I)
        {
            scalar d = Foam::sqrt(dx_*dy_);
            vector pos1 = currentIbpCoor_[I] + 1*d*IbpNorm_[I];
            vector pos2 = currentIbpCoor_[I] + 2*d*IbpNorm_[I];
            vector pos3 = currentIbpCoor_[I] + 3*d*IbpNorm_[I];
            vector pos4 = currentIbpCoor_[I] + 4*d*IbpNorm_[I];
            vector pos5 = currentIbpCoor_[I] + 5*d*IbpNorm_[I];
            label cellId1 = findCell_(pos1);
            label cellId2 = findCell_(pos2);
            label cellId3 = findCell_(pos3);
            label cellId4 = findCell_(pos4);
            label cellId5 = findCell_(pos5);
            IbpPressure_[I] = (5.0*d*p[cellId1] + 4*d*p[cellId2] +3*d*p[cellId3] + 2*d*p[cellId4] + 1*d*p[cellId5])/(1*d + 2.0*d + 3*d +4*d +5*d);
            IbpstressTensor_[I] = (5.0*d*stressTensor[cellId1] + 4*d*stressTensor[cellId2] 
                                + 3*d*stressTensor[cellId3] + 2*d*stressTensor[cellId4] + 1*d*stressTensor[cellId5])/(1*d + 2.0*d + 3*d +4*d +5*d);
        }*/

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
        
        std::string filename = std::to_string(index) + ".txt";
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
        IbpOutside_.clear();
        

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
        List<label> tmpOutside;
        
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
            tmpOutside.clear();

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
                                        vector n1 = mesh_.C()[cellId] - currentIbpCoor_[i];
                                        scalar tmpDot = n1 & IbpNorm_[i];
                                        if (tmpDot > 0)
                                        {
                                            tmpOutside.append(1); 
                                            tmpStencil.append(cellId);
                                            tmpDelta.append(delta);
                                            totDelta += delta;
                                        }
                                        else
                                        {
                                            tmpOutside.append(0); 
                                            tmpStencil.append(cellId);
                                            tmpDelta.append(delta);
                                        }
                                        
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }

            IbpStencil_.append(tmpStencil);
            IbpDelta_.append(tmpDelta);
            IbpOutside_.append(tmpOutside);
            IbpDeltaSum_[i] = totDelta;
            
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
        /*forAll(currentIbpCoor_, I)
        {
            scalar d = Foam::sqrt(dx_*dy_);
            vector pos1 = currentIbpCoor_[I] + d*IbpNorm_[I];
            vector pos2 = currentIbpCoor_[I] + 2*d*IbpNorm_[I];
            label cellId1 = findCell_(pos1);
            label cellId2 = findCell_(pos2);
            interpolatedIbpVel_[I] = (2.0*d*tU[cellId1] + d*tU[cellId2])/(d+2.0*d);
        }*/
        /*forAll(currentIbpCoor_, I)
        {
            scalar d = Foam::sqrt(dx_*dy_);
            vector pos1 = currentIbpCoor_[I] + d*IbpNorm_[I];
            vector pos2 = currentIbpCoor_[I] - d*IbpNorm_[I];
            label cellId1 = findCell_(pos1);
            label cellId2 = findCell_(pos2);
            interpolatedIbpVel_[I] = (tU[cellId1] + tU[cellId2])/(2.0);
        }*/

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
        /*forAll(currentIbpCoor_, I)
        {
            vector velocityDiff = desiredIbpVel_[I] - interpolatedIbpVel_[I];
            label cellId = findCell_(currentIbpCoor_[I]);
            f[cellId] = velocityDiff/dt/(dx_*dy_)*IbpVolume_[I];
        }*/
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
