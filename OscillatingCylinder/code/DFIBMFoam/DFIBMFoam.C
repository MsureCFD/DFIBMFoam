/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    pimpleFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

#include <memory>  // 包含智能指针的头文件
#include "IBM.H"
#include <nlohmann/json.hpp>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Read the fluid grid parameter in blockMesh and calculate the grid delta (uniform grid)
    const label  nx(readLabel(blockMeshDict.lookup("nx")));
    const label  ny(readLabel(blockMeshDict.lookup("ny")));
    const label  nz(readLabel(blockMeshDict.lookup("nz")));
    // 获取网格尺寸和网格划分
    const scalar Lx = mesh.bounds().max().x() - mesh.bounds().min().x();
    const scalar Ly = mesh.bounds().max().y() - mesh.bounds().min().y();
    const scalar Lz = mesh.bounds().max().z() - mesh.bounds().min().z();
    const scalar minX = mesh.bounds().min().x();
    const scalar minY = mesh.bounds().min().y();
    const scalar minZ = mesh.bounds().min().z();
    const scalar maxX = mesh.bounds().max().x();
    const scalar maxY = mesh.bounds().max().y();
    const scalar maxZ = mesh.bounds().max().z();
    // 计算单元格大小
    const scalar dx = Lx / nx;
    const scalar dy = Ly / ny;
    const scalar dz = Lz / nz;


    // read parameters of fish 
    std::ifstream i("circleProperties.json");
    nlohmann::json j;
    i >> j;
    vector centerCoordinate;
    int objNum = j["number"]; // number of object
    double diameter = j["diameter"]; 
    int nSection = j["nSection"];
    double xx = j["center of the circle"][0];
    double yy = j["center of the circle"][1];
    double zz = j["center of the circle"][2];
    centerCoordinate.x() = xx;
    centerCoordinate.y() = yy;
    centerCoordinate.z() = zz;


    
    std::vector<std::unique_ptr<IBM>> objs;
    scalar currentTime = runTime.value();
    for (label i = 0; i < objNum; i++)
    {
        auto obj = std::make_unique<IBM>(mesh, dx, dy, dz, nx, ny, nz, minX, minY, minZ, maxX, maxY, maxZ, objNum);
        obj->initCircle(currentTime, i, centerCoordinate, nSection, diameter);
        objs.push_back(std::move(obj));  // 将对象移动到向量中
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    label outputIndex = 0;
    for (label j = 0; j < objNum; j++)
    {
        objs[j]->writeVTK(j, outputIndex);
    }

    fileName outputFile("TEST.txt");
    OFstream os(outputFile);

    forAll(IBMf, I)
{
    IBMf[I] = vector::zero;
}


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        dimensionedScalar dt = runTime.deltaT();
        scalar ct = runTime.value();
        scalar ddt = runTime.deltaT().value();

        // update the velocity and coordinate of the fish
        for (label j = 0; j < objNum; j++)
        {
            objs[j]->updateIbpCoordinate(ct);
            objs[j]->calcDelta();
            objs[j]->updateIbpVelocity(ddt);
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        

        runTime.write();

        // output presurre force, shear force and total force on the solid surface
        volTensorField stressTensor = - nu * (fvc::grad(U) + T(fvc::grad(U)));

        for (label j = 0; j < objNum; j++)
        {
            objs[j]->writeForce(os, ct, p, stressTensor);
        }

        // output vtk 
        if (runTime.outputTime())
        {
            outputIndex += 1;
            for (label j = 0; j < objNum; j++)
            {
                objs[j]->writeVTK(j, outputIndex);
            }
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
