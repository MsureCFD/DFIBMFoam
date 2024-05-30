// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

/*Ustar = U;    

forAll(IBMf, I)
{
    IBMf[I] = vector::zero;
}
for (label i =0; i < 10; i++)
{
    tU = Ustar - dt*fvc::grad(p) + IBMf*dt;
    //tU = Ustar + IBMf*dt;

    for (label j = 0; j < objNum; j++)
    {
        objs[j]->interpolateVelocity(tU);
        objs[j]->distributeForce(IBMf, ddt);
    }

}*/

Ustar = U;    
forAll(IBMf, I)
{
    IBMf[I] = vector::zero;
}


for (label i =0; i < 10; i++)
{
    tU = Ustar - dt*fvc::grad(p) + IBMf*dt;
    //tU = U - dt*fvc::grad(p) + IBMf*dt;
    //tU = Ustar + IBMf*dt;

    for (label j = 0; j < objNum; j++)
    {
        objs[j]->interpolateVelocity(tU);
        objs[j]->distributeForce(IBMf, ddt);
    }

}

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U) + IBMf
);



fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);



if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}















