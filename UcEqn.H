    MRF.correctBoundaryVelocity(Uc);

    fvVectorMatrix UcEqn(Uc, rhoc.dimensions()*Uc.dimensions()*dimVol/dimTime);

    UcEqn = 
    (
        fvm::ddt(alphac, rhoc, Uc) + fvm::div(alphaRhoPhic, Uc)
      + MRF.DDt(alphac*rhoc, Uc)
      + turbulence->divDevRhoReff(Uc)
     ==
        alphac*rhoc()*g
      + bioBed.SU(Uc)
      + fvOptions(rhoc, Uc)
    );

    UcEqn.relax();

    fvOptions.constrain(UcEqn);

    if (pimple.momentumPredictor())
    {
        solve(UcEqn == -alphac*fvc::grad(p));

        fvOptions.correct(Uc);
        Kc = 0.5*magSqr(Uc);
    }
