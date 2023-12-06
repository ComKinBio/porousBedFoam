{
    volScalarField& he = thermo.he();

    const volScalarField alphaEff(turbulence->alphaEff());

    fvScalarMatrix EEqn
    (
        fvm::ddt(alphac, rhoc, he) + mvConvection->fvmDiv(alphaRhoPhic, he)
      + fvc::ddt(alphac, rhoc, Kc) + fvc::div(alphaRhoPhic, Kc)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(alphaRhoPhic/fvc::interpolate(rhoc), Uc),
                p,
                "div(phiv,p)"
            )
          : -dpdt*alphac
        )
      - fvm::laplacian                
        (
            fvc::interpolate(alphac)
           *fvc::interpolate(alphaEff),
            he
        )
     ==
        alphac*rhoc*(Uc&g)
      + alphac*Qdot
      + bioBed.Sh(he)
      + radiation->Sh(thermo, he)
      + fvOptions(rhoc, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
