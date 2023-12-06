Info<< "\nConstructing fuel bed" << endl;

basicGravityBio2DBed bioBed
(
    "bioBed",
    mesh,
    rhoc,
    Uc,
    g,
    slgThermo
);
