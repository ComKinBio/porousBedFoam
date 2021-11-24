Info<< "\nConstructing fuel bed" << endl;

basicGravityBio2DBed bioBed
(
    "bed2d",
    mesh,
    rhoc,
    Uc,
    g,
    slgThermo
);
