
label O2Id = slgThermo.carrierId("O2");
label H2Id = slgThermo.carrierId("H2");
label H2OId = slgThermo.carrierId("H2O");


volScalarField& GasYO2 = composition.Y(O2Id);
volScalarField& GasYH2 = composition.Y(H2Id);
volScalarField& GasYH2O = composition.Y(H2OId);

// Wi kg/kmol
scalar WO2 = composition.Wi(O2Id);
scalar WH2 = composition.Wi(H2Id);
scalar WH2O = composition.Wi(H2OId);

