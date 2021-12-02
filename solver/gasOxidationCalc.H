
// H2 + 0.5O2 = H2O
// k = 1e11*exp(10000/T)C_H2*C_O2^0.5, kmol/m3/s
// reaction heat = 142919 kJ/kg

forAll(GasYH2,i)
{
    
    if (GasYH2[i] > small)
    {
        scalar C_H2 = rhoc[i]*GasYH2[i]/WH2;
        scalar C_O2 = rhoc[i]*GasYO2[i]/WO2;

        scalar kH2 = 1e11*Foam::exp(-10000.0/T[i])*C_H2*Foam::sqrt(C_O2);
        
        scalar dMoleH2 = min(C_O2*2, min(C_H2, kH2*runTime.deltaTValue()));
       
        scalar dYH2 = dMoleH2*WH2/rhoc[i];
        scalar dYH2O = dMoleH2*WH2O/rhoc[i];
        scalar dYO2 = dMoleH2*0.5*WO2/rhoc[i];

        GasYH2[i] = max(0.0, GasYH2[i] - dYH2);
        GasYO2[i] = max(0.0, GasYO2[i] - dYO2);
        GasYH2O[i] = max(0.0, GasYH2O[i] + dYH2O);
        
        scalar H2ReactionHeat = dMoleH2*WH2*142919000/runTime.deltaTValue();

        Qdot[i] += H2ReactionHeat;

    }
}
