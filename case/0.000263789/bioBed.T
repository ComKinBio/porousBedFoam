/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0.000263789";
    object      bioBed.T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   nonuniform List<scalar> 9(338.0086864 338.0092183 338.0092666 338.0088057 338.0093097 1200 338.0089832 1200 1200);

boundaryField
{
    top
    {
        type            calculated;
        value           uniform 0;
    }
    bottom
    {
        type            calculated;
        value           uniform 0;
    }
    walls
    {
        type            calculated;
        value           uniform 0;
    }
}


// ************************************************************************* //
