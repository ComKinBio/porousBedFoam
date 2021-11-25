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
    object      bioBed.Cp;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 -1 0 0 0];

internalField   nonuniform List<scalar> 9(1493.534498 1493.534498 1493.534498 1493.534498 1493.534498 0 1493.534498 0 0);

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
