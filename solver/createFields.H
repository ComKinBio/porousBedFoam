#include "createRDeltaT.H"

#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField rhoc
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

// dynamic pressure field - used externally (optional)
volScalarField pDyn
(
    IOobject
    (
        "pDyn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimPressure, 0.0)
);

// velocity
Info<< "\nReading field Uc\n" << endl;
volVectorField Uc
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phic
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Uc) & mesh.Sf()
);

#include "createMRF.H"
#include "createBed.H"

const word bedName = bioBed.bedName();

mesh.setFluxRequired(p.name());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot",dimEnergy/dimVolume/dimTime, 0.0)
);

volScalarField alphac
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alphac",dimless, 0.0)
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy Kc\n" << endl;
volScalarField Kc("K", 0.5*magSqr(Uc));

scalar alphacMin = bioBed.alphaMin();

// Update alphac from the particle locations
alphac = max(1.0 - bioBed.beta(), alphacMin);
alphac.correctBoundaryConditions();

surfaceScalarField alphacf("alphaf", fvc::interpolate(alphac));

surfaceScalarField rhocPhic
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(rhoc)*phic
);

surfaceScalarField alphaRhoPhic
(
    IOobject
    (
        "alphaRhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphacf*rhocPhic
);


Info<< "Creating turbulence model\n" << endl;
autoPtr<bedTurbulenceModel> turbulence
(
    bedTurbulenceModel::New
    (
        alphac,
        rhoc,
        Uc,
        alphaRhoPhic,
        phic,
        thermo
    )
);

Info<< "Creating combustion model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> combustion
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

#include "createRadiationModel.H"
#include "createFvOptions.H"

IOdictionary radiationDict
(
    IOobject
    (
        "radiationProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

label radiationsolverFreq = radiationDict.lookupOrDefault<label>("solverFreq", 1);
radiationsolverFreq = max(1, radiationsolverFreq);

volScalarField bedap
(
    IOobject
    (
        bedName + ".ap",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimLength, 0)
);

volScalarField bedEp
(
    IOobject
    (
        bedName + ".Ep",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
);

volScalarField bedsigmap
(
    IOobject
    (
        bedName + ".sigmap",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimLength, 0)
);
