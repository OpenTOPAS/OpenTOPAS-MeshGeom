// Component for TsMRCP
#include "TsParameterManager.hh"
#include "G4Box.hh"
#include "TsMRCPParameterization.hh"

#include "TsMRCP.hh"

using namespace std;

TsMRCP::TsMRCP(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
		 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
				TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    // Define tetrahedral input either by (PhantomDirectory, Age, Sex) or by (PhantomDirectory, NodeFile, MaterialFile, EleFile)

    if (fPm->ParameterExists(GetFullParmName("PhantomDirectory")) &&
        fPm->ParameterExists(GetFullParmName("Age")) &&
        fPm->ParameterExists(GetFullParmName("Sex")) ) {
        G4cout << "Importing an ICRP tetrahedral mesh" << G4endl;

	    fPhantomDirectory = fPm->GetStringParameter(GetFullParmName("PhantomDirectory"));
        // Get parameters to identify the desired phantom
        if ( fPm->ParameterExists(GetFullParmName("Age")) )
        {
            fAge = fPm->GetStringParameter(GetFullParmName("Age"));
            if (fAge.length() > 2)
            {
                G4cerr << "Age for peadiatric phantoms need to be specified as a 2-digit string; for adults with 'A'." << G4endl;
                fPm->AbortSession(1);
            }
        }
        else
        {
            G4cerr << "Age needs to be specified to read the corresponding phantom." << G4endl;
            fPm->AbortSession(1);
        }

        if ( fPm->ParameterExists(GetFullParmName("Sex")))
        {
            fSex = fPm->GetStringParameter(GetFullParmName("Sex"));
            if (fSex != "F" && fSex != "M")
            {
                G4cerr << "Sex needs to be specified as F or M." << G4endl;
                fPm->AbortSession(1);
            }
        }
        else
        {
            G4cerr << "Sex needs to be specified to read the corresponding phantom." << G4endl;
            fPm->AbortSession(1);
        }

        G4String age = fAge;
        G4String phantomName = "MRCP_" + fAge + fSex;
        // Get data using TsTETModelImport
        fTetData = new TsTETModelImport(fPhantomDirectory, phantomName);

    }
    else if (fPm->ParameterExists(GetFullParmName("PhantomDirectory")) &&
             fPm->ParameterExists(GetFullParmName("NodeFile")) &&
             fPm->ParameterExists(GetFullParmName("MaterialFile")) &&
             fPm->ParameterExists(GetFullParmName("EleFile")) ) {
        G4cout << "Importing tetrahedral mess from file paths" << G4endl;
	    fPhantomDirectory = fPm->GetStringParameter(GetFullParmName("PhantomDirectory"));
        fNodeFile = fPm->GetStringParameter(GetFullParmName("NodeFile"));
        fMaterialFile = fPm->GetStringParameter(GetFullParmName("MaterialFile"));
        fEleFile = fPm->GetStringParameter(GetFullParmName("EleFile"));
        fTetData = new TsTETModelImport(fPhantomDirectory, fNodeFile, fMaterialFile, fEleFile);
    }
    else {
        G4cerr << "Must specify either (PhantomDirectory, Age, Sex) " << G4endl;
        G4cerr << "or (PhantomDirectory, NodeFile, MaterialFile, EleFile) " << G4endl;
        G4cerr << "to import tetrahedral geometry" << G4endl;
        fPm->AbortSession(1);
    }

	// Get variables for phantom information
	fPhantomSize = fTetData->GetPhantomSize();
	fPhantomBoxMin = fTetData->GetPhantomBoxMin();
	fPhantomBoxMax = fTetData->GetPhantomBoxMax();
	fNOfTetrahedrons = fTetData->GetNumTetrahedron();
}

TsMRCP::~TsMRCP() {}

G4VPhysicalVolume* TsMRCP::Construct()
{
	BeginConstruction();

	// Define the phantom container (buffer of 10% of the absolute max in each dimension)
    G4double xAbsMax = max(abs(fPhantomBoxMin.x()), abs(fPhantomBoxMax.x()));
    G4double yAbsMax = max(abs(fPhantomBoxMin.y()), abs(fPhantomBoxMax.y()));
    G4double zAbsMax = max(abs(fPhantomBoxMin.z()), abs(fPhantomBoxMax.z()));

	G4double totalHLX = xAbsMax + 0.1*xAbsMax;
	G4double totalHLY = yAbsMax + 0.1*yAbsMax;
	G4double totalHLZ = zAbsMax + 0.1*zAbsMax;

	G4Box* envelopeSolid = new G4Box(fName, totalHLX, totalHLY, totalHLZ);
    // Set envelop material to reduce number of warnings
    // of scorer called in non-parameterised volume
	G4String envelopeMaterial = "G4_AIR";

	fEnvelopeLog = CreateLogicalVolume(fName, envelopeMaterial, envelopeSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

	// Define the tetrahedral mesh phantom as a parameterized geometry
	G4VSolid* tetraSolid = new G4Tet("TetSolid", G4ThreeVector(), 
                                     G4ThreeVector(1.*cm,0,0),
                                     G4ThreeVector(0,1.*cm,0),
                                     G4ThreeVector(0,0,1.*cm));
	G4String material = "G4_Galactic";
	G4LogicalVolume* tetLogic = CreateLogicalVolume("TetLogic", material, tetraSolid);

	// Physical volume (phantom) constructed as parameterized geometry
    TsMRCPParameterization* fMRCPParam = new TsMRCPParameterization(fTetData);
	G4VPhysicalVolume* phantom = CreatePhysicalVolume("WholePhantom", tetLogic, fEnvelopePhys, kUndefined, fNOfTetrahedrons, fMRCPParam);

	InstantiateChildren();
    PrintPhantomInformation();

    // Initialize navigator in TsMRCPParameterization for querying material at a point
    fMRCPParam->InitializeNavigator(this->GetWorldName());

	return fEnvelopePhys;
}

void TsMRCP::PrintPhantomInformation()
{
	fTetData->PrintMaterialInformation();
}
