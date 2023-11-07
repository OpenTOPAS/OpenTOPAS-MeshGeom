#ifndef TsMRCP_hh
#define TsMRCP_hh

#include "TsVGeometryComponent.hh"
#include "TsTETModelImport.hh"


class TsMRCP : public TsVGeometryComponent
{
public:
	TsMRCP(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	virtual ~TsMRCP();
	G4VPhysicalVolume* Construct();
	void PrintPhantomInformation();

	std::map<G4int, G4Material*>	GetMaterialMap()	{ return fTetData->GetMaterialMap(); }
	std::map<G4int, G4String>		GetOrganNameMap() 	{ return fTetData->GetOrganNameMap(); }

private:
	G4String fPhantomDirectory;
	G4String fSex;
	G4String fAge;
	G4String fNodeFile;
	G4String fMaterialFile;
	G4String fEleFile;

    TsTETModelImport* fTetData;
    TsMRCPParameterization* fMRCPParam;
	G4ThreeVector fPhantomSize;
	G4ThreeVector fPhantomBoxMin, fPhantomBoxMax;
	G4int fNOfTetrahedrons;
};

#endif
