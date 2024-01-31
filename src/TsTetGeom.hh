#ifndef TsTetGeom_hh
#define TsTetGeom_hh

#include "TsVGeometryComponent.hh"
#include "TsTetModelImport.hh"
#include "TsTetGeomParameterization.hh"


class TsTetGeom : public TsVGeometryComponent
{
public:
	TsTetGeom(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	virtual ~TsTetGeom();
	G4VPhysicalVolume* Construct();
	void PrintPhantomInformation();

	std::map<G4int, G4Material*>	GetMaterialMap()	{ return fTetData->GetMaterialMap(); }
	std::map<G4int, G4String>		GetOrganNameMap() 	{ return fTetData->GetOrganNameMap(); }
    G4int GetDivisionCount(G4int dim);

private:
	G4String fPhantomDirectory;
	G4String fSex;
	G4String fAge;
	G4String fNodeFile;
	G4String fMaterialFile;
	G4String fEleFile;

    TsTetModelImport* fTetData;
    TsTetGeomParameterization* fTetGeomParam;
	G4ThreeVector fPhantomSize;
	G4ThreeVector fPhantomBoxMin, fPhantomBoxMax;
	G4int fNOfTetrahedrons;
};

#endif
