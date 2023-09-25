#ifndef TsMRCPParameterization_hh
#define TsMRCPParameterization_hh

#include "G4VPVParameterisation.hh"
#include "G4Tet.hh"
#include "TsTETModelImport.hh"

class TsMRCPParameterization : public G4VPVParameterisation
{
public:
	TsMRCPParameterization(TsTETModelImport* tetData);
	virtual ~TsMRCPParameterization();

	G4VSolid* ComputeSolid(const G4int copyNo, G4VPhysicalVolume* );
	void ComputeTransformation(const G4int, G4VPhysicalVolume*) const;
	G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* phy, const G4VTouchable*);
    G4Material* ComputeMaterial(const G4int copyNo);
    G4int       GetNumTetrahedron();
    G4Tet*      GetTetrahedron(const G4int copyNo);
    G4double      GetVolumeOfTet(const G4int copyNo);

private:
	TsTETModelImport* fTetData = NULL;
};

#endif
