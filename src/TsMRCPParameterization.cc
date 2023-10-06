// Extra Class for TsMRCP
#include "TsMRCPParameterization.hh"

TsMRCPParameterization::TsMRCPParameterization(TsTETModelImport* tetData)
	: G4VPVParameterisation(), fTetData(tetData) { }

TsMRCPParameterization::~TsMRCPParameterization() {}

G4VSolid* TsMRCPParameterization::ComputeSolid(const G4int copyNo, G4VPhysicalVolume*)
{
	return fTetData->GetTetrahedron(copyNo);
}

void TsMRCPParameterization::ComputeTransformation(const G4int, G4VPhysicalVolume*) const { }

G4Material* TsMRCPParameterization::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* phy, const G4VTouchable*)
{
	return fTetData->GetMaterial(fTetData->GetMaterialIndex(copyNo));
}

G4Material* TsMRCPParameterization::ComputeMaterial(const G4int copyNo)
{
    return fTetData->GetMaterial(fTetData->GetMaterialIndex(copyNo));
}

G4int TsMRCPParameterization::GetNumTetrahedron(){
    return fTetData->GetNumTetrahedron();
}

G4Tet* TsMRCPParameterization::GetTetrahedron(const G4int copyNo){
    return fTetData->GetTetrahedron(copyNo);
}

G4double TsMRCPParameterization::GetVolumeOfTet(const G4int copyNo){
    return fTetData->GetVolumeOfTet(copyNo);
}

std::vector<G4String> TsMRCPParameterization::GetMaterialNames(){
    std::vector<G4String> names;
    std::map<G4int, G4String> organNameMap = fTetData->GetOrganNameMap();
    for(auto &it:  organNameMap){
      names.push_back(it.second);
    }
    return names;
}
