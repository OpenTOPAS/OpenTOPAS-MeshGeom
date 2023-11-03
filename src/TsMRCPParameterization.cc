// Extra Class for TsMRCP
#include "TsMRCPParameterization.hh"
#include "G4TransportationManager.hh"


TsMRCPParameterization::TsMRCPParameterization(TsTETModelImport* tetData)
	: G4VPVParameterisation(), fTetData(tetData) {
    // Navigator is used to query the material at a given x, y, z location
    fNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
}

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

std::pair<G4ThreeVector, G4ThreeVector> TsMRCPParameterization::GetMaterialExtent(const G4String material){
    return fTetData->GetMaterialExtent(material);
}

const G4String TsMRCPParameterization::GetMaterialAtPoint(const G4ThreeVector point){
    G4VPhysicalVolume* volume = fNavigator->LocateGlobalPointAndSetup(point);
    return volume->GetLogicalVolume()->GetMaterial()->GetName();
}

std::vector<G4String> TsMRCPParameterization::GetMaterialNames(){
    std::vector<G4String> names;
    std::map<G4int, G4String> organNameMap = fTetData->GetOrganNameMap();
    for(auto &it:  organNameMap){
      names.push_back(it.second);
    }
    return names;
}
