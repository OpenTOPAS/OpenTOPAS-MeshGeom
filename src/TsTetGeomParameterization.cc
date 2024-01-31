// Extra Class for TsTetGeom
#include "TsTetGeomParameterization.hh"

#include "G4TransportationManager.hh"

TsTetGeomParameterization::TsTetGeomParameterization(TsTetModelImport* tetData)
	: G4VPVParameterisation(), fTetData(tetData) {
}

TsTetGeomParameterization::~TsTetGeomParameterization() {}

void TsTetGeomParameterization::InitializeNavigator(const G4String world){
    G4TransportationManager* transportationManager = G4TransportationManager::GetTransportationManager();
	fNavigator = transportationManager->GetNavigator(transportationManager->GetParallelWorld(world));
}

G4VSolid* TsTetGeomParameterization::ComputeSolid(const G4int copyNo, G4VPhysicalVolume*)
{
	return fTetData->GetTetrahedron(copyNo);
}

void TsTetGeomParameterization::ComputeTransformation(const G4int, G4VPhysicalVolume*) const { }

G4Material* TsTetGeomParameterization::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* phy, const G4VTouchable*)
{
	return fTetData->GetMaterial(fTetData->GetMaterialIndex(copyNo));
}

G4Material* TsTetGeomParameterization::ComputeMaterial(const G4int copyNo)
{
    return fTetData->GetMaterial(fTetData->GetMaterialIndex(copyNo));
}

G4int TsTetGeomParameterization::GetNumTetrahedron(){
    return fTetData->GetNumTetrahedron();
}

G4Tet* TsTetGeomParameterization::GetTetrahedron(const G4int copyNo){
    return fTetData->GetTetrahedron(copyNo);
}

G4double TsTetGeomParameterization::GetVolumeOfTet(const G4int copyNo){
    return fTetData->GetVolumeOfTet(copyNo);
}

std::pair<G4ThreeVector, G4ThreeVector> TsTetGeomParameterization::GetMaterialExtent(const G4String material){
    return fTetData->GetMaterialExtent(material);
}

std::vector<G4String> TsTetGeomParameterization::GetMaterialNames(){
    std::vector<G4String> names;
    std::map<G4int, G4String> organNameMap = fTetData->GetOrganNameMap();
    for(auto &it:  organNameMap){
      names.push_back(it.second);
    }
    return names;
}

G4String TsTetGeomParameterization::GetMaterialAtPoint(const G4ThreeVector point){
    // Method is for paramterized geometry, for non-parameterized use something like
    // physVol->GetVolume()->GetMaterial()->GetName()
    G4VPhysicalVolume* volume = fNavigator->LocateGlobalPointAndSetup(point);
    return fTetData->GetMaterial(fTetData->GetMaterialIndex(volume->GetCopyNo()))->GetName();
}
