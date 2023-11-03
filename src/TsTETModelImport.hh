#ifndef TsTETModelImport_hh
#define TsTETModelImport_hh

#include "G4SystemOfUnits.hh"
#include "G4Tet.hh"
#include "G4UIExecutive.hh"
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"

class TsTETModelImport
{
public:
	TsTETModelImport(G4String path, G4String phantomName);
	TsTETModelImport(G4String path, G4String nodeFile, G4String materialFile, G4String eleFile);
	~TsTETModelImport();

	// Get methods
	G4String		GetPhantomName()			{ return phantomName; }
	G4Material*		GetMaterial(G4int idx)		{ return materialMap[idx]; }
	G4int			GetNumTetrahedron()			{ return tetVector.size(); }
	G4int			GetMaterialIndex(G4int idx)	{ return materialVector[idx]; }
	G4Tet*			GetTetrahedron(G4int idx)	{ return tetVector[idx]; }
	G4double		GetVolume(G4int idx)		{ return volumeMap[idx]; }
	G4double		GetVolumeOfTet(G4int idx)	{ return tetVector[idx]->GetCubicVolume(); }
	std::map<G4int, G4double>	GetMassMap()	{ return massMap; }
	std::map<G4int, G4Colour>	GetColorMap()	{ return colorMap; }
	G4ThreeVector	GetPhantomSize()			{ return phantomSize; }
	G4ThreeVector	GetPhantomBoxMin()			{ return boundingBoxMin; }
	G4ThreeVector	GetPhantomBoxMax()			{ return boundingBoxMax; }
	void PrintMaterialInformation();

	std::map<G4int, G4Material*>	GetMaterialMap()	{ return materialMap; }
	std::map<G4int, G4String>		GetOrganNameMap() 	{ return organNameMap; }

private:
	// Private methods
	void DataRead(G4String, G4String);
	void MaterialRead(G4String);
	void ColorRead();

	G4String phantomDataPath;
	G4String phantomName;

	G4ThreeVector boundingBoxMin;
	G4ThreeVector boundingBoxMax;
	G4ThreeVector phantomSize;

	std::vector<G4ThreeVector> 	vertexVector;
	std::vector<G4Tet*>			tetVector;
	std::vector<G4int*>			eleVector;
	std::vector<G4int>			materialVector;
	std::map<G4int, G4int>		numTetMap;
	std::map<G4int, G4double>	volumeMap;
	std::map<G4int, G4double>	massMap;
	std::map<G4int, G4Colour>	colorMap;
    std::map<G4int, std::pair<G4ThreeVector, G4ThreeVector>> materialExtentMap;

	std::map<G4int, std::vector<std::pair<G4int, G4double>>>	materialIndexMap;
	std::vector<G4int>											materialIndex;
	std::map<G4int, G4Material*>								materialMap;
	std::map<G4int, G4double>									densityMap;
	std::map<G4int, G4String>									organNameMap;
};

#endif
