#ifndef TsMRCPScorer_hh
#define TsMRCPScorer_hh

#include "TsParameterManager.hh"
#include "TsScoringManager.hh"
#include "TsOutcomeModelList.hh"
#include "TsMRCP.hh"
#include "TsMRCPParameterization.hh"
#include "TsVBinnedScorer.hh"
#include "TsGeometryManager.hh"

#include "G4EmCalculator.hh"

class TsMRCPScorer : public TsVBinnedScorer
{
public:
	TsMRCPScorer(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
			G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    void GetAppropriatelyBinnedCopyOfComponent(G4String componentName);

	G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void RestoreResultsFromFile();

	~TsMRCPScorer();

protected:
    G4bool ScoreMaterialFlag(const G4String &mat);
    void BuildMaterialMap();
    std::map<G4String, G4bool> fMaterialMap;
    std::vector<G4String> AllICRPMaterials();

	G4Material* fReferenceMaterial;
    G4String* fICRPMaterials;
    G4int fNmaterials;
    G4String* fDVHsToSum;
    G4int fNDVHsToSum;

	G4EmCalculator fEmCalculator;
    TsMRCPParameterization* fmrcpParam;

    G4bool fUseMaterialFilter;
    G4bool fReportCVolHist;
    G4bool fReportDVolHist;
    G4bool fReportDoseByTet;
    G4bool fHistogramAutoMax;

    void Output(); // Overwrites parent definition
    std::vector<G4double> fVolumeHistogramVolumes;
    G4int fNUsedVolumes;
    G4double fTotalVolume;

    G4String fVHOutFileSpec1;
    G4String fVHOutFileSpec2;
    G4String fDoseElementsFileSpec1;

    void CalculateOneValue(G4int index);

    void PrintVHASCII(std::ostream& ofile);
    void PrintVHBinary(std::ostream& ofile);
    void PrintVHHeader(std::ostream& ofile);

private:
    G4int fCountInBin;
    G4double fSum;
};

#endif
