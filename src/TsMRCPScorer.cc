// Scorer for TsMRCPScorer
#include "TsMRCPScorer.hh"
#include "MeshGeomTools.hh"

#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include "G4String.hh"
#include "G4VPVParameterisation.hh"

TsMRCPScorer::TsMRCPScorer(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
		G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
		: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer), fEmCalculator(), fUseMaterialFilter(false), fReportDoseByTet(false), fNUsedVolumes(0), fTotalVolume(0.0), fHistogramAutoMax(true)
{
	SetUnit("Gy");

	G4String scoringMaterial = "G4_WATER";
	if (!fPm->ParameterExists(GetFullParmName("Material")))
		fReferenceMaterial = GetMaterial(scoringMaterial);
	else
	{
		G4cout << "ELSE" << G4endl;

		scoringMaterial = fPm->GetStringParameter(GetFullParmName("Material"));
		fReferenceMaterial = GetMaterial(scoringMaterial);
		if (!fReferenceMaterial) {
			G4cout << "Topas is exiting due to error in scoring setup." << G4endl;
			G4cout << "Unknown material, " << scoringMaterial << ", specified in: " << GetFullParmName("Material") << G4endl;
			fPm->AbortSession(1);
		}
	}

    G4String* reportValues;
	if (fPm->ParameterExists(GetFullParmName("Report"))) {
		reportValues = fPm->GetStringVector(GetFullParmName("Report"));
        G4int nValues = fPm->GetVectorLength(GetFullParmName("Report"));
        for (G4int i = 0; i < nValues; i++){
            if (reportValues[i] == "differentialvolumehistogram") {
                fReportDVolHist = true;
            } else if (reportValues[i] == "cumulativevolumehistogram") {
                fReportCVolHist = true;
            }
        }
	}


    if (fPm->ParameterExists(GetFullParmName("HistogramAutoMax"))){
        fHistogramAutoMax = fPm->GetBooleanParameter(GetFullParmName("HistogramAutoMax"));
    }

    if (fPm->ParameterExists(GetFullParmName("ReportDoseByTet"))){
        fReportDoseByTet = fPm->GetBooleanParameter(GetFullParmName("ReportDoseByTet"));
    }

    // Check files for writing
    G4int increment = 0;
    if (fReportCVolHist || fReportDVolHist){
        if (fOutputToBinary){
            fVHOutFileSpec1 = ConfirmCanOpen(fOutFileName+"_VolHist", ".bin", increment);
            fVHOutFileSpec2 = ConfirmCanOpen(fOutFileName+"_VolHist", ".binheader", increment);
        }
        else if (fOutputToCsv){
            fVHOutFileSpec1 = ConfirmCanOpen(fOutFileName+"_VolHist", ".csv", increment);
        }
    }
    if (fReportDoseByTet){
        // Binary output not implemented for dose by tet
        fDoseElementsFileSpec1 = ConfirmCanOpen(fOutFileName+"_DoseByTet", ".csv", increment);
    }

    // Get reference to MRCP paramterization. Assumes phantom name is hardcoded and that
    // we have only one phantom in the geometry
    G4String componentName = fPm->GetStringParameter(GetFullParmName("Component"));
    G4String volumeName = componentName + "/WholePhantom";
    TsVGeometryComponent* component = fGm->GetComponent(componentName);
    G4VPhysicalVolume* physVol = component->GetPhysicalVolume(volumeName);
    fmrcpParam = dynamic_cast<TsMRCPParameterization*>(physVol->GetParameterisation());

	// Get ICRP Material names to use as filters and precompute boolean map
	if (fPm->ParameterExists(GetFullParmName("ICRPMaterials"))){
		fICRPMaterials = fPm->GetStringVector(GetFullParmName("ICRPMaterials"));
		fNmaterials = fPm->GetVectorLength(GetFullParmName("ICRPMaterials"));
        fUseMaterialFilter = true;
        BuildMaterialMap();

        G4cout << "Scorer " << GetFullParmName("") << " is using the following " <<
            fNmaterials << " ICRP material(s): " << G4endl;
        for (G4int i=0; i<fNmaterials; i++){
            G4cout << "    " << i << " " << fICRPMaterials[i] << G4endl;
        }
        G4cout << G4endl;
    }

    G4bool restoreResultsFromFile = false;
    if (fPm->ParameterExists(GetFullParmName("RestoreResultsFromFile"))){
        restoreResultsFromFile = fPm->GetBooleanParameter(GetFullParmName("RestoreResultsFromFile"));
    }

    if (restoreResultsFromFile){
        G4cout << "Restoring results from files " << G4endl;
        // Summing DVH for restore results from file method
        if (fPm->ParameterExists(GetFullParmName("DVHsToSum"))){
            fDVHsToSum = fPm->GetStringVector(GetFullParmName("DVHsToSum"));
            fNDVHsToSum = fPm->GetVectorLength(GetFullParmName("DVHsToSum"));
            for (G4int i=0; i<fNDVHsToSum; i++){
                G4cout << " " << fDVHsToSum[i];
            }
            G4cout << G4endl;
        }
        else {
            G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
            G4cerr << "To restore results from file for the TsMRCP scorer DVHsToSum must be set." << G4endl;
            fPm->AbortSession(1);
        }
    }

    // // Test to check that the ids of the parameterisation match the IDs of fEvtMap/fFirstMomentMap
    // for (int i=0; i<10; i++){
    //     G4cout << G4endl;
    //     // Get vertices from parameterization
    //     G4Tet* tet = mrcpParam->GetTetrahedron(i);
    //     std::vector<G4ThreeVector> verts = tet->GetVertices();
    //     G4cout << "mrcpParam " << i << G4endl;
    //     G4cout << "param anchor " << verts[0][0] << " " << verts[0][1] << " " << verts[0][2] << G4endl;
    //     // Get vertices using the ComputeSolid method (this is the same index used by fEvtMap)
    //     tet = dynamic_cast<G4Tet*>(mrcpParam->ComputeSolid(i, physVol));
    //     verts = tet->GetVertices();
    //     G4cout << "computesolid anchor " << verts[0][0] << " " << verts[0][1] << " " << verts[0][2] << G4endl;
    // }
}

void TsMRCPScorer::GetAppropriatelyBinnedCopyOfComponent(G4String componentName)
{
    // Overwriting method of TsVScorer and TsVBinnedScorer to properly set number of divisions
    // Could possibly be accomplished through tscomponent::getdivisioncount method

	G4String componentNameLower = componentName;
#if GEANT4_VERSION_MAJOR >= 11
	G4StrUtil::to_lower(componentNameLower);
#else
	componentNameLower.toLower();
#endif
	if (componentNameLower == "world") {
		G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetName() << " is attempting to score in the World component." << G4endl;
		G4cerr << "This is the one component that can not have scorers." << G4endl;
		fPm->AbortSession(1);
	}

	G4String componentTypeString = "Ge/"+componentName+"/Type";
	G4String componentType = fPm->GetStringParameterWithoutMonitoring(componentTypeString);
#if GEANT4_VERSION_MAJOR >= 11
	G4StrUtil::to_lower(componentType);
#else
	componentType.toLower();
#endif

	fComponent = fGm->GetComponent(componentName);
	if (!fComponent) {
		G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetName() << " has unfound component:" << componentName << G4endl;
		fPm->AbortSession(1);
	}

	G4long testInLong = fmrcpParam->GetNumTetrahedron();

	if (testInLong > INT_MAX) {
		G4cerr << "Component " << GetName() << " has too many divisions." << G4endl;
		G4cerr << "Maximum allowed total number is " << INT_MAX << G4endl;
		G4cerr << "Number of Divisions was found to be: " << testInLong << G4endl;
		fPm->AbortSession(1);
	}

	fNDivisions = testInLong;
    fNBins = fNDivisions;
    G4cout << "***** Setting number of divisions " << fNDivisions << G4endl;

    // Ni/Nj/Nk used for output writing
    fNi = fNDivisions;
    fNj = 1;
    fNk = 1;

    // TODO: If there is an appropriately binned parallel copy of the component, use that instead
	// G4String nameWithCopyId = componentName + fGm->GetCopyIdFromBinning(fNi, fNj, fNk);
	// TsVGeometryComponent* parallelComponent = fGm->GetComponent(nameWithCopyId);
	// if (parallelComponent)
	// 	fComponent = parallelComponent;

	// Store component name including any copy Id so that if the geometry gets rebuilt, we can restore
	fComponentName = fComponent->GetNameWithCopyId();

	fDetector = fScm->GetDetector(fComponentName, this);
}

TsMRCPScorer::~TsMRCPScorer() {}

G4bool TsMRCPScorer::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive)
	{
		fSkippedWhileInactive++;
		return false;
	}

	// Find physical volume
	G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
	G4VPVParameterisation* physParam = physVol->GetParameterisation();

	ResolveSolid(aStep); // Sets fSolid
    // get index of G4Solid (could maybe be retrieved from fSolid instead
    G4int idx = ((G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable()))
		                               ->GetReplicaNumber(indexDepth);

    // G4Material* mat = fmrcpParam.ComputeMaterial(idx);
	G4double edep = aStep->GetTotalEnergyDeposit();
	if ( edep > 0. && ScoreMaterialFlag(aStep->GetPreStepPoint()->GetMaterial()->GetName()) )
	{
		G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();
		G4double dose = edep / (density * fSolid->GetCubicVolume()); // Change necessary from TsVScorer to avoid segmentation fault
		dose *= aStep->GetPreStepPoint()->GetWeight();

		G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
		G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();

		if (particle->GetPDGCharge() != 0)
		{
			// Convert to dose to material for charged particles from EM tables
			G4double materialStoppingPower = fEmCalculator.ComputeTotalDEDX(energy, particle, aStep->GetPreStepPoint()->GetMaterial());
			if (materialStoppingPower == 0.)
			{
				G4cout << "Scorer " << GetName() << " has omitted a hit that had energy: " << energy << " MeV." << G4endl;
				G4cout << "This hit could not be counted since the relevant material stopping power was unknown." << G4endl;
				return false;
			}
			G4double referenceMaterialStoppingPower = fEmCalculator.ComputeTotalDEDX(energy, particle, fReferenceMaterial);
			dose *= (density / fReferenceMaterial->GetDensity()) * (referenceMaterialStoppingPower / materialStoppingPower);
		}
		else
		{
			// Convert to dose to material for neutral particles from EM tables, assuming scaling factor from 100MeV protons
			G4double materialStoppingPower = fEmCalculator.ComputeTotalDEDX(100*MeV, G4Proton::ProtonDefinition(), aStep->GetPreStepPoint()->GetMaterial());
			if (materialStoppingPower == 0.)
			{
				G4cout << "Scorer " << GetName() << " has omitted a hit that had energy: " << energy << " MeV." << G4endl;
				G4cout << "This hit could not be counted since the relevant material stopping power was unknown." << G4endl;
				return false;
			}
			G4double referenceMaterialStoppingPower = fEmCalculator.ComputeTotalDEDX(100*MeV, G4Proton::ProtonDefinition(), fReferenceMaterial);
			dose *= (density / fReferenceMaterial->GetDensity()) * (referenceMaterialStoppingPower / materialStoppingPower);
		}
		AccumulateHit(aStep, dose, idx);
		return true;
	}
	return false;
}

void TsMRCPScorer::RestoreResultsFromFile()
{
    fRestoreResultsFromFile = true;
    MeshGeomTools::AbsDoseVolumeHistogram dvh;
    // Import CSVs
    if (fNDVHsToSum == 1){
        G4cout << "Restoring results from " << fDVHsToSum[0] << G4endl;
        dvh.loadCSV(fDVHsToSum[0]);
    }
    else {
        G4cout << "Summing DVHs: " << G4endl;
        for (G4int i=0; i<fNDVHsToSum; i++) {
            G4cout << "    " << fDVHsToSum[i] << G4endl;
            MeshGeomTools::AbsDoseVolumeHistogram newdvh(fDVHsToSum[i]);
            dvh = MeshGeomTools::sumDVHs(dvh, newdvh);
        }
    }

    // Write combined CSV to file
    std::stringstream header;
	header << "# Combined DVH for Scorer: " << GetNameWithSplitId() << G4endl;
	header << "# Combined from files:";
    for (G4int i=0; i<fNDVHsToSum; i++) {
        header << " " << fDVHsToSum[i];
    }
    header << G4endl;

    dvh.writeCSV(GetNameWithSplitId() + "_restored_DVH.csv", header.str());

    // Set appropriate values for outcome modeling
    G4double totalVol = 0;
    for (G4int i=0; i<std::size(dvh.fBinValues); i++){
        totalVol += dvh.fBinValues[i];
    }
    fVolumeHistogramVolumes.resize(dvh.fBinValues.size(), 0.0);
    fHistogramLowerValues.resize(dvh.fBinValues.size(), 0.0);
    for (G4int i=0; i < std::size(dvh.fBinValues); i++){
        fHistogramLowerValues[i] = dvh.fBinEdges[i];
        fVolumeHistogramVolumes[i] = dvh.fBinValues[i]/totalVol;
    }
    std::cout << std::endl;
}

void TsMRCPScorer::BuildMaterialMap(){
    // Precompute a boolean map over all materials in materials file
    for (auto &mat: fmrcpParam->GetMaterialNames()){
        fMaterialMap[mat] = false;
        for (int i=0; i<fNmaterials; i++){
            if (mat == fICRPMaterials[i]){
                fMaterialMap[mat] = true;
                break;
            }
        }
    }
}

G4bool TsMRCPScorer::ScoreMaterialFlag(const G4String &mat){
    if (fUseMaterialFilter){
        return fMaterialMap[mat];
    }
    else {
        return true;
    }
}

void TsMRCPScorer::Output() {
    // Skip histogram calculation if we have restored results from file
    if (fRestoreResultsFromFile){
        G4cout << "Using restored histogram values for outcome model." << G4endl;
        G4bool isDifVolHist = true;
        fOm->CalculateOutcome(fHistogramLowerValues, fVolumeHistogramVolumes,
                              isDifVolHist);
        return;
    }

    // Print sum of doses
    G4double totalSum = 0.0;
    for (G4int i=0; i<fNDivisions; i++){
        CalculateOneValue(i);
        totalSum += fSum;
    }
    G4cout << "Total sum of doses: " << totalSum << G4endl;


    // Handle volume histograms
    if (fReportCVolHist || fReportDVolHist) {
        // Recalculate histogram bin max value now that we have results
        G4double newMax = 0.0;
        if (fHistogramAutoMax) {
            for (int idx = 0; idx < fNDivisions; idx++) {
                CalculateOneValue(idx);
                if (fSum >= newMax) {
                    newMax = fSum;
                }
            }
            if (newMax > 0.0) {
                fHistogramMax = newMax;
                // Ensure bins above max so there is not overflow
                G4double binWidth = ( fHistogramMax - fHistogramMin ) / (fHistogramBins - 2);
                for (G4int i=0; i < fHistogramBins; i++){
                    fHistogramLowerValues[i] = fHistogramMin + i * binWidth;
                }
            }
        }
        fVolumeHistogramVolumes.resize(fHistogramBins, 0.0);

        G4double vol;
        for (G4int i = 0; i < fNDivisions; i++) {
            CalculateOneValue(i);
            if (ScoreMaterialFlag(fmrcpParam->ComputeMaterial(i)->GetName()) ) {
                vol = fmrcpParam->GetVolumeOfTet(i);
                fTotalVolume += vol;
                fNUsedVolumes++;
                TallyHistogramValue(fHistogramLowerValues, fVolumeHistogramVolumes, fSum, vol);
            }
        }

        // Normalize histogram
        for (G4int i=0; i < fHistogramBins; i++){
            fVolumeHistogramVolumes[i] /= fTotalVolume;
        }

        // Calculate cumulative volume histogram
        if (fReportCVolHist) {
            for (G4int i = fHistogramBins-1; i >= 1; i--) {
                fVolumeHistogramVolumes[i - 1] += fVolumeHistogramVolumes[i];
            }
        }

        if (fOutputToCsv) {
            std::ofstream ofile(fVHOutFileSpec1);
            if (ofile) {
                PrintVHHeader(ofile);
                PrintVHASCII(ofile);
                ofile.close();
            } else {
                G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
                G4cerr << "Output file: " << fVHOutFileSpec1 << " cannot be opened for Scorer name: " << GetName() << G4endl;
                fPm->AbortSession(1);
            }
        } else if (fOutputToBinary) {
            std::ofstream hfile(fVHOutFileSpec1);
            if (hfile) {
                PrintVHHeader(hfile);
                hfile << "# Binary file: " << fVHOutFileSpec2 << G4endl;
                hfile.close();
            } else {
                G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
                G4cerr << "Output file: " << fVHOutFileSpec1 << " cannot be opened for Scorer name: " << GetName() << G4endl;
                fPm->AbortSession(1);
            }

            std::ofstream ofile(fVHOutFileSpec2, std::ios::out | std::ios::binary);
            if (ofile) {
                PrintVHBinary(ofile);
                ofile.close();
            } else {
                G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
                G4cerr << "Output file: " << fVHOutFileSpec2 << " cannot be opened for Scorer name: " << GetName() << G4endl;
                fPm->AbortSession(1);
            }
        }
    }

    // Output full tet and dose information
    if (fReportDoseByTet){
        G4double vol;
        std::ofstream ofile(fOutFileName + "_DoseByTet.csv");
        ofile << "# TOPAS Version: " << fPm->GetTOPASVersion() << G4endl;
        ofile << "# Parameter File: " << fPm->GetTopParameterFileSpec() << G4endl;
        ofile << "# Results for Scorer: " << GetNameWithSplitId() << G4endl;
        if (fUseMaterialFilter){
            ofile << "# ICRP Materials scored:";
            for (G4int i=0; i<fNmaterials; i++){
                ofile << " " << fICRPMaterials[i];
            }
            ofile << G4endl;
        }
        ofile << "# Doses by individual tetrahedra" << G4endl;
        ofile << "# Index, Dose (Gy), Volume (mm^3), Counts" << G4endl;
        if (ofile) {
            ofile << std::setprecision(16); // for double value with 8 bytes
            for (G4int j = 0; j < fNDivisions; j++) {
                CalculateOneValue(j);
                vol = fmrcpParam->GetVolumeOfTet(j)/mm3;
                ofile << j << ", " << fSum << ", " << vol << "," << fCountInBin << G4endl;
            }
            ofile.close();
        } else {
            G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
            G4cerr << "Output file: " << fVHOutFileSpec1 << " cannot be opened for Scorer name: " << GetName() << G4endl;
            fPm->AbortSession(1);
        }
    }

    G4cout << "End of TsMRCPScorer output for Scorer: " << GetName() << G4endl;
}


void TsMRCPScorer::PrintVHASCII(std::ostream& ofile)
{
	ofile << std::setprecision(16); // for double value with 8 bytes
	for (int j = 0; j < fHistogramBins; j++)
		ofile << j << ", " << fHistogramLowerValues[j] << ", " << fVolumeHistogramVolumes[j] << G4endl;
	ofile << std::setprecision(6);
}


void TsMRCPScorer::PrintVHBinary(std::ostream& ofile)
{
	G4double* data = new G4double[fHistogramBins];
    for (int j = 0; j < fHistogramBins; j++){
        data[j] = fVolumeHistogramVolumes[j];
    }
	ofile.write( (char*) data, fHistogramBins*sizeof(G4double));
	delete[] data;
}

void TsMRCPScorer::PrintVHHeader(std::ostream& ofile)
{
    ofile << "# TOPAS Version: " << fPm->GetTOPASVersion() << G4endl;
    ofile << "# Parameter File: " << fPm->GetTopParameterFileSpec() << G4endl;
	ofile << "# Results for Scorer: " << GetNameWithSplitId() << G4endl;

	if (fReportCVolHist)
		ofile << "# Cumulative Volume Histogram over number of voxels: " << fNDivisions << G4endl;
    else if (fReportDVolHist)
		ofile << "# Differential Volume Histogram over number of voxels: " << fNDivisions << G4endl;

	ofile << "# BinNumber, LowerLimit of " << fQuantity;

	if (GetUnit()!="")
		ofile << " ( " << GetUnit() << " )";
	else if (fScm->AddUnitEvenIfItIsOne())
		ofile << " ( 1 )";

	ofile << ", Value" << G4endl;
}

void TsMRCPScorer::CalculateOneValue(G4int idx)
{
    fCountInBin = fCountMap[idx];
    fSum = fFirstMomentMap[idx] / GetUnitValue();
}
