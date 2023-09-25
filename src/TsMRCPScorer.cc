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
		: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer), fEmCalculator()
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

    // Current behavior is to always produce Absolute Differential DVH
    fReportAbsDVolHist = true;

    // Current behavior is to always produce dose to each tetrahedron
    fReportDoseByTet = true;

	// Get ICRP Material names to use as filters and precompute boolean map
	if (fPm->ParameterExists(GetFullParmName("ICRPMaterials"))){
		fICRPMaterials = fPm->GetStringVector(GetFullParmName("ICRPMaterials"));
		fNmaterials = fPm->GetVectorLength(GetFullParmName("ICRPMaterials"));
    }
    BuildMaterialMap();

    G4cout << "Scorer " << GetFullParmName("") << " is using the following " <<
        fNmaterials << " ICRP material(s): " << G4endl;
    for (G4int i=0; i<fNmaterials; i++){
	    G4cout << "    " << i << " " << fICRPMaterials[i] << G4endl;
    }
    G4cout << G4endl;

	if (!fPm->ParameterExists(GetFullParmName("UseBaseOutput")))
		fUseBaseOutput = false;
	else
    {
		fUseBaseOutput = fPm->GetBooleanParameter(GetFullParmName("UseBaseOutput"));
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

    // Get reference to MRCP paramterization. Assumes phantom name is hardcoded and that
    // we have only one phantom in the geometry
    G4String componentName = fPm->GetStringParameter(GetFullParmName("Component"));
    G4String volumeName = componentName + "/WholePhantom";
    TsVGeometryComponent* component = fGm->GetComponent(componentName);
    G4VPhysicalVolume* physVol = component->GetPhysicalVolume(volumeName);
    fmrcpParam = dynamic_cast<TsMRCPParameterization*>(physVol->GetParameterisation());

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
    fRelVolumeHistogramBinVols.resize(dvh.fBinValues.size(), 0.0);
    fHistogramLowerValues.resize(dvh.fBinValues.size(), 0.0);
    for (G4int i=0; i < std::size(dvh.fBinValues); i++){
        fHistogramLowerValues[i] = dvh.fBinEdges[i];
        fRelVolumeHistogramBinVols[i] = dvh.fBinValues[i]/totalVol;
    }
    std::cout << std::endl;
}

void TsMRCPScorer::BuildMaterialMap(){
    // Precompute a map over all ICRP materials for faster evaluation
    std::vector<G4String> allMaterials = AllICRPMaterials();
    for (auto &mat: allMaterials){
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
    return fMaterialMap[mat];
}

void TsMRCPScorer::Output() {
    // Skip histogram calculation if we have restored results from file
    if (fRestoreResultsFromFile){
        G4cout << "Using restored histogram values for outcome model." << G4endl;
        G4bool isDifVolHist = true;
        fOm->CalculateOutcome(fHistogramLowerValues, fRelVolumeHistogramBinVols,
                              isDifVolHist);
        return;
    }

    if (fUseBaseOutput) {
        G4cout << fUseBaseOutput << G4endl;
        G4cout << "Using the Output behavior of parent TsVBinnedScorer for TsMRCPScorer: " << GetName() << G4endl;
        TsVBinnedScorer::Output();
    }

    fTotalVolume = 0;
    fNUsedVolumes = 0;
    G4double vol = 0;
    // Calculate differential absolute dose volume histogram
    if (fReportAbsDVolHist) {
        G4cout << "Generating Absolute DVH " << G4endl;
        // Currently fNDivisions includes all tetrahedra of phantom, an if statement will skip
        // organs that do not match this scorer's name
        fAbsVolumeHistogramBinVols.resize(fHistogramBins, 0);
        for (G4int i = 0; i < fNDivisions; i++) {
            CalculateOneValue(i);
            if (fSum >= 0. && ScoreMaterialFlag(fmrcpParam->ComputeMaterial(i)->GetName()) ) {
                vol = fmrcpParam->GetVolumeOfTet(i);
                fTotalVolume += vol;
                fNUsedVolumes++;
                TallyHistogramValue(fHistogramLowerValues, fAbsVolumeHistogramBinVols, fSum, vol);
            }
        }
    }

    // Output to CSV
    std::ofstream ofile(fOutFileName + "_AbsDifVolHist.csv");
    ofile << "# TOPAS Version: " << fPm->GetTOPASVersion() << G4endl;
	ofile << "# Parameter File: " << fPm->GetTopParameterFileSpec() << G4endl;
	ofile << "# Results for Scorer: " << GetNameWithSplitId() << G4endl;
	ofile << "# ICRP Materials scored:";
    for (G4int i=0; i<fNmaterials; i++){
	    ofile << " " << fICRPMaterials[i];
    }
    ofile << G4endl;
    ofile << "# Absolute differential dose volume histogram over number of tetrahedra: " << fNUsedVolumes << G4endl;
    ofile << "# Total volume (mm^3): " << fTotalVolume << G4endl;
	ofile << "# BinNumber, LowerLimit of " << fQuantity;
	if (GetUnit()!="")
		ofile << " ( " << GetUnit() << " )";
	else if (fScm->AddUnitEvenIfItIsOne())
		ofile << " ( 1 )";
	ofile << ", Value" << G4endl;

    if (ofile) {
        ofile << std::setprecision(16); // for double value with 8 bytes
        for (G4int j = 0; j < fHistogramBins; j++){
            ofile << j << ", " << fHistogramLowerValues[j] << ", " << fAbsVolumeHistogramBinVols[j] << G4endl;
        }
        ofile.close();
    } else {
        G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
        G4cerr << "Output file: " << fVHOutFileSpec1 << " cannot be opened for Scorer name: " << GetName() << G4endl;
        fPm->AbortSession(1);
    }
    G4cout << "End of TsMRCPScorer output for Scorer: " << GetName() << G4endl;

    // Output full tet and dose information
    if (fReportDoseByTet){
        std::ofstream ofile(fOutFileName + "_DoseByTet.csv");
        ofile << "# TOPAS Version: " << fPm->GetTOPASVersion() << G4endl;
        ofile << "# Parameter File: " << fPm->GetTopParameterFileSpec() << G4endl;
        ofile << "# Results for Scorer: " << GetNameWithSplitId() << G4endl;
        ofile << "# ICRP Materials scored:";
        for (G4int i=0; i<fNmaterials; i++){
            ofile << " " << fICRPMaterials[i];
        }
        ofile << G4endl;
        ofile << "# Doses by individual tetrahedra" << G4endl;
        ofile << "# Index, Dose (Gy), Volume (mm^3), Counts" << G4endl;
        if (ofile) {
            ofile << std::setprecision(16); // for double value with 8 bytes
            for (G4int j = 0; j < fNDivisions; j++) {
                CalculateOneValue(j);
                if (fSum >= 0. && ScoreMaterialFlag(fmrcpParam->ComputeMaterial(j)->GetName()) ) {
                    vol = fmrcpParam->GetVolumeOfTet(j)/mm3;
                    ofile << j << ", " << fSum << ", " << vol << "," << fCountInBin << G4endl;
                }
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


std::vector<G4String> TsMRCPScorer::AllICRPMaterials(){
    return {
"Adrenal_left",
"Adrenal_right",
"ET1(0-8)",
"ET1(8-40)",
"ET1(40-50)",
"ET1(50-Surface)",
"ET2(-15-0)",
"ET2(0-40)",
"ET2(40-50)",
"ET2(50-55)",
"ET2(55-65)",
"ET2(65-Surface)",
"Oral_mucosa_tongue",
"Oral_mucosa_moth_floor",
"Oral_mucosa_lips_and_cheeks",
"Trachea",
"BB(-11--6)",
"BB(-6-0)",
"BB(0-10)",
"BB(10-35)",
"BB(35-40)",
"BB(40-50)",
"BB(50-60)",
"BB(60-70)",
"BB(70-surface)",
"Blood_in_large_arteries_head",
"Blood_in_large_veins_head",
"Blood_in_large_arteries_trunk",
"Blood_in_large_veins_trunk",
"Blood_in_large_arteries_arms",
"Blood_in_large_veins_arms",
"Blood_in_large_arteries_legs",
"Blood_in_large_veins_legs",
"Humeri_upper_cortical",
"Humeri_upper_spogiosa",
"Humeri_upper_medullary_cavity",
"Humeri_lower_cortical",
"Humeri_lower_spongiosa",
"Humeri_lower_medullary_cavity",
"Radii_cortical",
"Ulnae_cortical",
"Radii_spongiosa",
"Ulnae_spongiosa",
"Radii_medullary_cavity",
"Ulnae_medullary_cavity",
"Wrists_and_hand_bones_cortical",
"Wrists_and_hand_bones_spongiosa",
"Clavicles_cortical",
"Clavicles_spongiosa",
"Cranium_cortical",
"Cranium_spongiosa",
"Femora_upper_cortical",
"Femora_upper_spongiosa",
"Femora_upper_medullary_cavity",
"Femora_lower_cortical",
"Femora_lower_spongiosa",
"Femora_lower_medullary_cavity",
"Tibiae_cortical",
"Fibulae_cortical",
"Patellae_cortical",
"Tibiae_spongiosa",
"Fibulae_spongiosa",
"Patellae_spongiosa",
"Tibiae_medullary_cavity",
"Fibulae_medullary_cavity",
"Ankles_and_foot_cortical",
"Ankles_and_foot_spongiosa",
"Mandible_cortical",
"Mandible_spongiosa",
"Pelvis_cortical",
"Pelvis_spongiosa",
"Ribs_cortical",
"Ribs_spongiosa",
"Scapulae_cortical",
"Scapulae_spongiosa",
"Cervical_spine_cortical",
"Cervical_spine_spongiosa",
"Thoracic_spine_cortical",
"Thoracic_spine_spongiosa",
"Lumbar_spine_cortical",
"Lumbar_spine_spongiosa",
"Sacrum_cortical",
"Sacrum_spongiosa",
"Sternum_cortical",
"Sternum_spongiosa",
"Cartilage_trunk",
"Brain",
"Breast_left_adipose_tissue",
"Breast_left_glandular_tissue",
"Breast_right_adipose_tissue",
"Breast_right_glandular_tissue",
"Eye_lens_sensitive_left",
"Eye_lens_insensitive_left",
"Cornea_left",
"Aqueous_left",
"Vitreous_left",
"Eye_lens_sensitive_right",
"Eye_lens_insensitive_right",
"Cornea_right",
"Aqueous_right",
"Vitreous_right",
"Gall_bladder_wall",
"Gall_bladder_contents",
"Stomach_wall(0-60)",
"Stomach_wall(60-100)",
"Stomach_wall(100-300)",
"Stomach_wall(300-surface)",
"Stomach_contents",
"Small_intestine_wall(0-130)",
"Small_intestine_wall(130-150)",
"Small_intestine_wall(150-200)",
"Small_intestine_wall(200-surface)",
"Small_intestine_contents(-400-0)",
"Small_intestine_contents(centre--400)",
"Ascending_colon_wall(0-280)",
"Ascending_colon_wall(280-300)",
"Ascending_colon_wall(300-surface)",
"Ascending_colon_content",
"Transverse_colon_wall_right(0-280)",
"Transverse_colon_wall_right(280-300)",
"Transverse_colon_wall_right(300-surface)",
"Transverse_colon_contents_right",
"Transverse_colon_wall_left(0-280)",
"Transverse_colon_wall_left(280-300)",
"Transverse_colon_wall_left(300-surface)",
"Transverse_colon_content_left",
"Descending_colon_wall(0-280)",
"Descending_colon_wall(280-300)",
"Descending_colon_wall(300-surface)",
"Descending_colon_content",
"Sigmoid_colon_wall(0-280)",
"Sigmoid_colon_wall(280-300)",
"Sigmoid_colon_wall(300-surface)",
"Sigmoid_colon_contents",
"Rectum_wall(0-280)",
"Rectum_wall(280-300)",
"Rectum_wall(300-surface)",
"Rectum_contents",
"Heart_wall",
"Blood_in_heart_chamber",
"Kidney_left_cortex",
"Kidney_left_medulla",
"Kidney_left_pelvis",
"Kidney_right_cortex",
"Kidney_right_medulla",
"Kidney_right_pelvis",
"Liver",
"Lung(AI)_left",
"Lung(AI)_right",
"Lymphatic_nodes_ET",
"Lymphatic_nodes_thoracic",
"Lymphatic_nodes_head",
"Lymphatic_nodes_trunk",
"Lymphatic_nodes_arms",
"Lymphatic_nodes_legs",
"Muscle_head",
"Muscle_trunk",
"Muscle_arms",
"Muscle_legs",
"Oesophagus_wall(0-190)",
"Oesophagus_wall(190-200)",
"Oesophagus_wall(200-surface)",
"Oesophagus_contents",
"Pancreas",
"Pituitary_gland",
"Prostate",
"RST_head",
"RST_trunk",
"RST_arms",
"RST_legs",
"Salivary_glands_left",
"Salivary_glandss_right",
"Skin_head_insensitive",
"Skin_head_sensitive(40-100)",
"Skin_trunk_insensitive",
"Skin_trunk_sensitive(40-100)",
"Skin_arms_insensitive",
"Skin_arms_sensitive(40-100)",
"Skin_legs_insensitive",
"Skin_legs_sensitive(40-100)",
"Spinal_cord",
"Spleen",
"Teeth",
"Teeth_retention_region",
"Testis_left",
"Testis_right",
"Thymus",
"Thyroid",
"Tongue_upper(food)",
"Tongue_lower",
"Tonsils",
"Ureter_left",
"Ureter_right",
"Urinary_bladder_wall_insensitive",
"Urinary_bladder_wall_sensitive(86-193)",
"Urinary_bladder_content",
"Air_inside_body",
"nOvary_left",
"Ovary_right",
"Uterus",
"Small_intestine_contents(-500-0)",
"Small_intestine_contents(centre--500)",
"Urinary_bladder_wall_sensitive(99-212)",
"Cartilage_head",
"Urinary_bladder_wall_sensitive(71-238)",
"Skin_head_sensitive(50-100)",
"Skin_trunk_sensitive(50-100)",
"Skin_arms_sensitive(50-100)",
"Skin_legs_sensitive(50-100)",
"Urinary_bladder_wall_sensitive(116-238)",
"Urinary_bladder_wall_sensitive(111-227)",
"Cartilage_arms",
"Cartilage_legs",
"Urinary_bladder_wall_sensitive(54-232)",
};
}
