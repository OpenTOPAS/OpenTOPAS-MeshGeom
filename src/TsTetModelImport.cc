// Extra Class for TsTetGeom
#include "TsTetModelImport.hh"

TsTetModelImport::TsTetModelImport(G4String path, G4String phantomName)
{
    phantomDataPath = path;
    G4String eleFile = phantomName + ".ele";
    G4String nodeFile = phantomName + ".node";
    G4String materialFile = phantomName + ".material";
    MaterialRead(materialFile);
    DataRead(eleFile, nodeFile);
}

TsTetModelImport::TsTetModelImport(G4String path, G4String nodeFile, G4String materialFile, G4String eleFile)
{
    phantomDataPath = path;
    MaterialRead(materialFile);
    DataRead(eleFile, nodeFile);
}

TsTetModelImport::~TsTetModelImport() {}

void TsTetModelImport::DataRead(G4String eleFile, G4String nodeFile)
{
	G4String tempStr;
	G4int tempInt;

	// Read node file
	std::ifstream ifpNode;
	ifpNode.open((phantomDataPath + "/" + nodeFile).c_str());
	if (!ifpNode.is_open())
	{
		// Exception for the case there is no *.node file
		G4Exception("TsTetModelImport::DataRead", "", FatalErrorInArgument, G4String(" There is no " + phantomDataPath + "/" + nodeFile).c_str());
	}
	G4cout << "Opening TETGEN node (vertex points: x y z) file '" << nodeFile << "'" << G4endl;
    G4cout << (phantomDataPath + "/" + nodeFile).c_str() << G4endl;

	G4int numVertex;
	G4double xPos, yPos, zPos;
	G4double xMin(DBL_MAX), yMin(DBL_MAX), zMin(DBL_MAX);
	G4double xMax(DBL_MIN), yMax(DBL_MIN), zMax(DBL_MIN);

	ifpNode >> numVertex >> tempInt >> tempInt >> tempInt;

	for (G4int i = 0; i < numVertex; i++)
	{
		ifpNode >> tempInt >> xPos >> yPos >> zPos;

		// Set the unit
		xPos *= cm;
		yPos *= cm;
		zPos *= cm;

		// Save the node data as the form of std::vector<G4ThreeVector>
		vertexVector.push_back(G4ThreeVector(xPos, yPos, zPos));

		// Get the information of the bounding box of phantom
		if (xPos < xMin) xMin = xPos;
		if (xPos > xMax) xMax = xPos;
		if (yPos < yMin) yMin = yPos;
		if (yPos > yMax) yMax = yPos;
		if (zPos < zMin) zMin = zPos;
		if (zPos > zMax) zMax = zPos;
	}

	// Set the variables for the bounding box and phantom size
	boundingBoxMin = G4ThreeVector(xMin, yMin, zMin);
	boundingBoxMax = G4ThreeVector(xMax, yMax, zMax);
	phantomSize = G4ThreeVector(xMax - xMin, yMax - yMin, zMax - zMin);

	ifpNode.close();

	// Read ele file
	std::ifstream ifpEle;

	ifpEle.open((phantomDataPath + "/" + eleFile).c_str());
	if (!ifpEle.is_open())
	{
		// Exception for the case there is no *.ele file
		G4Exception("TsTetModelImport::DataRead", "", FatalErrorInArgument, G4String(" There is no " + eleFile).c_str());
	}
	G4cout << "Opening TETGEN elements (tetrahedron with node No.) file '" << eleFile << "'" << G4endl;
	G4cout << (phantomDataPath + "/" + eleFile).c_str() << G4endl;

	G4int numEle;
	ifpEle >> numEle >> tempInt >> tempInt;

    G4ThreeVector *minVec, *maxVec, *curPos;
	for (G4int i = 0; i < numEle; i++)
	{
		ifpEle >> tempInt;
		G4int* ele = new G4int[4];
		for (G4int j = 0; j < 4; j++)
		{
			ifpEle >> tempInt;
			ele[j] = tempInt;
		}
		eleVector.push_back(ele);
		ifpEle >> tempInt;
		materialVector.push_back(tempInt);

		// Save the element (tetrahedron) data as the form of std::vector<G4Tet*>
		tetVector.push_back(new G4Tet("Tet_Solid", vertexVector[ele[0]],  vertexVector[ele[1]],  vertexVector[ele[2]],  vertexVector[ele[3]]));

		// Calculate the total volume and the number of tetrahedrons for each organ
		std::map<G4int, G4double>::iterator findIter = volumeMap.find(materialVector[i]);
		if (findIter != volumeMap.end())
		{
			findIter->second += tetVector[i]->GetCubicVolume();
			numTetMap[materialVector[i]]++;
		}
		else
		{
			volumeMap[materialVector[i]] = tetVector[i]->GetCubicVolume();
			numTetMap[materialVector[i]] = 1;
		}

        // Update extent by looping over all vertices of each tet
        minVec = &(materialExtentMap[materialVector[i]].first);
        maxVec = &(materialExtentMap[materialVector[i]].second);
        for (G4int j=0; j<4; j++){
            curPos = &(vertexVector[ele[j]]);
            // G4cout << "Min: " << minVec->x() << " " << minVec->y() << " " << minVec->z() << G4endl;
            // G4cout << "Max: " << maxVec->x() << " " << maxVec->y() << " " << maxVec->z() << G4endl;
            // G4cout << "curPos: " << curPos->x() << " " << curPos->y() << " " << curPos->z() << G4endl;
            if (curPos->x() < minVec->x()) minVec->setX(curPos->x());
            if (curPos->x() > maxVec->x()) maxVec->setX(curPos->x());
            if (curPos->y() < minVec->y()) minVec->setY(curPos->y());
            if (curPos->y() > maxVec->y()) maxVec->setY(curPos->y());
            if (curPos->z() < minVec->z()) minVec->setZ(curPos->z());
            if (curPos->z() > maxVec->z()) maxVec->setZ(curPos->z());
        }
	}
	ifpEle.close();

    // Check that all materials from .ele file were defined in .material file
    for (G4int i = 0; i < (G4int)materialVector.size(); i++)
    {
        G4int idx = materialVector[i];
        if (!materialMap.count(idx)) {
            std::ostringstream oss;
            oss << "Material file did not define material number:" << idx << " found in .ele file";
            G4Exception("TsTetModelImport::DataRead", "", FatalErrorInArgument, G4String(oss.str()).c_str());
        }
    }
}

void TsTetModelImport::MaterialRead(G4String materialFile)
{
	// Read mateiral file (*.material)
	std::ifstream ifpMat;

	ifpMat.open((phantomDataPath + "/" + materialFile).c_str());
	if (!ifpMat.is_open())
	{
		// Exception for the case there is no *.material file
		G4Exception("TsTetModelImport::MaterialRead", "", FatalErrorInArgument, G4String(" There is no " + materialFile).c_str());
	}
	G4cout << "Opening material file '" << materialFile << "'" << G4endl;
	G4cout << (phantomDataPath + "/" + materialFile).c_str() << G4endl;

	char read_data[50];
	char* token;
	G4double zaid;
	G4double fraction;
	G4String materialName;
	G4double density;

	while(!ifpMat.eof())
	{
		ifpMat >> read_data;
		ifpMat >> materialName;
		ifpMat >> read_data;
		density = std::atof(read_data);
		ifpMat >> read_data;
		ifpMat >> read_data;
		if(read_data[0]!='m') continue;
		token = std::strtok(read_data, "m");
		G4int matID = std::atoi(token);
		materialIndex.push_back(matID);
		organNameMap[matID] = materialName;
        organNameToMatID[materialName] = matID;
		densityMap[matID] = density*g/cm3;

		for (G4int i = 0; ; i++)
		{
			ifpMat >> read_data;
			if(std::strcmp(read_data, "C") == 0 || ifpMat.eof()) break;

			zaid = (G4int)(std::atoi(read_data)/1000);
			ifpMat >> read_data;
			fraction = -1.0 * std::atof(read_data);
			materialIndexMap[matID].push_back(std::make_pair(G4int(zaid), fraction));
		}
	}
	ifpMat.close();

	// Construct materials for each organ
	G4NistManager* nistManager = G4NistManager::Instance();

	for (G4int i = 0; i < (G4int)materialIndex.size(); i++)
	{
		G4int idx = materialIndex[i];
		G4Material* mat = new G4Material(organNameMap[idx], densityMap[idx], G4int(materialIndexMap[idx].size()), kStateSolid, CLHEP::STP_Temperature, CLHEP::STP_Pressure);
		for (G4int j = 0; j < G4int(materialIndexMap[idx].size()); j++)
			mat->AddElement(nistManager->FindOrBuildElement(materialIndexMap[idx][j].first), materialIndexMap[idx][j].second);
		materialMap[idx] = mat;
		massMap[idx] = densityMap[idx] * volumeMap[idx];
	}

    // Initialize material extent map
	for (G4int i = 0; i < (G4int)materialIndex.size(); i++)
	{
		G4int idx = materialIndex[i];
        materialExtentMap[idx].first  = G4ThreeVector(DBL_MAX, DBL_MAX, DBL_MAX);
        materialExtentMap[idx].second = G4ThreeVector(DBL_MIN, DBL_MIN, DBL_MIN);
    }
}

void TsTetModelImport::ColorRead()
{
	// Read color data file
	std::ifstream ifpColor;

	ifpColor.open((phantomDataPath + "/" + "colour.dat").c_str());
	if (!ifpColor.is_open())
	{
		// Exception for the case there is no color.dat file
		G4Exception("TsTetModelImport::DataRead", "", FatalErrorInArgument, G4String("Color data file was not found").c_str());
	}
	G4cout << "Opening color data file 'colour.dat" << G4endl;

	G4int organID;
	G4double red, green, blue, alpha;
	while (ifpColor >> organID >> red >> green >> blue >> alpha)
		colorMap[organID] = G4Colour(red, green, blue, alpha);

	ifpColor.close();
}

void TsTetModelImport::PrintMaterialInformation()
{
	// Print the overall information for each organ
	G4cout << G4endl
		   << std::setw(9)	<< "Organ ID"
		   << std::setw(11)	<< "# of Tet"
		   << std::setw(11)	<< "vol [cm3]"
		   << std::setw(11)	<< "d [g/cm3]"
		   << std::setw(11)	<< "mass [g]"
		   << "\t" << "organ/tissue" << G4endl;
	G4cout << "-----------------------------------" << G4endl;

	std::map<G4int, G4Material*>::iterator matIter;
	G4cout << std::setiosflags(std::ios::fixed);
	G4cout.precision(3);
	for (matIter = materialMap.begin(); matIter != materialMap.end(); matIter++)
	{
		G4int idx = matIter->first;

		G4cout << std::setw(9)  << idx                         // organ ID
			   << std::setw(11) << numTetMap[idx]              // # of tetrahedrons
			   << std::setw(11) << volumeMap[idx]/cm3          // organ volume
			   << std::setw(11) << materialMap[idx]->GetDensity()/(g/cm3)      // organ density
			   << std::setw(11) << massMap[idx]/g              // organ mass
			   << "\t"<<materialMap[idx]->GetName() << G4endl; // organ name
	}

	// Print the extent of each organ 
	G4cout << G4endl
		   << std::setw(9)	<< "Organ ID"
		   << std::setw(11)	<< "Min X"
		   << std::setw(11)	<< "Min Y"
		   << std::setw(11)	<< "Min Z"
		   << std::setw(11)	<< "Max X"
		   << std::setw(11)	<< "Max Y"
		   << std::setw(11)	<< "Max Z"
		   << "\t" << "organ/tissue" << G4endl;
	G4cout << "-----------------------------------" << G4endl;

	G4cout << std::setiosflags(std::ios::fixed);
	G4cout.precision(3);
	for (matIter = materialMap.begin(); matIter != materialMap.end(); matIter++)
	{
		G4int idx = matIter->first;
		G4cout << std::setw(9)  << idx                         // organ ID
			   << std::setw(11) << materialExtentMap[idx].first.x()/cm
			   << std::setw(11) << materialExtentMap[idx].first.y()/cm
			   << std::setw(11) << materialExtentMap[idx].first.z()/cm
			   << std::setw(11) << materialExtentMap[idx].second.x()/cm
			   << std::setw(11) << materialExtentMap[idx].second.y()/cm
			   << std::setw(11) << materialExtentMap[idx].second.z()/cm
			   << "\t"<<materialMap[idx]->GetName() << G4endl; // organ name
	}
}

std::pair<G4ThreeVector, G4ThreeVector> TsTetModelImport::GetMaterialExtent(const G4String material) {
    return materialExtentMap[organNameToMatID[material]];
}
