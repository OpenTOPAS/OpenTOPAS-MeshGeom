// Extra Class for use by TsTetGeom
#include <MeshGeomTools.hh>

#include <algorithm>
#include <iterator>

namespace MeshGeomTools {

AbsDoseVolumeHistogram::AbsDoseVolumeHistogram()
{
    fBinEdges = {0.0, 1.0};
    fBinValues = {0.0};
};

AbsDoseVolumeHistogram::AbsDoseVolumeHistogram(G4String file) 
{
    loadCSV(file);
};

void AbsDoseVolumeHistogram::loadCSV(G4String file) 
{
    readCSV(file, fBinEdges, fBinValues);
    // TOPAS standard is for histograms to be lower(left) edges of bins
    // Remove last bin value if non-zero to remove overflow
    if(fBinValues[fBinValues.size()-1] == 0.0){
        fBinValues.pop_back();
    }
    else {
        G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
        G4cerr << "The histogram in " << file << " has overflow values that have no specified ";
        G4cerr << "upper limit."<< G4endl;
        exit(1);
    }
};


AbsDoseVolumeHistogram::~AbsDoseVolumeHistogram() {};

void AbsDoseVolumeHistogram::writeCSV(G4String file, G4String header)
{
    std::fstream fout;
	// Creates a new file, overwritting if exists
	fout.open(file, std::ios::out | std::ios::trunc);
    fout << header;
	// Insert data
	for (G4int i = 0; i<fBinValues.size(); i++)
	{
		fout << fBinEdges[i] << ", " << fBinEdges[i+1] << ", " << fBinValues[i] << "\n";
	}
	fout.close();
}

AbsDoseVolumeHistogram sumDVHs(AbsDoseVolumeHistogram hist1, AbsDoseVolumeHistogram hist2) 
{
    // Form union grid
    AbsDoseVolumeHistogram outHist;
    std::vector<G4double> result;
    std::set_union(hist1.fBinEdges.begin(), hist1.fBinEdges.end(), 
                   hist2.fBinEdges.begin(), hist2.fBinEdges.end(), 
                   std::back_inserter(outHist.fBinEdges));
    std::sort(outHist.fBinEdges.begin(), outHist.fBinEdges.end());

    outHist.fBinValues.resize(outHist.fBinEdges.size(), 0);

    std::vector<G4double> newVals1 = mapHistValsToGrid(hist1.fBinEdges, hist1.fBinValues, outHist.fBinEdges);
    std::vector<G4double> newVals2 = mapHistValsToGrid(hist2.fBinEdges, hist2.fBinValues, outHist.fBinEdges);
    for (G4int i=0; i<newVals1.size(); i++){
        newVals1[i] += newVals2[i];
    }
    outHist.fBinValues = newVals1;
    return outHist;
}

std::vector<G4double> mapHistValsToGrid(std::vector<G4double> &bins, std::vector<G4double> &vals,
                                        std::vector<G4double> &newbins)
{
    // Assuming a linear distribution across each bin, form the new DVH on the union 
    // grid of the old two
    // This function assumes that bins is a subset of newbins
    std::vector<G4double> newvals(newbins.size()-1, 0);
    G4int j=0;
    for (G4int i=0; i < vals.size(); i++){
        while (newbins[j] <= bins[i]){
            if (newbins[j+1] <= bins[i+1]){
                newvals[j] += vals[i] * (newbins[j+1] - newbins[j])/(bins[i+1] - bins[i]);
            }
            j++;
        }
    }
    return newvals;
};

extern "C" void testHistogramSumming(char* file1, char* file2, char* outfile)
{
    AbsDoseVolumeHistogram hist1(file1);
    AbsDoseVolumeHistogram hist2(file2);
    AbsDoseVolumeHistogram hist3 = sumDVHs(hist1, hist2);
    hist3.writeCSV(outfile);
}

} // end namespace
