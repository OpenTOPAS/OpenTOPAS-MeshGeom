#ifndef MeshGeomTools_hh
#define MeshGeomTools_hh

#include "globals.hh"

#include <iostream>
#include <fstream>
#include <vector>

namespace MeshGeomTools {

class AbsDoseVolumeHistogram 
{
public:
    AbsDoseVolumeHistogram();
    AbsDoseVolumeHistogram(G4String file);
    ~AbsDoseVolumeHistogram();

    // bin_edges is len(bin_values) + 1
    std::vector<G4double> fBinEdges;
    std::vector<G4double> fBinValues;
    void writeCSV(G4String file, G4String header="");
    void loadCSV(G4String file);
};

AbsDoseVolumeHistogram sumDVHs(AbsDoseVolumeHistogram hist1, AbsDoseVolumeHistogram hist2);
std::vector<G4double> mapHistValsToGrid(std::vector<G4double> &bins, std::vector<G4double> &values,
                                        std::vector<G4double> &newbins);

extern "C" void testHistogramSumming(char* file1, char* file2, char* outfile);

template<class T1, class T2>
void readCSV(G4String csv, std::vector<T1> &bins, std::vector<T2> &values)
{
    bins.clear();
    values.clear();

    // File pointer
    std::ifstream f;
    // Open csv file
    f.open(csv, std::ios::in);

    if (f.fail())
    {
        G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
        G4cerr << "Unable to open provided csv: " << csv << G4endl;
        exit(1);
    }
    G4int count = 0;
    G4String line, value;
    std::vector<G4String> row;

    // Read file
    while (getline(f, line))
    {
        if (line.substr(0,1) == "#") {G4cout << line << G4endl; continue;}
        row.clear();
        std::stringstream s(line);
        while (getline(s, value, ',')) { row.push_back(value); }
        bins.push_back(stod(row[1]));
        values.push_back(stod(row[2]));
        count++;
    }
    // Close file
    f.close();
}


} // end namespace
#endif
