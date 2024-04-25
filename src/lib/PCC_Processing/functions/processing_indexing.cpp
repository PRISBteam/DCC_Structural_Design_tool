///================================ A part of the PCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains top-down (based on the (k+1)-cell types) and bottom-up (based on the (k-1)-cell types) indexing of        **/
/**  k-cells, k={0,1,2,3} using the 'special' and 'induced' cell types assigned in the PCC Processing module.                     **/
///==============================================================================================================================///

/// Standard C++ libraries (STL):
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random> // Require C++ 11 and above
// #include <execution> // Require C++ 17 and above

// external libraries
#include <Eigen/Core>
#include <Eigen/SparseCore>

// local libraries
#include "../../PCC_Objects.h"
#include "../../PCC_Support_Functions.h" // It must be here - first in this list (!)

using namespace std; // standard namespace

typedef Eigen::SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'

extern std::vector<unsigned int> CellNumbs;
extern ofstream Out_logfile_stream;
extern std::vector<std::string> PCCpaths;
extern int dim;

#include "processing_indexing.h"
///------------------------------------------------------------------
/*!
 * @details
 * @param k_cell_type
 * @param Configuration_State
 * @return unsigned int Vector 'TopDownTypes(CellNumbs.at(k_cell_type))' of the (k+1) labels induced directly by the number of incident (k+1)-cell types
 */
std::vector<unsigned int> TopDown_cell_indexing(int k_cell_type, std::vector<std::vector<unsigned int>> &Configuration_State) {
// Output vector of the (k+1) labels induced directly by the number of incident (k+1)-cell types
    std::vector<unsigned int> TopDownTypes(CellNumbs.at(k_cell_type), 0); // CellNumbs.at(k_cell_type) is the number of k-cells

// Obtaining (k+1)-cells (coloumns) - k-cells (rows) Incidence matrix using the file paths.at(5 + (dim - 3))
    SpMat FES;
    if(k_cell_type != 3) {
        FES = SMatrixReader(PCCpaths.at((k_cell_type + 4) + (dim - 3)), CellNumbs.at(k_cell_type + (dim - 3)), CellNumbs.at((k_cell_type + 1) + (dim - 3))); // k-(k+1)) sparse incidence matrix; (k_cell_type + 4) is the corresponding paths to B_k incidence matrices
    }
    else{
        cout << "ERROR: TopDownTypes() function cannot be applied to 3-cells (!)" << endl;
        exit(1);
    }
// Creating sequence from the corresponding State Vector (configuration)
    std::vector<unsigned int> special_kp1_local_sequence;
    for (auto it = Configuration_State[k_cell_type+1].begin(); it != Configuration_State[k_cell_type+1].end(); ++it)
        if(*it > 0) {
            special_kp1_local_sequence.push_back(distance(Configuration_State[k_cell_type+1].begin(), it)); // add new element to the s_cells_sequence
        } // end if(it)

    for (auto kp1: special_kp1_local_sequence) // loop over all Special Faces
        for(int k = 0; k < CellNumbs.at(k_cell_type + (dim - 3)); ++k) // loop over all Edges
            if (FES.coeff(k, kp1) != 0) TopDownTypes.at(k)++;

    return TopDownTypes;
}