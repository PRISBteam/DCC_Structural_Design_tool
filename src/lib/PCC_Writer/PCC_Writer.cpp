///================================ PCC Writer module ===================================================================================///
///=====================================================================================================================================///
/** The function in this library perform a structured output of the calculation results **/
///===================================================================================================================================///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// Attached user defined C++ libraries:
#include <Eigen/SparseCore>

///-------------------------------------
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../PCC_Objects.h"
#include "../ini/ini_readers.h"


#include "functions/Writer_functions.h"
///-------------------------------------

using namespace std; // standard namespace

extern std::ofstream Out_logfile_stream;
extern std::string output_dir, source_path;

#include "PCC_Writer.h"
/// ---------------------------------------------------------------------------

/// # 1 # sequences only
void PCC_Writer(CellsDesign &new_cells_design) {
// Read PCC Writer specifications from the writer.ini file and output of the current configuration to the screen and .log file
    std::vector<int> writer_specifications; // vector<int> containing writer specifications and formats
    config_reader_writer(source_path, writer_specifications, Out_logfile_stream); // Read and output the initial configuration from the writer.ini file

    int output_counter = 0; // special counter for output numeration
    
    if (writer_specifications.at(0) == 1)
        PCC_CellSequences_Writer(new_cells_design, output_counter);

    return;
} /// END of PCC Writer module

/// # 2 # overloaded
void PCC_Writer(CellsDesign &new_cells_design, ProcessedComplex &pcc_processed) {
// Read PCC Writer specifications from the writer.ini file and output of the current configuration to the screen and .log file
    std::vector<int> writer_specifications; // vector<int> containing writer specifications and formats
    config_reader_writer(source_path, writer_specifications, Out_logfile_stream); // Read and output the initial configuration from the writer.ini file

    int output_counter = 0; // special counter for output numeration

    if (writer_specifications.at(0) == 1)
        PCC_CellSequences_Writer(new_cells_design, output_counter);
/**
    if (writer_specifications.at(2) == 1)
        PCC_Entropic_Writer(pcc_processed, output_counter);
    cout << "2" << endl;

    if (writer_specifications.at(5) == 1)
        PCC_AnalyticalFractions_Writer(pcc_processed, output_counter); // analytical solutions
    cout << "3" << endl;

    if(writer_specifications.at(7) == 1)
        PCC_Laplacians_Writer(pcc_processed, output_counter);
    cout << "4" << endl;
**/
/*
       char* stype = const_cast<char*>(P_type.c_str()); // Processing type variable
        char P_simulation_type = *stype;

  if (P_simulation_type == 'L') {
        Strips_distribution_Writer(output_counter); // strips distribution
        RWseries_Writer(RW_series_vector, output_counter); // vector of all strips

        /// agglomerations  output
        vector<agglomeration> agglomerations_vector; // vector containing aglomarations as its objects - only for L type simulations (!)
        Agglomerations_Writer(output_counter, RW_series_vector, mu_f, sigm_f);

    } // end of if(P_simulation_type == 'L')
*/
    return ;
} // END of PCC Writer module