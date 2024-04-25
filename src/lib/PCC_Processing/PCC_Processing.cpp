///======================================= PCC Processing module ================================================================================ ///
///============================================================================================================================================= ///
///* The interface use functions from Processing_<***>_functions.h C++ libraries to generate quasi-random or non-random processes of           *///
///* labelling of the k-cells ( k={0,1,2,3} ) of the pre-constructed polytopal cell complex (PCC) 'M' using the whole set of incidence        *///
///* and adjacency matrices                                                                                                                   *///
///* ---------------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2024 years as a module of the PCC Processing Design code (CPD code) *///
///* A part of the MATERiA codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years   *///
/// https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1                                                                  *///
///===================================================================================================================================== ///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// Attached user-defined C++ libraries:
// External
#include <Eigen/SparseCore>

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../ini/ini_readers.h"
#include "../PCC_Objects.h"

// Local
///---------------------------------------------------------
#include "functions/processing_assigned_labelling.h"
#include "functions/processing_induced_labelling.h"
#include "functions/processing_indexing.h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs;
extern ofstream Out_logfile_stream;
extern string source_path;
extern std::vector<std::string> PCCpaths;
extern int dim;

#include "PCC_Processing.h"
///* ========================================================= PCC PROCESSING FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
/*!
 * @details Employ functions from Processing_<***>_functions.h C++ libraries to generate quasi-random or non-random processes of labelling of the k-cells ( k={0,1,2,3} )
 * of the pre-constructed polytopal cell complex (PCC) 'M' using the whole set of incidence and adjacency matrices.
 * All the initial settings are written in 'processing.ini' file, including PCCpaths to the corresponding directories and execution types.
 * The key theoretical concepts are special (SCS), ordinary (OCS) k-cell sequences, state vectors (SVs), design vectors (DVs) and configurations.
 * @param configuration
 * @return CellDesign object
 */
CellDesign PCC_Processing(Config &configuration) {
/// Main output of the module 'special_cells_design' (CD) - class. In particular, it contains (1) special_nodes_sequence, (2) special_edges_sequence, (3) special_faces_sequence (in the 2D and 3D cases), and (4) special_polyhedrons_sequence (in the 3D case)
    CellDesign CD;

    std::vector<std::vector<unsigned int>> Configuration_sState = configuration.Get_Configuration_sState(); // definition of the local State Vectors with values from the object of Config class
    std::vector<std::vector<unsigned int>> Configuration_cState = configuration.Get_Configuration_iState(); // definition of the local State Vectors with values from the object of Config class

/// Processing vector variables
    std::vector <unsigned int> special_x_sequence; // variable sequence of 'special' (SCS) (different from 'ordinary' (OCS)) k-cells in order of their generation.
    std::vector <unsigned int> induced_x_sequence; // variable sequence of 'induced' (ICS) k-cells in order of their generation as the result of a kinetic process (historically they are also called "fractured" cells).
    std::vector<std::vector<unsigned int>> special_x_series; // vector of series of special k-cells (strips/chains)
    std::vector<std::vector<unsigned int>> induced_x_series; // vector of series of induced k-cells (like crack paths)

    double mu_f = 1.0, sigma_f = 0.0; // Only in the case of the 'L' PCC_Processing execution mode ('pp_mode') :: mean and dispersion for a lengthy defect sequences distribution.

    cout << "=========================================================================" << endl; Out_logfile_stream << "==============================================================================================================================================================" << endl;

/// Read configuration from processing.ini file :: the number of special cell types and calculation parameters.
    std::vector<vector<double>> max_sfractions_vectors(4), max_ifractions_vectors(4); // double vectors containing maximal values for Processing execution of fractions for [0][..] - nodes, [1][..] - edges, [2][..] - faces, [3][..] - polyhedrons
    std::vector<string> stype_vector(4), itype_vector(4); // special_types and induced_types vector of strings corresponding to the Processing execution type ON/OFF: {0,1} for all the possible 'k' values of k-cell read from the file 'config/processing.ini'; in both vectors: [0] - nodes, [1] - edges, [2] - faces, [3] - polyhedrons
    std::vector<double> pindex_vector(4); // supplementary index read from the 'config/processing.ini' file
    std::vector<string> sequence_source_paths(4); // k-sequence paths for the reading them from file(s) in the 'S' PCC_Processing execution mode (pp_mode) instead of the new Processing routes
    bool multiplexity = 0; // Bool variable indicating one {multiplexity = 0} or many {multiplexity = 1} labels can be assigned for each k-cell number

    // Reading of the configuration from the 'config/processing.ini' file
    config_reader_processing(source_path, sequence_source_paths, max_sfractions_vectors, max_ifractions_vectors, mu_f, sigma_f, stype_vector, itype_vector, pindex_vector, Out_logfile_stream); // void function

///* Cases for Processing types // cell type k :: 0 - nodes, 1 - edges, 2 - faces, 3 -polyhedrons - must coincide with the indexing of the CellNumbs.at(cell_type) vector (!) *///
///* ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- *///
/// Here '(dim - 3)' term is important for 1D and 2D cases (!) as in these cases there are no polyhedrons or even faces (2D) in the PCC
    for (int cell_type = (3 + (dim - 3)); cell_type >= 0; --cell_type) { // Loop over all types of k-cells in the PCC
        // clearence of vectors for each new cell type
        special_x_sequence.clear(); special_x_series.clear();

    /// I. Beginning of the processing of 'special' k-cells
    ///=======================================================
        if (stype_vector.at(cell_type) == "R" && max_sfractions_vectors.at(cell_type).size() > 0) { //  Random separate cells generation processing
            cout << "Random (R) mode processing in operation: cell_type : "s << cell_type << endl; Out_logfile_stream << "Random (R) mode processing in operation: cell_type : "s << cell_type << endl;

            multiplexity = (bool) pindex_vector.at(cell_type); // convert 'pindex' read from the 'config/processing.ini' file for the specific 'cell_type' to a bool variable 'multiplexity'.

            special_x_series = Processing_Random(cell_type, Configuration_sState, max_sfractions_vectors, multiplexity); // defined in the assigned labelling library

            /// Writing 'special_x_sequence' - each k-cell number appeared only once in the sequence. Configuration_sState and State Vectors show possible 'agglomeration' - several similar label per k-cell
            for (auto it = Configuration_sState[cell_type].begin(); it != Configuration_sState[cell_type].end(); ++it)
                if(*it > 0) {
                    special_x_sequence.push_back(distance(Configuration_sState[cell_type].begin(), it)); // add new element to the s_cells_sequence
                } // end if(it)

        } // End of 'R' type simulations (if)

        else if (stype_vector.at(cell_type) == "L" && max_sfractions_vectors.at(cell_type).size() > 0) { //  Random lengthy strips (chains) of cells generation processing
            cout << "Random strips/chains (L) mode processing in operation: cell_type : "s << cell_type << endl << endl; Out_logfile_stream << "Random strips/chains (L) mode processing in operation: cell_type : "s << cell_type << endl << endl;
            cout << "Average (mu) and dispersion (sigma): " << endl  << mu_f << "  " << sigma_f << endl;

            /// Obtaining the distribution of strip/chain lengths
            unsigned int bins_number = 50;
            std::vector<double> strip_lenghts_distribution = Log_normal_distribution(mu_f, sigma_f, bins_number); // double valued "continuous" distribution, where Log_normal_distribution() function is for obtaining strip_lenghts_distribution
            std::vector<unsigned int> cell_strip_distribution; // vector of positive integers containing "discrete" length distribution of special chains/strips of k-cells

            cell_strip_distribution.clear(); // clearing strip/chain length distribution vector for new 'cell_type' iterations
            for (auto  itr = strip_lenghts_distribution.begin(); itr != strip_lenghts_distribution.end(); ++itr) {
                cell_strip_distribution.push_back(max_sfractions_vectors[cell_type][0] * CellNumbs.at(cell_type) * (*itr) / (std::distance(strip_lenghts_distribution.begin(), itr) + 1.0)); // only for the first 'Xmax_fraction1' in the 'config/processing.ini' file values of max k-cell fractions
            }

            /// Random_Strips_Distribution() function call
            special_x_series = Processing_Random_Strips(cell_type, cell_strip_distribution, Configuration_sState, max_sfractions_vectors); // series of k-cells for each strip/chain

            /// Writing 'special_x_sequence' - each k-cell number appeared only once in the sequence. Configuration_sState and State Vectors show possible 'agglomeration' - several similar label per k-cell
            for (auto it = Configuration_sState[cell_type].begin(); it != Configuration_sState[cell_type].end(); ++it)
                if(*it > 0) {
                    special_x_sequence.push_back(distance(Configuration_sState[cell_type].begin(), it)); // add new element to the s_cells_sequence
                } // end if(it)

        } //End of 'L' type simulations (elseif)

        else if (stype_vector.at(cell_type) == "F" && max_sfractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
            // processing index :: 0 - direct special faces assignment;  1 - crystallographic ; 2 - configurational TJs-based entropy (deviatoric); //        if (pindex_vector.at(cell_type) == 0) { //        } else if (pindex_vector.at(cell_type) == 1) {
            cout << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
            Out_logfile_stream << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
///           if(cell_type == 2 + (dim - 3))             // cell type = 2 -> faces
///                special_x_sequence = Processing_maxFunctional(cell_type, Configuration_sState, max_fractions_vectors, pindex_vector.at(cell_type));
        } // End of 'F' type simulations (elseif)

        else if (stype_vector.at(cell_type) == "D" && max_sfractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
            cout << "Min (MAX-deviator) Functional processing in operation: cell_type : "s << cell_type << endl;
            Out_logfile_stream << "Min (MAX-deviator) Functional processing in operation: cell_type : "s << cell_type << endl;
///            if (max_fractions_vectors.at(cell_type).size() > 0)
///            special_x_sequence = Processing_minConfEntropy(2, Configuration_sState, max_fractions_vectors, pindex_vector.at(2));

        } // End of 'D' [S min] type simulations (elseif)

        else if (stype_vector.at(cell_type) == "S") { /// Reading structure from file
            char* kseq_sourcepath = const_cast<char*>(sequence_source_paths.at(cell_type).c_str());
            special_x_sequence = VectorIReader(kseq_sourcepath);

        // (!) Output +1 like in Neper, so he numbers should be modified back as -1
            for (auto it = special_x_sequence.begin(); it != special_x_sequence.end(); ++it)
                special_x_sequence.at(distance(special_x_sequence.begin(), it)) = *it - 1;

        // Cut up to max_fraction (!!)
            std::vector<unsigned int> temp_x_sequence = special_x_sequence; // temporarily new vector
            double total_max_sCell_fraction_processing = 0;
            for (int j = 0; j < max_sfractions_vectors[cell_type].size(); ++j)
                if(max_sfractions_vectors[cell_type][j] > 0)
                    total_max_sCell_fraction_processing += max_sfractions_vectors[cell_type][j];

            temp_x_sequence.clear();
            for (unsigned int k = 0; k < total_max_sCell_fraction_processing*CellNumbs.at(cell_type); ++k)
                temp_x_sequence.push_back(special_x_sequence.at(k));

            special_x_sequence.clear(); // putting everything back
            for(unsigned int p : temp_x_sequence)
                special_x_sequence.push_back(p);
            temp_x_sequence.clear();

        // Update of the corresponding Configuration State vector
            Configuration_sState[cell_type].clear();
            std::vector<int> State_vector(CellNumbs.at(cell_type), 0);
            for (unsigned int k_cell : special_x_sequence) // fill state vector from cells_sequence
                State_vector.at(k_cell) = 1;

            for (int var : State_vector)
                Configuration_sState[cell_type].push_back(var);
        } // End of 'S' [reading from file] type simulations (elseif)

        else if(max_sfractions_vectors[cell_type].size() > 0) cout << "ERROR [Processing] : unknown simulation type - please replace with 'R', 'L', 'F', 'D' or 'S'..!" << endl;

    ///* ONLY for 2-cells or 'grain boundaries' *///
        if (cell_type == (2 + (dim - 3))) { // '..+ (dim - 3)' because in the 2D case, boundaries become edges (!) and grains become faces of the corresponding 2D tessellation

            if (stype_vector.at(cell_type) == "Cm" &&
                max_sfractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
                cout << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
                Out_logfile_stream << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
///                if (max_sfractions_vectors.at(cell_type).size() > 0)
///                    special_x_sequence = Processing_maxF_crystallographic(2, Configuration_sState, max_sfractions_vectors, pindex_vector.at(2));
            } // end of 'Cm' type simulations (elseif)

            else if (stype_vector.at(cell_type) == "Cr" &&
                     max_sfractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
                cout << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
                Out_logfile_stream << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
///                if (max_sfractions_vectors.at(cell_type).size() > 0)
///                    special_x_sequence = Processing_maxP_crystallographic(2, Configuration_sState, max_fractions_vectors, pindex_vector.at(2));
                    //special_x_sequence = Processing_Random_crystallographic(2, Configuration_sState, max_fractions_vectors, pindex_vector.at(2));
            } // end of 'Cr' type simulations (elseif)

            /// II. Beginning of the processing of 'indexing' k-cells (TopDown and BottomUp)
            ///=======================================================

            if (special_x_sequence.size() < 2) { // ONLY if there was no processing and 'special' labelling of cells (!)
                cout << " Start TopDown_cell_indexing() based on the cell type k+1: " << cell_type << endl;

                std::vector<unsigned int> StateVector_indexing;
                StateVector_indexing = TopDown_cell_indexing(cell_type, Configuration_sState);

                /// 'special_x_sequence' - each k-cell number appeared only once in the sequence obtained by Configuration_sState
                for (auto it = StateVector_indexing.begin(); it != StateVector_indexing.end(); ++it)
                    if(*it > 0) {
                        special_x_sequence.push_back(distance(StateVector_indexing.begin(), it)); // add new element to the s_cells_sequence
                    } // end if(it)
                }

            /// III. Beginning of the processing of 'induced' (of 'fractured') k-cells
            ///=======================================================

            if (itype_vector.at(cell_type) == "Km" && max_ifractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
                cout << "Induced processing in operation: cell_type : "s << cell_type << endl;
                Out_logfile_stream << "Induced processing in operation: cell_type : "s << cell_type << endl;

                if (max_ifractions_vectors.at(cell_type).size() > 0)
                     induced_x_sequence = PCC_Kinematic_cracking(cell_type, special_x_sequence, Configuration_cState, max_ifractions_vectors);
            } // End of 'Km' type simulations (elseif)

            else if (itype_vector.at(cell_type) == "Kn" && max_ifractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
                cout << "Induced processing in operation: cell_type : "s << cell_type << endl; Out_logfile_stream << "Induced processing in operation: cell_type : "s << cell_type << endl;

///                if (max_ifractions_vectors.at(cell_type).size() > 0)
///                    induced_x_sequence = PCC_Kinetic_cracking(Configuration_sState, face_elastic_energies, large_crack);
            } // End of 'Km' type simulations (elseif)

            else if(max_ifractions_vectors[cell_type].size() > 0) cout << "ERROR [Processing] : unknown induced simulation type - please replace with 'Km' or 'Kn'..!" << endl;

        } // End of if (cell_type == 2)
// REPAIR    cout << "ctype_vector " << ctype_vector.at(cell_type + (3 - 3)) << "  " << max_cfractions_vectors[cell_type + (3 - 3)].size() << endl;

    /// Assigned sequences:
        CD.Set_special_sequence(special_x_sequence, cell_type); // (sequence, id)
        CD.Set_special_series(special_x_series, cell_type); // (series, id)
        CD.Set_special_configuration(Configuration_sState.at(cell_type), cell_type); // (configuration/design, id) - design vector with types

    /// Induced sequences:
        CD.Set_induced_sequence(induced_x_sequence, cell_type); // (sequence, ctype)

    } // END of for (int cell_type = 3; cell_type >= 0; --cell_type)

    if (configuration.Get_main_type() == "LIST"s) {
        cout << endl;
        Out_logfile_stream << endl;
        cout << "n-sequence size: " << CD.Get_n_special_sequence().size() << endl; Out_logfile_stream << "n-sequence size: " << CD.Get_n_special_sequence().size() << endl;
        cout << "e-sequence size: " << CD.Get_e_special_sequence().size() << endl; Out_logfile_stream << "e-sequence size: " << CD.Get_e_special_sequence().size() << endl;
        cout << "f-sequence size: " << CD.Get_f_special_sequence().size() << endl; Out_logfile_stream << "f-sequence size: " << CD.Get_f_special_sequence().size() << endl;
        cout << "p-sequence size: " << CD.Get_p_special_sequence().size() << endl; Out_logfile_stream << "p-sequence size: " << CD.Get_p_special_sequence().size() << endl;
        cout << endl; Out_logfile_stream << endl;

        cout << "n-induced-sequence size: " << CD.Get_n_induced_sequence().size() << endl; Out_logfile_stream << "n-induced-sequence size: " << CD.Get_n_induced_sequence().size() << endl;
        cout << "e-induced-sequence size: " << CD.Get_e_induced_sequence().size() << endl; Out_logfile_stream << "e-induced-sequence size: " << CD.Get_e_induced_sequence().size() << endl;
        cout << "f-induced-sequence size: " << CD.Get_f_induced_sequence().size() << endl; Out_logfile_stream << "f-induced-sequence size: " << CD.Get_f_induced_sequence().size() << endl;
        cout << "p-induced-sequence size: " << CD.Get_p_induced_sequence().size() << endl; Out_logfile_stream << "p-induced-sequence size: " << CD.Get_p_induced_sequence().size() << endl;
        cout << endl; Out_logfile_stream << endl;

        cout << "n-design vector size: " << CD.Get_n_design().size() << endl; Out_logfile_stream << "n-design vector size: " << CD.Get_n_design().size() << endl;
        cout << "e-design vector size: " << CD.Get_e_design().size() << endl; Out_logfile_stream << "e-design vector size: " << CD.Get_e_design().size() << endl;
        cout << "f-design vector size: " << CD.Get_f_design().size() << endl; Out_logfile_stream << "f-design vector size: " << CD.Get_f_design().size() << endl;
        cout << "p-design vector size: " << CD.Get_p_design().size() << endl; Out_logfile_stream << "p-design vector size: " << CD.Get_p_design().size() << endl;
        cout << endl; Out_logfile_stream << endl;
    } // End of if (configuration.Get_main_type() == "LIST"s)

    return CD;
} /// The END of PCC_Processing() function
