///******************************************************************************************************************************///
///************************   Polytopal Cell Complex (PCC) Processing Design :: (CPD code) (c)   *******************************///
///****************************************************************************************************************************///
///*                                        Version 4.0 | 15/04/2024                                                         *///
///**************************************************************************************************************************///
///************************************ Dr Elijah Borodin, Manchester, UK **************************************************///
///**************************************** Spring 2022 - Spring 2024  ****************************************************///
///***********************************************************************************************************************///
///*
///*    Code source:    https://github.com/PRISBteam/PCC_Processing_Design/
///*    Documentation:  https://prisbteam.github.io/
///*    PCC sources:    https://materia.team/
///*
///*  The project provides a reliable tool for 1. Obtaining, 2. Analysing and 3. Optimising of 'design vectors' as the sequences                         *///
///*  of k-cells containing in the k-skeletons, where k = {0,1,2,3}, of a Polytopal Cell Complex (PCC). Such PCCs can be created by external             *///
///*  codes based on the tessellation of 2D or 3D spaces by an agglomeration of polytopes (polygons in the 2D case or polyhedrons in 3D).                *///
///*  Graphs and networks (without loops) are considered as 1-complexes (1-PCCs) and also available for analysis similarly to the the 2D and 3D cases.   *///

///* Key terminology:                                                                                                                                                               *///
/// Tessellation's elements     ::   'nodes, 'edges', 'faces', 'polytopes' (with their measures - lengths, areas and volumes - and barycenter coordinates)                         ///
/// PCC's elements              ::   'k-cells' containing in 'k-skeletons', where k = {0,1,2,3}, with their degree fractions, and incident (k-1)-cells and (k+1)-cells.           ///
/// Material's elements         ::   'quadruple points', 'grain boundary junctions', 'grain boundaries', and 'grains' (with their orientations and barycenter coordinates often taken from EBSD or X-ray analysis)  ///

///* ----------------------------------------- *
///* Standard C++ (STL) libraries
///* ----------------------------------------- *
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <tuple>
#include <cmath>

///* ------------------------------------------------------------------------------- *
///* Attached user-defined C++ libraries (must be copied in the directory for STL):
///* ------------------------------------------------------------------------------- *
/// Eigen source: https://eigen.tuxfamily.org/ (2024)
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

/// Spectra source: https://spectralib.org/ (2024)
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>

/// Open MP library https://www.openmp.org/resources/openmp-compilers-tools/
// Included only in the parallelized version of the code

///------------------------------------------
using namespace std; // standard/STL namespace
using namespace Eigen; // Eigen library namespace
using namespace Spectra; // Spectra library namespace

/// Eigen library-based classes
typedef Triplet<double> Tr; // <Eigen> library class, which declares a triplet type with the nickname 'Tr' as the objects in the form T = T(i, j, value), where i and j are element's a(i,j) indices in the corresponding dense matrix and the third variable is its value
typedef SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'
typedef MatrixXd DMat; // <Eigen> library class, which declares a dense matrix type of doubles with the nickname 'DMat'

/// * ---------------------------------------------------------------------------------------------------------- *///
/// * ======================================== GLOBAL PROJECT VARIABLES ======================================== *///
/// * ------ Declaration of GLOBAL variables which can be seen in all the project modules and libraries -------- *///

/// Technical variables
std::string source_path = "../config/"s; char* sourcepath = const_cast<char*>(source_path.c_str()); // 'source_path' is a path to the directory reading from the 'config/main.ini' file
std::string main_type; // 'mode' from the config/main.ini file: 'LIST' for the execution one by one all the active (ON in the config file) project modules; 'TUTORIAL' as a specific educational mode; 'PERFORMANCE_TEST' as a special test for a computer performance and its ability to work with large PCCs, and the 'TASK' mode, where user-defined task scripts described in separate 'tasks/*.cpp' files are included with all the necessary modules and functions from the project's libraries.

std::vector<std::string> PCCpaths; // The vector containing the paths to all the PCC's matrices, measures and other supplementary data files
std::string source_dir, output_dir; // Input directory (if the initial configuration must be read from file) and output directory for the Writer module and the project log file as it is written in the 'config/main.ini' file
std::string sim_task; // Path to the corresponding 'tasks/*.cpp' file containing a 'simulation task' (for 'TASK' execution mode only) as it is written in the 'config/main.ini' file

/// Global 'log.txt' file output
std::ofstream Out_logfile_stream; // 'Processing_Design.log' file output of the entire computation process as a copy of the console output

/// PCC - related variables
int dim; // Tessellation dimension corresponding to the maximal 'k' in the PCC's k-cells: dim = 1 for graphs and networks, dim = 2 for the 2D plane polygonal tessellations, and dim = 3 for the 3D bulk volumetric tessellations, as it is specified in the 'config/main.ini' file.

/// Combinatorial
std::vector<unsigned int> CellNumbs; // the vector named CellNumbs containing the numbers of k-cells of different types 'k'. It is read from the 'number_of_cells.txt' file of a PCC.
// First line here is the number of nodes (0-cells), second - edges (1-cells), third - faces (2-cells) (in the 2D and 3D cases only), fourth - polyhedra (3-cells) (in the 3D case only)

/// Geometry
std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector; // vectors containing barycenter Cartesian coordinates of the corresponding tessellation's elements
// Global vectors of Cartesian coordinates for: (1) vertex coordinates, (2) barycentres of edges, (3) barycentres of faces and (4) barycentres of polyhedrons

/// Measures
std::vector<double> edge_lengths_vector, face_areas_vector, polyhedron_volumes_vector; // Global vectors of measures: edge lengths, face areas and polyhedra volumes

/// Time interval variables for different parts (modulus) of the CPD code
double Main_time = 0.0, S_time = 0.0, P_time = 0.0, C_time = 0.0, W_time = 0.0;

/// * ===================== MODULES and LIBRARIES ==============================* ///
///* =========================================================================== *///
// * A task for the future: include all the code modules and functions as a single C++ library here placed in the STL directory #include <cpd_lib>

/*! Various useful functions (mostly supplementary) are defined here */
#include "lib/PCC_Support_Functions.h" // It must be here - first in this list of libraries (!)

/*! Various set measures are defined here */
#include "lib/PCC_Measures.h"

/*! An Objects library contains classes of various objects related to the PCC's substructures and k-cells */
#include "lib/PCC_Objects.h"

/*! Section module calculates reduced PCC subcomplexes (parts of the initial PCC including plain cuts) inheriting reduced sequences of special cells and 'state vectors' of the original PCC */
/// #include "lib/PCC_Section/PCC_Subcomplex.h"

/*! Processing module assigned special IDs for the various elements (Nodes, Edges, Faces, Polytopes/Polyhedrons) of the space tessellation */
/* Output: module generates a design_sequences as the lists containing the sequences of k-cells possessing "special" IDs including
 * (1) ASSIGNED: k-Cells, k={0,1,2,3}, corresponding to different generation principles (random, maximum entropy,.. etc.),
 * (2) IMPOSED: m-Cells (where m < k) directly labelled based on the HIGH-ORDER k-Cell IDs
 * (3) INDUCED: i-Cells generated as a result of some KINETIC process. They are always DEPEND on the ASSIGNED design 's_cell_sequences' of special k-Cells and maybe also on the induced 's_induced_cell_sequences' of m-Cells */
#include "lib/PCC_Processing/PCC_Processing.h"

/*! The module provides vectors with the characteristics (entropic, spectral, etc.) representing evolution of state vectors as it given by the PCC_Processing module */
#include "lib/PCC_Characterisation/PCC_Characterisation.h"

/* Writer module performs formatted output of various data structures generated by other modules */
#include "lib/PCC_Writer/PCC_Writer.h"

/*!
 * @brief TUTORIAL :: educational course active in the 'TUTORIAL' CPD code execution mode, please see 'config/main.ini' file.
 * @param initial_configuration
 */
void tutorial(Config &initial_configuration);

/// PCC special structure-related variables :: see class 'Config' in Objects.h and Object.cpp files
Config initial_configuration, configuration; // Initial configuration as an object of the 'Config' class described in Objects.cpp project library

/*!
 * @brief PERFORMANCE_TEST :: testing mode of the CPD code execution, please see 'config/main.ini' file.
 * @param initial_configuration
 */
void performance_test(Config &initial_configuration);

///* ........................................................................................    Main    ................................................................ *///
//* (.h files) * @brief, @param and @return
//* (.cpp files) * @details (detailed descriptions)
/*!
* @details Implement the whole program execution according to the specifications written in 'config/_.ini' files.
* In particular, the LIST mode call the execution of the project modules one by one and compute the execution time for each of them.
* @param void
* @return 0 and the output to console and log file, if successful
*/
int main() {
    cout << "------------------------------------------------------------------------------------------------" << endl;

/// ========== Elapsing time Main =========== ///
    unsigned int Mn_time = clock();
    Main_time = (double) Mn_time;
    cout << endl; Out_logfile_stream << endl;
    cout << "Main execution time before modules is equal to  " << Main_time/ pow(10.0,6.0) <<  "  seconds" << endl;
    Out_logfile_stream << "Main execution time before modules is equal to  " << Main_time/ pow(10.0,6.0) <<  "  seconds" << endl;
    cout << endl;

/// Initial configuration reader and information output to the screen and into the 'Processing_Design.log' file
    initial_configuration.Read_config(); // Read_config() method of the class Config defined in Objects.h and described in Objects.cpp

/// Setting values of the global variables:: all the methods below are in the class Config defined in Objects.h and described in Objects.cpp
    std::vector<int> ConfigVector = initial_configuration.Get_ConfVector();
    dim = initial_configuration.Get_dim();
    source_dir = initial_configuration.Get_source_dir();
    output_dir = initial_configuration.Get_output_dir();
    PCCpaths = initial_configuration.Get_paths();
    main_type = initial_configuration.Get_main_type();
    sim_task = initial_configuration.Get_sim_task();
 /// ============================================================================== ///

/// Global log file output:
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
    Out_logfile_stream << "------------------------------------------------------------------------------------------------" << endl;

/// ========================================================================================================================================== ///
/// ================================================= PERFORMANCE_TEST MODE STARTS HERE ================================================================== ///
/// ========================================================================================================================================== ///
    if (main_type == "PERFORMANCE_TEST"s) { // testing mode 'PERFORMANCE_TEST' which should be further replaced in the 'config/main.ini' file with the 'PERFORMANCE_TEST' and then the 'TASK' or the 'LIST' modes
        cout << "==================================================================================================================================================" << endl; Out_logfile_stream << "==============================================================================================================================================================" << endl;
        cout << "\t\t\t\t\t\t\t\t\t\t[\tStart of the PCC Processing Design code \t]\t\t\t\t\t\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl; Out_logfile_stream << "\t\t\t\t\t[\tStart of the PCC Processing Design code\t]\t\t\t\t\t" << endl << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

        performance_test(initial_configuration); // output the 'performance_test.txt' file to the 'output_dir' showing the relative code execution times of the present server comparing with some reference execution times and suggest the preferable PCC sizes for various simulation tasks
} /// END of the SIMULATION MODE "PERFORMANCE_TEST" as specified in the config/main.ini file

/// ========================================================================================================================================== ///
/// ================================================= TUTORIAL MODE STARTS HERE ================================================================== ///
/// ========================================================================================================================================== ///
    else if (main_type == "TUTORIAL"s) { // TUTORIAL feature to facilitate the first acquaintance with the code: only in the 'TUTORIAL' execution type = the 'mode' variable in the config/main.ini file.
        cout << "==================================================================================================================================================" << endl; Out_logfile_stream << "==============================================================================================================================================================" << endl;
        cout << "\t\t\t\t\t\t\t\t\t\t[\tStart of the PCC Processing Design code \t]\t\t\t\t\t\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl; Out_logfile_stream << "\t\t\t\t\t[\tStart of the PCC Processing Design code\t]\t\t\t\t\t" << endl << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

        tutorial(initial_configuration);
    } /// END of the SIMULATION MODE "TUTORIAL" as specified in the config/main.ini file

/// ========================================================================================================================================== ///
/// ================================================= LIST MODE STARTS HERE ================================================================== ///
/// ========================================================================================================================================== ///
    else if ( main_type == "LIST"s ) { // In the LIST mode all the functions are calling one after another without additional loops and intermediate data output
    /// For all the more complicated simulation cases the TASK mode should be used - see it following next after the 'LIST' module.
        cout << "==================================================================================================================================================" << endl; Out_logfile_stream << "=======================================================================================================================================================================================================================================" << endl;
        cout << "\t\t\t\t\t\t\t\t\t\t[\tStart of the PCC Processing Design code \t]\t\t\t\t\t\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl; Out_logfile_stream << "\t\t\t\t\t[\tStart of the PCC Processing Design code\t]\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

    /// Initialisation of the current_configuration as equal to the initial_configuration
        configuration = initial_configuration;

    /// ====================== I. PCC Subcomplex module ======================
        if (ConfigVector.at(1) == 1) { // if the 'PCC_Section' parameter is switched 'ON' in the config/main.ini file
            cout << "-------------------------------------------------------------------------" << endl;             Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << " START of the PCC Subcomplex module " << endl; Out_logfile_stream << " START of the PCC Subcomplex module " << endl;

            ///* past module here */// Subcomplex(configuration); // Subcomplex(configuration) function

            // ================ Elapsing time for the Subcomplex module ================
            unsigned int Subcomplex_time = clock();
            S_time = (double) Subcomplex_time - Main_time;
            cout << "Section time is equal to  " << S_time/ pow(10.0,6.0) <<  "  seconds" << endl;
            cout << "-------------------------------------------------------" << endl;
            Out_logfile_stream << "Section time is equal to  " << S_time/ pow(10.0,6.0) <<  "  seconds" << endl;
            Out_logfile_stream << "-------------------------------------------------------" << endl;
        } // end if(SectionON)

    /// ====================== II. PCC Processing module ======================
        CellsDesign new_cells_design; // a class described in PCC_Objects.h contained (1) all special k-cell sequences and (2) all the design_<*>_vectors for all k-cells in the PCC

        if (ConfigVector.at(2) == 1) { // if the 'PCC_Processing' parameter is switched 'ON' in the config/main.ini file
            cout << "-------------------------------------------------------------------------" << endl;             Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << "START of the PCC Processing module " << endl; Out_logfile_stream << "START of the PCC Processing module " << endl;

            new_cells_design = PCC_Processing(configuration);

        // ================ Elapsing time for the Processing module ================
            unsigned int Processing_time = clock();
            P_time = (double) Processing_time - S_time - Main_time;
            cout << "Processing time is equal to  " << P_time/ pow(10.0,6.0) <<  "  seconds" << endl << endl; //cout << "-------------------------------------------------------------------------" << endl;
            Out_logfile_stream << "Processing time is equal to  " << P_time/ pow(10.0,6.0) <<  "  seconds" << endl << endl; //Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
        } // end if(ProcessingON)

    /// ====================== III. PCC Characterisation module ======================
        ProcessedComplex pcc_processed;  // a class described in PCC_Objects.h

        if (ConfigVector.at(3) == 1) { // if the 'PCC_Characterisation' parameter is switched 'ON' in the config/main.ini file
            cout << "-------------------------------------------------------------------------" << endl;             Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << "START of the PCC Characterisation module" << endl; Out_logfile_stream << "START of the PCC Characterisation module" << endl;
            cout << "=========================================================================" << endl; Out_logfile_stream << "==============================================================================================================================================================" << endl;

            pcc_processed = PCC_StructureCharacterisation(new_cells_design);

        // ===== Elapsing time for the Characterisation module ================
            unsigned int Characterisation_time = clock();
            C_time = (double) Characterisation_time - S_time - Main_time - P_time;
            cout << "Characterisation time is equal to  " << C_time/ pow(10.0,6.0) <<  "  seconds" << endl << endl; //cout << "-------------------------------------------------------------------------" << endl;
            Out_logfile_stream << "Characterisation time is equal to  " << C_time/ pow(10.0,6.0) <<  "  seconds" << endl << endl; //Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
        }// end if(CharacterisationON)

    /// ====================== IV. PCC Writer module ======================

        if (ConfigVector.at(6) == 1) { // if the 'PCC_Writer' parameter is switched 'ON' in the config/main.ini file
            cout << "-------------------------------------------------------------------------" << endl;             Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << "START of the PCC Writer module" << endl; Out_logfile_stream << "START of the PCC Writer module" << endl;
            cout << "=========================================================================" << endl; Out_logfile_stream << "==============================================================================================================================================================" << endl;

            PCC_Writer(new_cells_design, pcc_processed);

        // ================ Elapsing time for the Writer module ================
            unsigned int Writer_time = clock();
            W_time = (double) Writer_time - Main_time- C_time - S_time - P_time;
            cout << "Writer time is equal to  " << W_time/ pow(10.0,6.0) <<  "  seconds" << endl << endl;
        } // end if(WriterON)

    } /// END of the SIMULATION MODE "LIST" as specified in the config/main.ini file

/// ==========================================================================================================================================
/// ================================================= TASK MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
    else if ( main_type == "TASK"s ) { // In the TASK mode any piece of code using the project libraries can be included

        /// #include.. .cpp

    } /// END of the SIMULATION MODE "TASK" as specified in the config/main.ini file
/// ==========================================================================================================================================
    cout << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl << "\t\t\t\t\t\t\t\t\t\t[\tThe end of the PCC Processing Design\t]\t\t\t\t\t\t\t\t\t\t" << endl << "==================================================================================================================================================" << endl;
    Out_logfile_stream << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl << "\t\t\t\t\t\t\t\t\t\t[\tThe end of the PCC Processing Design\t]\t\t\t\t\t\t\t\t\t\t" << endl << "==================================================================================================================================================" << endl;

/// ================ Total CPD code Elapsing time ================ ///
    unsigned int end_time = clock();
    double fulltime = (double) end_time;
    cout << "Total " << dim << "D " << "runtime of the CPD code is equal to  " << fulltime / pow(10.0, 6.0) << "  seconds" << endl;
    Out_logfile_stream<< "Total " << dim << "D " << "runtime of the CPD code is equal to  " << fulltime / pow(10.0, 6.0) << "  seconds" << endl;
    cout << "-------------------------------------------------------------------------" << endl; Out_logfile_stream << "-------------------------------------------------------------------------" << endl;

    Out_logfile_stream.close(); // closing off-stream 'Processing_Design.log' file

    return 0; // success!
} /// The END of the Main function

/// ========================================= FUNCTIONS DEFINED IN MAIN MODULE =======================================///

/// ====================# 1 #========================== TUTORIAL ==================================================== ///
/*!
 * @details The TUTORIAL code execution mode is designed to make the first acquaintance with the code easier. It is simply a tour around
 * the typical output of the code finishing with the discussion of its user interface ('config/_.ini' files) and of how to
 * replace the TUTORIAL execution mode itself with the PERFORMANCE_TEST and then LIST or TASK execution modes.
 * @param initial_configuration
 */
void tutorial(Config &initial_configuration){ // teaching function
    cout << "Hello there! This is the TEACHING mode of the program execution (!), \n"
    "where the user passes tutorial helping they to understand how to work with the CPD code. \n"
    "It can be changed by replacing the 'mode' variable with 'list' or 'task' \n"
    "instead of the current 'tutorial' contained in the 'simulation_mode' in ../config/main.ini file (!)" << endl;
    cout << endl;
    cout << "Continue? Please type 'Y' for 'yes' or 'N' for 'no' " << endl;
    char if_continue = 'Y'; cin >> if_continue;
    if (if_continue == 'N') exit(0);
} /// END of the tutorial() function

/// ====================# 2 #====================== PERFORMANCE TEST ================================================ ///
/*!
 * @details The PERFORMANCE TEST code execution mode is aimed to allow user estimate the size (number of cells) of the PCCs which they can employ in calculation for different types of the Processing, Characterisation and Design tasks.
 * It writes the analysis data to the 'performance_test.txt' file contained in the 'output_dir' directory showing the relative code execution times of the present computer comparing to some reference execution times.
 * So the user can decide which computations can be done using the current facility for the reasonable time.
 * @param initial_configuration
 */
void performance_test(Config &initial_configuration) { // execution test function

    std::ofstream Out_performance_test_stream; // performance_test output stream to the output file

    /// Initialisation of the current_configuration = initial_configuration
    configuration = initial_configuration;

    // times
    double processing_execution_time = 0.0, full_processing_time = 0.0;
    unsigned int prev_time, new_time;

    int counter_max = 50000; //number of calculation series

    Out_performance_test_stream.open(output_dir + "performance_test.txt"s, ios::trunc); // creates 'performance_test.txt' file in the 'output_dir'
    Out_performance_test_stream.close();

    Out_performance_test_stream.open(output_dir + "performance_test.txt"s, ios::app); // open for writing 'performance_test.txt' file in the 'output_dir'
    Out_performance_test_stream << "------------------------------------------------------------------------------------------------" << endl;

    CellsDesign new_cells_design;
    for (int counter = 0; counter < counter_max; counter++ ) { // TEST LOOP
        prev_time = clock();
        std::vector<int> ConfigVector = configuration.Get_ConfVector();
        dim = configuration.Get_dim();
        source_dir = configuration.Get_source_dir();
        output_dir = configuration.Get_output_dir();
        PCCpaths = configuration.Get_paths();
        /// ---------------------------------------------------------------------- ///
        new_cells_design = PCC_Processing(configuration);

        /// ============= Elapsing time Processing ================ ///
        new_time = clock();
        processing_execution_time = (double) new_time - (double) prev_time;
        full_processing_time += processing_execution_time;
        //prev_time = processing_execution_time;

        // Output to the 'performance_test.txt' file in the 'output_dir
        Out_performance_test_stream << "Processing iteration " << counter + 1 << " tooks  " << processing_execution_time/ pow(10.0,6.0) <<  "  seconds" << endl << endl;

        cout << endl << "Processing iteration " << counter + 1 << " tooks  " << processing_execution_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------------------------" << endl;
        Out_logfile_stream << endl << "Processing iteration " << counter + 1 << " tooks  " << processing_execution_time/ pow(10.0,6.0) <<  "  seconds" << endl << "-------------------------------------------------------------------------" << endl;

    } // end of for (int counter = 0; counter < counter_max; counter++ ) loop

    // Output to the 'performance_test.txt' file in the 'output_dir
    Out_performance_test_stream << endl << "Full execution time for " <<  counter_max << " iterations is equal to  " << full_processing_time/ pow(10.0,6.0) <<  "  seconds" << endl;

    cout << "Full execution time for " <<  counter_max << " iterations is equal to  " << full_processing_time/ pow(10.0,6.0) <<  "  seconds" << endl;
    cout << "-------------------------------------------------------------------------" << endl; Out_logfile_stream << "Processing time is equal to  " << P_time/ pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------------------------" << endl;

    Out_performance_test_stream.close(); // close output to the 'performance_test.txt' file

} /// END of the performance_test() function

                                            /// *** H E A P *** ///