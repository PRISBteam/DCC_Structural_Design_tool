///******************************************************************************************************************************///
///*                                                                                                                           *///
///****************************************************************************************************************************///
///********************************** Dr Elijah Borodin, Manchester, UK, Spring 2022   ***************************************///
///**************************************************************************************************************************///
///************************   Polyhedral Cell Complex (PCC) Processing Design :: (CPD code) (c)   **************************///
///************************************************************************************************************************///
// Material :: quadruple points, triple junctions, grain boundaries, and grains
// Tessellation :: nodes, edges, faces, polytopes
// PCC :: k-cells with their measures (areas, volumes) and barycenter coordinates

///* ----------------------------------------- *
///* Standard C++ (STL) libraries
///* ----------------------------------------- *
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <tuple>
#include <numeric>
#include <algorithm>
// #include <thread> // Require C++ 17 and above
// #include <execution> // Require C++ 17 and above
#include <cmath>

///* ------------------------------------------------------------------------------- *
///* Attached user-defined C++ libraries (must be copied in the directory for STL):
///* ------------------------------------------------------------------------------- *
/// Eigen source: https://eigen.tuxfamily.org/ (2022)
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

/// Spectra source: https://spectralib.org/ (2022)
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>

/// Open MP library https://www.openmp.org/resources/openmp-compilers-tools/
// #include "/usr/local/Cellar/libomp/15.0.7/include/omp.h"
#include <omp.h>

/// Simple reader (library) for .ini files and specific DPD code-related functions for reading its .ini files
#include "lib/ini/ini.h"
#include "lib/ini/ini_readers.h"

// * For the future -- include the whole CPD modules and functions as a single C++ library here placed in the STL directory *//
//include { Processing + SupportFunctions + Objects + Characterisation + Writer + Section + Multiphysics + Kinetic modules }
///#include <CPD_libs>

///------------------------------------------
using namespace std; //Standard namespace
using namespace Eigen; //Eigen namespace
using namespace Spectra; //Spectra namespace

/// Eigen library based Triplets class containing objects in the form T = T(i, j, value), where i and j are element's a(i, j) indices in the corresponding dense matrix and the third variable is its value
typedef Triplet<double> Tr; // <Eigen library class> Declares a triplet's type with name - Tr
typedef SparseMatrix<double> SpMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles with name - SpMat
typedef MatrixXd DMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles with name - SpMat

/// * ------------------------------------------------------------------------------- -------------------------- *
/// * ============================================ GLOBAL VARIABLES ============================================ * ///
/// * ----------- Declaration of GLOBAL variables (can be used in all the project modules and libraries)-------- *
// Technical variables
/// File path to the configuration profile (MainPath is the path to the code's main directory)
string MainPath = "../";
string source_path = MainPath + "config/"s; char* sourcepath = const_cast<char*>(source_path.c_str()); // file path with main.ini

vector<char*> paths; // the vector with paths to all the PCC's matrices, sizes and other elements related to the input and output directories
string source_dir, output_dir; // input and output directories read from the file config_main.txt file
string sim_task; // path to the corresponding .cpp file (with "..") with a simulation task (form TASK main mode) as it is written in the main.ini file

// Global .log file output
ofstream Out_logfile_stream; //log.txt file output for the whole computation process

/// PCC related variables
// General
int dim0, dim; // problem dimension (2D or 3D)
std::vector<unsigned int> CellNumbs; // the vector named CellNumbs with the number of k-cells of different types; will be read from file

// PCC Complex Geometry
// Global vectors of Descartes Coordinates for: (1) vertices, (2) barycenters of edges, (3) barycenters of faces and (4) barycenters of polyhedrons
vector<tuple<double, double, double>> vertex_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector;

// Measures
std::vector<double> edge_lengths_vector, face_areas_vector, polyhedron_volumes_vector;

/// PCC special structure related variables
// Sequences
std::vector<unsigned int> special_n_Sequence, special_e_Sequence, special_f_Sequence, special_p_Sequence;
std::vector<vector<unsigned int>> special_cells_design(4); // is a list of special_<cells>_sequences as an output of the Processing module
// First line here is the sequence of special nodes Sequence_n_vector;
// Second line -- sequence of special edges Sequence_e_vector,
// Third -- sequence of special faces Sequence_f_vector;
// Forth -- sequence of special polyhedrons Sequence_p_vector.

// Configurations and  states
// State_Vector in the form : [Element index] - > [Type] = kind of a code related to the microstructure of PCC
// Configuration_State is a list of all state-vectors: from (1) State_p_vector to (4) State_n_vector
std::vector<int> State_p_vector, State_f_vector, State_e_vector, State_n_vector; // based on each sequences in the list, the State_<*>_vector's of special cells can be calculated, including State_n_vector, State_e_vector, State_f_vector and State_p_vector
std::vector<int> State_pfracture_vector, State_ffracture_vector, State_efracture_vector, State_nfracture_vector; // based on each sequences in the list, the State_<*>_vector's of special cells can be calculated, including State_n_vector, State_e_vector, State_f_vector and State_p_vector
std::vector<vector<int>> Configuration_State, Configuration_cState; //  is the list of all mentioned below State_<*>_vectors as a possible output of the Processing module
/* where n = "nodes", e = "edges", f = "faces" and p = "polyhedrons" */

/// * MODULES and LIBRARIES * ///
// #include all the main project's modules and libraries
// IMPORTANT! Each module here is just a library (*_lib.h) or simple function(s) definition without any additional code surrounded it

/* Various useful functions (it must be here - first in this list! ) */
#include "lib/PCC_SupportFunctions.h"

/* Various set measures */
#include "lib/measures.h"

/* Objects library contains classes of various objects related to PCC substructures and elements (it must be here - second in this list! ) */
#include "lib/PCC_Objects.h"

/* Section module calculates reduced PCC subcomplexes (including plain cuts) inheriting reduced sequences of special cells and State vectors of the original PCC */
//#include "lib/PCC_Section/PCC_Subcomplex.h"

/* Processing module assigned special IDs for the various elements (Vertices, Edges, Faces, Polyhedrons) of the space tessellation */
/* Output: module generates an s_sequences as the lists with the sequences of k-cells possessing "special" IDs including
 * (1) ASSIGNED: k-Cells, k= {0,1,2,3}, corresponding to different generation principles (random, maximum entropy, etc.),
 * (2) IMPOSED m-Cells (where m < k) by HIGH-ORDER k-Cells, and
 * (3) INDUCED i-Cells by some KINETIC process DEPENDENT on the assigned special k-Cells and imposed special m-Cells structures */
#include "lib/PCC_Processing/PCC_Processing.h"

/* Characterisation module calculates various structural characteristics of special substructures defined on the PCC elements */
#include "lib/PCC_Characterisation_/PCC_StructureCharacterisation.h"

/* Multiphysics module assign various physical quantities (energies, temperature, electrical conductivity) to the k-Cells of the PCC's subcomplexes */
//#include "lib/PCC_Multiphysics/PCC_Multiphysics.h"

/* Kinetic module assign some new values for the scalar or vector variables defined on the (Vertices, Edges, Faces, Polyhedrons) of the space tessellation (or new types of identifiers (IDs) of PCC's k-cells) */
/* As its output module currently generates one or several sequences of "fractured" cells' ID's; the new "fractured" State_cvector of faces (2-cells) can be generated */
//#include "lib/PCC_Kinetic/PCC_Kinetic.h"

/* Writer module perform formatted output of various data structures generated by other modules */
#include "lib/PCC_Writer/PCC_Writer.h"

/// PCC special "elements" related to Mutiphysics module
// #1# 0-cells "precipitates"
// #2# 1-cells "disclinations"
// #3# 2-cells "macro-cracks"
// std::vector <macrocrack> large_cracks_vector;
// #4# 3-cells "phases"

///* ........................................................................................    Main    ................................................................ *///
int main() {
    cout << "-------------------------------------------------------------------------------------" << endl;

/// Useful ofstreams cleaning for data output
    // # 1 # LogFile ofstream
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::trunc); // will close at the end of the main (!)
    Out_logfile_stream << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
    string main_type; /// SIMULATION MODE(LIST) or 'SIMULATION MODE(TASK) :: This define the global simulation mode: "LIST" for the "list" of modules implementing one by one (if ON) and "TASK" for the user-defined task scripts with modules and functions from the project's libraries
// ConfigVector (../config/main.ini) contains all the control variables needed for the program
    vector<int> ConfigVector = config_reader_main(source_path, source_dir, output_dir, main_type, Out_logfile_stream);
    dim0 = 3; //kind of a "standard" value
    dim = ConfigVector.at(0); // space dimension of the problem (dim = 2 or 3);

/// Below the file names with the sparse PCC matrices which must already exit in the input_dir and have the same names (!)
// They can be obtained by the PCC Generator tool (https://github.com/PRISBteam/Voronoi_PCC_Analyser) based on the Neper output (*.tess file) (https://neper.info/) or PCC Structure Generator tool (https://github.com/PRISBteam/PCC_Structure_Generator) for plasticity problems
//  PCC matrices
    string ssd0 = source_dir + "A0.txt"s, ssd1 = source_dir + "A1.txt"s, ssd2 = source_dir + "A2.txt"s, ssd3 = source_dir + "A3.txt"s, ssd4 = source_dir + "B1.txt"s, ssd5 = source_dir + "B2.txt"s, ssd6 = source_dir + "B3.txt"s;
//The next line just a technical procedure string to char arrays transformation needed to use them as the function arguments
    paths.push_back(const_cast<char*>(ssd0.c_str())); paths.push_back(const_cast<char*>(ssd1.c_str())); paths.push_back(const_cast<char*>(ssd2.c_str())); paths.push_back(const_cast<char*>(ssd3.c_str())); paths.push_back(const_cast<char*>(ssd4.c_str())); paths.push_back(const_cast<char*>(ssd5.c_str())); paths.push_back(const_cast<char*>(ssd6.c_str()));

//  PCC measures
//    std::vector<double> edge_lengths_vector;
    string ssd7 = source_dir + "face_areas.txt"s;
    paths.push_back(const_cast<char*>(ssd7.c_str()));
    string ssd8 = source_dir + "polyhedron_volumes.txt"s;
    paths.push_back(const_cast<char*>(ssd8.c_str()));

//  PCC geometry (barycenter coordinates)
//    vector<tuple<double, double, double>> vertex_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector;
    string ssd9 = source_dir + "polyhedron_seeds.txt"s; // grain (polyhedron) seeds
    paths.push_back(const_cast<char*>(ssd9.c_str()));
    string ssd10 = source_dir + "vertex_coordinates.txt"s; // vertex coordinates
    paths.push_back(const_cast<char*>(ssd10.c_str()));
    string ssd11 = source_dir + "face_normals.txt"s; // face normal vectors
    paths.push_back(const_cast<char*>(ssd11.c_str()));
///** add here edge_coordinates and face_coordinates

/// Vector with rhe numbers of PCC k-cells for k\in{0,1,2,3} from file
    if (is_file_exists(source_dir + "voro_Ncells.txt"s)) {
        string ncells = source_dir + "voro_Ncells.txt"s;
        char* number_of_cells = const_cast<char *>(ncells.c_str());
// :: CellNumbs vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Polyhedrons number
        CellNumbs = VectorIReader(number_of_cells); // VectorReader is a function from the PCC_SupportFunctions.h library; "number_of_cells" here if the PATH to file

// Special feature for the 2D case ("grains" -> 2-cells; "faces" -> 1-cells; "triple junctions (lines) -> 0-cells" in 2D):
/// WARNING: (Very important!) // It needs to make everything similar in the code for 2D and 3D cases with the same output of the PCC Generator tool (https://github.com/PRISBteam/Voronoi_PCC_Analyser) software (the same set of A0, A1..., B2 matrices)
 //       if (dim == 2) {
 /// it has been changed!           CellNumbs.push_back(CellNumbs.at(2)); CellNumbs.at(2) = CellNumbs.at(1); CellNumbs.at(1) = CellNumbs.at(0); CellNumbs.at(0) = NULL;
 //       }

/// CellNumbs output
cout << "=====================================================================================" << endl;  Out_logfile_stream << "==========================================================================================================================================================================" << endl;
        unsigned int t_length = 0, l_length = 0;
        for (int j : CellNumbs) cout << t_length++ << "-cells #\t" << j << endl;
        for (int j : CellNumbs) Out_logfile_stream << l_length++ << "-cells #\t" << j << endl;

    } // if voro_Ncells exists
    else cout << "ERROR: The file " << source_dir + "voro_Ncells.txt"s << " does not exists (!)" << endl;

/// Initial state::
    if(dim == 2) { // Does not exist in 2D: CellNumbs.at(3 + (dim0 - 3)) !!
        State_p_vector.resize(0,0); State_pfracture_vector.resize(0, 0);
    } else {
        State_p_vector.resize(CellNumbs.at(3 + (dim - 3)), 0); State_pfracture_vector.resize(CellNumbs.at(3 + (dim - 3)), 0);
    }
State_f_vector.resize(CellNumbs.at(2 + (dim0 - 3)), 0); State_e_vector.resize(CellNumbs.at(1 + (dim0 - 3)), 0); State_n_vector.resize(CellNumbs.at(0), 0);
State_ffracture_vector.resize(CellNumbs.at(2 + (dim0 - 3)), 0); State_efracture_vector.resize(CellNumbs.at(1 + (dim0 - 3)), 0); State_nfracture_vector.resize(CellNumbs.at(0), 0);

Configuration_State.push_back(State_n_vector); Configuration_State.push_back(State_e_vector); Configuration_State.push_back(State_f_vector); Configuration_State.push_back(State_p_vector); // the order here matters!
Configuration_cState.push_back(State_nfracture_vector); Configuration_cState.push_back(State_efracture_vector); Configuration_cState.push_back(State_ffracture_vector); Configuration_cState.push_back(State_pfracture_vector); // the order here matters!
// Configuration_State in 3D: [0] - nodes, [1] - edges, [2] - faces, [3] - polyhedrons

/// Output paths.vector to console and logfile out
    int npath = 0, lpath=0;
    cout << "_____________________________________________________________________________________" << endl;
    Out_logfile_stream << "_____________________________________________________________________________________" << endl;
    for (auto m : paths) cout <<"[" << npath++ << "]" << " paths:\t" << m << endl;
    for (auto p : paths) Out_logfile_stream <<"[" << lpath++ << "]" << " paths:\t" << p << endl;

/// ==========================================================================================================================================
/// ================================================= THE LIST MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
// In the LIST mode all the functions are calling one after another without loops
double S_time = 0.0, P_time = 0.0, C_time = 0.0, M_time = 0.0, K_time = 0.0, W_time = 0.0; // time interval variables for different modulus
CellsDesign new_cells_design; // an object (described in PCC_Objects.h) containing (1) all k-cell sequences and (2) all design_x_vectors for all k-cells

    if ( main_type == "LIST"s ) {
// I: PCC_Section.h module
        if (ConfigVector.at(1) == 1) { // if PCC_Section is SWITCH ON in the main.ini file
            cout << " START of the PCC Subcomplex module " << endl;
            Out_logfile_stream << " START of the PCC Subcomplex module " << endl;

     // Subcomplex();

            // ===== Elapsing time Subcomplex ================
            unsigned int Subcomplex_time = clock();
            S_time = (double) Subcomplex_time;
            cout << "Section time is equal to  " << S_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------" << endl;
            Out_logfile_stream << "Section time is equal to  " << S_time/ pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------" << endl;
        } // end if(SectionON)

/// II: PCC_Processing.h module
        if (ConfigVector.at(2) == 1) { // if PCC_Processing is SWITCH ON in the main.ini file
            cout << " START of the PCC Processing module " << endl;
            Out_logfile_stream << " START of the PCC Processing module " << endl;
            new_cells_design = PCC_Processing(Configuration_State, Configuration_cState);
/*
            // update global vectors
            if(new_cells_design.Get_p_sequence().size() > 0) {
                special_p_Sequence = new_cells_design.Get_p_sequence();
                special_cells_design.at(3) = special_p_Sequence;
            }
            if(new_cells_design.Get_f_sequence().size() > 0) {
                special_f_Sequence = new_cells_design.Get_f_sequence();
                special_cells_design.at(2) = special_f_Sequence;
            }
            if(new_cells_design.Get_e_sequence().size() > 0) {
                special_e_Sequence = new_cells_design.Get_e_sequence();
                special_cells_design.at(1) = special_e_Sequence;
            }
            if(new_cells_design.Get_n_sequence().size() > 0) {
                special_n_Sequence = new_cells_design.Get_n_sequence();
                special_cells_design.at(0) = special_n_Sequence;
            }
*/
// ===== Elapsing time Processing ================
            unsigned int Processing_time = clock();
            P_time = (double) Processing_time - S_time;
            cout << "Processing time is equal to  " << P_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------" << endl;
            Out_logfile_stream << "Processing time is equal to  " << P_time/ pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------" << endl;
        } // end if(ProcessingON)

/// III: DCC_Characterisation module
        ProcessedComplex pcc_processed;
        if (ConfigVector.at(3) == 1) {
            cout << "START of the PCC Structure Characterisation module" << endl;
            Out_logfile_stream << "START of the PCC Structure Characterisation module" << endl;

              pcc_processed = PCC_StructureCharacterisation(new_cells_design);

// ===== Elapsing time DCC_Characterisation ================
            unsigned int Characterisation_time = clock();
            C_time = (double) Characterisation_time - S_time - P_time;
            cout << "Characterisation time is equal to  " << C_time/pow(10.0,6.0) <<  "  seconds" << endl;
            cout << "-------------------------------------------------------" << endl;
        }// end if(CharacterisationON)

// IV: PCC_Multiphysics.h module
        if (ConfigVector.at(4) == 1) { // if PCC_Multiphysics is SWITCH ON in the main.ini file
            cout << "START of the PCC Multiphysics module" << endl;
            Out_logfile_stream << "START of the PCC Multiphysics module" << endl;

            // Multiphysics();

            // ===== Elapsing time Multiphysics ================
            unsigned int Multiphysics_time = clock();
            M_time = (double) Multiphysics_time - C_time - S_time - P_time;
            cout << "Section time is equal to  " << M_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------" << endl;
            Out_logfile_stream << "Section time is equal to  " << M_time/ pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------" << endl;
        } // end if(MultiphysicsON)

// V: DCC_Kinetic module
        if (ConfigVector.at(5) == 1) { // if PCC_Kinetic is SWITCH ON in the main.ini file
            cout << "START of the PCC Kinetic module" << endl;
            Out_logfile_stream << "START of the PCC Kinetic module" << endl;

            //kface_sequence = DCC_Kinetic(special_faces_sequence, K_type);

// ===== Elapsing time Kinetic ================
            unsigned int Kinetic_time = clock();
            K_time = (double) Kinetic_time - C_time - S_time - P_time - M_time;
            cout << "Kinetic time is equal to  " << K_time/ pow(10.0,6.0) <<  "  seconds" << endl;
            cout << "-------------------------------------------------------" << endl;
        }// end if(KineticON)

/// VI: PCC_Writing module
        if (ConfigVector.at(6) == 1) { // simply output ON/OFF for the PCC_Writer module on the screen
            cout << "START of the DCC Writer module" << endl;
            PCC_Writer(new_cells_design, pcc_processed);

// ===== Elapsing time Writer ================
            unsigned int Writer_time = clock();
            double W_time = (double) Writer_time - C_time - S_time - P_time - M_time - K_time;
            cout << "Writer time is equal to  " << W_time/ pow(10.0,6.0) <<  "  seconds" << endl;
        } // end if(WriterON)

    } /// end of the SIMULATION MODE "LIST" in the main.ini file
    else if ( main_type == "TASK"s ) {

/// ==========================================================================================================================================
/// ================================================= THE TASK MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
// #include sim_task (??? instead preprocessor related #include - call during the simulation time )

// # 1 # Simple sfaces processing
//#include "tasks/task_sFacesProcessing.cpp"

// # 2 # Processing design for Strips of Inclusions (with the parameters of average \mu and dispersion \sigma of the strip lengths distribution)
//#include "tasks/task_StripsProcessingDesign.cpp"

// # 3 # Macrocrack growth with multiple cracking simulations
//#include "tasks/task_macrocrack.cpp"

// # 4 # Several volumetric phases and imposed GB networks
        { /// *** GRAIN PHASE DESIGN TASK *** ////
cout << " START of the PCC Processing module " << endl;
            //processing
            new_cells_design = PCC_Processing(Configuration_State, Configuration_cState);

cout << " START of the PCC Characterisation  " << endl;
            //characterisation
            double omega = 0, sigma = 0, chi = 0, resistivity = 0;
            double g_Omega(std::vector<int> &grain_states_vector), gb_Sigma(std::vector<int> &face_states_vector), gb_Chi(std::vector<int> &face_states_vector), gb_Resistivity(std::vector<int> &face_states_vector, Eigen::SparseMatrix<double> const &WLMatrix);

            //only polyhedrons
            std::vector<int> grain_states_vector, face_states_vector(CellNumbs.at(2),0);

            if(new_cells_design.Get_p_sequence().size() > 1)
                for (auto pseq : new_cells_design.Get_p_design()) {
                    grain_states_vector.push_back(pseq);
                } // end for(auto pseq...)
//REPAIR for(auto p : grain_states_vector) cout << p << endl;

// new face types
SpMat GFS = SMatrixReader(paths.at(6), CellNumbs.at(2), CellNumbs.at(3)); // Faces-Grains sparse incidence matrix

vector<int> g_couple_types;
    for(unsigned int f = 0; f < CellNumbs.at(2); ++f) { // loop over all Faces
        g_couple_types.clear();
        for (unsigned int g = 0; g < CellNumbs.at(3); ++g) {// loop over all Grains in the PCC
            if (GFS.coeff(f, g) != 0) {
                g_couple_types.push_back(grain_states_vector.at(g));
            }
//REPAIR  cout << " g = " << g << " gsv = " << count grain_states_vector.size() << "  coeff " << GFS.coeff(f, g) << endl;
//REPAIR  cout << " f = " << f << " g = " << g << " g_couple_types size() = " << g_couple_types.size() << endl;
        }
        if (g_couple_types.size() < 2) g_couple_types.push_back(0);

                    if (g_couple_types.at(0) == 0 && g_couple_types.at(1) == 0) face_states_vector.at(f) = 0;
                    else if (g_couple_types.at(0) == 0 && g_couple_types.at(1) == 1 || g_couple_types.at(0) == 1 && g_couple_types.at(1) == 0) face_states_vector.at(f) = 1;
                    else if (g_couple_types.at(0) == 1 && g_couple_types.at(1) == 1) face_states_vector.at(f) = 2;
                    else if (g_couple_types.at(0) == 0 && g_couple_types.at(1) == 2 || g_couple_types.at(0) == 2 && g_couple_types.at(1) == 0) face_states_vector.at(f) = 3;
                    else if (g_couple_types.at(0) == 1 && g_couple_types.at(1) == 2 || g_couple_types.at(0) == 2 && g_couple_types.at(1) == 1) face_states_vector.at(f) = 4;
                    else if (g_couple_types.at(0) == 2 && g_couple_types.at(1) == 2) face_states_vector.at(f) = 5;
                    else if (g_couple_types.at(0) == 0 && g_couple_types.at(1) == 3 || g_couple_types.at(0) == 3 && g_couple_types.at(1) == 0) face_states_vector.at(f) = 6;
                    else if (g_couple_types.at(0) == 1 && g_couple_types.at(1) == 3 || g_couple_types.at(0) == 3 && g_couple_types.at(1) == 1) face_states_vector.at(f) = 7;
                    else if (g_couple_types.at(0) == 2 && g_couple_types.at(1) == 3 || g_couple_types.at(0) == 3 && g_couple_types.at(1) == 2) face_states_vector.at(f) = 8;
                    else if (g_couple_types.at(0) == 3 && g_couple_types.at(1) == 3) face_states_vector.at(f) = 9;
        } // end of for(unsigned int f = 0; f < CellNumbs.at(2); ++f)

        //    omega = g_Omega(grain_states_vector);
        //    sigma = gb_Sigma(face_states_vector);
        //    chi = gb_Chi(face_states_vector);

// Weighted Laplacian Matrix
SpMat WLMatrix = SMatrixReader(paths.at(5), (CellNumbs.at(1)), (CellNumbs.at(2))); //all Edges-Faces
// WLMatrix = FES* ; /// WLMatrix = FES * WEIGHT (U voltages for each type of GB)

/// cout << " befinning of the gb_Resistivity() computations " << endl;
/// resistivity = gb_Resistivity(face_states_vector, WLMatrix);

            cout << " START of the PCC Writer  " << endl;
            //writer
            int output_counter = 1;
            /// grain design
            PCC_CellSequences_Writer(new_cells_design, output_counter);

            /// face_design
            ofstream Out_face_design;
            string face_design_path = output_dir + "face_design.txt"s; // output directory
            Out_face_design.open(face_design_path, ios::trunc);
            for (auto fs: face_states_vector)
                Out_face_design << fs  << " ";
            cout << "(" << ++output_counter << ")  " << "Face design was successfully written in " << face_design_path << endl;
            Out_face_design.close();

            /// additional output
            ofstream Out_grain_measures;
            string grain_measures_path = output_dir + "grain_measures.txt"s; // output directory
            Out_grain_measures.open(grain_measures_path, ios::trunc);

            Out_grain_measures << omega << " " << sigma << " " << chi << " " << resistivity << endl;
            cout << "(" << ++output_counter << ")  " << "Grain measures were successfully written in " << grain_measures_path << endl;
            Out_grain_measures.close();


        } /// *** END of the GRAIN PHASE DESIGN TASK *** ////

    } /// end of the SIMULATION MODE "TASK" in the main.ini file
/// ==========================================================================================================================================

/// ================ Total Elapsing time ================ ///
    unsigned int end_time = clock();
    double fulltime = (double) end_time;
    cout << "-------------------------------------------------------------------------" << endl;
    Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
    cout << dim << "D " << "runtime is equal to  " << fulltime / pow(10.0, 6.0) << "  seconds" << endl;
    Out_logfile_stream << dim << "D " << "runtime is equal to  " << fulltime / pow(10.0, 6.0) << "  seconds" << endl;
    cout << "-------------------------------------------------------------------------" << endl << "\t\t\t\t\t[\tThe end of the PCC Processing\t]\t\t\t\t\t" << endl << "=========================================================================" << endl;
    Out_logfile_stream << "-------------------------------------------------------------------------" << endl << "\t\t\t\t\t[\tThe end of the PCC Processing\t]\t\t\t\t\t" << endl << "=========================================================================" << endl;

/// Closing log_file ofstream
    Out_logfile_stream.close();

    return 0;
} /// The END of Main function

/// ================================== FUNCTIONS DEFINED IN MAIN MODULE ==================================///

/// HEAP
