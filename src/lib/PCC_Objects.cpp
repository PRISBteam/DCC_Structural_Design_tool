/// Classes of special objects related to defect structures on a PCC elements
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace std; // standard namespace
using namespace Eigen; // standard namespace

extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector; // coordinate vectors defined globally
extern std::vector<unsigned int> CellNumbs; // number of cells in a PCC defined globally
extern std::string source_path, output_dir;
std::ofstream Out_local_logstream;

typedef Eigen::SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'

// local libraries
#include "PCC_Support_Functions.h" // It must be here - first in this list (!)

/// Tailored Reader for the main.ini file in the ../config/ subdirectory of the project based on the mINI C++ library (2018 Danijel Durakovic http://pulzed.com/, MIT license is included)
#include "ini/ini_readers.h"

#include "PCC_Objects.h"
///------------------------------------------------------------------


/// --------------------------------------------------------------------------------------------------- ///
/// ========== # 1 # ====================== CLASS 'Config' ============================================ ///
/// --------------------------------------------------------------------------------------------------- ///
/*!
 * @details Read the 'initial configuration' of the problem set in all the relevant '*.ini' files containing in the '\config' project directory using the functions from the 'ini_readers.cpp' project library (and only from there)
 *  Alternatively, Set_config() method allows to set all the values manually
 */
    int Config::Get_dim() const {
        return config_dim;
    }; //!@return dim

    std::string Config::Get_source_dir() const {
        return config_source_dir;
    }; //!@return source_dir

    std::string Config::Get_output_dir() const {
         return config_output_dir;
    }; //!@return output_dir

    std::string Config::Get_pcc_standard() const {
        return pcc_standard;
    }; //!@return pcc_standard

    std::vector<string> Config::Get_paths() const {
        return config_PCCpaths;
    }; //!@return PCCpaths

    std::vector<int> Config::Get_ConfVector() const {
        return config_ConfVector;
    }; //!@return ConfVector

    std::string Config::Get_main_type() const {
         return config_main_type;
    }; //!@return output_dir

    std::string Config::Get_sim_task() const {
        return config_sim_task;
    }; //!@return output_dir

  std::vector<std::vector<unsigned int>> Config::Get_Configuration_sState() const {
      return Configuration_sState;
    }; //!@return Configuration_sState

    std::vector<std::vector<unsigned int>> Config::Get_Configuration_iState() const {
        return Configuration_cState;
    }; //!@return Configuration_sState

     void Config::Read_config(){
        config_ConfVector = config_reader_main(source_path, config_source_dir, config_output_dir, pcc_standard, config_main_type);

        config_dim = config_ConfVector.at(0); // [0] space dimension of the problem (dim = 1, 2 or 3);
        if (config_dim != 1 && config_dim != 2 && config_dim != 3) {
            cout << "Wrong dimension ERROR! Please change 'dim' parameter in ../config/mail.ini file to 1, 2, or 3"
                 << endl;
            Out_local_logstream << "Wrong dimension ERROR! Please change 'dim' parameter in ../config/mail.ini file to 1, 2, or 3"
                    << endl;
            exit(1);
        }

  if (pcc_standard == "pcc1s") {
//  } // end of  if (pcc_standard == "pcc1s")

/// Several file PCCpaths to the sparse PCC's matrices which must already exit in the 'source_dir' and have the same names as below (!)
// They can be obtained by the PCC Generator tool (https://github.com/PRISBteam/Voronoi_PCC_Analyser) based on the Neper output (*.tess file) (https://neper.info/) or PCC Structure Generator tool (https://github.com/PRISBteam/PCC_Structure_Generator) for plasticity problems
//  PCC matrices, metrices and coordinates
        std::string ssd0, ssd1, ssd2, ssd3, ssd4, ssd5, ssd6, ssd7, ssd8, ssd9, ssd10, ssd11, ssd12, ssd13, ssd14; // PCC matrices

        if (is_file_exists(config_source_dir + "algebraic/A0.txt"s) && is_file_exists(config_source_dir + "algebraic/A1.txt"s) &&
            is_file_exists(config_source_dir + "algebraic/B1.txt"s) && config_dim >= 0)
            ssd0 = config_source_dir + "algebraic/A0.txt"s,
            ssd1 = config_source_dir + "algebraic/A1.txt"s,
            ssd4 = config_source_dir + "algebraic/B1.txt"s; // 1D
        else {
            ssd0 = "PCC reading ERROR: there is no A0 file in the PCC directory"s;
            ssd1 = "PCC reading ERROR: there is no A1 file in the PCC directory"s;
            ssd4 = "PCC reading ERROR: there is no B1 file in the PCC directory"s;
        }
        if (is_file_exists(config_source_dir + "/algebraic/A2.txt"s) && is_file_exists(config_source_dir + "algebraic/B2.txt"s) &&
                config_dim >= 2)
            ssd2 = config_source_dir + "algebraic/A2.txt"s,
            ssd5 = config_source_dir + "algebraic/B2.txt"s; // 2D
        else {
            ssd2 = "PCC reading ERROR: there is no A2 file in the PCC directory"s;
            ssd5 = "PCC reading ERROR: there is no B2 file in the PCC directory"s;
        }
        if (is_file_exists(config_source_dir + "/algebraic/A3.txt"s) && is_file_exists(config_source_dir + "algebraic/B3.txt"s) &&
                config_dim == 3)
            ssd3 = config_source_dir + "algebraic/A3.txt"s,
            ssd6 = config_source_dir + "algebraic/B3.txt"s; // 3D
        else {
            ssd3 = "PCC reading ERROR: there is no A3 file in the PCC directory"s;
            ssd6 = "PCC reading ERROR: there is no B3 file in the PCC directory"s;
        }
        if (config_dim > 3 || config_dim < 1) {
            cout
                    << "INPUT DATA ERROR (!) dim > 3 or < 1 as it specified in the ../config/main.ini file. Please, make it equal to 1, 2 or 3."
                    << endl;
            Out_local_logstream
                    << "INPUT DATA ERROR (!) dim > 3 or < 1 in the ../config/main.ini file. Please, make it equal to 1, 2 or 3."
                    << endl;
            exit(1);
        }
        //The next line just a technical procedure string to char arrays transformation needed to use them as the function arguments
        config_PCCpaths.push_back(const_cast<char *>(ssd0.c_str())); // A0
        config_PCCpaths.push_back(const_cast<char *>(ssd1.c_str())); // A1
        config_PCCpaths.push_back(const_cast<char *>(ssd2.c_str())); // A2
        config_PCCpaths.push_back(const_cast<char *>(ssd3.c_str())); // A3
        config_PCCpaths.push_back(const_cast<char *>(ssd4.c_str())); // B1
        config_PCCpaths.push_back(const_cast<char *>(ssd5.c_str())); // B2
        config_PCCpaths.push_back(const_cast<char *>(ssd6.c_str())); // B3

//  PCC measures
/// ---> edge lengths HERE (!) // see next 'ssd14 = source_dir + "edge_lengths.txt"s;' and replace with 'ssd6'
        if (is_file_exists(config_source_dir + "measures/face_areas.txt"s) && config_dim >= 2)
            ssd7 = config_source_dir + "measures/face_areas.txt"s; // 2D face_areas
        else
            ssd7 = "PCC reading ERROR: there is no 'face_areas' file in the PCC directory"s;
        config_PCCpaths.push_back(const_cast<char *>(ssd7.c_str()));
        if (is_file_exists(config_source_dir + "measures/polyhedron_volumes.txt"s) && config_dim == 3)
            ssd8 = config_source_dir + "measures/polyhedron_volumes.txt"s; // 3D polyhedron_volumes
        else
            ssd8 = "PCC reading ERROR: there is no 'polyhedron_volumes' file in the PCC directory"s;
        config_PCCpaths.push_back(const_cast<char *>(ssd8.c_str()));
        //    if (is_file_exists(source_dir + "measures/edge_lengths.txt"s) && dim >= 1) ssd14 = source_dir + "measures/edge_lengths.txt"s; // edge barycentres coordinates      PCCpaths.push_back(const_cast<char *>(ssd14.c_str()));

///  PCC geometry (barycenter coordinates)
// In the for of vector<tuple<double, double, double>> vertex_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector;
        if (is_file_exists(config_source_dir + "coordinates/polyhedron_seeds.txt"s) && config_dim == 3)
            ssd9 = config_source_dir + "coordinates/polyhedron_seeds.txt"s; // grain (polyhedron) seeds
        else
            ssd9 = "PCC reading ERROR: there is no 'polyhedron_seeds' file in the PCC directory"s;
        config_PCCpaths.push_back(const_cast<char *>(ssd9.c_str()));
        if (is_file_exists(config_source_dir + "coordinates/vertex_seeds.txt"s) && config_dim >= 1)
            ssd10 = config_source_dir + "coordinates/vertex_seeds.txt"s; // vertex coordinates
        else
            ssd10 = "PCC reading ERROR: there is no 'vertex_seeds' file in the PCC directory"s;
        config_PCCpaths.push_back(const_cast<char *>(ssd10.c_str()));
        if (is_file_exists(config_source_dir + "other/face_normals.txt"s) && config_dim >= 2)
            ssd11 = config_source_dir + "other/face_normals.txt"s; // face normal vectors
        else
            ssd11 = "PCC reading ERROR: there is no 'face_normals' file in the PCC directory"s;
        config_PCCpaths.push_back(const_cast<char *>(ssd11.c_str()));
        if (is_file_exists(config_source_dir + "coordinates/edge_seeds.txt"s) && config_dim >= 1)
            ssd12 = config_source_dir + "coordinates/edge_seeds.txt"s; // edge barycentres coordinates
        else
            ssd12 = "PCC reading ERROR: there is no 'edge_seeds' file in the PCC directory"s;
        config_PCCpaths.push_back(const_cast<char *>(ssd12.c_str()));
        if (is_file_exists(config_source_dir + "coordinates/face_seeds.txt"s) && config_dim >= 2)
            ssd13 = config_source_dir + "coordinates/face_seeds.txt"s; // face barycentres coordinates
        else
            ssd13 = "PCC reading ERROR: there is no 'face_seeds' file in the PCC directory"s;
        config_PCCpaths.push_back(const_cast<char *>(ssd13.c_str()));

/// Vector with rhe numbers of PCC k-cells for k\in{0,1,2,3} from file
        if (is_file_exists(config_source_dir + "/combinatorial/number_of_cells.txt"s)) {
            std::string ncells = config_source_dir + "/combinatorial/number_of_cells.txt"s;
            char *number_of_cells = const_cast<char *>(ncells.c_str());
            // :: CellNumbs vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Polyhedrons number
            CellNumbs = VectorIReader(
                    number_of_cells); // VectorReader is a function from the PCC_SupportFunctions.h library; "number_of_cells" here if the PATH to file
        } // if voro_Ncells exists
        else cout << "ERROR: The file " << config_source_dir + "/combinatorial/number_of_cells.txt"s
                  << " does not exists (!)" << endl;

/// CellNumbs output
        Out_local_logstream.open(config_output_dir + "Processing_Design.log"s, ios::app); // this *.log stream will be closed at the end of the main function
            cout << "=====================================================================================" << endl;
            Out_local_logstream
                    << "=========================================================================================================================================================================="
                    << endl;
            unsigned int t_length = 0;
            for (int cell_numb: CellNumbs) {
                cout << t_length << "-cells #\t" << cell_numb << endl;
                Out_local_logstream << t_length++ << "-cells #\t" << cell_numb << endl;
            } // end for (int cell_numb : CellNumbs)

        Configuration_sState = {State_p_vector, State_f_vector, State_e_vector, State_n_vector },
        Configuration_cState = {State_pfracture_vector, State_ffracture_vector, State_efracture_vector, State_nfracture_vector }; //  is the list of all mentioned below State_<*>_vectors and State_<*>fracture_vectors as the output of the Processing module // is the list of 'state vectors' analogous to the Configuration_sState but for 'cracked' (or induced) network of k-cells
/// Initial state::
        if (config_dim < 3) { // Does not exist in 2D: CellNumbs.at(3)
            State_p_vector.resize(0, 0);
            State_pfracture_vector.resize(0, 0);
            if (config_dim == 1)
                State_f_vector.resize(0, 0);
            State_ffracture_vector.resize(0, 0);
        } else {
            State_p_vector.resize(CellNumbs.at(3), 0);
            State_pfracture_vector.resize(CellNumbs.at(3), 0);
            State_f_vector.resize(CellNumbs.at(2), 0);
            State_ffracture_vector.resize(CellNumbs.at(2), 0);
        }
        State_e_vector.resize(CellNumbs.at(1), 0);
        State_n_vector.resize(CellNumbs.at(0), 0);
        State_efracture_vector.resize(CellNumbs.at(1), 0);
        State_nfracture_vector.resize(CellNumbs.at(0), 0);

        if (config_dim == 1) { // 1D case
            Configuration_sState.resize(2);
            Configuration_sState = {State_n_vector, State_e_vector};
            Configuration_cState.resize(2);
            Configuration_cState = {State_nfracture_vector, State_efracture_vector};
        }
        if (config_dim == 2) { // 2D case
            Configuration_sState.resize(3);
            Configuration_sState = {State_n_vector, State_e_vector, State_f_vector};
            Configuration_cState.resize(3);
            Configuration_cState = {State_nfracture_vector, State_efracture_vector, State_ffracture_vector};
        } else { // 3D case
            Configuration_sState.resize(4);
            Configuration_sState = {State_n_vector, State_e_vector, State_f_vector, State_p_vector};
            Configuration_cState.resize(4);
            Configuration_cState = {State_nfracture_vector, State_efracture_vector, State_ffracture_vector,
                                    State_pfracture_vector};
        } // Configuration_sState in 3D: [0] -> nodes, [1] -> edges, [2] -> faces, [3] -> polyhedrons

/// Output PCCpaths.vector to console and logfile out
        int npath = 0;
        cout << "_____________________________________________________________________________________" << endl;
        Out_local_logstream
                << "_____________________________________________________________________________________" << endl;
        for (auto path: config_PCCpaths) {
            cout << "[" << npath << "]" << " PCCpaths:\t" << path << endl;
            Out_local_logstream << "[" << npath++ << "]" << " PCCpaths:\t" << path << endl;
        }
        cout << endl; Out_local_logstream << endl;
  } // end of  if (pcc_standard == "pcc1s")
  else {
      cout << "ERROR in the reading PCC: please specify the correct 'pcc_standard' perameter in the config/main.ini file corresponding to the version of the PCC you use (please see technical documentation for more details, the first PCC standard has an ID 'pcc1s'" << endl;
      exit(1);
  }
        Out_local_logstream.close();
    }; // Read the 'initial configuration' of the problem set in all the relevant '_.ini' files containing in the '\config' project directory using the functions from the 'ini_readers.cpp' project library (and only from there)

    /// --------------------------------------- *** END of void Config::Read_config() method *** ------------------------------------------------ ///

    void Config::Set_config(const std::vector<int> &ConfigVector, const std::string &source_dir, int &dim, std::vector<char*> paths, std::vector<vector<int>> Configuration_State, std::vector<vector<int>> Configuration_cState){

    }; // manual setting of the configuration

// ========== END of class CONFIG functions description

/// --------------------------------------------------------------------------------------------------- ///
/// ========== # 2 # ====================== CLASS 'CellDesign' ============================================ ///
/// --------------------------------------------------------------------------------------------------- ///
/*!
 * @details Set a list of state_vectors corresponding to different dimensions 'k' of the k-cells in a PCC
 * @param pdesign
 * @param fdesign
 * @param edesign
 * @param ndesign
 */
void CellDesign::Set_designes(std::vector<unsigned int> pdesign, std::vector<unsigned int> fdesign, std::vector<unsigned int> edesign, std::vector<unsigned int> ndesign){
    p_special_design = pdesign; f_special_design = fdesign; e_special_design = edesign; n_special_design = ndesign;
}

void CellDesign::Set_induced_designs(std::vector<unsigned int> p_ind_design, std::vector<unsigned int> f_ind_design, std::vector<unsigned int> e_ind_design, std::vector<unsigned int> n_ind_design){
    p_induced_design = p_ind_design; f_induced_design = f_ind_design; e_induced_design = e_ind_design; n_induced_design = n_ind_design;
}

void CellDesign::Set_special_sequences(std::vector<unsigned int> psequence, std::vector<unsigned int> fsequence, std::vector<unsigned int> esequence, std::vector<unsigned int> nsequence){
    p_special_sequence = psequence; f_special_sequence = fsequence; e_special_sequence = esequence; n_special_sequence = nsequence;
    }

    void CellDesign::Set_induced_sequences(std::vector<unsigned int> p_ind_sequence, std::vector<unsigned int> f_ind_sequence, std::vector<unsigned int> e_ind_sequence, std::vector<unsigned int> n_ind_sequence){
        p_induced_sequence = p_ind_sequence; f_induced_sequence = f_ind_sequence; e_induced_sequence = e_ind_sequence; n_induced_sequence = n_ind_sequence;
    }

    void CellDesign::Set_special_sequence(std::vector<unsigned int> special_x_sequence, int cell_type){
        switch (cell_type) {
            case 3:
                p_special_sequence = special_x_sequence;
                break;
            case 2:
                f_special_sequence = special_x_sequence;
                break;
            case 1:
                e_special_sequence = special_x_sequence;
                break;
            case 0:
                n_special_sequence = special_x_sequence;
                break;
        } // end switch(cell_type)
    } // End of Set_special_sequence()

void CellDesign::Set_agglomeration_sequence(std::vector<Agglomeration> &agglomeration_x_sequence, int cell_type) {
    switch (cell_type) {
        case 3:
            p_agglomerations_map = agglomeration_x_sequence;
            break;
        case 2:
            f_agglomerations_map = agglomeration_x_sequence;
            break;
        case 1:
            e_agglomerations_map = agglomeration_x_sequence;
            break;
        case 0:
            n_agglomerations_map = agglomeration_x_sequence;
            break;
    } // end switch(cell_type)
} // End of Set_agglomeration_sequence()


void CellDesign::Set_induced_sequence(std::vector<unsigned int> induced_x_sequence, int cell_type){
        switch (cell_type) {
            case 3:
                p_induced_sequence = induced_x_sequence;
                break;
            case 2:
                f_induced_sequence = induced_x_sequence;
                break;
            case 1:
                e_induced_sequence = induced_x_sequence;
                break;
            case 0:
                n_induced_sequence = induced_x_sequence;
                break;
        } // end switch(cell_type)
    } // End of Set_induced_sequence()

    void CellDesign::Set_special_configuration(std::vector<unsigned int> special_x_configuration, int cell_type){
        switch (cell_type) {
            case 3:
                p_special_design = special_x_configuration;
                break;
            case 2:
                f_special_design = special_x_configuration;
                break;
            case 1:
                e_special_design = special_x_configuration;
                break;
            case 0:
                n_special_design = special_x_configuration;
                break;
        } // end switch(cell_type)
    } // End Set_special_configuration()

    void CellDesign::Set_special_series(std::vector<std::vector<unsigned int>> special_x_series, int cell_type){
        switch (cell_type) {
            case 3:
                p_special_series = special_x_series;
                break;
            case 2:
                f_special_series = special_x_series;
                break;
            case 1:
                e_special_series = special_x_series;
                break;
            case 0:
                n_special_series = special_x_series;
                break;
        } // end switch(cell_type)
    } // End Set_special_series()

    void CellDesign::Set_induced_design(std::vector<unsigned int> induced_x_design, int cell_type){
        switch (cell_type) {
            case 3:
                p_induced_design = induced_x_design;
                break;
            case 2:
                f_induced_design = induced_x_design;
                break;
            case 1:
                e_induced_design = induced_x_design;
                break;
            case 0:
                n_induced_design = induced_x_design;
                break;
        } // end switch(cell_type)
    } // End Set_induced_design()

void CellDesign::Set_induced_series(std::vector<std::vector<unsigned int>> induced_x_series, int cell_type){
    switch (cell_type) {
        case 3:
            p_induced_series = induced_x_series;
            break;
        case 2:
            f_induced_series = induced_x_series;
            break;
        case 1:
            e_induced_series = induced_x_series;
            break;
        case 0:
            n_induced_series = induced_x_series;
            break;
    } // end switch(ctype)
} // END Set_induced_series()

/// Get
    std::vector<unsigned int> CellDesign::Get_p_special_sequence(void) const {
        if (p_special_sequence.size() == 0) {
            cout << "WARNING: p_special_sequence did not set!" << endl;
            return {0};
        }
        else return p_special_sequence;
    }
    std::vector<unsigned int> CellDesign::Get_f_special_sequence(void) const {
        if (f_special_sequence.size() == 0) {
            cout << "WARNING: f_special_sequence did not set!" << endl;
            return {0};
        }
        else return f_special_sequence;
    }
    std::vector<unsigned int> CellDesign::Get_e_special_sequence(void) const {
        if (e_special_sequence.size() == 0) {
            cout << "WARNING: e_special_sequence did not set!" << endl;
            return {0};
        }
        else return e_special_sequence;
    }
    std::vector<unsigned int> CellDesign::Get_n_special_sequence(void) const {
        if (n_special_sequence.size() == 0) {
            cout << "WARNING: n_special_sequence did not set!" << endl;
            return {0};
        }
        else return n_special_sequence;
    }

    // Agglomerations
    std::vector<Agglomeration> CellDesign::Get_p_agglomeration_map(void) const {
        if (p_agglomerations_map.size() == 0) {
            cout << "WARNING: p_agglomerations_map did not set!" << endl;
            return {0};
        }
        else return p_agglomerations_map;
    }
std::vector<Agglomeration> CellDesign::Get_f_agglomeration_map(void) const {
    if (f_agglomerations_map.size() == 0) {
        cout << "WARNING: f_agglomerations_map did not set!" << endl;
        return {0};
    }
    else return f_agglomerations_map;
}
std::vector<Agglomeration> CellDesign::Get_e_agglomeration_map(void) const {
    if (e_agglomerations_map.size() == 0) {
        cout << "WARNING: e_agglomerations_map did not set!" << endl;
        return {0};
    }
    else return e_agglomerations_map;
}
std::vector<Agglomeration> CellDesign::Get_n_agglomeration_map(void) const {
    if (n_agglomerations_map.size() == 0) {
        cout << "WARNING: n_agglomerations_map did not set!" << endl;
        return {0};
    }
    else return n_agglomerations_map;
}
    // Special series

std::vector<std::vector<unsigned int>> CellDesign::Get_p_special_series(void) const {
    if (p_special_series.size() == 0) {
        cout << "WARNING: p_special_series did not set!" << endl;
        return {{0}};
    }
    else return p_special_series;
    }
std::vector<std::vector<unsigned int>> CellDesign::Get_f_special_series(void) const {
    if (f_special_series.size() == 0) {
        cout << "WARNING: f_special_series did not set!" << endl;
        return {{0}};
    }
    else return f_special_series;
    }
std::vector<std::vector<unsigned int>> CellDesign::Get_e_special_series(void) const {
    if (e_special_series.size() == 0) {
        cout << "WARNING: n_special_sequence did not set!" << endl;
        return {{0}};
    }
    else return e_special_series;
    }
std::vector<std::vector<unsigned int>> CellDesign::Get_n_special_series(void) const {
    if (n_special_series.size() == 0) {
        cout << "WARNING: n_special_sequence did not set!" << endl;
        return {{0}};
    }
    else return n_special_series;
    }

std::vector<unsigned int> CellDesign::Get_p_induced_sequence(void) const {
        if (p_induced_sequence.size() == 0) {
            cout << "WARNING: polyhedron induced sequence did not set!" << endl;
            return {0};
        }
        else return p_induced_sequence;
    }
    std::vector<unsigned int> CellDesign::Get_f_induced_sequence(void) const {
        if (f_induced_sequence.size() == 0) {
            cout << "WARNING: face induced sequence did not set!" << endl;
            return {0};
        }
        else return f_induced_sequence;
    }
    std::vector<unsigned int> CellDesign::Get_e_induced_sequence(void) const {
        if (e_induced_sequence.size() == 0) {
            cout << "WARNING: edge induced sequence did not set!" << endl;
            return {0};
        }
        else return e_induced_sequence;
    }
    std::vector<unsigned int> CellDesign::Get_n_induced_sequence(void) const {
        if (n_induced_sequence.size() == 0) {
            cout << "WARNING: node induced sequence did not set!" << endl;
            return {0};
        }
        else return n_induced_sequence;
    }
    std::vector<unsigned int> CellDesign::Get_p_design(void) const {
        if (p_special_design.size() == 0) {
            cout << "WARNING: p_special_design did not set!" << endl;
            return {0};
        }
        else return p_special_design;
    }
    std::vector<unsigned int> CellDesign::Get_f_design(void) const {
        if (f_special_design.size() == 0) {
            cout << "WARNING: f_special_design did not set!" << endl;
            return {0};
        } else return f_special_design;
    }
    std::vector<unsigned int> CellDesign::Get_e_design(void) const {
        if (e_special_design.size() == 0) {
            cout << "WARNING: e_special_design did not set!" << endl;
            return {0};
        } else return e_special_design;
    }
    std::vector<unsigned int> CellDesign::Get_n_design(void) const {
        if (n_special_design.size() == 0) {
            cout << "WARNING: n_special_design did not set!" << endl;
            return {0};
        }
        else return n_special_design;
    }
// ========== END of the class CELLS_DESIGN functions description

/// --------------------------------------------------------------------------------------------------- ///
/// ========== # 3 # ====================== CLASS 'Subcomplex' ============================================ ///
/// --------------------------------------------------------------------------------------------------- ///

Subcomplex::Subcomplex(unsigned int subcomplex_id_new) { // constructor 1, simple
    subcomplex_id = subcomplex_id_new;
}
    //2
    Subcomplex::Subcomplex(unsigned int subcomplex_id_new, std::vector <unsigned int> new_sub_grains_sequence) { // constructor 2, based on a sub_grains_sequence
        subcomplex_id = subcomplex_id_new;
        Set_grains_sequence(new_sub_grains_sequence);
    }

    /// Grains
    void Subcomplex::Set_grains_sequence(std::vector <unsigned int> new_sub_grains_sequence){
        sub_grains_sequence = new_sub_grains_sequence;
    }

    std::vector <unsigned int> Subcomplex::Get_grains_sequence(unsigned int subcomplex_id) const {
        if(sub_grains_sequence.size() != 0)
            return sub_grains_sequence;
        else return {0};
    }

    /// geometry
    void Subcomplex::Set_sub_grain_coordinates(std::vector<tuple<double, double, double>> new_sub_grain_coordinates){
        sub_grain_coordinates = new_sub_grain_coordinates;
    }
    std::vector<std::tuple<double, double, double>> Subcomplex::Get_sub_grain_coordinates(unsigned int subcomplex_id) const {
        return sub_grain_coordinates;
    }

    //2
//    std::vector <unsigned int>  Get_grains_sequence(unsigned int subcomplex_id, std::vector <unsigned int> new_sub_grains_sequence){
//        if(sub_grains_sequence.size() != 0) return sub_grains_sequence;
//        else { Set_grains_sequence(new_sub_grains_sequence); return sub_grains_sequence; } }

    /// Faces
    void Subcomplex::Set_faces_sequence(std::vector <unsigned int> new_sub_faces_sequence){
        sub_faces_sequence = new_sub_faces_sequence;
    }
    //1
    std::vector <unsigned int>  Subcomplex::Get_faces_sequence(unsigned int  subcomplex_id) const{
        if(sub_faces_sequence.size() != 0)
            return sub_faces_sequence;
        else return {0};
    }

    void Subcomplex::Set_common_faces_sequence(std::vector <unsigned int> new_common_faces_sequence){
        common_faces_sequence = new_common_faces_sequence; }
    std::vector <unsigned int> Subcomplex::Get_common_faces_sequence(unsigned int subcomplex_id) const {
        return common_faces_sequence; }

    void Subcomplex::Set_sfaces_sequence(std::vector <unsigned int> const &ssub_faces_sequence){
        s_sub_faces_sequence = ssub_faces_sequence;
    }
    std::vector <unsigned int> Subcomplex::Get_sfaces_sequence(unsigned int  subcomplex_id) const {
        if(s_sub_faces_sequence.size() > 0) return s_sub_faces_sequence;
        else {
            cout << "Caution! s_sub_faces_sequence = 0! Please add it to the corresponding subcomplex/crack"s << endl;
            return {0};
        }
    }

    void Subcomplex::Set_cfaces_sequence(std::vector <unsigned int> csub_faces_sequence){
        c_sub_faces_sequence = csub_faces_sequence;
    }
    std::vector <unsigned int> Subcomplex::Get_cfaces_sequence(unsigned int  subcomplex_id) const {
        return c_sub_faces_sequence;
    }

    /// geometry
    void Subcomplex::Set_common_faces_coordinates(std::vector<tuple<double, double, double>> new_common_faces_coordinates){
        common_faces_coordinates = new_common_faces_coordinates; }

    std::vector<tuple<double, double, double>> Subcomplex::Get_common_faces_coordinates(unsigned int subcomplex_id) const {
        return common_faces_coordinates; }

// ========== END of class SUBCOMPLEX functions description

/// --------------------------------------------------------------------------------------------------- ///
/// ========== # 4 # ====================== CLASS 'ProcessedComplex' ============================================ ///
/// --------------------------------------------------------------------------------------------------- ///

void ProcessedComplex::Set_design(CellDesign processed_pcc_design) {
        pcc_design = processed_pcc_design;
    }

// ========== END of the class PROCESSED_COMPLEX functions description

/// --------------------------------------------------------------------------------------------------- ///
/// ========== # 4 # ====================== CLASS 'Agglomeration' ============================================ ///
/// --------------------------------------------------------------------------------------------------- ///

//void Set_new_agglomeration(unsigned int AFace);
//void Set_agglomeration_type(std::string type);
//void Set_agglomeration_power(vector<vector<int>> const &RW_series_vector);
//void SetAvLength(vector<vector<int>> const &RW_series_vector); // Average length of strips related to this agglomeration

//unsigned int Get_face_agglomeration() const;
//int Get_agglomeration_power() const;
//int GetAPower(vector<vector<int>> const &RW_series_vector); /// overloaded /// BAD
//int GetAvLength() const; /// BAD
//int GetAvLength(vector<vector<int>> const &RW_series_vector)  /// overloaded /// BAD

Agglomeration::Agglomeration(unsigned int AFace) { // constructor 1 simple
        aface_number = AFace;
    }

Agglomeration::Agglomeration(unsigned int AFace, unsigned int AglPower) { // constructor 2 complex
        aface_number = AFace;
        apower = AglPower;
    }

    void Agglomeration::Set_new_agglomeration(unsigned int AFace) {
        aface_number = AFace;
    }

    void Agglomeration::Set_agglomeration_type(std::string type) {
        atype = type;
        if (atype == "rgo") {
            surface_energy = 0.2; // units [J/m^2]
            adhesion_energy = 0.4; // units [J/m^2]
        }
    }

    void Agglomeration::Set_agglomeration_power(vector<vector<unsigned int>> const &RW_series_vector) { /// Agglomerations power
        int AglPower = 0;
        //vector<int> aggl_vector(CellNumbs.at(2), 0); // state vector for agglomerations (with the size # cells initially filled with all 0s) : contains # of faces with agglomerations
        for (auto RWsfv : RW_series_vector)
            for (auto RWsf: RWsfv)
                if(RWsf == aface_number) AglPower += 1;

        apower = AglPower;
    }

    void Agglomeration::SetAvLength(vector<vector<unsigned int>> const &RW_series_vector) { // Average length of strips related to this agglomeration

        int ATotalLength = 0;
        for (auto RWsfv : RW_series_vector)
            for (auto RWsf: RWsfv)
                if(RWsf == aface_number) ATotalLength += RWsfv.size();

        a_average_strip_length = ATotalLength/ (double) Get_agglomeration_power(RW_series_vector);
    }

    unsigned int Agglomeration::Get_agglomeration_kcell_number() const {
        return aface_number;
    }

    int Agglomeration::Get_agglomeration_power() const {
        return apower;
    }

    int Agglomeration::Get_agglomeration_power(vector<vector<unsigned int>> const &RW_series_vector) {
        if (apower != 0) {
            return apower;
        } else
        {
            unsigned int AglPower = 0;
            for (auto RWsfv : RW_series_vector)
                for (auto RWsf: RWsfv)
                    if(RWsf == aface_number) AglPower += 1;

            apower = AglPower;
            return apower;
        }
    }

    int Agglomeration::GetAvLength() const {
        return a_average_strip_length;
    }

    int Agglomeration::GetAvLength(vector<vector<unsigned int>> const &RW_series_vector) {
        if (a_average_strip_length != 0) {
            return a_average_strip_length;
        } else
        {
            int ATotalLength = 0;
            for (auto RWsfv: RW_series_vector)
                for (auto RWsf: RWsfv)
                    if (RWsf == aface_number) ATotalLength += RWsfv.size();

            a_average_strip_length = ATotalLength / (double) Get_agglomeration_power(RW_series_vector);
            return a_average_strip_length;
        }
    }

// ========== END of the class AGGLOMERATION functions description


/// # 6 # The class of Polytopes in a PCC

    Polytope::Polytope(unsigned int grain_new_id) { // constructor 1 simple
        grain_id = grain_new_id;
    }

    void Polytope::Set_node_ids(unsigned int grain_id, SpMat const &GFS, SpMat const &FES, SpMat const &ENS) { // set the node ids
        /// GFS -> FES -> ENS
        if(node_ids.size() == 0) {
            for(unsigned int l = 0; l < CellNumbs.at(2); l++) {// over all Faces (l)
                if (GFS.coeff(l, grain_id) == 1) {
                    for (unsigned int j = 0; j < CellNumbs.at(1); j++) // over all Edges (j)
                        if (FES.coeff(j, l) == 1) { // at the chosen Face with ID = 'l'
                            for (unsigned int i = 0; i < CellNumbs.at(0); i++) // over all Nodes
                                if (ENS.coeff(i, j) == 1) node_ids.push_back(i); // at the chosen Face with ID = 'l'
                        } // end of if (FES.coeff(l, j) == 1)
                } // end of (GFS.coeff(m, l) == 1)
            } // end of for(unsigned int l = 0; l < CellNumbs.at(2); l++) - Faces

        }/// end of if(node_ids.size() == 0)

    } // end of Set_node_ids()

    void Polytope::Set_Faces_list(unsigned int grain_id, SpMat const &GFS) {
        for (unsigned int l = 0; l < CellNumbs.at(2); l++) // for each GB
            if (GFS.coeff(l, grain_id) == 1)
                Faces_list.push_back(l);
    } // end of Set_GBs_list()

    vector<unsigned int> Polytope::Get_Faces_list() const {
        if (Faces_list.size() > 0) return Faces_list;
        else { cout << "coution GBs_list.size() = 0! Please Set_GBs_list(unsigned int grain_id, SpMat const &GFS)  first!"s << endl;
            return {0};
        };
    } // end of Get_GBs_list()

    /// return - vector of all node (vertices) coordinates of a grain
    void Polytope::Set_node_coordinates(unsigned int grain_id) { // set the node ids from Tr = triplet list
//        for (auto  itr = node_ids.begin(); itr != node_ids.end(); ++itr)
        //if(find(node_ids.begin(), node_ids.end(), distance(node_ids.begin(), itr)) != node_ids.end())
//                node_coordinates.push_back(vertex_coordinates_vector.at(*itr)); // vector<unsigned int> node_ids;
        if(node_ids.size() > 0) {
            for (auto ids: node_ids) {
//REPAIR        cout << "vertex_coordinates_vector size " << vertex_coordinates_vector.size() << endl;
//REPAIR        cout << "ids: " << ids << " node_coordinates size: " << get<0>(vertex_coordinates_vector.at(ids)) << endl;
                node_coordinates.push_back(node_coordinates_vector.at(ids)); // vector<unsigned int> node_ids;
//REPAIR cout << "grain id: " << grain_id << " node_coordinates size: " << node_coordinates.size() << endl;
            }

        }
        else {
            cout << "Caution! The size of node_ids vector of the grain " << grain_id << " is 0 !" << endl;
            node_coordinates.push_back(make_tuple(0,0,0)); /// change after !!
        }
    } // the end of Set_node_coordinates

    std::vector<unsigned int> Polytope::Get_node_ids(unsigned int grain_id) const { // set the node ids
        return node_ids;
    }

std::vector<std::tuple<double, double, double>> Polytope::Get_node_coordinates(unsigned int grain_id) const { // set the node ids
        if (node_coordinates.size() != 0) {
            return node_coordinates;
        } else {
            throw std::invalid_argument( "Please call Set_node_coordinates(unsigned int grain_id, vector<tuple<double, double, double>> const &vertex_coordinates) method first!");
            return node_coordinates;
        }
    } // end of  Get_node_coordinates() method

    /// return - vector with two tuples : { x_min, y_min, z_min; x_max, y_max, z_max} of a grain witn number grain_id
    vector<tuple<double, double, double>> Polytope::Get_minmax_node_coordinates(unsigned int grain_id) const { // min and max {x,y,z} values of vertices for a grain
        vector<tuple<double, double, double>> minmax_tuple;
        ///Get_node_coordinates(grain_id)
        vector<tuple<double, double, double>> tup_node_coordinates = Get_node_coordinates(grain_id); // class Grain3D function Get_node_coordinates(grain_id)
//REPAIR        cout << "Xtup_node_coordinates " << get<0>(tup_node_coordinates.at(0)) << " Ytup_node_coordinates " <<get<1>(tup_node_coordinates.at(0)) << " Ztup_node_coordinates " << get<2>(tup_node_coordinates.at(0)) << endl;

        // separating in three parts
        vector<double> x_node_coordinates, y_node_coordinates, z_node_coordinates;
        for (auto itr = tup_node_coordinates.begin(); itr != tup_node_coordinates.end(); ++itr) {
            x_node_coordinates.push_back(get<0>(*itr));
            y_node_coordinates.push_back(get<1>(*itr));
            z_node_coordinates.push_back(get<2>(*itr));
        }
//REPAIR        cout << "x_node_coordinates " << (x_node_coordinates.at(0)) << " y_node_coordinates " << (y_node_coordinates.at(0)) << " z_node_coordinates " << (z_node_coordinates.at(0)) << endl;
//        std::for_each(tup_node_coordinates.begin(), tup_node_coordinates.end(), );

// Return iterators (!)
        auto itx = std::minmax_element(x_node_coordinates.begin(), x_node_coordinates.end());
        double xmin = *itx.first;
        double xmax = *itx.second;

        auto ity = std::minmax_element(y_node_coordinates.begin(), y_node_coordinates.end());
        double ymin = *ity.first;
        double ymax = *ity.second;

        auto itz = std::minmax_element(z_node_coordinates.begin(), z_node_coordinates.end());
        double zmin = *itz.first;
        double zmax = *itz.second;

        minmax_tuple = {make_tuple(xmin, ymin, zmin), make_tuple(xmax, ymax, zmax)};
//REPAIR !!        cout << "Xmin " << get<0>(minmax_tuple.at(0)) << " Ymin " <<get<1>(minmax_tuple.at(0)) << " Zmin " << get<2>(minmax_tuple.at(0)) << endl;

        return minmax_tuple;
    } // end of Get_minmax_node_coordinates()

/// ========== END of class grain3D functions description

/// # 6 # The class of a MACROCRACK

Macrocrack::Macrocrack(int crack_id_new, Subcomplex &half_plane_sub) : half_plane_subcomplex(0) { // constructor 3, with a subcomplex
        crack_id = crack_id_new;
        crack_length = half_plane_sub.sub_length; // crack length from the corresponding subcomplex size
        half_plane_subcomplex = half_plane_sub; //set subcomplex
    }

    double Macrocrack::Get_crack_length(int crack_id_new) const {
        return crack_length;
    }

    void Macrocrack::Set_real_crack_length(double sample_size) {
        real_crack_length = half_plane_subcomplex.sub_length * sample_size;
    }
    double Macrocrack::Get_real_crack_length() const {
        if (real_crack_length > 0.0) return real_crack_length;
        else {
            cout << "Caution! real_crack_length = 0! Please use Set_real_crack_length(double sample_size) before!"s << endl;
            return 0;
        }
    }

    void Macrocrack::Set_crack_plane() {
        crack_plane = {half_plane_subcomplex.a_n, half_plane_subcomplex.b_n, half_plane_subcomplex.c_n, half_plane_subcomplex.D_plane};
    }
    vector<double> Macrocrack::Get_crack_plane() const {
        if (crack_plane.size() > 0.0) return crack_plane;
        else {
            cout << "Caution! crack_plane.size() = 0! Please update first {half_plane_subcomplex.a_n, half_plane_subcomplex.b_n, half_plane_subcomplex.c_n, half_plane_subcomplex.D_plane} in the corresponding subcomplex!"s << endl;
            return {0};
        }
    }

    void Macrocrack::Set_multiple_cracking_energy(double total_energy) {
        multiple_cracking_energy = total_energy;
    }
    double Macrocrack::Get_multiple_cracking_energy() const {
        return multiple_cracking_energy;
    }

    std::vector <unsigned int> Macrocrack::Get_faces_sequence() const {
        return half_plane_subcomplex.Get_faces_sequence(crack_id); }

    std::vector <unsigned int> Macrocrack::Get_sfaces_sequence() const {
        return half_plane_subcomplex.Get_sfaces_sequence(0); }

    std::vector <tuple<double,double,double>> Macrocrack::Get_common_faces_coordinates(unsigned int  crack_id) const {
        return half_plane_subcomplex.Get_common_faces_coordinates(crack_id); }

/// ========== END of class MACROCRACK functions description

/// # V # The class of a CELLS_ENERGIES :: list of the energy_vectors corresponding to different dimensions 'k' of the k-cells in a PCC
//class CellEnergies
/// CellEnergies

/// ========== END of class CellEnergies functions description


/*
 /// # 0 # The class of Grain Boundaries in a PCC

class grain_boundary{

public:
    unsigned int GB_id; // grain boundary ID

    grain_boundary(unsigned int GB_new_id) { // class simple constructor
        GB_id = GB_new_id;
    }

/// GB state
    bool is_inclusion; // 0 - no, 1 - yes
    bool is_fractured; // 0 - no, 1 - yes
    bool is_agglomeration; // 0 - no, 1 - yes

/// Set values methods

    void Set_surface_energy(vector<double> GB_SE_vector ){
        surface_energy = GB_SE_vector.at(GB_id);
    }

    void Set_external_elastic_energy(vector<double> GB_EEE_vector) { /// EEE: external elastic energy vector
        external_elastic_energy = GB_EEE_vector.at(GB_id);
    }

    void Set_crack_interaction_energy(vector<double> GB_CIE_vector) {
        crack_interaction_energy = GB_CIE_vector.at(GB_id);
    }

    void Set_Bl_energy(vector<double> GB_BLE_vector) {
        Bl_energy = GB_BLE_vector.at(GB_id);
    }

    void Set_Cl_energy(vector<double> GB_CLE_vector){
        Cl_energy = GB_CLE_vector.at(GB_id);
    }

/// Get values methods
/*    double Get_surface_energy(unsigned int GB_id) const {
//        if(surface_energy != 0)
//            return surface_energy;
//        else {
//            Set_surface_energy(GB_SE_vector);
//            return surface_energy; } }
//    double Get_external_elastic_energy(unsigned int GB_id) const { /// EEE: external elastic energy vector
//        if(external_elastic_energy != 0)
//            return external_elastic_energy;
//        else {
//            Set_external_elastic_energy(GB_EEE_vector);
//            return external_elastic_energy; } }
//    double Get_crack_interaction_energy(unsigned int GB_id) const {
//        if(crack_interaction_energy != 0)
//            return crack_interaction_energy;
//        else {
//            Set_crack_interaction_energy(GB_CIE_vector);
//            return crack_interaction_energy; } }
//    double Get_Bl_energy(int GB_id) const {
//        if(Bl_energy != 0)
//            return Bl_energy;
//        else { Set_Bl_energy(GB_BLE_vector);
//            return Bl_energy; } }
//    double Get_Cl_energy(unsigned int GB_id) const {
//        if(Cl_energy != 0)
//            return Cl_energy;
//        else { Set_Cl_energy(GB_CLE_vector);
//            return Cl_energy; } }
//    double Get_total_energy() const {
        /// total_energy; //= surface_energy + external_elastic_energy + Bl_energy + Cl_energy;
//    }

private:

/// Combinatorial
    int GB_edges_number;  //equal to the number of neighbours
    int GB_nodes_number;

///Geometry
    double GB_area;
    double GB_perimeter;

    tuple<double,double,double> GB_barycentre_coordinates;

///Energies
    double surface_energy;
    double external_elastic_energy;
    double crack_interaction_energy;
    double Bl_energy;
    double Cl_energy;
    double total_energy; //= surface_energy + external_elastic_energy + Bl_energy + Cl_energy;
};

/// # 4 # The class of stress concentrators related to defects in one special face element of a PCC

class face_concentrator {
    double total_elastic_energy = 0;

private:
    //string conc_type; // bl or cl
    unsigned int conc_face_number = 0;
    unsigned int bl_index = 0;
    unsigned int cl_index = 0;
    double Bl_elastic_energy = 0;
    double Cl_elastic_energy = 0;
public:

    face_concentrator (unsigned int ConcFace) { // constructor 1 simple
        conc_face_number = ConcFace;
    }

    face_concentrator (unsigned int ConcFace, unsigned int Bl, unsigned int Cl) { // constructor 2 complex
        conc_face_number = ConcFace;
        bl_index = Bl;
        cl_index = Cl;
    }

    unsigned int GetBl_index()
    {
        return bl_index;
    }

    int GetCl_index()
    {
        return cl_index;
    }
}; // end of class

 */