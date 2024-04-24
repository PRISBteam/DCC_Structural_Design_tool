#ifndef PCC_PROCESSING_DESIGN_PCC_OBJECTS_H
#define PCC_PROCESSING_DESIGN_PCC_OBJECTS_H

#include <Eigen/SparseCore>

/// ==== # 0 # =============== Structure for the initial configuration  ========================= ///

/*!
 * @brief This class combine PCCpaths to directories and initial variables set in the config/_.ini files with the methods of their reading like
 * @public Get_config(), Set_config()
 * @param dim, source_dir, paths, ConfigVector
 * @return Configuration_sState, Configuration_cState
 */
class Config {
    private:
        int config_dim;
        std::string config_source_dir, config_output_dir; // Input and output directories as it is written in the 'config/main.ini' file
        std::string pcc_standard; // PCC standard as specified in the technical documentation for the project
        //std::vector<char*> config_PCCpaths;
        std::vector<std::string> config_PCCpaths;
        std::vector<int> config_ConfVector; // main module keys
        std::string config_main_type; // 'mode' from the config/main.ini file: 'LIST' (execution one by one all the active (ON) project modules), 'TUTORIAL' as a specific education mode, 'PERFORMANCE_TEST' or the 'TASK' mode :: This define the global simulation mode: 'LIST' for the "list" of modules implementing one by one (if ON) and 'TASK' for the user-defined task scripts with modules and functions included from the project's libraries
        std::string config_sim_task; // path to the corresponding *.cpp file containing a 'simulation task' (for 'TASK' execution mode only, not 'LIST') as it is written in the 'config/main.ini' file

    /// The list of all mentioned below State_<*>_vectors and State_<*>fracture_vectors as the output of the Processing module // is the list of 'state vectors' analogous to the Configuration_sState but for 'cracked' (or induced) network of k-cells
    /* where 'n' :: "nodes", 'e' :: "edges", 'f' :: "faces", and 'p' :: "polyhedrons" */
// State_Vector in the form : [Element index] - > [Element type], like [0, 0, 2, 1, 1, 0, 2, 4, 3, 3, 2, 0,... ...,2] containing all CellNumb.at(*) element types
        std::vector<unsigned int> State_p_vector, State_f_vector, State_e_vector, State_n_vector; // Normally the State_<*>_vector of special cells can be calculated based on the corresonding special_cell_sequences
        std::vector<unsigned int> State_pfracture_vector, State_ffracture_vector, State_efracture_vector, State_nfracture_vector; // separate vectors containing the other 'fractured' labels different from the 'special' ones. To be calculated based on the corresonding fractured_cell_sequences

/// Configuration_sState = { State_p_vector, State_f_vector, State_e_vector, State_n_vector } is a list of all 'state vectors': from (1) State_p_vector (on top, id = 0) to (4) State_n_vector (bottom, id = 3)
        std::vector<std::vector<unsigned int>> Configuration_sState;
        std::vector<std::vector<unsigned int>> Configuration_cState;

public:
    void Read_config(); // Read the 'initial configuration' of the problem set in all the relevant '*.ini' files containing in the '\config' project directory using the functions from the 'ini_readers.cpp' project library (and only from there)
    void Set_config(const std::vector<int> &ConfigVector, const std::string &source_dir, int &dim, std::vector<char*> paths, std::vector<std::vector<int>> Configuration_State, std::vector<std::vector<int>> Configuration_cState); // manual setting of the configuration

    int Get_dim(); //@return dim
    std::vector<int> Get_ConfVector(); //!@return ConfVector
    std::string Get_source_dir(); //!@return source_dir
    std::string Get_output_dir(); //!@return output_dir
    std::string Get_pcc_standard(); //!@return pcc_standard
    std::vector<std::string> Get_paths(); //!@return PCC PCCpaths
    std::string Get_main_type(); //!@return main_type
    std::string Get_sim_task(); //!@return sim_task path to the corresponding *.cpp file containing the task code

    std::vector<std::vector<unsigned int>> Get_Configuration_sState(); //!@return Configuration_sState
    std::vector<std::vector<unsigned int>> Get_Configuration_iState(); //!@return Configuration_iState
};
// ConfigVector (../config/main.ini) contains ALL the control variables needed for the program execution

/// ==== # 1 # =============== CellDesign class  ========================= ///
class CellDesign {
private:
    std::vector<unsigned int> p_special_design, f_special_design, e_special_design, n_special_design; // state vectors of special k-cells
    std::vector<unsigned int> p_induced_design, f_induced_design, e_induced_design, n_induced_design; // state vectors of induced k-cells

    std::vector<unsigned int> p_special_sequence, f_special_sequence, e_special_sequence, n_special_sequence; // sequences of UNIQUE special k-cell numbers
    std::vector<unsigned int> p_induced_sequence, f_induced_sequence, e_induced_sequence, n_induced_sequence; // sequences of UNIQUE induced k-cell numbers

    std::vector<std::vector<unsigned int>> p_special_series, f_special_series, e_special_series, n_special_series; // sequences of sequences of special k-cell numbers
    std::vector<std::vector<unsigned int>> p_induced_series, f_induced_series, e_induced_series, n_induced_series; // sequences of sequences of induced k-cell numbers


public:
    /// Set of variables
    CellDesign() {}; // constructor
    void Set_sequences(std::vector<unsigned int> psequence, std::vector<unsigned int> fsequence, std::vector<unsigned int> esequence, std::vector<unsigned int> nsequence);
    void Set_induced_sequences(std::vector<unsigned int> p_ind_sequence, std::vector<unsigned int> f_ind_sequence, std::vector<unsigned int> e_ind_sequence, std::vector<unsigned int> n_ind_sequence);
    void Set_designes(std::vector<unsigned int> pdesign, std::vector<unsigned int> fdesign, std::vector<unsigned int> edesign, std::vector<unsigned int> ndesign);
    void Set_induced_designs(std::vector<unsigned int> p_ind_design, std::vector<unsigned int> f_ind_design, std::vector<unsigned int> e_ind_design, std::vector<unsigned int> n_ind_design);
    void Set_sequence(std::vector<unsigned int> sequence, int ctype);
    void Set_induced_sequence(std::vector<unsigned int> ind_sequence, int ctype);
    void Set_design(std::vector<unsigned int> design, int ctype);
    void Set_induced_design(std::vector<unsigned int> ind_design, int ctype);

    // Get
    std::vector<unsigned int> Get_p_sequence(void);
    std::vector<unsigned int> Get_f_sequence(void);
    std::vector<unsigned int> Get_e_sequence(void);
    std::vector<unsigned int> Get_n_sequence(void);
    std::vector<unsigned int> Get_p_induced_sequence(void);
    std::vector<unsigned int> Get_f_induced_sequence(void);
    std::vector<unsigned int> Get_e_induced_sequence(void);
    std::vector<unsigned int> Get_n_induced_sequence(void);

    std::vector<unsigned int> Get_p_design(void);
    std::vector<unsigned int> Get_f_design(void);
    std::vector<unsigned int> Get_e_design(void);
    std::vector<unsigned int> Get_n_design(void);
};

/// ==== # 2 # =============== Processed Complex class  ========================= ///

class ProcessedComplex { // Essential for Characterisation module
// PCC processed with all its characteristics and design sequences

private:

public:
    /// Set variables
    CellDesign pcc_design;

    void Set_design(CellDesign processed_pcc_design);

    // Sequences of special k-cells
    std::vector<std::vector<unsigned int>> face_process_seq;
    std::vector<std::vector<int>> face_process_state;

    // Entropic analysis
    std::vector<double> e_entropy_mean_vector, e_entropy_skrew_vector, e_entropy_full_vector;
    std::vector<std::vector<double>> je_fractions_vector, de_fractions_vector;

    // Analytical solutions
    std::vector<std::vector<double>> j_analytical_rand_vector, d_analytical_rand_vector;
    std::vector<std::vector<double>> j_analytical_cryst_vector, d_analytical_cryst_vector;
    std::vector<std::tuple<double, double>> AnRandEntropies_vector, AnCrystEntropies_vector;

    // Laplacian lab
    std::vector<std::vector<double>> Betti_vector;
}; /// end of class ProcessedComplex

/// ==== # 3 # =============== Subcomplex class  ========================= ///

class Subcomplex {

private:
    /// 1. Combinatorics
    std::vector <unsigned int> sub_grains_sequence;
    std::vector <unsigned int> sub_faces_sequence;
    std::vector <unsigned int> sub_nodes_sequence;
    std::vector <unsigned int> common_faces_sequence;
    std::vector <unsigned int> s_sub_faces_sequence;
    std::vector <unsigned int> c_sub_faces_sequence;

    /// 2. Geometry
    //vector<tuple<double, double, double>> vertex_coordinates;
    std::vector<std::tuple<double, double, double>> common_faces_coordinates;
    std::vector<std::tuple<double, double, double>> sub_grain_coordinates;

public:
    unsigned int subcomplex_id;
    double sub_length;
    double a_n; double b_n; double c_n; double D_plane;
    std::vector<double> crack_plane = {a_n, b_n, c_n, D_plane};

    Subcomplex(unsigned int subcomplex_id_new);
    Subcomplex(unsigned int subcomplex_id_new, std::vector <unsigned int> new_sub_grains_sequence);

    /// Grains
    void Set_grains_sequence(std::vector <unsigned int> new_sub_grains_sequence);
    std::vector <unsigned int> Get_grains_sequence(unsigned int subcomplex_id);
    /// geometry
    void Set_sub_grain_coordinates(std::vector<std::tuple<double, double, double>> new_sub_grain_coordinates);
    std::vector<std::tuple<double, double, double>> Get_sub_grain_coordinates(unsigned int subcomplex_id);
    /// Faces
    void Set_faces_sequence(std::vector <unsigned int> new_sub_faces_sequence);
    std::vector <unsigned int>  Get_faces_sequence(unsigned int  subcomplex_id);
    void Set_common_faces_sequence(std::vector <unsigned int> new_common_faces_sequence);

    std::vector <unsigned int> Get_common_faces_sequence(unsigned int subcomplex_id);

    void Set_sfaces_sequence(std::vector <unsigned int> const &ssub_faces_sequence);
    std::vector <unsigned int> Get_sfaces_sequence(unsigned int  subcomplex_id);
    void Set_cfaces_sequence(std::vector <unsigned int> csub_faces_sequence);
    std::vector <unsigned int> Get_cfaces_sequence(unsigned int  subcomplex_id);

    /// geometry
    void Set_common_faces_coordinates(std::vector<std::tuple<double, double, double>> new_common_faces_coordinates);
    std::vector<std::tuple<double, double, double>> Get_common_faces_coordinates(unsigned int subcomplex_id);
};

/// ==== # 4 # =============== grain3D class  ========================= ///

class Polytope {

    std::vector<std::tuple<double, double, double>> minmax_node_coordinates; // a vecor containing two tuples: gmincoord{xmin,ymin,zmin},gmaxcoord{xmax,ymax,zmax}

private:
    // list of nodes
    std::vector<unsigned int> node_ids;

    // list of Faces
    std::vector<unsigned int> Faces_list;

    // list of triplets of nodes coordinated
    std::vector<std::tuple<double, double, double>> node_coordinates;

public:
    unsigned int grain_id;

    Polytope(unsigned int grain_new_id); // constructor 1

    void Set_node_ids(unsigned int grain_id, Eigen::SparseMatrix<double> const &GFS, Eigen::SparseMatrix<double> const &FES, Eigen::SparseMatrix<double> const &ENS);

    void Set_Faces_list(unsigned int grain_id, Eigen::SparseMatrix<double> const &GFS);

    std::vector<unsigned int> Get_Faces_list();

    /// return - vector of all node (vertices) coordinates of a grain
    void Set_node_coordinates(unsigned int grain_id);

    std::vector<unsigned int> Get_node_ids(unsigned int grain_id);

    std::vector<std::tuple<double, double, double>> Get_node_coordinates(unsigned int grain_id);

    /// return - vector with two tuples : { x_min, y_min, z_min; x_max, y_max, z_max} of a grain witn number grain_id
    std::vector<std::tuple<double, double, double>> Get_minmax_node_coordinates(unsigned int grain_id);

}; // end of class Polytope

/// ========== END of class Polytope functions description

/// # 6 # The class of a MACROCRACK
class Macrocrack {
    double total_fracture_energy = 0;
    Subcomplex half_plane_subcomplex; // geometry part

private:
    double a_n;
    double b_n;
    double c_n;
    double D_plane;
    double crack_length;
    double real_crack_length;
    std::vector<double> crack_plane = {a_n, b_n, c_n, D_plane};

public:
    int crack_id = 0;
    double surface_energy = 0;
    double bridging_energy = 0;
    double multiple_cracking_energy = 0;
    double stress_concentrators_energy = 0;

    Macrocrack(int crack_id_new, Subcomplex &half_plane_sub);

    double Get_crack_length(int crack_id_new);

    void Set_real_crack_length(double sample_size);

    double Get_real_crack_length();

    void Set_crack_plane();

    void Set_multiple_cracking_energy(double total_energy);

    double Get_multiple_cracking_energy();

    std::vector <unsigned int> Get_faces_sequence();

    std::vector <unsigned int> Get_sfaces_sequence();

    std::vector<double> Get_crack_plane();

    std::vector <std::tuple<double,double,double>> Get_common_faces_coordinates(unsigned int  crack_id);

}; // end of class MACROCRACK

/// ========== END of class Macrocrack functions description

/// # V # The class of a CELLS_ENERGIES :: list of the energy_vectors corresponding to different dimensions 'k' of the k-cells in a PCC
/// class CellEnergies {
///
/// }; // end of class CELLS_ENERGIES

/// ========== END of class CellEnergies functions description


#endif //PCC_PROCESSING_DESIGN_PCC_OBJECTS_H
