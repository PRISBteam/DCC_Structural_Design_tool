/// Classes of special objects related to defect structures on a PCC elements
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace std; // standard namespace
using namespace Eigen; // standard namespace


extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector; // coordinate vectors defined globally
extern std::vector<unsigned int> CellNumbs; // number of cells in a PCC defined globally

typedef Eigen::SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'

#include "PCC_Objects.h"

/// # 1 # The class of a CELLS_DESIGN :: list of the state_vectors corresponding to different dimensions 'k' of the k-cells in a PCC
//class CellsDesign
    void CellsDesign::Set_sequences(std::vector<unsigned int> psequence, std::vector<unsigned int> fsequence, std::vector<unsigned int> esequence, std::vector<unsigned int> nsequence){
        p_sequence = psequence; f_sequence = fsequence; e_sequence = esequence; n_sequence = nsequence;
    }
    void CellsDesign::Set_induced_sequences(std::vector<unsigned int> p_ind_sequence, std::vector<unsigned int> f_ind_sequence, std::vector<unsigned int> e_ind_sequence, std::vector<unsigned int> n_ind_sequence){
        p_induced_sequence = p_ind_sequence; f_induced_sequence = f_ind_sequence; e_induced_sequence = e_ind_sequence; n_induced_sequence = n_ind_sequence;
    }
    void CellsDesign::Set_designes(std::vector<int> pdesign, std::vector<int> fdesign, std::vector<int> edesign, std::vector<int> ndesign){
        p_design = pdesign; f_design = fdesign; e_design = edesign; n_design = ndesign;
    }
    void CellsDesign::Set_induced_designs(std::vector<int> p_ind_design, std::vector<int> f_ind_design, std::vector<int> e_ind_design, std::vector<int> n_ind_design){
        p_induced_design = p_ind_design; f_induced_design = f_ind_design; e_induced_design = e_ind_design; n_induced_design = n_ind_design;
    }
    void CellsDesign::Set_sequence(std::vector<unsigned int> sequence, int ctype){
        if (ctype == 3) p_sequence = sequence;
        else if (ctype == 2 ) f_sequence = sequence;
        else if (ctype == 1 ) e_sequence = sequence;
        else if (ctype == 0 ) n_sequence = sequence;
    }
    void CellsDesign::Set_induced_sequence(std::vector<unsigned int> ind_sequence, int ctype){
        if (ctype == 3) p_induced_sequence = ind_sequence;
        else if (ctype == 2 ) f_induced_sequence = ind_sequence;
        else if (ctype == 1 ) e_induced_sequence = ind_sequence;
        else if (ctype == 0 ) n_induced_sequence = ind_sequence;
    }
    void CellsDesign::Set_design(std::vector<int> design, int ctype){
        if (ctype == 3) p_design = design;
        else if (ctype == 2 ) f_design = design;
        else if (ctype == 1 ) e_design = design;
        else if (ctype == 0 ) n_design = design;
    }
    void CellsDesign::Set_induced_design(std::vector<int> ind_design, int ctype){
        if (ctype == 3) p_induced_design = ind_design;
        else if (ctype == 2 ) f_induced_design = ind_design;
        else if (ctype == 1 ) e_induced_design = ind_design;
        else if (ctype == 0 ) n_induced_design = ind_design;
    }
    // Get
    std::vector<unsigned int> CellsDesign::Get_p_sequence(void){
        if (p_sequence.size() == 0) {
            cout << "WARNING: p_sequence did not set!" << endl;
            return {0};
        }
        else return p_sequence;
    }
    std::vector<unsigned int> CellsDesign::Get_f_sequence(void){
        if (f_sequence.size() == 0) {
            cout << "WARNING: f_sequence did not set!" << endl;
            return {0};
        }
        else return f_sequence;
    }
    std::vector<unsigned int> CellsDesign::Get_e_sequence(void){
        if (e_sequence.size() == 0) {
            cout << "WARNING: e_sequence did not set!" << endl;
            return {0};
        }
        else return e_sequence;
    }
    std::vector<unsigned int> CellsDesign::Get_n_sequence(void){
        if (n_sequence.size() == 0) {
            cout << "WARNING: n_sequence did not set!" << endl;
            return {0};
        }
        else return n_sequence;
    }
    std::vector<unsigned int> CellsDesign::Get_p_induced_sequence(void) {
        if (p_induced_sequence.size() == 0) {
            cout << "WARNING: polyhedron induced sequence did not set!" << endl;
            return {0};
        }
        else return p_induced_sequence;
    }
    std::vector<unsigned int> CellsDesign::Get_f_induced_sequence(void){
        if (f_induced_sequence.size() == 0) {
            cout << "WARNING: face induced sequence did not set!" << endl;
            return {0};
        }
        else return f_induced_sequence;
    }
    std::vector<unsigned int> CellsDesign::Get_e_induced_sequence(void) {
        if (e_induced_sequence.size() == 0) {
            cout << "WARNING: edge induced sequence did not set!" << endl;
            return {0};
        }
        else return e_induced_sequence;
    }
    std::vector<unsigned int> CellsDesign::Get_n_induced_sequence(void) {
        if (n_induced_sequence.size() == 0) {
            cout << "WARNING: node induced sequence did not set!" << endl;
            return {0};
        }
        else return n_induced_sequence;
    }
    std::vector<int> CellsDesign::Get_p_design(void){
        if (p_design.size() == 0) {
            cout << "WARNING: p_design did not set!" << endl;
            return {0};
        }
        else return p_design;
    }
    std::vector<int> CellsDesign::Get_f_design(void){
        if (f_design.size() == 0) {
            cout << "WARNING: f_design did not set!" << endl;
            return {0};
        } else return f_design;
    }
    std::vector<int> CellsDesign::Get_e_design(void){
        if (e_design.size() == 0) {
            cout << "WARNING: e_design did not set!" << endl;
            return {0};
        } else return e_design;
    }
    std::vector<int> CellsDesign::Get_n_design(void){
        if (n_design.size() == 0) {
            cout << "WARNING: n_design did not set!" << endl;
            return {0};
        }
        else return n_design;
    }
/// ========== END of class CELLS_DESIGN functions description

/// # 3 # The class of SUBCOMPLEX

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
        sub_grains_sequence = new_sub_grains_sequence; }

    std::vector <unsigned int> Subcomplex::Get_grains_sequence(unsigned int subcomplex_id){
        if(sub_grains_sequence.size() != 0)
            return sub_grains_sequence;
        else return {0};
    }

    /// geometry
    void Subcomplex::Set_sub_grain_coordinates(std::vector<tuple<double, double, double>> new_sub_grain_coordinates){
        sub_grain_coordinates = new_sub_grain_coordinates; }
    std::vector<tuple<double, double, double>> Subcomplex::Get_sub_grain_coordinates(unsigned int subcomplex_id){
        return sub_grain_coordinates; }

    //2
//    std::vector <unsigned int>  Get_grains_sequence(unsigned int subcomplex_id, std::vector <unsigned int> new_sub_grains_sequence){
//        if(sub_grains_sequence.size() != 0) return sub_grains_sequence;
//        else { Set_grains_sequence(new_sub_grains_sequence); return sub_grains_sequence; } }

    /// Faces
    void Subcomplex::Set_faces_sequence(std::vector <unsigned int> new_sub_faces_sequence){
        sub_faces_sequence = new_sub_faces_sequence; }
    //1
    std::vector <unsigned int>  Subcomplex::Get_faces_sequence(unsigned int  subcomplex_id){
        if(sub_faces_sequence.size() != 0)
            return sub_faces_sequence;
        else return {0};
    }

    void Subcomplex::Set_common_faces_sequence(std::vector <unsigned int> new_common_faces_sequence){
        common_faces_sequence = new_common_faces_sequence; }
    std::vector <unsigned int> Subcomplex::Get_common_faces_sequence(unsigned int subcomplex_id){
        return common_faces_sequence; }

    void Subcomplex::Set_sfaces_sequence(std::vector <unsigned int> const &ssub_faces_sequence){
        s_sub_faces_sequence = ssub_faces_sequence;
    }
    std::vector <unsigned int> Subcomplex::Get_sfaces_sequence(unsigned int  subcomplex_id){
        if(s_sub_faces_sequence.size() > 0) return s_sub_faces_sequence;
        else {
            cout << "Caution! s_sub_faces_sequence = 0! Please add it to the corresponding subcomplex/crack"s << endl;
            return {0};
        }
    }

    void Subcomplex::Set_cfaces_sequence(std::vector <unsigned int> csub_faces_sequence){
        c_sub_faces_sequence = csub_faces_sequence;
    }
    std::vector <unsigned int> Subcomplex::Get_cfaces_sequence(unsigned int  subcomplex_id){
        return c_sub_faces_sequence;
    }

    /// geometry
    void Subcomplex::Set_common_faces_coordinates(std::vector<tuple<double, double, double>> new_common_faces_coordinates){
        common_faces_coordinates = new_common_faces_coordinates; }

    std::vector<tuple<double, double, double>> Subcomplex::Get_common_faces_coordinates(unsigned int subcomplex_id){
        return common_faces_coordinates; }

/// ========== END of class SUBCOMPLEX

/// # 4 # The class of a PROCESSED COMPLEX
    void ProcessedComplex::Set_design(CellsDesign processed_pcc_design) {
        pcc_design = processed_pcc_design;
    }

/// ========== END of class PROCESSED COMPLEX functions description


/// # 5 # The class of Polytopes in a PCC

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

    vector<unsigned int> Polytope::Get_Faces_list() {
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

    std::vector<unsigned int> Polytope::Get_node_ids(unsigned int grain_id) { // set the node ids
        return node_ids;
    }

std::vector<std::tuple<double, double, double>> Polytope::Get_node_coordinates(unsigned int grain_id) { // set the node ids
        if (node_coordinates.size() != 0) {
            return node_coordinates;
        } else if (node_coordinates.size() > 0) {
            Set_node_coordinates(grain_id);
            return node_coordinates;
        } else {
            throw std::invalid_argument(
                    "Please call Set_node_coordinates(unsigned int grain_id, vector<tuple<double, double, double>> const &vertex_coordinates) method first!");
            return node_coordinates;
        }
    } // end of  Get_node_coordinates() method

    /// return - vector with two tuples : { x_min, y_min, z_min; x_max, y_max, z_max} of a grain witn number grain_id
    vector<tuple<double, double, double>> Polytope::Get_minmax_node_coordinates(unsigned int grain_id) { // min and max {x,y,z} values of vertices for a grain
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

    double Macrocrack::Get_crack_length(int crack_id_new) {
        return crack_length;
    }

    void Macrocrack::Set_real_crack_length(double sample_size) {
        real_crack_length = half_plane_subcomplex.sub_length * sample_size;
    }
    double Macrocrack::Get_real_crack_length() {
        if (real_crack_length > 0.0) return real_crack_length;
        else {
            cout << "Caution! real_crack_length = 0! Please use Set_real_crack_length(double sample_size) before!"s << endl;
            return 0;
        }
    }

    void Macrocrack::Set_crack_plane() {
        crack_plane = {half_plane_subcomplex.a_n, half_plane_subcomplex.b_n, half_plane_subcomplex.c_n, half_plane_subcomplex.D_plane};
    }
    vector<double> Macrocrack::Get_crack_plane() {
        if (crack_plane.size() > 0.0) return crack_plane;
        else {
            cout << "Caution! crack_plane.size() = 0! Please update first {half_plane_subcomplex.a_n, half_plane_subcomplex.b_n, half_plane_subcomplex.c_n, half_plane_subcomplex.D_plane} in the corresponding subcomplex!"s << endl;
            return {0};
        }
    }

    void Macrocrack::Set_multiple_cracking_energy(double total_energy) {
        multiple_cracking_energy = total_energy;
    }
    double Macrocrack::Get_multiple_cracking_energy() {
        return multiple_cracking_energy;
    }

    std::vector <unsigned int> Macrocrack::Get_faces_sequence(){
        return half_plane_subcomplex.Get_faces_sequence(crack_id); }

    std::vector <unsigned int> Macrocrack::Get_sfaces_sequence(){
        return half_plane_subcomplex.Get_sfaces_sequence(0); }

    std::vector <tuple<double,double,double>> Macrocrack::Get_common_faces_coordinates(unsigned int  crack_id){
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
/*    double Get_surface_energy(unsigned int GB_id){
//        if(surface_energy != 0)
//            return surface_energy;
//        else {
//            Set_surface_energy(GB_SE_vector);
//            return surface_energy; } }
//    double Get_external_elastic_energy(unsigned int GB_id) { /// EEE: external elastic energy vector
//        if(external_elastic_energy != 0)
//            return external_elastic_energy;
//        else {
//            Set_external_elastic_energy(GB_EEE_vector);
//            return external_elastic_energy; } }
//    double Get_crack_interaction_energy(unsigned int GB_id) {
//        if(crack_interaction_energy != 0)
//            return crack_interaction_energy;
//        else {
//            Set_crack_interaction_energy(GB_CIE_vector);
//            return crack_interaction_energy; } }
//    double Get_Bl_energy(int GB_id){
//        if(Bl_energy != 0)
//            return Bl_energy;
//        else { Set_Bl_energy(GB_BLE_vector);
//            return Bl_energy; } }
//    double Get_Cl_energy(unsigned int GB_id){
//        if(Cl_energy != 0)
//            return Cl_energy;
//        else { Set_Cl_energy(GB_CLE_vector);
//            return Cl_energy; } }
//    double Get_total_energy(){
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

/// # 3 # The class of Agglomeration of defects in one special face element of a PCC
class agglomeration {
    double adhesion_energy = 0;
    double surface_energy = 0;

private:
    string atype; // like "rgo"
    unsigned int aface_number = 0;
    int          apower = 0;
    int          a_average_strip_length = 0;

public:
    agglomeration(unsigned int AFace) { // constructor 1 simple
        aface_number = AFace;
    }

    agglomeration(unsigned int AFace, int AglPower) { // constructor 2 complex
        aface_number = AFace;
        apower = AglPower;
    }

    void SetAFace(unsigned int AFace)
    {
        aface_number = AFace;
    }

    void SetAType(std::string type)
    {
        atype = type;
        if (atype == "rgo") {
            surface_energy = 0.2; // units [J/m^2]
            adhesion_energy = 0.4; // units [J/m^2]
        }

    }

    void SetAPower(vector<vector<int>> const &RW_series_vector) /// Agglomerations power
    {
        int AglPower = 0;
        //vector<int> aggl_vector(CellNumbs.at(2), 0); // state vector for agglomerations (with the size # cells initially filled with all 0s) : contains # of faces with agglomerations
        for (auto RWsfv : RW_series_vector)
            for (auto RWsf: RWsfv)
                if(RWsf == aface_number) AglPower += 1;

        apower = AglPower;
    }

    void SetAvLength(vector<vector<int>> const &RW_series_vector) /// Average length of strips related to this agglomeration
    {
        int ATotalLength = 0;
        for (auto RWsfv : RW_series_vector)
            for (auto RWsf: RWsfv)
                if(RWsf == aface_number) ATotalLength += RWsfv.size();

        a_average_strip_length = ATotalLength/ (double) GetAPower(RW_series_vector);
    }

    unsigned int GetAFace()
    {
        return aface_number;
    }

    int GetAPower()
    {
        return apower;
    }

    int GetAPower(vector<vector<int>> const &RW_series_vector)
    {
        if (apower != 0) {
            return apower;
        } else {
            int AglPower = 0;
            for (auto RWsfv : RW_series_vector)
                for (auto RWsf: RWsfv)
                    if(RWsf == aface_number) AglPower += 1;

            apower = AglPower;
            return apower;
        }
    }

    int GetAvLength()
    {
        return a_average_strip_length;
    }

    int GetAvLength(vector<vector<int>> const &RW_series_vector) {
        if (a_average_strip_length != 0) {
            return a_average_strip_length;
        } else {
            int ATotalLength = 0;
            for (auto RWsfv: RW_series_vector)
                for (auto RWsf: RWsfv)
                    if (RWsf == aface_number) ATotalLength += RWsfv.size();

            a_average_strip_length = ATotalLength / (double) GetAPower(RW_series_vector);
            return a_average_strip_length;
        }
    }

}; // end of class agglomeration

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