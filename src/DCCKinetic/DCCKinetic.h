///================================ DCC Kinetic module =============================================================///
///=======================================================================================================================///
/** The function in this library generate different quasi-random finetic processeson the elements of the pre-constructed ///
*   discrete sell complex (DCC)                                                                                         **/
///=====================================================================================================================///
/// 'W' for the 3D one-layer film, 'P' for the Ising-like model of Plasticity, 'F' for the Ising-like model of Fracture

// Triplets in the form T = T(i,j,value), where i and j element's indices in the corresponding dense matrix
typedef Triplet<double> Tr; // <Eigen library class> Declares a triplet's type name - Tr
typedef SparseMatrix<double> SpMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles name - SpMat
typedef tuple<double, double, double> Tup; // Eigen library class
typedef pair<double, double> Pr; // Eigen library class

/// Standard (STL) C++ libraries:
///------------------------------
///------------------------------
/// Attached user defined C++ libraries:
///-------------------------------------
///-------------------------------------
#include "Kinetic_Functions.h"
///-------------------------------------

std::vector <unsigned int> DCC_Kinetic(char* stype, std::vector<unsigned int> &s_faces_sequence, std::vector<char*> paths, char* input_folder, char* odir) {
    /// Specific functions
    std::vector <unsigned int> face_sequence;

/// Declaration of FUNCTIONS, see the function bodies at the end of file.
    Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols); // The function read any matrices from lists of triplets and create corresponding sparse matrices

//    std::vector<unsigned int> VectorReader(char* FilePath); // The function read any matrices from lists of triplets and create corresponding sparse matrices
//    std::vector<int> confCout(char* config, vector<int> const& configuration); // Read and Output configuration

/// Reading vector from the file "number_of_cells" the numbers of cells od different types ///
// ::      vector components: [0] - Nodes, [1] - Edges, [2] - Faces, [3] - Grains ::     ///
    // File path with the amounts of the different cells  (1st line for Nodes, 2nd line for Edges, 3rd line for Faces and (in 3D) 4th line for grains)
    string  ncells = input_folder + "number_of_cells.txt"s; char* number_of_cells = const_cast<char*>(ncells.c_str());
    std::vector<unsigned int> CellNumbs = VectorReader(number_of_cells);
    //Screen output for the numbers of cells in the DCC
    cout << "The number of different k-cells in the DCC:" << endl;
    unsigned int t_length = 0;
    for (int j : CellNumbs) cout << t_length++ << "-cells #\t" << j << endl;

////////////////////// Matrices initialisation part //////////////////////
    std::vector<Tr> SFaces_Triplet_list; // Probe vector of triplets

////////==================== Reading from files of all the adjacency matrices of the DCC ===================////
// AN - Nodes (Vertices or 0-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AE - Edges (Triple Junctions or 1-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AF - Faces (Grain Boundaries or 2-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AG - Grains (Volumes or 3-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MEN - Edges - Nodes (1-Cells to 0-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MFE - Faces - Edges (2-Cells to 1-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MGF - Grains - Faces (3-Cells to 2-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values // odir - source path
    /// Adjacency sparse matrix for nodes /// Adjacency sparse matrix for edges /// Adjacency sparse matrix for faces /// Adjacency sparse matrix for grains
//    SpMat ANS(CellNumbs.at(0), CellNumbs.at(0)), AES(CellNumbs.at(1), CellNumbs.at(1)),
//            AFS(CellNumbs.at(2), CellNumbs.at(2)), AGS(CellNumbs.at(3), CellNumbs.at(3));
SpMat AFS(CellNumbs.at(2), CellNumbs.at(2));
//    ANS = SMatrixReader(paths.at(0), (CellNumbs.at(0)), (CellNumbs.at(0))); //all Nodes
//    ANS = 0.5 * (ANS + SparseMatrix<double>(ANS.transpose())); // Full matrix instead of triagonal
//    AES = SMatrixReader(paths.at(1), (CellNumbs.at(1)), (CellNumbs.at(1))); //all Edges
//    AES = 0.5 * (AES + SparseMatrix<double>(AES.transpose())); // Full matrix instead of triagonal
    AFS = SMatrixReader(paths.at(2), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose())); // Full matrix instead of triagonal
//   AGS = SMatrixReader(paths.at(3), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
//    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose())); // Full matrix instead of triagonal
/// Incidence sparse matrix for Edges and Nodes /// Incidence sparse matrix for Faces and Edges /// Incidence sparse matrix for Grains and Faces
//    SpMat ENS(CellNumbs.at(0), CellNumbs.at(1)), FES(CellNumbs.at(1), CellNumbs.at(2)),
//            GFS(CellNumbs.at(2),CellNumbs.at(3));
//    ENS = SMatrixReader(paths.at(4), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges
SpMat   FES = SMatrixReader(paths.at(5), (CellNumbs.at(1)), (CellNumbs.at(2))); //all Edges-Faces
//    GFS = SMatrixReader(paths.at(6), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains

////////////////////////////////////// KINETIC PROCESSES ///////////////////////////////////////////////
    if (*stype == 'W') { // Wear
        double  ShearStress = 3.0*pow(10,8);
        /// Creation/Reading of grain orientations
        double Ori_angle = 0;
        vector<Tup> Grain_Orientations; // Vectors of triplets (in 2D {x,y,0}) for grain orientations
        for (unsigned int k = 0; k < CellNumbs.at(3); k++) {
            Ori_angle = rand() % 180; // Random generation of angle
            //if (Ori_angle < 90)
            Grain_Orientations.push_back(make_tuple(cos(Ori_angle),sin(Ori_angle),0));
           // cout << get<0>(Grain_Orientations.at(k)) << "\t" << get<1>(Grain_Orientations.at(k)) << "\t" <<get<2>(Grain_Orientations.at(k)) << "\t" << endl;
        }
        /// ======= Kinetic function ===================>>
        DCC_Kinetic_Wear(ShearStress, Grain_Orientations, FES, CellNumbs, input_folder, odir);

        /// ===== Grain orientations output =========>
        ofstream GrainOrientationsStream;

        GrainOrientationsStream.open(odir + "GrainOrientations.txt"s, ios::trunc);
        if (GrainOrientationsStream) {
            GrainOrientationsStream << "Grain Orientations" << endl;
            for (auto go : Grain_Orientations)
            GrainOrientationsStream << get<0>(go) << "\t" << get<1>(go) << endl; //<< "\t" <<get<2>(go)
            GrainOrientationsStream.close();
        } else cout << "Error: No such a directory for\t" << odir + "GrainOrientations.txt"s << endl;

    } ///End of 'Wear' type simulations

    else if (*stype == 'P') { // Plasticity
        vector<Tup> fraction_stress_temperature;

        /// ========= DCC_Kinetic main function call ==========>>

        fraction_stress_temperature = DCC_Kinetic_Plasticity(FES, CellNumbs, input_folder, odir);

        /// Analysis and output to file
        ofstream FSTStream;
        FSTStream.open(odir + "fraction_stress_temperature.txt"s, ios::trunc);
        if (FSTStream) {
            FSTStream << "(1) Plastic strain [%] (2) Yield strength [MPa] (3) Temperature [K])" << endl;
            for (auto go : fraction_stress_temperature)
                FSTStream << get<0>(go) * 100.0 << " \t" << get<1>(go)/pow(10,6) << " \t" << get<2>(go) << endl;
            FSTStream.close();

        } else cout << "Error: No such a directory for\t" << odir + "fraction_stress_temperature.txt"s << endl;
//        for (auto lk : fraction_stress_temperature) if (get<0>(lk)*CellNumbs.at(2)*Burgv/DCC_size > 0.002) { cout << "Nano-slips fraction\t" << get<0>(lk) << "\tYield strength\t" << get<1>(lk) << "\tTemperature\t" << get<2>(lk) << endl;              break; }

    } ///End of 'Plasticity' type simulations
    else if (*stype == 'F') { // Fracture
        face_sequence = DCC_Kinetic_cracking(s_faces_sequence, CellNumbs, AFS, FES, paths, input_folder, odir);

    } ///End of 'Fracture' type simulations
    else { cout << "ERROR [DCC_Kinetic] : unknown simulation type - please replace with 'W or P' " << endl; return face_sequence;}
    //Output streams

/// Closing and deleting

    return face_sequence;
} /// The end of DCC_Kinetic()

/// ================================== Related functions ==================================///
/* std::vector<int> confCout(char* config) {
    std::string line;
    std:vector<int> res;
    ifstream inConf(config);
    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n')) {
            if(line.at(0) == '#') res.push_back(1); // 1 and # means accept - the parameter will be calculated
            if(line.at(0) == '%') res.push_back(0); // 0 and % means ignore - the parameter will not be calculated
        }
    } else cout << "The file " << config << " cannot be read" << endl; // If something goes wrong

    cout << "External force effect:                                         "; if (res.at(1) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Internal gradient force (between adjacent grains) effect:      "; if (res.at(2) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Chi-factor (Xi) effect:                                        "; if (res.at(3) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "GBs and TJs types effect:                                      "; if (res.at(4) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain boundary dislocations effect:                            "; if (res.at(5) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;

    return res;
}
 */
/**

 ///=============================================================================================================================================////
/// ====================================================>  Random generation process   <========================================================////
///=============================================================================================================================================////
        std::vector<bool> SpecialCells(CellNumbs.at(2), 0), OrdinaryCells(CellNumbs.at(2),1); // New vectors initialised with 0 for SpecialCells and 1 for OrdinaryCells

        map<unsigned int,unsigned int> SpecialCellMap; // Mapping [k]->[l] from the set of all 2-Cells in the initial DCC to the set of newly generated special cells
        map<unsigned int, unsigned int>::iterator sit; // Special iterator for this map
        double ordinary_faces_fraction = 1.0, special_faces_fraction = 0.0;
        std::vector<unsigned int> SpecialCellNumbs(CellNumbs.at(2), 0); // Vector of the size equal to the total number of faces in DCC initialised with '0's
        for (unsigned int lit = 0; lit < SpecialCellNumbs.size(); lit++) SpecialCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

        /// Numerates newly created Faces during the random generation process
        long unsigned int numerator = 0;
        /// Maximal fraction for simulation loop max_sFaces_fraction = [0,1]
        double max_sFaces_fraction = 0.5;
        /// Step for output (structural analysis and data output will be performed after each output_step calculation steps (= number of newly converted elements))
        int output_step = 300;
///=============================================================================================================================================////
/// ================= Loop over all possible fractions of special cells [0,1] =======================>

/// The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)
        srand (time(NULL));

        do { // do{ ... }while(output_step) loop starting point
            int New2CellNumb = 0, NewFaceType = 1; // Only one possible Face type (binary model)

            New2CellNumb = rand() % (SpecialCellNumbs.size()-1); // Random generation of the boundary number in the range from 0 to SpecialCellNumbs.size()-1
            if (number_of_types > 1) NewFaceType = rand() % number_of_types; // Random chose for the chosen boundary to be assigned over all special types
            OrdinaryCells.at(SpecialCellNumbs.at(New2CellNumb)) = 0; // Replace the chosen element with 0 instead of 1 in the Original Faces vector

            SpecialCellMap[SpecialCellNumbs.at(New2CellNumb)] = numerator++; // Assign new elements to the map and increase numerator
            long unsigned int sit = SpecialCellMap[SpecialCellNumbs.at(New2CellNumb)]; // Just a useful variable for the number of newly converted Face (= numerator--)

            /// It is a key tricky point for the fast Random generation process: vector decreases after each turn from an ordinary to special face BUT we chose all the boundary one by one because their initial numeration stored in the SpecialCellNumbs[] vector !
            SpecialCellNumbs.erase(SpecialCellNumbs.begin() + New2CellNumb); // !!! Delete its element from the vector decreasing its size BUT

            /// Another key point: creation of a triplet list SFaces_Triplet_list for the sparse matrix of special faces
            for(unsigned int j = 0; j < CellNumbs.at(2); j++) //  Loop over all the Faces in the DCC
                if(j != sit && AFS.coeff(sit,j) == 1 && (SpecialCellMap.find(j) != SpecialCellMap.end())) // if (1) not the same element (2) we find an adjacent element (neighbour) (3)! this element already added in the special faces map (= was randomly chosen)
                    SFaces_Triplet_list.push_back(Tr(sit,SpecialCellMap[j],1)); // then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}

            // Special and Ordinary Faces fraction calculation
            unsigned int OCellAmount = std::count(OrdinaryCells.begin(), OrdinaryCells.end(), 1);
            ordinary_faces_fraction = OCellAmount / (double) CellNumbs.at(2);
            special_faces_fraction = 1.0 - ordinary_faces_fraction;

            ///=============================================================================================================================================////
            ///=================================== Characterisation module =============================///

            if( OCellAmount % output_step == 0 && ordinary_faces_fraction > 0.05) { // Periodicity of characterisation output

                /// Creation of the sparse adjacency (SAM_FacesGraph) matrix for special Faces /// Sparse adjacency matrix from the corresponding triplet list of special faces
                SpMat Id(numerator, numerator), SAM_FacesGraph(numerator, numerator), Face_Laplacian(numerator, numerator), Sym_Face_Laplacian(numerator, numerator), RW_Face_Laplacian(numerator, numerator);
                SAM_FacesGraph.setFromTriplets(SFaces_Triplet_list.begin(), SFaces_Triplet_list.end());

                /// Degree vector and matrix
                vector<unsigned int> SFace_degrees;
                for (unsigned int i = 0; i < numerator; i++) { //columns
                    unsigned int degree_Fcounter = 0;
                    for (unsigned int j = 0; j < numerator; j++) { // rows
                        if (SAM_FacesGraph.coeff(j, i) == 1 && i != j)
                            degree_Fcounter++;
                    }
                    SFace_degrees.push_back(degree_Fcounter);
                }
                //Creation of the S-Face degree matrix
                SpMat SFDegree(SFace_degrees.size(), SFace_degrees.size());
                SFDegree.setIdentity();
                unsigned int numerator_degree = 0;
                for (double num: SFace_degrees) {
                    SFDegree.coeffRef(numerator_degree,numerator_degree) = SFace_degrees.at(numerator_degree);
                    numerator_degree++;
                }

                if (configuration.at(1)) {  // External force effect
                }
                if (configuration.at(2)) {  // Internal gradient force effect
                }
                if (configuration.at(3)) { // Chi-factor - orientation effect
                }
                if (configuration.at(4)) { // GBs and TJs types effect
                }
                if (configuration.at(5)) { // Grain boundary dislocations effect
                }

//                cout << "Average Strain:" << MStrain << ";\t Size of S-Face Adjacency matrix \t" << SAM_FacesGraph.nonZeros() << endl;
            } // End of analysis and output iterator ( IF: iterator % X == 0 )
        }while(ordinary_faces_fraction > (1.0 - max_sFaces_fraction)); /// End of the Random generation process


 Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols) {

    Eigen::SparseMatrix<double> res(Rows,Cols);
    typedef Triplet<double> Tr; // Eigen library class
    std::vector<Tr> tripletList; // Probe vector of triplets

    int i = 0, j = 0, value = 0, t_length = 0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong
//Sparse AN matrix
    res.setFromTriplets(tripletList.begin(), tripletList.end());

//Remove all elements anf free the memory from the probe vector
    tripletList.clear();
    tripletList.shrink_to_fit();

    return res;
}

std::vector<unsigned int> VectorReader(char* FilePath) {
    std::vector<unsigned int> res;
    unsigned int i=0;
    ifstream inCellNumbers(FilePath);
    if (inCellNumbers.is_open()) { //If the file was successfully open, then
        while(inCellNumbers >> i) res.push_back(i);
    } else std:cout << "The file " << FilePath << " cannot be read" << endl; //If something goes wrong

    return res;
}
**/