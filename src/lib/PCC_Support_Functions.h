#ifndef PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H
#define PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H

#include <Eigen/SparseCore>

/// # 1 # Bool check if the file 'fileName' (path) exists in the directory
bool is_file_exists(const std::string fileName);

/// # 2 # Creates int std::vector from the file (containing in the directory) 'FilePath'
std::vector<unsigned int> VectorIReader(const char* FilePath); // creation int std::Vector from file

/// # 3 # Creates double std::vector from the file (containing in the directory) 'FilePath'
std::vector<double> VectorDReader(const char* FilePath);

/// # 4 # Creation Eigen::Sparse_Matrix from file
//Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols);
Eigen::SparseMatrix<double> SMatrixReader(std::string SMpath, unsigned int Rows, unsigned int Cols);

/// # 5 # Finding barycenter coordinates as a tuple<double, double, double> for a given 'facenumb' face
std::tuple<double, double, double> find_aGBseed(unsigned int facenumb, std::vector<char*> const paths, std::vector<unsigned int> & CellNumbs, std::vector<std::tuple<double, double, double>> & AllSeeds_coordinates);

template <typename TP> inline constexpr
int sign(TP x, std::false_type is_signed);
template <typename TP> inline constexpr
int sign(TP x, std::true_type is_signed);
template <typename TP> inline constexpr
int sign(TP x);

/*!
 * @brief  Face indexing based on the preliminsry indrxing of their edges (vector<double> TJsTypes)
 * @param face_number
 * @param FES
 * @param TJsTypes
 * @return
 */
std::vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, std::vector<double> const& TJsTypes);

/*!
 * @brief Node Types Calculation function
 * @param CellNumbs
 * @param s_faces_sequence
 * @param ENS
 * @return
 */
std::vector<double> NodesTypesCalc(std::vector<unsigned int> const &CellNumbs, std::vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const &ENS);

/*!
 * @brief  Edge Types Calculation function
 * @param CellNumbs
 * @param s_faces_sequence
 * @param FES
 * @return
 */
std::vector<double> EdgesTypesCalc(std::vector<unsigned int> const &CellNumbs, std::vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const &FES);


#endif //PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H