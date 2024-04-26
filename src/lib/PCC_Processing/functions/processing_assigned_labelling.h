#ifndef PCC_PROCESSING_DESIGN_PROCESSING_ASSIGNED_LABELLING_H
#define PCC_PROCESSING_DESIGN_PROCESSING_ASSIGNED_LABELLING_H

/*! ## 1 ##
 * @brief Basic Random generation machine for a new k-cell number and many other purposes: quasi-random choice of the element with # New2CellNumb from the list of numbers from 0 to OCellsNumb. OCellsNumb is the range of integer random numbers equal to the number of 'ordinary' k-cells.
 * @param OCellsNumb
 * @return random number
 */
unsigned int NewCellNumb_R(unsigned int OCellsNumb);

/*! ## 2 ##
 * @brief Basic Random Walks generation machine (with possible 'leaps') for a new k-cell chains of labels.
 * Leap_distance = 1 // for the future possibility of leaps: probability (hence frequency) and distance ((?)fraction of the complex size in grains).
 * @param iniFaceNumber
 * @param strip_length
 * @param Leap_frequency
 * @param Leap_distance
 * @return RW series of random numbers
 */
std::vector<unsigned int> NewCellsStrip_RW(int cell_type, unsigned int iniCellNumber, unsigned int strip_length, int Leap_frequency = 1, double Leap_distance = 1);

/*! ## 3 ##
 * @brief Log-normal distribution generator used for the RStrips_Distribution() processing function.
 * Its output vector contains the numbers of strips of each length starting with 1 expressed in # of faces.
 * Example: std::vector<int> strip_distribution = {2 4 5 27 8 6 3 1} means 2 strips of length 1 faces each, 4 strips of length 2 faces each,... , 1 strip of length 8 faces each.
 * @param mu_f
 * @param sigm_f
 * @param bins_number
 * @return strip_distribution_vector
 */
std::vector<double> Log_normal_distribution (double mu_f, double sigm_f, int bins_number);

/*! ## 4 ##
 * @brief The Random generation process function: quasi-random choice of the element with # New2CellNumb from the list of numbers {0 to OCellsNumb}
 * @param cell_type
 * @param Configuration_State
 * @param max_fractions_vectors
 * @return random state_vector of k-cells
 */
std::vector<std::vector<unsigned int>> Processing_Random(int cell_type, std::vector<std::vector<unsigned int>> &Configuration_State, std::vector<std::vector<double>> max_fractions_vectors, bool multiplexity);

/*! ## 5 ##
 * @brief Generate series special k-cells (series of strips or chains of k-cells) with the distribution of lengths taken from the 'face_strip_distribution' vector.
 * @param cell_type
 * @param cell_strip_distribution
 * @param Configuration_State
 * @param max_fractions_vectors
 * @return state_vector of k-cells
 */
std::vector<std::vector<unsigned int>> Processing_Random_Strips(int cell_type, std::vector<unsigned int> &cell_strip_distribution, std::vector<std::vector<unsigned int>> &Configuration_State, std::vector<std::vector<double>> max_fractions_vectors);

/*! ## 6 ##
 * @brief
 * @param cell_type
 * @param Configuration_State
 * @param max_fractions_vectors
 * @param p_index
 * @return state_vectors of k-cells
 */
std::vector<unsigned int> Processing_maxFunctional(int cell_type, std::vector<std::vector<int>> &Configuration_State, std::vector<std::vector<double>> const &max_fractions_vectors, bool multiplexity);

/*! ## 7 ##
 * @brief
 * @param cell_type
 * @param Configuration_State
 * @param max_fractions_vectors
 * @param p_index
 * @return
 */
std::vector<unsigned int> Processing_minConfEntropy(int cell_type, std::vector<std::vector<int>> &Configuration_State, std::vector<std::vector<double>> const &max_fractions_vectors, double p_index);

/*!
 * @brief
 * @param cell_type
 * @param Configuration_State
 * @param max_fractions_vectors
 * @param p_index
 * @return
 */
std::vector<unsigned int> Processing_maxF_crystallographic(int cell_type, std::vector<std::vector<int>> &Configuration_State, std::vector<std::vector<double>> const &max_fractions_vectors, double p_index);
// std::vector<unsigned int> Processing_maxP_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_sState, std::vector<vector<double>> const &max_fractions_vectors, double p_index);
// std::vector<unsigned int> Processing_Random_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_sState, std::vector<vector<double>> const &max_fractions_vectors, double p_index);

#endif //PCC_PROCESSING_DESIGN_PROCESSING_ASSIGNED_LABELLING_H
