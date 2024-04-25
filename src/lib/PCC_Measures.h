
std::vector<int> Edge_types_byFaces(std::vector<unsigned int> const &CellNumbs, std::vector<unsigned int> &special_face_sequence, std::vector<double> &j_fractions, std::vector<double> &d_fractions);

std::vector<double> j_fractions_vector(std::vector<int> const &TJsTypes); // based on Edges vector

std::vector<double> d_fractions_vector(std::vector<int> const &TJsTypes); // based on Edges vector

double Configuration_Entropy(std::vector<int> const &TJsTypes);

double Configuration_Entropy(std::vector<double> const &j_fractions); // 2nd  overloaded function

std::tuple<double, double> Configuration_Entropy_tuple(std::vector<double> const &j_fractions);

double Configuration_Entropy_change(std::vector<int> const &TJsTypes, std::vector<int> const &TJsTypesNew);

double Configuration_Entropy_change(std::vector<double> const &j_fractions, std::vector<int> const &TJsTypesNew); // 2nd overloaded function