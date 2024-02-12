///================================ DCC Edges types statistics, indices and configuration entropy =============================================================///
///==============================================================================================================================///
/** This subfile calculates the statistical characteristics of Edges incuding their Face indices and configurational entropy  **/
///==============================================================================================================================///

using namespace std; ///Standard namespace
vector<double> EdgesStat( std::vector<unsigned int> &s_faces_sequence, std::vector<unsigned int> const& CellNumbs, Eigen::SparseMatrix<double> const& FES, char* output_dir, double &Face_Entropy_Median, double &Face_Entropy_Skrew, double &informativeness, vector<double> &j_types_fractions)
{
    vector<double> TJsTypes(CellNumbs.at(1),0);
    map<unsigned int, unsigned int> res; // Here 100 is an arbitrary number of Edge types
    map<unsigned int, unsigned int>::iterator sit; // Special iterator for this map
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0;
    Face_Entropy_Median = 0; Face_Entropy_Skrew = 0;

    TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
//    for (auto sfe : TJsTypes) cout << sfe << "\t" ; cout << endl;

    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
    J0 = CellNumbs.at(1) - J1 - J2 - J3;
    Jall = CellNumbs.at(1);
// Conversion from numbers to fractions
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
    j0 = J0/Jall; j1 = J1/Jall; j2 = J2/Jall; j3 = J3/Jall;
    double j0s = j0, j1s = j1, j2s = j2, j3s = j3;
    /// using values with pow(10,-10) instead of 0s!
//    if (j0s == 0) j0s = pow(10,-30); if (j1s == 0) j1s = pow(10,-30); if (j2s == 0) j2s = pow(10,-30); if (j3s == 0) j3s = pow(10,-30); //Gives 0 in entropy!
    if (j0s != 0) j0s = j0* log2(j0); if (j1s != 0) j1s = j1* log2(j1); if (j2s != 0) j2s = j2* log2(j2); if (j3s != 0) j3s = j3* log2(j3); //Gives 0 in entropy!
    /// Configuration Entropy related with Faces
//    Configurational_Face_Entropy = - (j0s* log2(j0s) + j1s* log2(j1s) + j2s* log2(j2s) + j3s* log2(j3s));
    Configurational_Face_Entropy = - (j0s + j1s + j2s + j3s);

    /// Median part in the entropy decomposition
    j_types_fractions = {j0, j1, j2, j3}; /// using values with pow(10,-10) instead of 0s!
    if (j0s!=0 && j1s!=0 && j2s!=0 && j3s!=0) {
        Face_Entropy_Median = -(1.0 / j_types_fractions.size()) * log2(j0 * j1 * j2 * j3);
    } else Face_Entropy_Median = 0.0;

    /// Screw part (divergence from the uniform distribution -> S_max) in the entropy decomposition
    for (int j = 0; j < j_types_fractions.size(); j++)
        for (int i = 0; i < j; i++)
            if (j_types_fractions[i]!=0 && j_types_fractions[j]!=0) {
                Face_Entropy_Skrew +=
                        -(1.0 / j_types_fractions.size()) * (j_types_fractions[i] - j_types_fractions[j]) *
                        log2(j_types_fractions[i] / j_types_fractions[j]);
            } else Face_Entropy_Skrew += 0.0;

    /// Complexity/ Informativeness
    // Srand and Smax reading from files
    string input_filename_SstrRand = "Random_Entropy_100.txt"s, input_filename_SstrMAX = "Maximum_Entropy_100.txt"s;
    string input_RandomEntropy_dir = output_dir + input_filename_SstrRand, input_MAXEntropy_dir = output_dir + input_filename_SstrMAX;
    char* RandomEntropy_dir = const_cast<char*>(input_RandomEntropy_dir.c_str()); char* MAXEntropy_dir = const_cast<char*>(input_MAXEntropy_dir.c_str()); // From string to char for the passing folder path to a function
    // Format of the elements:: special Face fraction _ Conf entropy value _ Mean entropy value (log(p1*p2*..*pn))
    vector<vector<double>>  RandomEntropy = VectorVectors4Reader(RandomEntropy_dir);
    vector<vector<double>>  MAXEntropy = VectorVectors4Reader(MAXEntropy_dir);
    //for ( auto tpl : MAXEntropy) cout << tpl[0] << " " << tpl[1] << " " << tpl[2] << endl;

    /// Index of complexity:
    double RandomEntropy_p = 0, MAXEntropy_p = 0;
    vector<double> Delta_MAX;
    double sff =(double) s_faces_sequence.size()/ CellNumbs.at(2); // special face fraction
    for ( auto tpl : MAXEntropy)  Delta_MAX.push_back(abs(sff - tpl.at(0)));
    auto numb_Smax = std::min_element(std::begin(Delta_MAX), std::end(Delta_MAX)) - std::begin(Delta_MAX); // gives index of the max element
    Delta_MAX.clear();
    for ( auto tpl : RandomEntropy)  Delta_MAX.push_back(abs(sff - tpl.at(0)));
    auto numb_Srand = std::min_element(std::begin(Delta_MAX), std::end(Delta_MAX)) - std::begin(Delta_MAX); // gives index of the max element
    Delta_MAX.clear();
    /// Informativeness parameter
//REPAIR    cout << numb_Smax << " " << numb_Srand << endl;
//REPAIR    cout << Configurational_Face_Entropy << " " << RandomEntropy[numb_Srand][1] << " " << MAXEntropy[numb_Smax][1] << endl;
    informativeness =  ( Configurational_Face_Entropy - RandomEntropy[numb_Srand][1])/ ( MAXEntropy[numb_Smax][1] - RandomEntropy[numb_Srand][1]);
    if (informativeness > 1) informativeness = 1;
    /// ====== Data output =====================>
    /// Opening of the output streams
    string TJs_output_filename = "TJsLab_TJsTypes.txt"s, Entropy_output_filename = "TJsLab_ConTJsEntropy.txt"s,
            output_TJs_dir = output_dir + TJs_output_filename, output_Entropy_dir = output_dir + Entropy_output_filename;
    char* cTJs_dir = const_cast<char*>(output_TJs_dir.c_str()); char* cEntropy_dir = const_cast<char*>(output_Entropy_dir.c_str()); // From string to char for the passing folder path to a function

    ofstream OutTJsFile; OutTJsFile.open(cTJs_dir, ios::app);
    double special_faces_fraction =(double) s_faces_sequence.size()/ CellNumbs.at(2);
//    cout << Configurational_Face_Entropy << "\t" << Face_Entropy_Median << "\t"<< Face_Entropy_Median << endl;

    OutTJsFile << special_faces_fraction << "\t" << j0 << "\t" << j1 << "\t" << j2 << "\t" << j3 << Configurational_Face_Entropy << "\t" << Face_Entropy_Median << "\t\t" << Face_Entropy_Skrew << endl;
    OutTJsFile.close();

//    int b = 0;
//    for (unsigned int it : TJsTypes) res[b++] = it;

    return TJsTypes;
}

/*
 * ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////// Обработка матриц соседей: начльная статистика связей и матрицы смежности /////////////////////////////////////////////////////
	long lt=0; knn=1;
// Для каждой грани мы строим матрицу JEdgeNeigh[i][j] где на первом месте ее собственный J-тип а дальше список J-типов ее соседей
long knumb=0, J0Ncount=0, J00Ncount=0; 	knn=0;
	int **JEdgeNeigh;
	JEdgeNeigh = new int* [Edgenumb];
	for (int i = 0; i < Edgenumb; i++) JEdgeNeigh[i] = new int [100];	//снова закладываем до 100 гипотетических соседей
	int **JFaceNeigh;
	JFaceNeigh = new int* [Facenumb];
	for (int i = 0; i < Facenumb; i++) JFaceNeigh[i] = new int [100];	//снова закладываем до 100 гипотетических соседей

		for (int i = 0; i < Edgenumb; i++)
			for (int j = 0; j < 100; j++)
				JEdgeNeigh[i][j] = -6;
		for (int i = 0; i < Facenumb; i++)
			for (int j = 0; j < 100; j++)
				JFaceNeigh[i][j] = -1; //Все LAGBs!!!

//Сначала мы выясняем тип самой грани
for(int i = 0; i < Edgenumb; i++)	{
	for (int lk = 0; lk < Facenumb; lk++) if(MFE1[lk][i] == 1) J00Ncount++; //then ==2
		JEdgeNeigh[i][0] = J00Ncount;
			J00Ncount=0;
	lt=0; knn=1;


//Затем проходим всех ее соседей
	do{
		if(EdgeNeighbours[i][lt]>=0) knumb = EdgeNeighbours[i][lt]; else knumb=-1;
//------------------------------------------------------

		if(knumb>=0) for (int lk = 0; lk < Facenumb; lk++) if(MFE1[lk][knumb] == 1) J0Ncount++;  //then ==2

//-----------------------------------------------------		//cout <<i<< "   "<< knn << "   "<< JEdgeNeigh[0][1] << endl;
		JEdgeNeigh[i][knn++] = J0Ncount;

		J0Ncount=0;	lt++;
	}while(knumb>=0);
}

//Вывод в файл JEdgeNeigh.txt
			JEdgeN.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)JEdgeNeigh.txt", ios::trunc);
			for (int i = 0; i < Edgenumb; i++) {
				for (int j = 0; j < 100; j++) {
					if(JEdgeNeigh[i][j]>=0) JEdgeN << JEdgeNeigh[i][j] << "\t";
				}
			JEdgeN <<"Edge number="<<i<<"\n";
		}
        JEdgeN.close();

//***cout<<"initial (HAGBs)JEdgeNeigh.txt has been created"<<endl;

//Чистка файла под мощности тройных стыков с учетом соседей
	JEdgeN.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)TJspow.txt", ios::trunc);
	JEdgeN.close();


	int **JEN;
	JEN = new int* [30];
	for (int i = 0; i < 30; i++) JEN[i] = new int [30];
			for (int i = 0; i < 5; i++)
				for (int j = 0; j < 5; j++)
					JEN[i][j] = -10;

//Анализ матрицы JEdgeNeigh[i][j]
for (int i = 0; i < Edgenumb; i++) //для каждой грани
	for (int l = 0; l < 5; l++) //мы перебираем все варианты какой она может быть (с запасом до 30, хотя реально до 4-5 видов)
		if(JEdgeNeigh[i][0] == l) for (int j = 1; j < 100; j++) //и если она оказалась определенного типа, то мы перебираем всех ее соседей
										for (int k = 0; k < 5; k++)
											if(JEdgeNeigh[i][j] == k) JEN[l][k]++; //так что если сосед оказывается также определенного типа, то мы заносим их связь в матрицу JEN
//очевидно, что при таком алгоритме диагональные элементы учитываются дважды, то есть
for (int l = 0; l < 5; l++)
	for (int k = 0; k < 5; k++)
		if(k == l)  JEN[l][k]= 0.5*JEN[l][k];

//Вывод в файл начального распределения узлов по количеству соседей разного типа JEN.txt
//(сколько 1-узлов соединено с 1-узлами, 2 - узлами..., сколько 3-узлов соединено с 0 - узлами, 1 - узлами и тд)
			JENStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)JEN00.txt", ios::trunc);
			for (int i = 0; i < 5; i++) {
				for (int j = 0; j < 5; j++) {
					if(JEN[i][j]>=0) JENStream << JEN[i][j] << "\t";
//					if(JEN[i][j]>=0) JENStream << JEN[i][j]/Edgenumb << "\t"; //В долях /Edgenumb
				}
				JENStream <<"\n";
		}
        JENStream.close();

		JENStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)Jlk.txt", ios::trunc);  //			 JENStream << JEN[0][0] + JEN[1][1] + JEN[2][2] + JEN[3][3]<<"\t"<< JEN[0][1] + JEN[1][0] + JEN[0][2] + JEN[2][0] + JEN[0][3] + JEN[3][0] + JEN[1][2] + JEN[2][1] + JEN[1][3] + JEN[3][1] + JEN[2][3] + JEN[3][2]<< "\t"<< JEN[0][1] + JEN[1][0]<< "\t"<< JEN[0][2] + JEN[2][0]<< "\t"<< JEN[0][3] + JEN[0][3]<< "\t"<< JEN[1][2] + JEN[2][1]<< "\t"<< JEN[1][3] + JEN[3][1]<< "\t"<< JEN[2][3] + JEN[3][2];
		JENStream.close();

 *
 */