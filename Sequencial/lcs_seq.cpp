#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

/* DEBUG */

//#define DEBUG

/* RETURN CONSTANTS */

#define ERROR 1
#define SUCCESS 0

/* GENERAL CONSTANTS */

#define N_ARGUMENTS 2
#define ARGUMENT 1
#define PRINT_LIMIT 41

/* FUNCTIONS */

int validateArguments(int argc) {
	
	if(argc != N_ARGUMENTS) {
		#ifdef DEBUG
		printf("#ERROR - One and only one argument is allowed (file name).\n");
		#endif
		return ERROR;
	}
	return SUCCESS;
}

/* cost FUNCTION as provided by the teachers */
short cost(int x) {
	int i, n_iter = 20;
	double dcost = 0;
	for(i = 0; i < n_iter; i++) {
		dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
	}
	return (short) (dcost / n_iter + 0.1);
}

void lcsPrintMatrix(std::vector< std::vector<int> > matrix) {
	for(size_t i = 0; i < matrix.size(); i++) {
		for(size_t j = 0; j < matrix[i].size(); j++) {
			std::cout << "|" << matrix[i][j];
		}
		std::cout << "|" << std::endl;
	}
}

std::vector< std::vector<int> > lcsPopulateMatrix(std::string seq1, std::string seq2) {
	size_t rows = seq1.size()+1;
	size_t cols = seq2.size()+1;
	std::vector< std::vector<int> > matrix(rows, std::vector<int>(cols, 0));
	for(size_t i = 1; i < rows; i++) {
		for(size_t j = 1; j < cols; j++) {
			if(seq1[i-1] == seq2[j-1]) {
				matrix[i][j] = matrix[i-1][j-1] + cost(i);
			} else {
				matrix[i][j] = std::max(matrix[i-1][j], matrix[i][j-1]);
			}
		}
	}
	return matrix;
}

std::string lcsFindSubString(std::string seq1, std::string seq2,  std::vector< std::vector<int> > matrix) {
	int row = seq1.size(), col = seq2.size();
	std::string result = "";
	while(matrix[row][col] != 0) {
		if(seq1[row-1] == seq2[col-1]) {
			result.insert(result.begin(), seq1[row-1]);
			row--;
			col--;
		} else {
			if(matrix[row-1][col] > matrix[row][col-1]) {
				row--;
			} else {
				col--;
			}
		}
	}

 	return result;
}

/* MAIN */

int main(int argc, char* argv[]) {
	// this line allows for streams to be faster than buffer input
	std::ios_base::sync_with_stdio (false);

	if(validateArguments(argc) == ERROR) {
		#ifdef DEBUG
		std::cout << "#ERROR - Argument is not valid." << std::endl;
		#endif
		return ERROR;
	}

	std::string filename = argv[ARGUMENT];
	#ifdef DEBUG
	std::cout << "#META - The file name is: " << filename << std::endl;
	#endif

	std::ifstream inputFile;
  	inputFile.open(filename.c_str(), std::ifstream::in);

  	std::string seq1, seq2;
  	std::getline(inputFile, seq1); //reads the first line and ignores it (no need to read these values)
  	std::getline(inputFile, seq1);
  	std::getline(inputFile, seq2);
  	inputFile.close();

	#ifdef DEBUG
	std::cout << "#META - Sequence1: " << seq1 << std::endl;
	std::cout << "#META - Sequence2: " << seq2 << std::endl;
	#endif

	std::vector< std::vector<int> > matrix = lcsPopulateMatrix(seq1, seq2);
	std::string subString = lcsFindSubString(seq1, seq2, matrix);

	//size_t rows = seq1.size()+1;
	//size_t cols = seq2.size()+1;
	//std::vector< std::vector<int> > matrix(rows, std::vector<int>(cols, 0));

	#ifdef DEBUG
	lcsPrintMatrix(matrix);
	#endif
	
	std::cout << matrix[seq1.size()][seq2.size()] << std::endl;
	std::cout << subString << std::endl;

	return SUCCESS;
}
