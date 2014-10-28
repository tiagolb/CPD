#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <omp.h>

/* DEBUG */

//#define DEBUG
#define DEBUG_TIME

/* RETURN CONSTANTS */

#define ERROR 1
#define SUCCESS 0

/* GENERAL CONSTANTS */

#define N_ARGUMENTS 2
#define ARGUMENT 1
#define PRINT_LIMIT 41

#define FIRST_TRI_OFFSET 2
#define ZEROS_OFFSET 1

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

void lcsPrintMatrix(std::vector< std::vector<int> > & matrix) {
	std::cout << "line: " << matrix.size() << ", cols: " << matrix[0].size() << std::endl;
	for(size_t i = 0; i < matrix.size(); i++) {
		for(size_t j = 0; j < matrix[i].size(); j++) {
			std::cout << "|" << matrix[i][j];
		}
		std::cout << "|" << std::endl;
	}
}

void calcMatrixCell(size_t i, size_t j, std::vector< std::vector<int> > & matrix,
					std::string & seq1, std::string & seq2) {
	if(seq1[i-1] == seq2[j-1]) {
		matrix[i][j] = matrix[i-1][j-1] + cost(i);
	} else {
		matrix[i][j] = std::max(matrix[i-1][j], matrix[i][j-1]);
	}
}

std::vector< std::vector<int> > lcsPopulateMatrix_line(std::string & seq1, std::string & seq2) {
	size_t lines = seq1.size()+1;
	size_t cols = seq2.size()+1;
	std::vector< std::vector<int> > matrix(lines, std::vector<int>(cols, 0));
	
	#pragma omp parallel
	{
		// #ifdef DEBUG_TIME
		// #pragma omp single
		// std::cout << omp_get_num_threads() << std::endl;
		// #endif

		for(size_t diag = 0; diag < lines - FIRST_TRI_OFFSET; diag++) {
			#pragma omp for
			for(size_t line = diag + ZEROS_OFFSET; line > 0; line--) {
				size_t col = diag - line + FIRST_TRI_OFFSET;
				calcMatrixCell(line, col, matrix, seq1, seq2);
			}
		}

		for(size_t diag = 0; diag < cols - lines; diag++) {
			#pragma omp for
			for(size_t line = lines - ZEROS_OFFSET; line > 0; line--) {
				size_t col = lines - line + diag;
				calcMatrixCell(line, col, matrix, seq1, seq2);
			}
		}

		for(size_t diag = 0; diag < lines - ZEROS_OFFSET; diag++) {
			#pragma omp for
			for(size_t col = cols - lines + diag + ZEROS_OFFSET; col < cols; col++) {
				size_t line = cols - col + diag;
				calcMatrixCell(line, col, matrix, seq1, seq2);
			}
		}
	}

	return matrix;

}

std::vector< std::vector<int> > lcsPopulateMatrix_col(std::string & seq1, std::string & seq2) {
	size_t lines = seq1.size()+1;
	size_t cols = seq2.size()+1;
	std::vector< std::vector<int> > matrix(lines, std::vector<int>(cols, 0));
	
	#pragma omp parallel
	{
		// #ifdef DEBUG_TIME
		// #pragma omp single
		// std::cout << omp_get_num_threads() << std::endl;
		// #endif
		
		for(size_t diag = 0; diag < cols - FIRST_TRI_OFFSET; diag++) {
			#pragma omp for
			for(size_t col = diag + ZEROS_OFFSET; col > 0; col--) {
				size_t line = diag - col + FIRST_TRI_OFFSET;
				calcMatrixCell(line, col, matrix, seq1, seq2);
			}
		}

		for(size_t diag = 0; diag < lines - cols; diag++) {
			#pragma omp for
			for(size_t col = cols - ZEROS_OFFSET; col > 0; col--) {
				size_t line = cols - col + diag;
				calcMatrixCell(line, col, matrix, seq1, seq2);
			}
		}

		for(size_t diag = 0; diag < cols - ZEROS_OFFSET; diag++) {
			#pragma omp for
			for(size_t line = lines - cols + diag + ZEROS_OFFSET; line < lines; line++) {
				size_t col = lines - line + diag;
				calcMatrixCell(line, col, matrix, seq1, seq2);
			}
		}
	}

	return matrix;

}

/*
	duas funcoes, uma que ve a matrix pelas linhas quando
	seq1 > seq2 e outra que vê pelas colunas quando seq2 > seq1
*/

std::string lcsFindSubString(std::string & seq1, std::string & seq2,  std::vector< std::vector<int> > & matrix) {
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
	#ifdef DEBUG_TIME
	double start = omp_get_wtime();
	#endif
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

	std::vector< std::vector<int> > matrix;
	if(seq1.size() > seq2.size()) {
 		matrix = lcsPopulateMatrix_col(seq1, seq2);
	} else {
		matrix = lcsPopulateMatrix_line(seq1, seq2);
	}
	
	//lcsPrintMatrix(matrix);
	std::string subString = lcsFindSubString(seq1, seq2, matrix);
	
	std::cout << matrix[seq1.size()][seq2.size()] << std::endl;
	std::cout << subString << std::endl;

	#ifdef DEBUG_TIME
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;
	#endif

	return SUCCESS;
}
