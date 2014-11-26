#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include <mpi.h>

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

std::vector< std::vector<short> > lcsBuildMatrix(int id, int seq1_size, int seq2_size) {
	size_t rows = seq1_size + 1;
	size_t cols = (id == 0) ? seq2_size+1 : seq2_size;
	std::vector< std::vector<short> > matrix(rows, std::vector<short>(cols, 0));
	return matrix;
}

void calcMatrixCell(size_t i, size_t j, std::vector< std::vector<short> > & matrix,
					std::string & seq1, std::string & seq2) {
	if(seq1[i-1] == seq2[j-1]) {
		matrix[i][j] = matrix[i-1][j-1] + cost(i);
	} else {
		matrix[i][j] = std::max(matrix[i-1][j], matrix[i][j-1]);
	}
}

void lcsPopulateMatrixFirstCalc_col(std::vector< std::vector<short> > & matrix,
									std::string & seq1, std::string & seq2) {
}

void lcsPopulateMatrixMiddleCalc_col(std::vector< std::vector<short> > & matrix,
									std::string & seq1, std::string & seq2) {
}

void lcsPopulateMatrixLastCalc_col(std::vector< std::vector<short> > & matrix,
									std::string & seq1, std::string & seq2) {
}

void lcsPopulateMatrixFirstCalc_line(std::vector< std::vector<short> > & matrix,
									std::string & seq1, std::string & seq2) {
}

void lcsPopulateMatrixMiddleCalc_line(std::vector< std::vector<short> > & matrix,
									std::string & seq1, std::string & seq2) {
}

void lcsPopulateMatrixLastCalc_line(std::vector< std::vector<short> > & matrix,
									std::string & seq1, std::string & seq2) {
}

/* MAIN */

int main(int argc, char* argv[]) {
	#ifdef DEBUG_TIME
	double start = MPI_Wtime();
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

	MPI_Init(&argc, &argv);

	int id, p;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(id == 0) {
    	std::cout << "Este e' o numero de pc's: " << p << std::endl;
    	std::cout << "Este e' o numero de colunas: " << seq2.size() << std::endl;
    }

    int istart, iend;
    istart = id * seq2.size() / p;
    iend = (id+1) * seq2.size() / p;
    if(id == p-1) {
    	iend = seq2.size();
    }
    std::string substring = seq2.substr(istart, iend - istart);
    std::cout << "(" << id << ") minha string: " << substring << std::endl;

    std::vector< std::vector<short> > matrix = lcsBuildMatrix(id, seq1.size(), substring.size());
    std::cout << "(" << id << ") rows: " << matrix.size() << ", cols: " << matrix[0].size() << std::endl;

    if(id == 0) {
    	if(seq1.size() > substring.size()) {
 			lcsPopulateMatrixFirstCalc_col(matrix, seq1, substring);
		} else {
			lcsPopulateMatrixFirstCalc_line(matrix, seq1, substring);
		}

    } if (id == p-1) {
    	if(seq1.size() > substring.size()) {
 			lcsPopulateMatrixLastCalc_col(matrix, seq1, substring);
		} else {
			lcsPopulateMatrixLastCalc_line(matrix, seq1, substring);
		}
    } else {
    	if(seq1.size() > substring.size()) {
 			lcsPopulateMatrixMiddleCalc_col(matrix, seq1, substring);
		} else {
			lcsPopulateMatrixMiddleCalc_line(matrix, seq1, substring);
		}
    }

	#ifdef DEBUG_TIME
	double end = MPI_Wtime();
	std::cout << "time: " << end - start << std::endl;
	#endif

	MPI_Finalize();

	return SUCCESS;
}