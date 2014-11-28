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

#define FIRST_TRI_OFFSET 2
#define ZEROS_OFFSET 1
#define SEND_OFFSET 1

#define FIRST_TRI_MIDDLE_OFFSET 1
#define TAG_ID 841

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
	size_t cols = seq2_size+1;
	std::vector< std::vector<short> > matrix(rows, std::vector<short>(cols, 0));
	return matrix;
}

short lcsPopulateMatrixFirst(int id, int from, std::vector< std::vector<short> > & matrix,
					   std::string & seq1, std::string & seq2) {
	size_t rows = matrix.size();
	size_t cols = matrix[0].size();
	MPI_Request myRequest;
	for(size_t i = 1; i < rows; i++) {
		for(size_t j = 1; j < cols; j++) {
			if(seq1[i-1] == seq2[j-1]) {
				matrix[i][j] = matrix[i-1][j-1] + cost(i);
			} else {
				matrix[i][j] = std::max(matrix[i-1][j], matrix[i][j-1]);
			}
		}
		
		MPI_Isend(&matrix[i][cols-1], 1, MPI_SHORT, 1, i, MPI_COMM_WORLD, &myRequest);
		// como nao vamos escrever em matrix[i][j] de novo nao e' necessario MPI_Wait aqui
	}

	short result;
	MPI_Status status;
	MPI_Recv(&result, 1, MPI_SHORT, from, rows, MPI_COMM_WORLD, &status);
	return result;
}

void lcsPopulateMatrixMiddle(int id, std::vector< std::vector<short> > & matrix,
					   std::string & seq1, std::string & seq2) {
	int from = id - 1;
	int dest = id + 1;
	size_t rows = matrix.size();
	size_t cols = matrix[0].size();
	MPI_Status status;
	MPI_Request myRequest;
	for(size_t i = 1; i < rows; i++) {
		MPI_Recv(&matrix[i][0], 1, MPI_SHORT, from, i, MPI_COMM_WORLD, &status);
		
		for(size_t j = 1; j < cols; j++) {
			if(seq1[i-1] == seq2[j-1]) {
				matrix[i][j] = matrix[i-1][j-1] + cost(i);
			} else {
				matrix[i][j] = std::max(matrix[i-1][j], matrix[i][j-1]);
			}
		}

		MPI_Isend(&matrix[i][cols-1], 1, MPI_SHORT, dest, i, MPI_COMM_WORLD, &myRequest);
		// como nao vamos escrever em matrix[i][j] de novo nao e' necessario MPI_Wait aqui
	}
}

void lcsPopulateMatrixLast(int id, std::vector< std::vector<short> > & matrix,
					   std::string & seq1, std::string & seq2) {
	int from = id - 1;
	size_t rows = matrix.size();
	size_t cols = matrix[0].size();
	MPI_Status status;
	for(size_t i = 1; i < rows; i++) {

		MPI_Recv(&matrix[i][0], 1, MPI_SHORT, from, i, MPI_COMM_WORLD, &status);

		for(size_t j = 1; j < cols; j++) {
			if(seq1[i-1] == seq2[j-1]) {
				matrix[i][j] = matrix[i-1][j-1] + cost(i);
			} else {
				matrix[i][j] = std::max(matrix[i-1][j], matrix[i][j-1]);
			}
		}
	}

	MPI_Request myRequest;
	MPI_Isend(&matrix[rows-1][cols-1], 1, MPI_SHORT, 0, rows, MPI_COMM_WORLD, &myRequest);
	// como nao vamos escrever em matrix[i][j] de novo nao e' necessario MPI_Wait aqui
}

std::string lcsFindSubStringMaster(int processors, int id, std::string & seq1, std::string & seq2,
								   std::vector< std::vector<short> > & matrix) {

	MPI_Status status;
	std::string result = "";
	char * buffer = NULL;
	for(int p = processors-1; p > 0; p--) {
		int len = 0;
		MPI_Probe(p, p, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_CHAR, &len);
		buffer = new char[len];
		MPI_Recv(buffer, len, MPI_CHAR, p, p, MPI_COMM_WORLD, &status);
		result.insert(0, std::string(buffer));
		delete buffer;
	}

	int row, col = seq2.size();
	MPI_Recv(&row, 1, MPI_SHORT, id+1, id+1+TAG_ID, MPI_COMM_WORLD, &status);

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

void lcsFindSubStringSlave(int id, std::string & seq1, std::string & seq2,
						   std::vector< std::vector<short> > & matrix) {
	int row, col = seq2.size();
	MPI_Status status;
	MPI_Recv(&row, 1, MPI_SHORT, id+1, id+1+TAG_ID, MPI_COMM_WORLD, &status);
	
	std::string result = "";
	while(matrix[row][col] != 0 && col > 0) {
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
	MPI_Send(&row, 1, MPI_SHORT, id-1, id+TAG_ID, MPI_COMM_WORLD);
	
	char * buffer = new char[result.size()+1];
	result.copy(buffer, result.size(), 0);
	buffer[result.size()] = '\0';

	MPI_Send(buffer, result.size()+1, MPI_CHAR, 0, id, MPI_COMM_WORLD);
}

void lcsFindSubStringLastSlave(int id, std::string & seq1, std::string & seq2,
							   std::vector< std::vector<short> > & matrix) {
	int row = seq1.size(), col = seq2.size();
	std::string result = "";
	while(matrix[row][col] != 0 && col > 0) {
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
	MPI_Send(&row, 1, MPI_SHORT, id-1, id+TAG_ID, MPI_COMM_WORLD);
	
	char * buffer = new char[result.size()+1];
	result.copy(buffer, result.size(), 0);
	buffer[result.size()] = '\0';

	MPI_Send(buffer, result.size()+1, MPI_CHAR, 0, id, MPI_COMM_WORLD);
}

/* PRINT */

void lcsPrintMatrix(std::vector< std::vector<short> > & matrix) {
	for(size_t i = 0; i < matrix.size(); i++) {
		for(size_t j = 0; j < matrix[i].size(); j++) {
			std::cout << "|" << matrix[i][j];
		}
		std::cout << "|" << std::endl;
	}
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

    #ifdef DEBUG
    if(id == 0) {
    	std::cout << "Este e' o numero de pc's: " << p << std::endl;
    	std::cout << "Este e' o numero de colunas: " << seq2.size() << std::endl;
    }
    #endif

    std::string substring;
    char * buffer;
    MPI_Status status;

    if(id == 0) {
	    int first_start = id * seq2.size() / p;
	    int first_end = (id+1) * seq2.size() / p;
	    substring = seq2.substr(first_start, first_end - first_start);

	    int istart, iend;
	    for(int i = 1; i < p; i++) {
	    	istart = i * seq2.size() / p;
	    	iend = (i+1) * seq2.size() / p;
	    	if(i == p-1) {
	    		iend = seq2.size();
	    	}
	    	std::string sendSubstring = seq2.substr(istart, iend - istart);
	    	buffer = new char[sendSubstring.size()+1];
			sendSubstring.copy(buffer, sendSubstring.size(), 0);
			buffer[sendSubstring.size()] = '\0';
			MPI_Send(buffer, sendSubstring.size()+1, MPI_CHAR, i, i, MPI_COMM_WORLD);
	    } 
    } else {
    	int len = 0;
		MPI_Probe(0, id, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_CHAR, &len);
		buffer = new char[len];
		MPI_Recv(buffer, len, MPI_CHAR, 0, id, MPI_COMM_WORLD, &status);
		substring = std::string(buffer);
    }
    delete buffer;

    std::vector< std::vector<short> > matrix = lcsBuildMatrix(id, seq1.size(), substring.size());

    if(id == 0) {
    	short result = lcsPopulateMatrixFirst(id, p-1, matrix, seq1, substring);
    	std::string resultSubstring = lcsFindSubStringMaster(p, id, seq1, substring, matrix);

    	std::cout << result << std::endl;
    	std::cout << resultSubstring << std::endl;
    } else if(id == p-1) {
    	lcsPopulateMatrixLast(id, matrix, seq1, substring);
    	lcsFindSubStringLastSlave(id, seq1, substring, matrix);
    } else {
    	lcsPopulateMatrixMiddle(id, matrix, seq1, substring);
    	lcsFindSubStringSlave(id, seq1, substring, matrix);
    }

    // for(int i = 0; i < p; i++) {
    // 	MPI_Barrier(MPI_COMM_WORLD);
    // 	if(i == id) {
    // 		lcsPrintMatrix(matrix); 		
    // 	}
    // }

	MPI_Finalize();

	#ifdef DEBUG_TIME
	double end = MPI_Wtime();
	if(id == 0) {
		std::cout << "time: " << end - start << std::endl;	
	}
	#endif
	
	return SUCCESS;
}