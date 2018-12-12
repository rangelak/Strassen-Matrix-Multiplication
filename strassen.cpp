/************************************************************************************
********************** CS124: Programming assignment #2 *****************************
******************** Strassen Algorithm C++ Implementation **************************
*************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

// global variable for the dimension of the matrix
int dim;

// crossover point
int crossover;

/*******************************************************
********* Normal multiplication of two matrices. *******
*******************************************************/

void multiply (	std::vector< std::vector<int> > &A, 
				std::vector< std::vector<int> > &B, 
				std::vector< std::vector<int> > &C, int d)
{
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            for (int k = 0; k < d; k++)
            {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}


/*****************************************************
************** Functions for Strassen ****************
*****************************************************/

// function to split a matrix
void split (std::vector< std::vector<int> > &A, 
			std::vector< std::vector<int> > &B, 
			int row, int col, int d) 
{
    for (int i1 = 0, i2 = row; i1 < d; i1++, i2++)
    {
        for (int j1 = 0, j2 = col; j1 < d; j1++, j2++)
        {
            B[i1][j1] = A[i2][j2];
        }
    }
}

// join the 4 matrices in a resulting matrix
void join (	std::vector< std::vector<int> > &A, 
			std::vector< std::vector<int> > &B, 
			int row, int col, int d) 
{
    for (int i1 = 0, i2 = row; i1 < d; i1++, i2++)
    {
        for (int j1 = 0, j2 = col; j1 < d; j1++, j2++)
        {
        	B[i2][j2] = A[i1][j1];
        }
    }
}

// function to add two matrices. 
void add (	std::vector< std::vector<int> > &A, 
         	std::vector< std::vector<int> > &B, 
         	std::vector< std::vector<int> > &C, int d) 
{
    for (int i = 0; i < d; i++) 
    {
        for (int j = 0; j < d; j++) 
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

// function to subtract two matrices. 
void subtract (	std::vector< std::vector<int> > &A, 
         		std::vector< std::vector<int> > &B, 
         		std::vector< std::vector<int> > &C, int d) 
{
    for (int i = 0; i < d; i++) 
    {
        for (int j = 0; j < d; j++) 
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

// function that prints out the full matrix
void print_m (std::vector< std::vector<int> > m, int d)
{
	for (int i = 0; i < d; i++) 
	{
        for (int j = 0; j < d; j++) 
        {
             std::cout << "\t" << m[i][j];
        }
        std::cout << std::endl;
    }
}

// function to print diagonal of matrix
void print_diag (std::vector< std::vector<int> > m, int d)
{
	for (int i = 0; i < d; i++)
	{
		std::cout << m[i][i] << "\n";
	}
}

/********************************************************************
*********************** Strassen Algorithm **************************
********************************************************************/

void strassen (	std::vector< std::vector<int> > &A, 
              	std::vector< std::vector<int> > &B, 
              	std::vector< std::vector<int> > &C, int d)
{
	//if the dimension is not divisible by two, add padding
	if (d % 2 != 0)
	{
		d++;
		A.resize (d);
		B.resize (d);
		C.resize (d);
		
		// loop to resize the inside of the matrices
		for (int i = 0; i < d; i++)
		{
			A[i].resize (d);
			B[i].resize (d);
			C[i].resize (d);
		}
	}
	
	// if we are at the crossover point, call the normal multiplication
	if (d <= crossover)
	{
		multiply (A, B, C, d);
		return;
	}
	else
	{
		// setting the dimension of the new matrices
		int new_d = d/2;
        
        std::vector<int> inside (new_d);

        // initialize the submatices		
		std::vector< std::vector<int> > A11 (new_d, inside);
		std::vector< std::vector<int> > A12 (new_d, inside);
		std::vector< std::vector<int> > A21 (new_d, inside);
		std::vector< std::vector<int> > A22 (new_d, inside);
		std::vector< std::vector<int> > B11 (new_d, inside);
		std::vector< std::vector<int> > B12 (new_d, inside);
		std::vector< std::vector<int> > B21 (new_d, inside);
		std::vector< std::vector<int> > B22 (new_d, inside);
		std::vector< std::vector<int> > C11 (new_d, inside);
		std::vector< std::vector<int> > C12 (new_d, inside);
		std::vector< std::vector<int> > C21 (new_d, inside);
		std::vector< std::vector<int> > C22 (new_d, inside);

		// split matrices A and B in 4 submatrices
		split(A, A11, 0 , 0, new_d);
        split(A, A12, 0 , new_d, new_d);
        split(A, A21, new_d, 0, new_d);
        split(A, A22, new_d, new_d, new_d);
        split(B, B11, 0 , 0, new_d);
        split(B, B12, 0 , new_d, new_d);
        split(B, B21, new_d, 0, new_d);
        split(B, B22, new_d, new_d, new_d);

        // matrices to store results from arithmetic operations
		std::vector< std::vector<int> > result1 (new_d, inside);
		std::vector< std::vector<int> > result2 (new_d, inside);


        /************************************************
        *** 	M1 = (A11 + A22)(B11 + B22)
        *** 	M2 = (A21 + A22) B11
        *** 	M3 = A11 (B12 - B22)
        *** 	M4 = A22 (B21 - B11)
        *** 	M5 = (A11 + A12) B22
        *** 	M6 = (A21 - A11) (B11 + B12)
        *** 	M7 = (A12 - A22) (B21 + B22)
        ************************************************/

		// calculate M1
		add(A11, A22, result1, new_d);
		add(B11, B22, result2, new_d);
		std::vector< std::vector<int> > M1 (new_d, inside);
		strassen(result1, result2, M1, new_d);

		// calculate M2
		add(A21, A22, result1, new_d);
		std::vector< std::vector<int> > M2 (new_d, inside);
		strassen(result1, B11, M2, new_d);

		// calculate M3
		subtract(B12, B22, result2, new_d);
		std::vector< std::vector<int> > M3 (new_d, inside);
		strassen(A11, result2, M3, new_d);

		// calculate M4
		subtract(B21, B11, result2, new_d);
		std::vector< std::vector<int> > M4 (new_d, inside);
		strassen(A22, result2, M4, new_d);

		// calculate M5
		add(A11, A12, result1, new_d);
		std::vector< std::vector<int> > M5 (new_d, inside);
		strassen(result1, B22, M5, new_d);

		// calculate M6
		subtract(A21, A11, result1, new_d);
		add(B11, B12, result2, new_d);
		std::vector< std::vector<int> > M6 (new_d, inside);
		strassen(result1, result2, M6, new_d);

		// calculate M7
		subtract(A12, A22, result1, new_d);
		add(B21, B22, result2, new_d);
		std::vector< std::vector<int> > M7 (new_d, inside);
		strassen(result1, result2, M7, new_d);

		/**********************************************
        ***      C11 = M1 + M4 - M5 + M7
        ***      C12 = M3 + M5
        ***      C21 = M2 + M4
        ***      C22 = M1 - M2 + M3 + M6
        **********************************************/

		// calculating C11
       	add(M1, M4, result1, new_d);
       	add(result1, M7, result2, new_d);
       	subtract(result2, M5, C11, new_d);

       	// calculating C12
       	add(M3, M5, C12, new_d);

       	// calculating C21
       	add(M2, M4, C21, new_d);

       	// calculating C22
       	subtract(M1, M2, result1, new_d);
       	add(M3, M6, result2, new_d);
       	add(result1, result2, C22, new_d);

       	// add the resulting matrices in one matrix
       	join(C11, C, 0 , 0, new_d);
        join(C12, C, 0 , new_d, new_d);
        join(C21, C, new_d, 0, new_d);
        join(C22, C, new_d, new_d, new_d);
	}
}

/********************************************************************
************************* Main Function *****************************
********************************************************************/

int main (int argc, char* argv[])
{
	// ensure proper usage
	if (argc != 4)
	{
		std::cout << "Usage: ./strassen 0 dimension inputfile\n";
		return 1;
	}

	// ensure no negative dimensions of the array
	dim = atoi(argv[2]);
	if (dim < 1)
	{
		std::cout << "Dimension has to be greater than 0.\n";
		return 2;
	}

	// set the crossover
	crossover = 86;

	// vector for the inside of the matrices
	std::vector<int> inside (dim);

	// initialize matrices to be multiplied
	std::vector< std::vector<int> > A (dim, inside);
	std::vector< std::vector<int> > B (dim, inside);
	std::vector< std::vector<int> > C (dim, inside);


	// intialize how many lines we will be scanning
	int size_file = dim*dim*2;

	// open and read from the file
	std::ifstream infile(argv[3]);

	// check if file opened
	if (infile.is_open())
	{
		// each line in the file is a string
		std::string line;

		// iterate over lines of file and add them to matrices A,B
		for (int i = 0, j = 0; i < size_file; i++)
		{
			getline (infile, line);
			if (i < (dim*dim))
			{
				A[i/dim][i%dim] = std::stoi(line);
			}
			else 
			{
				B[j/dim][i%dim] = std::stoi(line);
				j++;
			}
		}
	}
	else
	{
		perror("Something went wrong; couldn't open file.\n");
		return 3;
	}

	// namespace for timing function
	using namespace std::chrono;
	
	// time functions
	/*    high_resolution_clock::time_point start1 = high_resolution_clock::now();	   
	**	  multiply(A, B, C, dim);
	**    high_resolution_clock::time_point end1 = high_resolution_clock::now();
	**    auto duration1 = duration_cast<microseconds>(end1 - start1).count();
	**    std::cout << "This took: "<< duration1 << " microseconds.\n";
	**/

	/*    high_resolution_clock::time_point start2 = high_resolution_clock::now();	   
	**	  strassen(A, B, C, dim);
	**	  high_resolution_clock::time_point end2 = high_resolution_clock::now();
	**    auto duration2 = duration_cast<microseconds>(end2 - start2).count();
	**    std::cout<< "This took: " << duration2 << " microseconds.\n";
	**/

	// run strassen
	strassen(A, B, C, dim);

    // print the matrix
    print_diag (C, dim);

   	return 0;
}

