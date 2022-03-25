#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#define m 10                                                                // Define the sparse matrix rows
#define n 6                                                                 // Define the sparse matrix columns
#define k 4                                                                 // Define the fat vector columns

//Function to create the sparse matrix
double create_sparseMatrices(int mrow, int mcol) {
	int i, j;
	int Spmatrix[m][n];                                                    //For storing the sparse matrix, used inside the function                  

	//Creating a sparse matrix
	for (i = 0; i < mrow; i++) {
		for (j = 0; j < mcol; j++) {
			if (rand() % 5 == 0) {                                        //Providing the sparsity of sparse matrix
				Spmatrix[i][j] = rand() % 100;                            //Generating random values for non zero entries  
			}
			else {
				Spmatrix[i][j] = 0;                                       //Giving zeroes to the remaining values
			}
			return Spmatrix[i][j];
		}
	}

}

double create_fatVector(int mcol, int vcol) {
	int i, j;
	int Fvec[n][k];                                                      //For storing the fat vector, used inside the function

	// Creating a fat vector
	if (vcol >= (0.5 * mcol))                                           //Specifying the limit of 50% on the number of columns with regards to number of rows
		for (i = 0; i < mcol; i++) {
			for (j = 0; j < vcol; j++) {
				Fvec[i][j] = rand() % 5;                                //Generating random numbers
				return Fvec[i][j];
				printf("%2d   ", Fvec[i][j]);
			}
			printf("\n");
		}
	else
		printf("This does not make a fat-vector, Please enter optimum value!");//Generating error message for wrong number of columns entered
}

 
//**********************************************************************************************************************
//Program starts
int main(int argc, char** argv) {
	int myrank, npes;                                                    //unique ID given to each processor, number of processors

	MPI_Init(&argc, &argv);                                              //Communicators
	MPI_Comm_size(MPI_COMM_WORLD, &npes);                                //Communicators
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);                              //Communicators
	int i, j;                                                            //Used for iterating
	int mat[m][n], FatVector[n][k], Result[m][k];

	double comm_start, comm_stop, comm1_start, comm1_stop, start, stop;  //For estimating the time
	double parallel_time, computation_time, recv_computation_time;       //For estimating the time taken
	//**********************************************************************************************************************


	printf("The Sparse matrix is:\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			mat[i][j] = create_sparseMatrices(m, n);					//Calling the function for generating sparse matrix
			printf("%2d   ", mat[i][j]);
		}
		printf("\n");
	}
	printf("The fat vector is:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < k; j++) {
			FatVector[i][j] = create_fatVector(n, k);					//Calling the function for generating fat vector
			printf("%2d   ", FatVector[i][j]);
		}
		printf("\n");
	}

	//**********************************************************************************************************************
	MPI_Barrier(MPI_COMM_WORLD);                                       //For stopping all other tasks before computing time
	start = MPI_Wtime();                                               //starting the computation here

	if (myrank == 0) {
		for (i = 0; i < m; i++) {									   //Iterating through the number of rows of sparse matrix
			for (j = 0; j < k; j++) {								   //Iterating through the number of columns of fat vector
				Result[i][j] = 0;                                      //Assigned zero to the resultant that is a row
				for (int l = 0; l < n; l++) {						   //Iterating through the number of columns of sparse matrix
					Result[i][j] += mat[i][l] * FatVector[l][j];
				}
			}
		}
	}
	stop = MPI_Wtime();                                                //For stopping all the processes of computation
	computation_time = (stop - start) * 1000;                           //Getting the Computation time

	//**********************************************************************************************************************
	if (myrank == 0) {
		printf("The Result matrix is:\n");
		for (i = 0; i < m; i++) {
			for (j = 0; j < k; j++) {
				printf("%d   ", Result[i][j]);                         //Printing the results
			}
			printf("\n");
		}
		printf("\n Time of execution: %f\n", computation_time);        //Computation time
	}
}