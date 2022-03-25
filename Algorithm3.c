#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
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

//***********************************************************************************************************************
//Program starts
int main(int argc, char** argv) {
	int myrank, npes;                                                    //unique ID given to each processor, number of processors

	MPI_Init(&argc, &argv);                                              //Communicators
	MPI_Comm_size(MPI_COMM_WORLD, &npes);                                //Communicators
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);                              //Communicators

	int mat[m][n], sparseMat[m][n];                                      //Sparse matrix of size (m x n), Recieving buffer for SPARSE MATRIX MPI_Scatterv size(m x n)
	int FatVector[n][k], Multiply[m][k];                                 // Fat vector of size(n x k), Resulting vector(m x k)
	int i, j, rem, size;                                                 //Used for iterating, used for storing the remainder from load balancing
	int load_array[npes],                                                //Used for load balancing storing number of rows for each processor
		rowToSend[npes],                                                 //Used for specifying number of rows to send for MPI_Scatterv
		displs[npes],                                                    //For providing the displacement needed during MPI_Scatterv
		counts[npes],                                                    //For number of rows to send of resultant in MPI_Gatherv
		vecDispls[npes];                                                 //For providing the diplacement needed during MPI_Gatherv
	double comm_start, comm_stop, comm1_start, comm1_stop, start, stop;  //For estimating the time
	double parallel_time, computation_time, recv_computation_time;       //For estimating the time taken
	int sum = 0;
	int vecSum = 0;
	size = m * n;
	int val[size], colIdx[size], rowIdx[size];                           //Indexes created for stroing value, row and column index
	
	//*************************************************************************************************************************

	if (myrank == 0) {                                                   //Running for 1 processor
		printf("The Sparse matrix is:\n");
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				mat[i][j] = create_sparseMatrices(m, n);                 //Calling the function for generating sparse matrix
				printf("%2d   ", mat[i][j]);                             //Printing sparse matrix
			}
			printf("\n");
		}
		printf("The fat vector is:\n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < k; j++) {
				FatVector[i][j] = create_fatVector(n, k);                //Calling the function for generating fat vector
				printf("%2d   ", FatVector[i][j]);                       //Printing fat vector
			}
			printf("\n");
		}
	}

	//*************************************************************************************************************************
	MPI_Barrier(MPI_COMM_WORLD);										 //For stopping all other tasks before communication time
	comm_start = MPI_Wtime();											 //Communication starts here
	MPI_Bcast(FatVector, n * k, MPI_INT, 0, MPI_COMM_WORLD);			 //Broadcasting the Fat vector on all processes

	rem = m % npes;                                                      //Getting the remaining rows after even distribution over processors
	for (j = 0; j < npes; j++) {
		load_array[j] = m / npes;                                        //Number of rows to rach each processor
		if (rem > 0) {                                                   //If there are any residual rows present to balance the load
			load_array[j] ++;                                            //Residual rows are added to the balanced load
			rem--;                                                       //count reduced as the load is getting balanced
		}
		rowToSend[j] = load_array[j] * n;                                //Rows to send or elements to senf for MPI_Scatterv
		displs[j] = sum;                                                 //Displacement for MPI_Scatterv
		sum += rowToSend[j];
	}

	int Result[load_array[myrank]][k];                                  //Initialing the Result matrix giving the size
	for (int r = 0; r < load_array[myrank]; ++r) {                     
		for (int c = 0; c < k; ++c) {
			Result[r][c] = 0;                                          //Initializing the result matrix by zeroes
		}
	}

	//Sparse matrix scattered on recieving buffer 'sparseMat'
	MPI_Scatterv(mat, rowToSend, displs, MPI_INT, sparseMat, rowToSend[myrank], MPI_INT, 0, MPI_COMM_WORLD);

	comm_stop = MPI_Wtime();                                            //First communication stops here

	//***********************************************************************************************************************

	MPI_Barrier(MPI_COMM_WORLD);                                     //For stopping all other tasks before computing time
	start = MPI_Wtime();                                             //starting the computation here

	int count = 0;                                                   //Initializing count as 0
	for (i = 0; i < load_array[myrank]; i++) {                       //number of rows reached to each processor
		for (j = 0; j < n; j++) {                                    //Iterating through the number of columns of sparse matrix
			if (sparseMat[i][j] != 0) {                              //Checking the non zeroes
				val[count] = sparseMat[i][j];                        //Creating the value index
				colIdx[count] = j;                                   //Creating the column index
				rowIdx[count] = i;                                   //Creating the row index
				count += 1;          
			}
		}
	}
	//***********************************************************************************************************************
	//Multiplication starting
	
	for (i = 0; i < load_array[myrank]; i++) {
		for (j = 0; j < count; j++) {
			if (rowIdx[j] == i) {                                     //Matching the row index
				for (int l = 0; l < k; l++) {                         //Iterating through the number of columns of fat vector
					Result[i][l] += val[j] * FatVector[colIdx[j]][l];
				}
			}
		}
	}

	stop = MPI_Wtime();                                               //For stopping all the processes of computation
//------------------------------------------------------------------------------------------------------------
	MPI_Barrier(MPI_COMM_WORLD);                                     //For stopping all other tasks before communication time
	comm1_start = MPI_Wtime();                                       //starting second communication here


	for (j = 0; j < npes; j++) {                                     //load balncing for MPI_Gatherv
		counts[j] = load_array[j] * k;                               //for giving the count of rows to MPI_Gatherv
		vecDispls[j] = vecSum;                                       //for displacemnt of MPI_Gtaherv
		vecSum += counts[j];
	}


	//Gathering all the resulatants from Result to combine into one matrix
	MPI_Gatherv(Result, counts[myrank], MPI_INT, Multiply, counts, vecDispls, MPI_INT, 0, MPI_COMM_WORLD);
	
	comm1_stop = MPI_Wtime();                                        //Stopping the second communication

	parallel_time = ((comm_stop - comm_start) + (comm1_stop - comm1_start)) * 1000; //Computing the communication time
	computation_time = (stop - start) * 1000;                                       //Getting the Computation time

	//Getting the maximum computation time
	MPI_Reduce(&computation_time, &recv_computation_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	//***********************************************************************************************************************

	if (myrank == 0) {
		printf("The Result matrix is:\n");                           //Printing the results
		for (i = 0; i < m; i++) {
			for (j = 0; j < k; j++) {
				printf("%2d  ", Multiply[i][j]);                     //Final result matrix
			}
			printf("\n");
		}
		printf("\n Time of computation: %f\n", recv_computation_time);  //Computation time
	}
	//***********************************************************************************************************************

	MPI_Finalize();                                                  //Finalizing the program

}
