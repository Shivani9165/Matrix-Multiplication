#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#define m 10                                                                // Define the sparse matrix rows
#define n 6                                                                 // Define the sparse matrix columns
#define k 4                                                                 // Define the fat vector columns

double create_sparseMatrices(int, int);
double create_fatVector(int, int);

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


	int mat[m][n], sparseMat[n][m],                                      //Sparse matrix of size (m x n), Recieving buffer for SPARSE MATRIX MPI_Scatterv size(m x n)
		trans[n][m];                                                     //transpose of matrix of size (n x m)
	int FatVector[n][k], fatVectorTemp[n][k], Result[m][k], Multiply[m][k];//Fat vector of size (n x k), Recieving buffer for FAT VECTOR MPI_Scatterv size(m xn)
	int i, j;                                                            //Used for iterating 
	int displs[npes],                                                    //Used for iterating
		vecDispls[npes], 
		sendCount[npes],                                                //Used for specifying number of rows to send in for MPI_Scatterv
		load_array[npes],                                               //Used for load balancing storing number of rows for each processor
		rowCount[npes];                                                 //For number of rows to fat vector in MPI_Scatterv
	double comm_start, comm_stop, comm1_start, comm1_stop, start, stop; //For estimating the time
	double parallel_time, computation_time, recv_computation_time;      //For estimating the time taken
	int rem = 0;                                                        //used for storing the remainder from load balancing
	int sum = 0; 
	int vecSum = 0;                                                     // used for storing the remainder from load balancing of vector

	
	//***************************************************************************************************************
	if (myrank == 0) {                                                  //Running for 1 processor
		printf("The Sparse matrix is:\n");                              //Printing sparse matrix
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				mat[i][j] = create_sparseMatrices(m, n);               //Calling the function for generating sparse matrix
				printf("%2d   ", mat[i][j]);
			}
			printf("\n");
		}
		printf("The fat vector is:\n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < k; j++) {
				FatVector[i][j] = create_fatVector(n, k);              //Calling the function for generating fat vector
				printf("%2d   ", FatVector[i][j]);                     //Printing fat vector
			}
			printf("\n");
		}
	}
	//**************************************************************************************************************
	//Declare and initialise the full array

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			trans[j][i] = mat[i][j];                                  //Generating transpose of matrix
		}
	}
	//*******************************************************************************************************************
	MPI_Barrier(MPI_COMM_WORLD);                                      //For stopping all other tasks before communication time
	comm_start = MPI_Wtime();                                         //Communication starts here

	//for Sparse matrix load-balancing (size = m * n)
	rem = n % npes;                                                   //intializing the remainder
	for (j = 0; j < npes; j++) {
		load_array[j] = n / npes;                                     //Number of columns to reach each processor
		if (rem > 0) {                                                //If there are any residual columns present to balance the load
			load_array[j] ++;                                         //Residual columns are added to the balanced load
			rem--;                                                    //count reduced as the load is getting balanced
		}
		sendCount[j] = load_array[j] * m;                             //Columns to send or elements to senf for MPI_Scatterv
		displs[j] = sum;                                              //Displacement for MPI_Scatterv
		sum += sendCount[j];
	}


	// for fat vector load-balancing (size = n * k)
	for (j = 0; j < npes; j++) {                                     
		rowCount[j] = load_array[j] * k;                             //load balncing for fat vector
		vecDispls[j] = vecSum;                                       //for givin the displacement needed for fat vector
		vecSum += rowCount[j];
	}

	//Sparse matrix scattered on recieving buffer 'sparseMat'
	MPI_Scatterv(trans, sendCount, displs, MPI_INT, sparseMat, sendCount[myrank], MPI_INT, 0, MPI_COMM_WORLD);

	//Fat vector scattered on recieving buffer 'fatVectorTemp'
	MPI_Scatterv(FatVector, rowCount, vecDispls, MPI_INT, fatVectorTemp, rowCount[myrank], MPI_INT, 0, MPI_COMM_WORLD);
	comm_stop = MPI_Wtime();                                       //First communication stops here
	//********************************************************************************************************************
	//Starting matrix multiplication here:

	MPI_Barrier(MPI_COMM_WORLD);                                     //First communication starts here
	start = MPI_Wtime();                                             //starting the computation here

	if (myrank < npes) {                                             //If the rank is less than the number of processors given
		for (i = 0; i < load_array[myrank]; i++) {                   //number of rows reached to each processor
			for (j = 0; j < k; j++) {                                //Iterating through the number of columns of fat vector
				for (int l = 0; l < m; l++) {                        //Iterating through the number of rows of sparse matrix
					if (i == 0)
						Result[l][j] = sparseMat[i][l] * fatVectorTemp[i][j];
					else
						Result[l][j] += sparseMat[i][l] * fatVectorTemp[i][j];
				}
			}
		}
	}
	stop = MPI_Wtime();                                              //For stopping all the processes of computation
	//**************************************************************************************************************************
	//Communicating again
	MPI_Barrier(MPI_COMM_WORLD);                                     //For stopping all other tasks before communication time
	comm1_start = MPI_Wtime();                                       //starting second communication here

	//Reducing the results obtained by Result and summing them by MPI_SUM, on resultant Multiply
	MPI_Reduce(Result, Multiply, m * k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	comm1_stop = MPI_Wtime();                                        //Stopping the second communication
	parallel_time = ((comm_stop - comm_start) + (comm1_stop - comm1_start)) * 1000; //Computing the communication time
	computation_time = (stop - start) * 1000;                                       //Getting the Computation time

	//Getting the maximum computation time
	MPI_Reduce(&computation_time, &recv_computation_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	//******************************************************************************************************************************
	if (myrank == 0) {
		printf("The Result matrix is:\n");
		for (i = 0; i < m; i++) {
			for (j = 0; j < k; j++) {
				printf("%d  ", Multiply[i][j]);                       //Printing the results
			}
			printf("\n");
		}
		printf("\n Time of computation: %f\n", recv_computation_time);  //Computation time
	}
	//******************************************************************************************************************************
	
	MPI_Finalize();                                                   //Finalizing the program
}