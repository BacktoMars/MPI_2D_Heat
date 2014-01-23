/* C Example */
#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

// use multiple MPI processes to do a 2D heat transfer simulation
// each process is used to deal with multiple columns.

//
//	-------------------------------------
//	|									|	
//	|									|
//	|NYPROB								|
//	|									|
//	|									|
//	-------------------------------------
//
//				NXPROB
//



//const int start = 0;
//const int end = 16;

#define NXPROB      8                 /* x dimension of problem grid */
#define NYPROB      6                 /* y dimension of problem grid */
#define STEPS       6                /* number of time steps */
#define MAXWORKER   8                  /* maximum number of worker tasks */
#define MINWORKER   3                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */

struct Parms { 
	float cx;
	float cy;
} parms = {0.1, 0.1};


void inidat(int nx, int ny, float *u);
void update(int start, int end, int nx, int ny, float *u1, float *u2);
void prtdat(int nx, int ny, float *u1, char *fnam);


int main (int argc, char* argv[])
{
	int rank, size;

	float *unew = (float *)malloc(sizeof(float)*NXPROB*NYPROB);
	float *uold = (float *)malloc(sizeof(float)*NXPROB*NYPROB);

	inidat(NXPROB, NYPROB, uold); // initialize the old temperature map
	inidat(NXPROB, NYPROB, unew);
	prtdat(NXPROB, NYPROB, uold, "par-t0");	

	int numWorkers;
	int colPerProcess;
	int left, right;
	int colOffset;


	MPI_Init (&argc, &argv);      /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */
	MPI_Status status;


	// define the column vector data type
	MPI_Datatype col_vector;
	MPI_Type_vector(NYPROB, 1, NXPROB, MPI_FLOAT, &col_vector);
	MPI_Type_commit(&col_vector);

	//	printf( "Hello world from process %d of %d\n", rank, size );

	// use rank 0 as the master process
	numWorkers = size - 1;
	colPerProcess = NXPROB/numWorkers;

	if (rank==0) // master node
	{
		printf("Distribute tasks to %d MPI processes\n", numWorkers);

		int j=0;

		int colOffset = 0;
		for (j=1; j<=numWorkers; j++)
		{
			if (j == 1)
				left = NONE;
			else
				left = j - 1;

			if (j == numWorkers)
				right = NONE;
			else
				right = j + 1;

			MPI_Send(&colOffset, 1, MPI_INT, j, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&colPerProcess, 1, MPI_INT, j, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&left, 1, MPI_INT, j, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&right, 1, MPI_INT, j, BEGIN, MPI_COMM_WORLD);

			int k=0;
			// send the column data to each MPI process
			for (k=0; k<NYPROB; k++)
				MPI_Send(uold+k*NXPROB+colOffset, colPerProcess, MPI_FLOAT, j, BEGIN, MPI_COMM_WORLD);

			colOffset = colOffset + colPerProcess;

		}
		// wait for results from workers
		for (j=1; j<=numWorkers; j++)
		{
			MPI_Recv(&colOffset,1, MPI_INT, j, DONE, MPI_COMM_WORLD, &status);
			MPI_Recv(&colPerProcess,1, MPI_INT, j, DONE, MPI_COMM_WORLD, &status);
			int k=0;
			for (k=0; k<NYPROB; k++){	
				MPI_Recv(uold+k*NXPROB+colOffset,colPerProcess, MPI_FLOAT, j, DONE, MPI_COMM_WORLD, &status);
				printf("received %d data from process %d\n", colPerProcess, j);
			}
		}
		char fileName[20];
		sprintf(fileName, "par-t%d", STEPS);
		prtdat(NXPROB, NYPROB, uold, fileName);
	}
	// worker processes
	if (rank != 0)
	{
		MPI_Recv(&colOffset, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
		MPI_Recv(&colPerProcess, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
		MPI_Recv(&left, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
		MPI_Recv(&right, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);

		int k=0;
		// get the column data from master process
		for (k=0; k<NYPROB; k++)
			MPI_Recv(uold+k*NXPROB+colOffset, colPerProcess, MPI_FLOAT, MASTER, BEGIN, MPI_COMM_WORLD, &status);

		int colStart = colOffset;
		int colEnd = colOffset + colPerProcess;

		printf("task=%d  start=%d  end=%d\n",rank, colStart, colEnd);	
		printf("left = %d right = %d \n", left, right);	
		int i, j;
		for (i=0; i<NYPROB; ++i)
		{	
			for (j=0; j<colPerProcess; ++j)
			{
				printf ( " %f ", uold[i*NXPROB + j + colStart] );
			}
			printf("\n");
		}


		int start, end;
		// get the start and end positions for each process
		if (left == NONE)	// leftmost block
			start = 1;
		else
			start = colStart;

		if (right == NONE) // rightmost block
			end = colEnd - 2;
		else
			end = colEnd - 1;

		int iter = 0;
		for (iter = 0; iter<STEPS; ++iter)
		{
			// interchange data for neighbour processes
			if (left != NONE)
			{
				MPI_Send(uold + colStart, 1, col_vector, left, 100, MPI_COMM_WORLD);
				MPI_Recv(uold + colStart - 1, 1, col_vector, left, 101, MPI_COMM_WORLD, &status);
			}
			if (right != NONE)
			{
				MPI_Send(uold + colEnd - 1, 1, col_vector, right, 101, MPI_COMM_WORLD );
				MPI_Recv(uold + colEnd, 1, col_vector, right, 100, MPI_COMM_WORLD, &status);
			}

			// update the local matrix
			int ix, iy;

			for (ix = start; ix <= end; ix++) 
			{	
				//printf("ix = %d\n", ix);	

				for (iy = 1; iy <= NYPROB-2; iy++)
				{ 
					unew[iy*NXPROB+ix] = uold[iy*NXPROB+ix]  + 
						parms.cx * (uold[(iy+1)*NXPROB+ix] +
								uold[(iy-1)*NXPROB+ix] - 
								2.0 * uold[iy*NXPROB+ix]) +
						parms.cy * ( uold[iy*NXPROB+ix+1] +
								uold[iy*NXPROB+ix-1] - 
								2.0 * uold[iy*NXPROB+ix]);
				}
			}

			// interchange uold and unew

			for (ix = start; ix <= end; ++ix)
				for (iy = 1; iy <= NYPROB - 2; ++iy)
				{	uold[iy*NXPROB + ix] = unew[iy*NXPROB + ix];
				}


			// send results back to master process
			/*
			   printf("updated u in process %d:\n", rank);
			   for (iy = 1; iy <= NYPROB-2; iy++)
			   {
			   for (ix = start; ix <= end; ix++) 
			   { 
			   printf("%8.1f", uold[iy*NXPROB+ix]);
			   if (ix != end) 
			   printf(" ");
			   else
			   printf("\n");
			   }	
			   }
			 */
			//MPI_Finalize();
		} // end iteration

		MPI_Send(&colOffset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
		MPI_Send(&colPerProcess, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
		for (k=0; k<NYPROB; k++)
			MPI_Send(uold+k*NXPROB+colOffset, colPerProcess, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);


	} // end rank!=0


	MPI_Finalize();

}

/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update(int start, int end, int nx, int ny, float *u1, float *u2)
{
	int ix, iy;
	for (ix = start; ix <= end; ix++) 
		for (iy = 1; iy <= ny-2; iy++) 
		{	//printf("ix = %d iy = %d\n", ix, iy);	
			u2[iy*nx+ix] = u1[iy*nx+ix]  + 
				parms.cx * ( u1[(iy+1)*nx+ix] +
						u1[(iy-1)*nx+ix] - 
						2.0 * u1[iy*nx+ix] ) +
				parms.cy * ( u1[iy*nx+ix+1] +
						u1[iy*nx+ix-1] - 
						2.0 * u1[iy*nx+ix]);
		}
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u) {
	/*	int ix, iy;

		for (ix = 0; ix <= nx-1; ix++) 
		for (iy = 0; iy <= ny-1; iy++)
	 *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
	 */
	int i,j;
	for (i=0; i<ny; ++i)
		for (j=0; j<nx; ++j)
		{
			u[i*nx+j] = (float)(i*(ny - i -1)*j*(nx - j -1));
		}
}

/**************************************************************************
 * subroutine prtdat
 **************************************************************************/
void prtdat(int nx, int ny, float *u1, char *fnam) {
	int i, j;
	FILE *fp;

	fp = fopen(fnam, "w");
	//	for (iy = ny-1; iy >= 0; iy--) {
	//		for (ix = 0; ix <= nx-1; ix++) {
	for (i=0; i < ny; i++) {
		for(j=0; j < nx; j++) {

			fprintf(fp, "%8.1f", u1[i*nx+j]);
			if (j != nx-1) 
				fprintf(fp, " ");
			else
				fprintf(fp, "\n");
		}
	}
	fclose(fp);
}
