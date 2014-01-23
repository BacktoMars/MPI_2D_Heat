/* C Example */
// sequential version
// boundaries are set to 0

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>


//const int start = 0;
//const int end = 16;

#define NXPROB      8                /* x dimension of problem grid */
#define NYPROB      8               /* y dimension of problem grid */
#define STEPS       1                /* number of time steps */
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



//void inidat(), prtdat(), update();

void inidat(int nx, int ny, float *u);
void update(int start, int end, int ny, float *u1, float *u2);
void prtdat(int nx, int ny, float *u1, char *fnam);

int main (int argc, char* argv[])
{
	int rank, size;

	printf("First the serial version!\n");
	//	float unew[NXPROB][NYPROB];
	//	float uold[NXPROB][NYPROB];

	float *unew = (float *)malloc(sizeof(float)*NXPROB*NYPROB);
	float *uold = (float *)malloc(sizeof(float)*NXPROB*NYPROB);

	inidat(NXPROB, NYPROB, uold); // initialize the old temperature map
	inidat(NXPROB, NYPROB, unew);
	prtdat(NXPROB, NYPROB, uold, "t0");	

	getchar();


	float *temp;

	int iternum = 0;
	for (iternum =0; iternum < STEPS; iternum++)
	{
		//printf("unew[0][1] = %f\n", *(unew + 1));
		//printf("uold[0][0] = %f\n", *(uold));
		update(1, NXPROB-2, NYPROB, uold, unew);


		int i=0;
		//prtdat(NXPROB, NYPROB, unew, stdout);	
		//for (i=0; i<20; i++)
		//{		
		//	printf("unew[%d][%d] = %f\n", i, i, *(unew+i*NYPROB+i));
		//	printf("uold[%d][%d] = %f\n", i, i, *(uold+i*NYPROB+i));
		//
		//}

		//printf("interchage two pointers.\n");
		// interchage uold and unew
		temp = uold;
		uold = unew;
		unew = temp;		

		//for (i=0; i<20; i++)
		//{		
		//	printf("unew[%d][%d] = %f\n", i, i, *(unew+i*NYPROB+i));
		//	printf("uold[%d][%d] = %f\n", i, i, *(uold+i*NYPROB+i));
		//
		//}
	}

	//char* outFileName = new char[100];
	char outFileName[20];
	sprintf(outFileName, "t%d", STEPS);
	prtdat(NXPROB, NYPROB, uold, outFileName);
	getchar();



//	MPI_Init (&argc, &argv);      /* starts MPI */
//	MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
//	MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */
//	printf( "Hello world from process %d of %d\n", rank, size );
//	MPI_Finalize();
	return 0;
}

/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update(int start, int end, int ny, float *u1, float *u2)
{
	int ix, iy;
	for (ix = start; ix <= end; ix++) 
		for (iy = 1; iy <= ny-2; iy++) 
{	printf("ix = %d iy = %d\n", ix, iy);	
		*(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
				parms.cx * (*(u1+(ix+1)*ny+iy) +
						*(u1+(ix-1)*ny+iy) - 
						2.0 * *(u1+ix*ny+iy)) +
				parms.cy * (*(u1+ix*ny+iy+1) +
						*(u1+ix*ny+iy-1) - 
						2.0 * *(u1+ix*ny+iy));
}
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u) {
	int ix, iy;

	for (ix = 0; ix <= nx-1; ix++) 
		for (iy = 0; iy <= ny-1; iy++)
			*(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}

/**************************************************************************
 * subroutine prtdat
 **************************************************************************/
void prtdat(int nx, int ny, float *u1, char *fnam) {
	int ix, iy;
	FILE *fp;

	fp = fopen(fnam, "w");
	for (iy = ny-1; iy >= 0; iy--) {
		for (ix = 0; ix <= nx-1; ix++) {
			fprintf(fp, "%8.1f", *(u1+ix*ny+iy));
			if (ix != nx-1) 
				fprintf(fp, " ");
			else
				fprintf(fp, "\n");
		}
	}
	fclose(fp);
}
