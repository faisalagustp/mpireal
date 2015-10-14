/*
**********************************************************************

                   Example 13 (vv_mult_blkstp_unf.c)

   Objective           : Vector_Vector Multiplication (Using Block Partition)

   Input               : Read files (vdata1.inp) for Vector_A and
                         (vdata2.inp) for Vector_B

   Output              : Result of Vector-Vector Multiplication on processor 0

   Necessary Condition : Number of Processes should be less than
                         or equal to 8. Vector size for Vectors A and
                         B should be properly striped. that is Vector
                         size should be properly divisible by
                         Number of processes used.

**********************************************************************
*/



#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define Epsilon 1.0E-10
int main(int argc, char** argv)
{

  int	   Numprocs, MyRank;
  int 	VectorSize, VectorSize_A, ScatterSize, VectorSize_B;
  int	   index, iproc;
  int	   Root = 0, ValidInput = 1, ValidOutput = 1;
  float  *Mybuffer_A, *Mybuffer_B,  MyFinalVector, FinalAnswer;
  float  CheckResultVector;
  FILE	*fp;
  int	   VectorA_FileStatus = 1, VectorB_FileStatus = 1;
  float  *Vector_A, *Vector_B;
  double start, end, startscatter, endscatter, startcomm, endcomm, commtime, totalcomm;


/* ........MPI Initialisation .......*/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &Numprocs);


  if(MyRank == Root) {

     //input of vector A
     VectorSize_A = 1024*1024*1024;
     Vector_A = (float *) malloc(VectorSize_A*sizeof(float));

     for(index=0;index < VectorSize_A; index++){
       Vector_A[index] = 1;
     }

     //input of vector A
     VectorSize_B = VectorSize_A;
     Vector_B = (float *) malloc(VectorSize_B*sizeof(float));

     for(index=0;index < VectorSize_B; index++){
       Vector_B[index] = 1;
     }



	  if(VectorSize_A != VectorSize_B)
			ValidInput = 0;
    start = MPI_Wtime();
	}

   MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(&VectorA_FileStatus, 1, MPI_INT, Root, MPI_COMM_WORLD);
   if(VectorA_FileStatus == 0) {
       if(MyRank == Root) printf("Can't open input file for Vector A\n");
	  MPI_Finalize();
		  exit(-1);
   }

	MPI_Bcast(&VectorB_FileStatus, 1, MPI_INT, Root, MPI_COMM_WORLD);
   if(VectorB_FileStatus == 0) {
     	  if(MyRank == Root) printf("Can't open input file for Vector B\n");
		  MPI_Finalize();
		  exit(-1);
	}

   MPI_Bcast(&ValidInput, 1, MPI_INT, Root, MPI_COMM_WORLD);
   if(ValidInput == 0) {
      if(MyRank == Root) printf("Size of both vectors is not equal ....\n");
		MPI_Finalize();
		exit(-1);
   }

   VectorSize = VectorSize_A;
   MPI_Bcast(&VectorSize, 1, MPI_INT, Root, MPI_COMM_WORLD);

   if(VectorSize < Numprocs) {
	 MPI_Finalize();
	if(MyRank == 0)
	 printf("VectorSize should be more than No of Processors ..... \n");
	 exit(0);
   }

   if(VectorSize % Numprocs != 0) {
	   MPI_Finalize();
	   if(MyRank == 0) {
	      printf("VectorSize Can not be Striped Evenly ..... \n");
	   }
	   exit(0);
   }

/* Scatter vector A and B */

   ScatterSize = VectorSize / Numprocs;
   Mybuffer_A = (float *)malloc(ScatterSize * sizeof(float));
   startscatter = MPI_Wtime();
   MPI_Scatter( Vector_A, ScatterSize, MPI_FLOAT,
					 Mybuffer_A, ScatterSize, MPI_FLOAT,
					 Root, MPI_COMM_WORLD);
   endscatter = MPI_Wtime();
   Mybuffer_B = (float *)malloc(ScatterSize * sizeof(float));
   MPI_Scatter( Vector_B, ScatterSize, MPI_FLOAT,Mybuffer_B,
					 ScatterSize, MPI_FLOAT,
					 Root, MPI_COMM_WORLD);


   commtime = (endscatter - startscatter);
   MPI_Reduce(&commtime, &totalcomm, 1, MPI_DOUBLE, MPI_MAX, Root, MPI_COMM_WORLD);

/* Calculate partial sum */

	MyFinalVector = 0.0;
	for(index = 0 ; index < ScatterSize ; index++)
		 MyFinalVector += (Mybuffer_A[index] * Mybuffer_B[index]);

/* Collective computation : Final answer on process 0 */
   MPI_Reduce(&MyFinalVector, &FinalAnswer, 1, MPI_FLOAT, MPI_SUM, Root, MPI_COMM_WORLD);



	if(MyRank == 0){
	   CheckResultVector = 0.0;

		for(index = 0 ; index < VectorSize ; index++)
	 	 CheckResultVector += (Vector_A[index] * Vector_B[index]);

		if(fabs((double)(FinalAnswer-CheckResultVector)) > Epsilon){
		   printf("Error %d\n",index);
			ValidOutput = 0;
		}
    end = MPI_Wtime();

    printf("Proc Size : %d\n",Numprocs);
    printf("Vector Size Size : %d\n",VectorSize_A);
    printf("Total Time : %f\n",end-start);
    printf("Communication Time : %f\n",totalcomm);
		/*
    if(ValidOutput)
		  printf("FinalResult = %f\n", FinalAnswer);
        else
		  printf("Result may be wrong\n");
    */
	   free(Vector_A);
	   free(Vector_B);
   }


   free(Mybuffer_A);
   free(Mybuffer_B);
   MPI_Finalize();
}
