
/*
*********************************************************************

                  Example 18 (mv_mult_blkstp.c)

   Objective           : Matrix_Vector multiplication
                        (Using block striped partitioning)

   Input               : Process 0 reads files (mdata.inp) for Matrix
                         and (vdata.inp) for Vector

   Output              : Process 0 prints the result of Matrix_Vector
                         Multiplication

   Necessary Condition : Number of rows of matrix should be less than
                         number of processors and perfectly divisible
                         by number of processors.Number of processors
                         should be less than or equal to 8

*********************************************************************
*/



#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv)
{

  int	   Numprocs, MyRank;
  int 	   NoofCols, NoofRows, VectorSize, ScatterSize;
  int	   index, irow, icol, iproc;
  int	   Root = 0, ValidOutput = 1;
  float    **Matrix, *Buffer, *Mybuffer, *Vector,
	   *MyFinalVector, *FinalVector;
  float    *CheckResultVector;
  FILE	   *fp;
  int	   MatrixFileStatus = 1, VectorFileStatus = 1;
  double    start, end, startcomm, endcomm, commtime, totalcomm;


  /* ........MPI Initialisation .......*/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &Numprocs);


  if(MyRank == 0)
  {
    //creating matrix
    start = MPI_Wtime();
    VectorSize = 1024*32;
    NoofCols = VectorSize;
    NoofRows = VectorSize;

    Matrix = (float **)malloc(NoofRows*sizeof(float *));
    for(irow=0 ;irow<NoofRows; irow++){
      Matrix[irow] = (float *)malloc(NoofCols*sizeof(float));
      for(icol=0; icol<NoofCols; icol++) {
 	       Matrix[irow][icol] = 1;
      }
    }

    Buffer = (float *)malloc(NoofRows*NoofCols*sizeof(float));
    index = 0;
    for(irow=0; irow< NoofRows; irow++){
      for(icol=0; icol< NoofCols; icol++) {
        Buffer[index] = Matrix[irow][icol];
        index++;
      }
    }

    //creating vector

    Vector = (float*)malloc(VectorSize*sizeof(float));
    for(index = 0; index<VectorSize; index++)
         Vector[index]=1;

   }  /* end  of if myrank = 0 */


   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Bcast (&MatrixFileStatus, 1, MPI_INT, Root, MPI_COMM_WORLD);
   if(MatrixFileStatus == 0) {
      if(MyRank == Root)
	 printf("Can't open input file for Matrix ..... \n");
      MPI_Finalize();
      exit(-1);
   }

   MPI_Bcast (&VectorFileStatus, 1, MPI_INT, Root, MPI_COMM_WORLD);
   if(VectorFileStatus == 0) {
      if(MyRank == Root)
	 printf("Can't open input file for Vector ..... \n");
      MPI_Finalize();
      exit(-1);
   }

   MPI_Bcast(&NoofRows, 1, MPI_INT, Root, MPI_COMM_WORLD);

   if(NoofRows < Numprocs) {
      if(MyRank == 0)
         printf("No of Rows should be more than No of Processors ... \n");
      MPI_Finalize();
      exit(0);
   }

   if(NoofRows % Numprocs != 0) {
      if(MyRank == 0)
	 printf("Matrix Can not be Striped Evenly ..... \n");
      MPI_Finalize();
      exit(0);
   }

   MPI_Bcast(&NoofCols, 1, MPI_INT, Root, MPI_COMM_WORLD);
   MPI_Bcast(&VectorSize, 1, MPI_INT, Root, MPI_COMM_WORLD);

   if(VectorSize != NoofCols){
      if(MyRank == 0){
         printf("Invalid input data..... \n");
	 printf("NoofCols should be equal to VectorSize\n");
      }
      MPI_Finalize();
      exit(0);
   }

   if(MyRank != 0)
      Vector = (float *)malloc(VectorSize*sizeof(float));
   startcomm = MPI_Wtime();
   MPI_Bcast(Vector, VectorSize, MPI_FLOAT, Root, MPI_COMM_WORLD);
   commtime = MPI_Wtime() - startcomm;

   ScatterSize = NoofRows / Numprocs;
   Mybuffer = (float *)malloc(ScatterSize * NoofCols * sizeof(float));
   startcomm = MPI_Wtime();
   MPI_Scatter( Buffer, ScatterSize * NoofCols, MPI_FLOAT, Mybuffer,
		ScatterSize * NoofCols, MPI_FLOAT, 0, MPI_COMM_WORLD);
   commtime += MPI_Wtime() - startcomm;

   MyFinalVector = (float *)malloc(ScatterSize*sizeof(float));

   for(irow = 0 ; irow < ScatterSize ; irow++) {
       MyFinalVector[irow] = 0;
       index = irow * NoofCols;
       for(icol = 0; icol < NoofCols; icol++)
	   MyFinalVector[irow] += (Mybuffer[index++] * Vector[icol]);
   }

   if( MyRank == 0)
      FinalVector = (float *)malloc(NoofRows*sizeof(float));
   startcomm = MPI_Wtime();
   MPI_Gather( MyFinalVector, ScatterSize, MPI_FLOAT, FinalVector,
	       ScatterSize, MPI_FLOAT, Root, MPI_COMM_WORLD);
   commtime += MPI_Wtime() - startcomm;

   MPI_Reduce(&commtime, &totalcomm, 1, MPI_DOUBLE, MPI_MAX, Root, MPI_COMM_WORLD);

   if( MyRank == 0) {
     end = MPI_Wtime();

     printf("Processor : %d\n",Numprocs);
     printf("Size : %d\n",VectorSize);
     printf("All Time : %f\n",end-start);
     printf("Comm Time : %f\n",totalcomm);
     //printf("Comm Time : %f\n",totalcomm);
     /*
      printf ("\n");
      printf(" --------------------------------------------------- \n");
      printf("Results of Gathering data  %d: \n", MyRank);
      printf("\n");

      for(index = 0; index < NoofRows; index++)
        printf(" FinalVector[%d] = %f \n", index, FinalVector[index]);
      printf(" --------------------------------------------------- \n");
      */
   }

   if(MyRank == 0){
      CheckResultVector = (float *)malloc(NoofRows*sizeof(float));
      for(irow = 0 ; irow < NoofRows ; irow++) {
	  CheckResultVector[irow] = 0;
	  for(icol = 0; icol < NoofCols; icol++){
	      CheckResultVector[irow] += (Matrix[irow][icol]*Vector[icol]);
	  }
	  if(fabs((double)(FinalVector[irow]-CheckResultVector[irow])) >
							    1.0E-10){
	     printf("Error %d\n",irow);
	     ValidOutput = 0;
	  }
       }
       free(CheckResultVector);
       if(ValidOutput)
	  printf("\n-------Correct Result------\n");
   }

   MPI_Finalize();

}
