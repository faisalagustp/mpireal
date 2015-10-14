/*
******************************************************************

                 Example 19 (mv_mult_checkerboard.c)

   Objective           : Matrix_Vector Multiplication
                         (Using Block CheckerBoard Partitioning)

   Input               : Process 0 reads files (mdata.inp) for Matrix
                         and (vdata.inp) for Vector

   Output              : Process 0 prints the result of Matrix_Vector
                         Multiplication

   Necessary Condition : Number of processors should be perfect
                         square and less than or equal to 8

******************************************************************
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char** argv)
{

  /* .......Variables Initialisation ......*/
  MPI_Status status;
  MPI_Comm row_comm, col_comm;
  MPI_Group MPI_GROUP_WORLD;

  int Numprocs, MyRank, root_p, Root = 0;
  int irow, icol, iproc, jproc, index, Proc_Id;
  int NoofRows, NoofCols, NoofRows_Bloc, NoofCols_Bloc;
  int Bloc_MatrixSize, Bloc_VectorSize, VectorSize;
  int Local_Index, Global_Row_Index, Global_Col_Index;

  float **Matrix, *Matrix_Array, *Bloc_Matrix, *Vector, *Bloc_Vector;
  float *FinalResult, *MyResult, *FinalVector;

  int *ranks, colsize, colrank, rowsize, rowrank, myrow;
  int		MatrixFileStatus = 1, VectorFileStatus = 1;
  double    start, end, startcomm, endcomm, commtime, totalcomm;


  FILE *fp;


  /* MPI Initialisation  */
  MPI_Init(&argc, &argv);
  MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &Numprocs);
  root_p = sqrt((double)Numprocs);
  if(Numprocs != (root_p * root_p)) {
     if( MyRank == 0)
        printf("Error : Number of processors should be perfect square \n");
     MPI_Finalize();
     exit(-1);
  }

  start = MPI_Wtime();


  if(MyRank == 0)
  {
    VectorSize = 1024*32;
    NoofCols = VectorSize;
    NoofRows = VectorSize;

    NoofRows_Bloc = NoofRows/root_p;
  	NoofCols_Bloc = NoofCols/root_p;

    Matrix = (float **)malloc(NoofRows*sizeof(float *));
    for(irow=0 ;irow<NoofRows; irow++){
      Matrix[irow] = (float *)malloc(NoofCols*sizeof(float));
      for(icol=0; icol<NoofCols; icol++) {
 	       Matrix[irow][icol] = 1;
      }
    }

    //creating vector

    Vector = (float*)malloc(VectorSize*sizeof(float));
    for(index = 0; index<VectorSize; index++)
         Vector[index]=1;


    /* .......Read the Input file ......*/
    if ((fp = fopen ("./data/mdata.inp", "r")) == NULL){
		 MatrixFileStatus = 0;
    }



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

   MPI_Bcast (&NoofRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast (&NoofCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast (&VectorSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

   if(NoofCols != VectorSize){
      if(MyRank == 0){
	 printf("Matrice and Vector Dimensions incompatible for Multiplication");
      }
      MPI_Finalize();
      exit(-1);
   }

   if(NoofRows % root_p != 0 || NoofCols % root_p != 0){
      if(MyRank == 0)
	 printf("Matrix can't be divided among processors equally");
      MPI_Finalize();
      exit(-1);
   }

   NoofRows_Bloc = NoofRows / root_p;
   NoofCols_Bloc = NoofCols / root_p;

   if(MyRank == 0){
      Matrix_Array = (float *)malloc(NoofRows*NoofCols*sizeof(float));
      /* Convert Matrix into 1-D array according to Checkerboard ......*/
      Local_Index = 0;
      for(iproc = 0; iproc < root_p; iproc++){
         for(jproc = 0; jproc < root_p; jproc++){
   	     Proc_Id = iproc * root_p + jproc;
	     for(irow = 0; irow < NoofRows_Bloc; irow++){
		 Global_Row_Index = iproc * NoofRows_Bloc + irow;
	         for (icol = 0; icol < NoofCols_Bloc; icol++){
		      Global_Col_Index = jproc * NoofCols_Bloc + icol;
	              Matrix_Array[Local_Index++] =
			       Matrix[Global_Row_Index][Global_Col_Index];
	         }
	     }
	 }
       }
   }


  /* Memory allocating for Bloc Matrix */
  Bloc_VectorSize = VectorSize / root_p;
  Bloc_MatrixSize = NoofRows_Bloc * NoofCols_Bloc;
  Bloc_Matrix = (float *) malloc (Bloc_MatrixSize * sizeof(float));

  startcomm = MPI_Wtime();
  MPI_Scatter(Matrix_Array, Bloc_MatrixSize, MPI_FLOAT, Bloc_Matrix,
	      Bloc_MatrixSize, MPI_FLOAT, 0, MPI_COMM_WORLD);
  commtime = MPI_Wtime() - startcomm;
  MPI_Barrier(MPI_COMM_WORLD);

  /* Creating groups of procesors row wise */
  myrow=MyRank/root_p;
  startcomm = MPI_Wtime();
  MPI_Comm_split(MPI_COMM_WORLD,myrow,MyRank,&row_comm);
  MPI_Comm_size(row_comm,&rowsize);
  MPI_Comm_rank(row_comm,&rowrank);
  commtime += MPI_Wtime() - startcomm;

  /* Creating groups of procesors column wise */
  myrow=MyRank%root_p;
  startcomm = MPI_Wtime();
  MPI_Comm_split(MPI_COMM_WORLD,myrow,MyRank,&col_comm);
  MPI_Comm_size(col_comm,&colsize);
  MPI_Comm_rank(col_comm,&colrank);
  commtime += MPI_Wtime() - startcomm;

  /* Scatter part of vector to all row master processors */
  Bloc_Vector = (float*) malloc(Bloc_VectorSize * sizeof(float));
  if(MyRank/root_p == 0){
     startcomm = MPI_Wtime();
     MPI_Scatter(Vector, Bloc_VectorSize, MPI_FLOAT, Bloc_Vector,
		 Bloc_VectorSize, MPI_FLOAT, 0,row_comm);
     commtime += MPI_Wtime() - startcomm;
  }

  /* Row master broadcasts its vector part to processors in its column */
  startcomm = MPI_Wtime();
  MPI_Bcast(Bloc_Vector, Bloc_VectorSize, MPI_FLOAT, 0, col_comm);
  commtime += MPI_Wtime() - startcomm;
  /* Multiplication done by all procs */

  MyResult   = (float*) malloc(NoofRows_Bloc * sizeof(float));
  index = 0;
  for(irow=0; irow < NoofRows_Bloc; irow++){
      MyResult[irow]=0;
      for(icol=0;icol< NoofCols_Bloc; icol++){
          MyResult[irow] += Bloc_Matrix[index++] * Bloc_Vector[icol];
      }
  }

   /* collect partial product from all procs on to master processor
		and add it to get final answer */
   if(MyRank == Root)
      FinalResult = (float *)malloc(NoofRows_Bloc*Numprocs*sizeof(float));

   startcomm = MPI_Wtime();
   MPI_Gather (MyResult, NoofRows_Bloc, MPI_FLOAT, FinalResult,
	       NoofRows_Bloc, MPI_FLOAT, 0, MPI_COMM_WORLD);
   commtime += MPI_Wtime() - startcomm;

   MPI_Reduce(&commtime, &totalcomm, 1, MPI_DOUBLE, MPI_MAX, Root, MPI_COMM_WORLD);

/*   if(MyRank == 0){
   for(irow=0; irow < NoofRows_Bloc*Numprocs; irow++)
       printf(" * FinalVector[%d] = %7.3f\n",irow,FinalResult[irow]);
   }
*/

   if(MyRank == 0){


      FinalVector = (float *) malloc(NoofRows * sizeof(float));
      index = 0;
      for(iproc=0; iproc<root_p; iproc++){
	  for(irow=0; irow<NoofRows_Bloc; irow++){
	      FinalVector[index]  = 0;
	      for(jproc=0; jproc<root_p; jproc++){
		  FinalVector[index] +=
		    FinalResult[iproc*root_p*NoofRows_Bloc +
				jproc*NoofRows_Bloc +irow];

	      }
        //printf(" FinalVector[%d] = %7.3f\n",index,FinalVector[index]);
	      index++;
	  }
      }


      end = MPI_Wtime();

      printf("Processor : %d\n",Numprocs);
      printf("Size : %d\n",VectorSize);
      printf("All Time : %f\n",end-start);
      printf("Comm Time : %f\n",commtime);
      //printf("Comm Time : %f\n",totalcomm);

   }


   /* Free the groups formed */
   MPI_Comm_free(&row_comm);
   MPI_Comm_free(&col_comm);


   MPI_Finalize();

}
