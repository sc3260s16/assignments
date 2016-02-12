//*****************************************************************************
//
// This program calculates the product of a square matrix with itself:
//
// B = A * A
//
// Please keep all code in this single file.
//
//
//*****************************************************************************
#include <stdio.h>
#include <stdlib.h>

//It is a good practice to sort this alphabetically (sorry!)

void readfile(FILE *ff, int size, int *matrix);

void mmult(int *matrixA, int *matrixB, int size);

void writefile(int *matrix, int size, FILE *ffo);

int main(int argc, char ** argv)
{
   
   // check command line arguments
   if ( argc != 3 ) {
      printf("This program computes the product of n x n matrix with itself\n");
      printf("Usage: ./matrix_multiply filename n\n");
      exit(0);
   }

   // TODO: parse input arguments
   int sizex=atoi(argv[2]);
   int sizeaux=sizex;
   sizex=sizex*sizex;

   FILE *ff;
   FILE *ffo;
   ff = fopen(argv[1],"r");

   if ( ff == NULL ) {
      printf("input_file.txt not opened, exiting...\n");
      exit(0);
   }

   // TODO: dynamically allocate space for matrix_A (input matrix) in 1d array
   int * matrix_A;  // declare input matrix
   matrix_A = (int*) malloc(sizex*sizeof(int));

   // TODO: call function to read data from file and copy into matrix_A
   readfile(ff,sizeaux, matrix_A);
   fclose(ff);
   
   ffo= fopen("out.dat", "w");

   // TODO: dynamically allocate space for matrix_B (output matrix) in 1d array
   int * matrix_B;  // declare output matrix
   matrix_B = (int*) malloc(sizex*sizeof(int));

   // TODO: call function to perform matrix multiplication ( matrix_B = matrix_A * matrix_A )
   
   mmult(matrix_A, matrix_B, sizeaux);

   // TODO: call function to write results (matrix_B) to stdout
   
   writefile(matrix_B,sizeaux, ffo);
   fclose(ffo);

   // TODO: free space allocated for matrix_A and matrix_B
   
   free(matrix_A);
   free(matrix_B);

   return 0;

}

void readfile(FILE *ff, int size, int *matrix){
   
   int i;  int j;

   for (j=0; j<size; j++){
      for (i=0; i<size; i++) {
         fscanf(ff,"%d\t", &matrix[size*j+i]);
      }
   }
}

void mmult(int *matrixA, int *matrixB, int size){
    
   int i, j, k;
   for(i=0; i<size*size; i++)
      matrixB[i]=0;

   for(i=0; i<size; i++){
      for(j=0; j<size; j++){
         for(k=0; k<size; k++){
            matrixB[j+i*size]+=matrixA[k+i*size]*matrixA[j+k*size];
         }
      }
   }
}

void writefile(int *matrix, int size, FILE *ffo){
    
   int i, j;

   for (j=0; j<size; j++){
      for (i=0; i<size; i++) {
          fprintf(ffo, "%d ", matrix[size*j+i]);
      }
      fprintf(ffo, "\n");
   }
}
