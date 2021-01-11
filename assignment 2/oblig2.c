#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>


void read_matrix_binaryformat (char* filename, double*** matrix,int* num_rows,
			       int* num_cols);
void write_matrix_binaryformat (char* filename, double** matrix,int num_rows, 
				int num_cols);
double* convert_to_vector(double **matrix,double *vector,int m,int n);
double** convert_to_matrix(double *vector,double **matrix,int m,int n);
double** matrixMultiply(double **matrixA,double **matrixB,int m,int n,int o);
void testResults(double **matrix,int m,int n);

int main(int argc, char *argv[]){

  /* matrix information */
  int a_m,a_n,b_m,b_n;
  int my_m1,my_m2,my_n1,my_n2; /* matrix partitioning */
  int num_rowsA,num_colsA,num_rowsB,num_colsB;
  char *filenameA,*filenameB,*filenameC;
  double **matrixA,**matrixB,**matrixC;

  double *myAvector,*myBvector,*aVector,*bVector,*cVector;
  double *myCvector;

  /* MPI Communications */
  int my_rank,num_procs;
  int sender,start_number1,end_number1; /* which process is sending to master */
  int start_number2,end_number2,start_number3;
  int number_to_receive, number_received;
  int number_of_data_to_send1,number_of_data_to_send2;

  
  //unsigned char *my_received_chars;
  MPI_Status status;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);  
  
  if (my_rank==0) {
    if (argc!=3){
      fprintf(stderr, "Use ./program <matrix A filename> <matrix B filename> \n");
      return -1;
    }
    // reads to filenames from the terminal
    filenameA=argv[1];
    filenameB=argv[2];

    read_matrix_binaryformat (filenameA,&matrixA,&num_rowsA,&num_colsA);
    read_matrix_binaryformat (filenameB,&matrixB,&num_rowsB,&num_colsB);

    printf("\nThe array A has %d rows and %d columns\n",num_rowsA,num_colsA);
    printf("The array B has %d rows and %d columns\n\n",num_rowsB,num_colsB);
    
    /* makes vectors containing all A and B values */
    aVector=convert_to_vector(matrixA,aVector,num_rowsA,num_colsA);
    bVector=convert_to_vector(matrixB,bVector,num_rowsB,num_colsB);
    cVector=malloc((num_rowsA*num_colsB)*sizeof(double));


    //printf("A first element is %f\n",matrixA[0][0]);
    //printf("aVector first element is %f\n",aVector[0]);
  }
  
  
  // Lets the processes know what sizes the arrays are 
  MPI_Bcast (&num_rowsA, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&num_colsA, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&num_rowsB, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&num_colsB, 1, MPI_INT, 0, MPI_COMM_WORLD);



  // divide the A and B matrix evenly among the MPI processes 
  my_m1 = num_rowsA/num_procs;
  my_m2 = num_rowsB/num_procs;
  my_n1 = num_colsA;
  my_n2 = num_colsB;
  // Lets the last process do the surplus part
  if (my_rank==(num_procs-1)){
    my_m1 = num_rowsA/num_procs + num_rowsA%num_procs;
    my_m2 = num_rowsB/num_procs + num_rowsB%num_procs;
  } 

  /* Allocates space for the matrices and arrays */
  myAvector = malloc((my_m1*my_n1)*sizeof(double));
  myBvector = malloc((my_m2*my_n2)*sizeof(double));

  int i=0;
  int *receive=malloc(num_procs*sizeof(int));
  int *disp=malloc(num_procs*sizeof(int));
  //int* sendCount=malloc(
  for (i=0;i<num_procs;i++){
    receive[i]=(num_rowsA/num_procs)*num_colsB;
    if (i==(num_procs-1)){
      receive[i]=(num_rowsA/num_procs + (num_rowsA%num_procs))*num_colsB;
    }
    disp[i]=i*(num_rowsA/num_procs)*num_colsB;
  }





  int i3=0;
  int *receive2=malloc(num_procs*sizeof(int));
  int *disp2=malloc(num_procs*sizeof(int));

  for (i3=0;i3<num_procs;i3++){
    receive2[i3]=(num_rowsB/num_procs)*num_colsB;
    if (i3==(num_procs-1)){
      receive2[i3]=(num_rowsB/num_procs + (num_rowsB%num_procs))*num_colsB;
    }
    disp2[i3]=i3*(num_rowsB/num_procs)*num_colsB;
  }




  if (my_rank==0){
    // assigns values for my part of matrix B
    int i2;
    for (i2=0;i2<(my_m2*my_n2);i2++){
      myBvector[i2]=bVector[i2];
    }


    /* divides the matrix data among the processes */
    int id=0;
    for (id=1;id<num_procs;id++){
      start_number1=id*(num_rowsA/num_procs)*num_colsA;
      end_number1=((id+1)*(num_rowsA/num_procs)*num_colsA)-1;
      start_number2=id*(num_rowsB/num_procs)*num_colsB;
      end_number2=((id+1)*(num_rowsB/num_procs)*num_colsB)-1;      
      if (id==(num_procs-1)){
	end_number1=(num_rowsA*num_colsA)-1;
	end_number2=(num_rowsB*num_colsB)-1;
      }
      number_of_data_to_send1=end_number1-start_number1+1;
      number_of_data_to_send2=end_number2-start_number2+1;
      
      MPI_Send(&number_of_data_to_send1,1,MPI_INT,id,1,MPI_COMM_WORLD);
      MPI_Send(&aVector[start_number1],number_of_data_to_send1,MPI_DOUBLE,
             id,1,MPI_COMM_WORLD);

      /*
      // sends entire B vector
      int send=num_rowsB*num_colsB;
      MPI_Send(&send,1,MPI_INT,id,1,MPI_COMM_WORLD);
      MPI_Send(&bVector[0],(num_rowsB*num_colsB),MPI_DOUBLE,
             id,1,MPI_COMM_WORLD);
      */


      // sends parts of B to all processes
      MPI_Send(&number_of_data_to_send2,1,MPI_INT,id,1,MPI_COMM_WORLD);
      MPI_Send(&bVector[start_number2],number_of_data_to_send2,MPI_DOUBLE,
             id,1,MPI_COMM_WORLD);
      
      



      printf("I have sent data to process %d\n",id);
    }
    MPI_Allgatherv(myBvector,my_m2*num_colsB,MPI_DOUBLE,bVector,receive2,disp2,MPI_DOUBLE,MPI_COMM_WORLD); 

  
    int i,counter=0;
    for (i=0;i<(my_m1*my_n2);i++){
      myAvector[counter]=aVector[counter];
      counter++;
    }
    double** my_Amatrix=convert_to_matrix(myAvector,my_Amatrix,my_m1,my_n1);
    double** my_Cmatrix=matrixMultiply(my_Amatrix,matrixB,my_m1,my_n1,num_colsB);

    myCvector=convert_to_vector(my_Cmatrix,myCvector,my_m1,num_colsB);

    /* make my a and b vectors */
    // adds the c values to the final vector
    //counter=0;
    for (i=0;i<(my_m1*num_colsB);i++){
      //cVector[i]=myCvector[i];
    }

    /*
    for (id=1;id<num_procs;id++){
       MPI_Recv(&number_to_receive,1,MPI_INT,id,1,MPI_COMM_WORLD,
      	       &status);
      sender=status.MPI_SOURCE;
      start_number3=id*(num_rowsA/num_procs)*num_colsB;
      MPI_Recv(&cVector[start_number3],number_to_receive,MPI_DOUBLE,sender,1,MPI_COMM_WORLD,&status);
      } */


  }
  else{
    /* I am a slave who receives data from master and does my own work */
    MPI_Recv(&number_received,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
    MPI_Recv(&myAvector[0],number_received,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&status);
    double** my_Amatrix=convert_to_matrix(myAvector,my_Amatrix,my_m1,my_n1);
    
    MPI_Recv(&number_received,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
    MPI_Recv(&myBvector[0],number_received,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&status);
 
    bVector=malloc((num_rowsB*num_colsB)*sizeof(double));
    MPI_Allgatherv(myBvector,my_m2*num_colsB,MPI_DOUBLE,bVector,receive2,disp2,MPI_DOUBLE,MPI_COMM_WORLD); 
    /*
    bVector=malloc((num_rowsB*num_colsB)*sizeof(double));
    MPI_Recv(&number_received,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
    MPI_Recv(&bVector[0],number_received,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&status);   
    */

    // converts bVector to matrix
    matrixB=convert_to_matrix(bVector,matrixB,num_rowsB,num_colsB);

    double** my_Cmatrix=matrixMultiply(my_Amatrix,matrixB,my_m1,my_n1,num_colsB);

    myCvector=convert_to_vector(my_Cmatrix,myCvector,my_m1,num_colsB);

    /*
    int send=my_m1*num_colsB;
    MPI_Send(&send,1,MPI_INT,0,1,MPI_COMM_WORLD);
    MPI_Send(&myCvector[0],send,MPI_DOUBLE,
             0,1,MPI_COMM_WORLD);   
    */

    /* need to get all parts of B, 0 already has it */
    /* can allocate entire B vector and receive my part into the vector*/

  }


  //number_received=my_m1*num_colsB;
  MPI_Gatherv(myCvector,my_m1*num_colsB,MPI_DOUBLE,cVector,receive,disp,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (my_rank==0){  
    // done computing matrix C
    matrixC=convert_to_matrix(cVector,matrixC,num_rowsA,num_colsB);    
    write_matrix_binaryformat("testfile.bin",matrixC,num_rowsA,num_colsB);
  
    //matrixC=matrixMultiply(matrixA,matrixB,num_rowsA,num_colsA,num_colsB);
    testResults(matrixC,num_rowsA,num_colsB);
  }
  // deallocate matrices


  MPI_Finalize();
}

double** matrixMultiply(double **matrixA,double **matrixB,int m,int n,int o){
  double **matrixC=(double**)malloc(m*sizeof(double*));
  int i=0;
  for (i=0;i<m;i++){
    matrixC[i]=(double*)malloc(o*sizeof(double));
  }
  int k,j;
  double sum;

#pragma omp parallel shared(matrixA,matrixB,matrixC) private(i,j,k) num_threads(4)
  { // start parallel region
#pragma omp for schedule(static)
  for(i=0;i<m;i++){
    for(j=0;j<o;j++) {  
      matrixC[i][j] = 0.0;  
      for(k=0;k<n;k++){  
	matrixC[i][j] = matrixC[i][j] + matrixA[i][k]*matrixB[k][j];  
      }
    }
  }
  
  } // end parallel region
  
  return matrixC;
}
void read_matrix_binaryformat (char* filename, double*** matrix,int* num_rows,
			       int* num_cols){
  int i;
  FILE* fp = fopen (filename,"rb");
  fread (num_rows, sizeof(int), 1, fp);
  fread (num_cols, sizeof(int), 1, fp);
  /* storage allocation of the matrix */
  *matrix = (double**)malloc((*num_rows)*sizeof(double*));
  (*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
  for (i=1; i<(*num_rows); i++)
    (*matrix)[i] = (*matrix)[i-1]+(*num_cols);
  /* read in the entire matrix */
  fread ((*matrix)[0], sizeof(double), (*num_rows)*(*num_cols), fp);
  fclose (fp);
}
void write_matrix_binaryformat (char* filename, double** matrix,int num_rows, 
				int num_cols){
  FILE *fp = fopen (filename,"wb");
  fwrite (&num_rows, sizeof(int), 1, fp);
  fwrite (&num_cols, sizeof(int), 1, fp);
  fwrite (matrix[0], sizeof(double), num_rows*num_cols, fp);
  fclose (fp);
}
double** convert_to_matrix(double *vector,double **matrix,int m,int n){
  matrix=(double**)malloc(m*sizeof(double*));
  int i=0;
  for (i=0;i<m;i++){
    matrix[i]=(double*)malloc(n*sizeof(double));
  }
  int j=0,k=0;
  int counter=0;
  for (j=0;j<m;j++){
    for (k=0;k<n;k++){
      matrix[j][k]=vector[counter];
      counter++;
    }
  }
  return matrix;
}
double* convert_to_vector(double **matrix,double *vector,int m,int n){
  //printf("m is %d and n is %d\n",m,n);
  vector = malloc((m*n)*sizeof(double));
  int i,j,counter=0;
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      vector[counter]=matrix[i][j];
      counter++;
    }
  }
  return vector;
}
void testResults(double **matrix,int m,int n){
  double **fasit;
  int m2,n2;
  char* name="small_matrix_c.bin";
  read_matrix_binaryformat(name,&fasit,&m2,&n2);
  int i,j;
  for (i=0;i<m;i++){
    for (j=0;j<2;j++){
      if (fasit[i][j]!=matrix[i][j]){
	printf("ooops %lf forskjellig fra fasit %lf \n",matrix[i][j],fasit[i][j]);
	printf("coordinates are i=%d j=%d\n",i,j);
      }
    }
  }
}
