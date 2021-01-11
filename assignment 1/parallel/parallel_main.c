#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* make use of two functions from the simplejpeg library */
void import_JPEG_file(const char *filename, unsigned char **image_chars,
		      int *image_height, int *image_width,
		      int *num_components);
void export_JPEG_file(const char *filename, unsigned char *image_chars,
		      int image_height, int image_width,
		      int num_components, int quality);

typedef struct{
  float** image_data; // a 2D array of floats
  int m; // pixels in x-direction
  int n; // pixels in y-direction ?? eller omvendt
} image;

/* note image struct have to be above in code */ 
void allocate_image(image *u, int m, int n);
void deallocate_image(image *u);

void convert_jpeg_to_image(const unsigned char* image_chars, image *u);
void convert_image_to_jpeg(const image *u, unsigned char* image_chars);
void iso_diffusion_denoising_parallel(image *u, image *u_bar, float kappa, int iters, int my_rank);

void add_to_whole(image *u,image *whole_image,int process_nr, int m, int n);

int main(int argc, char *argv[]){

  int m, n, c, iters;
  int my_m, my_n, my_rank, num_procs;
  float kappa;
  image u, u_bar, whole_image;
  unsigned char *image_chars, *my_image_chars;
  char *input_jpeg_filename, *output_jpeg_filename;

  /* MPI communications */
  int sender,start_number,end_number; /* which process is sending to master */
  int number_to_receive, number_received;
  int number_of_data_to_send;
  unsigned char *my_received_chars;
  MPI_Status status;

  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

  /* read from command line: kappa, iters, input_jpeg_filename, output_jpeg_filename */
  /* Exits and gives error if not enough arguments are given */

  /* use mpirun -np 4 parallel.main 100 0.1 infile.jpg outfile.jpg */
  /* to run with 4 processes, 100 iterations, kappa=0.1 */

  /* assigns values to variables */
  iters=atoi(argv[1]); /* converts string (chars) to int */
  kappa=atof(argv[2]); /* converts string to float (or double?) */
  input_jpeg_filename=argv[3];
  output_jpeg_filename=argv[4];

  if (my_rank==0) {
    if (argc!=5){
      fprintf(stderr, "Use ./program <iterasjoner> <kappa> <infile> <outfile>\n");
      return -1;
    }
    import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
    allocate_image (&whole_image, m, n);
    printf("\nThe picture is %d pixels high and %d wide\n",m,n);
  }

  /* Lets the processes know what m and n is */
  MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* divide the m x n pixels evenly among the MPI processes */
  my_m = m/num_procs;
  my_n = n;
  /* Lets the last process do the surplus pixels */
  if (my_rank==(num_procs-1)){
    my_m = m/num_procs + m%num_procs;
  } 

  /* In order to prevent lines in image we must have overlap */
  /* between bordering regions ????? */
  /* overlaps bottom own and top row the one below */

  if ((num_procs>1)){
    my_m=my_m+2; 
  }
  if (my_rank==(num_procs-1)){
    my_m=m/num_procs + m%num_procs+1;
  }
  if (my_rank==0){
    my_m=my_m+1;
  }

  allocate_image (&u, my_m, my_n);
  allocate_image (&u_bar, my_m, my_n);

  my_image_chars=malloc((my_m*my_n)*sizeof(unsigned char));

  /* each process asks process 0 for a partitioned region */
  /* of image_chars and copy the values into u */
  /* ... */
  if (my_rank==0){
    /* divides image data among the processes */
    int id=0;
    for (id=1;id<num_procs;id++){

      start_number=id*(m/num_procs)*n-n; /* upper row overlap */

      /* adds two rows for overlap ?? */
      end_number= ((id+1)*(m/num_procs)*n)-1+n; /* bottom row overlap */
      if (id==(num_procs-1)){
	end_number=(m*n)-1;
      }
      number_of_data_to_send=end_number-start_number+1;

      MPI_Send(&number_of_data_to_send,1,MPI_INT,id,1,MPI_COMM_WORLD);
      MPI_Send(&image_chars[start_number],number_of_data_to_send,MPI_UNSIGNED_CHAR,
	       id,1,MPI_COMM_WORLD);
      printf("I have sent data to process %d\n",id);
    }

    int i=0;
    for (i=0;i<(my_m*my_n);i++){
      my_image_chars[i]=image_chars[i];
    }

    convert_jpeg_to_image(my_image_chars, &u);
    iso_diffusion_denoising_parallel(&u, &u_bar, kappa, iters,my_rank);
    convert_image_to_jpeg(&u,image_chars);   

    id=0;

    for (id=1;id<num_procs;id++){
      /* printf("Trying to retrieve data from slave %d\n",id); */
      MPI_Recv(&number_to_receive,1,MPI_INT,id,1,MPI_COMM_WORLD,
      	       &status);
      sender=status.MPI_SOURCE;
      my_received_chars=(unsigned char *)malloc(number_to_receive*sizeof(unsigned char));

      MPI_Recv(&my_received_chars[0],number_to_receive,MPI_UNSIGNED_CHAR,sender,1,MPI_COMM_WORLD,&status);

      /* the process above has done your top row */
      /* starts at second row which is done by the process itself */
      int j=0;
      for (j=0;j<(number_to_receive-2*n);j++){
	image_chars[(id*(m/num_procs)*n)+j]=my_received_chars[j+n]; 
      }
      free(my_received_chars);
    }

    /* All work is done, lets make the image */
    /* overwrites work done by slaves. add my chars to image chars instead */
    /* for process 0 */
    /* u imagetojpeg instead?*/
    //convert_image_to_jpeg(&whole_image, image_chars);
    export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
    deallocate_image (&whole_image); //remove this??
   }
  else {
    /* I am a slave who receives data from my master and does my own work */
    /* when done I return data to master */
    MPI_Recv(&number_received,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
    my_image_chars=(unsigned char *)malloc(number_received*sizeof(unsigned char));
    MPI_Recv(&my_image_chars[0],number_received,MPI_UNSIGNED_CHAR,0,1,MPI_COMM_WORLD,&status);

    /* printf("received data from master %d\n",number_received); */
    convert_jpeg_to_image(my_image_chars, &u);
    iso_diffusion_denoising_parallel(&u, &u_bar, kappa, iters,my_rank);
    convert_image_to_jpeg(&u,my_image_chars); // now contains fixed data

    /* Easier to send data in one batch as unsigned chars since */
    /* they can be directly put into image chars and converted to jpeg */

    MPI_Send(&number_received,1,MPI_INT,0,1,MPI_COMM_WORLD);
    MPI_Send(&my_image_chars[0],number_received,MPI_UNSIGNED_CHAR,
	       0,1,MPI_COMM_WORLD);

    free(my_image_chars);
    printf("Process %d is done. Partial image returned\n",my_rank);
  }
  if (my_rank==0){
    printf("All processes done. Creating image\n");
  }

  deallocate_image (&u);
  deallocate_image (&u_bar);
  MPI_Finalize ();
  return 0;
}
void allocate_image(image *u, int m, int n){
  /* make double array from m and n (hva er x og y ?? */
  /* temp array sett temp array --> u sin array */
  /* set u sin m og n */
  u->m=m;
  u->n=n;
  u->image_data = (float**)malloc(m*sizeof(float*));
  int i=0;
  for (i=0;i<m;i++){
    u->image_data[i]=(float*)malloc(n*sizeof(float));
  }
}
void deallocate_image(image *u){
  /* free up all parts of the image */
  /* dynamically allocated arrays */
  int i=0;
  for (i=0;i<(u->m);i++){
    free(u->image_data[i]);
  }
  free(u->image_data);
}
void convert_jpeg_to_image(const unsigned char* image_chars, image *u){
  int counter=0;
  int i=0,j=0;
  for (i=0;i<(u->m);i++){
    for (j=0;j<(u->n);j++){
      /* can convert unsigned char directly to float??   */
      u->image_data[i][j]=(float)image_chars[counter];
      counter++;
    }
  }
}
void convert_image_to_jpeg(const image *u, unsigned char* image_chars){
  int i3=0,j3=0,counter2=0;
  for (i3=0;i3<(u->m);i3++){
    for (j3=0;j3<(u->n);j3++){
      image_chars[counter2]=(char)u->image_data[i3][j3];
      counter2++;
    }
  }
}
void iso_diffusion_denoising_parallel(image *u, image *u_bar, float kappa, int iters,int my_rank){
  int i=1,j=1;
  /* frame of boundary pixels around picture remain unchanged */
  int k=0; /* iterasjoner */
  for (k=0;k<iters;k++){
    for (i=1;i<(u->m-1);i++){ /* -1 due to border not included */
      for (j=1;j<(u->n-1);j++){
	u_bar->image_data[i][j]=u->image_data[i][j]+kappa*(u->image_data[i-1][j]+u->image_data[i][j-1]
							   -4*u->image_data[i][j]+u->image_data[i][j+1]
							   +u->image_data[i+1][j]);
      }
    }
    /* if i copy u_bar values to u at the end I can send u to char array when finished */
    /* no need to copy boundary pixels   ??  */
    
    int i2=0,j2=0;
    for (i2=1;i2<(u->m-1);i2++){
      for (j2=1;j2<(u->n-1);j2++){
	u->image_data[i2][j2]=u_bar->image_data[i2][j2];
      }
    } 

    /* prøvde å swappe pekere men programmet crasher */ 
    /* u->image_data=u_bar->image_data; */

    /* noe MPI info som trengs */
    MPI_Status status;
    int num_procs;
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
    
    /* prosessene sender overlap til hverandre etter hver iterasjon */
    if (my_rank == 0){
      /* sender min nest nederste til den under sin øverste */
      MPI_Send(&u->image_data[u->m-2][0],u->n,MPI_FLOAT,my_rank+1,1,MPI_COMM_WORLD);
      /* mottar den under sin nest øverste som min nederste */
      MPI_Recv(&u->image_data[u->m-1][0],u->n,MPI_FLOAT,my_rank+1,2,MPI_COMM_WORLD,&status);
    }
    if (my_rank==(num_procs-1)){
    /* mottar til min øverste fra den over sin nest nederste */
      MPI_Recv(&u->image_data[0][0],u->n,MPI_FLOAT,my_rank-1,1,MPI_COMM_WORLD,&status);
      /* sender min nest øverste til den over sin nederste */
      MPI_Send(&u->image_data[1][0],u->n,MPI_FLOAT,my_rank-1,2,MPI_COMM_WORLD);
    }
    if (my_rank>0 && my_rank<(num_procs-1)){
      /* mottar til min øverste fra den over sin nest nederste */
      MPI_Recv(&u->image_data[0][0],u->n,MPI_FLOAT,my_rank-1,1,MPI_COMM_WORLD,&status);
      /* sender min nest øverste til den over sin nederste */
      MPI_Send(&u->image_data[1][0],u->n,MPI_FLOAT,my_rank-1,2,MPI_COMM_WORLD);

      /* sender min nest nederste til den under sin øverste */
      MPI_Send(&u->image_data[u->m-2][0],u->n,MPI_FLOAT,my_rank+1,1,MPI_COMM_WORLD);
      /* mottar den under sin nest øverste som min nederste */
      MPI_Recv(&u->image_data[u->m-1][0],u->n,MPI_FLOAT,my_rank+1,2,MPI_COMM_WORLD,&status);
    }
  }
}
/* remove this?? */
/* cannot send image u in mpi */
void add_to_whole(image *u,image *whole_image,int process_nr,int m,int n){
  /* Where to start filling the whole image array */
  int counter=process_nr*m*n;
  int i=0,j=0;
  for (i=0;i<(u->m-2);i++){
      for (j=0;j<(u->n-2);j++){
	whole_image->image_data[counter+i][j]=u->image_data[i][j];
      }
  }
}
