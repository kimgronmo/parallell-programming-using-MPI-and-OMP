#include <stdio.h>
#include <stdlib.h>
//#include <string.h>

/* need to define functions used later in program */
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
void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters);




int main (int argc, char *argv[]) {
  /* write alle comments like this */
  /* defines variables */
  int m,n,c,iters;
  float kappa; 
  image u,u_bar;
  unsigned char *image_chars;
  char *input_jpeg_filename, *output_jpeg_filename;

  /* read info from command line */
  /* Exits and gives error if not enough arguments are given */
  if (argc!=5){
    fprintf(stderr, "Use ./program <iterasjoner> <kappa> <infile> <outfile>\n");
    return -1;
  }
  /* assigns values to variables */
  iters=atoi(argv[1]); /* converts string (chars) to int */
  kappa=atof(argv[2]); /* converts string to float (or double?) */
  input_jpeg_filename=argv[3];
  output_jpeg_filename=argv[4];

  /* test to see if read correctly
  printf("kappa is %f\n",kappa); 
  printf(input_jpeg_filename);
  printf(output_jpeg_filename);
  */

  import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
  printf("The picture is %d pixels high and %d pixels wide\n",m,n);
  printf("Starting denoising. Doing %d iterations: \n",iters);


  allocate_image(&u,m,n);
  allocate_image(&u_bar,m,n);

  convert_jpeg_to_image(image_chars, &u);

  iso_diffusion_denoising(&u, &u_bar,kappa,iters);

  /* since all values are copied back into u at the last iteration */
  /* I can send u instead of u_bar to be converted back */
  /* convert_image_to_jpeg(&u_bar, image_chars); */
  convert_image_to_jpeg(&u, image_chars);


  export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
  deallocate_image(&u);
  deallocate_image(&u_bar);
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
  //printf(u->image_data[0][1]);
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
void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters){
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
  }
}
