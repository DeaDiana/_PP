#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define COMMW MPI_COMM_WORLD
//#define HEIGHT 10
//#define WIDTH 10
#define EPS 0.01
#define T 0.01
#define DBG 0
#define CHECK 0
#define FROM_FILE 0
#define BIG_NUMBER 100000

int HEIGHT = 10;
int WIDTH = 10;

double norma(double*);

void generateMatrix(double** matr, int M_size, double hx_step, double hy_step, int Nx_number, int Ny_number)
{
    double a_coef = -2 * (1 / (hx_step * hx_step) + 1 / (hy_step * hy_step));
    double b_coef = 1 / (hx_step * hx_step);
    double c_coef = 1 / (hy_step * hy_step);

    printf("%d M_size, %lf hx_step, %lf hy_step, %d Nx_number, %d Ny_number\n", M_size, hx_step, hy_step, Nx_number, Ny_number);

    for(int i = 0; i < M_size; i++)
    {
        matr[i][i] = a_coef;
        if((i + 1) < M_size)
        {
            matr[i][i + 1] = b_coef;
            matr[i + 1][i] = b_coef;
        }
        if((i + Nx_number) < M_size)
        {
            matr[i][i + Nx_number] = c_coef;
            matr[i + Nx_number][i] = c_coef;
        }
    }

    printf("%lf\n", matr[3][4]);
}

void firstKindBoundaryCondition(double** matr, double *f_vector, int M_size, int Nx_number, int point_x, int point_y)
{
    int k = point_x * Nx_number + point_y;

    matr[k][k] = BIG_NUMBER;

    f_vector[k] = 1.0 * BIG_NUMBER;
}

int main(int argc, char *argv[]){

    FILE *file;
    int size = 0, rank = 0, i = 0, j = 0, iter = 0;
    int root = 0;
    int div_h = 0, mod_h = 0;
    double **partOfMatr = NULL, **matr = NULL;
    double *x_vector = NULL;
    double *loc_res = NULL, *result = NULL;
    int *loc_hghts = NULL, *starts_v = NULL;
    double *f_vector = NULL, *partOfF_vector = NULL;
    MPI_Status status;
    double prescission = 0.0;
    //==
    int Dx_size = 5;
    int Dy_size = 2;
    int Nx_number = 5;
    int Ny_number = 3;
    int M_size = Nx_number * Ny_number;
    double hx_step = Dx_size / (double)(Nx_number - 1);
    double hy_step = Dy_size / (double)(Ny_number - 1);
    HEIGHT = M_size;
    WIDTH = M_size;
    //==

    matr = (double **) malloc(sizeof(double*) * HEIGHT);
    x_vector = (double *) malloc(sizeof(double) * WIDTH);
    result = (double *) malloc(sizeof(double) * HEIGHT);
    f_vector = (double *) malloc(sizeof(double) * WIDTH);

    if((NULL == matr) || (NULL == x_vector) || (NULL == result) || (NULL == f_vector)){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
    }
    memset(f_vector, 0, sizeof(double) * (M_size));
    memset(x_vector, 0, sizeof(double) * (M_size));

    for(i = 0; i < HEIGHT; i++){
      matr[i] = (double*) malloc(sizeof(double) * WIDTH);
      memset(matr[i], 0, sizeof(double) * (M_size));
      if(NULL == matr[i]){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }
    }

    #if  FROM_FILE
    if(NULL == (file = fopen("matrix", "rb"))) {
       fprintf(stderr, "file with source matrix wasn't opened\n");
        exit(1); 
    }
                    //fill the source matrix and x_vector with numbers
    for(i = 0; i < HEIGHT; i++){
	for(j = 0; j < WIDTH; j++){
    	    fscanf(file, "%lf", &matr[i][j]);
	}
    }

    for(j = 0; j < WIDTH; j++){
        fscanf(file, "%lf", &f_vector[j]);
	x_vector[j] = f_vector[j];
    }
    fclose(file);
    #else
        generateMatrix(matr, M_size, hx_step, hy_step, Nx_number, Ny_number);
        firstKindBoundaryCondition(matr, f_vector, M_size, Nx_number, Nx_number / 2, Ny_number / 2);
    #endif

    for(int i = 0; i < M_size; i++)
    {
        for(int j = 0; j < M_size; j++)
        {
            printf("%.2lf ", matr[i][j]);
        }
        printf("% lf\n", f_vector[i]);
    }

    prescission = norma(f_vector) * EPS;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(COMMW, &size);
    MPI_Comm_rank(COMMW, &rank);

    loc_hghts = (int *) malloc(sizeof(int) * size); //heights of matrx's parts
    starts_v = (int *) malloc(sizeof(int) * (size + 1)); //the same navigation in vertical

    if((NULL == loc_hghts) || (NULL == starts_v)){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }

    div_h = HEIGHT / size;
    mod_h = HEIGHT % size;

    for(i = 0; i < size; i++){
    loc_hghts[i] = (i < mod_h) ? div_h + 1 : div_h; //first mod_h take +1 more row from the matrix
    starts_v[i] = (i <= mod_h) ? (div_h + 1) * i : (div_h + 1) * mod_h + div_h * (i - mod_h); // next process starts righter from previous
    }						//first mod_w shifts on steps (div_w + 1), the rest of other steps == (div_w + 1) and (div_w)
    starts_v[size] = HEIGHT; //bottom range

    partOfMatr = (double **) malloc(sizeof(double*) * loc_hghts[rank]);//allocating memory for working parts of matrix
    loc_res = (double *) malloc(sizeof(double) * (div_h + 1));	//and local result: size (div_h + 1) for any prcess (+1 - reserve)
								//so we take wich reserving (div_h + 1) for each process
    partOfF_vector = (double *) malloc(sizeof(double) * loc_hghts[rank]); //particular loc_hghts[rank] elements from source f_vector

    if((NULL == partOfMatr) || (NULL == loc_res) || (NULL == partOfF_vector)){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }

    for(i = 0; i < loc_hghts[rank]; i++){
      partOfMatr[i] = (double*) malloc(sizeof(double) * WIDTH);
      if(NULL == partOfMatr[i]){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }
    }

    while(1){

		    //fill numbers from the source matrix
    for(i = 0; i < loc_hghts[rank]; i++){
	for(j = 0; j < WIDTH; j++){
	    partOfMatr[i][j] = matr[i + starts_v[rank]][j];
	}
	partOfF_vector[i] = f_vector[i + starts_v[rank]]; //corect for square matrixs
    }

		    //production part of matrix on whole x_vector
    for(i = 0; i < loc_hghts[rank]; i++){
	for(j = 0; j < WIDTH; j++){
	    partOfMatr[i][j] *= x_vector[j];
        }
    }

    memset(loc_res, 0, sizeof(double) * (div_h + 1)); //length for local result with reserving 1 element
    for(i = 0; i < loc_hghts[rank]; i++){	//but fill only sense elements
	for(j = 0; j < WIDTH; j ++){
    	    loc_res[i] += partOfMatr[i][j]; //summing up for local result
	}
	loc_res[i] -= partOfF_vector[i];	// Axk - b

    }

    MPI_Allgatherv(loc_res, loc_hghts[rank], MPI_DOUBLE, //send particular loc_hghts[rank] - sometime less than div_h + 1
	    result, loc_hghts, starts_v, // loc_hghts contains lenghts of local results, starts contains their positions in result x_vector
	    MPI_DOUBLE, COMMW); // result will have exactly HEIGHT elements

    iter++;
    //chek norma of result <= EPS* || b || to finish algorithm or to continue
	if((norma(result) <= prescission) || (200 < iter)) {
	    break;	
	}
	else {
	    	for(i = 0; i < WIDTH; i++){
		    x_vector[i] = x_vector[i] - T * result[i];
		}
	}
    }

    if(root == rank){
        for(i = 0; i < WIDTH; i++){
	    printf("%f \n", x_vector[i]);
	}
	printf("\n the result is above \n");
    }

    for(i = 0; i < loc_hghts[rank]; i++){
      free(partOfMatr[i]);
    }
    free(partOfMatr);

    free(loc_hghts);
    free(starts_v);
    free(loc_res);
    free(partOfF_vector);

    MPI_Finalize();

    for(i = 0; i < HEIGHT; i++){
      free(matr[i]);
    }
    free(matr);
    free(x_vector);
    free(result);
    free(f_vector);

return 0;
}

double norma(double *v){

    double max = 0.0, compare = 0.0;
    int i = 0;

    for(i = 0; i < HEIGHT; i++){
	compare = (0 < v[i]) ? v[i] : -v[i];
        if(compare > max){
	     max = compare;
	}
    }

return max;
}