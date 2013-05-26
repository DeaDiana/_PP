/* PLATE_FILLING */
//#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BIG_NUMBER 100000

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

int main()
{
	int Dx_size = 5;
	int Dy_size = 2;
	int Nx_number = 5;
	int Ny_number = 3;
//	int k = i_index * Nx_number + j_index;
	int M_size = Nx_number * Ny_number;
	double hx_step = Dx_size / (double)(Nx_number - 1);
	double hy_step = Dy_size / (double)(Ny_number - 1);
	
	double **matr = (double **) malloc(sizeof(double*) * M_size);

	double *f_vector = (double *) malloc(sizeof(double) * M_size);
    double *x_vector = (double *) malloc(sizeof(double) * M_size);

    if((NULL == matr) || (NULL == f_vector) || (NULL == x_vector)){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
    }

    memset(f_vector, 0, sizeof(double) * (M_size));
    memset(x_vector, 0, sizeof(double) * (M_size));

    for(int i = 0; i < M_size; i++){
      matr[i] = (double*) malloc(sizeof(double) * M_size);
      memset(matr[i], 0, sizeof(double) * (M_size));
      if(NULL == matr[i]){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }
    }	

    generateMatrix(matr, M_size, hx_step, hy_step, Nx_number, Ny_number);
    firstKindBoundaryCondition(matr, f_vector, M_size, Nx_number, Nx_number / 2, Ny_number / 2);

    for(int i = 0; i < M_size; i++)
    {
    	for(int j = 0; j < M_size; j++)
    	{
    		printf("%.2lf ", matr[i][j]);
    	}
    	printf("% lf\n", f_vector[i]);
    }

	return 0;
}