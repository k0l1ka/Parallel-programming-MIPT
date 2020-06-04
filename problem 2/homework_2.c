#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>

int MAX_SIZE=50000;
double k=1.0;
long double T=0.1;
double START_VAL = 1;
double FRIDGE_VAL = 0;

double solution(double x, double t)
{
	double u, coef, sum, pi;
	int m;
	pi = 3.14159265358979323846;
	coef = 4.0/pi;
	sum = 0;
	for (m=0; m<20; m++)
	{
		sum += exp(-k*pi*pi*(2*m+1)*(2*m+1)*t)/(2*m+1)*sin(pi*(2*m+1)*x);
	}
	u = coef * sum;
	return u;
}

int main(int argc, char *argv[])
{
	int myrank, size, i, N, time_iter;
	MPI_Status Status;

	double begin, end, time;

	MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (myrank==0)
	{
		begin = MPI_Wtime();
	}

	long double h, tau;
	N = atoi(argv[1]);
	tau = atof(argv[2]);
	h = 1.0 / (double)N;	
	N += 1; // number of points
	if (tau == -1)
	{
		tau = 0.47 * h * h;
		T = 1e4 * tau;

	}
	

	time_iter = (int)((float)(T/tau));

	int batch, rest;
	batch = N / size;
	rest = N % size;

	int mem_grid[size];
	int st_grid[2*size];
	int temp = 0;

	for (i = 0; i < size; ++i)
	{
		int add = 0;
		if (rest != 0)
		{
			add = 1;
			rest--;
		}
		st_grid[2*i] = temp;
		st_grid[2*i+1] = st_grid[2*i] + batch + add;
		mem_grid[i] = batch + add;
		temp = st_grid[2*i+1];
	}

	int buf[2];
	MPI_Scatter (st_grid, 2, MPI_INT, buf, 2, MPI_INT, 0, MPI_COMM_WORLD);

	int a, b, m;
	a = buf[0];
	b = buf[1];
	m = b-a;

	double* grid = (double*) calloc(m+2, sizeof(double));

	for (i=0; i<m; i++)
	{
		grid[i+1] = START_VAL;
	}

	for (i=0; i<time_iter; i++)
	{
		int j;

		if ((a==0) && (b!=N))
		{
			grid[0] = FRIDGE_VAL;
			grid[1] = FRIDGE_VAL;
		}
		if ((a!=0) && (b==N))
		{
			grid[m] = FRIDGE_VAL;
			grid[m+1] = FRIDGE_VAL;
		}
		if ((a==0) && (b==N))
		{
			grid[0] = FRIDGE_VAL;
			grid[1] = FRIDGE_VAL;
			grid[m] = FRIDGE_VAL;
			grid[m+1] = FRIDGE_VAL;
		}

		if ((myrank % 2 == 0) && (myrank+1 < size)) 
		{
			MPI_Send(&grid[m], 1, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD);
			MPI_Recv(&grid[m+1], 1, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &Status);
		}
		if (myrank % 2 == 1)
		{
			MPI_Recv(&grid[0], 1, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &Status);
			MPI_Send(&grid[1], 1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
		}

		if ((myrank % 2 == 1) && (myrank+1 < size))
		{
			MPI_Send(&grid[m], 1, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD);
			MPI_Recv(&grid[m+1], 1, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &Status);
		}
		if ((myrank % 2 == 0) && (myrank > 0))
		{
			MPI_Recv(&grid[0], 1, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &Status);
			MPI_Send(&grid[1], 1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
		}

		double prev_val;
		prev_val = grid[0];
		for (j=1; j<m; j++)
		{
			double temp;
			temp = grid[j];
			grid[j] = grid[j] + tau * k * (grid[j+1] - 2 * grid[j] + prev_val) / (h * h);
			prev_val = temp;
		}

		grid[m] = grid[m] + tau * k * (grid[m+1] - 2 * grid[m] + prev_val) / (h * h);
		if (a==0)
		{
			grid[1] = 0;
		}
		if (b==N)
		{
			grid[m] = 0;
		}
	}
	if (myrank != 0)
	{
		MPI_Send(&grid[1], m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	if (myrank==0)
	{
		double parts[MAX_SIZE];
		int index;
		index = mem_grid[0];

		for (i=0; i<mem_grid[0]; i++)
		{
			parts[i] = grid[i+1];
		}

		for (i=1; i<size; i++)
		{
			MPI_Recv(&parts[index], mem_grid[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
			index += mem_grid[i];
		}
		
		end = MPI_Wtime();
	       	time = end - begin;

		int step;
		step = 0.1/h;
		
		for (i=0; i<11; i++)
		{
			double sol;
			int index;
			sol = solution(i*0.1, T);
			index = (int)((float)(i*step));
			printf("x = %f, grid_temp = %f, true_temp = %f, error = %f\n",
				 i*0.1, parts[index], sol, fabs(parts[index]-sol));
		}

		
		
		FILE * fp;
		fp = fopen("graph_data.txt", "a");
		fprintf(fp, "%d, %d, %f\n", N-1, size, time);
		fclose(fp);
	}

	MPI_Finalize();
	return 0;
}

