#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

// function calculus
double f(double x)
{
	return 1 / (1 + x*x);
}

// integral calculus
double integrate(float a, float b, float N)
{
	double s = 0;
	int i;
	for (i = a; i < b; i++)
        {
        	s += 0.5 / N * (f(i / N) + f((i + 1) / N));
        }
	return 4*s;
}

int main(int argc, char* argv[])
{
//       	float rest, step;
//	double S, sum, start1, start2, end1, end2, single_time, group_time;
        int myrank, size, i;
	float portion[2];
	float N;
	double S;
        MPI_Status Status;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&size);

//        printf("I am %d of %d\n", myrank, size);

	if (argc == 1)
                {
                        printf("Number of pieces N not entered!");
                        exit(0);
                }
        else
        	N = atoi(argv[1]);
//                printf("%f pieces of [0,1]\n", N);


        if (myrank == 0)
        {
	        double sum, start1, start2, end1, end2, single_time, group_time;
	        float rest, step;

		start2 = MPI_Wtime();

		rest = (int)N % size;
                step = (N - rest) / size;

//		printf("step = %f, rest = %f\n", step, rest);

		for(i = 1; i < size; i++)
		{
	        	portion[0] = rest + i * step; //start of the part for process number i
			portion[1] = portion[0] + step;  // end of that 
			MPI_Send(portion, 2, MPI_FLOAT, i, i, MPI_COMM_WORLD);
	      	}
//		printf("portions:%f, %f");

		//my part - first
		S = integrate(0, rest+step, N);
                printf("Process №1 part = %.10f\n", S);
	
		sum = S;
		for(i = 1; i < size; i++)
		{
	        	MPI_Recv(&S, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Status);
			sum += S;
			printf("Process №%d part = %.10f\n", i+1, S);	
		}	
		
		printf("Sum of this, I = %.10f\n", sum);

		end2 = MPI_Wtime();
		group_time = end2 - start2;
		
		//full integral
                start1 = MPI_Wtime();

	        S = integrate(0, N, N);
//                printf("I_0 = %f\n", S);

                end1 = MPI_Wtime();
		single_time = end1 - start1;

                printf("(by single process) I_0 = %.10f\n", S);
                printf("S(N=%.0f, p=%d) = %f\n", N, size, single_time / group_time);
 
	}

        if (myrank != 0)
        {
                MPI_Recv(portion, 2, MPI_FLOAT, 0, myrank, MPI_COMM_WORLD, &Status);
		
//		printf("%f %f\n", S, N);
		S = integrate(portion[0], portion[1], N);
//		printf("My piece S_%d = %f\n", myrank, S);

                MPI_Send(&S, 1, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
	
//               	printf("I am %d: %f-%f\n", myrank, portion[0],portion[1]);
        }

        MPI_Finalize();
        return 0;
}
 
