#include <pthread.h>
#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#include <semaphore.h>
#include <time.h>
#include <string.h>

int free_workers = 0;
double h, tau;
double L = 10.0;
double T = 4.0;
int num_threads;
int num_space_steps, num_time_steps;

static double *u;
static double *u_prev;
sem_t sem_count;
sem_t *sem_next;

double g(double x){
    if ((x <= 2.0) && (x >= 0.0)){
        return 2*x-x*x;
    }
    else{
        return 0.0;
    }
}

void* thread_work(void* param) {				//calculates function u at points in its portion during all time
    int i, j, batch, start_point_index, num_points_in_portion;
    start_point_index = ((int *)param)[0];
    num_points_in_portion = ((int *)param)[1] + 1; 		//+1 for extra point on the left
    
    double *local_u = malloc(num_points_in_portion * sizeof(double));

    for (i = 0; i < num_points_in_portion; i++){		//initialization of function 
        local_u[i] = g( 1.0 * (i + start_point_index - 1) * h );
    }

    
    int time;
    double prev, current;

    for (time = 0; time < num_time_steps; time++){
        if (time != 0){
            if (start_point_index > 0){					//the additional point from the left neighbour
                local_u[0] = u_prev[start_point_index-1];	
            }
            else{
                local_u[0] = 0;
            }
        }

        current = local_u[0];
        for (i = 1; i < num_points_in_portion; i++){ 	//one pass over points in portion
            prev = current;
            current = local_u[i];
            local_u[i] = current - tau/h * (current - prev);
        }

        
        for (i=1; i < num_points_in_portion; i++){     //each thread writes to its points
            u[start_point_index+i-1] = local_u[i];
        }
        

        sem_wait(&sem_count);				//all threads should end their work at current layers 
        free_workers++;
        sem_post(&sem_count);

        if (free_workers == num_threads){//code in if block is executed only once by last thread at this time step
            memcpy(u_prev, u, (num_space_steps+1)*sizeof(double));
            free_workers = 0;
            sem_post(&sem_next[time]);
        }

        sem_wait(&sem_next[time]);		//threads finish this iteration one by one after sem_post of the last thread
        sem_post(&sem_next[time]);
    }
    
}

int main(int argc, char *argv[]){
    num_threads = atoi(argv[1]);
    h = atof(argv[2]);
    tau = atof(argv[3]);

    num_time_steps = (int)(T/tau);
    num_space_steps = (int)(10.0/h);

    u = malloc((num_space_steps+1) * sizeof(double));
    u_prev = malloc((num_space_steps+1) * sizeof(double));
    sem_next = malloc(num_time_steps * sizeof(sem_t));
    
    pthread_t *pthr = malloc(num_threads * sizeof(pthread_t));

    int i, rc;
    struct timespec begin, end;
    double elapsed;

    int *buf[num_threads];						//start point number and portion size for each thread
    for (i = 0; i < num_threads; i++)
         buf[i] = (int *)malloc(2 * sizeof(int)); 

    clock_gettime(CLOCK_REALTIME, &begin);

    sem_init(&sem_count, 0, 1);
    for (i = 0; i < num_time_steps; i++){		//for all time steps semafor is locked
        sem_init(&sem_next[i], 0, 0);
    }

    int batch, rest, start_point_index; 
    batch = (num_space_steps+1) / num_threads;
    rest = (num_space_steps+1) % num_threads;
    start_point_index = 0;
    
    for (i = 0; i < num_threads; i++){			//data split among threads 			
        buf[i][0] = start_point_index;
        int add = 0;
        if (rest > 0){
            add = 1;
        }
        buf[i][1] = add + batch;
        start_point_index += buf[i][1];
        rest--;
    }

    for(i = 0; i < num_threads; i++){
        rc = pthread_create(&pthr[i], NULL, thread_work, (void*)buf[i]);
        if (rc) 
        	printf("ERROR; return code from pthread_create() is %d \n", rc);
    }

    for(i = 0; i < num_threads; i++){
        rc = pthread_join(pthr[i], NULL);
        if (rc) 
        	printf("ERROR; return code from pthread_join() is %d \n", rc);
    }

    for (i = 0; i < num_time_steps; i++){
        sem_destroy(&sem_next[i]);
    }
    sem_destroy(&sem_count);

    clock_gettime(CLOCK_REALTIME, &end);

    elapsed = end.tv_sec - begin.tv_sec;
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    //for (i = 0; i < num_space_steps+1; i++){
    //    printf("u(%.1f) = %f, exact_u = %f\n", i*h, u[i], g(i*h-1.0*num_time_steps*tau));
    //}    

    FILE * fp;
    fp = fopen("plot_data.txt", "a");
    fprintf(fp, "%d, %f\n", num_threads, elapsed);
    fclose(fp);

    free(pthr);
    return 0;
}


