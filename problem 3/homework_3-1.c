#include <pthread.h>
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <semaphore.h>
#include <time.h>

int inside_counter = 0;
int TOTAL_POINTS = 1e6;
int NUM_THREADS = 7;
float PI = 3.14159265358979323846;
sem_t sem;

void* monte_carlo(void* param) {
    int i, batch, seed;
    int inside_points = 0;
    seed = *(int*)param;
    batch = TOTAL_POINTS / NUM_THREADS;

    for (i = 0; i < batch; i++) {
        float x, y, z;
            
	x = PI * (rand_r(&seed) / ((double)RAND_MAX));
        y = 1.0 * (rand_r(&seed) /((double) RAND_MAX));
        z = PI * (rand_r(&seed) /((double) RAND_MAX));

	if ((x <= PI) && (y <= sin(x)) && (z <= x*y)){
            inside_points++;;
        }
    }
    sem_wait(&sem);
    inside_counter += inside_points;
    sem_post(&sem);
    // printf("Points inside %i/%i=%.3f\n", inside_points, batch, 1.0*inside_points/batch);
}

int main(int argc, char *argv[]){
    NUM_THREADS = atoi(argv[1]);
    pthread_t *pthr = malloc(NUM_THREADS * sizeof(pthread_t));
    int *buf = malloc(NUM_THREADS * sizeof(int));
    int i, err;
    struct timespec begin, end;
    double elapsed;

    clock_gettime(CLOCK_REALTIME, &begin);
    sem_init(&sem, 0, 1);
    for(i = 0; i < NUM_THREADS; i++){
        buf[i] = i;
        err = pthread_create(&pthr[i], NULL, monte_carlo, (void*)&buf[i]);
        if (err) 
        	printf("ERROR; return code from pthread_create() is %d \n", err);
    }

    for(i = 0; i < NUM_THREADS; i++){
        err = pthread_join(pthr[i], NULL);
        if (err)
        	printf("ERROR; return code from pthread_join() is %d \n", err);
    }
    sem_destroy(&sem);
    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = end.tv_sec - begin.tv_sec;
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    float volume;
    volume = PI * 1.0 * PI * inside_counter / TOTAL_POINTS;

    //FILE * fp;
    //fp = fopen("plot_data.txt", "a");
    //fprintf(fp, "%d, %f\n", NUM_THREADS, elapsed);
    //fclose(fp);
    printf("Integral = %f, Time = %f sec., Processes = %d\n", volume, elapsed, NUM_THREADS);
    free(buf);
    free(pthr);
    return 0;
}

// gcc code.c -lpthread -lrt -lm
