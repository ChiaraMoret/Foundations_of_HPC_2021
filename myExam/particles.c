#include <stdio.h>
#include <stdlib.h>
#include <time.h>

float* generate_particles(int N){

    int seed = time(NULL); // Random number generator seed (based on current time)
    srand(seed);

    //generate particles
    int Ndim = N*3;
    float* particles =(float*)malloc(Ndim*sizeof(float));
    int i;
    for (i = 0; i < Ndim; i++) {
        particles[i] = ((float)rand())/RAND_MAX;
        }
    return particles;
}

void main(int argc, char **argv)
{  
    //chech number of arguments
    if (argc!=2){
		printf("One argument required: number of particles to be generated\n");
		return -1;	
	}
    int Np = atoi(argv[1]);

     // if there is already a file, delete it
    remove("particles.bin");

    //open binary file 
    FILE *bf;
    bf = fopen("particles.bin", "wb");

    // the first digit is the number of particles 
    fwrite(&Np, 4, 1, bf);

    //generate particles
    float* particles = generate_particles(Np);

    //write paricles in file
    fwrite(particles, 4, Np*3 , bf);

    fclose(bf);
    free(particles);

}