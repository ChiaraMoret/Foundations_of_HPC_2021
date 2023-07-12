#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#define BLOCK_SIZE 100
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

double getCPUTime() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

float* generate_particles(int N){

    //generate particles
    int Ndim = N*3;
    float* particles =(float*)malloc(Ndim*sizeof(float));
    int seed = 100;
    unsigned short xsubi[3] = {seed, 0, 0};
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        unsigned short private_xsubi[3] = {xsubi[0] + id, xsubi[1], xsubi[2]};

        #pragma omp parallel for
        for (int i = 0; i < Ndim; i++) {
            particles[i] = (float)erand48(private_xsubi);
            }

    }

    
    
    
    return particles;
}

void box(int M, float R, int x_min, int x_max, float x, float y, float z, int* x_left, int* x_right, int* y_down, int* y_up, int* z_back, int* z_front) {
    //this is the box function 
    int x_l = (int)ceil((x - R) * M);
    int x_r = (int)floor((x + R) * M);
    
    int y_d = (int)ceil((y - R) * M);
    int y_u = (int)floor((y + R) * M);
    
    int z_b = (int)ceil((z - R) * M);
    int z_f = (int)floor((z + R) * M);

    *x_left = max(x_min, x_l);
    *x_right = min(x_r, x_max);
    *y_down = max(0, y_d);
    *y_up = min(y_u, M);
    *z_back = max(0, z_b);
    *z_front = min(z_f, M);
    

}
int main(int argc, char** argv)
{   //check number of arguments
    if (argc>4){
        printf("Too many input arguments\n");
        return -1;
	}
    if (argc<3){
        printf("At least 2 inputs arguments required: N grid number && R radius\n");
        return -1;
	}
    double start = getCPUTime();

    // some useful constants: N grid number; R radius; M used to compute distance between nodes in the grid
    const int N = atoi(argv[1]);
    const float R = atof(argv[2]);
    const int M = N-1;
	
    int Np, numFloats;
    float* particles = NULL;
    
	if (argc==3){
        srand(time(NULL));

       	//generate random Np  between 10 and 1k
       	Np = rand()% ((1000 + 1 - 10) + 10);
		numFloats =Np*3;

        printf("No input file has been provided. %i random particles will be generated.\n",Np);
        //generate particles
        particles = generate_particles(Np);
		
	}
	if (argc==4){
		//open binary file 
		FILE *bf;
		bf = fopen(argv[3],"rb");

		// check if file has been received correctly 
		if (bf == NULL) {       
			printf("Failed to open the file.\n");
			return -1;
		}
		// get number of particles
		fread(&Np, sizeof(int), 1, bf);

		numFloats = Np*3;

		// get array of particles
		particles = (float*)malloc(numFloats*sizeof(float));

        fseek(bf, sizeof(int), SEEK_SET);
        fread(particles, sizeof(float), numFloats, bf);

        // Close the file
        fclose(bf);
	}

    // compute array of node coordinates
    float grid_coordinates[N];

    for (int i = 0; i < N; i++) {
        grid_coordinates[i] = (float)i / M;
    }


    //prepare for blocking
    int numBlocks= N/BLOCK_SIZE;
    int blockReminder = N%BLOCK_SIZE;
    if (blockReminder!=0){
        numBlocks ++;
    } 
    
    int blockSizes[numBlocks] ;

     // ***produce the output binary file for the densities***
     // delete previous "densities.bin" file, if it exists
    remove("densities.bin");

    // create "densities.bin" file 
    FILE *densityf;
    densityf = fopen("densities.bin", "wb");

    // the first digit is the grid number
    fwrite(&N, 4, 1, densityf);

    //iterate over the blocks
    for(int blockNum = 0; blockNum < numBlocks; blockNum++){

        //We need the size of the block
        int blockSize = BLOCK_SIZE;
        if (blockNum == numBlocks -1 && blockReminder!=0){
            blockSize = blockReminder;
        }
        blockSizes[blockNum] = blockSize;

        //We need the indexes of coordinate x the block will deal with
        int x_min = (blockNum) * BLOCK_SIZE ;
        int x_max = x_min + blockSize -1;

        // density array, initialized to all zeroes
        float* densities = (float*)calloc(N * N * blockSize, sizeof(float));

        

        // ***update the density values***
        
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < numFloats; i += 3) {
            // particle coordinates
            float x_p = particles[i];
            float y_p = particles[i+1];
            float z_p = particles[i+2];

            // candidate nodes coordinates for the particle
            int x_left, x_right, y_down, y_up, z_back, z_front;
            box(M, R, x_min, x_max, x_p, y_p, z_p,  &x_left, &x_right, &y_down, &y_up, &z_back, &z_front);
            
            // update densities
            for (int x = x_left; x <= x_right; x++) {
                for (int y = y_down; y <= y_up; y++) {
                    for (int z = z_back; z <= z_front; z++) {
                        // compute (the square of) particle-node distance
                        float rho = powf(grid_coordinates[x] - x_p, 2) + powf(grid_coordinates[y] - y_p, 2) + powf(grid_coordinates[z] - z_p, 2);
                        // if the particle is in the radius of the candidate node, add 1 to the density for that node
                        if (rho <= R * R) {
                            int local_x = x-x_min;
                            #pragma omp atomic update
                            densities[local_x * N * N + y * N + z]++;
                            
                        }
                    }
                }
            }
            

        }

        // write densities in file
        fwrite(densities, 4, N*N*blockSize, densityf);
        //free the memory for the next block
        free(densities);
    
    }

    //close the densities.bin file
    fclose(densityf);

    // potentials array
    float* potentials = (float*)malloc(Np * 4 * sizeof(float));


    for (int i = 0; i < numFloats; i += 3) {
        // particle coordinates
        float x_p = particles[i];
        float y_p = particles[i+1];
        float z_p = particles[i+2];
        
        
        // find index of the paticle in the potentials array
        int n_p = (i/3) * 4;
        // save the particle's coordinates in the potentials array
        potentials[n_p]= x_p;
        potentials[n_p + 1]= y_p;
        potentials[n_p + 2]= z_p;
        
        // compute the potential using reduction
        float V = 0;
        #pragma omp parallel for reduction(+:V)
        for (int j=0; j<numFloats; j+=3){
            if (j!=i){
                float r = powf(particles[j]-x_p,2.) + powf(particles[j+1]-y_p,2.) + powf(particles[j+2]-z_p,2.);
                V += 1/sqrt(r);
            }
        }

        potentials[n_p + 3] = V;

    }
     //I no longer need the particle coordinates; I can free them
    free(particles);

    
	

	


	// ***produce the output binary file for the potentials***
	// delete previous "potentials.bin" file, if it exists
	remove("potentials.bin");

    //open binary file 
    FILE *potentialsf;
    potentialsf = fopen("potentials.bin", "wb");

    // the first digit is the number of particles
    fwrite(&Np, 4, 1, potentialsf);

    //write potentials in file
    fwrite(potentials, 4,4*Np, potentialsf);

    fclose(potentialsf);


 	
	free(potentials);

    // get elapsed time
	double end  = getCPUTime();
    

	printf("Elapsed time: %f seconds.\n",end - start);
    
    return 0;
}

