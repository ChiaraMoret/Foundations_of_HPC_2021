#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define BLOCK_SIZE 50


void box(int M, float R, int x_min, int x_max, float x, float y, float z, int* x_left, int* x_right, int* y_down, int* y_up, int* z_back, int* z_front) {
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

int main (int argc, char ** argv)
{   	//check number of arguments
    int myRank, nProcs, start, stop;
 	MPI_Init( &argc, &argv );                  	
 	MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
 	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	if (argc>4 && myRank==0){
		printf("Too many input arguments\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (argc<3&& myRank==0){
		printf("At least 2 inputs arguments required: N grid number && R radius\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	MPI_Barrier(MPI_COMM_WORLD);


	// some useful constants: N grid number; R radius; M used to compute distance between nodes in the grid
    const int N = atoi(argv[1]);
    const float R = atof(argv[2]);
    const int M = N-1;

    double begin = MPI_Wtime();

    int local_Np, Np, local_numFloats, numFloats;
	float* particles = NULL;
    if (argc==3){
        
        int seed = time(NULL)+myRank; 
        srand(seed);
		//the root decides the number of particles generated
		if (myRank == 0){
			//generate random Np  between 10 and 1k
			local_Np = rand()% (1000 + 1 - 10) + 10;
			local_Np /= nProcs;
			Np = local_Np*nProcs;

			printf("No input file has been provided. %i random particles will be generated.\n",local_Np*nProcs);
			
		}
        

		MPI_Bcast(&local_Np, 1, MPI_INT, 0, MPI_COMM_WORLD);
        

		//all processors recive the number of particles they should generate
		Np = local_Np*nProcs;
		local_numFloats = local_Np*3;
		numFloats=Np*3;


		//generate particles
        particles = (float*)malloc(Np*3*sizeof(float));
		
		float* local_particles = (float*)malloc(local_numFloats*sizeof(float));
		
		for (int i = 0; i < local_numFloats; i++) {
			local_particles[i] = ((float)rand())/RAND_MAX;
		}

		// all processors receive all particles	
		MPI_Allgather(local_particles, local_numFloats, MPI_FLOAT, particles, local_numFloats, MPI_FLOAT, MPI_COMM_WORLD);	
		
        free(local_particles);

		// define the indexes of particles which concern the processor
		start = myRank * local_Np;
		stop = (myRank + 1) * local_Np;
		
		start *=3;
		stop *=3;

	}
    

    if (argc==4){
		if(myRank==0){
			//open binary file 
			FILE *bf;
			bf = fopen(argv[3],"rb");

			// check if file has been recieved correctly 
			if (bf == NULL) {       
				printf("Failed to open the file.\n");
				MPI_Abort(MPI_COMM_WORLD, 1);
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

		MPI_Bcast(&Np, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if (myRank!=0){
			numFloats = Np*3;
			particles = (float*)malloc(numFloats*sizeof(float));
		}
		MPI_Bcast(particles, Np*3, MPI_FLOAT, 0, MPI_COMM_WORLD);

		// start and stop are used to denote the part of data of each processor's concern 
        local_Np = Np/nProcs;
		int reminder = Np%nProcs;                                                         
        
		if (myRank < reminder){
            local_Np ++;                                                           
			start = myRank * local_Np;                                          
        	stop = start + local_Np;
        	} else {                                                                          
        	start = myRank * local_Np + reminder;                                     
			stop = start + local_Np;
		}
        
		start *= 3;
		stop *= 3;

	}


    float grid_coordinates[N];                        
    for (int i=0; i<N; i++){
		grid_coordinates[i] = (float)i/M;
	}

	//prepare for blocking
    int numBlocks= N/BLOCK_SIZE;
    int blockReminder = N%BLOCK_SIZE;
    if (blockReminder!=0){
        numBlocks ++;
    } 
    
    int blockSizes[numBlocks] ;
	FILE *densityf;

	
	if (myRank==0){
		// ***produce the output binary file for the densities***
		// delete previous "densities.bin" file, if it exists
		remove("densities.bin");

		// create "densities.bin" file 
		//FILE *densityf;
		densityf = fopen("densities.bin", "wb");

		// the first digit is the grid number
		fwrite(&N, 4, 1, densityf);
	}
	

	
    //create local_potentials array (we'll update inside the first block)
    float* local_potentials = (float*)malloc(local_Np*4*sizeof(float));
	
	//we can iterate over the blocks
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


        // density arrays, initialized to all zeroes
		float* densities = NULL;
		if (myRank==0){
       		densities = (float*)calloc(N*N * blockSize, sizeof(float));
		}

		float* local_densities = (float*)calloc(N*N*blockSize,sizeof(float));
    

    


    
		//Each processor takes care of some particles
		for (int i=start; i< stop; i+=3){
			
			
			// particle coordinates
			float x_p = particles[i];
			float y_p = particles[i+1];
			float z_p = particles[i+2];
			

		
			// candidate nodes coordinates for the particle
			
			int x_left, x_right, y_down, y_up, z_back, z_front;
			box(M, R, x_min, x_max, x_p, y_p, z_p, &x_left, &x_right, &y_down, &y_up, &z_back, &z_front);


			//***We can update the local densities***
			for (int x=x_left; x<=x_right; x++){
				for (int y=y_down; y<=y_up; y++){
					for (int z=z_back; z<=z_front; z++){			
						// compute (the square of) particle-node distance
						float rho = powf(grid_coordinates[x]-x_p,2.) + powf(grid_coordinates[y]-y_p,2.) + powf(grid_coordinates[z]-z_p,2.);
						// if the particle is in the radius of the candidate node, add 1 to the density for that node
						if (rho<=R*R){
							int local_x = x-x_min;
							local_densities[local_x * N * N + y * N + z]++;
						}
					}
				}	
			}
			


		
			//***we can update the local potentials***
			//***we only update them in the first block***

			if (blockNum == 0){

				// find index of the paticle in the local_potentials array
				int n_p = ((i-start)/3) * 4;
				
				// fisrt  save the particle's coordinates
				local_potentials[n_p]= x_p;
				local_potentials[n_p + 1]= y_p;
				local_potentials[n_p + 2]= z_p;

				// compute the potential
				float V = 0;

				for (int j=0; j<numFloats; j+=3){
					if (j!=i){
						float r = powf(particles[j]-x_p,2.) + powf(particles[j+1]-y_p,2.) + powf(particles[j+2]-z_p,2.);
						V += 1/sqrt(r);
					}
				}
				local_potentials[n_p + 3] = V;
			}
		}
		
		//collect the local densities and sum them 
		MPI_Reduce(local_densities, densities, N*N*blockSize, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (myRank ==0){
			//the root updates the densities.bin file
			fwrite(densities, 4, N*N*blockSize, densityf);
		}

	//free memory for next block iteration
	free(local_densities);
	free(densities);
	}
	
	
	
    //***create file of potentials***
    
	// local_potentials have different lenghts
	// we need to use MPI_GatherV

	//lenght and displacement for the processor
    int p_floats=local_Np*4;
    int shift = (start/3)*4;

	//array of global_potentials, array of lenghts and array of displacements
    float* potentials = (float*)malloc(Np*4*sizeof(float));
    int* rcv_size = (int*)malloc( nProcs * sizeof(int)) ;
    int* rcv_shift = (int*)malloc( nProcs * sizeof(int)) ;

	//populate the arrays of lenghts and displacement
    MPI_Gather(&p_floats, 1, MPI_INT, rcv_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gather(&shift, 1, MPI_INT, rcv_shift, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
   
	// Concatenate the local potentials
    MPI_Gatherv(local_potentials, p_floats, MPI_FLOAT, potentials, rcv_size, rcv_shift, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if (myRank==0){
		//close the density file
		fclose(densityf);
	
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
    }
        
	//timing
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    double time = begin - end;
    

    if (myRank == 0) { /* use time on master node */
        printf("Runtime = %10.8f\n", end-start);
    }

    
    free(rcv_shift);
    free(rcv_size);
    free(potentials);
    free(local_potentials);
    free(particles);
    MPI_Finalize();
}