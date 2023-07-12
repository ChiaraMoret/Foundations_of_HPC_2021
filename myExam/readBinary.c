#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv)
{   
    if (argc!=2){
        printf("One argument required: the bianary file name");
        return -1;
    }

    // Open the binary file for reading in binary mode
    FILE *bf;
    bf = fopen(argv[1],"rb");
    if (bf == NULL) {
        printf("Failed to open the file.\n");
        return -1;
    }

    // Determine the size of the file
    fseek(bf, 0, SEEK_END);
    long fileSize = ftell(bf);
    fseek(bf, 0, SEEK_SET);

    // Calculate the number of floats in the file
    size_t numFloats = (fileSize - sizeof(int)) / sizeof(float);

    // Read the int 
    int points;
    fread(&points, sizeof(int), 1, bf);
    //check the sizes match

    if (pow(points,3) != numFloats && !strcmp(argv[1],"densities.bin")){
        printf("I expected %i density points, instead there were %zu", points*points*points, numFloats);
        return -1;
    }
    else if (points*4 != numFloats && !strcmp(argv[1],"potentials.bin")){
        printf("I expected %i potential points, instead there were %zu", points, numFloats/4);
        return -1;
    }
    else if (points*3 != numFloats ){
        printf("I expected %i particle points, instead there were %zu", points, numFloats/3);
        return -1;
    }
    
    //Read the floats into an array
    float* floats = (float*)malloc(numFloats * sizeof(float));
    fseek(bf, sizeof(int), SEEK_SET);
    fread(floats, sizeof(float), numFloats, bf);

    // Close the file
    fclose(bf);

    // Print the # of points and print the floats

    // case 1: particles
    if (!strcmp(argv[1],"particles.bin")){
        printf("There are %i particles: \n", points); 
        int j=0;
        for (size_t i = 0; i < numFloats-3; i+=3) {
            printf("[%f,%f,%f], ", floats[i],floats[i+1],floats[i+2]);
            j++;
            if (j%3==0){
                printf("\n");
            }
            
        }
    
    printf("[%f,%f,%f]",floats[numFloats-3],floats[numFloats-2],floats[numFloats-1]);
    printf("\n");
    }

    //case 2: densities.bin
    if (!strcmp(argv[1],"densities.bin")){
        printf("Grid number is %i, for a total of %i nodes. Here are the densities: \n ", points,points*points*points);
        printf("[%f ", floats[0]);
        int i=1;
        while(i < numFloats) {
            printf(" %f ", floats[i]);
            i++;
            if (i%(points*points)==0){
                printf("]\n");
                if(i!=(numFloats)){
                printf("[");
            }
            }
            if (i%points==0 && i%(points*points)!=0){
                printf("\n");
            }
            
        }
    }
    // case 3: potentials.bin
    if (!strcmp(argv[1],"potentials.bin")){
        
        printf("There are %i particles. Here are their coordinates and their potentials: \n", points);
        int j=0;
        for (int i = 0; i < numFloats-4; i+=4) {
                printf("[%f,%f,%f; %f ], ", floats[i],floats[i+1],floats[i+2],floats[i+3]);
                j++;
                if (j%3==0){
                    printf("\n");
                }
            }
        printf("[%f,%f,%f; %f ] \n", floats[numFloats-4],floats[numFloats-3],floats[numFloats-2],floats[numFloats-1]);
    }


    // Free the dynamically allocated memory
    free(floats);
}