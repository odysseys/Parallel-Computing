/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         cuda_kmeans.cu  (CUDA version)                            */
/*   Description:  Implementation of simple k-means clustering algorithm     */
/*                 This program takes an array of N data objects, each with  */
/*                 M coordinates and performs a k-means clustering given a   */
/*                 user-provided value of the number of clusters (K). The    */
/*                 clustering results are saved in 2 arrays:                 */
/*                 1. a returned array of size [K][N] indicating the center  */
/*                    coordinates of K clusters                              */
/*                 2. membership[N] stores the cluster center ids, each      */
/*                    corresponding to the cluster a data object is assigned */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department, Northwestern University                        */
/*            email: wkliao@ece.northwestern.edu                             */
/*   Copyright, 2005, Wei-keng Liao                                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Copyright (c) 2005 Wei-keng Liao
// Copyright (c) 2011 Serban Giuroiu
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include "kmeans.h"

static inline int nextPowerOfTwo(int n) {
    n--;

    n = n >>  1 | n;
    n = n >>  2 | n;
    n = n >>  4 | n;
    n = n >>  8 | n;
    n = n >> 16 | n;
//  n = n >> 32 | n;    //  For 64-bit ints

    return ++n;
}

__inline static
float euclid_dist_2_seq(int    numdims,  /* no. dimensions */
                    float *coord1,   /* [numdims] */
                    float *coord2)   /* [numdims] */
{
    int i;
    float ans=0.0;

    #pragma unroll 2
    for (i=0; i<numdims; i++)
        ans += (coord1[i]-coord2[i]) * (coord1[i]-coord2[i]);

    return(ans);
}

__inline static
int find_nearest_cluster_seq(int     numClusters, /* no. clusters */
                         int     numCoords,   /* no. coordinates */
                         float  *object,      /* [numCoords] */
                         float **clusters)    /* [numClusters][numCoords] */
{
    int   index, i;
    float dist, min_dist;

    /* find the cluster id that has min distance to object */
    index    = 0;
    min_dist = euclid_dist_2_seq(numCoords, object, clusters[0]);

    for (i=1; i<numClusters; i++) {
        dist = euclid_dist_2_seq(numCoords, object, clusters[i]);
        /* no need square root */
        if (dist < min_dist) { /* find the min and its array index */
            min_dist = dist;
            index    = i;
        }
    }
    return(index);
}

/*----< euclid_dist_2() >----------------------------------------------------*/
/* square of Euclid distance between two multi-dimensional points            */
__host__ __device__ inline static
float euclid_dist_2(int    numCoords,
                    int    numObjs,
                    int    numClusters,
                    float *objects,     // [numCoords][numObjs]
                    float *clusters,    // [numCoords][numClusters]
                    int    objectId,
                    int    clusterId)
{
    
    float ans=0.0;

    int i;

    #pragma unroll 4
   
    for (i = 0; i < numCoords; i++) {
        ans += (objects[numObjs * i + objectId] - clusters[numClusters * i + clusterId]) *
               (objects[numObjs * i + objectId] - clusters[numClusters * i + clusterId]);
    }

    return(ans);
}

/*----< find_nearest_cluster() >---------------------------------------------*/
__global__ static
void find_nearest_cluster(int numCoords,
                          int numObjs,
                          int numClusters,
                          float *objects,           //  [numCoords][numObjs]
                          float *deviceClusters,    //  [numCoords][numClusters]
                          int *membership,          //  [numObjs]
                          int *intermediates)
{
    extern __shared__ char sharedMemory[];

    //  The type chosen for membershipChanged must be large enough to support
    //  reductions! There are blockDim.x elements, one for each thread in the
    //  block.
    unsigned char *membershipChanged = (unsigned char *)sharedMemory;
    float *clusters = (float *)(sharedMemory + blockDim.x);

    membershipChanged[threadIdx.x] = 0;

    //  BEWARE: We can overrun our shared memory here if there are too many
    //  clusters or too many coordinates!
    for (int i = threadIdx.x; i < numClusters; i += blockDim.x) {
        for (int j = 0; j < numCoords; j++) {
            clusters[numClusters * j + i] = deviceClusters[numClusters * j + i];
        }
    }
    __syncthreads();

    int objectId = blockDim.x * blockIdx.x + threadIdx.x;

    if (objectId < numObjs) {
        int   index, i;
        float dist, min_dist;

        /* find the cluster id that has min distance to object */
        index    = 0;
        min_dist = euclid_dist_2(numCoords, numObjs, numClusters,
                                 objects, clusters, objectId, 0);


        #pragma unroll 4
        for (i=1; i<numClusters; i++) {
            dist = euclid_dist_2(numCoords, numObjs, numClusters,
                                 objects, clusters, objectId, i);
            /* no need square root */
            if (dist < min_dist) { /* find the min and its array index */
                min_dist = dist;
                index    = i;
            }
        }

        if (membership[objectId] != index) {
            membershipChanged[threadIdx.x] = 1;
            //membership[objectId] = index;
        }

        /* assign the membership to object objectId */
        membership[objectId] = index;

        __syncthreads();    //  For membershipChanged[]

        //  blockDim.x *must* be a power of two!
        #pragma unroll 4
        for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
            if (threadIdx.x < s) {
                membershipChanged[threadIdx.x] +=
                    membershipChanged[threadIdx.x + s];
            }
            __syncthreads();
        }

        if (threadIdx.x == 0) {
            intermediates[blockIdx.x] = membershipChanged[0];
        }
    }
}

__global__ static
void compute_delta(int *deviceIntermediates,
                   int numIntermediates,    //  The actual number of intermediates
                   int numIntermediates2)   //  The next power of two
{
    //  The number of elements in this array should be equal to
    //  numIntermediates2, the number of threads launched. It *must* be a power
    //  of two!
    extern __shared__ unsigned int intermediates[];

    //  Copy global intermediate values into shared memory.
    intermediates[threadIdx.x] =
        (threadIdx.x < numIntermediates) ? deviceIntermediates[threadIdx.x] : 0;

    __syncthreads();

    //  numIntermediates2 *must* be a power of two!
    #pragma unroll 4
    for (unsigned int s = numIntermediates2 / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            intermediates[threadIdx.x] += intermediates[threadIdx.x + s];
        }
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        deviceIntermediates[0] = intermediates[0];
    }
}

//numCoords, numObjs, numClusters, deviceNewClusters, deviceNewClusterSizes, deviceMembership
__global__ static
void update_new_cluster(int numCoords,
                        int numObjs,
                        int numClusters,
                        float ** deviceNewClusters,
                        int *deviceNewClusterSize,
                        int *deviceMembership, 
                        float *deviceObjects) {

    int index;
    for(int i = 0; i<numObjs; i++) {
        index = deviceMembership[i];
        deviceNewClusterSize[index]++;
        for(int j=0;j<numCoords;j++) {
            deviceNewClusters[j][index] += deviceObjects[j * numObjs + i];
        }
    }

    for (int i=0; i<numClusters; i++) {
        for (int j=0; j<numCoords; j++) {
            if (deviceNewClusterSize[i] > 0) {
                deviceNewClusters[j][index] =  deviceNewClusters[j][index] / deviceNewClusterSize[i];
            }
        }
    }
}

__global__ static
void updateClusterSize_Sum(int numObjs,
                           int numCoords,
                           int numClusters,
                           int *deviceMembership,
                           float *deviceObjects,
                           int *deviceNewClusterSize,
                           float *deviceNewClusters)//<<<numClusterBlocks, numThreadsPerClusterBlock>>>
{
    if((blockIdx.x * blockDim.x + threadIdx.x) < numObjs)
    {
        int index = deviceMembership[blockIdx.x * blockDim.x + threadIdx.x];
        
        atomicAdd(&deviceNewClusterSize[index], 1);
        //must use exclusive thread process to prevent race condition
        for(int i = 0; i < numCoords; i++)
            atomicAdd(&deviceNewClusters[i * numClusters + index], deviceObjects[i * numObjs + (blockIdx.x * blockDim.x + threadIdx.x)]);
            //must use exclusive thread process to prevent race condition
    }

}

__global__ static
void updatedimClusters(int *deviceNewClusterSize,
                       float *deviceNewClusters,
                       float *deviceClusters)
{
    if(deviceNewClusterSize[blockIdx.x] > 0)//calcaulate new cluster centroids by each thread
        deviceClusters[threadIdx.x * gridDim.x + blockIdx.x] = deviceNewClusters[threadIdx.x * gridDim.x + blockIdx.x] / deviceNewClusterSize[blockIdx.x];
    deviceNewClusters[threadIdx.x * gridDim.x + blockIdx.x] = 0;
}


/*----< cuda_kmeans() >-------------------------------------------------------*/
//
//  ----------------------------------------
//  DATA LAYOUT
//
//  objects         [numObjs][numCoords]
//  clusters        [numClusters][numCoords]
//  dimObjects      [numCoords][numObjs]
//  dimClusters     [numCoords][numClusters]
//  newClusters     [numCoords][numClusters]
//  deviceObjects   [numCoords][numObjs]
//  deviceClusters  [numCoords][numClusters]
//  ----------------------------------------
//
/* return an array of cluster centers of size [numClusters][numCoords]       */
float** cuda_kmeans(float **objects,      /* in: [numObjs][numCoords] */
                   int     numCoords,    /* no. features */
                   int     numObjs,      /* no. objects */
                   int     numClusters,  /* no. clusters */
                   float   threshold,    /* % objects change membership */
                   int    *membership,   /* out: [numObjs] */
                   int    *loop_iterations)
{   
    // When numObjs > 10000, cuda version is more efficient than sequencial version
    if (numObjs > 5000) {
        int      i, j, loop=0;
        int     *newClusterSize; /* [numClusters]: no. objects assigned in each
                                    new cluster */
        int      delta;          /* % of objects change their clusters */
        float  **dimObjects;
        float  **clusters;       /* out: [numClusters][numCoords] */
        float  **dimClusters;
        float  **newClusters;    /* [numCoords][numClusters] */

        float *deviceObjects;
        float *deviceClusters;
        int *deviceMembership;
        int *deviceIntermediates;

        // self
        float *deviceNewClusters;
        int *deviceNewClusterSize;

        //  Copy objects given in [numObjs][numCoords] layout to new
        //  [numCoords][numObjs] layout
        malloc2D(dimObjects, numCoords, numObjs, float);

        /**
        ******************************************************************************************************************************
        */ 
        for (i = 0; i < numCoords; i++) {
            for (j = 0; j < numObjs; j++) {
                dimObjects[i][j] = objects[j][i];
            }
        }
        /**
        ******************************************************************************************************************************
        */

        /* pick first numClusters elements of objects[] as initial cluster centers*/
        malloc2D(dimClusters, numCoords, numClusters, float);
        /**
        ******************************************************************************************************************************
        */ 
        
        for (i = 0; i < numCoords; i++) {
            for (j = 0; j < numClusters; j++) {
                dimClusters[i][j] = dimObjects[i][j];
            }
        }
        /**
        ******************************************************************************************************************************
        */

        /* initialize membership[] */
        // loop unrolling
        #pragma unroll 4
        for (i=0; i<numObjs; i++) membership[i] = -1;

        /* need to initialize newClusterSize and newClusters[0] to all 0 */
        newClusterSize = (int*) calloc(numClusters, sizeof(int));
        assert(newClusterSize != NULL);

        malloc2D(newClusters, numCoords, numClusters, float);
        memset(newClusters[0], 0, numCoords * numClusters * sizeof(float));

        //  To support reduction, numThreadsPerClusterBlock *must* be a power of
        //  two, and it *must* be no larger than the number of bits that will
        //  fit into an unsigned char, the type used to keep track of membership
        //  changes in the kernel.
        const unsigned int numThreadsPerClusterBlock = 1024;
        // Original numThreadsPerClusterBlock is 128, it's not enough to compute numObjs > 131072 in test case 3 and 4
        // We just use maximum number of threads per block: 1024
        const unsigned int numClusterBlocks =
            (numObjs + numThreadsPerClusterBlock - 1) / numThreadsPerClusterBlock;
        const unsigned int clusterBlockSharedDataSize =
            numThreadsPerClusterBlock * sizeof(unsigned char) +
            numClusters * numCoords * sizeof(float);

        const unsigned int numReductionThreads =
            nextPowerOfTwo(numClusterBlocks);
        const unsigned int reductionBlockSharedDataSize =
            numReductionThreads * sizeof(unsigned int);

        checkCuda(cudaMalloc(&deviceObjects, numObjs*numCoords*sizeof(float)));
        checkCuda(cudaMalloc(&deviceClusters, numClusters*numCoords*sizeof(float)));
        checkCuda(cudaMalloc(&deviceMembership, numObjs*sizeof(int)));
        checkCuda(cudaMalloc(&deviceIntermediates, numReductionThreads*sizeof(unsigned int)));
        checkCuda(cudaMalloc(&deviceNewClusterSize, numClusters*sizeof(int)));
        cudaMemset(deviceNewClusterSize, 0, numClusters*sizeof(int));
        checkCuda(cudaMalloc(&deviceNewClusters, numClusters*numCoords*sizeof(float)));
        cudaMemset(deviceNewClusters, 0, numCoords * numClusters * sizeof(float));


        checkCuda(cudaMemcpy(deviceObjects, dimObjects[0],
                  numObjs*numCoords*sizeof(float), cudaMemcpyHostToDevice));
        checkCuda(cudaMemcpy(deviceMembership, membership,
                  numObjs*sizeof(int), cudaMemcpyHostToDevice));
        cudaMemcpy(deviceClusters, dimClusters[0],
                      numClusters*numCoords*sizeof(float), cudaMemcpyHostToDevice);

        do {
            find_nearest_cluster
                <<< numClusterBlocks, numThreadsPerClusterBlock, clusterBlockSharedDataSize >>>
                (numCoords, numObjs, numClusters,
                 deviceObjects, deviceClusters, deviceMembership, deviceIntermediates);

            cudaThreadSynchronize(); //checkLastCudaError();

            compute_delta <<< 1, numReductionThreads, reductionBlockSharedDataSize >>>
                (deviceIntermediates, numClusterBlocks, numReductionThreads);

            cudaThreadSynchronize(); //checkLastCudaError();

            int d;
            checkCuda(cudaMemcpy(&d, deviceIntermediates,
                      sizeof(int), cudaMemcpyDeviceToHost));
            delta = d;

            checkCuda(cudaMemcpy(membership, deviceMembership,
                      numObjs*sizeof(int), cudaMemcpyDeviceToHost));

            // ******************************************************************
            // checkCuda(cudaMalloc(&deviceNewClusters, numObjs*numCoords*sizeof(float)));
            // checkCuda(cudaMalloc(&deviceNewClusterSizes, numClusters*sizeof(int)));
            // update_new_cluster
            //     <<<numClusterBlocks, numThreadsPerClusterBlock>>>(numCoords, numObjs, numClusters, deviceNewClusters, 
            //         deviceNewClusterSizes, deviceMembership, deviceObjects);
            // cudaThreadSynchronize();

            // checkCuda(cudaMemcpy(dimClusters, deviceNewClusters, numObjs*numCoords*sizeof(float), cudaMemcpyDeviceToHost));

            updateClusterSize_Sum <<<numClusterBlocks, numThreadsPerClusterBlock>>>
                (numObjs, numCoords, numClusters, deviceMembership, deviceObjects, deviceNewClusterSize, deviceNewClusters);

            cudaThreadSynchronize();

            //  TODO: Flip the nesting order
            //  TODO: Change layout of newClusters to [numClusters][numCoords]
            /* average the sum and replace old cluster centers with newClusters */
            
            updatedimClusters <<<numClusters, numCoords>>>
                (deviceNewClusterSize, deviceNewClusters, deviceClusters);

            cudaThreadSynchronize();

            cudaMemset(deviceNewClusterSize, 0, numClusters*sizeof(int));
            //reset NewClusterSize to 0 for next iteration

        } while ((loop <= 1 || delta > (threshold * numObjs)) && loop++ < 500 );
        // must do at least 1 time

        *loop_iterations = loop + 1;

        /* allocate a 2D space for returning variable clusters[] (coordinates
           of cluster centers) */
        cudaMemcpy(dimClusters[0], deviceClusters,
                      numClusters*numCoords*sizeof(float), cudaMemcpyDeviceToHost);
        //copy back dimCluster for output
        malloc2D(clusters, numClusters, numCoords, float);
        

        #pragma unroll 4
        for (i = 0; i < numClusters; i++) {
            for (j = 0; j < numCoords; j++) {
                clusters[i][j] = dimClusters[j][i];
            }
        }

        checkCuda(cudaFree(deviceObjects));
        checkCuda(cudaFree(deviceClusters));
        checkCuda(cudaFree(deviceMembership));
        checkCuda(cudaFree(deviceIntermediates));
        checkCuda(cudaFree(deviceNewClusters));
        checkCuda(cudaFree(deviceNewClusterSize));

        free(dimObjects[0]);
        free(dimObjects);
        free(dimClusters[0]);
        free(dimClusters);
        free(newClusters[0]);
        free(newClusters);
        free(newClusterSize);

        return clusters;
    } else { // When numObjs <= 5000, sequencial version is more efficient than cuda version
        int      i, j, index, loop=0;
        int     *newClusterSize; /* [numClusters]: no. objects assigned in each
                                    new cluster */
        int      delta;          /* % of objects change their clusters */
        float  **clusters;       /* out: [numClusters][numCoords] */
        float  **newClusters;    /* [numClusters][numCoords] */

        /* allocate a 2D space for returning variable clusters[] (coordinates
           of cluster centers) */
        clusters    = (float**) malloc(numClusters *             sizeof(float*));
        assert(clusters != NULL);
        clusters[0] = (float*)  malloc(numClusters * numCoords * sizeof(float));
        assert(clusters[0] != NULL);
        for (i=1; i<numClusters; i++)
            clusters[i] = clusters[i-1] + numCoords;

        /* pick first numClusters elements of objects[] as initial cluster centers*/
        for (i=0; i<numClusters; i++)
            for (j=0; j<numCoords; j++)
                clusters[i][j] = objects[i][j];

        /* initialize membership[] */
        for (i=0; i<numObjs; i++) membership[i] = -1;

        /* need to initialize newClusterSize and newClusters[0] to all 0 */
        newClusterSize = (int*) calloc(numClusters, sizeof(int));
        assert(newClusterSize != NULL);

        newClusters    = (float**) malloc(numClusters *            sizeof(float*));
        assert(newClusters != NULL);
        newClusters[0] = (float*)  calloc(numClusters * numCoords, sizeof(float));
        assert(newClusters[0] != NULL);
        for (i=1; i<numClusters; i++)
            newClusters[i] = newClusters[i-1] + numCoords;

        do {
            delta = 0.0;
            for (i=0; i<numObjs; i++) {
                /* find the array index of nestest cluster center */
                index = find_nearest_cluster_seq(numClusters, numCoords, objects[i],
                                             clusters);

                /* if membership changes, increase delta by 1 */
                if (membership[i] != index) delta += 1.0;

                /* assign the membership to object i */
                membership[i] = index;

                /* update new cluster centers : sum of objects located within */
                newClusterSize[index]++;
                for (j=0; j<numCoords; j++)
                    newClusters[index][j] += objects[i][j];
            }

            /* average the sum and replace old cluster centers with newClusters */
            for (i=0; i<numClusters; i++) {
                for (j=0; j<numCoords; j++) {
                    if (newClusterSize[i] > 0)
                        clusters[i][j] = newClusters[i][j] / newClusterSize[i];
                    newClusters[i][j] = 0.0;   /* set back to 0 */
                }
                newClusterSize[i] = 0;   /* set back to 0 */
            }
        } while (delta > (threshold * numObjs) && loop++ < 500);

        *loop_iterations = loop + 1;

        free(newClusters[0]);
        free(newClusters);
        free(newClusterSize);

        return clusters;
    }
}

