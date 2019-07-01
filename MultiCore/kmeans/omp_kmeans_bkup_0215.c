/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         kmeans_clustering.c  (OpenMP version)                     */
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

#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include "kmeans.h"
#include <string.h>
#include <x86intrin.h>

/*----< euclid_dist_2() >----------------------------------------------------*/
/* square of Euclid distance between two multi-dimensional points            */
__inline static
float euclid_dist_2(int    numdims,  /* no. dimensions */
                    float *coord1,   /* [numdims] */
                    float *coord2)   /* [numdims] */
{
    int i;
    float diff1, diff2, diff3, diff4;
    float dist1, dist2, dist3, dist4 = 0;
    
/*
    float ans=0.0;
     //#pragma omp parallel for reduction(+: ans)
    for (i=0; i<numdims; i++)
        ans += (coord1[i]-coord2[i]) * (coord1[i]-coord2[i]);
    return ans;
*/

    //SIMD
    __m128 vcoord1;
    __m128 vcoord2;
    __m128 vmul = _mm_setzero_ps();
    __m128 vsub;
    float accum_ans[4];

    for(i = 0; i < numdims; i = i + 4)
    {
/*   
        if((i + 4) > numdims)
        {
            if(numdims % 4 == 3)
            {
                vcoord1 = _mm_setr_ps(coord1[i], coord1[i+1], coord1[i+2], 0.0f);
                vcoord2 = _mm_setr_ps(coord2[i], coord2[i+1], coord2[i+2], 0.0f);
            }
            else if(numdims % 4 == 2)
            {
                vcoord1 = _mm_setr_ps(coord1[i], coord1[i+1], 0.0f, 0.0f);
                vcoord2 = _mm_setr_ps(coord2[i], coord2[i+1], 0.0f, 0.0f);
            }
            else if(numdims % 4 == 1)
            {
                vcoord1 = _mm_setr_ps(coord1[i], 0.0f, 0.0f, 0.0f);
                vcoord2 = _mm_setr_ps(coord2[i], 0.0f, 0.0f, 0.0f);
            }
        }
        else*/
        {
            vcoord1 = _mm_load_ps(&coord1[i]);
            vcoord2 = _mm_load_ps(&coord2[i]);
            //vcoord1 = _mm_setr_ps(coord1[i], coord1[i+1], coord1[i+2], coord1[i+3]);
            //vcoord2 = _mm_setr_ps(coord2[i], coord2[i+1], coord2[i+2], coord2[i+3]);
        }
        vsub = _mm_sub_ps(vcoord1, vcoord2);
        vmul = _mm_add_ps(vmul, _mm_mul_ps(vsub, vsub));
    }
    _mm_store_ps(&accum_ans[0], vmul);
    return (accum_ans[0] + accum_ans[1] + accum_ans[2] + accum_ans[3]);

/*
    //SIMD-like
    //#pragma omp parallel reduction(+: dist1, dist2, dist3, dist4)
    for (i = 0; i < numdims; i = i + 4)
    {
      diff1 = coord1[i    ] - coord2[i    ];
      diff2 = coord1[i + 1] - coord2[i + 1];
      diff3 = coord1[i + 2] - coord2[i + 2];
      diff4 = coord1[i + 3] - coord2[i + 3];
      dist1 += diff1 * diff1;
      dist2 += diff2 * diff2;
      dist3 += diff3 * diff3;
      dist4 += diff4 * diff4;
    }

    return(dist1 + dist2 + dist3 + dist4);
*/
}

/*----< find_nearest_cluster() >---------------------------------------------*/
__inline static
int find_nearest_cluster(int     numClusters, /* no. clusters */
                         int     numCoords,   /* no. coordinates */
                         float  *object,      /* [numCoords] */
                         float **clusters)    /* [numClusters][numCoords] */
{
    int   index, i;
    float dist, min_dist;

    /* find the cluster id that has min distance to object */
    index    = 0;
    min_dist = euclid_dist_2(numCoords, object, clusters[0]);

    for (i=1; i<numClusters; i++) {
        dist = euclid_dist_2(numCoords, object, clusters[i]);
        /* no need square root */
        if (dist < min_dist) { /* find the min and its array index */
            min_dist = dist;
            index    = i;
        }
    }
    return(index);
}


/*----< get next cerntroid index at initial >------------------------------------*/
__inline static
int get_next_centroid_index(int numObjs, float *nearest_dist)
{
  int i;
  float dist_sum_sq = 0;
  float *prob_each;
  float *prob_sum;
  double rand_num = (double) rand() / (RAND_MAX + 1.0);
  
  prob_each = (float*) malloc(numObjs * sizeof(float));
  prob_sum  = (float*) calloc(numObjs, sizeof(float));
  
  #pragma omp parallel for reduction(+: dist_sum_sq)
  for (i = 0; i < numObjs; i++)
    dist_sum_sq += nearest_dist[i];
  
  #pragma omp parallel for
  for (i = 0; i < numObjs; i++)
    prob_each[i] = nearest_dist[i] / dist_sum_sq;
  
  prob_sum[0] = prob_each[0];
  for (i = 1; i < numObjs; i++)
    prob_sum[i] = prob_sum[i-1] + prob_each[i];
  
  //printf("RAND_MAX = %d", RAND_MAX);
  for (i = 0; i < numObjs; i++)
    if (rand_num < prob_sum[i])
    {//printf("next centroid idx is: %d", i);
      free(prob_each);
      free(prob_sum);
      return i;
    }

}

/*----< kmeans_clustering() >------------------------------------------------*/
/* return an array of cluster centers of size [numClusters][numCoords]       */
float** omp_kmeans(int     is_perform_atomic, /* in: */
                   float **objects,           /* in: [numObjs][numCoords] */
                   int     numCoords,         /* no. coordinates */
                   int     numObjs,           /* no. objects */
                   int     numClusters,       /* no. clusters */
                   float   threshold,         /* % objects change membership */
                   int    *membership)        /* out: [numObjs] */
{

    int      i, j, k, index, loop=0;
    int     *newClusterSize; /* [numClusters]: no. objects assigned in each
                                new cluster */
    float    delta;          /* % of objects change their clusters */
    float  **clusters;       /* out: [numClusters][numCoords] */
    float  **newClusters;    /* [numClusters][numCoords] */
    double   timing;
    double timing2;
    int      nthreads;             /* no. threads */
    int    **local_newClusterSize; /* [nthreads][numClusters] */
    float ***local_newClusters;    /* [nthreads][numClusters][numCoords] */
    //srand(time(NULL));
    //srand(123456);
    int idx, next_centroid;
    float threshold_mul = threshold * numObjs;
    float newClusterSize_inverse;
    int reallo_numCoords = numCoords;;
    float **reallo_objs;


    //reallocate the objects if numCoords is not multiple of 4
    if (numCoords % 4 != 0)
    {
        reallo_numCoords = numCoords - (numCoords % 4) + 4;
        reallo_objs =  (float**) malloc(numObjs * sizeof(float*));
        assert(reallo_objs != NULL);
        reallo_objs[0] = (float*) calloc(numObjs * reallo_numCoords, sizeof(float));
        assert(reallo_objs[0] != NULL);
        memcpy(reallo_objs[0], objects[0], numCoords * sizeof(float));
        for (i=1; i<numObjs; i++)
        {
            reallo_objs[i] = reallo_objs[i-1] + reallo_numCoords;
            memcpy(reallo_objs[i], objects[i], numCoords * sizeof(float));
        }
    
        objects = reallo_objs;
        numCoords = reallo_numCoords;
    }
    

    nthreads = omp_get_max_threads();

    /* allocate a 2D space for returning variable clusters[] (coordinates
       of cluster centers) */
    clusters    = (float**) malloc(numClusters *             sizeof(float*));
    assert(clusters != NULL);
    clusters[0] = (float*)  malloc(numClusters * numCoords * sizeof(float));
    assert(clusters[0] != NULL);
    for (i=1; i<numClusters; i++)
        clusters[i] = clusters[i-1] + numCoords;

    /* pick first numClusters elements of objects[] as initial cluster centers*/
    /* Using kmeans++ to get initial cluster centroids*/

/*
    float *nearest_dist;
    nearest_dist = (float*) malloc(numObjs * sizeof(float));
    #pragma omp parallel for private(j)
    for (j = 0; j < reallo_numCoords; j++)
    {
        clusters[0][j] = objects[0][j];
    }

    for(i = 1; i < numClusters; i++)
    {
        #pragma omp parallel for private(idx)
        for(k = 0; k < numObjs; k++)
        {
            idx = find_nearest_cluster(i, reallo_numCoords, objects[k], clusters);
            nearest_dist[k] = euclid_dist_2(reallo_numCoords, objects[k], clusters[idx]);
        }
        next_centroid = get_next_centroid_index(numObjs, nearest_dist);
        //next_centroid = 1;
        #pragma omp parallel for
        for(j = 0; j < reallo_numCoords; j++)
        {
            clusters[i][j] = objects[next_centroid][j];
        }

    }
*/


    
    for (i=0; i<numClusters; i++)
        for (j=0; j<numCoords; j++)
            clusters[i][j] = objects[i][j];
    

    /* initialize membership[] */
    #pragma omp parallel for
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

    if (!is_perform_atomic) {
        /* each thread calculates new centers using a private space,
           ithen thread 0 does an array reduction on them. This approach
           should be faster */
        local_newClusterSize    = (int**) malloc(nthreads * sizeof(int*));
        assert(local_newClusterSize != NULL);
        local_newClusterSize[0] = (int*)  calloc(nthreads*numClusters,
                                                 sizeof(int));
        assert(local_newClusterSize[0] != NULL);
        for (i=1; i<nthreads; i++)
            local_newClusterSize[i] = local_newClusterSize[i-1]+numClusters;

        /* local_newClusters is a 3D array */
        local_newClusters    =(float***)malloc(nthreads * sizeof(float**));
        assert(local_newClusters != NULL);
        local_newClusters[0] =(float**) malloc(nthreads * numClusters *
                                               sizeof(float*));
        assert(local_newClusters[0] != NULL);
        for (i=1; i<nthreads; i++)
            local_newClusters[i] = local_newClusters[i-1] + numClusters;
        for (i=0; i<nthreads; i++) {
            for (j=0; j<numClusters; j++) {
                local_newClusters[i][j] = (float*)calloc(numCoords,
                                                         sizeof(float));
                assert(local_newClusters[i][j] != NULL);
            }
        }
    }
timing2=omp_get_wtime();
    if (_debug) timing = omp_get_wtime();
    do {
        delta = 0.0;

        if (is_perform_atomic) {
            #pragma omp parallel for \
                    private(i,j,index) \
                    firstprivate(numObjs,numClusters,numCoords) \
                    shared(objects,clusters,membership,newClusters,newClusterSize) \
                    schedule(static) \
                    reduction(+:delta)
            for (i=0; i<numObjs; i++) {
                /* find the array index of nestest cluster center */
                index = find_nearest_cluster(numClusters, numCoords, objects[i],
                                             clusters);

                /* if membership changes, increase delta by 1 */
                if (membership[i] != index) delta += 1.0;

                /* assign the membership to object i */
                membership[i] = index;

                /* update new cluster centers : sum of objects located within */
                #pragma omp atomic
                newClusterSize[index]++;
                for (j=0; j<numCoords; j++)
                    #pragma omp atomic
                    newClusters[index][j] += objects[i][j];
            }
        }
        else {
            #pragma omp parallel \
                    shared(objects,clusters,membership,local_newClusters,local_newClusterSize)  
            {
                int tid = omp_get_thread_num();
                #pragma omp for \
                            private(i,j,index) \
                            firstprivate(numObjs,numClusters,numCoords) \
                            schedule(static) \
                            reduction(+:delta)
                for (i=0; i<numObjs; i++) {
                    /* find the array index of nestest cluster center */
                    index = find_nearest_cluster(numClusters, numCoords,
                                                 objects[i], clusters);

                    /* if membership changes, increase delta by 1 */
                    if (membership[i] != index) delta += 1.0;

                    /* assign the membership to object i */
                    membership[i] = index;

                    /* update new cluster centers : sum of all objects located
                       within (average will be performed later) */
                    local_newClusterSize[tid][index]++;
                    for (j=0; j<numCoords; j++)
                        local_newClusters[tid][index][j] += objects[i][j];
                }
            } /* end of #pragma omp parallel */

            /* let the main thread perform the array reduction */
            //for (i=0; i<numClusters; i++) {
                //for (j=0; j<nthreads; j++) {
              for (i = 0; i < nthreads; i++) {
                for (j = 0; j < numClusters; j++) {
                    newClusterSize[j] += local_newClusterSize[i][j];
                    local_newClusterSize[i][j] = 0.0;
                    for (k=0; k<numCoords; k++) {
                        newClusters[j][k] += local_newClusters[i][j][k];
                        local_newClusters[i][j][k] = 0.0;
                    }
                }
            }
        }
//timing2=omp_get_wtime()-timing2;
//printf("the time for expectation is %f\n", timing2);
	timing2=omp_get_wtime();
        /* average the sum and replace old cluster centers with newClusters */
        //#pragma omp parallel for private(i, j)
        for (i=0; i<numClusters; i++) {
            for (j=0; j<numCoords; j++) {
                if (newClusterSize[i] > 1) {
                    newClusterSize_inverse = 1.0 / newClusterSize[i];
                    clusters[i][j] = newClusters[i][j] * newClusterSize_inverse;
                }
                newClusters[i][j] = 0.0;   /* set back to 0 */
            }
            newClusterSize[i] = 0;   /* set back to 0 */
        }
//            timing2=omp_get_wtime()-timing2;
//printf("the time for maximization is %f\n", timing2);
        //delta /= numObjs;
    } while (delta > threshold_mul && loop++ < 500);

    if (_debug) {
        timing = omp_get_wtime() - timing;
        printf("nloops = %2d (T = %7.4f)",loop,timing);
    }

    if (!is_perform_atomic) {
        free(local_newClusterSize[0]);
        free(local_newClusterSize);

        for (i=0; i<nthreads; i++)
            for (j=0; j<numClusters; j++)
                free(local_newClusters[i][j]);
        free(local_newClusters[0]);
        free(local_newClusters);
    }
    free(newClusters[0]);
    free(newClusters);
    free(newClusterSize);
    //free(nearest_dist);
    //free(reallo_objs[0]);
    //free(reallo_objs);
    return clusters;
}

