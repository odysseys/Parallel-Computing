#include<cutil_inline.h>  
#include<iostream>  

using namespace std; 

#define N 32

// Kernel definition
__global__ void MatAdd(float A[N], float B[N], float* C)  
{  
    int i = blockIdx.x * blockDim.x + threadIdx.x; //get thread index by built-in variables  
    if (i < N)  
        C[i] = A[i] + B[i];  
}       

int main()  
{  
    float A[N],B[N]; // host variable  
    float *dA, *dB; // device variable, to have same value with A,B  
    float *device_res, *host_res; // device and host result, to be device and host variable respectively  
  
    // initialize host variable  
    memset(A,0,sizeof(A));  
    memset(B,0,sizeof(B));  
    A[0] = 1;  
    B[0] = 2;  
  
  
    // allocate for device variable and set value to them  
    cudaMalloc((void**) &dA,N*sizeof(float));  
    cudaMalloc((void**) &dB,N*sizeof(float));  
    cudaMemcpy(dA, A, N*sizeof(float),cudaMemcpyHostToDevice);  
    cudaMemcpy(dB, B, N*sizeof(float),cudaMemcpyHostToDevice);  
  
    //malloc for host and device variable  
    host_res = (float*) malloc(N*sizeof(float));  
    cudaMalloc((void**)&device_res, N*sizeof(float));  
  
    // Kernel invocation  
    int threadsPerBlock = 16;  
    int numBlocks = N/threadsPerBlock;   
    MatAdd<<<numBlocks, threadsPerBlock>>>(dA, dB, device_res);  
  
    cudaMemcpy(host_res, device_res, N*sizeof(float),cudaMemcpyDeviceToHost); //copy from device to host  
      
    // validate  
    int i;  
    float sum = 0;  
    for(i=0;i<N;i++)  
        sum += host_res[i];  
    cout<<sum<<endl;  
  
    //free variables  
    cudaFree(dA);  
    cudaFree(dB);  
  
    cudaFree(device_res);  
    free(host_res);  
} 