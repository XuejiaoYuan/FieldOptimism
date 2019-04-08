


#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <chrono>
#include <stdio.h>
#include <string>
#include <iostream>
#include "DataStructure\Timer.h"
using namespace std;




class A {
	int m;
public:
	A() :m(1) {}
	__device__ __host__ int getM() { return m; }
	void* operator new(size_t len) {
		void* ptr;
		cudaMallocManaged(&ptr, len);
		return ptr;
	}
};

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);
void addWithCuda2(int* c, const int*a, const int *b, unsigned int size);

__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

int main()
{
	Timer t;
	t.resetStart();
    const int arraySize = 5;
    const int A[arraySize] = { 1, 2, 3, 4, 5 };
    const int B[arraySize] = { 10, 20, 30, 40, 50 };
	int C[arraySize] = { 0 };

	addWithCuda(C, A, B, arraySize);

	t.printDuration("cudaMalloc");
	//if (cudaStatus != cudaSuccess) {
	//    fprintf(stderr, "addWithCuda failed!");
	//    return 1;
	//}
	
	t.resetStart();
	int *a, *b;
	cudaMallocManaged(&a, arraySize * sizeof(int));
	cudaMallocManaged(&b, arraySize * sizeof(int));
	int *c = NULL;
	cudaMallocManaged(&c, arraySize * sizeof(int));

	for (int i = 0; i < arraySize; ++i) {
		a[i] = i;
		b[i] = 10 * i;
	}
	for(int i=0; i<10000; ++i)
		addWithCuda2(c, a, b, arraySize);
	t.printDuration("cudaMallocManaged");

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    //cudaStatus = cudaDeviceReset();
    //if (cudaStatus != cudaSuccess) {
    //    fprintf(stderr, "cudaDeviceReset failed!");
    //    return 1;
    //}

    return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
	Timer t;
	t.resetStart();
	cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }
	
    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }
	t.printDuration("malloc memcpy");

	t.resetStart();
    // Launch a kernel on the GPU with one thread for each element.
	for (int i = 0; i < 10000; ++i) {
		addKernel << <1, size >> >(dev_c, dev_a, dev_b);
		cudaStatus = cudaDeviceSynchronize();
	}
	t.printDuration("Ñ­»·");

    // Check for any errors launching the kernel
    //cudaStatus = cudaGetLastError();
    //if (cudaStatus != cudaSuccess) {
    //    fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
    //    goto Error;
    //}
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}


void addWithCuda2(int *c, const int*a, const int *b, unsigned int size) {
	//cudaMallocManaged(&b, size * sizeof(int));
	//cudaMallocManaged(&c, size * sizeof(int));

	addKernel << <1, size >> > (c, b, a);
	cudaDeviceSynchronize();
}