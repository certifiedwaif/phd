// first_gpu_program.cu

#include <stdio.h>
#include <cuda_runtime.h>

#define N 65000

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

__global__ void add(int *a, int *b, int *c)
{
	int tid = blockIdx.x;
	if (tid < N)
		c[tid] = a[tid] + b[tid];
}

int main(int argc, char **argv)
{
	int a[N], b[N], c[N];
	int *dev_a, *dev_b, *dev_c;

	// Allocate memory on the GPU
	HANDLE_ERROR(cudaMalloc((void**) &dev_a, N * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**) &dev_b, N * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**) &dev_c, N * sizeof(int)));

	// Fill the arrays 'a' and ''b' on the CPU
	for (int i = 0; i < N; i++) {
		a[i] = -i;
		b[i] = i * i;
	}

	// Copy the arrays 'a' and 'b' to the GPU
	HANDLE_ERROR(cudaMemcpy(dev_a, a, N * sizeof(int),
					cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_b, b, N * sizeof(int),
					cudaMemcpyHostToDevice));

	add<<<N, 1>>>(dev_a, dev_b, dev_c);

	// Copy array 'c' back from the GPU to the CPU
	HANDLE_ERROR(cudaMemcpy(c, dev_c, N * sizeof(int),
					cudaMemcpyDeviceToHost));

	// Display the results
	for (int i = 0; i < N; i++) {
		if (i % 256 == 0)
			printf("%d + %d = %d\n", a[i], b[i], c[i]);
	}

	// Free the memory allocated on the GPU
	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_c);

	return 0;
}

