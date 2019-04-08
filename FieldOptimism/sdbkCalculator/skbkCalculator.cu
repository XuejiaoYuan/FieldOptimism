#include "sdbkCalculator.cuh"

void SdBkCalc::discreteRayCasting() {
	int nThreads = 512;
	dim3 nBlocks;
	//GeometryFunc::setThreadsBlocks(nBlocks, nThreads, )
}

__device__ __host__ void SdBkCalc::calcIntersection3DDDA(Heliostat* helio, const Vector3d&dir, int*& esimate_grid) {

}

__device__ __host__ void SdBkCalc::_helio_calc(int index) {
	
}

