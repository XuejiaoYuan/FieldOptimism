#pragma once
#include "../Common/CommonFunc.h"

template<typename T>
class DeviceList {
public:
	T* d_list;
	T* backup_list;
	int listLength;
	int backupSize;
	int listSize;
	DeviceList(int _list_size = DEVICE_LIST_SIZE) :d_list(nullptr), backup_list(nullptr), listLength(0), listSize(_list_size) {
		if (_list_size <= 0)
			listSize = DEVICE_LIST_SIZE;
		//cudaMallocManaged(&d_list, listSize * sizeof(T));
		cudaMalloc((void**)&d_list, listSize * sizeof(T));
	};
	void mallocBackupSpace(int _backup_size = DEVICE_BACKUP_SIZE) {
		if (!backup_list) {
			//cudaMallocManaged(&backup_list, _backup_size * sizeof(T));
			cudaMalloc((void**)&backup_list, _backup_size * sizeof(T));
			backupSize = _backup_size;
		}
	};
	__host__ __device__ void append(T& data) {
		if (listLength < listSize)
			d_list[listLength++] = data;
	}
	~DeviceList() {
		cudaFree(d_list);
		if (backup_list) cudaFree(backup_list);

		d_list = nullptr;
		backup_list = nullptr;
	}
	void* operator new(size_t size) {
		void * ptr;
		//cudaMallocManaged(&ptr, size);
		cudaMalloc((void**)&ptr, size);
		return ptr;
	}
	void operator delete(void* ptr) {
		cudaFree(ptr);
		ptr = nullptr;
	}
	void* operator new[](size_t size) {
		void* ptr;
		//cudaMallocManaged(&ptr, size);
		cudaMalloc((void**)&ptr, size);
		return ptr;
	}
		void operator delete[](void* ptr) {
		cudaFree(ptr);
		ptr = nullptr;
	}
};
