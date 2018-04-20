/*
 * TwoDimArray.h
 *
 *  Created on: Mar 5, 2018
 *      Author: chuckjia
 */

#ifndef INCLUDE_TWODIMARRAY_H_
#define INCLUDE_TWODIMARRAY_H_

#include <assert.h>
#include <string>

#include "debug.h"

template<class T> class TwoDimArray {
public:

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Constructors
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	TwoDimArray() {
		data_slices_ = 0;
		data_area_ = 0;
		nrow_ = ncol_ = 0;
	}

	TwoDimArray(int _nrow, int _ncol) {
		nrow_ = _nrow;
		ncol_ = _ncol;

		data_slices_ = 0;
		data_area_ = 0;

		initializeStorage();
	}

	TwoDimArray(int _nrow, int _ncol, const T *array) {
		nrow_ = _nrow;
		ncol_ = _ncol;

		data_slices_ = 0;
		data_area_ = 0;

		initializeStorage();

		memcpy(data_area_, array, nrow_ * ncol_ * sizeof(T));
	}

	TwoDimArray(const TwoDimArray<T> &other) {
		assert(this != &other);

		data_slices_ = 0;
		data_area_ = 0;

		*this = other;
	}


	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Destructor
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	~TwoDimArray() {
		deallocateStorage();
	}


	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Assignment Operator
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	TwoDimArray<T> &operator=(const TwoDimArray<T> &other) {
		if(this == &other)
			return *this;

		if(!data_slices_ || nrow_ != other.nrow() || ncol_ != other.ncol()) {
			nrow_ = other.nrow();
			ncol_ = other.ncol();

			initializeStorage();
		}

		memcpy(data_area_, other.data_area_, nrow_ * ncol_ * sizeof(T));

		return *this;
	}


	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Element Access
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	inline T *operator[](int row) const { return data_slices_[row]; }


	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Informational
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	inline int nrow() const { return nrow_; }
	inline int ncol() const { return ncol_; }
	inline T *dataPtr() const { return data_area_; }
	inline T **rowPtrs() const { return data_slices_; }


protected:
	void deallocateStorage() {
		if(data_slices_) {
			delete[] data_slices_;
			delete[] data_area_;

			data_slices_ = 0;
			data_area_ = 0;
		}
	}

	void initializeStorage() {
		if(data_slices_)
			deallocateStorage();

		if(nrow_ > 0) {
			data_slices_ = new T *[nrow_];
			data_area_ = new T[nrow_ * ncol_];
		} else {
			data_slices_ = 0;
			data_area_ = 0;
		}

		T *curr_ptr = data_area_;
		for(int i = 0; i < nrow_; i++, curr_ptr += ncol_)
			data_slices_[i] = curr_ptr;
	}

	// Size of the matrix
	int nrow_, ncol_;

	// Array of pointers to the beginning of each matrix row in data_area
	T **data_slices_;

	// Pointer to actual matrix data
	T *data_area_;
};


#endif /* INCLUDE_TWODIMARRAY_H_ */

