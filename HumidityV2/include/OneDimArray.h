/*
 * OneDimArray.h
 *
 *  Created on: Apr 1, 2018
 *      Author: chuckjia
 */

#ifndef INCLUDE_ONEDIMARRAY_H_
#define INCLUDE_ONEDIMARRAY_H_

template<class T> class OneDimArray {

public:
	OneDimArray() {
		len_ = 0;
		data_area_ = 0;
	}

	OneDimArray(int _len) {
		len_ = _len;
		data_area_ = 0;
		initializeStorage();
	}

	OneDimArray(int _len, T* arr) {
		len_ = _len;

		data_area_ = 0;
		initializeStorage();

		memcpy(data_area_, arr, len_ * sizeof(T));
	}

	OneDimArray(const OneDimArray<T> &other) {
		assert(this != &other);

		data_area_ = 0;

		*this = other;
	}

	~OneDimArray() {
		deallocateStorage();
	}

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Assignment Operator
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	OneDimArray<T> &operator=(const OneDimArray<T> &other) {
		if (this == &other)
			return *this;

		if (!data_area_ || len_ != other.len()) {
			len_ = other.len();
			initializeStorage();
		}

		memcpy(data_area_, other.dataPtr(), len_ * sizeof(T));
		return *this;
	}

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Element Access
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	inline T &operator[](int i) const { return data_area_[i]; }

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Informational
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	inline int len() const { return len_; }
	inline T* dataPtr() const { return data_area_; }

protected:
	void deallocateStorage() {
		if (data_area_)
			delete[] data_area_;
		data_area_ = 0;
	}

	void initializeStorage() {
		if (data_area_)
			deallocateStorage();

		if (len_ > 0)
			data_area_ = new T[len_];
		else
			data_area_ = 0;
	}

	int len_;
	T* data_area_;
};



#endif /* INCLUDE_ONEDIMARRAY_H_ */
