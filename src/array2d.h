///////////////////////////////////////////////////////////////////////////////
///
///	\file	 array2d.h
///	\author  Paul Ullrich
///	\version March 7, 2024
///
///	<remarks>
///		Copyright 2024 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the BSD-3
///		License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _tmp_array2d_h_
#define _tmp_array2d_h_

namespace tmp {

///////////////////////////////////////////////////////////////////////////////

template <typename T>
class array2d {
	private:
		///	<summary>
		///		Size in the y direction.
		///	</summary>
		size_t m_size_y;

		///	<summary>
		///		Size in the x direction.
		///	</summary>
		size_t m_size_x;

		///	<summary>
		///		Pointer to contiguous memory block.
		///	</summary>
		T * m_vec;

	public:
		///	<summary>
		///		Default constructor.
		///	</summary>
		array2d() :
			m_size_x(0),
			m_size_y(0),
			m_vec(nullptr)
		{ }

		///	<summary>
		///		Constructor.
		///	</summary>
		array2d(
			size_t a_size_y,
			size_t a_size_x
		) :
			m_size_x(a_size_x),
			m_size_y(a_size_y)
		{
			m_vec = new T[m_size_y * m_size_x];
		}

		///	<summary>
		///		Constructor with initializer.
		///	</summary>
		array2d(
			size_t a_size_y,
			size_t a_size_x,
			const T & t
		) :
			m_size_x(a_size_x),
			m_size_y(a_size_y)
		{
			m_vec = new T[m_size_y * m_size_x];

			std::fill(m_vec, m_vec + (m_size_y * m_size_x), t);
		}

		///	<summary>
		///		Destructor.
		///	</summary>
		~array2d() {
			delete[] m_vec;
		}

	public:
		///	<summary>
		///		Resize.
		///	</summary>
		void resize(
			size_t a_size_y,
			size_t a_size_x
		) {
			if (m_vec != nullptr) {
				delete[] m_vec;
			}

			m_size_y = a_size_y;
			m_size_x = a_size_x;

			m_vec = new T[m_size_y * m_size_x];
		}

		///	<summary>
		///		Resize with initializer.
		///	</summary>
		void resize(
			size_t a_size_y,
			size_t a_size_x,
			const T & t
		) {
			if (m_vec != nullptr) {
				delete[] m_vec;
			}

			m_size_y = a_size_y;
			m_size_x = a_size_x;

			m_vec = new T[m_size_y * m_size_x];

			std::fill(m_vec, m_vec + (m_size_y * m_size_x), t);
		}

		///	<summary>
		///		Clear.
		///	</summary>
		void clear() {
			if (m_vec != nullptr) {
				delete[] m_vec;
			}

			m_size_y = 0;
			m_size_x = 0;
			m_vec = nullptr;
		}

	public:
		///	<summary>
		///		Get the y size.
		///	</summary>
		size_t size_y() const {
			return m_size_y;
		}

		///	<summary>
		///		Get the x size.
		///	</summary>
		size_t size_x() const {
			return m_size_x;
		}

		///	<summary>
		///		Accessor.
		///	</summary>
		T & operator()(size_t j, size_t i) {
			return m_vec[j * m_size_x + i];
		}

		///	<summary>
		///		Const accessor.
		///	</summary>
		const T & operator()(size_t j, size_t i) const {
			return m_vec[j * m_size_x + i];
		}

	public:
		///	<summary>
		///		Cast to pointer.
		///	</summary>
		operator T*() {
			return m_vec;
		}

		///	<summary>
		///		Cast to const pointer.
		///	</summary>
		operator T const *() const {
			return m_vec;
		}

		///	<summary>
		///		Access the 1D data structure.
		///	</summary>
		T* data1d() {
			return m_vec;
		}

		///	<summary>
		///		Access the 1D data structure.
		///	</summary>
		T const * data1d() const {
			return m_vec;
		}
};

///////////////////////////////////////////////////////////////////////////////

} // namespace tmp

#endif //_tmp_array2d_h_

