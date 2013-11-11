/*
  Copyright (c) 2012,2013 Tobias Brink

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef TERSOFF_NDARRAY_HPP
#define TERSOFF_NDARRAY_HPP

#include <cstdlib>
#include <cmath>
#include <stdexcept>

namespace mytest {
  // Some common integer powers. /////////////////////////////////////////
  template<typename T>
  T pow2(const T& b) {
    return b * b;
  }
  template<typename T>
  T pow3(const T& b) {
    return b * b * b;
  }

  // 3D vectors. /////////////////////////////////////////////////////////

  template<typename T> class Vec3D;

  //! Addition of Vec3D with scalar.
  template<typename T>
  Vec3D<T> operator+(Vec3D<T> lhs, const T& rhs) {
    lhs += rhs;
    return lhs;
  }

  //! Addition of Vec3D with Vec3D.
  template<typename T>
  Vec3D<T> operator+(Vec3D<T> lhs, const Vec3D<T>& rhs) {
    lhs += rhs;
    return lhs;
  }

  //! Add scalar to Vec3D.
  template<typename T>
  Vec3D<T> operator+(const T& lhs, Vec3D<T> rhs) {
    rhs += lhs;
    return rhs;
  }

  //! Subtract scalar from Vec3D.
  template<typename T>
  Vec3D<T> operator-(Vec3D<T> lhs, const T& rhs) {
    lhs -= rhs;
    return lhs;
  }

  //! Subtract Vec3D from Vec3D.
  template<typename T>
  Vec3D<T> operator-(Vec3D<T> lhs, const Vec3D<T>& rhs) {
    lhs -= rhs;
    return lhs;
  }

  //! Subtract Vec3D from scalar.
  template<typename T>
  Vec3D<T> operator-(const T& lhs, const Vec3D<T>& rhs) {
    return -rhs + lhs;
  }

  //! Multiplication of Vec3D with scalar.
  template<typename T>
  Vec3D<T> operator*(Vec3D<T> lhs, const T& rhs) {
    lhs *= rhs;
    return lhs;
  }

  //! Multiplication of Vec3D with Vec3D.
  template<typename T>
  Vec3D<T> operator*(Vec3D<T> lhs, const Vec3D<T>& rhs) {
    lhs *= rhs;
    return lhs;
  }

  //! Multiply scalar with Vec3D.
  template<typename T>
  Vec3D<T> operator*(const T& lhs, Vec3D<T> rhs) {
    rhs *= lhs;
    return rhs;
  }

  //! Devide Vec3D by scalar.
  template<typename T>
  Vec3D<T> operator/(Vec3D<T> lhs, const T& rhs) {
    lhs /= rhs;
    return lhs;
  }

  //! Devide Vec3D by Vec3D.
  template<typename T>
  Vec3D<T> operator/(Vec3D<T> lhs, const Vec3D<T>& rhs) {
    lhs /= rhs;
    return lhs;
  }

  //! Devide scalar by Vec3D.
  template<typename T>
  Vec3D<T> operator/(const T& lhs, Vec3D<T> rhs) {
    rhs.x = lhs / rhs.x;
    rhs.y = lhs / rhs.y;
    rhs.z = lhs / rhs.z;
    return rhs;
  }

  //! Vec3D dot product.
  template<typename T>
  T dot(const Vec3D<T>& lhs, const Vec3D<T>& rhs) {
    return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
  }

  /*!
    A 3D vector (x,y,z).
  */
  template<typename T>
  class Vec3D {
    friend T dot<T>(const Vec3D<T>&, const Vec3D<T>&);
    friend Vec3D<T> operator/<T>(const T&, Vec3D<T>);
  public:
    Vec3D();
    Vec3D(const T& x, const T& y, const T& z);
    Vec3D(const Vec3D<T>& rhs);

    T& operator[](int);
    const T& operator[](int i) const;

    Vec3D<T>& operator+=(const Vec3D<T>& rhs);
    Vec3D<T>& operator+=(const T& rhs);

    Vec3D<T>& operator-=(const Vec3D<T>& rhs);
    Vec3D<T>& operator-=(const T& rhs);
    Vec3D<T> operator-() const;

    Vec3D<T>& operator*=(const Vec3D<T>& rhs);
    Vec3D<T>& operator*=(const T& rhs);

    Vec3D<T>& operator/=(const Vec3D<T>& rhs);
    Vec3D<T>& operator/=(const T& rhs);

    T abs() const;
    T abssq() const;
  private:
    T x,y,z;
  };

  template<typename T> inline
  Vec3D<T>::Vec3D() { }

  template<typename T> inline
  Vec3D<T>::Vec3D(const T& x, const T& y, const T& z)
    : x(x), y(y), z(z)
  { }

  template<typename T> inline
  Vec3D<T>::Vec3D(const Vec3D<T>& rhs)
    : x(rhs.x), y(rhs.y), z(rhs.z)
  { }

  template<typename T> inline
  T& Vec3D<T>::operator[](int i) {
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      throw std::out_of_range("i must be in the range [0;2]");
    }
  }

  template<typename T> inline
  const T& Vec3D<T>::operator[](int i) const {
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
        throw std::out_of_range("i must be in the range [0;2]");
    }
  }

  template<typename T> inline
  Vec3D<T>& Vec3D<T>::operator+=(const Vec3D<T>& rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  }

  template<typename T> inline
  Vec3D<T>& Vec3D<T>::operator+=(const T& rhs) {
    x += rhs;
    y += rhs;
    z += rhs;
    return *this;
  }

  template<typename T> inline
  Vec3D<T>& Vec3D<T>::operator-=(const Vec3D<T>& rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  }

  template<typename T> inline
  Vec3D<T>& Vec3D<T>::operator-=(const T& rhs) {
    x -= rhs;
    y -= rhs;
    z -= rhs;
    return *this;
  }

  template<typename T> inline
  Vec3D<T> Vec3D<T>::operator-() const {
    return Vec3D<T>(-this->x, -this->y, -this->z);
  }

  template<typename T> inline
  Vec3D<T>& Vec3D<T>::operator*=(const Vec3D<T>& rhs) {
    x *= rhs.x;
    y *= rhs.y;
    z *= rhs.z;
    return *this;
  }

  template<typename T> inline
  Vec3D<T>& Vec3D<T>::operator*=(const T& rhs) {
    x *= rhs;
    y *= rhs;
    z *= rhs;
    return *this;
  }

  template<typename T> inline
  Vec3D<T>& Vec3D<T>::operator/=(const Vec3D<T>& rhs) {
    x /= rhs.x;
    y /= rhs.y;
    z /= rhs.z;
    return *this;
  }

  template<typename T> inline
  Vec3D<T>& Vec3D<T>::operator/=(const T& rhs) {
    x /= rhs;
    y /= rhs;
    z /= rhs;
    return *this;
  }

  template<typename T> inline
  T Vec3D<T>::abs() const {
    return std::sqrt(x*x + y*y + z*z);
  }

  template<typename T> inline
  T Vec3D<T>::abssq() const {
    return x*x + y*y + z*z;
  }


  // n-dimensional arrays. ///////////////////////////////////////////////

  /*!
    This wraps a pointer to type T to be used as a 2D array.

    @todo: move inline member functions out of the class definition like above.
  */
  template<typename T>
  class Array2D {
  public:
    //! External memory management.
    explicit Array2D(T* data, int extent_outer, int extent_inner)
      : data(data), extent_outer(extent_outer), extent_inner(extent_inner),
        keep_data(true)
    {}

    //! Manage data by the class.
    explicit Array2D(int extent_outer, int extent_inner)
      : data(static_cast<T*>(std::malloc(extent_outer*extent_inner*sizeof(T)))),
        extent_outer(extent_outer), extent_inner(extent_inner),
        keep_data(false)
    {}

    ~Array2D() {
      if (!keep_data)
        std::free(data);
    }

    //! Return size of array in dimension e (0 is the first dimension).
    int extent(int e) const {
      switch(e) {
      case 0:
        return extent_outer;
      case 1:
        return extent_inner;
      default:
        throw std::out_of_range("e must be in the range [0;1]");
      }
    }

    T& operator()(int i, int j) {
      //TODO: range check??
      return data[i*extent_inner + j];
    }
    const T& operator()(int i, int j) const {
      //TODO: range check??
      return data[i*extent_inner + j];
    }

  private:
    T* data;
    int extent_outer, extent_inner;
    bool keep_data;

    // Do not use!
    Array2D& operator=(const Array2D&) {}
  };


  /*!
    This wraps a pointer to type T to be used as a 2D array.

    @todo: move inline member functions out of the class definition like above.
  */
  template<typename T>
  class Array3D {
  public:
    explicit Array3D(T* data, int extent0, int extent1, int extent2)
      : data(data), extent0(extent0), extent1(extent1), extent2(extent2),
        fac_i(extent1*extent2),
        keep_data(true)
    {}

    explicit Array3D(int extent0, int extent1, int extent2)
      : data(static_cast<T*>(malloc(extent0*extent1*extent2*sizeof(T)))),
        extent0(extent0), extent1(extent1), extent2(extent2),
        fac_i(extent1*extent2),
        keep_data(false)
    {}

    ~Array3D() {
      if (!keep_data)
        std::free(data);
    }

    //! Return size of array in dimension e (0 is the first dimension).
    int extent(int e) const {
      switch(e) {
      case 0:
        return extent0;
      case 1:
        return extent1;
      case 2:
        return extent2;
      default:
        throw std::out_of_range("e must be in the range [0;2]");
      }
    }

    T& operator()(int i, int j, int k) {
      //TODO: range check??
      return data[i*fac_i + j*extent2 + k];
    }
    const T& operator()(int i, int j, int k) const {
      //TODO: range check??
      return data[i*fac_i + j*extent2 + k];
    }

    void operator=(const T& rhs) {
      for (int i = 0; i != extent0; ++i)
        for (int j = 0; j != extent1; ++j)
          for (int k = 0; k != extent2; ++k)
            (*this)(i,j,k) = rhs;
    }

    bool any() const {
      for (int i = 0; i != extent0; ++i)
        for (int j = 0; j != extent1; ++j)
          for (int k = 0; k != extent2; ++k)
            if ((*this)(i,j,k))
              return true;
      return false;
    }

    bool all() const {
      for (int i = 0; i != extent0; ++i)
        for (int j = 0; j != extent1; ++j)
          for (int k = 0; k != extent2; ++k)
            if (!(*this)(i,j,k))
              return false;
      return true;
    }

  private:
    T* data;
    int extent0, extent1, extent2;
    int fac_i;
    bool keep_data;

    // Do not use!
    Array3D& operator=(const Array3D&) {}
  };


  /*!
    This wraps a pointer to type T to be used as a 2D array.

    @todo: move inline member functions out of the class definition like above.
  */
  template<typename T>
  class Array4D {
  public:
    explicit Array4D(T* data,
                     int extent0, int extent1, int extent2, int extent3)
      : data(data),
        extent0(extent0), extent1(extent1), extent2(extent2), extent3(extent3),
        fac_i(extent1*extent2*extent3),
        fac_j(extent2*extent3),
        keep_data(true)
    {}

    explicit Array4D(int extent0, int extent1, int extent2, int extent3)
      : data(static_cast<T*>(malloc(extent0*extent1*extent2*extent3*sizeof(T)))),
        extent0(extent0), extent1(extent1), extent2(extent2), extent3(extent3),
        fac_i(extent1*extent2*extent3),
        fac_j(extent2*extent3),
        keep_data(false)
    {}

    ~Array4D() {
      if (!keep_data)
        std::free(data);
    }

    //! Return size of array in dimension e (0 is the first dimension).
    int extent(int e) const {
      switch(e) {
      case 0:
        return extent0;
      case 1:
        return extent1;
      case 2:
        return extent2;
      case 3:
        return extent3;
      default:
        throw std::out_of_range("e must be in the range [0;3]");
      }
    }

    T& operator()(int i, int j, int k, int l) {
      //TODO: range check??
      return data[i*fac_i + j*fac_j + k*extent3 + l];
    }
    const T& operator()(int i, int j, int k, int l) const {
      //TODO: range check??
      return data[i*fac_i + j*fac_j + k*extent3 + l];
    }

    void operator=(const T& rhs) {
      for (int i = 0; i != extent0; ++i)
        for (int j = 0; j != extent1; ++j)
          for (int k = 0; k != extent2; ++k)
            for (int l = 0; l != extent3; ++l)
              (*this)(i,j,k,l) = rhs;
    }

    bool any() const {
      for (int i = 0; i != extent0; ++i)
        for (int j = 0; j != extent1; ++j)
          for (int k = 0; k != extent2; ++k)
            for (int l = 0; l != extent3; ++l)
              if ((*this)(i,j,k,l))
                return true;
      return false;
    }

    bool all() const {
      for (int i = 0; i != extent0; ++i)
        for (int j = 0; j != extent1; ++j)
          for (int k = 0; k != extent2; ++k)
            for (int l = 0; l != extent3; ++l)
              if (!(*this)(i,j,k,l))
                return false;
      return true;
    }

  private:
    T* data;
    int extent0, extent1, extent2, extent3;
    int fac_i, fac_j;
    bool keep_data;

    // Do not use!
    Array4D& operator=(const Array4D&) {}
  };


}

#endif /* TERSOFF_NDARRAY_HPP */
