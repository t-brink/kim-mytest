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

  //! Vec3D cross product.
  template<typename T>
  Vec3D<T> cross(const Vec3D<T>& lhs, const Vec3D<T>& rhs) {
    return Vec3D<T>(lhs.y * rhs.z - lhs.z * rhs.y,
                    lhs.z * rhs.x - lhs.x * rhs.z,
                    lhs.x * rhs.y - lhs.y * rhs.x);
  }

  /*!
    A 3D vector (x,y,z).
  */
  template<typename T>
  class Vec3D {
    friend T dot<T>(const Vec3D<T>&, const Vec3D<T>&);
    friend Vec3D<T> cross<T>(const Vec3D<T>&, const Vec3D<T>&);
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


  // Symmetric 3*3 matrix in Voigt notation //////////////////////////////

  template<typename T> class Voigt6;

  // Multiply Voigt6 with scalar.
  template<typename T>
  Voigt6<T> operator*(Voigt6<T> lhs, const T& rhs) {
    lhs *= rhs;
    return lhs;
  }

  // Multiply scalar with Voigt6.
  template<typename T>
  Voigt6<T> operator*(const T& lhs, Voigt6<T> rhs) {
    rhs *= lhs;
    return rhs;
  }

  // Divide Voigt6 by scalar.
  template<typename T>
  Voigt6<T> operator/(Voigt6<T> lhs, const T& rhs) {
    lhs /= rhs;
    return lhs;
  }

  // Uses {xx,yy,zz,yz,xz,xy} to store data.  Data is stored
  // continuously.
  template<typename T>
  class Voigt6 {
  public:
    Voigt6(T xx_, T yy_, T zz_, T yz_, T xz_, T xy_)
      : data(static_cast<T*>(std::malloc(6*sizeof(T)))), keep_data(false),
        xx(data[0]), yy(data[1]), zz(data[2]),
        yz(data[3]), xz(data[4]), xy(data[5])
    {
      data[0] = xx_; data[1] = yy_; data[2] = zz_;
      data[3] = yz_; data[4] = xz_; data[5] = xy_;
    }

    // Uninitialized.
    Voigt6()
      : data(static_cast<T*>(std::malloc(6*sizeof(T)))), keep_data(false),
        xx(data[0]), yy(data[1]), zz(data[2]),
        yz(data[3]), xz(data[4]), xy(data[5])
    {}

    // Wrapper to externally managed memory.
    explicit Voigt6(T* d)
      : data(d), keep_data(true),
        xx(data[0]), yy(data[1]), zz(data[2]),
        yz(data[3]), xz(data[4]), xy(data[5])
    {}

    // Copy constructor.
    Voigt6(const Voigt6<T>& other)
      : data(static_cast<T*>(std::malloc(6*sizeof(T)))), keep_data(false),
        xx(data[0]), yy(data[1]), zz(data[2]),
        yz(data[3]), xz(data[4]), xy(data[5])
    {
      data[0] = other.data[0]; data[1] = other.data[1]; data[2] = other.data[2];
      data[3] = other.data[3]; data[4] = other.data[4]; data[5] = other.data[5];
    }

    // Move constructor.
    Voigt6(Voigt6<T>&& other)
      : xx(other.data[0]), yy(other.data[1]), zz(other.data[2]),
        yz(other.data[3]), xz(other.data[4]), xy(other.data[5])
    {
      data = other.data;
      keep_data = other.keep_data;
      other.keep_data = true; // Keep from double-freeing memory.
    }

    ~Voigt6() {
      if (!keep_data)
        std::free(data);
    }

    T& operator()(int i) {
      //TODO: range check??
      return data[i];
    }
    const T& operator()(int i) const {
      //TODO: range check??
      return data[i];
    }

    T& operator()(int i, int j) {
      switch (i) {
      case 0:
        switch (j) {
        case 0: return data[0]; // xx
        case 1: return data[5]; // xy
        case 2: return data[4]; // xz
        default: throw std::runtime_error("out of bounds");
        }
      case 1:
        switch (j) {
        case 0: return data[5]; // xy
        case 1: return data[1]; // yy
        case 2: return data[3]; // yz
        default: throw std::runtime_error("out of bounds");
        }
      case 2:
        switch (j) {
        case 0: return data[4]; // xz
        case 1: return data[3]; // yz
        case 2: return data[2]; // zz
        default: throw std::runtime_error("out of bounds");
        }
      default:
        throw std::runtime_error("out of bounds");
      }
    }
    const T& operator()(int i, int j) const {
      switch (i) {
      case 0:
        switch (j) {
        case 0: return data[0]; // xx
        case 1: return data[5]; // xy
        case 2: return data[4]; // xz
        default: throw std::runtime_error("out of bounds");
        }
      case 1:
        switch (j) {
        case 0: return data[5]; // xy
        case 1: return data[1]; // yy
        case 2: return data[3]; // yz
        default: throw std::runtime_error("out of bounds");
        }
      case 2:
        switch (j) {
        case 0: return data[4]; // xz
        case 1: return data[3]; // yz
        case 2: return data[2]; // zz
        default: throw std::runtime_error("out of bounds");
        }
      default:
        throw std::runtime_error("out of bounds");
      }
    }

    Voigt6<T>& operator*=(const T& rhs) {
      data[0] *= rhs; data[1] *= rhs; data[2] *= rhs;
      data[3] *= rhs; data[4] *= rhs; data[5] *= rhs;
      return *this;
    }

    Voigt6<T>& operator/=(const T& rhs) {
      data[0] /= rhs; data[1] /= rhs; data[2] /= rhs;
      data[3] /= rhs; data[4] /= rhs; data[5] /= rhs;
      return *this;
    }

  private:
    T* data;
    bool keep_data;
  public:
    const T& xx, yy, zz, yz, xz, xy;
  };


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

    Array2D<T>& operator*=(const T& rhs) {
      const int size = extent_outer * extent_inner;
      for (int i = 0; i != size; ++i)
        data[i] *= rhs;
      return *this;
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
