/*
  Copyright (c) 2017 Tobias Brink

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

#include "core.hpp"
#include "ndarray.hpp"

#include <vector>
#include <cmath>
#include <random>

using namespace std;
using namespace mytest;

double dist_sq(const Vec3D<double>& u,
               const Vec3D<double>& v,
               double width, double hx,
               double height, double hy,
               double depth, double hz,
               bool bcx, bool bcy, bool bcz) {
  double dx = u[0] - v[0];
  if (bcx) {
    if (dx > hx) {
      dx -= width;
    } else if (dx < -hx) {
      dx += width;
    }
  }
  double dy = u[1] - v[1];
  if (bcy) {
    if (dy > hy) {
      dy -= height;
    } else if (dy < -hy) {
      dy += height;
    }
  }
  double dz = u[2] - v[2];
  if (bcz) {
    if (dz > hz) {
      dz -= depth;
    } else if (dz < -hz) {
      dz += depth;
    }
  }
  return dx*dx + dy*dy + dz*dz;
}

bool in_neighborhood(const Vec3D<double>& p,
                     int i, int j, int k,
                     int len_i, int len_j, int len_k,
                     const Array3D< Vec3D<double> >& grid,
                     const Array3D<bool>& grid_filled,
                     double min_dist_sq,
                     double width, double hx,
                     double height, double hy,
                     double depth, double hz,
                     bool bcx, bool bcy, bool bcz) {
  const int imin = bcx ? i - 2 : max(0, i - 2);
  const int imax = bcx ? i + 3 : min(len_i, i + 3);
  const int jmin = bcy ? j - 2 : max(0, j - 2);
  const int jmax = bcy ? j + 3 : min(len_j, j + 3);
  const int kmin = bcz ? k - 2 : max(0, k - 2);
  const int kmax = bcz ? k + 3 : min(len_k, k + 3);
  for (int ii = imin; ii < imax; ++ii) {
    int iii = ii;
    if (bcx) { iii %= len_i; iii = iii < 0 ? iii + len_i : iii; }
    for (int jj = jmin; jj < jmax; ++jj) {
      int jjj = jj;
      if (bcy) { jjj %= len_j; jjj = jjj < 0 ? jjj + len_j : jjj; }
      for (int kk = kmin; kk < kmax; ++kk) {
        int kkk = kk;
        if (bcz) { kkk %= len_k; kkk = kkk < 0 ? kkk + len_k : kkk; }
        if (!grid_filled(iii,jjj,kkk))
          continue;
        const double d =
          dist_sq(grid(iii,jjj,kkk), p,
                  width, hx, height, hy, depth, hz,
                  bcx, bcy, bcz);
        if (d <= min_dist_sq)
          return true;
      }
    }
  }
  return false;
}

/* The output vector is a flattened natoms*3 array. */
void poisson_disk_sampling(double width, double height, double depth,
                           double min_dist,
                           vector<double>& positions, // <-- output
                           mt19937& rng) {
  const double cell_size = min_dist / sqrt(3);

  const unsigned len_i = ceil(width / cell_size);
  const unsigned len_j = ceil(height / cell_size);
  const unsigned len_k = ceil(depth / cell_size);

  Array3D< Vec3D<double> > grid(len_i, len_j, len_k);
  Array3D<bool> grid_filled(len_i, len_j, len_k);
  grid_filled = false;

  const double min_dist_sq = min_dist * min_dist;

  /*
  const double limit_xmin = 0;
  const double limit_xmax = width;
  const double limit_ymin = 0;
  const double limit_ymax = height;
  const double limit_zmin = 0;
  const double limit_zmax = depth;
  */
  const double hx = width / 2;
  const double hy = height / 2;
  const double hz = depth / 2;

  vector< Vec3D<double> > process_list;

  uniform_real_distribution<double> distrib_x(0, width);
  uniform_real_distribution<double> distrib_y(0, height);
  uniform_real_distribution<double> distrib_z(0, depth);

  const Vec3D<double> first_point(distrib_x(rng),
                                  distrib_y(rng),
                                  distrib_z(rng));

  process_list.push_back(first_point);
  positions.push_back(first_point[0]);
  positions.push_back(first_point[1]);
  positions.push_back(first_point[2]);
  const unsigned first_i = floor(first_point[0] / cell_size);
  const unsigned first_j = floor(first_point[1] / cell_size);
  const unsigned first_k = floor(first_point[2] / cell_size);
  grid(first_i, first_j, first_k) = first_point;
  grid_filled(first_i, first_j, first_k) = true;

  uniform_real_distribution<double> distrib_radius(min_dist, 2*min_dist);
  uniform_real_distribution<double> distrib_angle(0, 2 * M_PI);

  while (process_list.size() > 0) {
    // Pop a point.
    uniform_int_distribution<unsigned> distrib_idx(0, process_list.size()-1);
    const unsigned idx = distrib_idx(rng);
    const Vec3D<double> point = process_list[idx];
    for (unsigned i = idx; i < process_list.size()-1; ++i)
      process_list[i] = process_list[i+1];
    process_list.pop_back();

    for (unsigned unused_ctr = 0; unused_ctr < 500; ++unused_ctr) {
      const double radius = distrib_radius(rng);
      const double phi = distrib_angle(rng);
      const double rho = distrib_angle(rng);
      Vec3D<double> new_point(point[0] + radius * sin(phi) + cos(rho),
                              point[1] + radius * sin(phi) + sin(rho),
                              point[2] + radius * cos(phi));
      // Enforce PBCs.
      new_point[0] = fmod(new_point[0], width);
      new_point[0] = new_point[0] < 0 ? new_point[0] + width : new_point[0];
      new_point[1] = fmod(new_point[1], height);
      new_point[1] = new_point[1] < 0 ? new_point[1] + height : new_point[1];
      new_point[2] = fmod(new_point[2], depth);
      new_point[2] = new_point[2] < 0 ? new_point[2] + depth : new_point[2];
      // Get grid point.
      const unsigned i = floor(new_point[0] / cell_size);
      const unsigned j = floor(new_point[1] / cell_size);
      const unsigned k = floor(new_point[2] / cell_size);
      if (0 <= i && i < len_i &&
          0 <= j && j < len_j &&
          0 <= k && k < len_k &&
          !grid_filled(i,j,k) &&
          !in_neighborhood(new_point, i, j, k, len_i, len_j, len_k,
                           grid, grid_filled, min_dist_sq,
                           width, hx, height, hy, depth, hz,
                           true, true, true)) {
        process_list.push_back(new_point);
        positions.push_back(new_point[0]);
        positions.push_back(new_point[1]);
        positions.push_back(new_point[2]);
        grid(i, j, k) = new_point;
        grid_filled(i, j, k) = true;
      }
    }
  }
}

std::unique_ptr<Box> Box::random_box(double a, double b, double c,
                                     bool periodic_a,
                                     bool periodic_b,
                                     bool periodic_c,
                                     double min_dist,
                                     const std::string& atomtype,
                                     KIMNeigh neighmode,
                                     const std::string& name,
                                     std::mt19937& rng) {
  vector<double> positions;
  poisson_disk_sampling(a, b, c, min_dist, positions, rng);
  unsigned natoms = positions.size() / 3;
  auto coords = make_unique< Array2D<double> >(natoms, 3);
  auto types = make_unique< Array1DInit<string> >(natoms);
  for (unsigned i = 0; i < natoms; ++i) {
    (*coords)(i, 0) = positions[i*3 + 0];
    (*coords)(i, 1) = positions[i*3 + 1];
    (*coords)(i, 2) = positions[i*3 + 2];
    (*types)(i) = atomtype;
  }
  return make_unique<Box>(Vec3D<double>(a, 0, 0),
                          Vec3D<double>(0, b, 0),
                          Vec3D<double>(0, 0, c),
                          periodic_a, periodic_b, periodic_c,
                          move(coords), move(types),
                          neighmode, name);
}
