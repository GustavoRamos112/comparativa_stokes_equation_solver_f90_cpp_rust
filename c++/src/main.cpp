//* --------------------------------------------------------------------
//*  Main program - C++ translation of Rust version
//* --------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <array>
#include <tuple>
#include <algorithm>
#include <print>

namespace fs = std::filesystem;

constexpr size_t NNODES = 6;
constexpr size_t NQUAD = 3;

//* --------------------------------------------------------------------
//*  Index helpers (same as Rust phi_idx / psi_idx)
//* --------------------------------------------------------------------
inline size_t phi_idx(size_t it, size_t iquad, size_t iq, size_t k) {
  return ((it * NQUAD + iquad) * NNODES + iq) * 3 + k;
}
inline size_t psi_idx(size_t it, size_t iquad, size_t iq) {
  return (it * NQUAD + iquad) * NNODES + iq;
}

//* --------------------------------------------------------------------
//*  Bump struct (mirrors bump_struct::Bump)
//* --------------------------------------------------------------------
struct Bump {
  size_t nx, ny, mx, my;
  size_t maxnew, maxsec;
  size_t maxrow, nelemn, maxeqn, n_points;

  double anew, anext, aold, aprof;
  int iwrite;
  bool _long;
  size_t nband, neqn, nlband, nrow;
  size_t numnew, numsec;
  double reynld, rjpnew, rjpold, tolnew, tolsec;
  double xbleft, xbrite, xlngth, xprof, ylngth, ypert;

  std::vector<double> area;
  std::vector<double> dcda;
  std::vector<double> a;
  std::vector<double> f;
  std::vector<double> g;
  std::vector<double> gr;
  std::vector<int> iline;
  std::vector<int> iline_inv;
  std::vector<int> indx;
  std::vector<int> insc;
  std::vector<int> isotri;
  std::vector<size_t> node;
  std::vector<double> phi;
  std::vector<double> psi;
  std::vector<double> r;
  std::vector<double> res;
  std::vector<double> sens;
  std::vector<double> uprof;
  std::vector<double> xc;
  std::vector<double> xm;
  std::vector<double> yc;
  std::vector<double> ym;

  bool save_times;

  Bump(size_t nx_, size_t ny_, size_t mx_, size_t my_)
    : nx(nx_), ny(ny_), mx(mx_), my(my_),
      maxnew(4), maxsec(10),
      maxrow(27 * ny_),
      nelemn(2 * (nx_ - 1) * (ny_ - 1)),
      maxeqn(2 * mx_ * my_ + nx_ * ny_),
      n_points(mx_ * my_),
      anew(0.0), anext(0.3), aold(0.0), aprof(0.25),
      iwrite(10), _long(false),
      nband(0), neqn(0), nlband(0), nrow(0),
      numnew(0), numsec(0),
      reynld(1.0), rjpnew(0.0), rjpold(0.0),
      tolnew(0.0001), tolsec(0.0001),
      xbleft(1.0), xbrite(3.0), xlngth(10.0),
      xprof(4.0), ylngth(3.0), ypert(0.0),
      area(nelemn, 0.0),
      dcda(my_, 0.0),
      a(1, 0.0),
      f(1, 0.0),
      g(1, 0.0),
      gr(my_ * my_, 0.0),
      iline(my_, 0),
      iline_inv(),
      indx(n_points * 2, 0),
      insc(n_points, 0),
      isotri(nelemn, 0),
      node(nelemn * NNODES, 0),
      phi(nelemn * NQUAD * NNODES * 3, 0.0),
      psi(nelemn * NQUAD * NNODES, 0.0),
      r(my_, 0.0),
      res(1, 0.0),
      sens(1, 0.0),
      uprof(my_, 0.0),
      xc(n_points, 0.0),
      xm(nelemn * NQUAD, 0.0),
      yc(n_points, 0.0),
      ym(nelemn * NQUAD, 0.0),
      save_times(true)
  {}
};

// Forward declarations
void timestamp();
double bsp(
  size_t it, size_t iq, size_t id, const std::vector<size_t> &node,
  const std::vector<double> &xc, double xq,
  const std::vector<double> &yc, double yq
);
void daxpy_v(
  int n, double da, const std::vector<double> &dx, int incx,
  std::vector<double> &dy, int incy
);
void daxpy_m(
  int n, double da, std::vector<double> &abd, size_t stride,
  size_t col_x, size_t start_row_x, int incx,
  size_t col_y, size_t start_row_y, int incy
);
void daxpy(
  int n, double da, const std::vector<double> &abd_x, size_t stride,
  size_t col_x, size_t start_row_x, int incx,
  std::vector<double> &b_y, size_t start_idx_y, int incy
);
double ddot(
  int n, const std::vector<double> &abd, size_t stride, size_t col,
  size_t row_start, int incx, const std::vector<double> &b,
  size_t start_b, int incy
);
int dgbfa(
  std::vector<double> &abd, size_t lda, size_t n, size_t ml, size_t mu,
  std::vector<int> &ipvt
);
void dgbsl(
  std::vector<double> &abd, size_t lda, size_t n, size_t ml, size_t mu,
  const std::vector<int> &ipvt, std::vector<double> &b, int job
);
void dscal_v(int n, double sa, std::vector<double> &x, int incx);
void dscal_m(
  size_t n, double sa, std::vector<double> &abd, size_t stride,
  size_t col, size_t start_row, int incx
);
std::string file_name_inc(const std::string& file_name);
std::vector<double> getg(const std::vector<double> &f, const std::vector<int> &iline, size_t my);
void gram(Bump& bump);
size_t i4_modp(int i, int j);
size_t i4_wrap(size_t ival, size_t ilo, size_t ihi);
size_t idamax_v(int n, const std::vector<double> &dx, int incx);
size_t idamax_m(
  int n, const std::vector<double> &abd, size_t stride,
  size_t col, size_t start_row, int incx
);
int igetl(int i, const std::vector<int> &iline);
void linsys(Bump& bump, int itype);
void nstoke(Bump& bump);
std::tuple<double, double, double> qbf(
  double xq, double yq, size_t it, size_t inn,
  const std::vector<size_t> &node,
  const std::vector<double> &xc,
  const std::vector<double> &yc
);
double refbsp(double xq, double yq, size_t iq);
std::tuple<double, double, double> refqbf(
  double x, double y, size_t inn,
  double etax, double etay,
  double xix, double xiy
);
void resid(Bump& bump);
void setban(Bump& bump);
void setbas(Bump& bump);
void setgrd(int ibump, Bump& bump);
void setlin(Bump& bump);
void setqud(Bump& bump);
void setxy(Bump& bump);
double ubdry(size_t iuk, double yy);
double ubump(Bump& bump, size_t ip, size_t iqq, size_t it, size_t iukk);
std::tuple<std::array<double, 2>, std::array<double, 2>, std::array<double, 2>>
  _ubump_uval(
    Bump& bump, size_t it, double xix, double xiy, double xq, double yq,
    double _det, double etax, double etay
  );
std::tuple<std::array<double, 2>, std::array<double, 2>, std::array<double, 2>>
  uval(
    double etax, double etay, Bump& bump, size_t it,
    double xix, double xiy, double xq, double yq
  );
std::tuple<double, double, double, double, double>
  trans(
    size_t it, const std::vector<size_t> &node,
    const std::vector<double> &xc, double xq,
    const std::vector<double> &yc, double yq
  );
void uv_write(Bump& bump, const std::vector<double> &f, std::ofstream& file);
void xy_write(Bump& bump, std::ofstream& file);

//* --------------------------------------------------------------------
//*  Main program
//* --------------------------------------------------------------------
int main() {
  std::string uv_dir = "data/uv";
  std::string uv_file = uv_dir + "/uv_000.txt";
  std::string xy_dir = "data/xy";
  std::string xy_file = xy_dir + "\\xy_000.txt";

  size_t nx = 21;
  size_t ny = 7;
  size_t mx = 2 * nx - 1;
  size_t my = 2 * ny - 1;

  Bump bump(nx, ny, mx, my);

  if (!fs::exists(uv_dir)) {
    fs::create_directories(uv_dir);
  }
  if (!fs::exists(xy_dir)) {
    fs::create_directories(xy_dir);
  }

  auto inicio = std::chrono::high_resolution_clock::now();
  timestamp();
  std::println();
  std::println("BUMP");
  std::println("  C++ version");
  std::println("  Control problem for channel flow over a bump.");
  std::println();
  std::println("  The bump will be generated with a height of {}", bump.aprof);
  std::println();
  std::println("  NX = {}", nx);
  std::println("  NY = {}", ny);
  std::println("  Number of elements = {}", bump.nelemn);
  std::println("  Reynolds number =  {}", bump.reynld);
  std::println("  Secant tolerance = {}", bump.tolsec);
  std::println("  Newton tolerance = {}", bump.tolnew);

  int ibump = 2;
  setgrd(ibump, bump);

  bump.f.resize(bump.neqn, 0.0);
  bump.g.resize(bump.neqn, 0.0);
  bump.res.resize(bump.neqn, 0.0);
  bump.sens.resize(bump.neqn, 0.0);

  bump.ypert = bump.aprof;
  setxy(bump);
  setqud(bump);
  setbas(bump);
  setlin(bump);
  setban(bump);

  bump.a.resize(bump.nrow * bump.neqn, 0.0);

  nstoke(bump);
  resid(bump);

  bump.uprof = getg(bump.g, bump.iline, bump.my);

  if (1 <= bump.iwrite) {
    std::println();
    std::println("Velocity profile:");
    std::println();
    for (size_t i = 0; i < bump.my; i++) {
      std::print("{}\t", bump.uprof[i]);
      if ((i + 1) % 5 == 0) {
        std::println();
      }
    }
    std::println();
  }

  gram(bump);

  xy_file = file_name_inc(xy_file);
  {
    std::ofstream f_xy(xy_file);
    xy_write(bump, f_xy);
  }

  uv_file = file_name_inc(uv_file);
  {
    std::ofstream f_uv(uv_file);
    uv_write(bump, bump.f, f_uv);
  }

  for (size_t i = 0; i < bump.neqn; i++) {
    bump.g[i] = 0.0;
  }

  size_t iter = 1;
  while (true) {
    if (iter > bump.maxsec) {
      break;
    }
    std::println();
    std::println("Secant iteration {}", iter);

    bump.numsec += 1;

    bump.ypert = bump.anew;
    setxy(bump);
    setqud(bump);
    setbas(bump);
    nstoke(bump);

    bump.uprof = getg(bump.g, bump.iline, bump.my);

    if (1 <= bump.iwrite) {
      std::println();
      std::println("Velocity profile:");
      std::println();
      for (size_t i = 0; i < bump.my; i++) {
        std::print("{}\t", bump.uprof[i]);
        if ((i + 1) % 5 == 0) {
          std::println();
        }
      }
      std::println();
    }

    linsys(bump, -2);

    bump.dcda = getg(bump.sens, bump.iline, bump.my);

    if (2 <= bump.iwrite) {
      std::println();
      std::println("Sensitivities:");
      std::println();
      for (size_t i = 0; i < bump.my; i++) {
        std::print("{}\t", bump.dcda[i]);
        if ((i + 1) % 5 == 0) {
          std::println();
        }
      }
      std::println();
    }

    bump.rjpnew = 0.0;
    for (size_t i = 0; i < bump.my; i++) {
      double temp = -bump.r[i];
      for (size_t j = 0; j < bump.my; j++) {
        temp += bump.gr[i * bump.my + j] * bump.uprof[j];
      }
      bump.rjpnew += 2.0 * bump.dcda[i] * temp;
    }

    xy_file = file_name_inc(xy_file);
    {
      std::ofstream f_xy(xy_file);
      xy_write(bump, f_xy);
    }

    uv_file = file_name_inc(uv_file);
    {
      std::ofstream f_uv(uv_file);
      uv_write(bump, bump.f, f_uv);
    }

    std::println();
    std::println("  Parameter = {}, J prime = {}", bump.anew, bump.rjpnew);

    if (1 < iter) {
      double denom = bump.rjpnew - bump.rjpold;
      bump.anext = (std::abs(denom) > 1e-30)
        ? bump.aold - bump.rjpold * (bump.anew - bump.aold) / denom
        : bump.anew;
    }

    bump.aold = bump.anew;
    bump.anew = bump.anext;
    bump.rjpold = bump.rjpnew;

    double test = (bump.anew != 0.0)
      ? std::abs(bump.anew - bump.aold) / bump.anew
      : 0.0;

    std::println("  New value of parameter = {}", bump.anew);
    std::println("  Convergence test = {}", test);

    if (std::abs(bump.anew - bump.aold) <= std::abs(bump.anew) * bump.tolsec && 1 < iter) {
      std::println("Secant iteration converged.");
      break;
    }
    iter += 1;
  }
  if (bump.maxsec < iter) {
    std::println("  Secant iteration failed to converge.");
  }

  auto duracion = std::chrono::high_resolution_clock::now() - inicio;
  double duracion_sec = std::chrono::duration<double>(duracion).count();

  std::println();
  std::println("  Total execution time = {} s", duracion_sec);
  std::println("  Number of secant steps = {}", bump.numsec);
  std::println("  Number of Newton steps = {}", bump.numnew);
  std::println();
  std::println("BUMP:");
  std::println("  Normal end of execution.");
  std::println();
  timestamp();

  if (bump.save_times) {
    std::ofstream file("times.txt", std::ios::app);
    file << duracion_sec << "\n";
    file.flush();
  }

  return 0;
}

//* --------------------------------------------------------------------
//*  TIMESTAMP
//* --------------------------------------------------------------------
void timestamp()
{
  auto ahora = std::chrono::system_clock::now();
  std::println("{0:%Y-%m-%d %H:%M:%S}", ahora);
}

//* --------------------------------------------------------------------
//*  BSP - linear basis function for pressure
//* --------------------------------------------------------------------
double bsp(
  size_t it, size_t iq, size_t id, const std::vector<size_t> &node,
  const std::vector<double> &xc, double xq,
  const std::vector<double> &yc, double yq
) {
  size_t l1 = iq;
  size_t l2 = i4_wrap(iq + 1, 0, 2);
  size_t l3 = i4_wrap(iq + 2, 0, 2);

  size_t g1 = node[it * NNODES + l1];
  size_t g2 = node[it * NNODES + l2];
  size_t g3 = node[it * NNODES + l3];

  double d = (xc[g2] - xc[g1]) * (yc[g3] - yc[g1])
           - (xc[g3] - xc[g1]) * (yc[g2] - yc[g1]);

  if (id == 0) {
    return 1.0 + ((yc[g2] - yc[g3]) * (xq - xc[g1]) + (xc[g3] - xc[g2]) * (yq - yc[g1])) / d;
  } else if (id == 1) {
    return (yc[g2] - yc[g3]) / d;
  } else if (id == 2) {
    return (xc[g3] - xc[g2]) / d;
  } else {
    std::println();
    std::println("BSP - Fatal error!");
    std::println("  Illegal local index value for linear basis.");
    std::println("  Legal values are 1, 2 or 3.");
    std::println("  The input value was ID = {}", id + 1);
    std::exit(1);
  }
}

//* --------------------------------------------------------------------
//*  DAXPY - constant times a vector plus a vector
//* --------------------------------------------------------------------
void daxpy_v(
  int n, double da, const std::vector<double> &dx, int incx,
  std::vector<double> &dy, int incy
) {
  if (n <= 0 || da == 0.0) return;
  int ix, iy;

  if (incx != 1 || incy != 1) {
    ix = (incx >= 0) ? 0 : (-n + 1) * incx;
    iy = (incy >= 0) ? 0 : (-n + 1) * incy;
    for (int i = 0; i < n; i++) {
      dy[iy] += da * dx[ix];
      ix += incx;
      iy += incy;
    }
  } else {
    int m = n % 4;
    for (int i = 0; i < m; i++) {
      dy[i] += da * dx[i];
    }
    for (int i = m; i < n; i += 4) {
      dy[i]     += da * dx[i];
      dy[i + 1] += da * dx[i + 1];
      dy[i + 2] += da * dx[i + 2];
      dy[i + 3] += da * dx[i + 3];
    }
  }
}

void daxpy_m(
  int n, double da, std::vector<double> &abd, size_t stride,
  size_t col_x, size_t start_row_x, int incx,
  size_t col_y, size_t start_row_y, int incy
) {
  if (n <= 0 || da == 0.0) return;
  for (int i = 0; i < n; i++) {
    size_t idx_x = (start_row_x + (i * incx)) * stride + col_x;
    size_t idx_y = (start_row_y + (i * incy)) * stride + col_y;
    abd[idx_y] += da * abd[idx_x];
  }
}

void daxpy(
  int n, double da, const std::vector<double> &abd_x, size_t stride,
  size_t col_x, size_t start_row_x, int incx,
  std::vector<double> &b_y, size_t start_idx_y, int incy
) {
  if (n <= 0 || da == 0.0) return;
  for (int i = 0; i < n; i++) {
    size_t idx_x = (start_row_x + (i * incx)) * stride + col_x;
    size_t idx_y = start_idx_y + (i * incy);
    b_y[idx_y] += da * abd_x[idx_x];
  }
}

//* --------------------------------------------------------------------
//*  DDOT - dot product of two vectors
//* --------------------------------------------------------------------
double ddot(
  int n, const std::vector<double> &abd, size_t stride, size_t col,
  size_t row_start, int incx, const std::vector<double> &b,
  size_t start_b, int incy
) {
  if (n <= 0) return 0.0;
  double dtemp = 0.0;
  if (incx != 1 || incy != 1) {
    int ix = (int)row_start;
    int iy = (int)start_b;
    for (int i = 0; i < n; i++) {
      dtemp += abd[(size_t)ix * stride + col] * b[(size_t)iy];
      ix += incx;
      iy += incy;
    }
  } else {
    int m = n % 5;
    for (int i = 0; i < m; i++) {
      dtemp += abd[(row_start + i) * stride + col] * b[start_b + i];
    }
    for (int i = m; i < n; i += 5) {
      dtemp += abd[(row_start + i) * stride + col] * b[start_b + i]
             + abd[(row_start + i + 1) * stride + col] * b[start_b + i + 1]
             + abd[(row_start + i + 2) * stride + col] * b[start_b + i + 2]
             + abd[(row_start + i + 3) * stride + col] * b[start_b + i + 3]
             + abd[(row_start + i + 4) * stride + col] * b[start_b + i + 4];
    }
  }
  return dtemp;
}

//* --------------------------------------------------------------------
//*  DGBFA - factor a real band matrix by elimination
//* --------------------------------------------------------------------
int dgbfa(
  std::vector<double> &abd, size_t lda, size_t n, size_t ml, size_t mu,
  std::vector<int> &ipvt
) {
  size_t m = ml + mu + 1;
  int info = 0;
  size_t j0 = mu + 1;
  size_t j1 = std::min(n, m) - 1;

  for (size_t jz = j0; jz < j1; jz++) {
    size_t i0 = m - jz;
    for (size_t i = i0; i < ml; i++) {
      abd[i * n + jz] = 0.0;
    }
  }

  size_t jz = j1;
  size_t ju = 0;

  for (size_t k = 0; k < n - 1; k++) {
    jz += 1;
    if (jz < n) {
      for (size_t i = 0; i < ml; i++) {
        abd[i * n + jz] = 0.0;
      }
    }

    size_t lm = std::min(ml, n - k - 1);
    size_t l = idamax_m((int)(lm + 1), abd, n, k, m - 1, 1) + (m - 1);
    ipvt[k] = (int)(l + k - (m - 1));

    if (abd[l * n + k] == 0.0) {
      info = (int)k;
    } else {
      if (l != m - 1) {
        double t = abd[l * n + k];
        abd[l * n + k] = abd[(m - 1) * n + k];
        abd[(m - 1) * n + k] = t;
      }

      double t = -1.0 / abd[(m - 1) * n + k];
      dscal_m(lm, t, abd, n, k, m, 1);

      ju = std::max(ju, mu + (size_t)ipvt[k]);
      ju = std::min(ju, n - 1);
      size_t mm = m - 1;
      size_t l_cur = l;

      for (size_t j = k + 1; j < ju + 1; j++) {
        l_cur = l_cur - 1;
        mm = mm - 1;
        double t2 = abd[l_cur * n + j];
        if (l_cur != mm) {
          abd[l_cur * n + j] = abd[mm * n + j];
          abd[mm * n + j] = t2;
        }
        daxpy_m((int)lm, t2, abd, n, k, m, 1, j, mm + 1, 1);
      }
    }
  }

  ipvt[n - 1] = (int)(n - 1);

  if (abd[(m - 1) * n + (n - 1)] == 0.0) {
    info = (int)(n - 1);
  }

  return info;
}

//* --------------------------------------------------------------------
//*  DGBSL - solve a real banded system factored by DGBFA
//* --------------------------------------------------------------------
void dgbsl(
  std::vector<double> &abd, size_t lda, size_t n, size_t ml, size_t mu,
  const std::vector<int> &ipvt, std::vector<double> &b, int job
) {
  size_t m = mu + ml + 1;
  if (job == 0) {
    if (0 < ml) {
      for (size_t k = 0; k < n - 1; k++) {
        size_t lm = std::min(ml, n - k - 1);
        size_t l = (size_t)ipvt[k];
        double t = b[l];
        if (l != k) {
          b[l] = b[k];
          b[k] = t;
        }
        daxpy((int)lm, t, abd, n, k, m, 1, b, k + 1, 1);
      }
    }

    for (size_t k = n; k-- > 0; ) {
      b[k] /= abd[(m - 1) * n + k];
      size_t lm = std::min(k, m - 1);
      size_t la = m - 1 - lm;
      size_t lb = k - lm;
      double t = -b[k];
      daxpy((int)lm, t, abd, n, k, la, 1, b, lb, 1);
    }
  } else {
    for (size_t k = 0; k < n; k++) {
      size_t lm = std::min(k, m - 1);
      size_t la = m - 1 - lm;
      size_t lb = k - lm;
      double t = ddot((int)lm, abd, n, k, la, 1, b, lb, 1);
      b[k] = (b[k] - t) / abd[(m - 1) * n + k];
    }

    if (0 < ml) {
      for (size_t k = n - 1; k-- > 0; ) {
        size_t lm = std::min(ml, n - k - 1);
        b[k] += ddot((int)lm, abd, n, k, m, 1, b, k + 1, 1);
        size_t l = (size_t)ipvt[k];
        if (l != k) {
          double t = b[l];
          b[l] = b[k];
          b[k] = t;
        }
      }
    }
  }
}

//* --------------------------------------------------------------------
//*  DSCAL - scale a vector by a constant
//* --------------------------------------------------------------------
void dscal_v(int n, double sa, std::vector<double> &x, int incx) {
  if (n <= 0) return;
  if (incx == 1) {
    int m = n % 5;
    for (int i = 0; i < m; i++) {
      x[i] *= sa;
    }
    for (int i = m; i < n; i += 5) {
      x[i]     *= sa;
      x[i + 1] *= sa;
      x[i + 2] *= sa;
      x[i + 3] *= sa;
      x[i + 4] *= sa;
    }
  } else {
    int ix = (incx >= 0) ? 0 : (-n + 1) * incx;
    for (int i = 0; i < n; i++) {
      x[ix] *= sa;
      ix += incx;
    }
  }
}

void dscal_m(
  size_t n, double sa, std::vector<double> &abd, size_t stride,
  size_t col, size_t start_row, int incx
) {
  if (n == 0) return;
  for (size_t i = 0; i < n; i++) {
    abd[(start_row + (i * incx)) * stride + col] *= sa;
  }
}

//* --------------------------------------------------------------------
//*  FILE_NAME_INC - increment a partially numeric filename
//* --------------------------------------------------------------------
std::string file_name_inc(const std::string& file_name) {
  if (file_name.empty()) {
    std::println();
    std::println("FILE_NAME_INC - Fatal error!");
    std::println("The input string is empty.");
    std::exit(1);
  }
  std::string s = file_name;
  int change = 0;
  for (int i = (int)s.length() - 1; i >= 0; i--) {
    if (s[i] >= '0' && s[i] <= '9') {
      change += 1;
      char digit = s[i] - '0';
      digit += 1;
      if (digit == 10) {
        s[i] = '0';
      } else {
        s[i] = (char)(digit + '0');
        return s;
      }
    }
  }
  if (change == 0) {
    return std::string();
  } else {
    return s;
  }
}

//* --------------------------------------------------------------------
//*  GETG - extract values of a quantity along the profile line
//* --------------------------------------------------------------------
std::vector<double> getg(const std::vector<double> &f, const std::vector<int> &iline, size_t my) {
  std::vector<double> u(my, 0.0);
  for (size_t i = 0; i < my; i++) {
    int j = iline[i];
    u[i] = (j <= 0) ? 0.0 : f[(size_t)(j - 1)];
  }
  return u;
}

//* --------------------------------------------------------------------
//*  GRAM - compute the Gram matrix and R vector
//* --------------------------------------------------------------------
void gram(Bump& bump) {
  double wt[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  double yq_gauss[3] = {-0.7745966692, 0.0, 0.7745966692};

  for (size_t it = 0; it < bump.my; it++) {
    bump.r[it] = 0.0;
  }
  for (size_t it = 0; it < bump.my; it++) {
    for (size_t j = 0; j < bump.my; j++) {
      bump.gr[it * bump.my + j] = 0.0;
    }
  }

  for (size_t it = 0; it < bump.nelemn; it++) {
    size_t k = bump.node[it * NNODES];
    size_t kk = bump.node[it * NNODES + 1];

    if (std::abs(bump.xc[k] - bump.xprof) > 1.0e-4 ||
        std::abs(bump.xc[kk] - bump.xprof) > 1.0e-4) {
      continue;
    }

    for (size_t iquad = 0; iquad < 3; iquad++) {
      double bma2 = (bump.yc[kk] - bump.yc[k]) / 2.0;
      double ar = bma2 * wt[iquad];
      double x = bump.xprof;
      double y = bump.yc[k] + bma2 * (yq_gauss[iquad] + 1.0);

      double uiqdpt = 0.0;
      for (size_t iq = 0; iq < NNODES; iq++) {
        if (iq == 0 || iq == 1 || iq == 3) {
          auto [bb, _bx, _by] = qbf(x, y, it, iq, bump.node, bump.xc, bump.yc);
          size_t ip = bump.node[it * NNODES + iq];
          int iun = bump.indx[ip * 2];
          if (0 < iun) {
            int ii = igetl(iun, bump.iline);
            uiqdpt += bb * bump.uprof[(size_t)(ii - 1)];
          } else if (iun == -1) {
            double ubc = ubdry(1, bump.yc[ip]);
            uiqdpt += bb * ubc;
          }
        }
      }

      for (size_t iq = 0; iq < NNODES; iq++) {
        if (iq == 0 || iq == 1 || iq == 3) {
          size_t ip = bump.node[it * NNODES + iq];
          auto [bb, _bx, _by] = qbf(x, y, it, iq, bump.node, bump.xc, bump.yc);
          int i_val = bump.indx[ip * 2];
          if (0 < i_val) {
            int ii = igetl(i_val, bump.iline);
            bump.r[(size_t)(ii - 1)] += bb * uiqdpt * ar;
            for (size_t iqq = 0; iqq < NNODES; iqq++) {
              if (iqq == 0 || iqq == 1 || iqq == 3) {
                size_t ipp = bump.node[it * NNODES + iqq];
                auto [bbb, _bbx_, _bby_] = qbf(x, y, it, iqq, bump.node, bump.xc, bump.yc);
                int j_val = bump.indx[ipp * 2];
                if (j_val != 0) {
                  int jj = igetl(j_val, bump.iline);
                  bump.gr[((size_t)(ii - 1)) * bump.my + ((size_t)(jj - 1))] += bb * bbb * ar;
                }
              }
            }
          }
        }
      }
    }
  }

  if (3 <= bump.iwrite) {
    std::println();
    std::println("Gram matrix:");
    std::println();
    for (size_t i = 0; i < bump.my; i++) {
      for (size_t j = 0; j < bump.my; j++) {
        std::println("{} {} {}", i + 1, j + 1, bump.gr[i * bump.my + j]);
      }
    }
    std::println();
    std::println("R vector:");
    std::println();
    for (size_t i = 0; i < bump.my; i++) {
      std::println("{}", bump.r[i]);
    }
  }
}

//* --------------------------------------------------------------------
//*  I4_MODP - nonnegative remainder of integer division
//* --------------------------------------------------------------------
size_t i4_modp(int i, int j) {
  if (j == 0) {
    std::println();
    std::println("I4_MODP - Fatal error!");
    std::println("  Illegal divisor J = {}", j);
    std::exit(1);
  }
  int value = i % j;
  if (value < 0) {
    value += std::abs(j);
  }
  return (size_t)value;
}

//* --------------------------------------------------------------------
//*  I4_WRAP - force an integer between limits by wrapping
//* --------------------------------------------------------------------
size_t i4_wrap(size_t ival, size_t ilo, size_t ihi) {
  size_t jlo = std::min(ilo, ihi);
  size_t jhi = std::max(ilo, ihi);
  int wide = (int)(jhi - jlo + 1);
  if (wide == 1) {
    return jlo;
  }
  return jlo + i4_modp((int)(ival - jlo), wide);
}

//* --------------------------------------------------------------------
//*  IDAMAX - index of the vector element of maximum absolute value
//* --------------------------------------------------------------------
size_t idamax_v(int n, const std::vector<double> &dx, int incx) {
  if (n < 1 || incx <= 0) {
    return 0;
  }
  if (n == 1) {
    return 0;
  }
  if (incx == 1) {
    double dmax = std::abs(dx[0]);
    size_t idamax_val = 0;
    for (int i = 1; i < n; i++) {
      double val = std::abs(dx[i]);
      if (dmax < val) {
        idamax_val = (size_t)i;
        dmax = val;
      }
    }
    return idamax_val;
  } else {
    int ix = 0;
    double dmax = std::abs(dx[0]);
    size_t idamax_val = 0;
    ix += incx;
    for (int i = 1; i < n; i++) {
      double val = std::abs(dx[ix]);
      if (dmax < val) {
        idamax_val = (size_t)i;
        dmax = val;
      }
      ix += incx;
    }
    return idamax_val;
  }
}

size_t idamax_m(
  int n, const std::vector<double> &abd, size_t stride,
  size_t col, size_t start_row, int incx
) {
  if (n < 1 || incx <= 0) {
    return 0;
  }
  double dmax = std::abs(abd[start_row * stride + col]);
  size_t idamax_val = 0;
  for (int i = 1; i < n; i++) {
    double current_val = std::abs(abd[(start_row + (i * incx)) * stride + col]);
    if (current_val > dmax) {
      idamax_val = (size_t)i;
      dmax = current_val;
    }
  }
  return idamax_val;
}

//* --------------------------------------------------------------------
//*  IGETL - get the local unknown number along the profile line
//* --------------------------------------------------------------------
int igetl(int i, const std::vector<int> &iline) {
  for (size_t j = 0; j < iline.size(); j++) {
    if (iline[j] == i) {
      return (int)(j + 1);
    }
  }
  return -1;
}

//* --------------------------------------------------------------------
//*  LINSYS - solve the linearized Navier Stokes equation
//* --------------------------------------------------------------------
void linsys(Bump& bump, int itype) {
  int ioff = (int)(bump.nlband + bump.nlband + 1);
  double visc = 1.0 / bump.reynld;
  bool use_sens = (itype == -2);
  size_t neqn = bump.neqn;
  size_t nrow = bump.nrow;

  std::vector<double> rhs(neqn, 0.0);
  std::fill(bump.a.begin(), bump.a.begin() + nrow * neqn, 0.0);

  for (size_t it = 0; it < bump.nelemn; it++) {
    size_t elem_nodes[6] = {
      bump.node[it * NNODES],
      bump.node[it * NNODES + 1],
      bump.node[it * NNODES + 2],
      bump.node[it * NNODES + 3],
      bump.node[it * NNODES + 4],
      bump.node[it * NNODES + 5]
    };
    int elem_indx_u[6] = {
      bump.indx[elem_nodes[0] * 2],
      bump.indx[elem_nodes[1] * 2],
      bump.indx[elem_nodes[2] * 2],
      bump.indx[elem_nodes[3] * 2],
      bump.indx[elem_nodes[4] * 2],
      bump.indx[elem_nodes[5] * 2]
    };
    int elem_indx_v[6] = {
      bump.indx[elem_nodes[0] * 2 + 1],
      bump.indx[elem_nodes[1] * 2 + 1],
      bump.indx[elem_nodes[2] * 2 + 1],
      bump.indx[elem_nodes[3] * 2 + 1],
      bump.indx[elem_nodes[4] * 2 + 1],
      bump.indx[elem_nodes[5] * 2 + 1]
    };
    int elem_insc[6] = {
      bump.insc[elem_nodes[0]],
      bump.insc[elem_nodes[1]],
      bump.insc[elem_nodes[2]],
      bump.insc[elem_nodes[3]],
      bump.insc[elem_nodes[4]],
      bump.insc[elem_nodes[5]]
    };
    int isotri = bump.isotri[it];
    double area = bump.area[it];
    double ar = area / 3.0;

    for (size_t iquad = 0; iquad < NQUAD; iquad++) {
      double yq = bump.ym[it * NQUAD + iquad];
      double xq = bump.xm[it * NQUAD + iquad];

      double det, etax, etay, xix, xiy;

      if (isotri == 1) {
        std::tie(det, etax, etay, xix, xiy) = trans(it, bump.node, bump.xc, xq, bump.yc, yq);
        ar = det * area / 3.0;
      } else {
        det = 0.0; etax = 0.0; etay = 0.0; xix = 0.0; xiy = 0.0;
      }

      auto [un, unx, uny] = uval(etax, etay, bump, it, xix, xiy, xq, yq);

      std::array<double, NNODES> phi_bb;
      std::array<double, NNODES> phi_bx;
      std::array<double, NNODES> phi_by;
      std::array<double, NNODES> psi_bbl;
      for (size_t i = 0; i < NNODES; i++) {
        phi_bb[i] = bump.phi[phi_idx(it, iquad, i, 0)];
        phi_bx[i] = bump.phi[phi_idx(it, iquad, i, 1)];
        phi_by[i] = bump.phi[phi_idx(it, iquad, i, 2)];
        psi_bbl[i] = bump.psi[psi_idx(it, iquad, i)];
      }

      for (size_t iq = 0; iq < NNODES; iq++) {
        double bb = phi_bb[iq];
        double bx = phi_bx[iq];
        double by = phi_by[iq];
        double bbl = psi_bbl[iq];
        int ihor = elem_indx_u[iq];
        int iver = elem_indx_v[iq];
        int iprs = elem_insc[iq];

        if (0 < ihor) {
          rhs[(size_t)(ihor - 1)] += ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
        }
        if (0 < iver) {
          rhs[(size_t)(iver - 1)] += ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
        }

        for (size_t iqq = 0; iqq < NNODES; iqq++) {
          double bbb = phi_bb[iqq];
          double bbx = phi_bx[iqq];
          double bby = phi_by[iqq];
          double bbbl = psi_bbl[iqq];
          int ju = elem_indx_u[iqq];
          int jv = elem_indx_v[iqq];
          int jp = elem_insc[iqq];
          size_t ipp = elem_nodes[iqq];

          if (0 < ju) {
            if (0 < ihor) {
              int iuse = ihor - ju + ioff;
              size_t row_off = (size_t)(iuse - 1) * neqn;
              size_t col = (size_t)(ju - 1);
              bump.a[row_off + col] += ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
            }
            if (0 < iver) {
              int iuse = iver - ju + ioff;
              size_t row_off = (size_t)(iuse - 1) * neqn;
              size_t col = (size_t)(ju - 1);
              bump.a[row_off + col] += ar * bb * bbb * unx[1];
            }
            if (0 < iprs) {
              int iuse = iprs - ju + ioff;
              size_t row_off = (size_t)(iuse - 1) * neqn;
              size_t col = (size_t)(ju - 1);
              bump.a[row_off + col] += ar * bbx * bbl;
            }
          } else if (ju == itype) {
            double ubc = (ju == -1) ? ubdry(1, bump.yc[ipp]) : ubump(bump, ipp, iqq, it, 1);

            if (0 < ihor) {
              double aij = ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
              rhs[(size_t)(ihor - 1)] -= ubc * aij;
            }
            if (0 < iver) {
              double aij = ar * bb * bbb * unx[1];
              rhs[(size_t)(iver - 1)] -= ubc * aij;
            }
            if (0 < iprs) {
              double aij = ar * bbx * bbl;
              rhs[(size_t)(iprs - 1)] -= ubc * aij;
            }
          }

          if (0 < jv) {
            if (0 < ihor) {
              int iuse = ihor - jv + ioff;
              size_t row_off = (size_t)(iuse - 1) * neqn;
              size_t col = (size_t)(jv - 1);
              bump.a[row_off + col] += ar * bb * bbb * uny[0];
            }
            if (0 < iver) {
              int iuse = iver - jv + ioff;
              size_t row_off = (size_t)(iuse - 1) * neqn;
              size_t col = (size_t)(jv - 1);
              bump.a[row_off + col] += ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
            }
            if (0 < iprs) {
              int iuse = iprs - jv + ioff;
              size_t row_off = (size_t)(iuse - 1) * neqn;
              size_t col = (size_t)(jv - 1);
              bump.a[row_off + col] += ar * bby * bbl;
            }
          } else if (jv == itype) {
            double ubc = (jv == -1) ? ubdry(2, bump.yc[ipp]) : ubump(bump, ipp, iqq, it, 2);

            if (0 < ihor) {
              double aij = ar * bb * bbb * uny[0];
              rhs[(size_t)(ihor - 1)] -= ubc * aij;
            }
            if (0 < iver) {
              double aij = ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
              rhs[(size_t)(iver - 1)] -= ubc * aij;
            }
            if (0 < iprs) {
              double aij = ar * bby * bbl;
              rhs[(size_t)(iprs - 1)] -= ubc * aij;
            }
          }

          if (0 < jp) {
            if (0 < ihor) {
              int iuse = ihor - jp + ioff;
              size_t row_off = (size_t)(iuse - 1) * neqn;
              size_t col = (size_t)(jp - 1);
              bump.a[row_off + col] -= ar * bx * bbbl;
            }
            if (0 < iver) {
              int iuse = iver - jp + ioff;
              size_t row_off = (size_t)(iuse - 1) * neqn;
              size_t col = (size_t)(jp - 1);
              bump.a[row_off + col] -= ar * by * bbbl;
            }
          }
        }
      }
    }
  }

  rhs[neqn - 1] = 0.0;
  for (size_t j_1based = neqn - bump.nlband; j_1based < neqn; j_1based++) {
    size_t j = j_1based - 1;
    int i = (int)(neqn - j_1based + (size_t)ioff);
    bump.a[(size_t)(i - 1) * neqn + j] = 0.0;
  }
  bump.a[(size_t)(ioff - 1) * neqn + (neqn - 1)] = 1.0;

  std::vector<int> ipvt(neqn, 0);
  int info = dgbfa(bump.a, bump.maxrow, neqn, bump.nlband, bump.nlband, ipvt);

  if (info != 0) {
    std::println();
    std::println("LINSYS - fatal error!");
    std::println("DGBFA returns INFO = {}", info + 1);
    std::exit(1);
  }

  std::vector<double> &actual_rhs = use_sens ? bump.sens : bump.f;
  actual_rhs = rhs;
  dgbsl(bump.a, bump.maxrow, neqn, bump.nlband, bump.nlband, ipvt, actual_rhs, 0);
}

//* --------------------------------------------------------------------
//*  NSTOKE - solve Navier Stokes using Newton iteration
//* --------------------------------------------------------------------
void nstoke(Bump& bump) {
  for (size_t iter = 0; iter < bump.maxnew; iter++) {
    bump.numnew += 1;

    linsys(bump, -1);

    for (size_t i = 0; i < bump.neqn; i++) {
      bump.g[i] -= bump.f[i];
    }

    size_t imax = idamax_v((int)bump.neqn, bump.g, 1);
    double diff = std::abs(bump.g[imax]);
    std::println("NSTOKE: Iteration {}, MaxNorm(diff) = {}", iter + 1, diff);

    for (size_t i = 0; i < bump.neqn; i++) {
      bump.g[i] = bump.f[i];
    }

    if (diff <= bump.tolnew) {
      std::println("NSTOKE converged.");
      return;
    }

    if (iter == bump.maxnew - 1) {
      std::println("NSTOKE failed!");
      std::exit(1);
    }
  }
}

//* --------------------------------------------------------------------
//*  QBF - evaluate quadratic basis functions
//* --------------------------------------------------------------------
std::tuple<double, double, double> qbf(
  double xq, double yq, size_t it, size_t inn,
  const std::vector<size_t> &node,
  const std::vector<double> &xc,
  const std::vector<double> &yc
) {
  size_t in1, in2, in3;
  size_t i1, i2, i3;
  double d, c, t, s, bb, bx, by;

  if (inn <= 2) {
    in1 = inn;
    in2 = (inn + 1) % 3;
    in3 = (inn + 2) % 3;
    i1 = node[it * NNODES + in1];
    i2 = node[it * NNODES + in2];
    i3 = node[it * NNODES + in3];
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    t = 1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d;
    bb = t * (2.0 * t - 1.0);
    bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d;
    by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d;
  } else {
    size_t inn_local = inn - 3;
    in1 = inn_local;
    in2 = (inn_local + 1) % 3;
    in3 = (inn_local + 2) % 3;
    i1 = node[it * NNODES + in1];
    i2 = node[it * NNODES + in2];
    i3 = node[it * NNODES + in3];
    size_t j1 = i2;
    size_t j2 = i3;
    size_t j3 = i1;
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    c = (xc[j2] - xc[j1]) * (yc[j3] - yc[j1]) - (xc[j3] - xc[j1]) * (yc[j2] - yc[j1]);
    t = 1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d;
    s = 1.0 + ((yc[j2] - yc[j3]) * (xq - xc[j1]) + (xc[j3] - xc[j2]) * (yq - yc[j1])) / c;
    bb = 4.0 * s * t;
    bx = 4.0 * (t * (yc[j2] - yc[j3]) / c + s * (yc[i2] - yc[i3]) / d);
    by = 4.0 * (t * (xc[j3] - xc[j2]) / c + s * (xc[i3] - xc[i2]) / d);
  }
  return {bb, bx, by};
}

//* --------------------------------------------------------------------
//*  REFBSP - evaluate linear basis functions in a reference triangle
//* --------------------------------------------------------------------
double refbsp(double xq, double yq, size_t iq) {
  if (iq == 0) return 1.0 - xq;
  else if (iq == 1) return yq;
  else if (iq == 2) return xq - yq;
  return 0.0;
}

//* --------------------------------------------------------------------
//*  REFQBF - evaluate quadratic basis functions on reference triangle
//* --------------------------------------------------------------------
std::tuple<double, double, double> refqbf(
  double x, double y, size_t inn,
  double etax, double etay,
  double xix, double xiy
) {
  double bb, tbx, tby, bx, by;

  if (inn == 0) {
    bb = 1.0 - 3.0 * x + 2.0 * x * x;
    tbx = -3.0 + 4.0 * x;
    tby = 0.0;
  } else if (inn == 1) {
    bb = -y + 2.0 * y * y;
    tbx = 0.0;
    tby = -1.0 + 4.0 * y;
  } else if (inn == 2) {
    bb = -x + 2.0 * x * x + y - 4.0 * x * y + 2.0 * y * y;
    tbx = -1.0 + 4.0 * x - 4.0 * y;
    tby = 1.0 - 4.0 * x + 4.0 * y;
  } else if (inn == 3) {
    bb = 4.0 * y - 4.0 * x * y;
    tbx = -4.0 * y;
    tby = 4.0 - 4.0 * x;
  } else if (inn == 4) {
    bb = 4.0 * x * y - 4.0 * y * y;
    tbx = 4.0 * y;
    tby = 4.0 * x - 8.0 * y;
  } else if (inn == 5) {
    bb = 4.0 * x - 4.0 * x * x - 4.0 * y + 4.0 * x * y;
    tbx = 4.0 - 8.0 * x + 4.0 * y;
    tby = -4.0 + 4.0 * x;
  } else {
    std::println("REFQBF - Illegal value of IN = {}", inn + 1);
    std::exit(1);
  }

  bx = tbx * xix + tby * etax;
  by = tbx * xiy + tby * etay;

  return {bb, bx, by};
}

//* --------------------------------------------------------------------
//*  RESID - compute the residual
//* --------------------------------------------------------------------
void resid(Bump& bump) {
  double visc = 1.0 / bump.reynld;

  for (size_t it = 0; it < bump.neqn; it++) {
    bump.res[it] = 0.0;
  }

  for (size_t it = 0; it < bump.nelemn; it++) {
    double ar = bump.area[it] / 3.0;

    for (size_t iquad = 0; iquad < NQUAD; iquad++) {
      double yq = bump.ym[it * NQUAD + iquad];
      double xq = bump.xm[it * NQUAD + iquad];

      double det = 0.0, etax = 0.0, etay = 0.0, xix = 0.0, xiy = 0.0;

      if (bump.isotri[it] == 1) {
        std::tie(det, etax, etay, xix, xiy) = trans(it, bump.node, bump.xc, xq, bump.yc, yq);
        ar = det * bump.area[it] / 3.0;
      }

      auto [un, unx, uny] = uval(etax, etay, bump, it, xix, xiy, xq, yq);

      for (size_t iq = 0; iq < NNODES; iq++) {
        size_t ip = bump.node[it * NNODES + iq];
        double bb = bump.phi[phi_idx(it, iquad, iq, 0)];
        double bx = bump.phi[phi_idx(it, iquad, iq, 1)];
        double by = bump.phi[phi_idx(it, iquad, iq, 2)];
        double bbl = bump.psi[psi_idx(it, iquad, iq)];
        int ihor = bump.indx[ip * 2];
        int iver = bump.indx[ip * 2 + 1];
        int iprs = bump.insc[ip];

        if (0 < ihor) {
          bump.res[(size_t)(ihor - 1)] -= ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
        }
        if (0 < iver) {
          bump.res[(size_t)(iver - 1)] -= ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
        }

        for (size_t iqq = 0; iqq < NNODES; iqq++) {
          size_t ipp = bump.node[it * NNODES + iqq];
          double bbb = bump.phi[phi_idx(it, iquad, iqq, 0)];
          double bbx = bump.phi[phi_idx(it, iquad, iqq, 1)];
          double bby = bump.phi[phi_idx(it, iquad, iqq, 2)];
          double bbbl = bump.psi[psi_idx(it, iquad, iqq)];
          int ju = bump.indx[ipp * 2];
          int jv = bump.indx[ipp * 2 + 1];
          int jp = bump.insc[ipp];

          if (0 < ju) {
            if (0 < ihor) {
              bump.res[(size_t)(ihor - 1)] += ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]))
                * bump.g[(size_t)(ju - 1)];
            }
            if (0 < iver) {
              bump.res[(size_t)(iver - 1)] += ar * bb * bbb * unx[1] * bump.g[(size_t)(ju - 1)];
            }
            if (0 < iprs) {
              bump.res[(size_t)(iprs - 1)] += ar * bbx * bbl * bump.g[(size_t)(ju - 1)];
            }
          } else if (ju == -1) {
            double ubc = ubdry(1, bump.yc[ipp]);
            if (0 < ihor) {
              double aij = ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
              bump.res[(size_t)(ihor - 1)] += ubc * aij;
            }
            if (0 < iver) {
              double aij = ar * bb * bbb * unx[1];
              bump.res[(size_t)(iver - 1)] += ubc * aij;
            }
            if (0 < iprs) {
              double aij = ar * bbx * bbl;
              bump.res[(size_t)(iprs - 1)] += ubc * aij;
            }
          } else if (ju == -2) {
            double ubc = ubump(bump, ipp, iqq, it, 1);
            if (0 < ihor) {
              double aij = ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
              bump.res[(size_t)(ihor - 1)] += ubc * aij;
            }
            if (0 < iver) {
              double aij = ar * bb * bbb * unx[1];
              bump.res[(size_t)(iver - 1)] += ubc * aij;
            }
            if (0 < iprs) {
              double aij = ar * bbx * bbl;
              bump.res[(size_t)(iprs - 1)] += ubc * aij;
            }
          }

          if (0 < jv) {
            if (0 < ihor) {
              bump.res[(size_t)(ihor - 1)] += ar * bb * bbb * uny[0] * bump.g[(size_t)(jv - 1)];
            }
            if (0 < iver) {
              bump.res[(size_t)(iver - 1)] += ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]))
                * bump.g[(size_t)(jv - 1)];
            }
            if (0 < iprs) {
              bump.res[(size_t)(iprs - 1)] += ar * bby * bbl * bump.g[(size_t)(jv - 1)];
            }
          } else if (jv == -1) {
            double ubc = ubdry(2, bump.yc[ipp]);
            if (0 < ihor) {
              bump.res[(size_t)(ihor - 1)] += ubc * ar * bb * bbb * uny[0];
            }
            if (0 < iver) {
              double aij = ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
              bump.res[(size_t)(iver - 1)] += ubc * aij;
            }
            if (0 < iprs) {
              bump.res[(size_t)(iprs - 1)] += ubc * ar * bby * bbl;
            }
          } else if (jv == -2) {
            double ubc = ubump(bump, ipp, iqq, it, 2);
            if (0 < ihor) {
              bump.res[(size_t)(ihor - 1)] += ubc * ar * bb * bbb * uny[0];
            }
            if (0 < iver) {
              double aij = ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
              bump.res[(size_t)(iver - 1)] += ubc * aij;
            }
            if (0 < iprs) {
              bump.res[(size_t)(iprs - 1)] += ubc * ar * bby * bbl;
            }
          }

          if (0 < jp) {
            if (0 < ihor) {
              bump.res[(size_t)(ihor - 1)] -= ar * bx * bbbl * bump.g[(size_t)(jp - 1)];
            }
            if (0 < iver) {
              bump.res[(size_t)(iver - 1)] -= ar * by * bbbl * bump.g[(size_t)(jp - 1)];
            }
          }
        }
      }
    }
  }

  bump.res[bump.neqn - 1] = bump.g[bump.neqn - 1];

  double rmax = 0.0;
  size_t imax = 0;
  size_t ibad = 0;

  for (size_t i = 0; i < bump.neqn; i++) {
    double test = std::abs(bump.res[i]);
    if (rmax < test) {
      rmax = test;
      imax = i;
    }
    if (1.0e-3 < test) {
      ibad += 1;
    }
  }

  if (1 <= bump.iwrite) {
    std::println();
    std::println("RESIDUAL INFORMATION:");
    std::println();
    std::println("Worst residual is number {}", imax + 1);
    std::println("of magnitude {}", rmax);
    std::println();
    std::println("Number of \"bad\" residuals is {} out of {}", ibad, bump.neqn);
    std::println();
  }

  if (2 <= bump.iwrite) {
    std::println("Raw residuals:");
    std::println();
    size_t idx = 0;
    for (size_t j = 0; j < bump.n_points; j++) {
      if (0 < bump.indx[j * 2]) {
        if (std::abs(bump.res[idx]) <= 1.0e-3) {
          std::println(" U {} {} {}", idx + 1, j + 1, bump.res[idx]);
        } else {
          std::println("*U {} {} {}", idx + 1, j + 1, bump.res[idx]);
        }
        idx += 1;
      }
      if (0 < bump.indx[j * 2 + 1]) {
        if (std::abs(bump.res[idx]) <= 1.0e-3) {
          std::println(" V {} {} {}", idx + 1, j + 1, bump.res[idx]);
        } else {
          std::println("*V {} {} {}", idx + 1, j + 1, bump.res[idx]);
        }
        idx += 1;
      }
      if (0 < bump.insc[j]) {
        if (std::abs(bump.res[idx]) <= 1.0e-3) {
          std::println(" P {} {} {}", idx + 1, j + 1, bump.res[idx]);
        } else {
          std::println("*P {} {} {}", idx + 1, j + 1, bump.res[idx]);
        }
        idx += 1;
      }
    }
  }
}

//* --------------------------------------------------------------------
//*  SETBAN - compute the half band width
//* --------------------------------------------------------------------
void setban(Bump& bump) {
  bump.nlband = 0;

  for (size_t it = 0; it < bump.nelemn; it++) {
    for (size_t iq = 0; iq < NNODES; iq++) {
      size_t ip = bump.node[it * NNODES + iq];
      for (size_t iuk = 0; iuk < 3; iuk++) {
        int i_val = (iuk == 2) ? bump.insc[ip] : bump.indx[ip * 2 + iuk];
        if (0 < i_val) {
          for (size_t iqq = 0; iqq < NNODES; iqq++) {
            size_t ipp = bump.node[it * NNODES + iqq];
            for (size_t iukk = 0; iukk < 3; iukk++) {
              int j_val = (iukk == 2) ? bump.insc[ipp] : bump.indx[ipp * 2 + iukk];
              if (0 < j_val && j_val > i_val) {
                size_t diff = (size_t)(j_val - i_val);
                if (bump.nlband < diff) {
                  bump.nlband = diff;
                }
              }
            }
          }
        }
      }
    }
  }

  bump.nband = bump.nlband + bump.nlband + 1;
  bump.nrow = bump.nlband + bump.nlband + bump.nlband + 1;

  std::println();
  std::println("SETBAN:");
  std::println();
  std::println("  Lower bandwidth = {}", bump.nlband);
  std::println("  Total bandwidth = {}", bump.nband);
  std::println("  Required matrix rows = {}", bump.nrow);

  if (bump.maxrow < bump.nrow) {
    std::println("SETBAN - NROW is too large!");
    std::println("The maximum allowed is {}", bump.maxrow);
    std::println("This problem requires NROW = {}", bump.nrow);
    std::exit(1);
  }
}

//* --------------------------------------------------------------------
//*  SETBAS - evaluate basis functions at each integration point
//* --------------------------------------------------------------------
void setbas(Bump& bump) {
  for (size_t it = 0; it < bump.nelemn; it++) {
    for (size_t j = 0; j < NQUAD; j++) {
      double xq = bump.xm[it * NQUAD + j];
      double yq = bump.ym[it * NQUAD + j];
      auto [det, etax, etay, xix, xiy] = trans(it, bump.node, bump.xc, xq, bump.yc, yq);

      for (size_t iq = 0; iq < NNODES; iq++) {
        double bb, bx, by;
        if (bump.isotri[it] == 0) {
          bump.psi[psi_idx(it, j, iq)] = bsp(it, iq, 0, bump.node, bump.xc, xq, bump.yc, yq);
          std::tie(bb, bx, by) = qbf(xq, yq, it, iq, bump.node, bump.xc, bump.yc);
        } else {
          std::tie(bb, bx, by) = refqbf(xq, yq, iq, etax, etay, xix, xiy);
          bump.psi[psi_idx(it, j, iq)] = refbsp(xq, yq, iq);
        }
        bump.phi[phi_idx(it, j, iq, 0)] = bb;
        bump.phi[phi_idx(it, j, iq, 1)] = bx;
        bump.phi[phi_idx(it, j, iq, 2)] = by;
      }
    }
  }
}

//* --------------------------------------------------------------------
//*  SETGRD - set up the geometric grid
//* --------------------------------------------------------------------
void setgrd(int ibump, Bump& bump) {
  std::println();
  std::println("SETGRD:");
  std::println();

  if (bump.ny < bump.nx) {
    bump._long = true;
    std::println("Using vertical ordering.");
  } else {
    bump._long = false;
    std::println("Using horizontal ordering.");
  }

  if (ibump == 0) {
    std::println("No isoparametric elements will be used.");
  } else if (ibump == 1) {
    std::println("Isoparametric elements directly on bump.");
  } else if (ibump == 2) {
    std::println("All elements above bump are isoparametric.");
  } else if (ibump == 3) {
    std::println("All elements are isoparametric.");
  } else {
    std::println("Unexpected value of IBUMP = {}", ibump);
    std::exit(1);
  }

  size_t nbleft = (size_t)std::round(bump.xbleft * (double)(bump.mx - 1) / bump.xlngth);
  size_t nbrite = (size_t)std::round(bump.xbrite * (double)(bump.mx - 1) / bump.xlngth);
  std::println("Bump extends from {} at node {}", bump.xbleft, nbleft + 1);
  std::println("               to {} at node {}", bump.xbrite, nbrite + 1);

  bump.neqn = 0;
  size_t ielemn = 0;

  for (size_t ip = 0; ip < bump.n_points; ip++) {
    size_t ic, jc;
    if (bump._long) {
      ic = ip / bump.my;
      jc = ip % bump.my;
    } else {
      ic = ip % bump.mx;
      jc = ip / bump.mx;
    }
    size_t icnt = (ic + 1) % 2;
    size_t jcnt = (jc + 1) % 2;

    if ((icnt == 1 && jcnt == 1) && (ic != bump.mx - 1) && (jc != bump.my - 1)) {
      size_t ip1, ip2;
      if (bump._long) {
        ip1 = ip + bump.my;
        ip2 = ip + bump.my + bump.my;

        bump.node[ielemn * NNODES] = ip;
        bump.node[ielemn * NNODES + 1] = ip + 2;
        bump.node[ielemn * NNODES + 2] = ip2 + 2;
        bump.node[ielemn * NNODES + 3] = ip + 1;
        bump.node[ielemn * NNODES + 4] = ip1 + 2;
        bump.node[ielemn * NNODES + 5] = ip1 + 1;

        if (ibump == 0) {
          bump.isotri[ielemn] = 0;
        } else if (ibump == 1) {
          bump.isotri[ielemn] = 0;
        } else if (ibump == 2) {
          bump.isotri[ielemn] = (nbleft <= ic && ic < nbrite) ? 1 : 0;
        } else {
          bump.isotri[ielemn] = 1;
        }

        ielemn += 1;

        bump.node[ielemn * NNODES] = ip;
        bump.node[ielemn * NNODES + 1] = ip2 + 2;
        bump.node[ielemn * NNODES + 2] = ip2;
        bump.node[ielemn * NNODES + 3] = ip1 + 1;
        bump.node[ielemn * NNODES + 4] = ip2 + 1;
        bump.node[ielemn * NNODES + 5] = ip1;

        if (ibump == 0) {
          bump.isotri[ielemn] = 0;
        } else if (ibump == 1) {
          bump.isotri[ielemn] = (jc == 0 && nbleft <= ic && ic < nbrite) ? 1 : 0;
        } else if (ibump == 2) {
          bump.isotri[ielemn] = (nbleft <= ic && ic < nbrite) ? 1 : 0;
        } else {
          bump.isotri[ielemn] = 1;
        }

        ielemn += 1;
      } else {
        ip1 = ip + bump.mx;
        ip2 = ip + bump.mx + bump.mx;

        bump.node[ielemn * NNODES] = ip;
        bump.node[ielemn * NNODES + 1] = ip2;
        bump.node[ielemn * NNODES + 2] = ip2 + 2;
        bump.node[ielemn * NNODES + 3] = ip1;
        bump.node[ielemn * NNODES + 4] = ip2 + 1;
        bump.node[ielemn * NNODES + 5] = ip1 + 1;

        if (ibump == 0) {
          bump.isotri[ielemn] = 0;
        } else if (ibump == 1) {
          bump.isotri[ielemn] = 0;
        } else if (ibump == 2) {
          bump.isotri[ielemn] = (nbleft <= ic && ic < nbrite) ? 1 : 0;
        } else {
          bump.isotri[ielemn] = 1;
        }

        ielemn += 1;

        bump.node[ielemn * NNODES] = ip;
        bump.node[ielemn * NNODES + 1] = ip2 + 2;
        bump.node[ielemn * NNODES + 2] = ip + 2;
        bump.node[ielemn * NNODES + 3] = ip1 + 1;
        bump.node[ielemn * NNODES + 4] = ip1 + 2;
        bump.node[ielemn * NNODES + 5] = ip + 1;

        if (ibump == 0) {
          bump.isotri[ielemn] = 0;
        } else if (ibump == 1) {
          bump.isotri[ielemn] = (jc == 0 && nbleft <= ic && ic < nbrite) ? 1 : 0;
        } else if (ibump == 2) {
          bump.isotri[ielemn] = (nbleft <= ic && ic < nbrite) ? 1 : 0;
        } else {
          bump.isotri[ielemn] = 1;
        }

        ielemn += 1;
      }
    }

    if (ic == 0 && 0 < jc && jc < bump.my - 1) {
      bump.indx[ip * 2] = -1;
      bump.indx[ip * 2 + 1] = -1;
    } else if (ic == bump.mx - 1 && 0 < jc && jc < bump.my - 1) {
      bump.neqn += 1;
      bump.indx[ip * 2] = (int)bump.neqn;
      bump.indx[ip * 2 + 1] = 0;
    } else if (jc == 0 && ielemn > 0 && bump.isotri[ielemn - 1] == 1) {
      bump.indx[ip * 2] = -2;
      bump.indx[ip * 2 + 1] = -2;
    } else if (ic == 0 || ic == bump.mx - 1 || jc == 0 || jc == bump.my - 1) {
      bump.indx[ip * 2] = 0;
      bump.indx[ip * 2 + 1] = 0;
    } else {
      bump.neqn += 2;
      bump.indx[ip * 2] = (int)(bump.neqn - 1);
      bump.indx[ip * 2 + 1] = (int)bump.neqn;
    }

    if (jcnt == 1 && icnt == 1) {
      bump.neqn += 1;
      bump.insc[ip] = (int)bump.neqn;
    } else {
      bump.insc[ip] = 0;
    }
  }

  if (1 <= bump.iwrite) {
    std::println();
    std::println("     I     indx 1, indx 2, insc");
    std::println();
    for (size_t i = 0; i < bump.n_points; i++) {
      std::println("{}\t{}\t{}\t{}", i + 1, bump.indx[i * 2], bump.indx[i * 2 + 1], bump.insc[i]);
    }
    std::println();
    std::println("Isoparametric triangles:");
    std::println();
    for (size_t i = 0; i < bump.nelemn; i++) {
      if (bump.isotri[i] == 1) {
        std::println("{}", i + 1);
      }
    }
    std::println();
    std::println("   IT   node(IT,*)");
    std::println();
    for (size_t it = 0; it < bump.nelemn; it++) {
      std::print("{}\t", it + 1);
      for (size_t i = 0; i < 6; i++) {
        std::print("{}\t", bump.node[it * NNODES + i] + 1);
      }
      std::println();
    }
  }

  std::println();
  std::println("SETGRD: Number of unknowns = {}", bump.neqn);

  if (bump.maxeqn < bump.neqn) {
    std::println("SETGRD - Too many unknowns!");
    std::println("The maximum allowed is MAXEQN = {}", bump.maxeqn);
    std::println("This problem requires neqn = {}", bump.neqn);
    std::exit(1);
  }
}

//* --------------------------------------------------------------------
//*  SETLIN - determine unknown numbers along the profile line
//* --------------------------------------------------------------------
void setlin(Bump& bump) {
  size_t itemp = (size_t)std::round(
    2.0 * (double)(bump.nx - 1) * bump.xprof / bump.xlngth);

  size_t nodex0 = bump._long
    ? itemp * (2 * bump.ny - 1)
    : itemp;

  std::println();
  std::println("SETLIN:");
  std::println();
  std::println("  Profile generated at X = {}", bump.xprof);
  std::println("  which is above node  = {}", nodex0 + 1);

  for (size_t i = 0; i < bump.my; i++) {
    size_t ip = bump._long
      ? nodex0 + i
      : nodex0 + bump.mx * i;
    bump.iline[i] = bump.indx[ip * 2];
  }

  bump.iline_inv.assign(bump.neqn + 1, -1);
  for (size_t j = 0; j < bump.my; j++) {
    int val = bump.iline[j];
    if (val > 0) {
      bump.iline_inv[(size_t)val] = (int)(j + 1);
    }
  }

  if (1 <= bump.iwrite) {
    std::println();
    std::println("  Indices of unknowns along the profile line:");
    std::println();
    for (size_t i = 0; i < bump.my; i += 5) {
      for (size_t j = i; j < std::min(i + 5, bump.my); j++) {
        std::print("{:5}", bump.iline[j]);
      }
      std::println();
    }
  }
}

//* --------------------------------------------------------------------
//*  SETQUD - set midpoint quadrature rule information
//* --------------------------------------------------------------------
void setqud(Bump& bump) {
  for (size_t it = 0; it < bump.nelemn; it++) {
    size_t ip1 = bump.node[it * NNODES];
    size_t ip2 = bump.node[it * NNODES + 1];
    size_t ip3 = bump.node[it * NNODES + 2];
    double x1 = bump.xc[ip1];
    double x2 = bump.xc[ip2];
    double x3 = bump.xc[ip3];
    double y1 = bump.yc[ip1];
    double y2 = bump.yc[ip2];
    double y3 = bump.yc[ip3];

    if (bump.isotri[it] == 0) {
      bump.xm[it * NQUAD]     = 0.5 * (x1 + x2);
      bump.xm[it * NQUAD + 1] = 0.5 * (x2 + x3);
      bump.xm[it * NQUAD + 2] = 0.5 * (x3 + x1);
      bump.ym[it * NQUAD]     = 0.5 * (y1 + y2);
      bump.ym[it * NQUAD + 1] = 0.5 * (y2 + y3);
      bump.ym[it * NQUAD + 2] = 0.5 * (y3 + y1);
      bump.area[it] = 0.5 * std::abs(
        (y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3));
    } else {
      bump.xm[it * NQUAD]     = 0.5;
      bump.ym[it * NQUAD]     = 0.5;
      bump.xm[it * NQUAD + 1] = 1.0;
      bump.ym[it * NQUAD + 1] = 0.5;
      bump.xm[it * NQUAD + 2] = 0.5;
      bump.ym[it * NQUAD + 2] = 0.0;
      bump.area[it] = 0.5;
    }
  }

  if (3 <= bump.iwrite) {
    std::println();
    std::println("SETQUD: Element Areas and Quadrature points:");
    std::println();
    for (size_t i = 0; i < bump.nelemn; i++) {
      std::println("{}\t{}", i + 1, bump.area[i]);
      for (size_t j = 0; j < NQUAD; j++) {
        std::println("{}\t{}\t{}\t{}", i + 1, j + 1, bump.xm[i * NQUAD + j], bump.ym[i * NQUAD + j]);
      }
    }
  }
}

//* --------------------------------------------------------------------
//*  SETXY - set the grid coordinates based on the parameter value
//* --------------------------------------------------------------------
void setxy(Bump& bump) {
  for (size_t ip = 0; ip < bump.n_points; ip++) {
    size_t ic, jc;
    if (bump._long) {
      ic = ip / bump.my;
      jc = ip % bump.my;
    } else {
      ic = ip % bump.mx;
      jc = ip / bump.mx;
    }
    bump.xc[ip] = (double)ic * bump.xlngth / (double)(bump.nx + bump.nx - 2);

    double ybot = -bump.ypert * (bump.xc[ip] - 3.0) * (bump.xc[ip] - 1.0);
    double ylo = std::max(0.0, ybot);

    bump.yc[ip] = ((double)(bump.my - 1 - jc) * ylo + (double)jc * bump.ylngth)
                  / (double)(bump.ny + bump.ny - 2);
  }

  if (2 <= bump.iwrite) {
    std::println();
    std::println("SETXY:");
    std::println();
    std::println("     I     XC     YC");
    std::println();
    for (size_t i = 0; i < bump.n_points; i++) {
      std::println("{}\t{}\t{}", i + 1, bump.xc[i], bump.yc[i]);
    }
  }
}

//* --------------------------------------------------------------------
//*  TRANS - calculate the element transformation mapping
//* --------------------------------------------------------------------
std::tuple<double, double, double, double, double> trans(
  size_t it, const std::vector<size_t> &node,
  const std::vector<double> &xc, double xq,
  const std::vector<double> &yc, double yq
) {
  size_t i1 = node[it * NNODES];
  size_t i2 = node[it * NNODES + 1];
  size_t i3 = node[it * NNODES + 2];
  size_t i4 = node[it * NNODES + 3];
  size_t i5 = node[it * NNODES + 4];
  size_t i6 = node[it * NNODES + 5];

  double x1 = xc[i1]; double y1 = yc[i1];
  double x2 = xc[i2]; double y2 = yc[i2];
  double x3 = xc[i3]; double y3 = yc[i3];
  double x4 = xc[i4]; double y4 = yc[i4];
  double x5 = xc[i5]; double y5 = yc[i5];
  double x6 = xc[i6]; double y6 = yc[i6];

  double a1 = 2.0 * x3 - 4.0 * x6 + 2.0 * x1;
  double b1 = -4.0 * x3 - 4.0 * x4 + 4.0 * x5 + 4.0 * x6;
  double c1 = 2.0 * x2 + 2.0 * x3 - 4.0 * x5;
  double d1 = -3.0 * x1 - x3 + 4.0 * x6;
  double e1 = -x2 + x3 + 4.0 * x4 - 4.0 * x6;

  double a2 = 2.0 * y3 - 4.0 * y6 + 2.0 * y1;
  double b2 = -4.0 * y3 - 4.0 * y4 + 4.0 * y5 + 4.0 * y6;
  double c2 = 2.0 * y2 + 2.0 * y3 - 4.0 * y5;
  double d2 = -3.0 * y1 - y3 + 4.0 * y6;
  double e2 = -y2 + y3 + 4.0 * y4 - 4.0 * y6;

  double dxdxi = 2.0 * a1 * xq + b1 * yq + d1;
  double dxdeta = b1 * xq + 2.0 * c1 * yq + e1;
  double dydxi = 2.0 * a2 * xq + b2 * yq + d2;
  double dydeta = b2 * xq + 2.0 * c2 * yq + e2;

  double det = (2.0 * a1 * b2 - 2.0 * a2 * b1) * xq * xq
    + (4.0 * a1 * c2 - 4.0 * a2 * c1) * xq * yq
    + (2.0 * b1 * c2 - 2.0 * b2 * c1) * yq * yq
    + (2.0 * a1 * e2 + b2 * d1 - b1 * d2 - 2.0 * a2 * e1) * xq
    + (2.0 * c2 * d1 + b1 * e2 - b2 * e1 - 2.0 * c1 * d2) * yq
    + d1 * e2 - d2 * e1;

  double eps = 1e-30;
  if (std::abs(det) < eps) {
    det = eps;
  }

  double xix   =  dydeta / det;
  double xiy   = -dxdeta / det;
  double etax  = -dydxi  / det;
  double etay  =  dxdxi  / det;

  return {det, etax, etay, xix, xiy};
}

//* --------------------------------------------------------------------
//*  UBDRY - parabolic inflow boundary condition
//* --------------------------------------------------------------------
double ubdry(size_t iuk, double yy) {
  if (iuk == 1) {
    return (-2.0 * yy + 6.0) * yy / 9.0;
  } else {
    return 0.0;
  }
}

//* --------------------------------------------------------------------
//*  UBUMP - sensitivity dU/dA on the bump
//* --------------------------------------------------------------------
double ubump(Bump& bump, size_t ip, size_t iqq, size_t it, size_t iukk) {
  double det = 1.0;
  double etax = 0.0;
  double etay = 0.0;
  double xix = 1.0;
  double xiy = 1.0;
  double xq, yq;

  if (bump.isotri[it] == 0) {
    xq = bump.xc[ip];
    yq = bump.yc[ip];
  } else {
    switch (iqq) {
      case 0: xq = 0.0; yq = 0.0; break;
      case 1: xq = 1.0; yq = 1.0; break;
      case 2: xq = 1.0; yq = 0.0; break;
      case 3: xq = 0.5; yq = 0.5; break;
      case 4: xq = 1.0; yq = 0.5; break;
      case 5: xq = 0.5; yq = 0.0; break;
      default: xq = 0.0; yq = 0.0;
    }
    std::tie(det, etax, etay, xix, xiy) = trans(it, bump.node, bump.xc, xq, bump.yc, yq);
  }

  auto [_un, _unx, uny] = _ubump_uval(bump, it, xix, xiy, xq, yq, det, etax, etay);

  if (iukk == 1) {
    return -uny[0] * (bump.xc[ip] - 1.0) * (bump.xc[ip] - 3.0);
  } else if (iukk == 2) {
    return -uny[1] * (bump.xc[ip] - 1.0) * (bump.xc[ip] - 3.0);
  } else {
    std::println("UBUMP called for iukk = {}", iukk);
    std::exit(1);
  }
}

std::tuple<std::array<double, 2>, std::array<double, 2>, std::array<double, 2>> _ubump_uval(
  Bump& bump, size_t it, double xix, double xiy, double xq, double yq,
  double _det, double etax, double etay
) {
  std::array<double, 2> un = {0.0, 0.0};
  std::array<double, 2> unx = {0.0, 0.0};
  std::array<double, 2> uny = {0.0, 0.0};

  for (size_t iq = 0; iq < NNODES; iq++) {
    auto [bb, bx, by] = (bump.isotri[it] == 1)
      ? refqbf(xq, yq, iq, etax, etay, xix, xiy)
      : qbf(xq, yq, it, iq, bump.node, bump.xc, bump.yc);
    size_t ip_local = bump.node[it * NNODES + iq];

    for (size_t iuk = 0; iuk < 2; iuk++) {
      int iun = bump.indx[ip_local * 2 + iuk];
      if (0 < iun) {
        un[iuk]  += bb  * bump.g[(size_t)(iun - 1)];
        unx[iuk] += bx  * bump.g[(size_t)(iun - 1)];
        uny[iuk] += by  * bump.g[(size_t)(iun - 1)];
      } else if (iun == -1) {
        double ubc = ubdry(iuk + 1, bump.yc[ip_local]);
        un[iuk]  += bb  * ubc;
        unx[iuk] += bx  * ubc;
        uny[iuk] += by  * ubc;
      }
    }
  }
  return {un, unx, uny};
}

//* --------------------------------------------------------------------
//*  UVAL - evaluate velocities at a given quadrature point
//* --------------------------------------------------------------------
std::tuple<std::array<double, 2>, std::array<double, 2>, std::array<double, 2>> uval(
  double etax, double etay, Bump& bump, size_t it,
  double xix, double xiy, double xq, double yq
) {
  std::array<double, 2> un = {0.0, 0.0};
  std::array<double, 2> unx = {0.0, 0.0};
  std::array<double, 2> uny = {0.0, 0.0};

  for (size_t iq = 0; iq < NNODES; iq++) {
    auto [bb, bx, by] = (bump.isotri[it] == 1)
      ? refqbf(xq, yq, iq, etax, etay, xix, xiy)
      : qbf(xq, yq, it, iq, bump.node, bump.xc, bump.yc);
    size_t ip_local = bump.node[it * NNODES + iq];

    for (size_t iuk = 0; iuk < 2; iuk++) {
      int iun = bump.indx[ip_local * 2 + iuk];
      if (0 < iun) {
        un[iuk]  += bb  * bump.g[(size_t)(iun - 1)];
        unx[iuk] += bx  * bump.g[(size_t)(iun - 1)];
        uny[iuk] += by  * bump.g[(size_t)(iun - 1)];
      } else if (iun == -1) {
        double ubc = ubdry(iuk + 1, bump.yc[ip_local]);
        un[iuk]  += bb  * ubc;
        unx[iuk] += bx  * ubc;
        uny[iuk] += by  * ubc;
      }
    }
  }
  return {un, unx, uny};
}

//* --------------------------------------------------------------------
//*  UV_WRITE - write a velocity file
//* --------------------------------------------------------------------
void uv_write(Bump& bump, const std::vector<double> &f, std::ofstream& file) {
  for (size_t ip = 0; ip < bump.n_points; ip++) {
    int k_u = bump.indx[ip * 2];
    double u;
    if (k_u < 0) {
      u = ubdry(1, bump.yc[ip]);
    } else if (k_u == 0) {
      u = 0.0;
    } else {
      u = f[(size_t)(k_u - 1)];
    }

    int k_v = bump.indx[ip * 2 + 1];
    double v;
    if (k_v < 0) {
      v = ubdry(2, bump.yc[ip]);
    } else if (k_v == 0) {
      v = 0.0;
    } else {
      v = f[(size_t)(k_v - 1)];
    }

    file << "  " << std::scientific << std::setw(14) << u
         << "  " << std::scientific << std::setw(14) << v << "\n";
  }
}

//* --------------------------------------------------------------------
//*  XY_WRITE - write node coordinate data
//* --------------------------------------------------------------------
void xy_write(Bump& bump, std::ofstream& file) {
  for (size_t ip = 0; ip < bump.n_points; ip++) {
    file << "  " << std::scientific << std::setw(14) << bump.xc[ip]
         << "  " << std::scientific << std::setw(14) << bump.yc[ip] << "\n";
  }
}
