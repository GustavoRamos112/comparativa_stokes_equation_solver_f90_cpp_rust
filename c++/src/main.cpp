//* --------------------------------------------------------------------
//*  Main program
//* --------------------------------------------------------------------

#include <vector>
#include <string>
#include <tuple>

#include <chrono>

#include <format>
#include <print>

#include <fstream>
#include <iomanip>

#include <cstdlib>
#include <cmath>

#include <array>

#include <filesystem>
namespace fs = std::filesystem;

const bool save_times = true;

//* -------------------------------------------------------------------
//* Declaracion de funciones
void timestamp();
std::tuple<bool, int> setgrd(
  int ibump, 
  std::vector<std::vector<int>> &indx, 
  std::vector<int> &insc,
  std::vector<int> &isotri, 
  int iwrite, bool _long, int maxeqn,
  int mx, int my, int nelemn,
  int neqn, int nnodes,
  std::vector<std::vector<int>> &node, 
  int n_points, int nx, int ny,
  double xbleft, double xbrite, double xlngth
);
void setxy(
  int iwrite, bool _long, int mx, int my, int n_points, 
  int nx, int ny, 
  std::vector<double> &xc, double xlngth, 
  std::vector<double> &yc, double ylngth, double ypert
);
void setqud(
  std::vector<double> &area, std::vector<int> &isotri, int iwrite,
  int nelemn, int nnodes, std::vector<std::vector<int>> &node,
  int n_points, int nquad,
  std::vector<double> &xc, std::vector<std::vector<double>> &xm,
  std::vector<double> &yc, std::vector<std::vector<double>> &ym
);
void setbas(
  const std::vector<int> &isotri, int nelemn, int nnodes, 
  std::vector<std::vector<int>> &node, int n_points, int nquad, 
  std::vector<std::vector<std::vector<std::vector<double>>>> &phi, 
  std::vector<std::vector<std::vector<double>>> &psi, 
  std::vector<double> &xc, 
  std::vector<std::vector<double>> &xm, 
  std::vector<double> &yc, 
  std::vector<std::vector<double>> &ym
);
std::tuple<double, double, double, double, double> trans(
  int it, int nelemn, int nnodes, 
  const std::vector<std::vector<int>> &node, 
  int n_points, 
  const std::vector<double> &xc, double xq, 
  const std::vector<double> &yc, double yq
);
double bsp(
  int it, int iq, int id, int nelemn, int nnodes, 
  const std::vector<std::vector<int>> &node, int n_points, 
  const std::vector<double> &xc, double xq, 
  const std::vector<double> &yc, double yq
);
inline int i4_modp(int i, int j);
inline int i4_wrap(int ival, int ilo, int ihi);
inline double ubdry(int iuk, double yy);
inline double refbsp(double xq, double yq, int iq);
inline int igetl(int i, std::vector<int> &iline, int my);
std::tuple<double, double, double> refqbf(
  double x, double y, int inn, 
  double etax, double etay, 
  double xix, double xiy
);
std::tuple<double, double, double> qbf(
  double xq, double yq, int it, 
  int inn, int nelemn, int nnodes, 
  const std::vector<std::vector<int>> &node, int n_points, 
  const std::vector<double> &xc, const std::vector<double> &yc
);
void setlin(
  std::vector<int> &iline, 
  std::vector<std::vector<int>> &indx, 
  int iwrite, bool _long, int mx, int my, 
  int n_points, int nx, int ny, 
  double xlngth, double xprof
);
std::tuple<int, int, int> setban(
  std::vector<std::vector<int>> &indx,
  std::vector<int> &insc, int maxrow,
  int nband, int nelemn, int nlband,
  int nnodes, std::vector<std::vector<int>> &node, 
  int n_points, int nrow
);
int nstoke(
  std::vector<std::vector<double>> &a,
  std::vector<double> &area,
  std::vector<double> &f, std::vector<double> &g,
  std::vector<std::vector<int>> &indx,
  std::vector<int> &insc, std::vector<int> &isotri,
  int maxnew, int maxrow, int nband,
  int nelemn, int neqn, int nlband,
  int nnodes, std::vector<std::vector<int>> &node,
  int n_points, int nquad, int nrow,
  int numnew,
  std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
  std::vector<std::vector<std::vector<double>>> &psi,
  double reynld, double tolnew,
  std::vector<double> &xc, std::vector<std::vector<double>> &xm,
  std::vector<double> &yc, std::vector<std::vector<double>> &ym
);
int idamax_v(int n, std::vector<double> &dx, int incx);
int idamax_m(int n, const std::vector<std::vector<double>>& abd, int col, int start_row, int incx);
void dscal_v(
  int n, double sa, 
  std::vector<double>& x, int incx
);
void dscal_m(
  int n, double sa, 
  std::vector<std::vector<double>>& abd, 
  int col, int start_row, int incx
);
void daxpy_v(
  int n, double da, 
  std::vector<double>& dx, int incx, 
  std::vector<double>& dy, int incy
);
void daxpy_m(
  int n, double da, 
  std::vector<std::vector<double>> &abd_x, 
  int col_x, int start_row_x, int incx,
  std::vector<std::vector<double>>& abd_y, 
  int col_y, int start_row_y, int incy
);
void daxpy(
  int n, double da, 
  const std::vector<std::vector<double>>& abd_x, int col_x, int start_row_x, int incx,
  std::vector<double>& b_y, int start_idx_y, int incy
);
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> uval(
  double etax,
  double etay,
  const std::vector<double> &g,
  const std::vector<std::vector<int>> &indx,
  const std::vector<int> &isotri,
  int it,
  int nelemn,
  int neqn,
  int nnodes,
  const std::vector<std::vector<int>> &node,
  int n_points,
  const std::vector<double> &xc,
  double xix,
  double xiy,
  double xq,
  const std::vector<double> &yc,
  double yq
);
double ubump(
  const std::vector<double> &g,
  const std::vector<std::vector<int>> &indx,
  int ip, int iqq,
  const std::vector<int> &isotri,
  int it, int iukk, int nelemn,
  int neqn, int nnodes,
  const std::vector<std::vector<int>> &node, int n_points,
  const std::vector<double> &xc, const std::vector<double> &yc
);
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> _ubump_uval(
  const std::vector<double> &g,
  const std::vector<std::vector<int>> &indx,
  const std::vector<int> &isotri,
  int it, int nelemn, int neqn, int nnodes,
  const std::vector<std::vector<int>> &node,
  int n_points, const std::vector<double> &xc,
  double xix, double xiy, double xq,
  const std::vector<double> &yc,
  double yq, double det, double etax, double etay
);
double ddot(
  int n, 
  const std::vector<std::vector<double>> &abd, 
  int col, int row_start, int incx, 
  const std::vector<double> &b, int start_b, int incy
);
void resid(
  std::vector<double> &area, std::vector<double> &g,
  std::vector<std::vector<int>> &indx,
  std::vector<int> &insc, std::vector<int> &isotri,
  int iwrite, int nelemn, int neqn,
  int nnodes, std::vector<std::vector<int>> &node,
  int n_points, int nquad,
  std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
  std::vector<std::vector<std::vector<double>>> &psi, 
  std::vector<double> &res, double reynld, 
  std::vector<double> &xc, std::vector<std::vector<double>> &xm, 
  std::vector<double> &yc, std::vector<std::vector<double>> &ym
);
std::vector<double> getg(
  const std::vector<double> &f, const std::vector<int> &iline, int my, int neqn
);
void gram(
  std::vector<std::vector<double>> &gr, 
  std::vector<int> &iline, 
  std::vector<std::vector<int>> &indx, 
  int iwrite, int my, int nelemn, 
  int nnodes, std::vector<std::vector<int>> &node, 
  int n_points, std::vector<double> &r, 
  std::vector<double> &uprof, 
  std::vector<double> &xc, double xprof, 
  std::vector<double> &yc
);
std::string file_name_inc(std::string file_name);
void xy_write(
  std::ofstream& xy_file_obj, int n_points, 
  std::vector<double>& xc, std::vector<double>& yc
);
void uv_write(
  const std::vector<double>& f, 
  const std::vector<std::vector<int>>& indx, 
  std::ofstream& uv_file_obj, 
  int neqn, int n_points, 
  const std::vector<double>& yc
);
void linsys(
  std::vector<std::vector<double>> &a, 
  std::vector<double> &area,
  std::vector<double> &f, std::vector<double> &g,
  std::vector<std::vector<int>> &indx,
  std::vector<int> &insc, std::vector<int> &isotri,
  int itype, int maxrow, int nband,
  int nelemn, int neqn, int nlband,
  int nnodes, std::vector<std::vector<int>> &node,
  int n_points, int nquad, int nrow,
  std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
  std::vector<std::vector<std::vector<double>>> &psi, double reynld,
  std::vector<double> &xc, std::vector<std::vector<double>> &xm,
  std::vector<double> &yc, std::vector<std::vector<double>> &ym
);
//* -------------------------------------------------------------------

int main(void)
{
  constexpr int maxnew = 4;
  constexpr int maxsec = 10;
  constexpr int nx = 21;
  constexpr int ny = 7;
  constexpr int maxrow = 27 * ny;
  constexpr int nelemn = 2 * (nx - 1) * (ny - 1);
  constexpr int mx = 2 * nx - 1;
  constexpr int my = 2 * ny - 1;
  constexpr int maxeqn = 2 * mx * my + nx * ny;
  constexpr int n_points = mx * my;
  constexpr int nnodes = 6;
  constexpr int nquad = 3;

  std::vector<std::vector<double>> a;
  double anew = 0.0;
  double anext = 0.3;
  double aold = 0.0;
  double aprof = 0.25;
  std::vector<double> area(nelemn, 0.0);
  std::vector<double> dcda(my, 0.0);
  std::vector<double> f;
  std::vector<double> g;
  std::vector<std::vector<double>> gr(my, std::vector<double>(my, 0.0));
  std::vector<int> iline(my, 0);
  std::vector<std::vector<int>> indx(n_points, std::vector<int>(2, 0));
  std::vector<int> insc(n_points, 0);
  std::vector<int> isotri(nelemn, 0);
  int iwrite = 0;
  bool _long = false;
  int nband = 0;
  int neqn = 0;
  int nlband = 0;
  std::vector<std::vector<int>> node(nelemn, std::vector<int>(nnodes, 0));
  int nrow = 0;
  int numnew = 0;
  int numsec = 0;
  // phi = np.zeros((nelemn, nquad, nnodes, 3), dtype=np.float64)
  std::vector<std::vector<std::vector<std::vector<double>>>> phi(
    nelemn, std::vector<std::vector<std::vector<double>>>(
        nquad, std::vector<std::vector<double>>(
          nnodes, std::vector<double>(3, 0.0)
        )
    )
  );
  // psi = np.zeros((nelemn, nquad, nnodes), dtype=np.float64)
  std::vector<std::vector<std::vector<double>>> psi(
    nelemn, std::vector<std::vector<double>>(
      nquad, std::vector<double>(nnodes, 0.0)
    )
  );
  std::vector<double> r(my, 0.0);
  std::vector<double> res;
  double reynld = 1.0;
  double rjpnew = 0.0;
  double rjpold = 0.0;
  std::vector<double> sens;
  double tolnew = 0.0001;
  double tolsec = 0.0001;
  std::vector<double> uprof(my, 0.0);
  std::string uv_dir = "data/uv";
  std::string uv_file = std::format("{}/uv_000.txt", uv_dir);
  double xbleft = 1.0;
  double xbrite = 3.0;
  std::vector<double> xc(n_points, 0.0);
  double xlngth = 10.0;
  std::vector<std::vector<double>> xm(nelemn, std::vector<double>(nquad, 0.0));
  double xprof = 4.0;
  std::string xy_dir = "data/xy";
  std::string xy_file = std::format("{}\\xy_000.txt", xy_dir);
  std::vector<double> yc(n_points, 0.0);
  double ylngth = 3.0;
  std::vector<std::vector<double>> ym(nelemn, std::vector<double>(nquad, 0.0));
  double ypert = 0.0;
  int itype;
  double temp, denom, test;

  
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
  std::println("  The bump will be generated with a height of {}", aprof);
  std::println();
  std::println("  NX = {}", nx);
  std::println("  NY = {}", ny);
  std::println("  Number of elements = {}", nelemn);
  std::println("  Reynolds number =  {}", reynld);
  std::println("  Secant tolerance = {}", tolsec);
  std::println("  Newton tolerance = {}", tolnew);

  int ibump = 2;
  std::tie(_long, neqn) = setgrd(
    ibump, indx, insc,
    isotri, iwrite, _long,
    maxeqn, mx, my,
    nelemn, neqn, nnodes,
    node, n_points, nx,
    ny, xbleft, xbrite, xlngth
  );

  f.resize(neqn, 0.0);
  g.resize(neqn, 0.0);
  res.resize(neqn, 0.0);
  sens.resize(neqn, 0.0);


  ypert = aprof;
  setxy(iwrite, _long, mx, my, n_points, nx, ny, xc, xlngth, yc, ylngth, ypert);

  setqud(area, isotri, iwrite, nelemn, nnodes, node, n_points, nquad, xc, xm, yc, ym);

  setbas(isotri, nelemn, nnodes, node, n_points, nquad, phi, psi, xc, xm, yc, ym);

  setlin(iline, indx, iwrite, _long, mx, my, n_points, nx, ny, xlngth, xprof);

  std::tie(nband, nlband, nrow) = setban(
    indx, insc, maxrow, nband, nelemn, nlband, nnodes, node, n_points, nrow
  );

  a.resize(nrow, std::vector<double>(neqn, 0.0));

  numnew = nstoke(
    a, area, f, g, indx, insc, isotri,
    maxnew, maxrow, nband,
    nelemn, neqn, nlband,
    nnodes, node, n_points, nquad, nrow,
    numnew, phi, psi, reynld, tolnew,
    xc, xm, yc, ym
  );
  
  resid(
    area, g, indx, insc, isotri,
    iwrite, nelemn, neqn, nnodes, node, n_points,
    nquad, phi, psi, res, reynld, xc,
    xm, yc, ym
  );

  uprof = getg(g, iline, my, neqn);

  if (1 <= iwrite) {
    std::println();
    std::println("Velocity profile:");
    std::println();
    for (int i = 0; i < my; i++) {
      std::print("{}\t", uprof[i]);
      if ((i + 1) % 5 == 0)
        std::println();
    }
    std::println();
  }

  gram(
    gr, iline, indx, iwrite, my, nelemn, 
    nnodes, node, n_points, r, uprof, xc, xprof, yc
  );

  xy_file = file_name_inc(xy_file);

  std::ofstream f_xy(xy_file); // Abre el archivo en modo escritura
  if (f_xy.is_open()) {
    xy_write(f_xy, n_points, xc, yc);
    f_xy.close();
  }

  uv_file = file_name_inc(uv_file);
  //uv_unit = get_unit()
  std::ofstream f_uv(uv_file);
  if (f_uv.is_open()) {
    // Llamada a uv_write con los parámetros correctos
    uv_write(f, indx, f_uv, neqn, n_points, yc);
    f_uv.close();
  }

  //g[:neqn] = 0.0
  for(int i = 0; i < neqn; i++) {
    g[i] = 0.0;
  }
  
  int iter;
  for (iter = 1; iter <= maxsec; iter++) {
    std::println();
    std::println("Secant iteration {}", iter);

    numsec += 1;

    ypert = anew;
    setxy(iwrite, _long, mx, my, n_points, nx, ny, xc, xlngth, yc, ylngth, ypert);

    setqud(area, isotri, iwrite, nelemn, nnodes, node, n_points, nquad, xc, xm, yc, ym);

    setbas(isotri, nelemn, nnodes, node, n_points, nquad, phi, psi, xc, xm, yc, ym);

    numnew = nstoke(
      a, area, f, g, indx, insc,
      isotri, maxnew, maxrow, nband,
      nelemn, neqn, nlband, nnodes, node, 
      n_points, nquad, nrow, numnew, phi,
      psi, reynld, tolnew, xc,
      xm, yc, ym
    );

    uprof = getg(g, iline, my, neqn);

    if (1 <= iwrite) {
      std::println();
      std::println("Velocity profile:");
      std::println();
      for (int i = 0; i < my; i++) {
        std::print("{}\t", uprof[i]);
        if ((i + 1) % 5 == 0) std::println();
      }
      std::println();
    }

    itype = -2;

    linsys(
      a, area, sens, g, indx, insc, isotri, itype, maxrow,
      nband, nelemn, neqn, nlband, nnodes, node,
      n_points, nquad, nrow, phi, psi,
      reynld, xc, xm, yc, ym
    );

    dcda = getg(sens, iline, my, neqn);

    if (2 <= iwrite) {
      std::println();
      std::println("Sensitivities:");
      std::println();
      for (int i = 0; i < my; i++) {
        std::print("{}\t", dcda[i]);
        if ((i + 1) % 5 == 0) {
          std::println();
        }
      }
      std::println();
    }

    rjpnew = 0.0;
    for (int i = 0; i < my; i++) {
      temp = -r[i];
      for (int j = 0; j < my; j++) {
        temp += gr[i][j] * uprof[j];
      }
      rjpnew += 2.0 * dcda[i] * temp;
    }

    xy_file = file_name_inc(xy_file);

    std::ofstream f_xy(xy_file); // Abre el archivo en modo escritura
    if (f_xy.is_open()) {
      xy_write(f_xy, n_points, xc, yc);
      f_xy.close();
    }

    uv_file = file_name_inc(uv_file);
    //uv_unit = get_unit()
    std::ofstream f_uv(uv_file);
    if (f_uv.is_open()) {
      // Llamada a uv_write con los parámetros correctos
      uv_write(f, indx, f_uv, neqn, n_points, yc);
      f_uv.close();
    }

    std::println();
    std::println("  Parameter = {}, J prime = {}", anew, rjpnew);

    if (1 < iter) {
      denom = rjpnew - rjpold;
      if (std::abs(denom) > 1e-30) 
        anext = aold - rjpold * (anew - aold) / denom;
      else anext = anew;
    }

    aold = anew;
    anew = anext;
    rjpold = rjpnew;

    if (anew != 0.0) test = std::abs(anew - aold) / anew;
    else test = 0.0;

    std::println("  New value of parameter = {}", anew);
    std::println("  Convergence test = {}", test);

    if ((std::abs(anew - aold) <= std::abs(anew) * tolsec) and (1 < iter)) {
      std::println("Secant iteration converged.");
      break;
    }
  }
  if (maxsec < iter)
    std::println("  Secant iteration failed to converge.");

  auto fin = std::chrono::high_resolution_clock::now();
  auto duracion = std::chrono::duration_cast<std::chrono::milliseconds>(fin - inicio);

  
  std::println();
  std::println("  Total execution time = {} ms.", duracion.count());
  std::println("  Number of secant steps = {}", numsec);
  std::println("  Number of Newton steps = {}", numnew);
  std::println();
  std::println("BUMP:");
  std::println("  Normal end of execution.");
  std::println();
  timestamp();
  
  if (save_times) {
    std::ofstream archivo("times.txt", std::ios_base::app);
    archivo << duracion.count() << std::endl;
    archivo.close();
  }
  return 0;
}

//* --------------------------------------------------------------------
//*  BSP - linear basis function for pressure
//* --------------------------------------------------------------------
double bsp(
  int it, int iq, int id, int nelemn, int nnodes, 
  const std::vector<std::vector<int>> &node, int n_points, 
  const std::vector<double> &xc, double xq, 
  const std::vector<double> &yc, double yq
) {
  int l1 = iq;
  int l2 = i4_wrap(iq + 1, 0, 2);  // 0-indexed: 1->2, 2->0, 3->1
  int l3 = i4_wrap(iq + 2, 0, 2);

  int g1 = node[it][l1];
  int g2 = node[it][l2];
  int g3 = node[it][l3];

  double d = (
    (xc[g2] - xc[g1]) * (yc[g3] - yc[g1]) 
    - (xc[g3] - xc[g1]) * (yc[g2] - yc[g1])
  );

  if (id == 0)  // 0-indexed: id == 1 -> 0
    return (
      1.0 + ((yc[g2] - yc[g3]) * (xq - xc[g1]) + (xc[g3] - xc[g2]) * (yq - yc[g1])) / d
    );
  else if (id == 1)
    return (yc[g2] - yc[g3]) / d;
  else if (id == 2)
    return (xc[g3] - xc[g2]) / d;
  else {
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
  int n, double da, 
  std::vector<double>& dx, int incx, 
  std::vector<double>& dy, int incy
) {
  if (n <= 0 or da == 0.0) return;
  int ix, iy, m;

  if (incx != 1 or incy != 1) {
    if (incx >= 0) ix = 0;
    else ix = (-n + 1) * incx;
    
    if (incy >= 0) iy = 0;
    else iy = (-n + 1) * incy;

    for (int i = 0; i < n; i++) {
      dy[iy] += da * dx[ix];
      ix += incx;
      iy += incy;
    }
  } else {
    m = n % 4;
    for (int i = 0; i < m; i++) {
      dy[i] += da * dx[i];
    }
    for (int i = m; i < n; i += 4) {
      dy[i] += da * dx[i];
      dy[i + 1] += da * dx[i + 1];
      dy[i + 2] += da * dx[i + 2];
      dy[i + 3] += da * dx[i + 3];
    }
  }
}
void daxpy_m(
  int n, double da, 
  std::vector<std::vector<double>> &abd_x, 
  int col_x, int start_row_x, int incx,
  std::vector<std::vector<double>>& abd_y, 
  int col_y, int start_row_y, int incy
) {
  if (n <= 0 || da == 0.0) return;

  for (int i = 0; i < n; ++i) {
    abd_y[start_row_y + (i * incy)][col_y] += da * abd_x[start_row_x + (i * incx)][col_x];
  }
}
void daxpy(
  int n, double da, 
  const std::vector<std::vector<double>>& abd_x, int col_x, int start_row_x, int incx,
  std::vector<double>& b_y, int start_idx_y, int incy
) {
    
    if (n <= 0 || da == 0.0) return;

    for (int i = 0; i < n; ++i) {
        b_y[start_idx_y + (i * incy)] += da * abd_x[start_row_x + (i * incx)][col_x];
    }
}

//* --------------------------------------------------------------------
//*  DDOT - dot product of two vectors
//* --------------------------------------------------------------------
double ddot(
  int n, 
  const std::vector<std::vector<double>> &abd, 
  int col, int row_start, int incx, 
  const std::vector<double> &b, int start_b, int incy
) {
    
  if (n <= 0) return 0.0;

  double dtemp = 0.0;

  // Caso cuando los incrementos no son unitarios
  if (incx != 1 || incy != 1) {
    int ix = row_start;
    int iy = start_b;
    for (int i = 0; i < n; ++i) {
      dtemp += abd[ix][col] * b[iy];
      ix += incx;
      iy += incy;
    }
  } 
  // Caso con incremento unitario (Optimización manual)
  else {
    int m = n % 5;
    // Procesar el resto
    for (int i = 0; i < m; ++i) {
      dtemp += abd[row_start + i][col] * b[start_b + i];
    }
    // Procesar en bloques de 5
    for (int i = m; i < n; i += 5) {
      dtemp += (
        abd[row_start + i][col] * b[start_b + i]
        + abd[row_start + i + 1][col] * b[start_b + i + 1]
        + abd[row_start + i + 2][col] * b[start_b + i + 2]
        + abd[row_start + i + 3][col] * b[start_b + i + 3]
        + abd[row_start + i + 4][col] * b[start_b + i + 4]
      );
    }
  }

  return dtemp;
}

//* --------------------------------------------------------------------
//*  DGBFA - factor a real band matrix by elimination
//* --------------------------------------------------------------------
std::tuple<int, std::vector<int>> dgbfa(
  std::vector<std::vector<double>> &abd, int lda, 
  int n, int ml, int mu, 
  std::vector<int> &ipvt
) {
  int m = ml + mu + 1;
  int info = 0;

  int j0 = mu + 1;  // 0-indexed adjustment
  int j1 = std::min(n, m) - 1;
  int i0;
  double t;

  for (int jz = j0; jz < j1; jz++) {
    i0 = m - jz;
    for (int i = i0; i < ml; i++) {
      abd[i][jz] = 0.0;
    }
  }

  int jz = j1;
  int ju = 0, lm, l, mm;

  for (int k = 0; k < n - 1; k++) {
    jz += 1;
    if (jz < n)
      for (int i = 0; i < ml; i++) {
        abd[i][jz] = 0.0;
      }

    lm = std::min(ml, n - k - 1);
    l = idamax_m(lm + 1, abd, k, m - 1, 1) + (m - 1);
    ipvt[k] = l + k - (m - 1);

    if (abd[l][k] == 0.0)
      info = k;
    else {
      if (l != m - 1) {
        double t = abd[l][k];
        abd[l][k] = abd[m - 1][k];
        abd[m - 1][k] = t;
      }

      t = -1.0 / abd[m - 1][k];
      dscal_m(lm, t, abd, k, m, 1);

      ju = std::min(std::max(ju, mu + ipvt[k]), n - 1);
      mm = m - 1;

      for (int j = k + 1; j < ju + 1; j++) {
        l -= 1;
        mm -= 1;
        double t = abd[l][j];
        if (l != mm) {
          abd[l][j] = abd[mm][j];
          abd[mm][j] = t;
        }
        daxpy_m(lm, t, abd, k, m, 1, abd, j, mm + 1, 1);
      }
    }
  }

  ipvt[n - 1] = n - 1;

  if (abd[m - 1][n - 1] == 0.0)
    info = n - 1;

  return std::tuple {info, ipvt};
}

//* --------------------------------------------------------------------
//*  DGBSL - solve a real banded system factored by DGBFA
//* --------------------------------------------------------------------
void dgbsl(
  std::vector<std::vector<double>> &abd, int lda,
  int n, int ml, int mu,
  std::vector<int> &ipvt, 
  std::vector<double> &b, int job
) {
  int m = mu + ml + 1;
  double t;
  int l, lm, la, lb;

  if (job == 0) {
    if (0 < ml) {
      for (int k = 0; k < n - 1; k++) {
        lm = std::min(ml, n - k - 1);
        l = ipvt[k];
        t = b[l];
        if (l != k) {
          b[l] = b[k];
          b[k] = t;
        }
        daxpy(lm, t, abd, k, m, 1, b, k + 1, 1);
      }
    }

    for (int k = n - 1; k >= 0; k--) {
      b[k] /= abd[m - 1][k];
      lm = std::min(k, m - 1);
      la = m - 1 - lm;
      lb = k - lm;
      t = -b[k];
      daxpy(lm, t, abd, k, la, 1, b, lb, 1);
    }
  } else {
    for (int k = 0; k < n; k++) {
      lm = std::min(k, m - 1);
      la = m - 1 - lm;
      lb = k - lm;
      t = ddot(lm, abd, k, la, 1, b, lb, 1);
      b[k] = (b[k] - t) / abd[m - 1][k];
    }

    if (0 < ml) {
      for (int k = n - 2; k >= 0; k--) {
        lm = std::min(ml, n - k - 1);
        b[k] += ddot(lm, abd, k, m, 1, b, k + 1, 1);
        l = ipvt[k];
        if (l != k) {
          t = b[l];
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
void dscal_v(
  int n, double sa, 
  std::vector<double>& x, int incx
) {
  int ix, m;
  if (n <= 0) return;

  if (incx == 1) {
    m = n % 5;
    for (int i = 0; i < m; i++) {
      x[i] *= sa;
    }
    for (int i = m; i < n; i += 5) {
      x[i] *= sa;
      x[i + 1] *= sa;
      x[i + 2] *= sa;
      x[i + 3] *= sa;
      x[i + 4] *= sa;
    }
  } else {
    if (incx >= 0) ix = 0;
    else ix = (-n + 1) * incx;
    for (int i = 0; i < n; i++) {
      x[ix] *= sa;
      ix += incx;
    }
  }
}
void dscal_m(int n, double sa, std::vector<std::vector<double>>& abd, int col, int start_row, int incx) {
  if (n <= 0) return;

  for (int i = 0; i < n; ++i) {
    abd[start_row + (i * incx)][col] *= sa;
  }
}
//* --------------------------------------------------------------------
//*  FILE_NAME_INC - increment a partially numeric filename
//* --------------------------------------------------------------------
std::string file_name_inc(std::string file_name) {
  if (file_name.empty()) {
    std::println("\nFILE_NAME_INC - Fatal error!");
    std::println("The input string is empty.");
    std::exit(1);
  }

  int change = 0;
  // Recorremos el string desde el final hacia el inicio
  for (int i = file_name.length() - 1; i >= 0; --i) {
    if (file_name[i] >= '0' && file_name[i] <= '9') {
      change++;
      int digit = file_name[i] - '0';
      digit++;
      
      if (digit == 10) {
        digit = 0;
        file_name[i] = '0';
      } else {
        file_name[i] = static_cast<char>(digit + '0');
        return file_name; // Retorno temprano si no hay acarreo
      }
  }
  }

  return (change == 0) ? "" : file_name;
}

//* --------------------------------------------------------------------
//*  GETG - extract values of a quantity along the profile line
//* --------------------------------------------------------------------
std::vector<double> getg(
  const std::vector<double> &f, const std::vector<int> &iline, int my, int neqn
) {
  std::vector<double> u(my, 0.0);
  int j;
  for (int i = 0; i < my; i++) {
    j = iline[i];
    if (j < 0)  // 0-indexed: j <= 0 means j == -1 (Fortran: j <= 0)
      u[i] = 0.0;
    else if (j == 0) u[i] = 0.0;
    else u[i] = f[j - 1];  // Convert 1-based index from iline to 0-based
  }
  return u;
}

//* --------------------------------------------------------------------
//*  GRAM - compute the Gram matrix and R vector
//* --------------------------------------------------------------------
void gram(
  std::vector<std::vector<double>> &gr, 
  std::vector<int> &iline, 
  std::vector<std::vector<int>> &indx, 
  int iwrite, int my, int nelemn, 
  int nnodes, std::vector<std::vector<int>> &node, 
  int n_points, std::vector<double> &r, 
  std::vector<double> &uprof, 
  std::vector<double> &xc, double xprof, 
  std::vector<double> &yc
) {
  constexpr std::array<double, 3> wt = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  constexpr std::array<double, 3> yq_gauss = {-0.7745966692, 0.0, 0.7745966692};

  for (int it = 0; it < static_cast<int>(r.size()); it++) {
    r[it] = 0.0;
  }

  for (int it = 0; it < static_cast<int>(gr.size()); it++) {
    for (int j = 0; j < static_cast<int>(gr[it].size()); j++) {
      gr[it][j] = 0.0;
    }
  }
  for (int it = 0; it < nelemn; it++) {
    int k = node[it][0];
    int kk = node[it][1];

    if ((std::abs(xc[k] - xprof) > 1.0e-4) or (std::abs(xc[kk] - xprof) > 1.0e-4))
      continue;

    for (int iquad = 0; iquad < 3; iquad++) {
      double bma2 = (yc[kk] - yc[k]) / 2.0;
      double ar = bma2 * wt[iquad];
      double x = xprof;
      double y = yc[k] + bma2 * (yq_gauss[iquad] + 1.0);

      double uiqdpt = 0.0;
      for (int iq = 0; iq < nnodes; iq++) {
        if ((iq == 0) or (iq == 1) or (iq == 3)) {
          auto [bb, bx, by] = qbf(x, y, it, iq, nelemn, nnodes, node, n_points, xc, yc);
          int ip = node[it][iq];
          int iun = indx[ip][0];
          if (0 < iun) {
            int ii = igetl(iun, iline, my);
            uiqdpt += bb * uprof[ii - 1];
          } else if (iun == -1) {
            double ubc = ubdry(1, yc[ip]);
            uiqdpt += bb * ubc;
          }
        }
      }

      for (int iq = 0; iq < nnodes; iq++) {
        if ((iq == 0) or (iq == 1) or (iq == 3)) {
          int ip = node[it][iq];
          auto [bb, bx, by] = qbf(x, y, it, iq, nelemn, nnodes, node, n_points, xc, yc);
          int i_val = indx[ip][0];
          if (0 < i_val) {
            int ii = igetl(i_val, iline, my);
            r[ii - 1] += bb * uiqdpt * ar;
            for (int iqq = 0; iqq < nnodes; iqq++) {
              if (iqq == 0 or iqq == 1 or iqq == 3) {
                int ipp = node[it][iqq];
                auto [bbb, bbx_, bby_] = qbf(
                  x, y, it, iqq, nelemn, nnodes, node, n_points, xc, yc
                );
                int j_val = indx[ipp][0];
                if (j_val != 0) {
                  int jj = igetl(j_val, iline, my);
                  gr[ii - 1][jj - 1] += bb * bbb * ar;
                }
              }
            }
          }
        }
      }
    }
  }
  if (3 <= iwrite) {
    std::println();
    std::println("Gram matrix:");
    std::println();
    for (int i = 0; i < my; i++) {
      for (int j = 0; j < my; j++) {
        std::println("{} {} {}", i + 1, j + 1, gr[i][j]);
      }
    }
    std::println();
    std::println("R vector:");
    std::println();
    for (int i = 0; i < my; i++) {
      std::println("{}", r[i]);
    }
  }
}

//* --------------------------------------------------------------------
//*  I4_MODP - nonnegative remainder of integer division
//* --------------------------------------------------------------------
inline int i4_modp(int i, int j) {
  if (j == 0) {
    std::println();
    std::println("I4_MODP - Fatal error!");
    std::println("  Illegal divisor J = {}", j);
    std::exit(1);
  }
  int value = i % j;
  if (value < 0)
    value += std::abs(j);
  return value;
}

//* --------------------------------------------------------------------
//*  I4_WRAP - force an integer between limits by wrapping
//* --------------------------------------------------------------------
inline int i4_wrap(int ival, int ilo, int ihi) {
  int jlo = std::min(ilo, ihi);
  int jhi = std::max(ilo, ihi);
  int wide = jhi - jlo + 1;
  if (wide == 1)
    return jlo;
  return jlo + i4_modp(ival - jlo, wide);
}

//* --------------------------------------------------------------------
//*  IDAMAX - index of the vector element of maximum absolute value
//* --------------------------------------------------------------------
int idamax_v(int n, std::vector<double> &dx, int incx) {
  double dmax;
  int ix;

  if (n < 1 or incx <= 0)
    return -1;

  int idamax_val = 0;

  if (n == 1)
    return 0;

  if (incx == 1) {
    dmax = std::abs(dx[0]);
    for (int i = 1; i < n; i++) {
      if (dmax < std::abs(dx[i])) {
        idamax_val = i;
        dmax = std::abs(dx[i]);
      }
    }
  } else {
    ix = 0;
    dmax = std::abs(dx[0]);
    ix += incx;
    for (int i = 1; i < n; i++) {
      if (dmax < std::abs(dx[ix])) {
        idamax_val = i;
        dmax = std::abs(dx[ix]);
      }
      ix += incx;
    }
  }

  return idamax_val;
}

int idamax_m(int n, const std::vector<std::vector<double>>& abd, int col, int start_row, int incx) {
  if (n < 1 || incx <= 0) return -1;

  int idamax_val = 0;
  double dmax = std::abs(abd[start_row][col]);
  
  // Recorremos desde start_row hasta start_row + n
  for (int i = 1; i < n; i++) {
    double current_val = std::abs(abd[start_row + (i*incx)][col]);
    if (current_val > dmax) {
      idamax_val = i;
      dmax = current_val;
    }
  }
  return idamax_val;
}

//* --------------------------------------------------------------------
//*  IGETL - get the local unknown number along the profile line
//* --------------------------------------------------------------------
inline int igetl(int i, std::vector<int> &iline, int my) {
  for (int j = 0; j < my; j++) {
    if (iline[j] == i) return j + 1; // 1-indexed
  }
  return -1;
}

//* --------------------------------------------------------------------
//*  LINSYS - solve the linearized Navier Stokes equation
//* --------------------------------------------------------------------
void linsys(
  std::vector<std::vector<double>> &a, 
  std::vector<double> &area,
  std::vector<double> &f, std::vector<double> &g,
  std::vector<std::vector<int>> &indx,
  std::vector<int> &insc, std::vector<int> &isotri,
  int itype, int maxrow, int nband,
  int nelemn, int neqn, int nlband,
  int nnodes, std::vector<std::vector<int>> &node,
  int n_points, int nquad, int nrow,
  std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
  std::vector<std::vector<std::vector<double>>> &psi, double reynld,
  std::vector<double> &xc, std::vector<std::vector<double>> &xm,
  std::vector<double> &yc, std::vector<std::vector<double>> &ym
) {
  int ioff = (
    nlband + nlband + 1
  );  // This is 1-indexed offset to main diagonal row in band storage
  double visc = 1.0 / reynld;
  int ip, ihor, iver, iprs;
  int ipp;
  double bbb, bbx, bby;
  int ju, jv, jp, iuse, info;

  for (int i = 0; i < neqn; i++) {
    f[i] = 0.0;
  }

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < neqn; j++) {
      a[i][j] = 0.0;
    }
  }

  std::vector<double> un(2, 0.0);
  std::vector<double> unx(2, 0.0);
  std::vector<double> uny(2, 0.0);
  double det, etax, etay, xix, xiy;
  double ar, yq, xq, bbl, bbbl, ubc, aij;
  double bb, bx, by;

  for (int it = 0; it < nelemn; it++) {
    ar = area[it] / 3.0;

    for (int iquad = 0; iquad < nquad; iquad++) {
      yq = ym[it][iquad];
      xq = xm[it][iquad];

      det = 0.0;
      etax = 0.0;
      etay = 0.0;
      xix = 0.0;
      xiy = 0.0;
      ubc = 0.0;

      if (isotri[it] == 1) {
        std::tie(det, etax, etay, xix, xiy) = trans(
          it, nelemn, nnodes, node, n_points, xc, xq, yc, yq
        );
        ar = det * area[it] / 3.0;
      }

      std::tie(un, unx, uny) = uval(
        etax, etay, g, indx, isotri, it, nelemn,
        neqn, nnodes, node, n_points, xc, xix,
        xiy, xq, yc, yq
      );

      for (int iq = 0; iq < nnodes; iq++) {
        ip = node[it][iq];
        bb = phi[it][iquad][iq][0];
        bx = phi[it][iquad][iq][1];
        by = phi[it][iquad][iq][2];
        bbl = psi[it][iquad][iq];
        ihor = indx[ip][0];
        iver = indx[ip][1];
        iprs = insc[ip];

        if (0 < ihor)
          f[ihor - 1] += ar * bb * (un[0] * unx[0] + un[1] * uny[0]);

        if (0 < iver)
          f[iver - 1] += ar * bb * (un[0] * unx[1] + un[1] * uny[1]);

        for (int iqq = 0; iqq < nnodes; iqq++) {
          ipp = node[it][iqq];
          bbb = phi[it][iquad][iqq][0];
          bbx = phi[it][iquad][iqq][1];
          bby = phi[it][iquad][iqq][2];
          bbbl = psi[it][iquad][iqq];
          ju = indx[ipp][0];
          jv = indx[ipp][1];
          jp = insc[ipp];

          if (0 < ju) {
            if (0 < ihor) {
              iuse = ihor - ju + ioff;
              a[iuse - 1][ju - 1] += ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
              );
            }

            if (0 < iver) {
              iuse = iver - ju + ioff;
              a[iuse - 1][ju - 1] += ar * bb * bbb * unx[1];
            }

            if (0 < iprs) {
              iuse = iprs - ju + ioff;
              a[iuse - 1][ju - 1] += ar * bbx * bbl;
            }
          } else if (ju == itype) {
            if (ju == -1)
              ubc = ubdry(1, yc[ipp]);
            else if (ju == -2) {
              ubc = ubump(
                g,
                indx,
                ipp, iqq,
                isotri,
                it, 1, nelemn,
                neqn, nnodes,
                node, n_points,
                xc, yc
              );
            }

            if (0 < ihor) {
              aij = ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
              );
              f[ihor - 1] -= ubc * aij;
            }

            if (0 < iver) {
              aij = ar * bb * bbb * unx[1];
              f[iver - 1] -= ubc * aij;
            }

            if (0 < iprs) {
              aij = ar * bbx * bbl;
              f[iprs - 1] -= ubc * aij;
            }
          }

          if (0 < jv) {
            if (0 < ihor) {
              iuse = ihor - jv + ioff;
              a[iuse - 1][jv - 1] += ar * bb * bbb * uny[0];
            }

            if (0 < iver) {
              iuse = iver - jv + ioff;
              a[iuse - 1][jv - 1] += ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
              );
            }

            if (0 < iprs) {
              iuse = iprs - jv + ioff;
              a[iuse - 1][jv - 1] += ar * bby * bbl;
            }

          } else if (jv == itype) {
            if (jv == -1)
              ubc = ubdry(2, yc[ipp]);
            else if (jv == -2) {
              ubc = ubump(
                g,
                indx,
                ipp,
                iqq,
                isotri,
                it,
                2,
                nelemn,
                neqn,
                nnodes,
                node,
                n_points,
                xc,
                yc
              );
            }

            if (0 < ihor) {
              aij = ar * bb * bbb * uny[0];
              f[ihor - 1] -= ubc * aij;
            }

            if (0 < iver) {
              aij = ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
              );
              f[iver - 1] -= ubc * aij;
            }

            if (0 < iprs) {
              aij = ar * bby * bbl;
              f[iprs - 1] -= ubc * aij;
            }
          }
          if (0 < jp) {
            if (0 < ihor) {
              iuse = ihor - jp + ioff;
              a[iuse - 1][jp - 1] -= ar * bx * bbbl;
            }

            if (0 < iver) {
              iuse = iver - jp + ioff;
              a[iuse - 1][jp - 1] -= ar * by * bbbl;
            }
          }
        }
      }
    }
  }

  f[neqn - 1] = 0.0;
  int i, j;
  for (int j_1based = neqn - nlband; j_1based < neqn; j_1based++) {
    j = j_1based - 1;
    i = neqn - j_1based + ioff;  // 1-indexed row in band storage
    a[i - 1][j] = 0.0;
  }
  a[ioff - 1][neqn - 1] = 1.0;  // nband == ioff, set diagonal of last eqn to 1

  // Save right-hand side in f, solve in-place
  // DGBFA and DGBSL work on A in LINPACK band format

  // We need to call dgbfa with the right parameters
  // a is (maxrow, neqn) but the band part is (nrow, neqn)
  // In LINPACK, LDA = maxrow, N = neqn, ML = nlband, MU = nlband
  //ipvt = np.zeros(neqn, dtype=np.int32)
  std::vector<int> ipvt(neqn, 0);
  std::tie(info, ipvt) = dgbfa(a, maxrow, neqn, nlband, nlband, ipvt);

  if (info != 0) {
    std::println();
    std::println("LINSYS - fatal error!");
    std::println("DGBFA returns INFO = {}", info+1);  // +1 for 1-indexed
    std::exit(1);
  }

  int job = 0;
  dgbsl(a, maxrow, neqn, nlband, nlband, ipvt, f, job);
}


//* --------------------------------------------------------------------
//*  NSTOKE - solve Navier Stokes using Newton iteration
//* --------------------------------------------------------------------
int nstoke(
  std::vector<std::vector<double>> &a,
  std::vector<double> &area,
  std::vector<double> &f, std::vector<double> &g,
  std::vector<std::vector<int>> &indx,
  std::vector<int> &insc, std::vector<int> &isotri,
  int maxnew, int maxrow, int nband,
  int nelemn, int neqn, int nlband,
  int nnodes, std::vector<std::vector<int>> &node,
  int n_points, int nquad, int nrow, int numnew,
  std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
  std::vector<std::vector<std::vector<double>>> &psi,
  double reynld, double tolnew,
  std::vector<double> &xc, std::vector<std::vector<double>> &xm,
  std::vector<double> &yc, std::vector<std::vector<double>> &ym
) {
  int itype, imax;
  double diff;

  for (int iter = 0; iter < maxnew; iter++) {
    numnew++;

    itype = -1;
    linsys(
      a, area,
      f, g,
      indx,
      insc,
      isotri,
      itype, maxrow, nband,
      nelemn, neqn, nlband,
      nnodes, node,
      n_points, nquad, nrow,
      phi,
      psi, reynld,
      xc, xm,
      yc, ym
    );

    for (int i = 0; i < neqn; i++) {
      g[i] -= f[i];
    }

    imax = idamax_v(neqn, g, 1);
    
    diff = std::abs(g[imax]);
    std::println("NSTOKE: Iteration {}, MaxNorm(diff) = {}", iter+1, diff);

    for (int i = 0; i < neqn; i++) {
      g[i] = f[i];
    }

    if (diff <= tolnew) {
      std::println("NSTOKE converged.");
      return numnew;
    }

    if (iter == maxnew - 1) {
      std::println("NSTOKE failed!");
      std::exit(1);
    }
  }

  return numnew;
}

//// --------------------------------------------------------------------
////  PVAL - compute pressure table at all nodes
//// --------------------------------------------------------------------
//def pval(g, insc, _long, mx, my, nelemn, neqn, nnodes, node, n_points):
//  press = np.zeros((mx, my), dtype=np.float64)
//
//  for it in range(nelemn):
//    for iq in range(3):
//      ip = node[it, iq]
//      ivar = insc[ip]
//      if _long:
//        i = (ip) // my
//        j = (ip) % my
//      else:
//        i = (ip) % mx
//        j = (ip) // mx
//
//      if 0 < ivar:
//        press[i, j] = g[ivar - 1]
//
//  for i in range(1, mx - 1, 2):
//    for j in range(0, my, 2):
//      press[i, j] = 0.5 * (press[i - 1, j] + press[i + 1, j])
//
//  for j in range(1, my - 1, 2):
//    for i in range(0, mx, 2):
//      press[i, j] = 0.5 * (press[i, j - 1] + press[i, j + 1])
//
//  for j in range(1, my - 1, 2):
//    for i in range(1, mx - 1, 2):
//      press[i, j] = 0.5 * (press[i - 1, j - 1] + press[i + 1, j + 1])
//
//  return press

//* --------------------------------------------------------------------
//*  QBF - evaluate quadratic basis functions
//* --------------------------------------------------------------------
std::tuple<double, double, double> qbf(
  double xq, double yq, int it, 
  int inn, int nelemn, int nnodes, 
  const std::vector<std::vector<int>> &node, int n_points, 
  const std::vector<double> &xc, const std::vector<double> &yc
) {
  int in1, in2, in3, i1, i2, i3, inn_local;
  int j1, j2, j3;
  double d, c, t, s, bb, bx, by;

  if (inn <= 2) {  // 0-indexed: in=0,1,2 corresponds to Fortran in=1,2,3
    in1 = inn;
    in2 = (inn + 1) % 3;
    in3 = (inn + 2) % 3;
    i1 = node[it][in1];
    i2 = node[it][in2];
    i3 = node[it][in3];
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    t = (
      1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d
    );
    bb = t * (2.0 * t - 1.0);
    bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d;
    by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d;
  } else {
    inn_local = inn - 3;
    in1 = inn_local;
    in2 = (inn_local + 1) % 3;
    in3 = (inn_local + 2) % 3;
    i1 = node[it][in1];
    i2 = node[it][in2];
    i3 = node[it][in3];
    j1 = i2;
    j2 = i3;
    j3 = i1;
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    c = (xc[j2] - xc[j1]) * (yc[j3] - yc[j1]) - (xc[j3] - xc[j1]) * (yc[j2] - yc[j1]);
    t = (
      1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d
    );
    s = (
      1.0 + ((yc[j2] - yc[j3]) * (xq - xc[j1]) + (xc[j3] - xc[j2]) * (yq - yc[j1])) / c
    );
    bb = 4.0 * s * t;
    bx = 4.0 * (t * (yc[j2] - yc[j3]) / c + s * (yc[i2] - yc[i3]) / d);
    by = 4.0 * (t * (xc[j3] - xc[j2]) / c + s * (xc[i3] - xc[i2]) / d);
  }
  return std::tuple {bb, bx, by};
}

//* --------------------------------------------------------------------
//*  REFBSP - evaluate linear basis functions in a reference triangle
//* --------------------------------------------------------------------
inline double refbsp(double xq, double yq, int iq){
  // 0-indexed: iq=0 -> psi1, iq=1 -> psi2, iq=2 -> psi3
  if (iq == 0)
    return 1.0 - xq;
  else if (iq == 1)
    return yq;
  else if (iq == 2)
    return xq - yq;
  return 0.0;
}

//* --------------------------------------------------------------------
//*  REFQBF - evaluate quadratic basis functions on reference triangle
//* --------------------------------------------------------------------
std::tuple<double, double, double> refqbf(
  double x, double y, int inn, 
  double etax, double etay, 
  double xix, double xiy
) {
  // 0-indexed: inn=0..5 corresponds to Fortran in=1..6
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

  return std::tuple {bb, bx, by};
}

//* --------------------------------------------------------------------
//*  RESID - compute the residual
//* --------------------------------------------------------------------
void resid(
  std::vector<double> &area, std::vector<double> &g,
  std::vector<std::vector<int>> &indx,
  std::vector<int> &insc, std::vector<int> &isotri,
  int iwrite, int nelemn, int neqn,
  int nnodes, std::vector<std::vector<int>> &node,
  int n_points, int nquad,
  std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
  std::vector<std::vector<std::vector<double>>> &psi, 
  std::vector<double> &res, double reynld, 
  std::vector<double> &xc, std::vector<std::vector<double>> &xm, 
  std::vector<double> &yc, std::vector<std::vector<double>> &ym
) {
  int itype = -1;
  double visc = 1.0 / reynld;

  for (int it = 0; it < neqn; it++) {
    res[it] = 0.0;
  }

  double ubc = 0.0;
  double ar, yq, xq, det, etax, etay, xix, xiy, aij;
  int ip, iprs, ihor, iver, ipp, ju, jv, jp;
  double bb, bx, by, bbl, bbb, bbx, bby, bbbl;

  std::vector<double> un(2, 0.0);
  std::vector<double> unx(2, 0.0);
  std::vector<double> uny(2, 0.0);

  for (int it = 0; it < nelemn; it++) {
    ar = area[it] / 3.0;

    for (int iquad = 0; iquad < nquad; iquad++) {
      yq = ym[it][iquad];
      xq = xm[it][iquad];

      det = 0.0;
      etax = 0.0;
      etay = 0.0;
      xix = 0.0;
      xiy = 0.0;

      if (isotri[it] == 1) {
        std::tie(det, etax, etay, xix, xiy) = trans(
          it, nelemn, nnodes, node, n_points, xc, xq, yc, yq
        );
        ar = det * area[it] / 3.0;
      }

      std::tie(un, unx, uny) = uval(
        etax, etay, g, indx, isotri, it,
        nelemn, neqn, nnodes, node, n_points, xc,
        xix, xiy, xq, yc, yq
      );

      for (int iq = 0; iq < nnodes; iq++) {
        ip = node[it][iq];
        bb = phi[it][iquad][iq][0];
        bx = phi[it][iquad][iq][1];
        by = phi[it][iquad][iq][2];
        bbl = psi[it][iquad][iq];
        iprs = insc[ip];
        ihor = indx[ip][0];
        iver = indx[ip][1];

        if (0 < ihor)
          res[ihor - 1] -= ar * bb * (un[0] * unx[0] + un[1] * uny[0]);

        if (0 < iver)
          res[iver - 1] -= ar * bb * (un[0] * unx[1] + un[1] * uny[1]);

        for (int iqq = 0; iqq < nnodes; iqq++) {
          ipp = node[it][iqq];
          bbb = phi[it][iquad][iqq][0];
          bbx = phi[it][iquad][iqq][1];
          bby = phi[it][iquad][iqq][2];
          bbbl = psi[it][iquad][iqq];
          ju = indx[ipp][0];
          jv = indx[ipp][1];
          jp = insc[ipp];

          if (0 < ju) {
            if (0 < ihor) {
              res[ihor - 1] += (
                ar
                * (
                  visc * (by * bby + bx * bbx)
                  + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
                )
                * g[ju - 1]
              );
            }

            if (0 < iver)
              res[iver - 1] += ar * bb * bbb * unx[1] * g[ju - 1];

            if (0 < iprs)
              res[iprs - 1] += ar * bbx * bbl * g[ju - 1];

          } else if (ju == itype) {
            if (ju == -2) {
              ubc = ubump(
                g, indx, ipp, iqq, isotri, it, 1,
                nelemn, neqn, nnodes, node, n_points,
                xc, yc
              );
            } else if (ju == -1)
              ubc = ubdry(1, yc[ipp]);

            if (0 < ihor) {
              aij = ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
              );
              res[ihor - 1] += ubc * aij;
            }

            if (0 < iver) {
              aij = ar * bb * bbb * unx[1];
              res[iver - 1] += ubc * aij;
            }

            if (0 < iprs) {
              aij = ar * bbx * bbl;
              res[iprs - 1] += ubc * aij;
            }
          }
          if (0 < jv) {
            if (0 < ihor)
              res[ihor - 1] += ar * bb * bbb * uny[0] * g[jv - 1];

            if (0 < iver) {
              res[iver - 1] += (
                ar
                * (
                  visc * (by * bby + bx * bbx)
                  + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
                )
                * g[jv - 1]
              );
            }

            if (0 < iprs)
              res[iprs - 1] += ar * bby * bbl * g[jv - 1];

          } else if (jv == itype) {
            if (jv == -2) {
              ubc = ubump(
                g, indx, ipp, iqq, isotri, it, 2,
                nelemn, neqn, nnodes, node, n_points, 
                xc, yc
              );
            } else if (jv == -1) ubc = ubdry(2, yc[ipp]);

            if (0 < ihor) {
              aij = ar * bb * bbb * uny[0];
              res[ihor - 1] += ubc * aij;
            }

            if (0 < iver) {
              aij = ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
              );
              res[iver - 1] += ubc * aij;
            }

            if (0 < iprs) {
              aij = ar * bby * bbl;
              res[iprs - 1] += ubc * aij;
            }
          }
          if (0 < jp) {
            if (0 < ihor)
              res[ihor - 1] -= ar * bx * bbbl * g[jp - 1];

            if (0 < iver)
              res[iver - 1] -= ar * by * bbbl * g[jp - 1];
          }
        }
      }
    }
  }

  res[neqn - 1] = g[neqn - 1];

  double rmax = 0.0;
  int imax = 0;
  int ibad = 0;
  double test;

  for (int i = 0; i < neqn; i++) {
    test = std::abs(res[i]);
    if (rmax < test) {
      rmax = test;
      imax = i;
    }
    if (1.0e-3 < test)
      ibad += 1;
  }

  if (1 <= iwrite) {
    std::println();
    std::println("RESIDUAL INFORMATION:");
    std::println();
    std::println("Worst residual is number {}", imax+1);
    std::println("of magnitude {}", rmax);
    std::println();
    std::println("Number of \"bad\" residuals is {} out of {}", ibad, neqn);
    std::println();
  }

  if (2 <= iwrite) {
    std::println("Raw residuals:");
    std::println();
    int idx = 0;
    for (int j = 0; j < n_points; j++) {
      if (0 < indx[j][0]) {
        if (std::abs(res[idx]) <= 1.0e-3)
          std::println(" U {} {} {}", idx+1, j+1, res[idx]);
        else
          std::println("*U {} {} {}", idx+1, j+1, res[idx]);
        idx += 1;
      }

      if (0 < indx[j][1]) {
        if (std::abs(res[idx]) <= 1.0e-3)
          std::println(" V {} {} {}", idx+1, j+1, res[idx]);
        else
          std::println("*V {} {} {}", idx+1, j+1, res[idx]);
        idx += 1;
      }

      if (0 < insc[j]) {
        if (std::abs(res[idx]) <= 1.0e-3)
          std::println(" P {} {} {}", idx+1, j+1, res[idx]);
        else
          std::println("*P {} {} {}", idx+1, j+1, res[idx]);
        idx += 1;
      }
    }
  }
}


//* --------------------------------------------------------------------
//*  SETBAN - compute the half band width
//* --------------------------------------------------------------------
std::tuple<int, int, int> setban(
  std::vector<std::vector<int>> &indx, 
  std::vector<int> &insc, int maxrow, 
  int nband, int nelemn, int nlband, 
  int nnodes, std::vector<std::vector<int>> &node, 
  int n_points, int nrow
) {
  int ip, i_val, ipp, j_val;
  nlband = 0;

  for (int it = 0; it < nelemn; it++) {
    for (int iq = 0; iq < nnodes; iq++) {
      ip = node[it][iq];
      for (int iuk = 0; iuk < 3; iuk++) {
        if (iuk == 2)
          i_val = insc[ip];
        else
          i_val = indx[ip][iuk];
        if (0 < i_val) {
          for (int iqq = 0; iqq < nnodes; iqq++) {
            ipp = node[it][iqq];
            for (int iukk = 0; iukk < 3; iukk++) {
              if (iukk == 2)
                j_val = insc[ipp];
              else
                j_val = indx[ipp][iukk];
              if (0 < j_val)
                nlband = std::max(nlband, j_val - i_val);
            }
          }
        }
      }
    }
  }

  nband = nlband + nlband + 1;
  nrow = nlband + nlband + nlband + 1;

  std::println();
  std::println("SETBAN:");
  std::println();
  std::println("  Lower bandwidth = {}", nlband);
  std::println("  Total bandwidth = {}", nband);
  std::println("  Required matrix rows = {}", nrow);

  if (maxrow < nrow) {
    std::println("SETBAN - NROW is too large!");
    std::println("The maximum allowed is {}", maxrow);
    std::exit(1);
  }

  return std::tuple{nband, nlband, nrow};
}

//* --------------------------------------------------------------------
//*  SETBAS - evaluate basis functions at each integration point
//* --------------------------------------------------------------------
void setbas(
  const std::vector<int> &isotri, int nelemn, int nnodes, 
  std::vector<std::vector<int>> &node, int n_points, int nquad, 
  std::vector<std::vector<std::vector<std::vector<double>>>> &phi, 
  std::vector<std::vector<std::vector<double>>> &psi, 
  std::vector<double> &xc, 
  std::vector<std::vector<double>> &xm, 
  std::vector<double> &yc, 
  std::vector<std::vector<double>> &ym
) {
  double xq, yq;
  double det, etax, etay, xix, xiy;
  double bb, bx, by;

  for (int it = 0; it < nelemn; it++) {
    for (int j = 0; j < nquad; j++) {
      xq = xm[it][j];
      yq = ym[it][j];
      std::tie(det, etax, etay, xix, xiy) = trans(
        it, nelemn, nnodes, node, n_points, xc, xq, yc, yq
      );

      for (int iq = 0; iq < nnodes; iq++) {
        if (isotri[it] == 0) {
          psi[it][j][iq] = bsp(
            it, iq, 0, nelemn, nnodes, node, n_points, xc, xq, yc, yq
          );
          std::tie(bb, bx, by) = qbf(xq, yq, it, iq, nelemn, nnodes, node, n_points, xc, yc);
        } else {
          std::tie(bb, bx, by) = refqbf(xq, yq, iq, etax, etay, xix, xiy);
          psi[it][j][iq] = refbsp(xq, yq, iq);
        }
        phi[it][j][iq][0] = bb;
        phi[it][j][iq][1] = bx;
        phi[it][j][iq][2] = by;
      }
    }
  }
}

//* --------------------------------------------------------------------
//*  SETGRD - set up the geometric grid
//* --------------------------------------------------------------------
std::tuple<bool, int> setgrd(
  int ibump, 
  std::vector<std::vector<int>> &indx, 
  std::vector<int> &insc,
  std::vector<int> &isotri, 
  int iwrite, bool _long, int maxeqn,
  int mx, int my, int nelemn,
  int neqn, int nnodes,
  std::vector<std::vector<int>> &node, 
  int n_points, int nx, int ny,
  double xbleft, double xbrite, double xlngth
) {
  std::println();
  std::println("SETGRD:");
  std::println();

  if (ny < nx){
    _long = true;
    std::println("Using vertical ordering.");
  } else{
    _long = false;
    std::println("Using horizontal ordering.");
  }

  if (ibump == 0)
    std::println("No isoparametric elements will be used.");
  else if (ibump == 1)
    std::println("Isoparametric elements directly on bump.");
  else if (ibump == 2)
  std::println("All elements above bump are isoparametric.");
  else if (ibump == 3)
  std::println("All elements are isoparametric.");
  else {
    std::println("Unexpected value of IBUMP = {}", ibump);
    std::exit(EXIT_FAILURE);
  }

  int nbleft = round(xbleft * (mx - 1) / xlngth);  // 0-indexed
  int nbrite = round(xbrite * (mx - 1) / xlngth);
  std::println("Bump extends from {} at node {}", xbleft, nbleft + 1);
  std::println("               to {} at node {}", xbrite, nbrite + 1);

  neqn = 0;
  int ielemn = 0;

  for (int ip = 0; ip < n_points; ip++) {
    int ic, jc;
    if (_long) {
      ic = ip / my;
      jc = ip % my;
    } else {
      ic = ip % mx;
      jc = ip / mx;
    }
    int icnt = (ic + 1) % 2;
    int jcnt = (jc + 1) % 2;

    if ((icnt == 1 and jcnt == 1) and (ic != mx - 1) and (jc != my - 1)) {
      int ip1, ip2;
      if (_long) {
        ip1 = ip + my;
        ip2 = ip + my + my;

        node[ielemn][0] = ip;
        node[ielemn][1] = ip + 2;
        node[ielemn][2] = ip2 + 2;
        node[ielemn][3] = ip + 1;
        node[ielemn][4] = ip1 + 2;
        node[ielemn][5] = ip1 + 1;

        if (ibump == 0)
          isotri[ielemn] = 0;
        else if (ibump == 1)
          isotri[ielemn] = 0;
        else if (ibump == 2)
          isotri[ielemn] = (nbleft <= ic && ic < nbrite) ? 1 : 0;
        else
          isotri[ielemn] = 1;

        ielemn += 1;

        node[ielemn][0] = ip;
        node[ielemn][1] = ip2 + 2;
        node[ielemn][2] = ip2;
        node[ielemn][3] = ip1 + 1;
        node[ielemn][4] = ip2 + 1;
        node[ielemn][5] = ip1;

        if (ibump == 0)
          isotri[ielemn] = 0;
        else if (ibump == 1)
          isotri[ielemn] = (jc == 0 and nbleft <= ic  and ic < nbrite)? 1 : 0;
        else if (ibump == 2)
          isotri[ielemn] = (nbleft <= ic and ic < nbrite)? 1 : 0;
        else
          isotri[ielemn] = 1;

        ielemn += 1;
      } else {
        ip1 = ip + mx;
        ip2 = ip + mx + mx;

        node[ielemn][0] = ip;
        node[ielemn][1] = ip2;
        node[ielemn][2] = ip2 + 2;
        node[ielemn][3] = ip1;
        node[ielemn][4] = ip2 + 1;
        node[ielemn][5] = ip1 + 1;

        if (ibump == 0)
          isotri[ielemn] = 0;
        else if (ibump == 1)
          isotri[ielemn] = 0;
        else if (ibump == 2)
          isotri[ielemn] = (nbleft <= ic and ic < nbrite)? 1 : 0;
        else
          isotri[ielemn] = 1;

        ielemn += 1;

        node[ielemn][0] = ip;
        node[ielemn][1] = ip2 + 2;
        node[ielemn][2] = ip + 2;
        node[ielemn][3] = ip1 + 1;
        node[ielemn][4] = ip1 + 2;
        node[ielemn][5] = ip + 1;

        if (ibump == 0)
          isotri[ielemn] = 0;
        else if (ibump == 1)
          isotri[ielemn] = (jc == 0 and nbleft <= ic and ic < nbrite)? 1 : 0;
        else if (ibump == 2)
          isotri[ielemn] = (nbleft <= ic and ic < nbrite)? 1 : 0;
        else
          isotri[ielemn] = 1;

        ielemn += 1;
      }
    }

    if (ic == 0 and 0 < jc and jc < my - 1) {
      indx[ip][0] = -1;
      indx[ip][1] = -1;
    } else if (ic == mx - 1 and 0 < jc and jc < my - 1) {
      neqn += 1;
      indx[ip][0] = neqn;
      indx[ip][1] = 0;
    } else if (jc == 0 and ielemn > 0 and isotri[ielemn - 1] == 1) {
      indx[ip][0] = -2;
      indx[ip][1] = -2;
    } else if (ic == 0 or ic == mx - 1 or jc == 0 or jc == my - 1) {
      indx[ip][0] = 0;
      indx[ip][1] = 0;
    } else {
      neqn += 2;
      indx[ip][0] = neqn - 1;
      indx[ip][1] = neqn;
    }
    
    if (jcnt == 1 and icnt == 1) {
      neqn += 1;
      insc[ip] = neqn;
    } else
      insc[ip] = 0;
  }
  if (1 <= iwrite) {
    std::println();
    std::println("     I     INDX 1, INDX 2, INSC");
    std::println();
    for (int i = 0; i < n_points; i++) {
      std::println(
        "{}\t{}\t{}\t{}", 
        i+1, indx[i][0], indx[i][1], insc[i]
      );
    }
    std::println();
    std::println("Isoparametric triangles:");
    std::println();
    for (int i = 0; i < nelemn; i++){
      if (isotri[i] == 1)
        std::println("{}", i + 1);
    }
    std::println();
    std::println("   IT   NODE(IT,*)");
    std::println();
    for (int it = 0; it < nelemn; it++){
      std::print("{}\t", it+1);
      for (int i = 0; i < 6; i++) {
        std::print("{}\t", node[it][i]+1);
      }
      std::println();
    }
  }

  std::println();
  std::println("SETGRD: Number of unknowns = {}", neqn);

  if (maxeqn < neqn) {
    std::println("SETGRD - Too many unknowns!");
    std::println("The maximum allowed is MAXEQN = {}", maxeqn);
    std::println("This problem requires NEQN = {}", neqn);
    std::exit(1);
  }
  return std::tuple {_long, neqn};
}

//* --------------------------------------------------------------------
//*  SETLIN - determine unknown numbers along the profile line
//* --------------------------------------------------------------------
void setlin(
  std::vector<int> &iline, 
  std::vector<std::vector<int>> &indx, 
  int iwrite, bool _long, int mx, int my, 
  int n_points, int nx, int ny, 
  double xlngth, double xprof
) {
  int itemp, nodex0, ip;

  itemp = std::round(2.0 * (nx - 1) * xprof / xlngth);

  if (_long)
    nodex0 = itemp * (2 * ny - 1);
  else
    nodex0 = itemp;

  std::println();
  std::println("SETLIN:");
  std::println();
  std::println("  Profile generated at X = {}", xprof);
  std::println("  which is above node  = {}", nodex0 + 1);

  for (int i = 0; i < my; i++) {
    if (_long)
      ip = nodex0 + i;
    else
      ip = nodex0 + mx * i;
    iline[i] = indx[ip][0];
  }

  if (1 <= iwrite) {
    std::println();
    std::println("  Indices of unknowns along the profile line:");
    std::println();
    for (int i = 0; i < my; i += 5) {
      for (int j = i; j < std::min(i + 5, my); ++j) {
        std::print("{:5d}", iline[j]);
      }
      std::println("");
    }
  }
}

//* --------------------------------------------------------------------
//*  SETQUD - set midpoint quadrature rule information
//* --------------------------------------------------------------------
void setqud(
  std::vector<double> &area, std::vector<int> &isotri, int iwrite,
  int nelemn, int nnodes, std::vector<std::vector<int>> &node,
  int n_points, int nquad,
  std::vector<double> &xc, std::vector<std::vector<double>> &xm,
  std::vector<double> &yc, std::vector<std::vector<double>> &ym
) {
  for (int it = 0; it < nelemn; it++) {
    int ip1 = node[it][0];
    int ip2 = node[it][1];
    int ip3 = node[it][2];
    double x1 = xc[ip1];
    double x2 = xc[ip2];
    double x3 = xc[ip3];
    double y1 = yc[ip1];
    double y2 = yc[ip2];
    double y3 = yc[ip3];

    if (isotri[it] == 0) {
      xm[it][0] = 0.5 * (x1 + x2);
      xm[it][1] = 0.5 * (x2 + x3);
      xm[it][2] = 0.5 * (x3 + x1);
      ym[it][0] = 0.5 * (y1 + y2);
      ym[it][1] = 0.5 * (y2 + y3);
      ym[it][2] = 0.5 * (y3 + y1);
      area[it] = 0.5 * std::abs(
        (y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3)
      );
    } else {
      xm[it][0] = 0.5;
      ym[it][0] = 0.5;
      xm[it][1] = 1.0;
      ym[it][1] = 0.5;
      xm[it][2] = 0.5;
      ym[it][2] = 0.0;
      area[it] = 0.5;
    }
  }
  if (3 <= iwrite) {
    std::println();
    std::println("SETQUD: Element Areas and Quadrature points:");
    std::println();
    for (int i = 0; i < nelemn; i++) {
      std::println("{}\t{}", i + 1, area[i]);
      for (int j =0; j < nquad; j++) {
        std::println("{}\t{}\t{}\t{}", i + 1, j + 1, xm[i][j], ym[i][j]);
      }
    }
  }
}

//* --------------------------------------------------------------------
//*  SETXY - set the grid coordinates based on the parameter value
//* --------------------------------------------------------------------
void setxy(
  int iwrite, bool _long, int mx, int my, int n_points, 
  int nx, int ny, 
  std::vector<double> &xc, double xlngth, 
  std::vector<double> &yc, double ylngth, double ypert
) {
  for (int ip = 0; ip < n_points; ip++){
    int ic, jc;
    if (_long) {
      ic = ip / my;
      jc = ip % my;
    } else {
      ic = ip % mx;
      jc = ip / mx;
    }
    xc[ip] = ic * xlngth / (2 * nx - 2);

    double ybot = -ypert * (xc[ip] - 3.0) * (xc[ip] - 1.0);
    double ylo = std::max(0.0, ybot);

    yc[ip] = ((my - 1 - jc) * ylo + jc * ylngth) / (2 * ny - 2);
  }
  if (2 <= iwrite) {
    std::println();
    std::println("SETXY:");
    std::println();
    std::println("     I     XC     YC");
    std::println();
    for (int i = 0; i < n_points; i++) {
      std::println("{}\t{}\t{}", i+1, xc[i], yc[i]);
    }
  }
}

//* --------------------------------------------------------------------
//*  TIMESTAMP - std::println current date and time
//* --------------------------------------------------------------------
void timestamp()
{
  auto ahora = std::chrono::system_clock::now();
  //? Formato: Año-Mes-Dia Hora:Minuto:Segundo
  //? %F -> %Y-%m-%d, %T -> %H:%M:%S
  std::println("{0:%Y-%m-%d %H:%M:%S}", ahora);
}

//* --------------------------------------------------------------------
//*  TRANS - calculate the element transformation mapping
//* --------------------------------------------------------------------
std::tuple<double, double, double, double, double> trans(
  int it, int nelemn, int nnodes, 
  const std::vector<std::vector<int>> &node, 
  int n_points, 
  const std::vector<double> &xc, double xq, 
  const std::vector<double> &yc, double yq
) {
  int i1 = node[it][0];
  int i2 = node[it][1];
  int i3 = node[it][2];
  int i4 = node[it][3];
  int i5 = node[it][4];
  int i6 = node[it][5];

  double x1 = xc[i1];
  double y1 = yc[i1];
  double x2 = xc[i2];
  double y2 = yc[i2];
  double x3 = xc[i3];
  double y3 = yc[i3];
  double x4 = xc[i4];
  double y4 = yc[i4];
  double x5 = xc[i5];
  double y5 = yc[i5];
  double x6 = xc[i6];
  double y6 = yc[i6];

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

  double det = (
    (2.0 * a1 * b2 - 2.0 * a2 * b1) * xq * xq
    + (4.0 * a1 * c2 - 4.0 * a2 * c1) * xq * yq
    + (2.0 * b1 * c2 - 2.0 * b2 * c1) * yq * yq
    + (2.0 * a1 * e2 + b2 * d1 - b1 * d2 - 2.0 * a2 * e1) * xq
    + (2.0 * c2 * d1 + b1 * e2 - b2 * e1 - 2.0 * c1 * d2) * yq
    + d1 * e2
    - d2 * e1
  );

  constexpr double eps = 1e-30;
  if (std::abs(det) < eps)
    det = eps;

  double xix = dydeta / det;
  double xiy = -dxdeta / det;
  double etax = -dydxi / det;
  double etay = dxdxi / det;

  return std::tuple {det, etax, etay, xix, xiy};
}

//* --------------------------------------------------------------------
//*  UBDRY - parabolic inflow boundary condition
//* --------------------------------------------------------------------
inline double ubdry(int iuk, double yy) {
  if (iuk == 1)
    return (-2.0 * yy + 6.0) * yy / 9.0;
  else
    return 0.0;
}

//* --------------------------------------------------------------------
//*  UBUMP - sensitivity dU/dA on the bump
//* --------------------------------------------------------------------
double ubump(
  const std::vector<double> &g,
  const std::vector<std::vector<int>> &indx,
  int ip, int iqq,
  const std::vector<int> &isotri,
  int it, int iukk, int nelemn,
  int neqn, int nnodes,
  const std::vector<std::vector<int>> &node, int n_points,
  const std::vector<double> &xc, const std::vector<double> &yc
) {
  double det = 1.0;
  double etax = 0.0;
  double etay = 0.0;
  double xix = 1.0;
  double xiy = 1.0;
  double xq = 0.0;
  double yq = 0.0;

  if (isotri[it] == 0) {
    xq = xc[ip];
    yq = yc[ip];
  } else {
    if (iqq == 0) {
      xq = 0.0;
      yq = 0.0;
    } else if (iqq == 1) {
      xq = 1.0;
      yq = 1.0;
    } else if (iqq == 2) {
      xq = 1.0;
      yq = 0.0;
    } else if (iqq == 3) {
      xq = 0.5;
      yq = 0.5;
    } else if (iqq == 4) {
      xq = 1.0;
      yq = 0.5;
    } else if (iqq == 5) {
      xq = 0.5;
      yq = 0.0;
    }
    
    std::tie(det, etax, etay, xix, xiy) = trans(
      it, nelemn, nnodes, node, n_points, xc, xq, yc, yq
    );
  }

  auto [un, unx, uny] = _ubump_uval(
    g, indx, isotri, it, nelemn, 
    neqn, nnodes, node, n_points,
    xc, xix, xiy, xq, yc, yq, 
    det, etax, etay
  );

  if (iukk == 1)
    return -uny[0] * (xc[ip] - 1.0) * (xc[ip] - 3.0);
  else if (iukk == 2)
    return -uny[1] * (xc[ip] - 1.0) * (xc[ip] - 3.0);
  else {
    std::println("UBUMP called for iukk = {}", iukk);
    std::exit(1);
  }
}
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> _ubump_uval(
  const std::vector<double> &g,
  const std::vector<std::vector<int>> &indx,
  const std::vector<int> &isotri,
  int it, int nelemn, int neqn, int nnodes,
  const std::vector<std::vector<int>> &node,
  int n_points, const std::vector<double> &xc,
  double xix, double xiy, double xq,
  const std::vector<double> &yc,
  double yq, double det, double etax, double etay
) {
  std::vector<double> un(2, 0.0);
  std::vector<double> unx(2, 0.0);
  std::vector<double> uny(2, 0.0);
  
  double bb, bx, by, ubc;

  for (int iq = 0; iq < nnodes; iq++) {
    if (isotri[it] == 1)
      std::tie(bb, bx, by) = refqbf(xq, yq, iq, etax, etay, xix, xiy);
    else
      std::tie(bb, bx, by) = qbf(xq, yq, it, iq, nelemn, nnodes, node, n_points, xc, yc);
    int ip_local = node[it][iq];

    for (int iuk = 0; iuk < 2; iuk++) {
      int iun = indx[ip_local][iuk];
      if (0 < iun) {
        un[iuk] += bb * g[iun - 1];
        unx[iuk] += bx * g[iun - 1];
        uny[iuk] += by * g[iun - 1];
      } else if (iun == -1) {
        ubc = ubdry(iuk + 1, yc[ip_local]);
        un[iuk] += bb * ubc;
        unx[iuk] += bx * ubc;
        uny[iuk] += by * ubc;
      }
    }
  }

  return std::tuple {un, unx, uny};
}

//* --------------------------------------------------------------------
//*  UVAL - evaluate velocities at a given quadrature point
//* --------------------------------------------------------------------
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> uval(
  double etax,
  double etay,
  const std::vector<double> &g,
  const std::vector<std::vector<int>> &indx,
  const std::vector<int> &isotri,
  int it,
  int nelemn,
  int neqn,
  int nnodes,
  const std::vector<std::vector<int>> &node,
  int n_points,
  const std::vector<double> &xc,
  double xix,
  double xiy,
  double xq,
  const std::vector<double> &yc,
  double yq
) {
  std::vector<double> un(2, 0.0);
  std::vector<double> unx(2, 0.0);
  std::vector<double> uny(2, 0.0);
  double bb, bx, by, ubc;

  for (int iq = 0; iq < nnodes; iq++) {
    if (isotri[it] == 1)
      std::tie(bb, bx, by) = refqbf(xq, yq, iq, etax, etay, xix, xiy);
    else
      std::tie(bb, bx, by) = qbf(xq, yq, it, iq, nelemn, nnodes, node, n_points, xc, yc);
    
    int ip_local = node[it][iq];

    for (int iuk = 0; iuk < 2; iuk++) {
      int iun = indx[ip_local][iuk];
      if (0 < iun) {
        un[iuk] += bb * g[iun - 1];
        unx[iuk] += bx * g[iun - 1];
        uny[iuk] += by * g[iun - 1];
      } else if (iun == -1) {
        ubc = ubdry(iuk + 1, yc[ip_local]);
        un[iuk] += bb * ubc;
        unx[iuk] += bx * ubc;
        uny[iuk] += by * ubc;
      }
    }
  }

  return std::tuple {un, unx, uny};
}

//* --------------------------------------------------------------------
//*  UV_WRITE - write a velocity file
//* --------------------------------------------------------------------
void uv_write(
  const std::vector<double>& f, 
  const std::vector<std::vector<int>>& indx, 
  std::ofstream& uv_file_obj, 
  int neqn, int n_points, 
  const std::vector<double>& yc
) {
  
  uv_file_obj << std::scientific << std::setprecision(6);
  double u, v;

  for (int ip = 0; ip < n_points; ++ip) {
    // Lógica para 'u'
    int k_u = indx[ip][0];

    if (k_u < 0) u = ubdry(1, yc[ip]);
    else if (k_u == 0) u = 0.0;
    else u = f[k_u - 1];

    // Lógica para 'v'
    int k_v = indx[ip][1];

    if (k_v < 0) v = ubdry(2, yc[ip]);
    else if (k_v == 0) v = 0.0;
    else v = f[k_v - 1];

    uv_file_obj << "  " << std::setw(14) << u 
                << "  " << std::setw(14) << v << "\n";
  }
}

//* --------------------------------------------------------------------
//*  XY_WRITE - write node coordinate data
//* --------------------------------------------------------------------
void xy_write(
  std::ofstream& xy_file_obj, int n_points, 
  std::vector<double>& xc, std::vector<double>& yc
) {
  xy_file_obj << std::scientific << std::setprecision(6);
  
  for (int ip = 0; ip < n_points; ++ip) {
    // El formato :14.6e en Python se traduce así:
    xy_file_obj << "  " << std::setw(14) << xc[ip] 
                << "  " << std::setw(14) << yc[ip] << "\n";
  }
}
