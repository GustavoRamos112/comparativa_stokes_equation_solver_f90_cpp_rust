
use chrono::Local;

use crate::program::bump_struct;

pub fn phi_idx(it: usize, iquad: usize, iq: usize, k: usize, nquad: usize, nnodes: usize) -> usize {
  return ((it * nquad + iquad) * nnodes + iq) * 3 + k
}
pub fn psi_idx(it: usize, iquad: usize, iq: usize, nquad: usize, nnodes: usize) -> usize {
  return (it * nquad + iquad) * nnodes + iq
}

//* --------------------------------------------------------------------
//*  BSP - linear basis function for pressure
//* --------------------------------------------------------------------
pub fn bsp(
  it: usize, iq: usize, id: usize,
  node: &Vec<usize>,
  xc: &Vec<f64>, xq: f64,
  yc: &Vec<f64>, yq: f64,
  nnodes: &usize
) -> f64 {
  let l1 = iq;
  let l2 = i4_wrap(iq + 1, 0, 2);
  let l3 = i4_wrap(iq + 2, 0, 2);

  let g1: usize = node[it * nnodes + l1];
  let g2: usize = node[it * nnodes + l2];
  let g3: usize = node[it * nnodes + l3];

  let d: f64 =
    (xc[g2] - xc[g1]) * (yc[g3] - yc[g1]) 
    - (xc[g3] - xc[g1]) * (yc[g2] - yc[g1]);

  if id == 0 {
    return 1.0 + ((yc[g2] - yc[g3]) * (xq - xc[g1]) + (xc[g3] - xc[g2]) * (yq - yc[g1])) / d;
  }
  else if id == 1 {
    return (yc[g2] - yc[g3]) / d;
  } else if id == 2 {
    return (xc[g3] - xc[g2]) / d;
  } else {
    println!();
    println!("BSP - Fatal error!");
    println!("  Illegal local index value for linear basis.");
    println!("  Legal values are 1, 2 or 3.");
    println!("  The input value was ID = {}", id + 1);
    std::process::exit(1);
  }
}

////* --------------------------------------------------------------------
////*  DAXPY - constant times a vector plus a vector
////* --------------------------------------------------------------------
//void daxpy_v(
//  int n, double da, 
//  std::vector<double>& dx, int incx, 
//  std::vector<double>& dy, int incy
//) {
//  if (n <= 0 or da == 0.0) return;
//  int ix, iy, m;
//
//  if (incx != 1 or incy != 1) {
//    if (incx >= 0) ix = 0;
//    else ix = (-n + 1) * incx;
//    
//    if (incy >= 0) iy = 0;
//    else iy = (-n + 1) * incy;
//
//    for (int i = 0; i < n; i++) {
//      dy[iy] += da * dx[ix];
//      ix += incx;
//      iy += incy;
//    }
//  } else {
//    m = n % 4;
//    for (int i = 0; i < m; i++) {
//      dy[i] += da * dx[i];
//    }
//    for (int i = m; i < n; i += 4) {
//      dy[i] += da * dx[i];
//      dy[i + 1] += da * dx[i + 1];
//      dy[i + 2] += da * dx[i + 2];
//      dy[i + 3] += da * dx[i + 3];
//    }
//  }
//}
//void daxpy_m(
//  int n, double da, 
//  std::vector<double> &abd,
//  int stride, int col_x, int start_row_x, int incx,
//  int col_y, int start_row_y, int incy
//) {
//  if (n <= 0 || da == 0.0) return;
//
//  for (int i = 0; i < n; ++i) {
//    abd[(start_row_y + (i * incy)) * stride + col_y] += da * abd[(start_row_x + (i * incx)) * stride + col_x];
//  }
//}
//void daxpy(
//  int n, double da, 
//  const std::vector<double>& abd_x, int stride, int col_x, int start_row_x, int incx,
//  std::vector<double>& b_y, int start_idx_y, int incy
//) {
//    
//    if (n <= 0 || da == 0.0) return;
//
//    for (int i = 0; i < n; ++i) {
//        b_y[start_idx_y + (i * incy)] += da * abd_x[(start_row_x + (i * incx)) * stride + col_x];
//    }
//}
//
////* --------------------------------------------------------------------
////*  DDOT - dot product of two vectors
////* --------------------------------------------------------------------
//double ddot(
//  int n, 
//  const std::vector<double> &abd,
//  int stride, int col, int row_start, int incx, 
//  const std::vector<double> &b, int start_b, int incy
//) {
//    
//  if (n <= 0) return 0.0;
//
//  double dtemp = 0.0;
//
//  if (incx != 1 || incy != 1) {
//    int ix = row_start;
//    int iy = start_b;
//    for (int i = 0; i < n; ++i) {
//      dtemp += abd[ix * stride + col] * b[iy];
//      ix += incx;
//      iy += incy;
//    }
//  } 
//  else {
//    int m = n % 5;
//    for (int i = 0; i < m; ++i) {
//      dtemp += abd[(row_start + i) * stride + col] * b[start_b + i];
//    }
//    for (int i = m; i < n; i += 5) {
//      dtemp += (
//        abd[(row_start + i) * stride + col] * b[start_b + i]
//        + abd[(row_start + i + 1) * stride + col] * b[start_b + i + 1]
//        + abd[(row_start + i + 2) * stride + col] * b[start_b + i + 2]
//        + abd[(row_start + i + 3) * stride + col] * b[start_b + i + 3]
//        + abd[(row_start + i + 4) * stride + col] * b[start_b + i + 4]
//      );
//    }
//  }
//
//  return dtemp;
//}
//
////* --------------------------------------------------------------------
////*  DGBFA - factor a real band matrix by elimination
////* --------------------------------------------------------------------
//std::tuple<int, std::vector<int>> dgbfa(
//  std::vector<double> &abd, int lda, 
//  int n, int ml, int mu, 
//  std::vector<int> &ipvt
//) {
//  int m = ml + mu + 1;
//  int info = 0;
//
//  int j0 = mu + 1;
//  int j1 = min(n, m) - 1;
//  int i0;
//
//  for (int jz = j0; jz < j1; jz++) {
//    i0 = m - jz;
//    for (int i = i0; i < ml; i++) {
//      abd[i * n + jz] = 0.0;
//    }
//  }
//
//  int jz = j1;
//  int ju = 0, lm, l, mm;
//
//  double t;
//  for (int k = 0; k < n - 1; k++) {
//    jz += 1;
//    if (jz < n)
//      for (int i = 0; i < ml; i++) {
//        abd[i * n + jz] = 0.0;
//      }
//
//    lm = min(ml, n - k - 1);
//    l = idamax_m(lm + 1, abd, n, k, m - 1, 1) + (m - 1);
//    ipvt[k] = l + k - (m - 1);
//
//    if (abd[l * n + k] == 0.0)
//      info = k;
//    else {
//      if (l != m - 1) {
//        t = abd[l * n + k];
//        abd[l * n + k] = abd[(m - 1) * n + k];
//        abd[(m - 1) * n + k] = t;
//      }
//
//      t = -1.0 / abd[(m - 1) * n + k];
//      dscal_m(lm, t, abd, n, k, m, 1);
//
//      ju = min(max(ju, mu + ipvt[k]), n - 1);
//      mm = m - 1;
//
//      for (int j = k + 1; j < ju + 1; j++) {
//        l -= 1;
//        mm -= 1;
//        t = abd[l * n + j];
//        if (l != mm) {
//          abd[l * n + j] = abd[mm * n + j];
//          abd[mm * n + j] = t;
//        }
//        daxpy_m(lm, t, abd, n, k, m, 1, j, mm + 1, 1);
//      }
//    }
//  }
//
//  ipvt[n - 1] = n - 1;
//
//  if (abd[(m - 1) * n + (n - 1)] == 0.0)
//    info = n - 1;
//
//  return std::tuple {info, ipvt};
//}
//
////* --------------------------------------------------------------------
////*  DGBSL - solve a real banded system factored by DGBFA
////* --------------------------------------------------------------------
//void dgbsl(
//  std::vector<double> &abd, int lda,
//  int n, int ml, int mu,
//  std::vector<int> &ipvt, 
//  std::vector<double> &b, int job
//) {
//  int m = mu + ml + 1;
//  double t;
//  int l, lm, la, lb;
//
//  if (job == 0) {
//    if (0 < ml) {
//      for (int k = 0; k < n - 1; k++) {
//        lm = min(ml, n - k - 1);
//        l = ipvt[k];
//        t = b[l];
//        if (l != k) {
//          b[l] = b[k];
//          b[k] = t;
//        }
//        daxpy(lm, t, abd, n, k, m, 1, b, k + 1, 1);
//      }
//    }
//
//    for (int k = n - 1; k >= 0; k--) {
//      b[k] /= abd[(m - 1) * n + k];
//      lm = min(k, m - 1);
//      la = m - 1 - lm;
//      lb = k - lm;
//      t = -b[k];
//      daxpy(lm, t, abd, n, k, la, 1, b, lb, 1);
//    }
//  } else {
//    for (int k = 0; k < n; k++) {
//      lm = min(k, m - 1);
//      la = m - 1 - lm;
//      lb = k - lm;
//      t = ddot(lm, abd, n, k, la, 1, b, lb, 1);
//      b[k] = (b[k] - t) / abd[(m - 1) * n + k];
//    }
//
//    if (0 < ml) {
//      for (int k = n - 2; k >= 0; k--) {
//        lm = min(ml, n - k - 1);
//        b[k] += ddot(lm, abd, n, k, m, 1, b, k + 1, 1);
//        l = ipvt[k];
//        if (l != k) {
//          t = b[l];
//          b[l] = b[k];
//          b[k] = t;
//        }
//      }
//    }
//  }
//}
//
////* --------------------------------------------------------------------
////*  DSCAL - scale a vector by a constant
////* --------------------------------------------------------------------
//void dscal_v(
//  int n, double sa, 
//  std::vector<double>& x, int incx
//) {
//  int ix, m;
//  if (n <= 0) return;
//
//  if (incx == 1) {
//    m = n % 5;
//    for (int i = 0; i < m; i++) {
//      x[i] *= sa;
//    }
//    for (int i = m; i < n; i += 5) {
//      x[i] *= sa;
//      x[i + 1] *= sa;
//      x[i + 2] *= sa;
//      x[i + 3] *= sa;
//      x[i + 4] *= sa;
//    }
//  } else {
//    if (incx >= 0) ix = 0;
//    else ix = (-n + 1) * incx;
//    for (int i = 0; i < n; i++) {
//      x[ix] *= sa;
//      ix += incx;
//    }
//  }
//}
//void dscal_m(int n, double sa, std::vector<double>& abd, int stride, int col, int start_row, int incx) {
//  if (n <= 0) return;
//
//  for (int i = 0; i < n; ++i) {
//    abd[(start_row + (i * incx)) * stride + col] *= sa;
//  }
//}
////* --------------------------------------------------------------------
////*  FILE_NAME_INC - increment a partially numeric filename
////* --------------------------------------------------------------------
//std::string file_name_inc(std::string file_name) {
//  if (file_name.empty()) {
//    println!("\nFILE_NAME_INC - Fatal error!");
//    println!("The input string is empty.");
//    std::process::exit(1);
//  }
//
//  int change = 0;
//  for (int i = file_name.length() - 1; i >= 0; --i) {
//    if (file_name[i] >= '0' && file_name[i] <= '9') {
//      change++;
//      int digit = file_name[i] - '0';
//      digit++;
//      
//      if (digit == 10) {
//        digit = 0;
//        file_name[i] = '0';
//      } else {
//        file_name[i] = static_cast<char>(digit + '0');
//        return file_name;
//      }
//  }
//  }
//
//  return (change == 0) ? "" : file_name;
//}
//
////* --------------------------------------------------------------------
////*  GETG - extract values of a quantity along the profile line
////* --------------------------------------------------------------------
//std::array<double, my> getg(
//  const std::vector<double> &f, const std::array<int, my> &iline, int bump.neqn
//) {
//  std::array<double, my> u{};
//  int j;
//  for (int i = 0; i < my; i++) {
//    j = iline[i];
//    if (j < 0)
//      u[i] = 0.0;
//    else if (j == 0) u[i] = 0.0;
//    else u[i] = f[j - 1];
//  }
//  return u;
//}
//
////* --------------------------------------------------------------------
////*  GRAM - compute the Gram matrix && R vector
////* --------------------------------------------------------------------
//void gram(
//  std::array<double, my * my> &gr,
//  std::array<int, my> &iline,
//  std::array<int, bump.n_points * 2> &indx,
//  int bump.iwrite,
//  std::array<int, bump.nelemn * bump.nnodes> &node,
//  std::array<double, my> &r,
//  std::array<double, my> &uprof,
//  std::array<double, bump.n_points> &xc, double xprof,
//  std::array<double, bump.n_points> &yc
//) {
//  constexpr std::array<double, 3> wt = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
//  constexpr std::array<double, 3> yq_gauss = {-0.7745966692, 0.0, 0.7745966692};
//
//  for (int it = 0; it < my; it++) {
//    r[it] = 0.0;
//  }
//
//  for (int it = 0; it < my; it++) {
//    for (int j = 0; j < my; j++) {
//      gr[it * my + j] = 0.0;
//    }
//  }
//  for (int it = 0; it < bump.nelemn; it++) {
//    int k = bump.node[it * bump.nnodes];
//    int kk = bump.node[it * bump.nnodes + 1];
//
//    if ((std::abs(bump.xc[k] - xprof) > 1.0e-4) or (std::abs(bump.xc[kk] - xprof) > 1.0e-4))
//      continue;
//
//    for (int iquad = 0; iquad < 3; iquad++) {
//      double bma2 = (yc[kk] - bump.yc[k]) / 2.0;
//      double ar = bma2 * wt[iquad];
//      double x = xprof;
//      double y = bump.yc[k] + bma2 * (yq_gauss[iquad] + 1.0);
//
//      double uiqdpt = 0.0;
//      for (int iq = 0; iq < bump.nnodes; iq++) {
//        if ((iq == 0) or (iq == 1) or (iq == 3)) {
//          auto [bb, bx, by] = qbf(x, y, it, iq, bump.node, xc, yc);
//          int ip = bump.node[it * bump.nnodes + iq];
//          int iun = bump.indx[ip * 2];
//          if (0 < iun) {
//            int ii = igetl(iun, iline);
//            uiqdpt += bb * uprof[ii - 1];
//          } else if (iun == -1) {
//            double ubc = ubdry(1, bump.yc[ip]);
//            uiqdpt += bb * ubc;
//          }
//        }
//      }
//
//      for (int iq = 0; iq < bump.nnodes; iq++) {
//        if ((iq == 0) or (iq == 1) or (iq == 3)) {
//          int ip = bump.node[it * bump.nnodes + iq];
//          auto [bb, bx, by] = qbf(x, y, it, iq, bump.node, xc, yc);
//          int i_val = bump.indx[ip * 2];
//          if (0 < i_val) {
//            int ii = igetl(i_val, iline);
//            r[ii - 1] += bb * uiqdpt * ar;
//            for (int iqq = 0; iqq < bump.nnodes; iqq++) {
//              if (iqq == 0 or iqq == 1 or iqq == 3) {
//                int ipp = bump.node[it * bump.nnodes + iqq];
//                auto [bbb, bbx_, bby_] = qbf(
//                  x, y, it, iqq, bump.node, xc, yc
//                );
//                int j_val = bump.indx[ipp * 2];
//                if (j_val != 0) {
//                  int jj = igetl(j_val, iline);
//                  gr[(ii - 1) * my + (jj - 1)] += bb * bbb * ar;
//                }
//              }
//            }
//          }
//        }
//      }
//    }
//  }
//  if (3 <= bump.iwrite) {
//    println!();
//    println!("Gram matrix:");
//    println!();
//    for (int i = 0; i < my; i++) {
//      for (int j = 0; j < my; j++) {
//        println!("{} {} {}", i + 1, j + 1, gr[i * my + j]);
//      }
//    }
//    println!();
//    println!("R vector:");
//    println!();
//    for (int i = 0; i < my; i++) {
//      println!("{}", r[i]);
//    }
//  }
//}

//* --------------------------------------------------------------------
//*  I4_MODP - nonnegative remainder of integer division
//* --------------------------------------------------------------------
pub fn i4_modp(i: i32, j: i32) -> usize {
  if j == 0 {
    println!();
    println!("I4_MODP - Fatal error!");
    println!("  Illegal divisor J = {}", j);
    std::process::exit(1);
  }
  let mut value: i32 = (i as i32) % j;
  if value < 0 {
    value += j.abs();
  }
  return value as usize;
}

//* --------------------------------------------------------------------
//*  I4_WRAP - force an integer between limits by wrapping
//* --------------------------------------------------------------------
pub fn i4_wrap(ival: usize, ilo: usize, ihi: usize) -> usize {
  let jlo: usize = ilo.min(ihi);
  let jhi: usize = ilo.max(ihi);
  let wide: i32 = (jhi - jlo + 1) as i32;
  if wide == 1 {
    return jlo;
  }
  return jlo + i4_modp((ival - jlo) as i32, wide);
}

////* --------------------------------------------------------------------
////*  IDAMAX - index of the vector element of maximum absolute value
////* --------------------------------------------------------------------
//int idamax_v(int n, std::vector<double> &dx, int incx) {
//  double dmax;
//  int ix;
//
//  if (n < 1 or incx <= 0)
//    return -1;
//
//  int idamax_val = 0;
//
//  if (n == 1)
//    return 0;
//
//  if (incx == 1) {
//    dmax = std::abs(dx[0]);
//    for (int i = 1; i < n; i++) {
//      if (dmax < std::abs(dx[i])) {
//        idamax_val = i;
//        dmax = std::abs(dx[i]);
//      }
//    }
//  } else {
//    ix = 0;
//    dmax = std::abs(dx[0]);
//    ix += incx;
//    for (int i = 1; i < n; i++) {
//      if (dmax < std::abs(dx[ix])) {
//        idamax_val = i;
//        dmax = std::abs(dx[ix]);
//      }
//      ix += incx;
//    }
//  }
//
//  return idamax_val;
//}
//
//int idamax_m(int n, const std::vector<double>& abd, int stride, int col, int start_row, int incx) {
//  if (n < 1 || incx <= 0) return -1;
//
//  int idamax_val = 0;
//  double dmax = std::abs(abd[start_row * stride + col]);
//  
//  for (int i = 1; i < n; i++) {
//    double current_val = std::abs(abd[(start_row + (i * incx)) * stride + col]);
//    if (current_val > dmax) {
//      idamax_val = i;
//      dmax = current_val;
//    }
//  }
//  return idamax_val;
//}
//
////* --------------------------------------------------------------------
////*  IGETL - get the local unknown number along the profile line
////* --------------------------------------------------------------------
//inline int igetl(int i, const std::array<int, my> &iline) {
//  for (int j = 0; j < my; j++) {
//    if (iline[j] == i) return j + 1;
//  }
//  return -1;
//}
//
////* --------------------------------------------------------------------
////*  LINSYS - solve the linearized Navier Stokes equation
////* --------------------------------------------------------------------
//void linsys(
//  std::vector<double> &a,
//  std::array<double, bump.nelemn> &area,
//  std::vector<double> &f, std::vector<double> &g,
//  std::array<int, bump.n_points * 2> &indx,
//  std::array<int, bump.n_points> &insc, std::array<int, bump.nelemn> &isotri,
//  int itype, int maxrow, int nband,
//  int bump.neqn, int nlband,
//  std::array<int, bump.nelemn * bump.nnodes> &node,
//  int nrow,
//  std::array<double, bump.nelemn * bump.nquad * bump.nnodes * 3> &phi,
//  std::array<double, bump.nelemn * bump.nquad * bump.nnodes> &psi, double reynld,
//  std::array<double, bump.n_points> &xc, std::array<double, bump.nelemn * bump.nquad> &xm,
//  std::array<double, bump.n_points> &yc, std::array<double, bump.nelemn * bump.nquad> &ym
//) {
//  int ioff = (
//    nlband + nlband + 1
//  );
//  double visc = 1.0 / reynld;
//  int ip, ihor, iver, iprs;
//  int ipp;
//  double bbb, bbx, bby;
//  int ju, jv, jp, iuse, info;
//
//  for (int i = 0; i < bump.neqn; i++) {
//    f[i] = 0.0;
//  }
//
//  for (int i = 0; i < nrow; i++) {
//    for (int j = 0; j < bump.neqn; j++) {
//      a[i * bump.neqn + j] = 0.0;
//    }
//  }
//
//  std::array<double, 2> un{};
//  std::array<double, 2> unx{};
//  std::array<double, 2> uny{};
//  double det, etax, etay, xix, xiy;
//  double ar, yq, xq, bbl, bbbl, ubc, aij;
//  double bb, bx, by;
//
//  for (int it = 0; it < bump.nelemn; it++) {
//    ar = bump.area[it] / 3.0;
//
//    for (int iquad = 0; iquad < bump.nquad; iquad++) {
//      yq = bump.ym[it * bump.nquad + iquad];
//      xq = bump.xm[it * bump.nquad + iquad];
//
//      det = 0.0;
//      etax = 0.0;
//      etay = 0.0;
//      xix = 0.0;
//      xiy = 0.0;
//      ubc = 0.0;
//
//      if (bump.isotri[it] == 1) {
//        std::tie(det, etax, etay, xix, xiy) = trans(
//          it, bump.node, xc, xq, yc, yq
//        );
//        ar = det * bump.area[it] / 3.0;
//      }
//
//      std::tie(un, unx, uny) = uval(
//        etax, etay, g, bump.indx, bump.isotri, it,
//        bump.neqn, bump.node, xc, xix,
//        xiy, xq, yc, yq
//      );
//
//      for (int iq = 0; iq < bump.nnodes; iq++) {
//        ip = bump.node[it * bump.nnodes + iq];
//        bb = bump.phi[phi_idx(it, iquad, iq, 0)];
//        bx = bump.phi[phi_idx(it, iquad, iq, 1)];
//        by = bump.phi[phi_idx(it, iquad, iq, 2)];
//        bbl = psi[psi_idx(it, iquad, iq)];
//        ihor = bump.indx[ip * 2];
//        iver = bump.indx[ip * 2 + 1];
//        iprs = bump.insc[ip];
//
//        if (0 < ihor)
//          f[ihor - 1] += ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
//
//        if (0 < iver)
//          f[iver - 1] += ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
//
//        for (int iqq = 0; iqq < bump.nnodes; iqq++) {
//          ipp = bump.node[it * bump.nnodes + iqq];
//          bbb = bump.phi[phi_idx(it, iquad, iqq, 0)];
//          bbx = bump.phi[phi_idx(it, iquad, iqq, 1)];
//          bby = bump.phi[phi_idx(it, iquad, iqq, 2)];
//          bbbl = psi[psi_idx(it, iquad, iqq)];
//          ju = bump.indx[ipp * 2];
//          jv = bump.indx[ipp * 2 + 1];
//          jp = bump.insc[ipp];
//
//          if (0 < ju) {
//            if (0 < ihor) {
//              iuse = ihor - ju + ioff;
//              a[(iuse - 1) * bump.neqn + (ju - 1)] += ar * (
//                visc * (by * bby + bx * bbx)
//                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
//              );
//            }
//
//            if (0 < iver) {
//              iuse = iver - ju + ioff;
//              a[(iuse - 1) * bump.neqn + (ju - 1)] += ar * bb * bbb * unx[1];
//            }
//
//            if (0 < iprs) {
//              iuse = iprs - ju + ioff;
//              a[(iuse - 1) * bump.neqn + (ju - 1)] += ar * bbx * bbl;
//            }
//          } else if (ju == itype) {
//            if (ju == -1)
//              ubc = ubdry(1, bump.yc[ipp]);
//            else if (ju == -2) {
//              ubc = ubump(
//                g,
//                bump.indx,
//                ipp, iqq,
//                bump.isotri,
//                it, 1,
//                bump.neqn,
//                bump.node,
//                xc, yc
//              );
//            }
//
//            if (0 < ihor) {
//              aij = ar * (
//                visc * (by * bby + bx * bbx)
//                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
//              );
//              f[ihor - 1] -= ubc * aij;
//            }
//
//            if (0 < iver) {
//              aij = ar * bb * bbb * unx[1];
//              f[iver - 1] -= ubc * aij;
//            }
//
//            if (0 < iprs) {
//              aij = ar * bbx * bbl;
//              f[iprs - 1] -= ubc * aij;
//            }
//          }
//
//          if (0 < jv) {
//            if (0 < ihor) {
//              iuse = ihor - jv + ioff;
//              a[(iuse - 1) * bump.neqn + (jv - 1)] += ar * bb * bbb * uny[0];
//            }
//
//            if (0 < iver) {
//              iuse = iver - jv + ioff;
//              a[(iuse - 1) * bump.neqn + (jv - 1)] += ar * (
//                visc * (by * bby + bx * bbx)
//                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
//              );
//            }
//
//            if (0 < iprs) {
//              iuse = iprs - jv + ioff;
//              a[(iuse - 1) * bump.neqn + (jv - 1)] += ar * bby * bbl;
//            }
//
//          } else if (jv == itype) {
//            if (jv == -1)
//              ubc = ubdry(2, bump.yc[ipp]);
//            else if (jv == -2) {
//              ubc = ubump(
//                g,
//                bump.indx,
//                ipp,
//                iqq,
//                bump.isotri,
//                it,
//                2,
//                bump.neqn,
//                bump.node,
//                xc,
//                yc
//              );
//            }
//
//            if (0 < ihor) {
//              aij = ar * bb * bbb * uny[0];
//              f[ihor - 1] -= ubc * aij;
//            }
//
//            if (0 < iver) {
//              aij = ar * (
//                visc * (by * bby + bx * bbx)
//                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
//              );
//              f[iver - 1] -= ubc * aij;
//            }
//
//            if (0 < iprs) {
//              aij = ar * bby * bbl;
//              f[iprs - 1] -= ubc * aij;
//            }
//          }
//          if (0 < jp) {
//            if (0 < ihor) {
//              iuse = ihor - jp + ioff;
//              a[(iuse - 1) * bump.neqn + (jp - 1)] -= ar * bx * bbbl;
//            }
//
//            if (0 < iver) {
//              iuse = iver - jp + ioff;
//              a[(iuse - 1) * bump.neqn + (jp - 1)] -= ar * by * bbbl;
//            }
//          }
//        }
//      }
//    }
//  }
//
//  f[bump.neqn - 1] = 0.0;
//  int i, j;
//  for (int j_1based = bump.neqn - nlband; j_1based < bump.neqn; j_1based++) {
//    j = j_1based - 1;
//    i = bump.neqn - j_1based + ioff;
//    a[(i - 1) * bump.neqn + j] = 0.0;
//  }
//  a[(ioff - 1) * bump.neqn + (bump.neqn - 1)] = 1.0;
//
//  std::vector<int> ipvt(bump.neqn, 0);
//  std::tie(info, ipvt) = dgbfa(a, maxrow, bump.neqn, nlband, nlband, ipvt);
//
//  if (info != 0) {
//    println!();
//    println!("LINSYS - fatal error!");
//    println!("DGBFA returns INFO = {}", info+1);
//    std::process::exit(1);
//  }
//
//  int job = 0;
//  dgbsl(a, maxrow, bump.neqn, nlband, nlband, ipvt, f, job);
//}
//
//
////* --------------------------------------------------------------------
////*  NSTOKE - solve Navier Stokes using Newton iteration
////* --------------------------------------------------------------------
//int nstoke(
//  std::vector<double> &a,
//  std::array<double, bump.nelemn> &area,
//  std::vector<double> &f, std::vector<double> &g,
//  std::array<int, bump.n_points * 2> &indx,
//  std::array<int, bump.n_points> &insc, std::array<int, bump.nelemn> &isotri,
//  int maxrow, int nband,
//  int bump.neqn, int nlband,
//  std::array<int, bump.nelemn * bump.nnodes> &node,
//  int nrow, int numnew,
//  std::array<double, bump.nelemn * bump.nquad * bump.nnodes * 3> &phi,
//  std::array<double, bump.nelemn * bump.nquad * bump.nnodes> &psi,
//  double reynld, double tolnew,
//  std::array<double, bump.n_points> &xc, std::array<double, bump.nelemn * bump.nquad> &xm,
//  std::array<double, bump.n_points> &yc, std::array<double, bump.nelemn * bump.nquad> &ym
//) {
//  int itype, imax;
//  double diff;
//
//  for (int iter = 0; iter < maxnew; iter++) {
//    numnew++;
//
//    itype = -1;
//    linsys(
//      a, area,
//      f, g,
//      bump.indx,
//      bump.insc,
//      bump.isotri,
//      itype, maxrow, nband,
//      bump.neqn, nlband,
//      bump.node,
//      nrow,
//      bump.phi,
//      psi, reynld,
//      xc, xm,
//      yc, ym
//    );
//
//    for (int i = 0; i < bump.neqn; i++) {
//      g[i] -= f[i];
//    }
//
//    imax = idamax_v(bump.neqn, g, 1);
//    
//    diff = std::abs(g[imax]);
//    println!("NSTOKE: Iteration {}, MaxNorm(diff) = {}", iter+1, diff);
//
//    for (int i = 0; i < bump.neqn; i++) {
//      g[i] = f[i];
//    }
//
//    if (diff <= tolnew) {
//      println!("NSTOKE converged.");
//      return numnew;
//    }
//
//    if (iter == maxnew - 1) {
//      println!("NSTOKE failed!");
//      std::process::exit(1);
//    }
//  }
//
//  return numnew;
//}

//* --------------------------------------------------------------------
//*  QBF - evaluate quadratic basis functions
//* --------------------------------------------------------------------
pub fn qbf(
  xq: f64, yq: f64, it: usize,
  inn: usize,
  node: &Vec<usize>,
  xc: &Vec<f64>, yc: &Vec<f64>, 
  nnodes: &usize
) -> (f64, f64, f64) {
  let in1;
  let in2;
  let in3;
  let i1;
  let i2;
  let i3;
  let inn_local;

  let j1;
  let j2;
  let j3;
  
  let d;
  let c;
  let t;
  let s;
  let bb;
  let bx;
  let by;

  if inn <= 2 {
    in1 = inn;
    in2 = (inn + 1) % 3;
    in3 = (inn + 2) % 3;
    i1 = node[it * nnodes + in1];
    i2 = node[it * nnodes + in2];
    i3 = node[it * nnodes + in3];
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    t =
      1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d;
    bb = t * (2.0 * t - 1.0);
    bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d;
    by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d;
  } else {
    inn_local = inn - 3;
    in1 = inn_local;
    in2 = (inn_local + 1) % 3;
    in3 = (inn_local + 2) % 3;
    i1 = node[it * nnodes + in1];
    i2 = node[it * nnodes + in2];
    i3 = node[it * nnodes + in3];
    j1 = i2;
    j2 = i3;
    j3 = i1;
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    c = (xc[j2] - xc[j1]) * (yc[j3] - yc[j1]) - (xc[j3] - xc[j1]) * (yc[j2] - yc[j1]);
    t = 1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d;
    s = 1.0 + ((yc[j2] - yc[j3]) * (xq - xc[j1]) + (xc[j3] - xc[j2]) * (yq - yc[j1])) / c;
    bb = 4.0 * s * t;
    bx = 4.0 * (t * (yc[j2] - yc[j3]) / c + s * (yc[i2] - yc[i3]) / d);
    by = 4.0 * (t * (xc[j3] - xc[j2]) / c + s * (xc[i3] - xc[i2]) / d);
  }
  return (bb, bx, by);
}

//* --------------------------------------------------------------------
//*  REFBSP - evaluate linear basis functions in a reference triangle
//* --------------------------------------------------------------------
pub fn refbsp(xq: f64, yq: f64, iq: usize) -> f64{
  if iq == 0 {
    return 1.0 - xq;
  } else if iq == 1 {
    return yq;
  } else if iq == 2 {
    return xq - yq;
  }
  return 0.0;
}

//* --------------------------------------------------------------------
//*  REFQBF - evaluate quadratic basis functions on reference triangle
//* --------------------------------------------------------------------
pub fn refqbf(
  x: f64, y: f64, inn: usize, 
  etax: f64, etay: f64, 
  xix: f64, xiy: f64
) -> (f64, f64, f64) {
  let bb: f64;
  let tbx: f64;
  let tby: f64;
  let bx: f64;
  let by: f64;

  if inn == 0 {
    bb = 1.0 - 3.0 * x + 2.0 * x * x;
    tbx = -3.0 + 4.0 * x;
    tby = 0.0;
  } else if inn == 1 {
    bb = -y + 2.0 * y * y;
    tbx = 0.0;
    tby = -1.0 + 4.0 * y;
  } else if inn == 2 {
    bb = -x + 2.0 * x * x + y - 4.0 * x * y + 2.0 * y * y;
    tbx = -1.0 + 4.0 * x - 4.0 * y;
    tby = 1.0 - 4.0 * x + 4.0 * y;
  } else if inn == 3 {
    bb = 4.0 * y - 4.0 * x * y;
    tbx = -4.0 * y;
    tby = 4.0 - 4.0 * x;
  } else if inn == 4 {
    bb = 4.0 * x * y - 4.0 * y * y;
    tbx = 4.0 * y;
    tby = 4.0 * x - 8.0 * y;
  } else if inn == 5 {
    bb = 4.0 * x - 4.0 * x * x - 4.0 * y + 4.0 * x * y;
    tbx = 4.0 - 8.0 * x + 4.0 * y;
    tby = -4.0 + 4.0 * x;
  } else {
    println!("REFQBF - Illegal value of IN = {}", inn + 1);
    std::process::exit(1);
  }

  bx = tbx * xix + tby * etax;
  by = tbx * xiy + tby * etay;

  return (bb, bx, by)
}

////* --------------------------------------------------------------------
////*  RESID - compute the residual
////* --------------------------------------------------------------------
//void resid(
//  std::array<double, bump.nelemn> &area, std::vector<double> &g,
//  std::array<int, bump.n_points * 2> &indx,
//  std::array<int, bump.n_points> &insc, std::array<int, bump.nelemn> &isotri,
//  int bump.iwrite, int bump.neqn,
//  std::array<int, bump.nelemn * bump.nnodes> &node,
//  std::array<double, bump.nelemn * bump.nquad * bump.nnodes * 3> &phi,
//  std::array<double, bump.nelemn * bump.nquad * bump.nnodes> &psi,
//  std::vector<double> &res, double reynld,
//  std::array<double, bump.n_points> &xc, std::array<double, bump.nelemn * bump.nquad> &xm,
//  std::array<double, bump.n_points> &yc, std::array<double, bump.nelemn * bump.nquad> &ym
//) {
//  int itype = -1;
//  double visc = 1.0 / reynld;
//
//  for (int it = 0; it < bump.neqn; it++) {
//    res[it] = 0.0;
//  }
//
//  double ubc = 0.0;
//  double ar, yq, xq, det, etax, etay, xix, xiy, aij;
//  int ip, iprs, ihor, iver, ipp, ju, jv, jp;
//  double bb, bx, by, bbl, bbb, bbx, bby, bbbl;
//
//  std::array<double, 2> un{};
//  std::array<double, 2> unx{};
//  std::array<double, 2> uny{};
//
//  for (int it = 0; it < bump.nelemn; it++) {
//    ar = bump.area[it] / 3.0;
//
//    for (int iquad = 0; iquad < bump.nquad; iquad++) {
//      yq = bump.ym[it * bump.nquad + iquad];
//      xq = bump.xm[it * bump.nquad + iquad];
//
//      det = 0.0;
//      etax = 0.0;
//      etay = 0.0;
//      xix = 0.0;
//      xiy = 0.0;
//
//      if (bump.isotri[it] == 1) {
//        std::tie(det, etax, etay, xix, xiy) = trans(
//          it, bump.node, xc, xq, yc, yq
//        );
//        ar = det * bump.area[it] / 3.0;
//      }
//
//      std::tie(un, unx, uny) = uval(
//        etax, etay, g, bump.indx, bump.isotri, it,
//        bump.neqn, bump.node, xc,
//        xix, xiy, xq, yc, yq
//      );
//
//      for (int iq = 0; iq < bump.nnodes; iq++) {
//        ip = bump.node[it * bump.nnodes + iq];
//        bb = bump.phi[phi_idx(it, iquad, iq, 0)];
//        bx = bump.phi[phi_idx(it, iquad, iq, 1)];
//        by = bump.phi[phi_idx(it, iquad, iq, 2)];
//        bbl = psi[psi_idx(it, iquad, iq)];
//        iprs = bump.insc[ip];
//        ihor = bump.indx[ip * 2];
//        iver = bump.indx[ip * 2 + 1];
//
//        if (0 < ihor)
//          res[ihor - 1] -= ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
//
//        if (0 < iver)
//          res[iver - 1] -= ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
//
//        for (int iqq = 0; iqq < bump.nnodes; iqq++) {
//          ipp = bump.node[it * bump.nnodes + iqq];
//          bbb = bump.phi[phi_idx(it, iquad, iqq, 0)];
//          bbx = bump.phi[phi_idx(it, iquad, iqq, 1)];
//          bby = bump.phi[phi_idx(it, iquad, iqq, 2)];
//          bbbl = psi[psi_idx(it, iquad, iqq)];
//          ju = bump.indx[ipp * 2];
//          jv = bump.indx[ipp * 2 + 1];
//          jp = bump.insc[ipp];
//
//          if (0 < ju) {
//            if (0 < ihor) {
//              res[ihor - 1] += (
//                ar
//                * (
//                  visc * (by * bby + bx * bbx)
//                  + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
//                )
//                * g[ju - 1]
//              );
//            }
//
//            if (0 < iver)
//              res[iver - 1] += ar * bb * bbb * unx[1] * g[ju - 1];
//
//            if (0 < iprs)
//              res[iprs - 1] += ar * bbx * bbl * g[ju - 1];
//
//          } else if (ju == itype) {
//            if (ju == -2) {
//              ubc = ubump(
//                g, bump.indx, ipp, iqq, bump.isotri, it, 1,
//                bump.neqn, bump.node,
//                xc, yc
//              );
//            } else if (ju == -1)
//              ubc = ubdry(1, bump.yc[ipp]);
//
//            if (0 < ihor) {
//              aij = ar * (
//                visc * (by * bby + bx * bbx)
//                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
//              );
//              res[ihor - 1] += ubc * aij;
//            }
//
//            if (0 < iver) {
//              aij = ar * bb * bbb * unx[1];
//              res[iver - 1] += ubc * aij;
//            }
//
//            if (0 < iprs) {
//              aij = ar * bbx * bbl;
//              res[iprs - 1] += ubc * aij;
//            }
//          }
//          if (0 < jv) {
//            if (0 < ihor)
//              res[ihor - 1] += ar * bb * bbb * uny[0] * g[jv - 1];
//
//            if (0 < iver) {
//              res[iver - 1] += (
//                ar
//                * (
//                  visc * (by * bby + bx * bbx)
//                  + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
//                )
//                * g[jv - 1]
//              );
//            }
//
//            if (0 < iprs)
//              res[iprs - 1] += ar * bby * bbl * g[jv - 1];
//
//          } else if (jv == itype) {
//            if (jv == -2) {
//              ubc = ubump(
//                g, bump.indx, ipp, iqq, bump.isotri, it, 2,
//                bump.neqn, bump.node,
//                xc, yc
//              );
//            } else if (jv == -1) ubc = ubdry(2, bump.yc[ipp]);
//
//            if (0 < ihor) {
//              aij = ar * bb * bbb * uny[0];
//              res[ihor - 1] += ubc * aij;
//            }
//
//            if (0 < iver) {
//              aij = ar * (
//                visc * (by * bby + bx * bbx)
//                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
//              );
//              res[iver - 1] += ubc * aij;
//            }
//
//            if (0 < iprs) {
//              aij = ar * bby * bbl;
//              res[iprs - 1] += ubc * aij;
//            }
//          }
//          if (0 < jp) {
//            if (0 < ihor)
//              res[ihor - 1] -= ar * bx * bbbl * g[jp - 1];
//
//            if (0 < iver)
//              res[iver - 1] -= ar * by * bbbl * g[jp - 1];
//          }
//        }
//      }
//    }
//  }
//
//  res[bump.neqn - 1] = g[bump.neqn - 1];
//
//  double rmax = 0.0;
//  int imax = 0;
//  int ibad = 0;
//  double test;
//
//  for (int i = 0; i < bump.neqn; i++) {
//    test = std::abs(res[i]);
//    if (rmax < test) {
//      rmax = test;
//      imax = i;
//    }
//    if (1.0e-3 < test)
//      ibad += 1;
//  }
//
//  if (1 <= bump.iwrite) {
//    println!();
//    println!("RESIDUAL INFORMATION:");
//    println!();
//    println!("Worst residual is number {}", imax+1);
//    println!("of magnitude {}", rmax);
//    println!();
//    println!("Number of \"bad\" residuals is {} out of {}", ibad, bump.neqn);
//    println!();
//  }
//
//  if (2 <= bump.iwrite) {
//    println!("Raw residuals:");
//    println!();
//    int idx = 0;
//    for (int j = 0; j < bump.n_points; j++) {
//      if (0 < bump.indx[j * 2]) {
//        if (std::abs(res[idx]) <= 1.0e-3)
//          println!(" U {} {} {}", idx+1, j+1, res[idx]);
//        else
//          println!("*U {} {} {}", idx+1, j+1, res[idx]);
//        idx += 1;
//      }
//
//      if (0 < bump.indx[j * 2 + 1]) {
//        if (std::abs(res[idx]) <= 1.0e-3)
//          println!(" V {} {} {}", idx+1, j+1, res[idx]);
//        else
//          println!("*V {} {} {}", idx+1, j+1, res[idx]);
//        idx += 1;
//      }
//
//      if (0 < bump.insc[j]) {
//        if (std::abs(res[idx]) <= 1.0e-3)
//          println!(" P {} {} {}", idx+1, j+1, res[idx]);
//        else
//          println!("*P {} {} {}", idx+1, j+1, res[idx]);
//        idx += 1;
//      }
//    }
//  }
//}
//
//
////* --------------------------------------------------------------------
////*  SETBAN - compute the half band width
////* --------------------------------------------------------------------
//std::tuple<int, int, int> setban(
//  std::array<int, bump.n_points * 2> &indx,
//  std::array<int, bump.n_points> &insc, int maxrow,
//  int nband, int nlband,
//  std::array<int, bump.nelemn * bump.nnodes> &node,
//  int nrow
//) {
//  int ip, i_val, ipp, j_val;
//  nlband = 0;
//
//  for (int it = 0; it < bump.nelemn; it++) {
//    for (int iq = 0; iq < bump.nnodes; iq++) {
//      ip = bump.node[it * bump.nnodes + iq];
//      for (int iuk = 0; iuk < 3; iuk++) {
//        if (iuk == 2)
//          i_val = bump.insc[ip];
//        else
//          i_val = bump.indx[ip * 2 + iuk];
//        if (0 < i_val) {
//          for (int iqq = 0; iqq < bump.nnodes; iqq++) {
//            ipp = bump.node[it * bump.nnodes + iqq];
//            for (int iukk = 0; iukk < 3; iukk++) {
//              if (iukk == 2)
//                j_val = bump.insc[ipp];
//              else
//                j_val = bump.indx[ipp * 2 + iukk];
//              if (0 < j_val)
//                nlband = max(nlband, j_val - i_val);
//            }
//          }
//        }
//      }
//    }
//  }
//
//  nband = nlband + nlband + 1;
//  nrow = nlband + nlband + nlband + 1;
//
//  println!();
//  println!("SETBAN:");
//  println!();
//  println!("  Lower bandwidth = {}", nlband);
//  println!("  Total bandwidth = {}", nband);
//  println!("  Required matrix rows = {}", nrow);
//
//  if (maxrow < nrow) {
//    println!("SETBAN - NROW is too large!");
//    println!("The maximum allowed is {}", maxrow);
//    println!("This problem requires NROW = {}", nrow);
//    std::process::exit(1);
//  }
//
//  return std::tuple{nband, nlband, nrow};
//}

//* --------------------------------------------------------------------
//*  SETBAS - evaluate basis functions at each integration point
//* --------------------------------------------------------------------
pub fn setbas(
  bump: &mut bump_struct::Bump
) {
  let mut xq;
  let mut yq;
  let mut det;
  let mut etax;
  let mut etay;
  let mut xix;
  let mut xiy;
  let mut bb;
  let mut bx;
  let mut by;

  for it in 0..bump.nelemn {
    for j in 0..bump.nquad {
      xq = bump.xm[it * bump.nquad + j];
      yq = bump.ym[it * bump.nquad + j];
      (det, etax, etay, xix, xiy) = trans(
        it, &bump.node, &bump.xc, xq, &bump.yc, yq, &bump.nnodes
      );

      for iq in 0..bump.nnodes {
        if bump.isotri[it] == 0 {
          bump.psi[psi_idx(it, j, iq, bump.nquad, bump.nnodes)] = bsp(
            it, iq, 0, &bump.node, &bump.xc, xq, &bump.yc, yq, &bump.nnodes
          );
          (bb, bx, by) = qbf(xq, yq, it, iq, &bump.node, &bump.xc, &bump.yc, &bump.nnodes);
        } else {
          (bb, bx, by) = refqbf(xq, yq, iq, etax, etay, xix, xiy);
          bump.psi[psi_idx(it, j, iq, bump.nquad, bump.nnodes)] = refbsp(xq, yq, iq);
        }
        bump.phi[phi_idx(it, j, iq, 0, bump.nquad, bump.nnodes)] = bb;
        bump.phi[phi_idx(it, j, iq, 1, bump.nquad, bump.nnodes)] = bx;
        bump.phi[phi_idx(it, j, iq, 2, bump.nquad, bump.nnodes)] = by;
      }
    }
  }
}

//* --------------------------------------------------------------------
//*  SETGRD - set up the geometric grid
//* --------------------------------------------------------------------
pub fn setgrd(
  ibump: i32, bump: &mut bump_struct::Bump
)/*  -> (bool, usize) */ {
  println!();
  println!("SETGRD:");
  println!();

  if bump.ny < bump.nx {
    bump._long = true;
    println!("Using vertical ordering.");
  } else{
    bump._long = false;
    println!("Using horizontal ordering.");
  }

  if ibump == 0 {
    println!("No isoparametric elements will be used.");
  } else if ibump == 1 {
    println!("Isoparametric elements directly on bump.");
  } else if ibump == 2 {
    println!("All elements above bump are isoparametric.");
  } else if ibump == 3 {
    println!("All elements are isoparametric.");
  } else {
    println!("Unexpected value of IBUMP = {}", ibump);
    std::process::exit(1);
  }

  let nbleft: usize = (bump.xbleft * (bump.mx as f64 - 1.0) / bump.xlngth).round() as usize;
  let nbrite: usize = (bump.xbrite * (bump.mx as f64 - 1.0) / bump.xlngth).round() as usize;
  println!("Bump extends from {} at bump.node {}", bump.xbleft, nbleft + 1);
  println!("               to {} at bump.node {}", bump.xbrite, nbrite + 1);

  bump.neqn = 0;
  let mut ielemn = 0;

  for ip in 0..bump.n_points {
    let ic: usize; 
    let jc: usize;
    if bump._long {
      ic = ip / bump.my;
      jc = ip % bump.my;
    } else {
      ic = ip % bump.mx;
      jc = ip / bump.mx;
    }
    let icnt: usize = (ic + 1) % 2;
    let jcnt: usize = (jc + 1) % 2;

    if (icnt == 1 && jcnt == 1) && (ic != bump.mx - 1) && (jc != bump.my - 1) {
      let ip1: usize;
      let ip2: usize;
      if bump._long {
        ip1 = ip + bump.my;
        ip2 = ip + bump.my + bump.my;

        bump.node[ielemn *  bump.nnodes] = ip;
        bump.node[ielemn *  bump.nnodes + 1] = ip + 2;
        bump.node[ielemn *  bump.nnodes + 2] = ip2 + 2;
        bump.node[ielemn *  bump.nnodes + 3] = ip + 1;
        bump.node[ielemn *  bump.nnodes + 4] = ip1 + 2;
        bump.node[ielemn *  bump.nnodes + 5] = ip1 + 1;

        if ibump == 0 {
          bump.isotri[ielemn] = 0;
        } else if ibump == 1 {
          bump.isotri[ielemn] = 0;
        } else if ibump == 2 {
          bump.isotri[ielemn] = if nbleft <= ic && ic < nbrite {1} else {0};
        } else {
          bump.isotri[ielemn] = 1;
        }

        ielemn += 1;

        bump.node[ielemn *  bump.nnodes] = ip;
        bump.node[ielemn *  bump.nnodes + 1] = ip2 + 2;
        bump.node[ielemn *  bump.nnodes + 2] = ip2;
        bump.node[ielemn *  bump.nnodes + 3] = ip1 + 1;
        bump.node[ielemn *  bump.nnodes + 4] = ip2 + 1;
        bump.node[ielemn *  bump.nnodes + 5] = ip1;

        if ibump == 0 {
          bump.isotri[ielemn] = 0;
        } else if ibump == 1 {
          bump.isotri[ielemn] = if jc == 0 && nbleft <= ic  && ic < nbrite {1} else {0};
        } else if ibump == 2 {
          bump.isotri[ielemn] = if nbleft <= ic && ic < nbrite {1} else {0};
        } else {
          bump.isotri[ielemn] = 1;
        }

        ielemn += 1;
      } else {
        ip1 = ip + bump.mx;
        ip2 = ip + bump.mx + bump.mx;

        bump.node[ielemn * bump.nnodes] = ip;
        bump.node[ielemn * bump.nnodes + 1] = ip2;
        bump.node[ielemn * bump.nnodes + 2] = ip2 + 2;
        bump.node[ielemn * bump.nnodes + 3] = ip1;
        bump.node[ielemn * bump.nnodes + 4] = ip2 + 1;
        bump.node[ielemn * bump.nnodes + 5] = ip1 + 1;

        if ibump == 0 {
          bump.isotri[ielemn] = 0;
        } else if ibump == 1 {
          bump.isotri[ielemn] = 0;
        } else if ibump == 2 {
          bump.isotri[ielemn] = if nbleft <= ic && ic < nbrite {1} else {0};
        } else {
          bump.isotri[ielemn] = 1;
        }

        ielemn += 1;

        bump.node[ielemn * bump.nnodes] = ip;
        bump.node[ielemn * bump.nnodes + 1] = ip2 + 2;
        bump.node[ielemn * bump.nnodes + 2] = ip + 2;
        bump.node[ielemn * bump.nnodes + 3] = ip1 + 1;
        bump.node[ielemn * bump.nnodes + 4] = ip1 + 2;
        bump.node[ielemn * bump.nnodes + 5] = ip + 1;

        if ibump == 0 {
          bump.isotri[ielemn] = 0;
        } else if ibump == 1 {
          bump.isotri[ielemn] = if jc == 0 && nbleft <= ic && ic < nbrite {1} else {0};
        } else if ibump == 2 {
          bump.isotri[ielemn] = if nbleft <= ic && ic < nbrite {1} else {0};
        } else {
          bump.isotri[ielemn] = 1;
        }

        ielemn += 1;
      }
    }

    if ic == 0 && 0 < jc && jc < bump.my - 1 {
      bump.indx[ip * 2] = -1;
      bump.indx[ip * 2 + 1] = -1;
    } else if ic == bump.mx - 1 && 0 < jc && jc < bump.my - 1 {
      bump.neqn += 1;
      bump.indx[ip * 2] = bump.neqn as i32;
      bump.indx[ip * 2 + 1] = 0;
    } else if jc == 0 && ielemn > 0 && bump.isotri[ielemn - 1] == 1 {
      bump.indx[ip * 2] = -2;
      bump.indx[ip * 2 + 1] = -2;
    } else if ic == 0 || ic == bump.mx - 1 || jc == 0 || jc == bump.my - 1 {
      bump.indx[ip * 2] = 0;
      bump.indx[ip * 2 + 1] = 0;
    } else {
      bump.neqn += 2;
      bump.indx[ip * 2] = (bump.neqn - 1) as i32;
      bump.indx[ip * 2 + 1] = bump.neqn as i32;
    }
    
    if jcnt == 1 && icnt == 1 {
      bump.neqn += 1;
      bump.insc[ip] = bump.neqn as i32;
    } else {
      bump.insc[ip] = 0;
    }
  }
  if 1 <= bump.iwrite {
    println!();
    println!("     I     bump.indx 1, bump.indx 2, bump.insc");
    println!();
    for i in 0..bump.n_points {
      println!(
        "{}\t{}\t{}\t{}", 
        i+1, bump.indx[i * 2], bump.indx[i * 2 + 1], bump.insc[i]
      );
    }
    println!();
    println!("Isoparametric triangles:");
    println!();
    for i in 0..bump.nelemn {
      if bump.isotri[i] == 1 {
        println!("{}", i + 1);
      }
    }
    println!();
    println!("   IT   bump.node(IT,*)");
    println!();
    for it in 0..bump.nelemn {
      print!("{}\t", it+1);
      for i in 0..6 {
        print!("{}\t", bump.node[it * bump.nnodes + i]+1);
      }
      println!();
    }
  }

  println!();
  println!("SETGRD: Number of unknowns = {}", bump.neqn);

  if bump.maxeqn < bump.neqn {
    println!("SETGRD - Too many unknowns!");
    println!("The maximum allowed is MAXEQN = {}", bump.maxeqn);
    println!("This problem requires bump.neqn = {}", bump.neqn);
    std::process::exit(1);
  }
  //return (_long, bump.neqn);
}

////* --------------------------------------------------------------------
////*  SETLIN - determine unknown numbers along the profile line
////* --------------------------------------------------------------------
//void setlin(
//  std::array<int, my> &iline,
//  std::array<int, bump.n_points * 2> &indx,
//  int bump.iwrite, bool _long,
//  double xlngth, double xprof
//) {
//  int itemp, bump.nodex0, ip;
//
//  itemp = std::round(2.0 * (nx - 1) * xprof / xlngth);
//
//  if (_long)
//    bump.nodex0 = itemp * (2 * ny - 1);
//  else
//    bump.nodex0 = itemp;
//
//  println!();
//  println!("SETLIN:");
//  println!();
//  println!("  Profile generated at X = {}", xprof);
//  println!("  which is above bump.node  = {}", bump.nodex0 + 1);
//
//  for (int i = 0; i < my; i++) {
//    if (_long)
//      ip = bump.nodex0 + i;
//    else
//      ip = bump.nodex0 + bump.mx * i;
//    iline[i] = bump.indx[ip * 2];
//  }
//
//  if (1 <= bump.iwrite) {
//    println!();
//    println!("  Indices of unknowns along the profile line:");
//    println!();
//    for (int i = 0; i < my; i += 5) {
//      for (int j = i; j < min(i + 5, my); ++j) {
//        print!("{:5d}", iline[j]);
//      }
//      println!("");
//    }
//  }
//}

//* --------------------------------------------------------------------
//*  SETQUD - set midpoint quadrature rule information
//* --------------------------------------------------------------------
pub fn setqud(
  bump: &mut bump_struct::Bump
) {
  for it in 0..bump.nelemn {
    let ip1: usize = bump.node[it * bump.nnodes];
    let ip2: usize = bump.node[it * bump.nnodes + 1];
    let ip3: usize = bump.node[it * bump.nnodes + 2];
    let x1: f64 = bump.xc[ip1];
    let x2: f64 = bump.xc[ip2];
    let x3: f64 = bump.xc[ip3];
    let y1: f64 = bump.yc[ip1];
    let y2: f64 = bump.yc[ip2];
    let y3: f64 = bump.yc[ip3];

    if bump.isotri[it] == 0 {
      bump.xm[it * bump.nquad] = 0.5 * (x1 + x2);
      bump.xm[it * bump.nquad + 1] = 0.5 * (x2 + x3);
      bump.xm[it * bump.nquad + 2] = 0.5 * (x3 + x1);
      bump.ym[it * bump.nquad] = 0.5 * (y1 + y2);
      bump.ym[it * bump.nquad + 1] = 0.5 * (y2 + y3);
      bump.ym[it * bump.nquad + 2] = 0.5 * (y3 + y1);
      bump.area[it] = 0.5 * (
        (y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3)
      ).abs();
    } else {
      bump.xm[it * bump.nquad] = 0.5;
      bump.ym[it * bump.nquad] = 0.5;
      bump.xm[it * bump.nquad + 1] = 1.0;
      bump.ym[it * bump.nquad + 1] = 0.5;
      bump.xm[it * bump.nquad + 2] = 0.5;
      bump.ym[it * bump.nquad + 2] = 0.0;
      bump.area[it] = 0.5;
    }
  }
  if 3 <= bump.iwrite {
    println!();
    println!("SETQUD: Element Areas && Quadrature points:");
    println!();
    for i in 0..bump.nelemn {
      println!("{}\t{}", i + 1, bump.area[i]);
      for j in 0..bump.nquad {
        println!("{}\t{}\t{}\t{}", i + 1, j + 1, bump.xm[i * bump.nquad + j], bump.ym[i * bump.nquad + j]);
      }
    }
  }
}

//* --------------------------------------------------------------------
//*  SETXY - set the grid coordinates based on the parameter value
//* --------------------------------------------------------------------
pub fn setxy(
  bump: &mut bump_struct::Bump
) {
  for ip in 0..bump.n_points {
    let ic;
    let jc;
    if bump._long {
      ic = ip / bump.my;
      jc = ip % bump.my;
    } else {
      ic = ip % bump.mx;
      jc = ip / bump.mx;
    }
    bump.xc[ip] = (ic as f64) * bump.xlngth / ((bump.nx + bump.nx) as f64 - 2.0);

    let ybot: f64 = -bump.ypert * (bump.xc[ip] - 3.0) * (bump.xc[ip] - 1.0);
    let ylo: f64 = (0.0_f64).max(ybot);

    bump.yc[ip] = ((bump.my - 1 - jc) as f64 * ylo + (jc as f64) * bump.ylngth) / ((bump.ny + bump.ny) as f64 - 2.0);
  }
  if 2 <= bump.iwrite {
    println!();
    println!("SETXY:");
    println!();
    println!("     I     XC     YC");
    println!();
    for i in 0..bump.n_points {
      println!("{}\t{}\t{}", i+1, bump.xc[i], bump.yc[i]);
    }
  }
}

//* --------------------------------------------------------------------
//*  TIMESTAMP - println! current date && time
//* --------------------------------------------------------------------
pub fn timestamp() {
  let local_time: chrono::prelude::DateTime<Local> = Local::now();
  println!("Local time: {}", local_time.format("%Y-%m-%d %H:%M:%S"));
}

//* --------------------------------------------------------------------
//*  TRANS - calculate the element transformation mapping
//* --------------------------------------------------------------------
pub fn trans(
  it: usize,
  node: &Vec<usize>,
  xc: &Vec<f64>, xq: f64,
  yc: &Vec<f64>, yq: f64, 
  nnodes: &usize
) -> (f64, f64, f64, f64, f64) {
  let i1 = node[it * nnodes];
  let i2 = node[it * nnodes + 1];
  let i3 = node[it * nnodes + 2];
  let i4 = node[it * nnodes + 3];
  let i5 = node[it * nnodes + 4];
  let i6 = node[it * nnodes + 5];

  let x1 = xc[i1];
  let y1 = yc[i1];
  let x2 = xc[i2];
  let y2 = yc[i2];
  let x3 = xc[i3];
  let y3 = yc[i3];
  let x4 = xc[i4];
  let y4 = yc[i4];
  let x5 = xc[i5];
  let y5 = yc[i5];
  let x6 = xc[i6];
  let y6 = yc[i6];

  let a1 = 2.0 * x3 - 4.0 * x6 + 2.0 * x1;
  let b1 = -4.0 * x3 - 4.0 * x4 + 4.0 * x5 + 4.0 * x6;
  let c1 = 2.0 * x2 + 2.0 * x3 - 4.0 * x5;
  let d1 = -3.0 * x1 - x3 + 4.0 * x6;
  let e1 = -x2 + x3 + 4.0 * x4 - 4.0 * x6;

  let a2 = 2.0 * y3 - 4.0 * y6 + 2.0 * y1;
  let b2 = -4.0 * y3 - 4.0 * y4 + 4.0 * y5 + 4.0 * y6;
  let c2 = 2.0 * y2 + 2.0 * y3 - 4.0 * y5;
  let d2 = -3.0 * y1 - y3 + 4.0 * y6;
  let e2 = -y2 + y3 + 4.0 * y4 - 4.0 * y6;

  let dxdxi = 2.0 * a1 * xq + b1 * yq + d1;
  let dxdeta = b1 * xq + 2.0 * c1 * yq + e1;
  let dydxi = 2.0 * a2 * xq + b2 * yq + d2;
  let dydeta = b2 * xq + 2.0 * c2 * yq + e2;

  let mut det =
    (2.0 * a1 * b2 - 2.0 * a2 * b1) * xq * xq
    + (4.0 * a1 * c2 - 4.0 * a2 * c1) * xq * yq
    + (2.0 * b1 * c2 - 2.0 * b2 * c1) * yq * yq
    + (2.0 * a1 * e2 + b2 * d1 - b1 * d2 - 2.0 * a2 * e1) * xq
    + (2.0 * c2 * d1 + b1 * e2 - b2 * e1 - 2.0 * c1 * d2) * yq
    + d1 * e2
    - d2 * e1;

  let eps = 1e-30;
  if det.abs() < eps {
    det = eps;
  }

  let xix: f64 = dydeta / det;
  let xiy: f64 = -dxdeta / det;
  let etax: f64 = -dydxi / det;
  let etay: f64 = dxdxi / det;

  return (det, etax, etay, xix, xiy);
}

////* --------------------------------------------------------------------
////*  UBDRY - parabolic inflow boundary condition
////* --------------------------------------------------------------------
//inline double ubdry(int iuk, double yy) {
//  if (iuk == 1)
//    return (-2.0 * yy + 6.0) * yy / 9.0;
//  else
//    return 0.0;
//}
//
////* --------------------------------------------------------------------
////*  UBUMP - sensitivity dU/dA on the bump
////* --------------------------------------------------------------------
//double ubump(
//  const std::vector<double> &g,
//  const std::array<int, bump.n_points * 2> &indx,
//  int ip, int iqq,
//  const std::array<int, bump.nelemn> &isotri,
//  int it, int iukk,
//  int bump.neqn,
//  const std::array<int, bump.nelemn * bump.nnodes> &node,
//  const std::array<double, bump.n_points> &xc, const std::array<double, bump.n_points> &yc
//) {
//  double det = 1.0;
//  double etax = 0.0;
//  double etay = 0.0;
//  double xix = 1.0;
//  double xiy = 1.0;
//  double xq = 0.0;
//  double yq = 0.0;
//
//  if (bump.isotri[it] == 0) {
//    xq = bump.xc[ip];
//    yq = bump.yc[ip];
//  } else {
//    if (iqq == 0) {
//      xq = 0.0;
//      yq = 0.0;
//    } else if (iqq == 1) {
//      xq = 1.0;
//      yq = 1.0;
//    } else if (iqq == 2) {
//      xq = 1.0;
//      yq = 0.0;
//    } else if (iqq == 3) {
//      xq = 0.5;
//      yq = 0.5;
//    } else if (iqq == 4) {
//      xq = 1.0;
//      yq = 0.5;
//    } else if (iqq == 5) {
//      xq = 0.5;
//      yq = 0.0;
//    }
//    
//    std::tie(det, etax, etay, xix, xiy) = trans(
//      it, bump.node, xc, xq, yc, yq
//    );
//  }
//
//  auto [un, unx, uny] = _ubump_uval(
//    g, bump.indx, bump.isotri, it,
//    bump.neqn, bump.node,
//    xc, xix, xiy, xq, yc, yq, 
//    det, etax, etay
//  );
//
//  if (iukk == 1)
//    return -uny[0] * (bump.xc[ip] - 1.0) * (bump.xc[ip] - 3.0);
//  else if (iukk == 2)
//    return -uny[1] * (bump.xc[ip] - 1.0) * (bump.xc[ip] - 3.0);
//  else {
//    println!("UBUMP called for iukk = {}", iukk);
//    std::process::exit(1);
//  }
//}
//std::tuple<std::array<double, 2>, std::array<double, 2>, std::array<double, 2>> _ubump_uval(
//  const std::vector<double> &g,
//  const std::array<int, bump.n_points * 2> &indx,
//  const std::array<int, bump.nelemn> &isotri,
//  int it, int bump.neqn,
//  const std::array<int, bump.nelemn * bump.nnodes> &node,
//  const std::array<double, bump.n_points> &xc,
//  double xix, double xiy, double xq,
//  const std::array<double, bump.n_points> &yc,
//  double yq, double det, double etax, double etay
//) {
//  std::array<double, 2> un{};
//  std::array<double, 2> unx{};
//  std::array<double, 2> uny{};
//  
//  double bb, bx, by, ubc;
//
//  for (int iq = 0; iq < bump.nnodes; iq++) {
//    if (bump.isotri[it] == 1)
//      std::tie(bb, bx, by) = refqbf(xq, yq, iq, etax, etay, xix, xiy);
//    else
//      std::tie(bb, bx, by) = qbf(xq, yq, it, iq, bump.node, xc, yc);
//    int ip_local = bump.node[it * bump.nnodes + iq];
//
//    for (int iuk = 0; iuk < 2; iuk++) {
//      int iun = bump.indx[ip_local * 2 + iuk];
//      if (0 < iun) {
//        un[iuk] += bb * g[iun - 1];
//        unx[iuk] += bx * g[iun - 1];
//        uny[iuk] += by * g[iun - 1];
//      } else if (iun == -1) {
//        ubc = ubdry(iuk + 1, bump.yc[ip_local]);
//        un[iuk] += bb * ubc;
//        unx[iuk] += bx * ubc;
//        uny[iuk] += by * ubc;
//      }
//    }
//  }
//
//  return std::tuple {un, unx, uny};
//}
//
////* --------------------------------------------------------------------
////*  UVAL - evaluate velocities at a given quadrature point
////* --------------------------------------------------------------------
//std::tuple<std::array<double, 2>, std::array<double, 2>, std::array<double, 2>> uval(
//  double etax,
//  double etay,
//  const std::vector<double> &g,
//  const std::array<int, bump.n_points * 2> &indx,
//  const std::array<int, bump.nelemn> &isotri,
//  int it,
//  int bump.neqn,
//  const std::array<int, bump.nelemn * bump.nnodes> &node,
//  const std::array<double, bump.n_points> &xc,
//  double xix,
//  double xiy,
//  double xq,
//  const std::array<double, bump.n_points> &yc,
//  double yq
//) {
//  std::array<double, 2> un{};
//  std::array<double, 2> unx{};
//  std::array<double, 2> uny{};
//  double bb, bx, by, ubc;
//
//  for (int iq = 0; iq < bump.nnodes; iq++) {
//    if (bump.isotri[it] == 1)
//      std::tie(bb, bx, by) = refqbf(xq, yq, iq, etax, etay, xix, xiy);
//    else
//      std::tie(bb, bx, by) = qbf(xq, yq, it, iq, bump.node, xc, yc);
//    
//    int ip_local = bump.node[it * bump.nnodes + iq];
//
//    for (int iuk = 0; iuk < 2; iuk++) {
//      int iun = bump.indx[ip_local * 2 + iuk];
//      if (0 < iun) {
//        un[iuk] += bb * g[iun - 1];
//        unx[iuk] += bx * g[iun - 1];
//        uny[iuk] += by * g[iun - 1];
//      } else if (iun == -1) {
//        ubc = ubdry(iuk + 1, bump.yc[ip_local]);
//        un[iuk] += bb * ubc;
//        unx[iuk] += bx * ubc;
//        uny[iuk] += by * ubc;
//      }
//    }
//  }
//
//  return std::tuple {un, unx, uny};
//}
//
////* --------------------------------------------------------------------
////*  UV_WRITE - write a velocity file
////* --------------------------------------------------------------------
//void uv_write(
//  const std::vector<double>& f, 
//  const std::array<int, bump.n_points * 2>& bump.indx, 
//  std::ofstream& uv_file_obj, 
//  int bump.neqn,
//  const std::array<double, bump.n_points>& yc
//) {
//  
//  uv_file_obj << std::scientific << std::setprecision(6);
//  double u, v;
//
//  for (int ip = 0; ip < bump.n_points; ++ip) {
//    int k_u = bump.indx[ip * 2];
//
//    if (k_u < 0) u = ubdry(1, bump.yc[ip]);
//    else if (k_u == 0) u = 0.0;
//    else u = f[k_u - 1];
//
//    int k_v = bump.indx[ip * 2 + 1];
//
//    if (k_v < 0) v = ubdry(2, bump.yc[ip]);
//    else if (k_v == 0) v = 0.0;
//    else v = f[k_v - 1];
//
//    uv_file_obj << "  " << std::setw(14) << u 
//                << "  " << std::setw(14) << v << "\n";
//  }
//}
//
////* --------------------------------------------------------------------
////*  XY_WRITE - write bump.node coordinate data
////* --------------------------------------------------------------------
//void xy_write(
//  std::ofstream& xy_file_obj,
//  std::array<double, bump.n_points>& xc, std::array<double, bump.n_points>& yc
//) {
//  xy_file_obj << std::scientific << std::setprecision(6);
//  
//  for (int ip = 0; ip < bump.n_points; ++ip) {
//    xy_file_obj << "  " << std::setw(14) << bump.xc[ip] 
//                << "  " << std::setw(14) << bump.yc[ip] << "\n";
//  }
//}
