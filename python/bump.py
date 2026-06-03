"""
BUMP - Control problem for channel flow over a bump.

Translation from Fortran90 to Python (using only NumPy).

Original Fortran code by John Burkardt.
Discretization of Cost and Sensitivities in Shape Optimization,
in Computation and Control IV, Birkhaeuser, 1995.

Navier Stokes equations (primitive variable formulation):
  U dUdx + V dUdy + dPdx - mu*(ddU/dxdx + ddU/dydy) = F1
  U dVdx + V dVdy + dPdy - mu*(ddV/dxdx + ddV/dydy) = F2
  dUdx + dVdy = 0
"""

import numpy as np
import time
import datetime
import sys
import os


# --------------------------------------------------------------------
#  Main program
# --------------------------------------------------------------------
def main() -> int:
  maxnew = 4
  maxsec = 10
  nx = 21
  ny = 7
  maxrow = 27 * ny
  nelemn = 2 * (nx - 1) * (ny - 1)
  mx = 2 * nx - 1
  my = 2 * ny - 1
  maxeqn = 2 * mx * my + nx * ny
  n_points = mx * my
  nnodes = 6
  nquad = 3

  a: np.ndarray[tuple[int, int], np.dtype[np.float64]] = np.zeros(
    (maxrow, maxeqn), dtype=np.float64
  )
  anew = 0.0
  anext = 0.3
  aold = 0.0
  aprof = 0.25
  area = np.zeros(nelemn, dtype=np.float64)
  dcda = np.zeros(my, dtype=np.float64)
  f = np.zeros(maxeqn, dtype=np.float64)
  g = np.zeros(maxeqn, dtype=np.float64)
  gr = np.zeros((my, my), dtype=np.float64)
  iline = np.zeros(my, dtype=np.int32)
  indx = np.zeros((n_points, 2), dtype=np.int32)
  insc = np.zeros(n_points, dtype=np.int32)
  isotri = np.zeros(nelemn, dtype=np.int32)
  iwrite = 10
  long = False
  nband = 0
  neqn = 0
  nlband = 0
  node = np.zeros((nelemn, nnodes), dtype=np.int32)
  nrow = 0
  numnew = 0
  numsec = 0
  phi = np.zeros((nelemn, nquad, nnodes, 3), dtype=np.float64)
  psi = np.zeros((nelemn, nquad, nnodes), dtype=np.float64)
  r = np.zeros(my, dtype=np.float64)
  res = np.zeros(maxeqn, dtype=np.float64)
  reynld = 1.0
  rjpnew = 0.0
  rjpold = 0.0
  sens = np.zeros(maxeqn, dtype=np.float64)
  tolnew = 0.0001
  tolsec = 0.0001
  uprof = np.zeros(my, dtype=np.float64)
  uv_dir = "data\\uv"
  uv_file = f"{uv_dir}\\uv_000.txt"
  xbleft = 1.0
  xbrite = 3.0
  xc = np.zeros(n_points, dtype=np.float64)
  xlngth = 10.0
  xm = np.zeros((nelemn, nquad), dtype=np.float64)
  xprof = 4.0
  xy_dir = "data\\xy"
  xy_file = f"{xy_dir}\\xy_000.txt"
  yc = np.zeros(n_points, dtype=np.float64)
  ylngth = 3.0
  ym = np.zeros((nelemn, nquad), dtype=np.float64)
  ypert = 0.0

  cpu1 = time.time()

  if not os.path.exists(uv_dir):
    os.makedirs(uv_dir)
  if not os.path.exists(xy_dir):
    os.makedirs(xy_dir)

  timestamp()
  print()
  print("BUMP")
  print("  Python version")
  print("  Control problem for channel flow over a bump.")
  print()
  print(f"  The bump will be generated with a height of {aprof}")
  print()
  print(f"  NX = {nx}")
  print(f"  NY = {ny}")
  print(f"  Number of elements = {nelemn}")
  print(f"  Reynolds number =  {reynld}")
  print(f"  Secant tolerance = {tolsec}")
  print(f"  Newton tolerance = {tolnew}")

  ibump = 2
  long, neqn = setgrd(
    ibump,
    indx,
    insc,
    isotri,
    iwrite,
    long,
    maxeqn,
    mx,
    my,
    nelemn,
    neqn,
    nnodes,
    node,
    n_points,
    nx,
    ny,
    xbleft,
    xbrite,
    xlngth,
  )

  ypert = aprof
  setxy(iwrite, long, mx, my, n_points, nx, ny, xc, xlngth, yc, ylngth, ypert)

  setqud(area, isotri, iwrite, nelemn, nnodes, node, n_points, nquad, xc, xm, yc, ym)

  setbas(isotri, nelemn, nnodes, node, n_points, nquad, phi, psi, xc, xm, yc, ym)

  setlin(iline, indx, iwrite, long, mx, my, n_points, nx, ny, xlngth, xprof)

  nband, nlband, nrow = setban(
    indx, insc, maxrow, nband, nelemn, nlband, nnodes, node, n_points, nrow
  )

  g[:neqn] = 0.0

  numnew = nstoke(
    a,
    area,
    f,
    g,
    indx,
    insc,
    isotri,
    maxnew,
    maxrow,
    nband,
    nelemn,
    neqn,
    nlband,
    nnodes,
    node,
    n_points,
    nquad,
    nrow,
    numnew,
    phi,
    psi,
    reynld,
    tolnew,
    xc,
    xm,
    yc,
    ym,
  )

  resid(
    area,
    g,
    indx,
    insc,
    isotri,
    iwrite,
    nelemn,
    neqn,
    nnodes,
    node,
    n_points,
    nquad,
    phi,
    psi,
    res,
    reynld,
    xc,
    xm,
    yc,
    ym,
  )

  uprof = getg(g, iline, my, neqn)

  if 1 <= iwrite:
    print()
    print("Velocity profile:")
    print()
    for i in range(my):
      print(f"{uprof[i]:14.6e}", end="")
      if (i + 1) % 5 == 0:
        print()
    print()

  gram(
    gr, iline, indx, iwrite, my, nelemn, nnodes, node, n_points, r, uprof, xc, xprof, yc
  )

  xy_file = file_name_inc(xy_file)
  # xy_unit = get_unit()
  with open(xy_file, "w") as f_xy:
    xy_write(f_xy, n_points, xc, yc)

  uv_file = file_name_inc(uv_file)
  # uv_unit = get_unit()
  with open(uv_file, "w") as f_uv:
    uv_write(f, indx, f_uv, neqn, n_points, yc)

  g[:neqn] = 0.0

  for iter in range(1, maxsec + 1):
    print()
    print(f"Secant iteration {iter}")

    numsec += 1

    ypert = anew
    setxy(iwrite, long, mx, my, n_points, nx, ny, xc, xlngth, yc, ylngth, ypert)

    setqud(area, isotri, iwrite, nelemn, nnodes, node, n_points, nquad, xc, xm, yc, ym)

    setbas(isotri, nelemn, nnodes, node, n_points, nquad, phi, psi, xc, xm, yc, ym)

    numnew = nstoke(
      a,
      area,
      f,
      g,
      indx,
      insc,
      isotri,
      maxnew,
      maxrow,
      nband,
      nelemn,
      neqn,
      nlband,
      nnodes,
      node,
      n_points,
      nquad,
      nrow,
      numnew,
      phi,
      psi,
      reynld,
      tolnew,
      xc,
      xm,
      yc,
      ym,
    )

    uprof = getg(g, iline, my, neqn)

    if 1 <= iwrite:
      print()
      print("Velocity profile:")
      print()
      for i in range(my):
        print(f"{uprof[i]:14.6e}", end="")
        if (i + 1) % 5 == 0:
          print()
      print()

    itype = -2

    sens = linsys(
      a,
      area,
      sens,
      g,
      indx,
      insc,
      isotri,
      itype,
      maxrow,
      nband,
      nelemn,
      neqn,
      nlband,
      nnodes,
      node,
      n_points,
      nquad,
      nrow,
      phi,
      psi,
      reynld,
      xc,
      xm,
      yc,
      ym,
    )

    dcda = getg(sens, iline, my, neqn)

    if 2 <= iwrite:
      print()
      print("Sensitivities:")
      print()
      for i in range(my):
        print(f"{dcda[i]:14.6e}", end="")
        if (i + 1) % 5 == 0:
          print()
      print()

    rjpnew = 0.0
    for i in range(my):
      temp = -r[i]
      for j in range(my):
        temp += gr[i, j] * uprof[j]
      rjpnew += 2.0 * dcda[i] * temp

    xy_file = file_name_inc(xy_file)
    # xy_unit = get_unit()
    with open(xy_file, "w") as f_xy:
      xy_write(f_xy, n_points, xc, yc)

    uv_file = file_name_inc(uv_file)
    # uv_unit = get_unit()
    with open(uv_file, "w") as f_uv:
      uv_write(f, indx, f_uv, neqn, n_points, yc)

    print()
    print(f"  Parameter = {anew}, J prime = {rjpnew}")

    if 1 < iter:
      denom = rjpnew - rjpold
      if abs(denom) > 1e-30:
        anext = aold - rjpold * (anew - aold) / denom
      else:
        anext = anew

    aold = anew
    anew = anext
    rjpold = rjpnew

    if anew != 0.0:
      test = abs(anew - aold) / anew
    else:
      test = 0.0

    print(f"  New value of parameter = {anew}")
    print(f"  Convergence test = {test}")

    if abs(anew - aold) <= abs(anew) * tolsec and 1 < iter:
      print("Secant iteration converged.")
      break
    else:
      print("  Secant iteration failed to converge.")

  cpu2 = time.time()
  print()
  print(f"  Total execution time = {cpu2 - cpu1} seconds.")
  print(f"  Number of secant steps = {numsec}")
  print(f"  Number of Newton steps = {numnew}")
  print()
  print("BUMP:")
  print("  Normal end of execution.")
  print()
  timestamp()

  return 0


# --------------------------------------------------------------------
#  BSP - linear basis function for pressure
# --------------------------------------------------------------------
def bsp(it, iq, id, nelemn, nnodes, node, n_points, xc, xq, yc, yq):
  l1 = iq
  l2 = i4_wrap(iq + 1, 0, 2)  # 0-indexed: 1->2, 2->0, 3->1
  l3 = i4_wrap(iq + 2, 0, 2)

  g1 = node[it, l1]
  g2 = node[it, l2]
  g3 = node[it, l3]

  d = (xc[g2] - xc[g1]) * (yc[g3] - yc[g1]) - (xc[g3] - xc[g1]) * (yc[g2] - yc[g1])

  if id == 0:  # 0-indexed: id == 1 -> 0
    return (
      1.0 + ((yc[g2] - yc[g3]) * (xq - xc[g1]) + (xc[g3] - xc[g2]) * (yq - yc[g1])) / d
    )
  elif id == 1:
    return (yc[g2] - yc[g3]) / d
  elif id == 2:
    return (xc[g3] - xc[g2]) / d
  else:
    print()
    print("BSP - Fatal error!")
    print("  Illegal local index value for linear basis.")
    print("  Legal values are 1, 2 or 3.")
    print(f"  The input value was ID = {id + 1}")
    sys.exit(1)


# --------------------------------------------------------------------
#  DAXPY - constant times a vector plus a vector
# --------------------------------------------------------------------
def daxpy(n, da, dx, incx, dy, incy):
  if n <= 0 or da == 0.0:
    return

  if incx != 1 or incy != 1:
    if incx >= 0:
      ix = 0
    else:
      ix = (-n + 1) * incx
    if incy >= 0:
      iy = 0
    else:
      iy = (-n + 1) * incy
    for _ in range(n):
      dy[iy] += da * dx[ix]
      ix += incx
      iy += incy
  else:
    m = n % 4
    for i in range(m):
      dy[i] += da * dx[i]
    for i in range(m, n, 4):
      dy[i] += da * dx[i]
      dy[i + 1] += da * dx[i + 1]
      dy[i + 2] += da * dx[i + 2]
      dy[i + 3] += da * dx[i + 3]


# --------------------------------------------------------------------
#  DDOT - dot product of two vectors
# --------------------------------------------------------------------
def ddot(n, dx, incx, dy, incy):
  if n <= 0:
    return 0.0

  dtemp = 0.0

  if incx != 1 or incy != 1:
    if incx >= 0:
      ix = 0
    else:
      ix = (-n + 1) * incx
    if incy >= 0:
      iy = 0
    else:
      iy = (-n + 1) * incy
    for _ in range(n):
      dtemp += dx[ix] * dy[iy]
      ix += incx
      iy += incy
  else:
    m = n % 5
    for i in range(m):
      dtemp += dx[i] * dy[i]
    for i in range(m, n, 5):
      dtemp += (
        dx[i] * dy[i]
        + dx[i + 1] * dy[i + 1]
        + dx[i + 2] * dy[i + 2]
        + dx[i + 3] * dy[i + 3]
        + dx[i + 4] * dy[i + 4]
      )

  return dtemp


# --------------------------------------------------------------------
#  DGBFA - factor a real band matrix by elimination
# --------------------------------------------------------------------
def dgbfa(abd, lda, n, ml, mu, ipvt):
  m = ml + mu + 1
  info = 0

  j0 = mu + 1  # 0-indexed adjustment
  j1 = min(n, m) - 1

  for jz in range(j0, j1):
    i0 = m - jz
    for i in range(i0, ml):
      abd[i, jz] = 0.0

  jz = j1
  ju = 0

  for k in range(n - 1):
    jz += 1
    if jz < n:
      abd[0:ml, jz] = 0.0

    lm = min(ml, n - k - 1)
    l = idamax(lm + 1, abd[m - 1 :, k], 1) + (m - 1)
    ipvt[k] = l + k - (m - 1)

    if abd[l, k] == 0.0:
      info = k
    else:
      if l != m - 1:
        t = abd[l, k]
        abd[l, k] = abd[m - 1, k]
        abd[m - 1, k] = t

      t = -1.0 / abd[m - 1, k]
      dscal(lm, t, abd[m:, k], 1)

      ju = min(max(ju, mu + ipvt[k]), n - 1)
      mm = m - 1

      for j in range(k + 1, ju + 1):
        l -= 1
        mm -= 1
        t = abd[l, j]
        if l != mm:
          abd[l, j] = abd[mm, j]
          abd[mm, j] = t
        daxpy(lm, t, abd[m:, k], 1, abd[mm + 1 :, j], 1)

  ipvt[n - 1] = n - 1

  if abd[m - 1, n - 1] == 0.0:
    info = n - 1

  return info, ipvt


# --------------------------------------------------------------------
#  DGBSL - solve a real banded system factored by DGBFA
# --------------------------------------------------------------------
def dgbsl(abd, lda, n, ml, mu, ipvt, b, job):
  m = mu + ml + 1

  if job == 0:
    if 0 < ml:
      for k in range(n - 1):
        lm = min(ml, n - k - 1)
        l = ipvt[k]
        t = b[l]
        if l != k:
          b[l] = b[k]
          b[k] = t
        daxpy(lm, t, abd[m:, k], 1, b[k + 1 :], 1)

    for k in range(n - 1, -1, -1):
      b[k] /= abd[m - 1, k]
      lm = min(k, m - 1)
      la = m - 1 - lm
      lb = k - lm
      t = -b[k]
      daxpy(lm, t, abd[la:, k], 1, b[lb:], 1)
  else:
    for k in range(n):
      lm = min(k, m - 1)
      la = m - 1 - lm
      lb = k - lm
      t = ddot(lm, abd[la:, k], 1, b[lb:], 1)
      b[k] = (b[k] - t) / abd[m - 1, k]

    if 0 < ml:
      for k in range(n - 2, -1, -1):
        lm = min(ml, n - k - 1)
        b[k] += ddot(lm, abd[m:, k], 1, b[k + 1 :], 1)
        l = ipvt[k]
        if l != k:
          t = b[l]
          b[l] = b[k]
          b[k] = t

  return b


# --------------------------------------------------------------------
#  DSCAL - scale a vector by a constant
# --------------------------------------------------------------------
def dscal(n, sa, x, incx):
  if n <= 0:
    return

  if incx == 1:
    m = n % 5
    for i in range(m):
      x[i] *= sa
    for i in range(m, n, 5):
      x[i] *= sa
      x[i + 1] *= sa
      x[i + 2] *= sa
      x[i + 3] *= sa
      x[i + 4] *= sa
  else:
    if incx >= 0:
      ix = 0
    else:
      ix = (-n + 1) * incx
    for _ in range(n):
      x[ix] *= sa
      ix += incx


# --------------------------------------------------------------------
#  FILE_NAME_INC - increment a partially numeric filename
# --------------------------------------------------------------------
def file_name_inc(file_name: str) -> str:
  if not file_name:
    print()
    print("FILE_NAME_INC - Fatal error!")
    print("  The input string is empty.")
    sys.exit(1)

  name_list: list[str] = list(file_name)
  change = 0

  for i in range(len(name_list) - 1, -1, -1):
    c = name_list[i]
    if "0" <= c <= "9":
      change += 1
      digit = ord(c) - 48
      digit += 1
      if digit == 10:
        digit = 0
      c = chr(digit + 48)
      name_list[i] = c
      if c != "0":
        return "".join(name_list)

  if change == 0:
    return ""

  return "".join(name_list)


# --------------------------------------------------------------------
#  GET_UNIT - return a free FORTRAN unit number (mocked for Python)
# --------------------------------------------------------------------
# _next_unit = 10
# def get_unit():
#  global _next_unit
#  _next_unit += 1
#  return _next_unit


# --------------------------------------------------------------------
#  GETG - extract values of a quantity along the profile line
# --------------------------------------------------------------------
def getg(f, iline, my, neqn):
  u = np.zeros(my, dtype=np.float64)
  for i in range(my):
    j = iline[i]
    if j < 0:  # 0-indexed: j <= 0 means j == -1 (Fortran: j <= 0)
      u[i] = 0.0
    elif j == 0:
      u[i] = 0.0
    else:
      u[i] = f[j - 1]  # Convert 1-based index from iline to 0-based
  return u


# --------------------------------------------------------------------
#  GRAM - compute the Gram matrix and R vector
# --------------------------------------------------------------------
def gram(
  gr, iline, indx, iwrite, my, nelemn, nnodes, node, n_points, r, uprof, xc, xprof, yc
):
  wt = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0], dtype=np.float64)
  yq_gauss = np.array([-0.7745966692, 0.0, 0.7745966692], dtype=np.float64)

  r[:] = 0.0
  gr[:, :] = 0.0

  for it in range(nelemn):
    k = node[it, 0]
    kk = node[it, 1]

    if abs(xc[k] - xprof) > 1.0e-4 or abs(xc[kk] - xprof) > 1.0e-4:
      continue

    for iquad in range(3):
      bma2 = (yc[kk] - yc[k]) / 2.0
      ar = bma2 * wt[iquad]
      x = xprof
      y = yc[k] + bma2 * (yq_gauss[iquad] + 1.0)

      uiqdpt = 0.0
      for iq in range(nnodes):
        if iq in (0, 1, 3):
          bb, bx, by = qbf(x, y, it, iq, nelemn, nnodes, node, n_points, xc, yc)
          ip = node[it, iq]
          iun = indx[ip, 0]
          if 0 < iun:
            ii = igetl(iun, iline, my)
            uiqdpt += bb * uprof[ii - 1]
          elif iun == -1:
            ubc = ubdry(1, yc[ip])
            uiqdpt += bb * ubc

      for iq in range(nnodes):
        if iq in (0, 1, 3):
          ip = node[it, iq]
          bb, bx, by = qbf(x, y, it, iq, nelemn, nnodes, node, n_points, xc, yc)
          i_val = indx[ip, 0]
          if 0 < i_val:
            ii = igetl(i_val, iline, my)
            r[ii - 1] += bb * uiqdpt * ar
            for iqq in range(nnodes):
              if iqq in (0, 1, 3):
                ipp = node[it, iqq]
                bbb, bbx_, bby_ = qbf(
                  x, y, it, iqq, nelemn, nnodes, node, n_points, xc, yc
                )
                j_val = indx[ipp, 0]
                if j_val != 0:
                  jj = igetl(j_val, iline, my)
                  gr[ii - 1, jj - 1] += bb * bbb * ar

  if 3 <= iwrite:
    print()
    print("Gram matrix:")
    print()
    for i in range(my):
      for j in range(my):
        print(i + 1, j + 1, gr[i, j])
    print()
    print("R vector:")
    print()
    for i in range(my):
      print(r[i])


# --------------------------------------------------------------------
#  I4_MODP - nonnegative remainder of integer division
# --------------------------------------------------------------------
def i4_modp(i, j):
  if j == 0:
    print()
    print("I4_MODP - Fatal error!")
    print(f"  Illegal divisor J = {j}")
    sys.exit(1)
  value = i % j
  if value < 0:
    value += abs(j)
  return value


# --------------------------------------------------------------------
#  I4_WRAP - force an integer between limits by wrapping
# --------------------------------------------------------------------
def i4_wrap(ival, ilo, ihi):
  jlo = min(ilo, ihi)
  jhi = max(ilo, ihi)
  wide = jhi - jlo + 1
  if wide == 1:
    return jlo
  return jlo + i4_modp(ival - jlo, wide)


# --------------------------------------------------------------------
#  IDAMAX - index of the vector element of maximum absolute value
# --------------------------------------------------------------------
def idamax(n, dx, incx):
  if n < 1 or incx <= 0:
    return -1

  idamax_val = 0

  if n == 1:
    return 0

  if incx == 1:
    dmax = abs(dx[0])
    for i in range(1, n):
      if dmax < abs(dx[i]):
        idamax_val = i
        dmax = abs(dx[i])
  else:
    ix = 0
    dmax = abs(dx[0])
    ix += incx
    for i in range(1, n):
      if dmax < abs(dx[ix]):
        idamax_val = i
        dmax = abs(dx[ix])
      ix += incx

  return idamax_val


# --------------------------------------------------------------------
#  IGETL - get the local unknown number along the profile line
# --------------------------------------------------------------------
def igetl(i, iline, my):
  for j in range(my):
    if iline[j] == i:
      return j + 1  # 1-indexed
  return -1


# --------------------------------------------------------------------
#  LINSYS - solve the linearized Navier Stokes equation
# --------------------------------------------------------------------
def linsys(
  a,
  area,
  f,
  g,
  indx,
  insc,
  isotri,
  itype,
  maxrow,
  nband,
  nelemn,
  neqn,
  nlband,
  nnodes,
  node,
  n_points,
  nquad,
  nrow,
  phi,
  psi,
  reynld,
  xc,
  xm,
  yc,
  ym,
):
  ioff = (
    nlband + nlband + 1
  )  # This is 1-indexed offset to main diagonal row in band storage
  visc = 1.0 / reynld

  f[:neqn] = 0.0
  a[:nrow, :neqn] = 0.0

  un = np.zeros(2, dtype=np.float64)
  unx = np.zeros(2, dtype=np.float64)
  uny = np.zeros(2, dtype=np.float64)

  for it in range(nelemn):
    ar = area[it] / 3.0

    for iquad in range(nquad):
      yq = ym[it, iquad]
      xq = xm[it, iquad]

      det = 0.0
      etax = 0.0
      etay = 0.0
      xix = 0.0
      xiy = 0.0
      ubc = 0.0

      if isotri[it] == 1:
        det, etax, etay, xix, xiy = trans(
          it, nelemn, nnodes, node, n_points, xc, xq, yc, yq
        )
        ar = det * area[it] / 3.0

      uval_out = uval(
        etax,
        etay,
        g,
        indx,
        isotri,
        it,
        nelemn,
        neqn,
        nnodes,
        node,
        n_points,
        xc,
        xix,
        xiy,
        xq,
        yc,
        yq,
      )
      un[:], unx[:], uny[:] = uval_out

      for iq in range(nnodes):
        ip = node[it, iq]
        bb = phi[it, iquad, iq, 0]
        bx = phi[it, iquad, iq, 1]
        by = phi[it, iquad, iq, 2]
        bbl = psi[it, iquad, iq]
        ihor = indx[ip, 0]
        iver = indx[ip, 1]
        iprs = insc[ip]

        if 0 < ihor:
          f[ihor - 1] += ar * bb * (un[0] * unx[0] + un[1] * uny[0])

        if 0 < iver:
          f[iver - 1] += ar * bb * (un[0] * unx[1] + un[1] * uny[1])

        for iqq in range(nnodes):
          ipp = node[it, iqq]
          bbb = phi[it, iquad, iqq, 0]
          bbx = phi[it, iquad, iqq, 1]
          bby = phi[it, iquad, iqq, 2]
          bbbl = psi[it, iquad, iqq]
          ju = indx[ipp, 0]
          jv = indx[ipp, 1]
          jp = insc[ipp]

          if 0 < ju:
            if 0 < ihor:
              iuse = ihor - ju + ioff
              a[iuse - 1, ju - 1] += ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
              )

            if 0 < iver:
              iuse = iver - ju + ioff
              a[iuse - 1, ju - 1] += ar * bb * bbb * unx[1]

            if 0 < iprs:
              iuse = iprs - ju + ioff
              a[iuse - 1, ju - 1] += ar * bbx * bbl

          elif ju == itype:
            if ju == -1:
              ubc = ubdry(1, yc[ipp])
            elif ju == -2:
              ubc = ubump(
                g,
                indx,
                ipp,
                iqq,
                isotri,
                it,
                1,
                nelemn,
                neqn,
                nnodes,
                node,
                n_points,
                xc,
                yc,
              )

            if 0 < ihor:
              aij = ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
              )
              f[ihor - 1] -= ubc * aij

            if 0 < iver:
              aij = ar * bb * bbb * unx[1]
              f[iver - 1] -= ubc * aij

            if 0 < iprs:
              aij = ar * bbx * bbl
              f[iprs - 1] -= ubc * aij

          if 0 < jv:
            if 0 < ihor:
              iuse = ihor - jv + ioff
              a[iuse - 1, jv - 1] += ar * bb * bbb * uny[0]

            if 0 < iver:
              iuse = iver - jv + ioff
              a[iuse - 1, jv - 1] += ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
              )

            if 0 < iprs:
              iuse = iprs - jv + ioff
              a[iuse - 1, jv - 1] += ar * bby * bbl

          elif jv == itype:
            if jv == -1:
              ubc = ubdry(2, yc[ipp])
            elif jv == -2:
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
                yc,
              )

            if 0 < ihor:
              aij = ar * bb * bbb * uny[0]
              f[ihor - 1] -= ubc * aij

            if 0 < iver:
              aij = ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
              )
              f[iver - 1] -= ubc * aij

            if 0 < iprs:
              aij = ar * bby * bbl
              f[iprs - 1] -= ubc * aij

          if 0 < jp:
            if 0 < ihor:
              iuse = ihor - jp + ioff
              a[iuse - 1, jp - 1] -= ar * bx * bbbl

            if 0 < iver:
              iuse = iver - jp + ioff
              a[iuse - 1, jp - 1] -= ar * by * bbbl

  f[neqn - 1] = 0.0
  for j_1based in range(neqn - nlband, neqn):
    j = j_1based - 1
    i = neqn - j_1based + ioff  # 1-indexed row in band storage
    a[i - 1, j] = 0.0
  a[ioff - 1, neqn - 1] = 1.0  # nband == ioff, set diagonal of last eqn to 1

  # Save right-hand side in f, solve in-place
  # DGBFA and DGBSL work on A in LINPACK band format

  # We need to call dgbfa with the right parameters
  # a is (maxrow, neqn) but the band part is (nrow, neqn)
  # In LINPACK, LDA = maxrow, N = neqn, ML = nlband, MU = nlband
  ipvt = np.zeros(neqn, dtype=np.int32)
  info, ipvt = dgbfa(a, maxrow, neqn, nlband, nlband, ipvt)

  if info != 0:
    print()
    print("LINSYS - fatal error!")
    print(f"DGBFA returns INFO = {info + 1}")  # +1 for 1-indexed
    sys.exit(1)

  job = 0
  f = dgbsl(a, maxrow, neqn, nlband, nlband, ipvt, f, job)

  return f


# --------------------------------------------------------------------
#  NSTOKE - solve Navier Stokes using Newton iteration
# --------------------------------------------------------------------
def nstoke(
  a,
  area,
  f,
  g,
  indx,
  insc,
  isotri,
  maxnew,
  maxrow,
  nband,
  nelemn,
  neqn,
  nlband,
  nnodes,
  node,
  n_points,
  nquad,
  nrow,
  numnew,
  phi,
  psi,
  reynld,
  tolnew,
  xc,
  xm,
  yc,
  ym,
):
  for iter in range(maxnew):
    numnew += 1

    itype = -1
    f = linsys(
      a,
      area,
      f,
      g,
      indx,
      insc,
      isotri,
      itype,
      maxrow,
      nband,
      nelemn,
      neqn,
      nlband,
      nnodes,
      node,
      n_points,
      nquad,
      nrow,
      phi,
      psi,
      reynld,
      xc,
      xm,
      yc,
      ym,
    )

    g[:neqn] = g[:neqn] - f[:neqn]

    imax = idamax(neqn, g, 1)
    diff = abs(g[imax])
    print(f"NSTOKE: Iteration {iter + 1}, MaxNorm(diff) = {diff}")

    g[:neqn] = f[:neqn]

    if diff <= tolnew:
      print("NSTOKE converged.")
      return numnew

    if iter == maxnew - 1:
      print("NSTOKE failed!")
      sys.exit(1)

  return numnew


# --------------------------------------------------------------------
#  PVAL - compute pressure table at all nodes
# --------------------------------------------------------------------
def pval(g, insc, long, mx, my, nelemn, neqn, nnodes, node, n_points):
  press = np.zeros((mx, my), dtype=np.float64)

  for it in range(nelemn):
    for iq in range(3):
      ip = node[it, iq]
      ivar = insc[ip]
      if long:
        i = (ip) // my
        j = (ip) % my
      else:
        i = (ip) % mx
        j = (ip) // mx

      if 0 < ivar:
        press[i, j] = g[ivar - 1]

  for i in range(1, mx - 1, 2):
    for j in range(0, my, 2):
      press[i, j] = 0.5 * (press[i - 1, j] + press[i + 1, j])

  for j in range(1, my - 1, 2):
    for i in range(0, mx, 2):
      press[i, j] = 0.5 * (press[i, j - 1] + press[i, j + 1])

  for j in range(1, my - 1, 2):
    for i in range(1, mx - 1, 2):
      press[i, j] = 0.5 * (press[i - 1, j - 1] + press[i + 1, j + 1])

  return press


# --------------------------------------------------------------------
#  QBF - evaluate quadratic basis functions
# --------------------------------------------------------------------
def qbf(xq, yq, it, inn, nelemn, nnodes, node, n_points, xc, yc):
  if inn <= 2:  # 0-indexed: in=0,1,2 corresponds to Fortran in=1,2,3
    in1 = inn
    in2 = (inn + 1) % 3
    in3 = (inn + 2) % 3
    i1 = node[it, in1]
    i2 = node[it, in2]
    i3 = node[it, in3]
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1])
    t = (
      1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d
    )
    bb = t * (2.0 * t - 1.0)
    bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d
    by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d
  else:
    inn_local = inn - 3
    in1 = inn_local
    in2 = (inn_local + 1) % 3
    in3 = (inn_local + 2) % 3
    i1 = node[it, in1]
    i2 = node[it, in2]
    i3 = node[it, in3]
    j1 = i2
    j2 = i3
    j3 = i1
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1])
    c = (xc[j2] - xc[j1]) * (yc[j3] - yc[j1]) - (xc[j3] - xc[j1]) * (yc[j2] - yc[j1])
    t = (
      1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d
    )
    s = (
      1.0 + ((yc[j2] - yc[j3]) * (xq - xc[j1]) + (xc[j3] - xc[j2]) * (yq - yc[j1])) / c
    )
    bb = 4.0 * s * t
    bx = 4.0 * (t * (yc[j2] - yc[j3]) / c + s * (yc[i2] - yc[i3]) / d)
    by = 4.0 * (t * (xc[j3] - xc[j2]) / c + s * (xc[i3] - xc[i2]) / d)

  return bb, bx, by


# --------------------------------------------------------------------
#  REFBSP - evaluate linear basis functions in a reference triangle
# --------------------------------------------------------------------
def refbsp(xq, yq, iq):
  # 0-indexed: iq=0 -> psi1, iq=1 -> psi2, iq=2 -> psi3
  if iq == 0:
    return 1.0 - xq
  elif iq == 1:
    return yq
  elif iq == 2:
    return xq - yq
  return 0.0


# --------------------------------------------------------------------
#  REFQBF - evaluate quadratic basis functions on reference triangle
# --------------------------------------------------------------------
def refqbf(x, y, inn, etax, etay, xix, xiy):
  # 0-indexed: inn=0..5 corresponds to Fortran in=1..6
  if inn == 0:
    bb = 1.0 - 3.0 * x + 2.0 * x * x
    tbx = -3.0 + 4.0 * x
    tby = 0.0
  elif inn == 1:
    bb = -y + 2.0 * y * y
    tbx = 0.0
    tby = -1.0 + 4.0 * y
  elif inn == 2:
    bb = -x + 2.0 * x * x + y - 4.0 * x * y + 2.0 * y * y
    tbx = -1.0 + 4.0 * x - 4.0 * y
    tby = 1.0 - 4.0 * x + 4.0 * y
  elif inn == 3:
    bb = 4.0 * y - 4.0 * x * y
    tbx = -4.0 * y
    tby = 4.0 - 4.0 * x
  elif inn == 4:
    bb = 4.0 * x * y - 4.0 * y * y
    tbx = 4.0 * y
    tby = 4.0 * x - 8.0 * y
  elif inn == 5:
    bb = 4.0 * x - 4.0 * x * x - 4.0 * y + 4.0 * x * y
    tbx = 4.0 - 8.0 * x + 4.0 * y
    tby = -4.0 + 4.0 * x
  else:
    print(f"REFQBF - Illegal value of IN = {inn + 1}")
    sys.exit(1)

  bx = tbx * xix + tby * etax
  by = tbx * xiy + tby * etay

  return bb, bx, by


# --------------------------------------------------------------------
#  RESID - compute the residual
# --------------------------------------------------------------------
def resid(
  area,
  g,
  indx,
  insc,
  isotri,
  iwrite,
  nelemn,
  neqn,
  nnodes,
  node,
  n_points,
  nquad,
  phi,
  psi,
  res,
  reynld,
  xc,
  xm,
  yc,
  ym,
):
  itype = -1
  visc = 1.0 / reynld

  res[:neqn] = 0.0
  ubc = 0.0

  un = np.zeros(2, dtype=np.float64)
  unx = np.zeros(2, dtype=np.float64)
  uny = np.zeros(2, dtype=np.float64)

  for it in range(nelemn):
    ar = area[it] / 3.0

    for iquad in range(nquad):
      yq = ym[it, iquad]
      xq = xm[it, iquad]

      det = 0.0
      etax = 0.0
      etay = 0.0
      xix = 0.0
      xiy = 0.0

      if isotri[it] == 1:
        det, etax, etay, xix, xiy = trans(
          it, nelemn, nnodes, node, n_points, xc, xq, yc, yq
        )
        ar = det * area[it] / 3.0

      uval_out = uval(
        etax,
        etay,
        g,
        indx,
        isotri,
        it,
        nelemn,
        neqn,
        nnodes,
        node,
        n_points,
        xc,
        xix,
        xiy,
        xq,
        yc,
        yq,
      )
      un[:], unx[:], uny[:] = uval_out

      for iq in range(nnodes):
        ip = node[it, iq]
        bb = phi[it, iquad, iq, 0]
        bx = phi[it, iquad, iq, 1]
        by = phi[it, iquad, iq, 2]
        bbl = psi[it, iquad, iq]
        iprs = insc[ip]
        ihor = indx[ip, 0]
        iver = indx[ip, 1]

        if 0 < ihor:
          res[ihor - 1] -= ar * bb * (un[0] * unx[0] + un[1] * uny[0])

        if 0 < iver:
          res[iver - 1] -= ar * bb * (un[0] * unx[1] + un[1] * uny[1])

        for iqq in range(nnodes):
          ipp = node[it, iqq]
          bbb = phi[it, iquad, iqq, 0]
          bbx = phi[it, iquad, iqq, 1]
          bby = phi[it, iquad, iqq, 2]
          bbbl = psi[it, iquad, iqq]
          ju = indx[ipp, 0]
          jv = indx[ipp, 1]
          jp = insc[ipp]

          if 0 < ju:
            if 0 < ihor:
              res[ihor - 1] += (
                ar
                * (
                  visc * (by * bby + bx * bbx)
                  + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
                )
                * g[ju - 1]
              )

            if 0 < iver:
              res[iver - 1] += ar * bb * bbb * unx[1] * g[ju - 1]

            if 0 < iprs:
              res[iprs - 1] += ar * bbx * bbl * g[ju - 1]

          elif ju == itype:
            if ju == -2:
              ubc = ubump(
                g,
                indx,
                ipp,
                iqq,
                isotri,
                it,
                1,
                nelemn,
                neqn,
                nnodes,
                node,
                n_points,
                xc,
                yc,
              )
            elif ju == -1:
              ubc = ubdry(1, yc[ipp])

            if 0 < ihor:
              aij = ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1])
              )
              res[ihor - 1] += ubc * aij

            if 0 < iver:
              aij = ar * bb * bbb * unx[1]
              res[iver - 1] += ubc * aij

            if 0 < iprs:
              aij = ar * bbx * bbl
              res[iprs - 1] += ubc * aij

          if 0 < jv:
            if 0 < ihor:
              res[ihor - 1] += ar * bb * bbb * uny[0] * g[jv - 1]

            if 0 < iver:
              res[iver - 1] += (
                ar
                * (
                  visc * (by * bby + bx * bbx)
                  + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
                )
                * g[jv - 1]
              )

            if 0 < iprs:
              res[iprs - 1] += ar * bby * bbl * g[jv - 1]

          elif jv == itype:
            if jv == -2:
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
                yc,
              )
            elif jv == -1:
              ubc = ubdry(2, yc[ipp])

            if 0 < ihor:
              aij = ar * bb * bbb * uny[0]
              res[ihor - 1] += ubc * aij

            if 0 < iver:
              aij = ar * (
                visc * (by * bby + bx * bbx)
                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0])
              )
              res[iver - 1] += ubc * aij

            if 0 < iprs:
              aij = ar * bby * bbl
              res[iprs - 1] += ubc * aij

          if 0 < jp:
            if 0 < ihor:
              res[ihor - 1] -= ar * bx * bbbl * g[jp - 1]

            if 0 < iver:
              res[iver - 1] -= ar * by * bbbl * g[jp - 1]

  res[neqn - 1] = g[neqn - 1]

  rmax = 0.0
  imax = 0
  ibad = 0

  for i in range(neqn):
    test = abs(res[i])
    if rmax < test:
      rmax = test
      imax = i
    if 1.0e-3 < test:
      ibad += 1

  if 1 <= iwrite:
    print()
    print("RESIDUAL INFORMATION:")
    print()
    print(f"Worst residual is number {imax + 1}")
    print(f"of magnitude {rmax}")
    print()
    print(f'Number of "bad" residuals is {ibad} out of {neqn}')
    print()

  if 2 <= iwrite:
    print("Raw residuals:")
    print()
    idx = 0
    for j in range(n_points):
      if 0 < indx[j, 0]:
        if abs(res[idx]) <= 1.0e-3:
          print(f" U {idx + 1} {j + 1} {res[idx]:14.6e}")
        else:
          print(f"*U {idx + 1} {j + 1} {res[idx]:14.6e}")
        idx += 1

      if 0 < indx[j, 1]:
        if abs(res[idx]) <= 1.0e-3:
          print(f" V {idx + 1} {j + 1} {res[idx]:14.6e}")
        else:
          print(f"*V {idx + 1} {j + 1} {res[idx]:14.6e}")
        idx += 1

      if 0 < insc[j]:
        if abs(res[idx]) <= 1.0e-3:
          print(f" P {idx + 1} {j + 1} {res[idx]:14.6e}")
        else:
          print(f"*P {idx + 1} {j + 1} {res[idx]:14.6e}")
        idx += 1


# --------------------------------------------------------------------
#  SETBAN - compute the half band width
# --------------------------------------------------------------------
def setban(indx, insc, maxrow, nband, nelemn, nlband, nnodes, node, n_points, nrow):
  nlband = 0

  for it in range(nelemn):
    for iq in range(nnodes):
      ip = node[it, iq]
      for iuk in range(3):
        if iuk == 2:
          i_val = insc[ip]
        else:
          i_val = indx[ip, iuk]
        if 0 < i_val:
          for iqq in range(nnodes):
            ipp = node[it, iqq]
            for iukk in range(3):
              if iukk == 2:
                j_val = insc[ipp]
              else:
                j_val = indx[ipp, iukk]
              if 0 < j_val:
                nlband = max(nlband, j_val - i_val)

  nband = nlband + nlband + 1
  nrow = nlband + nlband + nlband + 1

  print()
  print("SETBAN:")
  print()
  print(f"  Lower bandwidth = {nlband}")
  print(f"  Total bandwidth = {nband}")
  print(f"  Required matrix rows = {nrow}")

  if maxrow < nrow:
    print("SETBAN - NROW is too large!")
    print(f"The maximum allowed is {maxrow}")
    sys.exit(1)

  return nband, nlband, nrow


# --------------------------------------------------------------------
#  SETBAS - evaluate basis functions at each integration point
# --------------------------------------------------------------------
def setbas(isotri, nelemn, nnodes, node, n_points, nquad, phi, psi, xc, xm, yc, ym):
  for it in range(nelemn):
    for j in range(nquad):
      xq = xm[it, j]
      yq = ym[it, j]
      det, etax, etay, xix, xiy = trans(
        it, nelemn, nnodes, node, n_points, xc, xq, yc, yq
      )
      for iq in range(nnodes):
        if isotri[it] == 0:
          psi[it, j, iq] = bsp(
            it, iq, 0, nelemn, nnodes, node, n_points, xc, xq, yc, yq
          )
          bb, bx, by = qbf(xq, yq, it, iq, nelemn, nnodes, node, n_points, xc, yc)
        else:
          bb, bx, by = refqbf(xq, yq, iq, etax, etay, xix, xiy)
          psi[it, j, iq] = refbsp(xq, yq, iq)
        phi[it, j, iq, 0] = bb
        phi[it, j, iq, 1] = bx
        phi[it, j, iq, 2] = by


# --------------------------------------------------------------------
#  SETGRD - set up the geometric grid
# --------------------------------------------------------------------
def setgrd(
  ibump,
  indx,
  insc,
  isotri,
  iwrite,
  long,
  maxeqn,
  mx,
  my,
  nelemn,
  neqn,
  nnodes,
  node,
  n_points,
  nx,
  ny,
  xbleft,
  xbrite,
  xlngth,
):
  print()
  print("SETGRD:")
  print()

  if ny < nx:
    long = True
    print("Using vertical ordering.")
  else:
    long = False
    print("Using horizontal ordering.")

  if ibump == 0:
    print("No isoparametric elements will be used.")
  elif ibump == 1:
    print("Isoparametric elements directly on bump.")
  elif ibump == 2:
    print("All elements above bump are isoparametric.")
  elif ibump == 3:
    print("All elements are isoparametric.")
  else:
    print(f"Unexpected value of IBUMP = {ibump}")
    sys.exit(1)

  nbleft = int(round(xbleft * (mx - 1) / xlngth))  # 0-indexed
  nbrite = int(round(xbrite * (mx - 1) / xlngth))
  print(f"Bump extends from {xbleft} at node {nbleft}")
  print(f"               to {xbrite} at node {nbrite}")

  neqn = 0
  ielemn = 0

  for ip in range(n_points):
    if long:
      ic = ip // my
      jc = ip % my
    else:
      ic = ip % mx
      jc = ip // mx

    icnt = (ic + 1) % 2
    jcnt = (jc + 1) % 2

    if (icnt == 1 and jcnt == 1) and (ic != mx - 1) and (jc != my - 1):
      if long:
        ip1 = ip + my
        ip2 = ip + my + my

        node[ielemn, 0] = ip
        node[ielemn, 1] = ip + 2
        node[ielemn, 2] = ip2 + 2
        node[ielemn, 3] = ip + 1
        node[ielemn, 4] = ip1 + 2
        node[ielemn, 5] = ip1 + 1

        if ibump == 0:
          isotri[ielemn] = 0
        elif ibump == 1:
          isotri[ielemn] = 0
        elif ibump == 2:
          isotri[ielemn] = 1 if (nbleft <= ic < nbrite) else 0
        else:
          isotri[ielemn] = 1

        ielemn += 1

        node[ielemn, 0] = ip
        node[ielemn, 1] = ip2 + 2
        node[ielemn, 2] = ip2
        node[ielemn, 3] = ip1 + 1
        node[ielemn, 4] = ip2 + 1
        node[ielemn, 5] = ip1

        if ibump == 0:
          isotri[ielemn] = 0
        elif ibump == 1:
          isotri[ielemn] = 1 if (jc == 0 and nbleft <= ic < nbrite) else 0
        elif ibump == 2:
          isotri[ielemn] = 1 if (nbleft <= ic < nbrite) else 0
        else:
          isotri[ielemn] = 1

        ielemn += 1

      else:
        ip1 = ip + mx
        ip2 = ip + mx + mx

        node[ielemn, 0] = ip
        node[ielemn, 1] = ip2
        node[ielemn, 2] = ip2 + 2
        node[ielemn, 3] = ip1
        node[ielemn, 4] = ip2 + 1
        node[ielemn, 5] = ip1 + 1

        if ibump == 0:
          isotri[ielemn] = 0
        elif ibump == 1:
          isotri[ielemn] = 0
        elif ibump == 2:
          isotri[ielemn] = 1 if (nbleft <= ic < nbrite) else 0
        else:
          isotri[ielemn] = 1

        ielemn += 1

        node[ielemn, 0] = ip
        node[ielemn, 1] = ip2 + 2
        node[ielemn, 2] = ip + 2
        node[ielemn, 3] = ip1 + 1
        node[ielemn, 4] = ip1 + 2
        node[ielemn, 5] = ip + 1

        if ibump == 0:
          isotri[ielemn] = 0
        elif ibump == 1:
          isotri[ielemn] = 1 if (jc == 0 and nbleft <= ic < nbrite) else 0
        elif ibump == 2:
          isotri[ielemn] = 1 if (nbleft <= ic < nbrite) else 0
        else:
          isotri[ielemn] = 1

        ielemn += 1

    if ic == 0 and 0 < jc < my - 1:
      indx[ip, 0] = -1
      indx[ip, 1] = -1

    elif ic == mx - 1 and 0 < jc < my - 1:
      neqn += 1
      indx[ip, 0] = neqn
      indx[ip, 1] = 0

    elif jc == 0 and ielemn > 0 and isotri[ielemn - 1] == 1:
      indx[ip, 0] = -2
      indx[ip, 1] = -2

    elif ic == 0 or ic == mx - 1 or jc == 0 or jc == my - 1:
      indx[ip, 0] = 0
      indx[ip, 1] = 0

    else:
      neqn += 2
      indx[ip, 0] = neqn - 1
      indx[ip, 1] = neqn

    if jcnt == 1 and icnt == 1:
      neqn += 1
      insc[ip] = neqn
    else:
      insc[ip] = 0

  if 1 <= iwrite:
    print()
    print("     I     INDX 1, INDX 2, INSC")
    print()
    for i in range(n_points):
      print(f"{i + 1:5d}{indx[i, 0]:5d}{indx[i, 1]:5d}{insc[i]:5d}")
    print()
    print("Isoparametric triangles:")
    print()
    for i in range(nelemn):
      if isotri[i] == 1:
        print(i + 1)
    print()
    print("   IT   NODE(IT,*)")
    print()
    for it in range(nelemn):
      print(f"{it + 1:6d}", end="")
      for i in range(6):
        print(f"{node[it, i] + 1:6d}", end="")
      print()

  print()
  print(f"SETGRD: Number of unknowns = {neqn}")

  if maxeqn < neqn:
    print("SETGRD - Too many unknowns!")
    print(f"The maximum allowed is MAXEQN = {maxeqn}")
    print(f"This problem requires NEQN = {neqn}")
    sys.exit(1)

  return long, neqn


# --------------------------------------------------------------------
#  SETLIN - determine unknown numbers along the profile line
# --------------------------------------------------------------------
def setlin(iline, indx, iwrite, long, mx, my, n_points, nx, ny, xlngth, xprof):
  itemp = int(round(2.0 * (nx - 1) * xprof / xlngth))

  if long:
    nodex0 = itemp * (2 * ny - 1)
  else:
    nodex0 = itemp

  print()
  print("SETLIN:")
  print()
  print(f"  Profile generated at X = {xprof}")
  print(f"  which is above node  = {nodex0 + 1}")

  for i in range(my):
    if long:
      ip = nodex0 + i
    else:
      ip = nodex0 + mx * i
    iline[i] = indx[ip, 0]

  if 1 <= iwrite:
    print()
    print("  Indices of unknowns along the profile line:")
    print()
    for i in range(0, my, 5):
      print(" ".join(f"{iline[j]:5d}" for j in range(i, min(i + 5, my))))


# --------------------------------------------------------------------
#  SETQUD - set midpoint quadrature rule information
# --------------------------------------------------------------------
def setqud(area, isotri, iwrite, nelemn, nnodes, node, n_points, nquad, xc, xm, yc, ym):
  for it in range(nelemn):
    ip1 = node[it, 0]
    ip2 = node[it, 1]
    ip3 = node[it, 2]
    x1 = xc[ip1]
    x2 = xc[ip2]
    x3 = xc[ip3]
    y1 = yc[ip1]
    y2 = yc[ip2]
    y3 = yc[ip3]

    if isotri[it] == 0:
      xm[it, 0] = 0.5 * (x1 + x2)
      xm[it, 1] = 0.5 * (x2 + x3)
      xm[it, 2] = 0.5 * (x3 + x1)
      ym[it, 0] = 0.5 * (y1 + y2)
      ym[it, 1] = 0.5 * (y2 + y3)
      ym[it, 2] = 0.5 * (y3 + y1)
      area[it] = 0.5 * abs(
        (y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3)
      )
    else:
      xm[it, 0] = 0.5
      ym[it, 0] = 0.5
      xm[it, 1] = 1.0
      ym[it, 1] = 0.5
      xm[it, 2] = 0.5
      ym[it, 2] = 0.0
      area[it] = 0.5

  if 3 <= iwrite:
    print()
    print("SETQUD: Element Areas and Quadrature points:")
    print()
    for i in range(nelemn):
      print(i + 1, area[i])
      for j in range(nquad):
        print(i + 1, j + 1, xm[i, j], ym[i, j])


# --------------------------------------------------------------------
#  SETXY - set the grid coordinates based on the parameter value
# --------------------------------------------------------------------
def setxy(iwrite, long, mx, my, n_points, nx, ny, xc, xlngth, yc, ylngth, ypert):
  for ip in range(n_points):
    if long:
      ic = ip // my
      jc = ip % my
    else:
      ic = ip % mx
      jc = ip // mx

    xc[ip] = ic * xlngth / (2 * nx - 2)

    ybot = -ypert * (xc[ip] - 3.0) * (xc[ip] - 1.0)
    ylo = max(0.0, ybot)

    yc[ip] = ((my - 1 - jc) * ylo + jc * ylngth) / (2 * ny - 2)

  if 2 <= iwrite:
    print()
    print("SETXY:")
    print()
    print("     I     XC     YC")
    print()
    for i in range(n_points):
      print(f"{i + 1:5d}{xc[i]:12.5f}{yc[i]:12.5f}")


# --------------------------------------------------------------------
#  TIMESTAMP - print current date and time
# --------------------------------------------------------------------
def timestamp():
  now = datetime.datetime.now()
  print(now.strftime("%d %B %Y   %I:%M:%S.%f %p"))


# --------------------------------------------------------------------
#  TRANS - calculate the element transformation mapping
# --------------------------------------------------------------------
def trans(
  it, nelemn, nnodes, node, n_points, xc, xq, yc, yq
) -> tuple[float, float, float, float, float]:
  i1 = node[it, 0]
  i2 = node[it, 1]
  i3 = node[it, 2]
  i4 = node[it, 3]
  i5 = node[it, 4]
  i6 = node[it, 5]

  x1 = xc[i1]
  y1 = yc[i1]
  x2 = xc[i2]
  y2 = yc[i2]
  x3 = xc[i3]
  y3 = yc[i3]
  x4 = xc[i4]
  y4 = yc[i4]
  x5 = xc[i5]
  y5 = yc[i5]
  x6 = xc[i6]
  y6 = yc[i6]

  a1 = 2.0 * x3 - 4.0 * x6 + 2.0 * x1
  b1 = -4.0 * x3 - 4.0 * x4 + 4.0 * x5 + 4.0 * x6
  c1 = 2.0 * x2 + 2.0 * x3 - 4.0 * x5
  d1 = -3.0 * x1 - x3 + 4.0 * x6
  e1 = -x2 + x3 + 4.0 * x4 - 4.0 * x6

  a2 = 2.0 * y3 - 4.0 * y6 + 2.0 * y1
  b2 = -4.0 * y3 - 4.0 * y4 + 4.0 * y5 + 4.0 * y6
  c2 = 2.0 * y2 + 2.0 * y3 - 4.0 * y5
  d2 = -3.0 * y1 - y3 + 4.0 * y6
  e2 = -y2 + y3 + 4.0 * y4 - 4.0 * y6

  dxdxi = 2.0 * a1 * xq + b1 * yq + d1
  dxdeta = b1 * xq + 2.0 * c1 * yq + e1
  dydxi = 2.0 * a2 * xq + b2 * yq + d2
  dydeta = b2 * xq + 2.0 * c2 * yq + e2

  det = (
    (2.0 * a1 * b2 - 2.0 * a2 * b1) * xq * xq
    + (4.0 * a1 * c2 - 4.0 * a2 * c1) * xq * yq
    + (2.0 * b1 * c2 - 2.0 * b2 * c1) * yq * yq
    + (2.0 * a1 * e2 + b2 * d1 - b1 * d2 - 2.0 * a2 * e1) * xq
    + (2.0 * c2 * d1 + b1 * e2 - b2 * e1 - 2.0 * c1 * d2) * yq
    + d1 * e2
    - d2 * e1
  )

  eps = 1e-30
  if abs(det) < eps:
    det = eps

  xix = dydeta / det
  xiy = -dxdeta / det
  etax = -dydxi / det
  etay = dxdxi / det

  return det, etax, etay, xix, xiy


# --------------------------------------------------------------------
#  UBDRY - parabolic inflow boundary condition
# --------------------------------------------------------------------
def ubdry(iuk, yy):
  if iuk == 1:
    return (-2.0 * yy + 6.0) * yy / 9.0
  else:
    return 0.0


# --------------------------------------------------------------------
#  UBUMP - sensitivity dU/dA on the bump
# --------------------------------------------------------------------
def ubump(
  g, indx, ip, iqq, isotri, it, iukk, nelemn, neqn, nnodes, node, n_points, xc, yc
):
  det = 1.0
  etax = 0.0
  etay = 0.0
  xix = 1.0
  xiy = 1.0
  xq = 0.0
  yq = 0.0

  if isotri[it] == 0:
    xq = xc[ip]
    yq = yc[ip]
  else:
    if iqq == 0:
      xq = 0.0
      yq = 0.0
    elif iqq == 1:
      xq = 1.0
      yq = 1.0
    elif iqq == 2:
      xq = 1.0
      yq = 0.0
    elif iqq == 3:
      xq = 0.5
      yq = 0.5
    elif iqq == 4:
      xq = 1.0
      yq = 0.5
    elif iqq == 5:
      xq = 0.5
      yq = 0.0

    det, etax, etay, xix, xiy = trans(
      it, nelemn, nnodes, node, n_points, xc, xq, yc, yq
    )

  un = np.zeros(2, dtype=np.float64)
  unx = np.zeros(2, dtype=np.float64)
  uny = np.zeros(2, dtype=np.float64)

  un_out, unx_out, uny_out = _ubump_uval(
    g,
    indx,
    isotri,
    it,
    nelemn,
    neqn,
    nnodes,
    node,
    n_points,
    xc,
    xix,
    xiy,
    xq,
    yc,
    yq,
    det,
    etax,
    etay,
  )
  un[:], unx[:], uny[:] = un_out, unx_out, uny_out

  if iukk == 1:
    return -uny[0] * (xc[ip] - 1.0) * (xc[ip] - 3.0)
  elif iukk == 2:
    return -uny[1] * (xc[ip] - 1.0) * (xc[ip] - 3.0)
  else:
    print(f"UBUMP called for iukk = {iukk}")
    sys.exit(1)


def _ubump_uval(
  g,
  indx,
  isotri,
  it,
  nelemn,
  neqn,
  nnodes,
  node,
  n_points,
  xc,
  xix,
  xiy,
  xq,
  yc,
  yq,
  det,
  etax,
  etay,
):
  un = np.zeros(2, dtype=np.float64)
  unx = np.zeros(2, dtype=np.float64)
  uny = np.zeros(2, dtype=np.float64)

  for iq in range(nnodes):
    if isotri[it] == 1:
      bb, bx, by = refqbf(xq, yq, iq, etax, etay, xix, xiy)
    else:
      bb, bx, by = qbf(xq, yq, it, iq, nelemn, nnodes, node, n_points, xc, yc)
    ip_local = node[it, iq]

    for iuk in range(2):
      iun = indx[ip_local, iuk]
      if 0 < iun:
        un[iuk] += bb * g[iun - 1]
        unx[iuk] += bx * g[iun - 1]
        uny[iuk] += by * g[iun - 1]
      elif iun == -1:
        ubc = ubdry(iuk + 1, yc[ip_local])
        un[iuk] += bb * ubc
        unx[iuk] += bx * ubc
        uny[iuk] += by * ubc

  return un, unx, uny


# --------------------------------------------------------------------
#  UVAL - evaluate velocities at a given quadrature point
# --------------------------------------------------------------------
def uval(
  etax,
  etay,
  g,
  indx,
  isotri,
  it,
  nelemn,
  neqn,
  nnodes,
  node,
  n_points,
  xc,
  xix,
  xiy,
  xq,
  yc,
  yq,
):
  un = np.zeros(2, dtype=np.float64)
  unx = np.zeros(2, dtype=np.float64)
  uny = np.zeros(2, dtype=np.float64)

  for iq in range(nnodes):
    if isotri[it] == 1:
      bb, bx, by = refqbf(xq, yq, iq, etax, etay, xix, xiy)
    else:
      bb, bx, by = qbf(xq, yq, it, iq, nelemn, nnodes, node, n_points, xc, yc)
    ip_local = node[it, iq]

    for iuk in range(2):
      iun = indx[ip_local, iuk]
      if 0 < iun:
        un[iuk] += bb * g[iun - 1]
        unx[iuk] += bx * g[iun - 1]
        uny[iuk] += by * g[iun - 1]
      elif iun == -1:
        ubc = ubdry(iuk + 1, yc[ip_local])
        un[iuk] += bb * ubc
        unx[iuk] += bx * ubc
        uny[iuk] += by * ubc

  return un, unx, uny


# --------------------------------------------------------------------
#  UV_WRITE - write a velocity file
# --------------------------------------------------------------------
def uv_write(f, indx, uv_file_obj, neqn, n_points, yc):
  for ip in range(n_points):
    k = indx[ip, 0]
    if k < 0:
      u = ubdry(1, yc[ip])
    elif k == 0:
      u = 0.0
    else:
      u = f[k - 1]

    k = indx[ip, 1]
    if k < 0:
      v = ubdry(2, yc[ip])
    elif k == 0:
      v = 0.0
    else:
      v = f[k - 1]

    uv_file_obj.write(f"  {u:14.6e}  {v:14.6e}\n")


# --------------------------------------------------------------------
#  XY_WRITE - write node coordinate data
# --------------------------------------------------------------------
def xy_write(xy_file_obj, n_points: int, xc: np.ndarray, yc: np.ndarray) -> None:
  for ip in range(n_points):
    xy_file_obj.write(f"  {xc[ip]:14.6e}  {yc[ip]:14.6e}\n")


# --------------------------------------------------------------------
#  Run the program
# --------------------------------------------------------------------
if __name__ == "__main__":
  main()
