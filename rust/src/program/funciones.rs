#![allow(dead_code)]

use chrono::Local;

use crate::program::bump_struct;

const NNODES: usize = 6;
const NQUAD: usize = 3;

#[inline(always)]
pub fn phi_idx(it: usize, iquad: usize, iq: usize, k: usize) -> usize {
    ((it * NQUAD + iquad) * NNODES + iq) * 3 + k
}
#[inline(always)]
pub fn psi_idx(it: usize, iquad: usize, iq: usize) -> usize {
    (it * NQUAD + iquad) * NNODES + iq
}

//* --------------------------------------------------------------------
//*  BSP - linear basis function for pressure
//* --------------------------------------------------------------------
pub fn bsp(
    it: usize,
    iq: usize,
    id: usize,
    node: &[usize],
    xc: &[f64],
    xq: f64,
    yc: &[f64],
    yq: f64,
) -> f64 {
    let l1 = iq;
    let l2 = i4_wrap(iq + 1, 0, 2);
    let l3 = i4_wrap(iq + 2, 0, 2);

    let g1: usize = node[it * NNODES + l1];
    let g2: usize = node[it * NNODES + l2];
    let g3: usize = node[it * NNODES + l3];

    let d: f64 = (xc[g2] - xc[g1]) * (yc[g3] - yc[g1]) - (xc[g3] - xc[g1]) * (yc[g2] - yc[g1]);

    if id == 0 {
        return 1.0 + ((yc[g2] - yc[g3]) * (xq - xc[g1]) + (xc[g3] - xc[g2]) * (yq - yc[g1])) / d;
    } else if id == 1 {
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

//* --------------------------------------------------------------------
//*  DAXPY - constant times a vector plus a vector
//* --------------------------------------------------------------------
#[allow(dead_code)]
pub fn daxpy_v(n: i32, da: f64, dx: &[f64], incx: i32, dy: &mut [f64], incy: i32) {
    if n <= 0 || da == 0.0 {
        return;
    }
    let (mut ix, mut iy);

    if incx != 1 || incy != 1 {
        ix = if incx >= 0 { 0 } else { (-n + 1) * incx };
        iy = if incy >= 0 { 0 } else { (-n + 1) * incy };
        for _ in 0..n {
            dy[iy as usize] += da * dx[ix as usize];
            ix += incx;
            iy += incy;
        }
    } else {
        let m = (n % 4) as usize;
        let n = n as usize;
        for i in 0..m {
            dy[i] += da * dx[i];
        }
        for i in (m..n).step_by(4) {
            dy[i] += da * dx[i];
            dy[i + 1] += da * dx[i + 1];
            dy[i + 2] += da * dx[i + 2];
            dy[i + 3] += da * dx[i + 3];
        }
    }
}
pub fn daxpy_m(
    n: i32,
    da: f64,
    abd: &mut [f64],
    stride: usize,
    col_x: usize,
    start_row_x: usize,
    incx: i32,
    col_y: usize,
    start_row_y: usize,
    incy: i32,
) {
    if n <= 0 || da == 0.0 {
        return;
    }
    for i in 0..(n as usize) {
        let idx_x = (start_row_x + (i as i32 * incx) as usize) * stride + col_x;
        let idx_y = (start_row_y + (i as i32 * incy) as usize) * stride + col_y;
        abd[idx_y] += da * abd[idx_x];
    }
}
pub fn daxpy(
    n: i32,
    da: f64,
    abd_x: &[f64],
    stride: usize,
    col_x: usize,
    start_row_x: usize,
    incx: i32,
    b_y: &mut [f64],
    start_idx_y: usize,
    incy: i32,
) {
    if n <= 0 || da == 0.0 {
        return;
    }
    for i in 0..(n as usize) {
        let idx_x = (start_row_x + (i as i32 * incx) as usize) * stride + col_x;
        let idx_y = start_idx_y + (i as i32 * incy) as usize;
        b_y[idx_y] += da * abd_x[idx_x];
    }
}
//
//* --------------------------------------------------------------------
//*  DDOT - dot product of two vectors
//* --------------------------------------------------------------------
pub fn ddot(
    n: i32,
    abd: &[f64],
    stride: usize,
    col: usize,
    row_start: usize,
    incx: i32,
    b: &[f64],
    start_b: usize,
    incy: i32,
) -> f64 {
    if n <= 0 {
        return 0.0;
    }
    let mut dtemp = 0.0;
    if incx != 1 || incy != 1 {
        let mut ix = row_start as i32;
        let mut iy = start_b as i32;
        for _ in 0..n {
            dtemp += abd[ix as usize * stride + col] * b[iy as usize];
            ix += incx;
            iy += incy;
        }
    } else {
        let m = (n % 5) as usize;
        let n = n as usize;
        for i in 0..m {
            dtemp += abd[(row_start + i) * stride + col] * b[start_b + i];
        }
        for i in (m..n).step_by(5) {
            dtemp += abd[(row_start + i) * stride + col] * b[start_b + i]
                + abd[(row_start + i + 1) * stride + col] * b[start_b + i + 1]
                + abd[(row_start + i + 2) * stride + col] * b[start_b + i + 2]
                + abd[(row_start + i + 3) * stride + col] * b[start_b + i + 3]
                + abd[(row_start + i + 4) * stride + col] * b[start_b + i + 4];
        }
    }
    dtemp
}
//
//* --------------------------------------------------------------------
//*  DGBFA - factor a real band matrix by elimination
//* --------------------------------------------------------------------
#[allow(unused_variables)]
pub fn dgbfa(
    abd: &mut [f64],
    _lda: usize,
    n: usize,
    ml: usize,
    mu: usize,
    ipvt: &mut [i32],
) -> i32 {
    let m = ml + mu + 1;
    let mut info: i32 = 0;
    let j0 = mu + 1;
    let j1 = n.min(m) - 1;

    for jz in j0..j1 {
        let i0 = m - jz;
        for i in i0..ml {
            abd[i * n + jz] = 0.0;
        }
    }

    let mut jz = j1;
    let mut ju = 0;

    for k in 0..n - 1 {
        jz += 1;
        if jz < n {
            for i in 0..ml {
                abd[i * n + jz] = 0.0;
            }
        }

        let lm = ml.min(n - k - 1);
        let l = idamax_m((lm + 1) as i32, abd, n, k, m - 1, 1) + (m - 1);
        ipvt[k] = (l + k - (m - 1)) as i32;

        if abd[l * n + k] == 0.0 {
            info = k as i32;
        } else {
            if l != m - 1 {
                let t = abd[l * n + k];
                abd[l * n + k] = abd[(m - 1) * n + k];
                abd[(m - 1) * n + k] = t;
            }

            let t = -1.0 / abd[(m - 1) * n + k];
            dscal_m(lm, t, abd, n, k, m, 1);

            ju = ju.max(mu + ipvt[k] as usize).min(n - 1);
            let mut mm = m - 1;
            let mut l_cur = l;

            for j in k + 1..ju + 1 {
                l_cur = l_cur.wrapping_sub(1);
                mm = mm.wrapping_sub(1);
                let t = abd[l_cur * n + j];
                if l_cur != mm {
                    abd[l_cur * n + j] = abd[mm * n + j];
                    abd[mm * n + j] = t;
                }
                daxpy_m(lm as i32, t, abd, n, k, m, 1, j, mm + 1, 1);
            }
        }
    }

    ipvt[n - 1] = (n - 1) as i32;

    if abd[(m - 1) * n + (n - 1)] == 0.0 {
        info = (n - 1) as i32;
    }

    info
}

//* --------------------------------------------------------------------
//*  DGBSL - solve a real banded system factored by DGBFA
//* --------------------------------------------------------------------
pub fn dgbsl(
    abd: &mut [f64],
    lda: usize,
    n: usize,
    ml: usize,
    mu: usize,
    ipvt: &[i32],
    b: &mut [f64],
    job: i32,
) {
    let m = mu + ml + 1;
    if job == 0 {
        if 0 < ml {
            for k in 0..n - 1 {
                let lm = ml.min(n - k - 1);
                let l = ipvt[k] as usize;
                let t = b[l];
                if l != k {
                    b[l] = b[k];
                    b[k] = t;
                }
                daxpy(lm as i32, t, abd, n, k, m, 1, b, k + 1, 1);
            }
        }

        for k in (0..n).rev() {
            b[k] /= abd[(m - 1) * n + k];
            let lm = k.min(m - 1);
            let la = m - 1 - lm;
            let lb = k - lm;
            let t = -b[k];
            daxpy(lm as i32, t, abd, n, k, la, 1, b, lb, 1);
        }
    } else {
        for k in 0..n {
            let lm = k.min(m - 1);
            let la = m - 1 - lm;
            let lb = k - lm;
            let t = ddot(lm as i32, abd, n, k, la, 1, b, lb, 1);
            b[k] = (b[k] - t) / abd[(m - 1) * n + k];
        }

        if 0 < ml {
            for k in (0..n - 1).rev() {
                let lm = ml.min(n - k - 1);
                b[k] += ddot(lm as i32, abd, n, k, m, 1, b, k + 1, 1);
                let l = ipvt[k] as usize;
                if l != k {
                    let t = b[l];
                    b[l] = b[k];
                    b[k] = t;
                }
            }
        }
    }
}
//
//* --------------------------------------------------------------------
//*  DSCAL - scale a vector by a constant
//* --------------------------------------------------------------------
pub fn dscal_v(n: i32, sa: f64, x: &mut [f64], incx: i32) {
    if n <= 0 {
        return;
    }
    if incx == 1 {
        let m = (n % 5) as usize;
        let n = n as usize;
        for i in 0..m {
            x[i] *= sa;
        }
        for i in (m..n).step_by(5) {
            x[i] *= sa;
            x[i + 1] *= sa;
            x[i + 2] *= sa;
            x[i + 3] *= sa;
            x[i + 4] *= sa;
        }
    } else {
        let mut ix = if incx >= 0 { 0 } else { (-n + 1) * incx };
        for _ in 0..n {
            x[ix as usize] *= sa;
            ix += incx;
        }
    }
}
pub fn dscal_m(
    n: usize,
    sa: f64,
    abd: &mut [f64],
    stride: usize,
    col: usize,
    start_row: usize,
    incx: i32,
) {
    if n == 0 {
        return;
    }
    for i in 0..n {
        abd[(start_row + (i as i32 * incx) as usize) * stride + col] *= sa;
    }
}
//* --------------------------------------------------------------------
//*  FILE_NAME_INC - increment a partially numeric filename
//* --------------------------------------------------------------------
pub fn file_name_inc(file_name: &str) -> String {
    if file_name.is_empty() {
        println!("\nFILE_NAME_INC - Fatal error!");
        println!("The input string is empty.");
        std::process::exit(1);
    }
    let mut chars: Vec<char> = file_name.chars().collect();
    let mut change = 0;
    for i in (0..chars.len()).rev() {
        if chars[i].is_ascii_digit() {
            change += 1;
            let mut digit = chars[i] as u8 - b'0';
            digit += 1;
            if digit == 10 {
                chars[i] = '0';
            } else {
                chars[i] = (digit + b'0') as char;
                return chars.into_iter().collect();
            }
        }
    }
    if change == 0 {
        String::new()
    } else {
        chars.into_iter().collect()
    }
}
//
//* --------------------------------------------------------------------
//*  GETG - extract values of a quantity along the profile line
//* --------------------------------------------------------------------
pub fn getg(f: &[f64], iline: &[i32], my: usize) -> Vec<f64> {
    let mut u = vec![0.0; my];
    for i in 0..my {
        let j = iline[i];
        u[i] = if j <= 0 { 0.0 } else { f[(j - 1) as usize] };
    }
    u
}
//
//* --------------------------------------------------------------------
//*  GRAM - compute the Gram matrix and R vector
//* --------------------------------------------------------------------
pub fn gram(bump: &mut bump_struct::Bump) {
    let wt: [f64; 3] = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
    let yq_gauss: [f64; 3] = [-0.7745966692, 0.0, 0.7745966692];

    for it in 0..bump.my {
        bump.r[it] = 0.0;
    }
    for it in 0..bump.my {
        for j in 0..bump.my {
            bump.gr[it * bump.my + j] = 0.0;
        }
    }

    for it in 0..bump.nelemn {
        let k = bump.node[it * NNODES];
        let kk = bump.node[it * NNODES + 1];

        if (bump.xc[k] - bump.xprof).abs() > 1.0e-4 || (bump.xc[kk] - bump.xprof).abs() > 1.0e-4 {
            continue;
        }

        for iquad in 0..3 {
            let bma2 = (bump.yc[kk] - bump.yc[k]) / 2.0;
            let ar = bma2 * wt[iquad];
            let x = bump.xprof;
            let y = bump.yc[k] + bma2 * (yq_gauss[iquad] + 1.0);

            let mut uiqdpt = 0.0;
            for iq in 0..NNODES {
                if iq == 0 || iq == 1 || iq == 3 {
                    let (bb, _bx, _by) = qbf(x, y, it, iq, &bump.node, &bump.xc, &bump.yc);
                    let ip = bump.node[it * NNODES + iq];
                    let iun = bump.indx[ip * 2];
                    if 0 < iun {
                        let ii = bump.iline_inv[iun as usize];
                        uiqdpt += bb * bump.uprof[(ii - 1) as usize];
                    } else if iun == -1 {
                        let ubc = ubdry(1, bump.yc[ip]);
                        uiqdpt += bb * ubc;
                    }
                }
            }

            for iq in 0..NNODES {
                if iq == 0 || iq == 1 || iq == 3 {
                    let ip = bump.node[it * NNODES + iq];
                    let (bb, _bx, _by) = qbf(x, y, it, iq, &bump.node, &bump.xc, &bump.yc);
                    let i_val = bump.indx[ip * 2];
                    if 0 < i_val {
                        let ii = bump.iline_inv[i_val as usize];
                        bump.r[(ii - 1) as usize] += bb * uiqdpt * ar;
                        for iqq in 0..NNODES {
                            if iqq == 0 || iqq == 1 || iqq == 3 {
                                let ipp = bump.node[it * NNODES + iqq];
                                let (bbb, _bbx_, _bby_) =
                                    qbf(x, y, it, iqq, &bump.node, &bump.xc, &bump.yc);
                                let j_val = bump.indx[ipp * 2];
                                if j_val != 0 {
                                    let jj = bump.iline_inv[j_val as usize];
                                    bump.gr[((ii - 1) as usize) * bump.my + ((jj - 1) as usize)] +=
                                        bb * bbb * ar;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if 3 <= bump.iwrite {
        println!();
        println!("Gram matrix:");
        println!();
        for i in 0..bump.my {
            for j in 0..bump.my {
                println!("{} {} {}", i + 1, j + 1, bump.gr[i * bump.my + j]);
            }
        }
        println!();
        println!("R vector:");
        println!();
        for i in 0..bump.my {
            println!("{}", bump.r[i]);
        }
    }
}

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

//* --------------------------------------------------------------------
//*  IDAMAX - index of the vector element of maximum absolute value
//* --------------------------------------------------------------------
pub fn idamax_v(n: i32, dx: &[f64], incx: i32) -> usize {
    if n < 1 || incx <= 0 {
        return 0;
    }
    if n == 1 {
        return 0;
    }
    if incx == 1 {
        let mut dmax = dx[0].abs();
        let mut idamax_val = 0;
        let n = n as usize;
        for i in 1..n {
            let val = dx[i].abs();
            if dmax < val {
                idamax_val = i;
                dmax = val;
            }
        }
        idamax_val
    } else {
        let mut ix = 0i32;
        let mut dmax = dx[0].abs();
        let mut idamax_val = 0;
        ix += incx;
        for i in 1..n {
            let val = dx[ix as usize].abs();
            if dmax < val {
                idamax_val = i;
                dmax = val;
            }
            ix += incx;
        }
        idamax_val as usize
    }
}

pub fn idamax_m(
    n: i32,
    abd: &[f64],
    stride: usize,
    col: usize,
    start_row: usize,
    incx: i32,
) -> usize {
    if n < 1 || incx <= 0 {
        return 0;
    }
    let mut dmax = abd[start_row * stride + col].abs();
    let mut idamax_val: usize = 0;
    for i in 1..(n as usize) {
        let current_val = abd[(start_row + (i as i32 * incx) as usize) * stride + col].abs();
        if current_val > dmax {
            idamax_val = i;
            dmax = current_val;
        }
    }
    idamax_val
}
//
//* --------------------------------------------------------------------
//*  IGETL - get the local unknown number along the profile line
//* --------------------------------------------------------------------
pub fn igetl(i: i32, iline: &[i32]) -> i32 {
    for (j, &val) in iline.iter().enumerate() {
        if val == i {
            return (j + 1) as i32;
        }
    }
    -1
}
//
//* --------------------------------------------------------------------
//*  LINSYS - solve the linearized Navier Stokes equation
//* --------------------------------------------------------------------
pub fn linsys(bump: &mut bump_struct::Bump, itype: i32) {
    let ioff = bump.nlband + bump.nlband + 1;
    let visc = 1.0 / bump.reynld;
    let use_sens = itype == -2;
    let neqn = bump.neqn;
    let nrow = bump.nrow;

    let mut rhs = vec![0.0; neqn];
    bump.a[..nrow * neqn].fill(0.0);

    let nelemn = bump.nelemn;
    let node_v = &bump.node;
    let indx_v = &bump.indx;
    let insc_v = &bump.insc;
    let xc_v = &bump.xc;
    let yc_v = &bump.yc;
    let phi_v = &bump.phi;
    let psi_v = &bump.psi;
    let area_v = &bump.area;
    let xm_v = &bump.xm;
    let ym_v = &bump.ym;
    let isotri_v = &bump.isotri;

    for it in 0..nelemn {
        let elem_nodes = [
            node_v[it * NNODES],
            node_v[it * NNODES + 1],
            node_v[it * NNODES + 2],
            node_v[it * NNODES + 3],
            node_v[it * NNODES + 4],
            node_v[it * NNODES + 5],
        ];
        let elem_indx_u = [
            indx_v[elem_nodes[0] * 2],
            indx_v[elem_nodes[1] * 2],
            indx_v[elem_nodes[2] * 2],
            indx_v[elem_nodes[3] * 2],
            indx_v[elem_nodes[4] * 2],
            indx_v[elem_nodes[5] * 2],
        ];
        let elem_indx_v = [
            indx_v[elem_nodes[0] * 2 + 1],
            indx_v[elem_nodes[1] * 2 + 1],
            indx_v[elem_nodes[2] * 2 + 1],
            indx_v[elem_nodes[3] * 2 + 1],
            indx_v[elem_nodes[4] * 2 + 1],
            indx_v[elem_nodes[5] * 2 + 1],
        ];
        let elem_insc = [
            insc_v[elem_nodes[0]],
            insc_v[elem_nodes[1]],
            insc_v[elem_nodes[2]],
            insc_v[elem_nodes[3]],
            insc_v[elem_nodes[4]],
            insc_v[elem_nodes[5]],
        ];
        let isotri = isotri_v[it];
        let area = area_v[it];
        let mut ar = area / 3.0;

        for iquad in 0..NQUAD {
            let yq = ym_v[it * NQUAD + iquad];
            let xq = xm_v[it * NQUAD + iquad];

            let (det, etax, etay, xix, xiy);

            if isotri == 1 {
                (det, etax, etay, xix, xiy) = trans(it, node_v, xc_v, xq, yc_v, yq);
                ar = det * area / 3.0;
            } else {
                (det, etax, etay, xix, xiy) = (0.0, 0.0, 0.0, 0.0, 0.0);
            }

            let (un, unx, uny) = uval(etax, etay, bump, it, xix, xiy, xq, yq);

            let phi_bb: [f64; NNODES] = [
                phi_v[phi_idx(it, iquad, 0, 0)],
                phi_v[phi_idx(it, iquad, 1, 0)],
                phi_v[phi_idx(it, iquad, 2, 0)],
                phi_v[phi_idx(it, iquad, 3, 0)],
                phi_v[phi_idx(it, iquad, 4, 0)],
                phi_v[phi_idx(it, iquad, 5, 0)],
            ];
            let phi_bx: [f64; NNODES] = [
                phi_v[phi_idx(it, iquad, 0, 1)],
                phi_v[phi_idx(it, iquad, 1, 1)],
                phi_v[phi_idx(it, iquad, 2, 1)],
                phi_v[phi_idx(it, iquad, 3, 1)],
                phi_v[phi_idx(it, iquad, 4, 1)],
                phi_v[phi_idx(it, iquad, 5, 1)],
            ];
            let phi_by: [f64; NNODES] = [
                phi_v[phi_idx(it, iquad, 0, 2)],
                phi_v[phi_idx(it, iquad, 1, 2)],
                phi_v[phi_idx(it, iquad, 2, 2)],
                phi_v[phi_idx(it, iquad, 3, 2)],
                phi_v[phi_idx(it, iquad, 4, 2)],
                phi_v[phi_idx(it, iquad, 5, 2)],
            ];
            let psi_bbl: [f64; NNODES] = [
                psi_v[psi_idx(it, iquad, 0)],
                psi_v[psi_idx(it, iquad, 1)],
                psi_v[psi_idx(it, iquad, 2)],
                psi_v[psi_idx(it, iquad, 3)],
                psi_v[psi_idx(it, iquad, 4)],
                psi_v[psi_idx(it, iquad, 5)],
            ];

            for iq in 0..NNODES {
                let bb = phi_bb[iq];
                let bx = phi_bx[iq];
                let by = phi_by[iq];
                let bbl = psi_bbl[iq];
                let ihor = elem_indx_u[iq];
                let iver = elem_indx_v[iq];
                let iprs = elem_insc[iq];

                if 0 < ihor {
                    rhs[(ihor - 1) as usize] += ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
                }
                if 0 < iver {
                    rhs[(iver - 1) as usize] += ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
                }

                for iqq in 0..NNODES {
                    let bbb = phi_bb[iqq];
                    let bbx = phi_bx[iqq];
                    let bby = phi_by[iqq];
                    let bbbl = psi_bbl[iqq];
                    let ju = elem_indx_u[iqq];
                    let jv = elem_indx_v[iqq];
                    let jp = elem_insc[iqq];
                    let ipp = elem_nodes[iqq];

                    if 0 < ju {
                        if 0 < ihor {
                            let iuse = (ihor - ju + ioff as i32) as usize;
                            let row_off = (iuse - 1) * neqn;
                            let col = (ju - 1) as usize;
                            bump.a[row_off + col] += ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
                        }
                        if 0 < iver {
                            let iuse = (iver - ju + ioff as i32) as usize;
                            let row_off = (iuse - 1) * neqn;
                            let col = (ju - 1) as usize;
                            bump.a[row_off + col] += ar * bb * bbb * unx[1];
                        }
                        if 0 < iprs {
                            let iuse = (iprs - ju + ioff as i32) as usize;
                            let row_off = (iuse - 1) * neqn;
                            let col = (ju - 1) as usize;
                            bump.a[row_off + col] += ar * bbx * bbl;
                        }
                    } else if ju == itype {
                        let ubc = if ju == -1 {
                            ubdry(1, yc_v[ipp])
                        } else {
                            ubump(bump, ipp, iqq, it, 1)
                        };

                        if 0 < ihor {
                            let aij = ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
                            rhs[(ihor - 1) as usize] -= ubc * aij;
                        }
                        if 0 < iver {
                            let aij = ar * bb * bbb * unx[1];
                            rhs[(iver - 1) as usize] -= ubc * aij;
                        }
                        if 0 < iprs {
                            let aij = ar * bbx * bbl;
                            rhs[(iprs - 1) as usize] -= ubc * aij;
                        }
                    }

                    if 0 < jv {
                        if 0 < ihor {
                            let iuse = (ihor - jv + ioff as i32) as usize;
                            let row_off = (iuse - 1) * neqn;
                            let col = (jv - 1) as usize;
                            bump.a[row_off + col] += ar * bb * bbb * uny[0];
                        }
                        if 0 < iver {
                            let iuse = (iver - jv + ioff as i32) as usize;
                            let row_off = (iuse - 1) * neqn;
                            let col = (jv - 1) as usize;
                            bump.a[row_off + col] += ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
                        }
                        if 0 < iprs {
                            let iuse = (iprs - jv + ioff as i32) as usize;
                            let row_off = (iuse - 1) * neqn;
                            let col = (jv - 1) as usize;
                            bump.a[row_off + col] += ar * bby * bbl;
                        }
                    } else if jv == itype {
                        let ubc = if jv == -1 {
                            ubdry(2, yc_v[ipp])
                        } else {
                            ubump(bump, ipp, iqq, it, 2)
                        };

                        if 0 < ihor {
                            let aij = ar * bb * bbb * uny[0];
                            rhs[(ihor - 1) as usize] -= ubc * aij;
                        }
                        if 0 < iver {
                            let aij = ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
                            rhs[(iver - 1) as usize] -= ubc * aij;
                        }
                        if 0 < iprs {
                            let aij = ar * bby * bbl;
                            rhs[(iprs - 1) as usize] -= ubc * aij;
                        }
                    }

                    if 0 < jp {
                        if 0 < ihor {
                            let iuse = (ihor - jp + ioff as i32) as usize;
                            let row_off = (iuse - 1) * neqn;
                            let col = (jp - 1) as usize;
                            bump.a[row_off + col] -= ar * bx * bbbl;
                        }
                        if 0 < iver {
                            let iuse = (iver - jp + ioff as i32) as usize;
                            let row_off = (iuse - 1) * neqn;
                            let col = (jp - 1) as usize;
                            bump.a[row_off + col] -= ar * by * bbbl;
                        }
                    }
                }
            }
        }
    }

    rhs[neqn - 1] = 0.0;
    for j_1based in (neqn - bump.nlband)..neqn {
        let j = j_1based - 1;
        let i = neqn - j_1based + ioff;
        bump.a[(i - 1) * neqn + j] = 0.0;
    }
    bump.a[(ioff - 1) * neqn + (neqn - 1)] = 1.0;

    let mut ipvt = vec![0i32; neqn];
    let info = dgbfa(
        &mut bump.a,
        bump.maxrow,
        neqn,
        bump.nlband,
        bump.nlband,
        &mut ipvt,
    );

    if info != 0 {
        println!();
        println!("LINSYS - fatal error!");
        println!("DGBFA returns INFO = {}", info + 1);
        std::process::exit(1);
    }

    let actual_rhs: &mut [f64] = if use_sens {
        &mut bump.sens
    } else {
        &mut bump.f
    };
    actual_rhs.copy_from_slice(&rhs);
    dgbsl(
        &mut bump.a,
        bump.maxrow,
        neqn,
        bump.nlband,
        bump.nlband,
        &ipvt,
        actual_rhs,
        0,
    );
}

//* --------------------------------------------------------------------
//*  NSTOKE - solve Navier Stokes using Newton iteration
//* --------------------------------------------------------------------
pub fn nstoke(bump: &mut bump_struct::Bump) {
    for iter in 0..bump.maxnew {
        bump.numnew += 1;

        linsys(bump, -1);

        for i in 0..bump.neqn {
            bump.g[i] -= bump.f[i];
        }

        let imax = idamax_v(bump.neqn as i32, &bump.g, 1);
        let diff = bump.g[imax].abs();
        println!("NSTOKE: Iteration {}, MaxNorm(diff) = {}", iter + 1, diff);

        for i in 0..bump.neqn {
            bump.g[i] = bump.f[i];
        }

        if diff <= bump.tolnew {
            println!("NSTOKE converged.");
            return;
        }

        if iter == bump.maxnew - 1 {
            println!("NSTOKE failed!");
            std::process::exit(1);
        }
    }
}

//* --------------------------------------------------------------------
//*  QBF - evaluate quadratic basis functions
//* --------------------------------------------------------------------
pub fn qbf(
    xq: f64,
    yq: f64,
    it: usize,
    inn: usize,
    node: &[usize],
    xc: &[f64],
    yc: &[f64],
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
        i1 = node[it * NNODES + in1];
        i2 = node[it * NNODES + in2];
        i3 = node[it * NNODES + in3];
        d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
        t = 1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d;
        bb = t * (2.0 * t - 1.0);
        bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d;
        by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d;
    } else {
        inn_local = inn - 3;
        in1 = inn_local;
        in2 = (inn_local + 1) % 3;
        in3 = (inn_local + 2) % 3;
        i1 = node[it * NNODES + in1];
        i2 = node[it * NNODES + in2];
        i3 = node[it * NNODES + in3];
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
pub fn refbsp(xq: f64, yq: f64, iq: usize) -> f64 {
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
    x: f64,
    y: f64,
    inn: usize,
    etax: f64,
    etay: f64,
    xix: f64,
    xiy: f64,
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

    return (bb, bx, by);
}

//* --------------------------------------------------------------------
//*  RESID - compute the residual
//* --------------------------------------------------------------------
pub fn resid(bump: &mut bump_struct::Bump) {
    let visc = 1.0 / bump.reynld;

    for it in 0..bump.neqn {
        bump.res[it] = 0.0;
    }

    for it in 0..bump.nelemn {
        let mut ar = bump.area[it] / 3.0;

        for iquad in 0..NQUAD {
            let yq = bump.ym[it * NQUAD + iquad];
            let xq = bump.xm[it * NQUAD + iquad];

            let (det, etax, etay, xix, xiy) = if bump.isotri[it] == 1 {
                let result = trans(it, &bump.node, &bump.xc, xq, &bump.yc, yq);
                ar = result.0 * bump.area[it] / 3.0;
                result
            } else {
                (0.0, 0.0, 0.0, 0.0, 0.0)
            };

            let (un, unx, uny) = uval(etax, etay, bump, it, xix, xiy, xq, yq);

            for iq in 0..NNODES {
                let ip = bump.node[it * NNODES + iq];
                let bb = bump.phi[phi_idx(it, iquad, iq, 0)];
                let bx = bump.phi[phi_idx(it, iquad, iq, 1)];
                let by = bump.phi[phi_idx(it, iquad, iq, 2)];
                let bbl = bump.psi[psi_idx(it, iquad, iq)];
                let ihor = bump.indx[ip * 2];
                let iver = bump.indx[ip * 2 + 1];
                let iprs = bump.insc[ip];

                if 0 < ihor {
                    bump.res[(ihor - 1) as usize] -= ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
                }
                if 0 < iver {
                    bump.res[(iver - 1) as usize] -= ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
                }

                for iqq in 0..NNODES {
                    let ipp = bump.node[it * NNODES + iqq];
                    let bbb = bump.phi[phi_idx(it, iquad, iqq, 0)];
                    let bbx = bump.phi[phi_idx(it, iquad, iqq, 1)];
                    let bby = bump.phi[phi_idx(it, iquad, iqq, 2)];
                    let bbbl = bump.psi[psi_idx(it, iquad, iqq)];
                    let ju = bump.indx[ipp * 2];
                    let jv = bump.indx[ipp * 2 + 1];
                    let jp = bump.insc[ipp];

                    if 0 < ju {
                        if 0 < ihor {
                            bump.res[(ihor - 1) as usize] += ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]))
                                * bump.g[(ju - 1) as usize];
                        }
                        if 0 < iver {
                            bump.res[(iver - 1) as usize] +=
                                ar * bb * bbb * unx[1] * bump.g[(ju - 1) as usize];
                        }
                        if 0 < iprs {
                            bump.res[(iprs - 1) as usize] +=
                                ar * bbx * bbl * bump.g[(ju - 1) as usize];
                        }
                    } else if ju == -1 {
                        let ubc = ubdry(1, bump.yc[ipp]);
                        if 0 < ihor {
                            let aij = ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
                            bump.res[(ihor - 1) as usize] += ubc * aij;
                        }
                        if 0 < iver {
                            let aij = ar * bb * bbb * unx[1];
                            bump.res[(iver - 1) as usize] += ubc * aij;
                        }
                        if 0 < iprs {
                            let aij = ar * bbx * bbl;
                            bump.res[(iprs - 1) as usize] += ubc * aij;
                        }
                    } else if ju == -2 {
                        let ubc = ubump(bump, ipp, iqq, it, 1);
                        if 0 < ihor {
                            let aij = ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
                            bump.res[(ihor - 1) as usize] += ubc * aij;
                        }
                        if 0 < iver {
                            let aij = ar * bb * bbb * unx[1];
                            bump.res[(iver - 1) as usize] += ubc * aij;
                        }
                        if 0 < iprs {
                            let aij = ar * bbx * bbl;
                            bump.res[(iprs - 1) as usize] += ubc * aij;
                        }
                    }

                    if 0 < jv {
                        if 0 < ihor {
                            bump.res[(ihor - 1) as usize] +=
                                ar * bb * bbb * uny[0] * bump.g[(jv - 1) as usize];
                        }
                        if 0 < iver {
                            bump.res[(iver - 1) as usize] += ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]))
                                * bump.g[(jv - 1) as usize];
                        }
                        if 0 < iprs {
                            bump.res[(iprs - 1) as usize] +=
                                ar * bby * bbl * bump.g[(jv - 1) as usize];
                        }
                    } else if jv == -1 {
                        let ubc = ubdry(2, bump.yc[ipp]);
                        if 0 < ihor {
                            bump.res[(ihor - 1) as usize] += ubc * ar * bb * bbb * uny[0];
                        }
                        if 0 < iver {
                            let aij = ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
                            bump.res[(iver - 1) as usize] += ubc * aij;
                        }
                        if 0 < iprs {
                            bump.res[(iprs - 1) as usize] += ubc * ar * bby * bbl;
                        }
                    } else if jv == -2 {
                        let ubc = ubump(bump, ipp, iqq, it, 2);
                        if 0 < ihor {
                            bump.res[(ihor - 1) as usize] += ubc * ar * bb * bbb * uny[0];
                        }
                        if 0 < iver {
                            let aij = ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
                            bump.res[(iver - 1) as usize] += ubc * aij;
                        }
                        if 0 < iprs {
                            bump.res[(iprs - 1) as usize] += ubc * ar * bby * bbl;
                        }
                    }

                    if 0 < jp {
                        if 0 < ihor {
                            bump.res[(ihor - 1) as usize] -=
                                ar * bx * bbbl * bump.g[(jp - 1) as usize];
                        }
                        if 0 < iver {
                            bump.res[(iver - 1) as usize] -=
                                ar * by * bbbl * bump.g[(jp - 1) as usize];
                        }
                    }
                }
            }
        }
    }

    bump.res[bump.neqn - 1] = bump.g[bump.neqn - 1];

    let mut rmax = 0.0;
    let mut imax = 0;
    let mut ibad = 0;

    for i in 0..bump.neqn {
        let test = bump.res[i].abs();
        if rmax < test {
            rmax = test;
            imax = i;
        }
        if 1.0e-3 < test {
            ibad += 1;
        }
    }

    if 1 <= bump.iwrite {
        println!();
        println!("RESIDUAL INFORMATION:");
        println!();
        println!("Worst residual is number {}", imax + 1);
        println!("of magnitude {}", rmax);
        println!();
        println!(
            "Number of \"bad\" residuals is {} out of {}",
            ibad, bump.neqn
        );
        println!();
    }

    if 2 <= bump.iwrite {
        println!("Raw residuals:");
        println!();
        let mut idx = 0;
        for j in 0..bump.n_points {
            if 0 < bump.indx[j * 2] {
                if bump.res[idx].abs() <= 1.0e-3 {
                    println!(" U {} {} {}", idx + 1, j + 1, bump.res[idx]);
                } else {
                    println!("*U {} {} {}", idx + 1, j + 1, bump.res[idx]);
                }
                idx += 1;
            }
            if 0 < bump.indx[j * 2 + 1] {
                if bump.res[idx].abs() <= 1.0e-3 {
                    println!(" V {} {} {}", idx + 1, j + 1, bump.res[idx]);
                } else {
                    println!("*V {} {} {}", idx + 1, j + 1, bump.res[idx]);
                }
                idx += 1;
            }
            if 0 < bump.insc[j] {
                if bump.res[idx].abs() <= 1.0e-3 {
                    println!(" P {} {} {}", idx + 1, j + 1, bump.res[idx]);
                } else {
                    println!("*P {} {} {}", idx + 1, j + 1, bump.res[idx]);
                }
                idx += 1;
            }
        }
    }
}
//
//
//* --------------------------------------------------------------------
//*  SETBAN - compute the half band width
//* --------------------------------------------------------------------
pub fn setban(bump: &mut bump_struct::Bump) {
    bump.nlband = 0;

    for it in 0..bump.nelemn {
        for iq in 0..NNODES {
            let ip = bump.node[it * NNODES + iq];
            for iuk in 0..3 {
                let i_val = if iuk == 2 {
                    bump.insc[ip]
                } else {
                    bump.indx[ip * 2 + iuk]
                };
                if 0 < i_val {
                    for iqq in 0..NNODES {
                        let ipp = bump.node[it * NNODES + iqq];
                        for iukk in 0..3 {
                            let j_val = if iukk == 2 {
                                bump.insc[ipp]
                            } else {
                                bump.indx[ipp * 2 + iukk]
                            };
                            if 0 < j_val && j_val > i_val {
                                let diff = (j_val - i_val) as usize;
                                if bump.nlband < diff {
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

    println!();
    println!("SETBAN:");
    println!();
    println!("  Lower bandwidth = {}", bump.nlband);
    println!("  Total bandwidth = {}", bump.nband);
    println!("  Required matrix rows = {}", bump.nrow);

    if bump.maxrow < bump.nrow {
        println!("SETBAN - NROW is too large!");
        println!("The maximum allowed is {}", bump.maxrow);
        println!("This problem requires NROW = {}", bump.nrow);
        std::process::exit(1);
    }
}

//* --------------------------------------------------------------------
//*  SETBAS - evaluate basis functions at each integration point
//* --------------------------------------------------------------------
pub fn setbas(bump: &mut bump_struct::Bump) {
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
        for j in 0..NQUAD {
            xq = bump.xm[it * NQUAD + j];
            yq = bump.ym[it * NQUAD + j];
            (det, etax, etay, xix, xiy) = trans(it, &bump.node, &bump.xc, xq, &bump.yc, yq);

            for iq in 0..NNODES {
                if bump.isotri[it] == 0 {
                    bump.psi[psi_idx(it, j, iq)] =
                        bsp(it, iq, 0, &bump.node, &bump.xc, xq, &bump.yc, yq);
                    (bb, bx, by) = qbf(xq, yq, it, iq, &bump.node, &bump.xc, &bump.yc);
                } else {
                    (bb, bx, by) = refqbf(xq, yq, iq, etax, etay, xix, xiy);
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
pub fn setgrd(ibump: i32, bump: &mut bump_struct::Bump) /*  -> (bool, usize) */
{
    println!();
    println!("SETGRD:");
    println!();

    if bump.ny < bump.nx {
        bump._long = true;
        println!("Using vertical ordering.");
    } else {
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
    println!(
        "Bump extends from {} at bump.node {}",
        bump.xbleft,
        nbleft + 1
    );
    println!(
        "               to {} at bump.node {}",
        bump.xbrite,
        nbrite + 1
    );

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

                bump.node[ielemn * NNODES] = ip;
                bump.node[ielemn * NNODES + 1] = ip + 2;
                bump.node[ielemn * NNODES + 2] = ip2 + 2;
                bump.node[ielemn * NNODES + 3] = ip + 1;
                bump.node[ielemn * NNODES + 4] = ip1 + 2;
                bump.node[ielemn * NNODES + 5] = ip1 + 1;

                if ibump == 0 {
                    bump.isotri[ielemn] = 0;
                } else if ibump == 1 {
                    bump.isotri[ielemn] = 0;
                } else if ibump == 2 {
                    bump.isotri[ielemn] = if nbleft <= ic && ic < nbrite { 1 } else { 0 };
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

                if ibump == 0 {
                    bump.isotri[ielemn] = 0;
                } else if ibump == 1 {
                    bump.isotri[ielemn] = if jc == 0 && nbleft <= ic && ic < nbrite {
                        1
                    } else {
                        0
                    };
                } else if ibump == 2 {
                    bump.isotri[ielemn] = if nbleft <= ic && ic < nbrite { 1 } else { 0 };
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

                if ibump == 0 {
                    bump.isotri[ielemn] = 0;
                } else if ibump == 1 {
                    bump.isotri[ielemn] = 0;
                } else if ibump == 2 {
                    bump.isotri[ielemn] = if nbleft <= ic && ic < nbrite { 1 } else { 0 };
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

                if ibump == 0 {
                    bump.isotri[ielemn] = 0;
                } else if ibump == 1 {
                    bump.isotri[ielemn] = if jc == 0 && nbleft <= ic && ic < nbrite {
                        1
                    } else {
                        0
                    };
                } else if ibump == 2 {
                    bump.isotri[ielemn] = if nbleft <= ic && ic < nbrite { 1 } else { 0 };
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
                i + 1,
                bump.indx[i * 2],
                bump.indx[i * 2 + 1],
                bump.insc[i]
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
            print!("{}\t", it + 1);
            for i in 0..6 {
                print!("{}\t", bump.node[it * NNODES + i] + 1);
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

//* --------------------------------------------------------------------
//*  SETLIN - determine unknown numbers along the profile line
//* --------------------------------------------------------------------
pub fn setlin(bump: &mut bump_struct::Bump) {
    let itemp = (2.0 * (bump.nx as f64 - 1.0) * bump.xprof / bump.xlngth).round() as usize;

    let nodex0 = if bump._long {
        itemp * (2 * bump.ny - 1)
    } else {
        itemp
    };

    println!();
    println!("SETLIN:");
    println!();
    println!("  Profile generated at X = {}", bump.xprof);
    println!("  which is above node  = {}", nodex0 + 1);

    for i in 0..bump.my {
        let ip = if bump._long {
            nodex0 + i
        } else {
            nodex0 + bump.mx * i
        };
        bump.iline[i] = bump.indx[ip * 2];
    }

    bump.iline_inv = vec![-1; bump.neqn + 1];
    for (j, &val) in bump.iline.iter().enumerate() {
        if val > 0 {
            bump.iline_inv[val as usize] = (j + 1) as i32;
        }
    }

    if 1 <= bump.iwrite {
        println!();
        println!("  Indices of unknowns along the profile line:");
        println!();
        for i in (0..bump.my).step_by(5) {
            for j in i..((i + 5).min(bump.my)) {
                print!("{:5}", bump.iline[j]);
            }
            println!();
        }
    }
}

//* --------------------------------------------------------------------
//*  SETQUD - set midpoint quadrature rule information
//* --------------------------------------------------------------------
pub fn setqud(bump: &mut bump_struct::Bump) {
    for it in 0..bump.nelemn {
        let ip1: usize = bump.node[it * NNODES];
        let ip2: usize = bump.node[it * NNODES + 1];
        let ip3: usize = bump.node[it * NNODES + 2];
        let x1: f64 = bump.xc[ip1];
        let x2: f64 = bump.xc[ip2];
        let x3: f64 = bump.xc[ip3];
        let y1: f64 = bump.yc[ip1];
        let y2: f64 = bump.yc[ip2];
        let y3: f64 = bump.yc[ip3];

        if bump.isotri[it] == 0 {
            bump.xm[it * NQUAD] = 0.5 * (x1 + x2);
            bump.xm[it * NQUAD + 1] = 0.5 * (x2 + x3);
            bump.xm[it * NQUAD + 2] = 0.5 * (x3 + x1);
            bump.ym[it * NQUAD] = 0.5 * (y1 + y2);
            bump.ym[it * NQUAD + 1] = 0.5 * (y2 + y3);
            bump.ym[it * NQUAD + 2] = 0.5 * (y3 + y1);
            bump.area[it] =
                0.5 * ((y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3)).abs();
        } else {
            bump.xm[it * NQUAD] = 0.5;
            bump.ym[it * NQUAD] = 0.5;
            bump.xm[it * NQUAD + 1] = 1.0;
            bump.ym[it * NQUAD + 1] = 0.5;
            bump.xm[it * NQUAD + 2] = 0.5;
            bump.ym[it * NQUAD + 2] = 0.0;
            bump.area[it] = 0.5;
        }
    }
    if 3 <= bump.iwrite {
        println!();
        println!("SETQUD: Element Areas && Quadrature points:");
        println!();
        for i in 0..bump.nelemn {
            println!("{}\t{}", i + 1, bump.area[i]);
            for j in 0..NQUAD {
                println!(
                    "{}\t{}\t{}\t{}",
                    i + 1,
                    j + 1,
                    bump.xm[i * NQUAD + j],
                    bump.ym[i * NQUAD + j]
                );
            }
        }
    }
}

//* --------------------------------------------------------------------
//*  SETXY - set the grid coordinates based on the parameter value
//* --------------------------------------------------------------------
pub fn setxy(bump: &mut bump_struct::Bump) {
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

        bump.yc[ip] = ((bump.my - 1 - jc) as f64 * ylo + (jc as f64) * bump.ylngth)
            / ((bump.ny + bump.ny) as f64 - 2.0);
    }
    if 2 <= bump.iwrite {
        println!();
        println!("SETXY:");
        println!();
        println!("     I     XC     YC");
        println!();
        for i in 0..bump.n_points {
            println!("{}\t{}\t{}", i + 1, bump.xc[i], bump.yc[i]);
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
    node: &[usize],
    xc: &[f64],
    xq: f64,
    yc: &[f64],
    yq: f64,
) -> (f64, f64, f64, f64, f64) {
    let i1 = node[it * NNODES];
    let i2 = node[it * NNODES + 1];
    let i3 = node[it * NNODES + 2];
    let i4 = node[it * NNODES + 3];
    let i5 = node[it * NNODES + 4];
    let i6 = node[it * NNODES + 5];

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

    let mut det = (2.0 * a1 * b2 - 2.0 * a2 * b1) * xq * xq
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

//* --------------------------------------------------------------------
//*  UBDRY - parabolic inflow boundary condition
//* --------------------------------------------------------------------
pub fn ubdry(iuk: usize, yy: f64) -> f64 {
    if iuk == 1 {
        (-2.0 * yy + 6.0) * yy / 9.0
    } else {
        0.0
    }
}
//
//* --------------------------------------------------------------------
//*  UBUMP - sensitivity dU/dA on the bump
//* --------------------------------------------------------------------
pub fn ubump(bump: &bump_struct::Bump, ip: usize, iqq: usize, it: usize, iukk: usize) -> f64 {
    let mut det = 1.0;
    let mut etax = 0.0;
    let mut etay = 0.0;
    let mut xix = 1.0;
    let mut xiy = 1.0;
    let (xq, yq);

    if bump.isotri[it] == 0 {
        xq = bump.xc[ip];
        yq = bump.yc[ip];
    } else {
        (xq, yq) = match iqq {
            0 => (0.0, 0.0),
            1 => (1.0, 1.0),
            2 => (1.0, 0.0),
            3 => (0.5, 0.5),
            4 => (1.0, 0.5),
            5 => (0.5, 0.0),
            _ => (0.0, 0.0),
        };
        (det, etax, etay, xix, xiy) = trans(it, &bump.node, &bump.xc, xq, &bump.yc, yq);
    }

    let (_un, _unx, uny) = _ubump_uval(bump, it, xix, xiy, xq, yq, det, etax, etay);

    if iukk == 1 {
        -uny[0] * (bump.xc[ip] - 1.0) * (bump.xc[ip] - 3.0)
    } else if iukk == 2 {
        -uny[1] * (bump.xc[ip] - 1.0) * (bump.xc[ip] - 3.0)
    } else {
        println!("UBUMP called for iukk = {}", iukk);
        std::process::exit(1);
    }
}

pub fn _ubump_uval(
    bump: &bump_struct::Bump,
    it: usize,
    xix: f64,
    xiy: f64,
    xq: f64,
    yq: f64,
    _det: f64,
    etax: f64,
    etay: f64,
) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let mut un = [0.0; 2];
    let mut unx = [0.0; 2];
    let mut uny = [0.0; 2];

    for iq in 0..NNODES {
        let (bb, bx, by) = if bump.isotri[it] == 1 {
            refqbf(xq, yq, iq, etax, etay, xix, xiy)
        } else {
            qbf(xq, yq, it, iq, &bump.node, &bump.xc, &bump.yc)
        };
        let ip_local = bump.node[it * NNODES + iq];

        for iuk in 0..2 {
            let iun = bump.indx[ip_local * 2 + iuk];
            if 0 < iun {
                un[iuk] += bb * bump.g[(iun - 1) as usize];
                unx[iuk] += bx * bump.g[(iun - 1) as usize];
                uny[iuk] += by * bump.g[(iun - 1) as usize];
            } else if iun == -1 {
                let ubc = ubdry(iuk + 1, bump.yc[ip_local]);
                un[iuk] += bb * ubc;
                unx[iuk] += bx * ubc;
                uny[iuk] += by * ubc;
            }
        }
    }
    (un, unx, uny)
}
//
//* --------------------------------------------------------------------
//*  UVAL - evaluate velocities at a given quadrature point
//* --------------------------------------------------------------------
pub fn uval(
    etax: f64,
    etay: f64,
    bump: &bump_struct::Bump,
    it: usize,
    xix: f64,
    xiy: f64,
    xq: f64,
    yq: f64,
) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let mut un = [0.0; 2];
    let mut unx = [0.0; 2];
    let mut uny = [0.0; 2];

    for iq in 0..NNODES {
        let (bb, bx, by) = if bump.isotri[it] == 1 {
            refqbf(xq, yq, iq, etax, etay, xix, xiy)
        } else {
            qbf(xq, yq, it, iq, &bump.node, &bump.xc, &bump.yc)
        };
        let ip_local = bump.node[it * NNODES + iq];

        for iuk in 0..2 {
            let iun = bump.indx[ip_local * 2 + iuk];
            if 0 < iun {
                un[iuk] += bb * bump.g[(iun - 1) as usize];
                unx[iuk] += bx * bump.g[(iun - 1) as usize];
                uny[iuk] += by * bump.g[(iun - 1) as usize];
            } else if iun == -1 {
                let ubc = ubdry(iuk + 1, bump.yc[ip_local]);
                un[iuk] += bb * ubc;
                unx[iuk] += bx * ubc;
                uny[iuk] += by * ubc;
            }
        }
    }
    (un, unx, uny)
}
//
//* --------------------------------------------------------------------
//*  UV_WRITE - write a velocity file
//* --------------------------------------------------------------------
pub fn uv_write(bump: &bump_struct::Bump, f: &[f64], file: &mut dyn std::io::Write) {
    for ip in 0..bump.n_points {
        let k_u = bump.indx[ip * 2];
        let u = if k_u < 0 {
            ubdry(1, bump.yc[ip])
        } else if k_u == 0 {
            0.0
        } else {
            f[(k_u - 1) as usize]
        };

        let k_v = bump.indx[ip * 2 + 1];
        let v = if k_v < 0 {
            ubdry(2, bump.yc[ip])
        } else if k_v == 0 {
            0.0
        } else {
            f[(k_v - 1) as usize]
        };

        writeln!(file, "  {:14e}  {:14e}", u, v).unwrap();
    }
}

//* --------------------------------------------------------------------
//*  XY_WRITE - write node coordinate data
//* --------------------------------------------------------------------
pub fn xy_write(bump: &bump_struct::Bump, file: &mut dyn std::io::Write) {
    for ip in 0..bump.n_points {
        writeln!(file, "  {:14e}  {:14e}", bump.xc[ip], bump.yc[ip]).unwrap();
    }
}
