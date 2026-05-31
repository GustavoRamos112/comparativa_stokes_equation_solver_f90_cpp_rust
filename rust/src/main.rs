//* --------------------------------------------------------------------
//*  Main program
//* --------------------------------------------------------------------

use std::path::Path;
use std::time::Instant;
use std::fs::OpenOptions;
use std::io::Write;

mod program;
use program::*;


fn main() {
  let mut itype: i32;
  let mut temp: f64;
  let mut denom: f64;
  let mut test: f64;
  
  let uv_dir: &str = "data/uv";
  let uv_file: String = format!("{}/uv_000.txt", uv_dir);
  let xy_dir: &str = "data/xy";
  let xy_file: String = format!("{}\\xy_000.txt", xy_dir);
  
  let uv_path: &Path = Path::new(&uv_dir);
  let xy_path: &Path = Path::new(&xy_dir);
  
  let nx: usize = 21;
  let ny: usize = 7;
  let mx: usize = 2 * nx - 1;
  let my: usize = 2 * ny - 1;
  let nnodes: usize = 6;
  let nquad: usize = 3;

  let mut bump: bump_struct::Bump = bump_struct::Bump::new(nx, ny, mx, my, nnodes, nquad);

  if !uv_path.exists() {
    std::fs::create_dir_all(uv_path).unwrap();
  }

  if !xy_path.exists() {
    std::fs::create_dir_all(xy_path).unwrap();
  }

  let inicio: Instant = Instant::now();
  funciones::timestamp();
  println!();
  println!("BUMP");
  println!("  C++ version");
  println!("  Control problem for channel flow over a bump.");
  println!();
  println!("  The bump will be generated with a height of {}", bump.aprof);
  println!();
  println!("  NX = {}", nx);
  println!("  NY = {}", ny);
  println!("  Number of elements = {}", bump.nelemn);
  println!("  Reynolds number =  {}", bump.reynld);
  println!("  Secant tolerance = {}", bump.tolsec);
  println!("  Newton tolerance = {}", bump.tolnew);

  let mut ibump: i32 = 2;
  /* (bump._long, bump.neqn) = */ funciones::setgrd(
    ibump, &mut bump
  );

  bump.f.resize(bump.neqn, 0.0);
  bump.g.resize(bump.neqn, 0.0);
  bump.res.resize(bump.neqn, 0.0);
  bump.sens.resize(bump.neqn, 0.0);


  bump.ypert = bump.aprof;
  funciones::setxy(
    &mut bump
  );

  funciones::setqud(
    &mut bump
  );

  funciones::setbas(&mut bump);

  //setlin(iline, indx, iwrite, _long, xlngth, xprof);

//  tie(nband, nlband, nrow) = setban(
//    indx, insc, maxrow, nband, nlband, node, nrow
//  );
//
//  a.resize(nrow * neqn, 0.0);
//
//  numnew = nstoke(
//    a, area, f, g, indx, insc, isotri,
//    maxrow, nband,
//    neqn, nlband,
//    node, nrow, numnew,
//    phi, psi, reynld, tolnew,
//    xc, xm, yc, ym
//  );
//  
//  resid(
//    area, g, indx, insc, isotri,
//    iwrite, neqn, node,
//    phi, psi, res, reynld,
//    xc, xm, yc, ym
//  );
//
//  uprof = getg(g, iline, neqn);
//
//  if (1 <= iwrite) {
//    println!();
//    println!("Velocity profile:");
//    println!();
//    for (int i = 0; i < my; i++) {
//      print("{}\t", uprof[i]);
//      if ((i + 1) % 5 == 0)
//        println!();
//    }
//    println!();
//  }
//
//  gram(
//    gr, iline, indx, iwrite, node, r, uprof, xc, xprof, yc
//  );
//
//  xy_file = file_name_inc(xy_file);
//
//  ofstream f_xy(xy_file);
//  if (f_xy.is_open()) {
//    xy_write(f_xy, xc, yc);
//    f_xy.close();
//  }
//
//  uv_file = file_name_inc(uv_file);
//  ofstream f_uv(uv_file);
//  if (f_uv.is_open()) {
//    uv_write(f, indx, f_uv, neqn, yc);
//    f_uv.close();
//  }
//
//  for(int i = 0; i < neqn; i++) {
//    g[i] = 0.0;
//  }
//  
//  int iter;
//  for (iter = 1; iter <= maxsec; iter++) {
//    println!();
//    println!("Secant iteration {}", iter);
//
//    numsec += 1;
//
//    ypert = anew;
//    setxy(iwrite, _long, xc, xlngth, yc, ylngth, ypert);
//
//    setqud(area, isotri, iwrite, node, xc, xm, yc, ym);
//
//    setbas(isotri, node, phi, psi, xc, xm, yc, ym);
//
//    numnew = nstoke(
//      a, area, f, g, indx, insc,
//      isotri, maxrow, nband,
//      neqn, nlband, node,
//      nrow, numnew, phi,
//      psi, reynld, tolnew, xc,
//      xm, yc, ym
//    );
//
//    uprof = getg(g, iline, neqn);
//
//    if (1 <= iwrite) {
//      println!();
//      println!("Velocity profile:");
//      println!();
//      for (int i = 0; i < my; i++) {
//        print("{}\t", uprof[i]);
//        if ((i + 1) % 5 == 0) println!();
//      }
//      println!();
//    }
//
//    itype = -2;
//
//    linsys(
//      a, area, sens, g, indx, insc, isotri, itype, maxrow,
//      nband, neqn, nlband, node,
//      nrow, phi, psi,
//      reynld, xc, xm, yc, ym
//    );
//
//    dcda = getg(sens, iline, neqn);
//
//    if (2 <= iwrite) {
//      println!();
//      println!("Sensitivities:");
//      println!();
//      for (int i = 0; i < my; i++) {
//        print("{}\t", dcda[i]);
//        if ((i + 1) % 5 == 0) {
//          println!();
//        }
//      }
//      println!();
//    }
//
//    rjpnew = 0.0;
//    for (int i = 0; i < my; i++) {
//      temp = -r[i];
//      for (int j = 0; j < my; j++) {
//        temp += gr[i * my + j] * uprof[j];
//      }
//      rjpnew += 2.0 * dcda[i] * temp;
//    }
//
//    xy_file = file_name_inc(xy_file);
//
//    ofstream f_xy(xy_file);
//    if (f_xy.is_open()) {
//      xy_write(f_xy, xc, yc);
//      f_xy.close();
//    }
//
//    uv_file = file_name_inc(uv_file);
//    ofstream f_uv(uv_file);
//    if (f_uv.is_open()) {
//      uv_write(f, indx, f_uv, neqn, yc);
//      f_uv.close();
//    }
//
//    println!();
//    println!("  Parameter = {}, J prime = {}", anew, rjpnew);
//
//    if (1 < iter) {
//      denom = rjpnew - rjpold;
//      if (abs(denom) > 1e-30) 
//        anext = aold - rjpold * (anew - aold) / denom;
//      else anext = anew;
//    }
//
//    aold = anew;
//    anew = anext;
//    rjpold = rjpnew;
//
//    if (anew != 0.0) test = abs(anew - aold) / anew;
//    else test = 0.0;
//
//    println!("  New value of parameter = {}", anew);
//    println!("  Convergence test = {}", test);
//
//    if ((abs(anew - aold) <= abs(anew) * tolsec) and (1 < iter)) {
//      println!("Secant iteration converged.");
//      break;
//    }
//  }
//  if (maxsec < iter)
//    println!("  Secant iteration failed to converge.");
  
  let duracion: std::time::Duration = inicio.elapsed();
  
  println!();
  println!("  Total execution time = {:?}", duracion);
  println!("  Number of secant steps = {}", bump.numsec);
  println!("  Number of Newton steps = {}", bump.numnew);
  println!();
  println!("BUMP:");
  println!("  Normal end of execution.");
  println!();
  funciones::timestamp();
  
  if bump.save_times {
    let mut file = OpenOptions::new()
      .append(true)
      .create(true) // Creates the file if it's missing
      .open("archive.txt").unwrap();
    write!(file, "{}", duracion.as_secs()).unwrap();
  }
}