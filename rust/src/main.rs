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
  let uv_dir: &str = "data/uv";
  let mut uv_file: String = format!("{}/uv_000.txt", uv_dir);
  let xy_dir: &str = "data/xy";
  let mut xy_file: String = format!("{}\\xy_000.txt", xy_dir);

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
  println!("  Rust version");
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

  let ibump: i32 = 2;
  funciones::setgrd(ibump, &mut bump);

  bump.f.resize(bump.neqn, 0.0);
  bump.g.resize(bump.neqn, 0.0);
  bump.res.resize(bump.neqn, 0.0);
  bump.sens.resize(bump.neqn, 0.0);

  bump.ypert = bump.aprof;
  funciones::setxy(&mut bump);
  funciones::setqud(&mut bump);
  funciones::setbas(&mut bump);
  funciones::setlin(&mut bump);
  funciones::setban(&mut bump);

  bump.a.resize(bump.nrow * bump.neqn, 0.0);

  funciones::nstoke(&mut bump);
  funciones::resid(&mut bump);

  bump.uprof = funciones::getg(&bump.g, &bump.iline, bump.my);

  if 1 <= bump.iwrite {
    println!();
    println!("Velocity profile:");
    println!();
    for i in 0..bump.my {
      print!("{}\t", bump.uprof[i]);
      if (i + 1) % 5 == 0 { println!(); }
    }
    println!();
  }

  funciones::gram(&mut bump);

  xy_file = funciones::file_name_inc(&xy_file);
  {
    let mut f_xy = std::fs::File::create(&xy_file).unwrap();
    funciones::xy_write(&bump, &mut f_xy);
  }

  uv_file = funciones::file_name_inc(&uv_file);
  {
    let mut f_uv = std::fs::File::create(&uv_file).unwrap();
    funciones::uv_write(&bump, &bump.f, &mut f_uv);
  }

  for i in 0..bump.neqn {
    bump.g[i] = 0.0;
  }

  let mut iter = 1;
  loop {
    if iter > bump.maxsec { break; }
    println!();
    println!("Secant iteration {}", iter);

    bump.numsec += 1;

    bump.ypert = bump.anew;
    funciones::setxy(&mut bump);
    funciones::setqud(&mut bump);
    funciones::setbas(&mut bump);
    funciones::nstoke(&mut bump);

    bump.uprof = funciones::getg(&bump.g, &bump.iline, bump.my);

    if 1 <= bump.iwrite {
      println!();
      println!("Velocity profile:");
      println!();
      for i in 0..bump.my {
        print!("{}\t", bump.uprof[i]);
        if (i + 1) % 5 == 0 { println!(); }
      }
      println!();
    }

    funciones::linsys(&mut bump, -2);

    bump.dcda = funciones::getg(&bump.sens, &bump.iline, bump.my);

    if 2 <= bump.iwrite {
      println!();
      println!("Sensitivities:");
      println!();
      for i in 0..bump.my {
        print!("{}\t", bump.dcda[i]);
        if (i + 1) % 5 == 0 { println!(); }
      }
      println!();
    }

    bump.rjpnew = 0.0;
    for i in 0..bump.my {
      let mut temp = -bump.r[i];
      for j in 0..bump.my {
        temp += bump.gr[i * bump.my + j] * bump.uprof[j];
      }
      bump.rjpnew += 2.0 * bump.dcda[i] * temp;
    }

    xy_file = funciones::file_name_inc(&xy_file);
    {
      let mut f_xy = std::fs::File::create(&xy_file).unwrap();
      funciones::xy_write(&bump, &mut f_xy);
    }

    uv_file = funciones::file_name_inc(&uv_file);
    {
      let mut f_uv = std::fs::File::create(&uv_file).unwrap();
      funciones::uv_write(&bump, &bump.f, &mut f_uv);
    }

    println!();
    println!("  Parameter = {}, J prime = {}", bump.anew, bump.rjpnew);

    if 1 < iter {
      let denom = bump.rjpnew - bump.rjpold;
      bump.anext = if denom.abs() > 1e-30 {
        bump.aold - bump.rjpold * (bump.anew - bump.aold) / denom
      } else {
        bump.anew
      };
    }

    bump.aold = bump.anew;
    bump.anew = bump.anext;
    bump.rjpold = bump.rjpnew;

    let test = if bump.anew != 0.0 { (bump.anew - bump.aold).abs() / bump.anew } else { 0.0 };

    println!("  New value of parameter = {}", bump.anew);
    println!("  Convergence test = {}", test);

    if (bump.anew - bump.aold).abs() <= bump.anew.abs() * bump.tolsec && 1 < iter {
      println!("Secant iteration converged.");
      break;
    }
    iter += 1;
  }
  if bump.maxsec < iter {
    println!("  Secant iteration failed to converge.");
  }

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
      .open("times.txt").unwrap();
    write!(file, "{}\n", duracion.as_secs_f64()).unwrap();
    file.flush().unwrap();
  }
}