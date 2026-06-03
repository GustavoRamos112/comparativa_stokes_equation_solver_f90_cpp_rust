# Traduccion de un codigo f90 legacy a rust, c++ y python


## Descripcion de las traduciones
### Codigo original: 

````
fortran\bump.f90
````

Compilado con (en el directorio fortran):
````
flang fortran\bump.f90 -O3 -o .\fortran\build\bump.exe
````

### Codigo python

````
python\bump.py
````

Para ejecutar:
````
python bump.py
````


### Codigo C++

````
c++\src\main.cpp
````

Compilado con (en el directorio c++):
````
clang++ src\main.cpp -std=c++23 -o3 -o build\main.exe
````

### Codigo Rust

````
rust\src\main.rs
````

Para ejecutar:
````
cargo run --release
````

# Traduccion de un codigo f90 legacy a rust, c++ y python


## Descripcion de las traduciones
### Codigo original: 

````
fortran\bump.f90
````

Compilado con (en el directorio fortran):
````
flang fortran\bump.f90 -O3 -o .\fortran\build\bump.exe
````

### Codigo python

````
python\bump.py
````

Para ejecutar:
````
python bump.py
````


### Codigo C++

````
c++\src\main.cpp
````

Compilado con (en el directorio c++):
````
clang++ src\main.cpp -std=c++23 -o3 -o build\main.exe
````

### Codigo Rust

````
rust\src\main.rs
````

Para ejecutar:
````
cargo run --release
````

## Toma del tiempo

A excepcion de f90, rust y c++ guardaran el tiempo de ejecucion en un archivo aparte predefinido como tiempos.txt en cada uno, un script de python se encargara de correr el programa n veces y promediara el tiempo que se guarda en el archivo ya mencionado.

### Caso Fortran 90
En el caso de fortran 90 Python mide el timepo de ejcucion de esta forma:
````python
start: float = time.perf_counter()
subprocess.run(
  comando, 
  stdout=subprocess.DEVNULL, 
  stderr=subprocess.DEVNULL
)
end: float = time.perf_counter()
tiempos.append(end - start)
````

y se promedia el tiempo medido por python.

> Nota: la diferencia con el tiempo dado por el propio programa es despreciable (1%~2% de diferencia o incluso menos).

## Resultados:

````
Repeticiones:  10
+--------------+--------------+---------------+
|   Tiempo F90 |   Tiempo Cpp |   Tiempo Rust |
+==============+==============+===============+
|     0.710615 |     0.807117 |      0.606199 |
+--------------+--------------+---------------+
|     0.354249 |     1.20506  |      0.592848 |
+--------------+--------------+---------------+
|     0.373087 |     1.06755  |      0.583329 |
+--------------+--------------+---------------+
|     0.312353 |     0.715874 |      1.06763  |
+--------------+--------------+---------------+
|     0.331267 |     0.713371 |      1.12133  |
+--------------+--------------+---------------+
|     0.334456 |     0.682206 |      0.716348 |
+--------------+--------------+---------------+
|     0.341999 |     0.696173 |      0.579489 |
+--------------+--------------+---------------+
|     0.323372 |     0.683231 |      0.556674 |
+--------------+--------------+---------------+
|     0.322331 |     0.677896 |      0.545364 |
+--------------+--------------+---------------+
|     0.407791 |     0.689118 |      0.537743 |
+--------------+--------------+---------------+
+------------+----------+
|  Lenguaje  |    Media |
+============+==========+
|    f90     | 0.381152 |
+------------+----------+
|    c++     | 0.79376  |
+------------+----------+
|    rust    | 0.690695 |
+------------+----------+
````