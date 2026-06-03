# Traducción de un código F90 legacy a Rust, C++ y Python

## Descripción de las traducciones

### Código original

```
fortran/bump.f90
```

Compilado con (en el directorio `fortran`):

```bash
flang fortran/bump.f90 -O3 -o ./fortran/build/bump.exe
```

### Código Python

```
python/bump.py
```

Para ejecutar:

```bash
python bump.py
```

### Código C++

```
c++/src/main.cpp
```

Compilado con (en el directorio `c++`):

```bash
clang++ src/main.cpp -std=c++23 -O3 -o build/main.exe
```

### Código Rust

```
rust/src/main.rs
```

Para ejecutar:

```bash
cargo run --release
```

## Toma de tiempos

A excepción de F90, Rust y C++ guardarán el tiempo de ejecución en un archivo aparte predefinido (`tiempos.txt`) en cada uno de sus directorios. Un script de Python se encargará de ejecutar el programa *n* veces y promediará el tiempo que se guarda en el archivo ya mencionado.

### Caso Fortran 90

En el caso de Fortran 90, Python mide el tiempo de ejecución de esta forma:

```python
start: float = time.perf_counter()
subprocess.run(
    comando,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL
)
end: float = time.perf_counter()
tiempos.append(end - start)
```

y se promedia el tiempo medido por Python.

> **Nota:** la diferencia con el tiempo dado por el propio programa es despreciable (1 % – 2 % de diferencia o incluso menos).

## Resultados

```
Repeticiones: 10

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
|  Lenguaje  |   Media  |
+============+==========+
|    f90     | 0.381152 |
+------------+----------+
|    c++     | 0.79376  |
+------------+----------+
|    rust    | 0.690695 |
+------------+----------+
```


## Requerimentos necesarios:

### Rust

````bash
chrono
````

### Python
````
pip install tabulate
````

#### Para su version de bump
````
pip install numpy
````