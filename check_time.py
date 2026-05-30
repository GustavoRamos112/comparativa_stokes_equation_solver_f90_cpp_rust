
import time
import subprocess
import os
# pyrefly: ignore [untyped-import]
from tabulate import tabulate

def check_times_file (dir_times: str) -> None:
  if os.path.exists(dir_times):
    os.remove(dir_times)

  if not os.path.exists(dir_times):
    with open(dir_times, "w") as archivo:
      archivo.write("")

def medir_tiempo(directorio: str, comando: list, repeticiones: int) -> tuple[list[float], float]:
    # Cambiar de directorio de forma segura
    original_dir: str = os.getcwd()
    os.chdir(directorio)
    
    # 1. Warm-up (Lanzar una vez para cargar en caché)
    subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    tiempos: list[float] = []
    for _ in range(repeticiones):
        # 2. Usar perf_counter para alta precisión
        start: float = time.perf_counter()
        subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        end: float = time.perf_counter()
        tiempos.append(end - start)
    
    os.chdir(original_dir)
    
    # Retornar el promedio
    return (tiempos, sum(tiempos) / len(tiempos))

if __name__ == "__main__":
  _dir = "c++"
  _command: list[str] = ["build/main.exe"]
  _repetitions = 20
  
  archivo_times: str = "times.txt"
  dir_times: str = f"{_dir}/{archivo_times}"

  check_times_file(dir_times)  

  tiempos_python: list[float] = []
  media_python: float = 0.0

  tiempos_python, media_python = medir_tiempo(
    directorio=_dir, comando=_command, 
    repeticiones=_repetitions
  )

  with open(dir_times, "r") as archivo:
    tiempos_programa: list[float] = [float(l)/1000 for l in archivo.readlines()]

  tiempos_programa.pop(0)

  
  tiempos_generales: list[list[float]] = [
    [tiempos_python[i], tiempos_programa[i]] 
    for i in range(len(tiempos_python))
  ]

  
  media_programa: float = sum(tiempos_programa) / len(tiempos_programa)

  print("Repeticiones: ", _repetitions)

  print(tabulate(
    tiempos_generales, 
    headers=["Tiempo Python", "Tiempo comando"], 
    tablefmt="grid")
  )
  print(tabulate(
    [
      ["Python", media_python],
      ["Programa", media_programa]
    ], 
    headers=["Medida", "Tiempo"], 
    tablefmt="grid")
  )
  
  

