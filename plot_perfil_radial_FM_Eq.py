# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:14:06 2026

@author: francisco ruelas
"""
import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
import glob

# 1. Configuración
nx = 400
variable = "Density"  # Cambia a "Velocity_X" o "rho"

# 2. Busca y ordena todos los archivos .vtk automáticamente
archivos_vtk = sorted(glob.glob("*.vtk"))

if not archivos_vtk:
    print("No se encontraron archivos. Verifica la ruta o el prefijo.")
    exit()

print(f"Graficando {len(archivos_vtk)} pasos temporales...")

# 3. Prepara la gráfica
plt.figure(figsize=(12, 7))

# 4. Ciclo maestro para graficar todos los tiempos superpuestos
for i, archivo in enumerate(archivos_vtk):
    if i < 120:
        malla = pv.read(archivo)

        # Extraer el eje X (radio físico) y los datos
        coords = malla.points[:nx]
        r_phys = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2 + coords[:, 2]**2)
        datos = malla.point_data[variable][:nx]

        # El secreto visual: alpha=0.1 hace que cada línea sea casi transparente.
        # Al sumarse 100 líneas, resaltan los modos de oscilación.
        # Solo le ponemos 'label' a la primera línea para que la leyenda no se llene.
        etiqueta = "Evolución temporal" if i == 0 else None
        plt.plot(r_phys, datos, color='teal', alpha=0.5, label=etiqueta)
    else:
        pass

# 5. Formato limpio y científico
plt.xlabel("Radio (r)", fontsize=14)
plt.ylabel(variable, fontsize=14)
# plt.yscale('log')
plt.title("Envolvente de Oscilación del Toro Ecuatorial", fontsize=16)

# Opcional: Zoom al borde interno del disco donde ocurre la respiración
plt.ylim(1e-11, 2.0)
plt.xlim(1, 40)

plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
