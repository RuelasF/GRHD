# -*- coding: utf-8 -*-
"""
Created on Mon May 11 16:03:12 2026

@author: francisco ruelas
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import matplotlib.colors as colors

# ====================================================================
# 1. PARÁMETROS FÍSICOS Y ESCALADO (10 M_sol a 10 kpc)
# ====================================================================
M_BH = 10.0
M_sol_g = 1.98892e33
G = 6.6743e-8
c = 2.9979e10
kpc_cm = 3.085677581e21

T_to_sec = (G * (M_BH * M_sol_g)) / (c**3)
L_to_cm = (G * (M_BH * M_sol_g)) / (c**2)

D_extraccion_cm = 1000.0 * L_to_cm
D_observador_cm = 10000.0 * kpc_cm

factor_distancia = D_extraccion_cm / D_observador_cm
factor_geometrico = 2.0

# ====================================================================
# 2. CARGA DE DATOS
# ====================================================================
datos = np.loadtxt('FishboneMoncrief_Equatorial_test/weno5_hlle/GW_signal.dat')

tiempo_num = datos[:, 0]
h_plus_num = datos[:, 1]

# Aislar régimen no lineal (T > 1000)
mascara = tiempo_num >= 100.0
tiempo_sec = tiempo_num[mascara] * T_to_sec

# Normalizamos para que el tiempo inicie en 0 en la gráfica
tiempo_sec = tiempo_sec - tiempo_sec[0]

h_plus_fisico = h_plus_num[mascara] * factor_distancia * factor_geometrico

dt_sec = np.mean(np.diff(tiempo_sec))
fs_hz = 1.0 / dt_sec

# Limpieza básica
h_plus_limpio = signal.detrend(h_plus_fisico, type='linear')
ventana = signal.windows.tukey(len(tiempo_sec), alpha=0.1)
h_plus_limpio = h_plus_limpio * ventana

# ====================================================================
# 3. ESPECTROGRAMA (Mapa de Tiempo-Frecuencia)
# ====================================================================
# Ajustar la resolución: Ventanas de ~0.04 segundos
nperseg_puntos = int(0.04 * fs_hz)
# Un 'overlap' del 95% da una imagen extremadamente fluida y continua
noverlap_puntos = int(nperseg_puntos * 0.95)

f, t_spec, Sxx = signal.spectrogram(h_plus_limpio, fs=fs_hz, window='hann',
                                    nperseg=nperseg_puntos, noverlap=noverlap_puntos,
                                    scaling='density')

# ====================================================================
# 4. VISUALIZACIÓN CIENTÍFICA (Estilo LIGO)
# ====================================================================
plt.figure(figsize=(12, 6))

# Usamos LogNorm para resaltar el pico principal sobre el ruido numérico de fondo
# El colormap 'viridis' o 'magma' son los estándares actuales en astrofísica
mesh = plt.pcolormesh(t_spec, f, Sxx,
                      norm=colors.LogNorm(vmin=np.max(Sxx)*1e-4, vmax=np.max(Sxx)),
                      cmap='magma', shading='gouraud')

# Barra de energía
cbar = plt.colorbar(mesh, pad=0.02)
cbar.set_label(r'Densidad de Potencia (Strain$^2$/Hz)', fontsize=12)

# Frecuencia principal teórica para referencia visual
f_dominante = 180.0
plt.axhline(f_dominante, color='cyan', linestyle='--', linewidth=1.5, alpha=0.7,
            label=f'Frecuencia Principal ({f_dominante} Hz)')

plt.title(rf'Espectrograma de Radiación Gravitacional (Inestabilidad Hidrodinámica)' +
          '\n'+rf'$M_{{BH}} = {M_BH} M_\odot$, D = 10 kpc', fontsize=15)
plt.ylabel('Frecuencia [Hz]', fontsize=13)
plt.xlabel('Tiempo [s]', fontsize=13)

# Enfocamos en la banda auditiva (LIGO) y en escala logarítmica
plt.ylim(1, 20000)
plt.yscale('log')

plt.legend(loc='upper right', fontsize=11)
plt.tight_layout()
plt.show()
