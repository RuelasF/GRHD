# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:48:31 2026

@author: francisco ruelas
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Intenta importar gwpy para descargar la curva observacional real
try:
    from gwpy.timeseries import TimeSeries
    has_gwpy = True
except ImportError:
    has_gwpy = False
    print("ADVERTENCIA: gwpy no está instalado. Ejecuta 'pip install gwpy'")

# ====================================================================
# 1. PARÁMETROS FÍSICOS Y ESCALADO (Agujero Negro 10 M_sol a 10 kpc)
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
# GW_signal.dat ya contiene la proyección face-on para un observador sobre +z.
factor_geometrico = 1.0

# ====================================================================
# 2. CARGA DE DATOS Y CONVERSIÓN DE TU SIMULACIÓN
# ====================================================================
datos = np.loadtxt('FishboneMoncrief_Equatorial_test/weno5_hlle/GW_signal.dat')

tiempo_num = datos[:, 0]
h_plus_num = datos[:, 1]
h_cross_num = datos[:, 2]

# Aislar régimen no lineal (T > 1000)
mascara = tiempo_num >= 1000.0
tiempo_sec = tiempo_num[mascara] * T_to_sec
h_plus_fisico = h_plus_num[mascara] * factor_distancia * factor_geometrico
h_cross_fisico = h_cross_num[mascara] * factor_distancia * factor_geometrico

# La extracción se realiza cada cierto número de pasos adaptativos. Welch
# requiere muestreo uniforme, por lo que interpolamos antes del análisis.
if len(tiempo_sec) < 2 or np.any(np.diff(tiempo_sec) <= 0.0):
    raise ValueError('La serie GW necesita al menos dos tiempos estrictamente crecientes.')
tiempo_uniforme = np.linspace(tiempo_sec[0], tiempo_sec[-1], len(tiempo_sec))
h_plus_fisico = np.interp(tiempo_uniforme, tiempo_sec, h_plus_fisico)
h_cross_fisico = np.interp(tiempo_uniforme, tiempo_sec, h_cross_fisico)
tiempo_sec = tiempo_uniforme
dt_sec = tiempo_sec[1] - tiempo_sec[0]
fs_hz = 1.0 / dt_sec

# ====================================================================
# 3. PROCESAMIENTO: DETREND, WINDOWING Y WELCH
# ====================================================================
h_plus_detrend = signal.detrend(h_plus_fisico, type='linear')
h_cross_detrend = signal.detrend(h_cross_fisico, type='linear')

ventana = signal.windows.tukey(len(tiempo_sec), alpha=0.1)
h_plus_limpio = h_plus_detrend * ventana
h_cross_limpio = h_cross_detrend * ventana

# Resolución optimizada para bajas frecuencias
nperseg_puntos = min(int(0.05 / dt_sec), len(tiempo_sec))

frecuencias_hz, psd_plus = signal.welch(
    h_plus_limpio, fs_hz, window='hann', nperseg=nperseg_puntos)
frecuencias_hz, psd_cross = signal.welch(
    h_cross_limpio, fs_hz, window='hann', nperseg=nperseg_puntos)

# CONVERSIÓN A ASD (Amplitude Spectral Density)
asd_plus = np.sqrt(psd_plus)
asd_cross = np.sqrt(psd_cross)

idx_max = np.argmax(asd_plus[1:]) + 1
f_dominante_hz = frecuencias_hz[idx_max]

# ====================================================================
# 4. CURVA DE SENSIBILIDAD EMPÍRICA DE LIGO O4 (GWPY)
# ====================================================================
# ====================================================================
# 4. CURVA DE SENSIBILIDAD EMPÍRICA DE LIGO (GWPY - O3)
# ====================================================================
if has_gwpy:
    print("Descargando 1024 segundos de datos observacionales estables de LIGO...")
    # Tiempo GPS garantizado de la campaña O3b (Detector estable y "locked")
    gps_start = 1251332000
    gps_end = gps_start + 1024

    data_empirica = TimeSeries.fetch_open_data('H1', gps_start, gps_end, cache=True)

    # ASD con resolución de 0.25 Hz para que se vean bien las líneas espectrales
    asd_ligo_emp = data_empirica.asd(fftlength=4, overlap=2, method="welch")

    # Extraer valores y limpiar posibles NaNs/Ceros del instrumento
    freqs_emp = asd_ligo_emp.frequencies.value
    vals_emp = asd_ligo_emp.value
    valido = (freqs_emp > 0) & (vals_emp > 0)

# ====================================================================
# 5. VISUALIZACIÓN OBSERVACIONAL FINAL
# ====================================================================
plt.figure(figsize=(11, 7))

plt.loglog(frecuencias_hz, asd_plus, label=r'Simulación GRHD ($h_+$)',
           color='darkblue', linewidth=1.5)
plt.loglog(frecuencias_hz, asd_cross, label=r'Simulación GRHD ($h_\times$)',
           color='darkred', linewidth=1.5, alpha=0.8)

if has_gwpy:
    plt.loglog(freqs_emp[valido], vals_emp[valido],
               label='Ruido Empírico LIGO (Hanford)', color='black', linewidth=1.0, alpha=0.7)

plt.axvline(f_dominante_hz, color='gray', linestyle=':', linewidth=2,
            label=f'Frecuencia Principal: {f_dominante_hz:.1f} Hz')

plt.title(rf'Strain de Onda Gravitacional vs Ruido Observacional LIGO' +
          '\n'+rf'($D = 10$ kpc, $M_{{BH}} = {M_BH} M_\odot$)', fontsize=15)
plt.xlabel('Frecuencia [Hz]', fontsize=13)
plt.ylabel(r'Densidad Espectral de Amplitud (ASD) [$1/\sqrt{\mathrm{Hz}}$]', fontsize=13)

# Ajuste maestro para la ventana auditiva y amplitudes de LIGO
plt.xlim(10, 2000)
plt.ylim(1e-25, 1e-14)  # Abarca desde el ruido cuántico hasta tu inestabilidad masiva

plt.grid(True, which="major", ls="-", alpha=0.5)
plt.grid(True, which="minor", ls=":", alpha=0.3)
plt.legend(fontsize=11, loc='upper right')
plt.tight_layout()
plt.show()
