# Manual de usuario y desarrollo de GRHD2

Este manual describe cómo compilar, configurar, ejecutar, verificar y extender
el código contenido en este repositorio. Su fuente de verdad es el código
Fortran actual y los archivos versionados en `par/` y `tests/`.

> **Estado físico:** GRHD para un fluido perfecto sobre una métrica analítica
> fija. No es GRMHD, no evoluciona el espacio-tiempo y no incluye autogravedad
> del fluido. La salida `GW_signal.dat` es un proxy cuadrupolar Finn--Evans de
> campo débil, no un *strain* observable ni una extracción gauge-invariant.

## 1. Requisitos

El flujo probado utiliza:

- GNU Fortran (`gfortran`);
- GNU Make;
- OpenMP, suministrado por `gfortran`;
- Python 3 para los scripts de análisis;
- NumPy y Matplotlib para las gráficas;
- VisIt o ParaView, opcionales, para inspeccionar los VTK.

En Debian/Ubuntu puede instalarse la base con:

```bash
sudo apt update
sudo apt install gfortran make python3 python3-numpy python3-matplotlib
```

La instalación de paquetes es una operación del sistema y no forma parte del
`Makefile`.

## 2. Compilación rápida

Desde la raíz del repositorio:

```bash
make -j
```

El ejecutable resultante es `grhd2`. Las banderas de producción actuales son
`-O3`, OpenMP, optimización específica de la máquina y LTO. Para una auditoría
de memoria o punto flotante conviene reemplazar temporalmente `FFLAGS` por las
banderas de depuración comentadas en el `Makefile`.

Para eliminar objetos, módulos, ejecutable y binarios de prueba:

```bash
make clean
```

## 3. Primer cálculo

Compruebe primero el archivo de parámetros sin iniciar la simulación:

```bash
./grhd2 --check par/fishbone_equatorial.par
```

Si la validación termina correctamente, ejecute:

```bash
export OMP_NUM_THREADS=8
./grhd2 par/fishbone_equatorial.par
```

`OMP_NUM_THREADS` debe ajustarse al equipo. El programa exige exactamente un
archivo `.par`; ya no existe una ejecución válida sin argumentos.

El ejemplo genera un toro de Fishbone--Moncrief ecuatorial en Kerr--Schild con
`a=0.9`, malla `400 x 1 x 200`, WENO5--HLLE y radio logarítmico. La salida se
guarda, por defecto, en:

```text
FM_KS_a0.9_data/weno5_hlle/
```

Aunque se evolucionan las direcciones radial y azimutal, una sola celda polar
no resuelve la estructura vertical. Esta configuración debe describirse como
**modelo ecuatorial 2.5D**, no como simulación 3D.

## 4. Archivos de parámetros

### 4.1 Sintaxis

Cada línea activa tiene la forma:

```text
nombre = valor       # comentario opcional
```

Las claves y valores simbólicos no distinguen mayúsculas de minúsculas. Las
rutas conservan su escritura. Una clave desconocida, duplicada, incompatible
con el problema o con valor no admisible detiene el programa con un mensaje que
indica archivo y línea.

Para dominios esféricos, `theta_min`, `theta_max`, `phi_min` y `phi_max` se
escriben como múltiplos de pi:

```text
theta_min = 0
theta_max = 1       # pi
phi_min   = 0
phi_max   = 2       # 2*pi
```

Se admiten cocientes, por ejemplo `theta_min = 1/2`.

### 4.2 Problemas disponibles

| `problem` | Descripción | Geometría usual | Estado recomendado de uso |
|---|---|---|---|
| `sod` | Tubo de choque relativista | Minkowski cartesiana 1D | Prueba física pendiente de reporte actualizado |
| `strong_shock` | Choque fuerte | Minkowski cartesiana 1D | Prueba física pendiente de reporte actualizado |
| `shu_osher` | Interacción choque--onda | Minkowski cartesiana 1D | Validación cuantitativa pendiente |
| `kelvin_helmholtz` | Inestabilidad de cizalla | Minkowski cartesiana 2D | Estudio de convergencia pendiente |
| `axisymmetric_jet` | Jet relativista axisimétrico | Minkowski cilíndrica | Validación sistemática pendiente |
| `convergence` | Onda sinusoidal advectada | Minkowski cartesiana 1D | Campaña de convergencia disponible |
| `michel` | Acreción esférica estacionaria | EF esférica | Solución analítica implementada; convergencia final pendiente |
| `dust` | Acreción de polvo | EF esférica | Exploratorio |
| `offaxis_blast` | Explosión fuera del eje | EF esférica | Exploratorio; proxy GW |
| `fishbone_equatorial` | Toro FM en `r-phi` | EF o KS esférica | Caso científico principal 2.5D |
| `fishbone_sagittal` | Corte FM en `r-theta` | EF esférica | Caso auxiliar |
| `fishbone_3d` | Toro FM con tres direcciones activas | EF o KS esférica | Implementado; requiere validar resolución |

Cada problema carga un *preset* desde `initialization.f90`; el archivo `.par`
sobrescribe después los valores especificados. Es aconsejable declarar de
forma explícita todos los parámetros relevantes para una corrida científica.

### 4.3 Métricas y geometrías

Valores admitidos:

```text
metric   = minkowski | eddington_finkelstein | kerr_schild
geometry = cartesian | cylindrical | spherical
```

Compatibilidad actual:

- Minkowski: cartesiana o cilíndrica;
- Eddington--Finkelstein: esférica, Schwarzschild y `spin = 0`;
- Kerr--Schild: esférica y `|spin| <= bh_mass`.

En coordenada radial logarítmica se usa `x = ln(r)` y se requiere `r_min > 0`:

```text
logarithmic_r = true
```

### 4.4 Métodos numéricos

```text
reconstruction = godunov | tvd | weno3 | mp5 | weno5
tvd_limiter    = minmod | superbee | mc
riemann_solver = hlle | hllc
shock_sensor   = true | false
```

La profundidad de celdas fantasma se selecciona automáticamente: una para
Godunov, dos para TVD/WENO3 y tres para MP5/WENO5. `HLLD` no es una opción
válida: sólo existe un esqueleto que aborta si se invoca internamente.

HLLE es el *baseline* robusto. HLLC está implementado para GRHD mediante un
marco ortonormal local de cada cara y vuelve a HLLE ante estados degenerados o
no admisibles.

El sensor de choques puede degradar localmente la reconstrucción hacia métodos
más robustos. Desactivarlo permite comparar reconstructores puros, pero no es
una decisión universalmente más estable.

### 4.5 Malla, tiempo y salida

Parámetros generales:

```text
nx = 400
ny = 1
nz = 200
r_min = 1.2
r_max = 40.0
theta_min = 0
theta_max = 1
phi_min = 0
phi_max = 2
logarithmic_r = true

final_time = 1000.0
cfl = 0.4
save_interval = 10.0
checkpoint_interval = 100.0

output_prefix = FM_KS_a0.9
output_folder = FM_KS_a0.9_data
vtk_mapping = physical
```

`vtk_mapping` acepta:

- `physical`: mapeo cartesiano físico que conserva la torsión Kerr--Schild;
- `untwisted`: representación oblata sin la rotación radial, sólo para
  visualización.

Esta opción no cambia el estado ni las coordenadas de evolución.

### 4.6 Diagnósticos y perturbaciones

```text
extract_gw = false
extract_mdot = false
ppi_diagnostics = true
diagnostic_stride = 100

perturbation = none | pressure_noise | density_noise | density_mode
perturbation_seed = 3435
perturbation_time = 1000.0
perturbation_amplitude = 0.01
perturbation_mode = 1.0
```

`pressure_noise` y `density_noise` usan una semilla reproducible. `density_mode`
inyecta un modo azimutal coherente cuyo número es `perturbation_mode`. La
amplitud debe satisfacer `0 <= perturbation_amplitude < 1`.

Una perturbación solicitada se aplica al comienzo del primer paso para el cual
el tiempo acumulado ya alcanzó `perturbation_time`; el instante efectivo se
registra en `perturbation_events.dat`.

### 4.7 Parámetros específicos por problema

El listado completo se mantiene en `par/README.txt`. Los prefijos son:

| Problema | Claves específicas |
|---|---|
| Sod | `sod_interface`, densidades y presiones izquierda/derecha |
| Choque fuerte | `strong_interface`, densidades y presiones izquierda/derecha |
| Shu--Osher | `shu_interface`, estados izquierdo/derecho, amplitud y número de onda |
| Kelvin--Helmholtz | `khi_half_width`, densidades, velocidades, presión, amplitud y número de onda |
| Jet | estado ambiente, geometría de la boquilla y estado del jet |
| Convergencia | densidad, amplitud, velocidad, presión y número de onda |
| Michel | radio y densidad críticos |
| Polvo | `dust_accretion_constant` |
| Explosión fuera del eje | centro, ancho, fondo y amplitudes |
| Fishbone--Moncrief | `fm_inner_radius`, `fm_pressure_max_radius`, `fm_polytropic_constant` |

Los parámetros de una condición inicial distinta a `problem` son rechazados;
esto evita configuraciones híbridas accidentales.

## 5. Reinicio desde checkpoint

Los checkpoints actuales usan el identificador `GRHDCP_V3` y contienen malla,
física, métodos, estado conservado, tiempo, paso y estado de la perturbación.

Ejemplo mínimo:

```text
restart = true
restart_file = FM_KS_a0.9_data/checkpoint_weno5_hlle_00500.rst

final_time = 1000.0
cfl = 0.4
save_interval = 10.0
checkpoint_interval = 100.0
output_prefix = FM_KS_a0.9_restart
output_folder = FM_KS_a0.9_restart_data
vtk_mapping = physical

shock_sensor = false
extract_gw = false
extract_mdot = false
ppi_diagnostics = true
diagnostic_stride = 100
perturbation = none
perturbation_seed = 3435
perturbation_time = 1000.0
perturbation_amplitude = 0.01
perturbation_mode = 1.0
```

En un reinicio sólo pueden cambiarse controles de tiempo, salida, diagnósticos,
perturbación y sensor. La malla, métrica, EOS y arquitectura numérica se leen
del checkpoint y el lector rechaza intentos de modificarlas.

Antes de reiniciar, valide conjuntamente parámetros y checkpoint:

```bash
./grhd2 --check reinicio.par
```

Los checkpoints son volcados binarios Fortran; deben tratarse como ligados a la
versión del formato, precisión, compilador/arquitectura y forma del estado.
Conserve siempre el código y manifiesto que los generaron.

## 6. Salidas

La ruta normal es:

```text
output_folder/scheme_name/
```

donde `scheme_name` combina reconstructor y Riemann, por ejemplo
`weno5_hlle`.

| Archivo | Contenido |
|---|---|
| `run_manifest.txt` | Configuración efectiva y metadatos físicos |
| `PREFIX_SCHEME_NNNN.vtk` | Estado para VisIt/ParaView |
| `checkpoint_SCHEME_TTTTT.rst` | Reinicio binario |
| `m_dot.dat` | `t`, tasa de acreción |
| `ppi_modes.dat` | masa y modos complejos `m=1,...,4` |
| `global_diagnostics.dat` | masas, momento angular, extremos, `v^2` y conteos |
| `perturbation_events.dat` | instante y parámetros de cada inyección |
| `GW_signal.dat` | `t`, `h_plus`, `h_cross` del proxy Finn--Evans |
| `convergence_*.dat` | normas para onda advectada o Michel |

El manifiesto se añade al reiniciar, delimitando cada sesión con `[run]` y
`[/run]`.

### 6.1 Interpretación de diagnósticos

La tasa de acreción se integra en la primera celda exterior al horizonte con
la velocidad coordenada efectiva `alpha*v^r - beta^r`.

Los modos PPI se definen a partir de la masa azimutal:

```text
C_m = integral D exp(-i m phi) d^3x
A_m = |C_m| / C_0
```

Se reportan partes real e imaginaria, amplitud normalizada y fase para
`m=1,...,4`.

En `global_diagnostics.dat`, `N_invalid=0` es necesario pero no suficiente para
validar una corrida. También deben examinarse conservación, sensibilidad a
resolución, floors, fronteras y método.

Con `ny=1`, las masas y amplitudes integradas heredan una cuadratura polar de
una sola celda. No deben interpretarse como integrales tridimensionales
resueltas.

## 7. Visualización y análisis

Para crear perfiles de Fishbone--Moncrief:

```bash
python3 plot_fm_profiles.py
```

Con un patrón explícito:

```bash
python3 plot_fm_profiles.py \
  'FM_KS_a0.9_data/weno5_hlle/*.vtk' \
  --field Density --frames 100 \
  --output figures/fm_a09_density.png --no-show
```

Consulte todas las opciones con:

```bash
python3 plot_fm_profiles.py --help
python3 plot_advected_wave_convergence.py --help
```

`gw_fm.py` y `gw_qtrans_fm.py` son utilidades históricas de análisis. Antes de
usarlas en resultados de tesis debe comprobarse que sus convenciones coincidan
con el manifiesto y formato de la campaña concreta.

## 8. Pruebas disponibles

Ejecute la regresión completa versionada con:

```bash
make test-hllc
make test-gw
make test-fm
make test-polar
make test-cfl
```

Cobertura actual:

- `test-hllc`: estados uniformes, contactos, métrica general, atmósfera,
  invariancia de escala y salto relativista fuerte;
- `test-gw`: radio del horizonte, transformaciones KS, cinemática cartesiana y
  tensor de esfuerzo Finn--Evans;
- `test-fm`: propiedades puntuales del dato inicial Fishbone--Moncrief en EF y
  KS;
- `test-polar`: periodicidad azimutal, cruce transpolar y paridades;
- `test-cfl`: condición CFL aditiva en 1D, 2D y 3D con *shift*.

Estas pruebas no sustituyen campañas de convergencia ni validación física.
Antes de publicar resultados deben registrarse compilador, banderas, commit,
archivo `.par`, número de hilos y salida completa de las pruebas.

## 9. Arquitectura del código

| Archivo | Responsabilidad |
|---|---|
| `variables.f90` | Estado global, índices, constantes y configuración |
| `parameters.f90` | Lectura, tipado y validación de `.par` |
| `metrics.f90` | Interfaz 3+1, métricas, derivadas y Christoffel |
| `conditions.f90` | Datos iniciales, fronteras y perturbaciones |
| `equations.f90` | Valencia, flujos físicos, fuentes y C2P |
| `fluxes.f90` | Velocidades características, HLLE y HLLC |
| `reconstruction.f90` | Godunov, TVD, WENO3, MP5 y WENO5-Z |
| `evolution.f90` | Reconstrucción por dirección y suma del RHS |
| `output.f90` | VTK, checkpoints, manifiesto y diagnósticos |
| `initialization.f90` | Presets, malla y caché geométrica |
| `grhd2.f90` | Arranque, SSP-RK3 y bucle principal |

Flujo de una corrida:

```text
leer cabecera .par
  -> cargar preset
  -> aplicar .par y validar
  -> seleccionar métrica y métodos
  -> construir malla y caché geométrica
  -> inicializar primitivas y conservadas
  -> calcular CFL
  -> repetir SSP-RK3
       -> fronteras
       -> reconstrucción por cada dirección activa
       -> Riemann en caras
       -> sumar divergencias y fuentes en un RHS
       -> recuperar primitivas
  -> diagnósticos, VTK y checkpoints
```

El algoritmo reconstruye tiras por dirección, pero no aplica operadores
temporales separados: las contribuciones espaciales se suman antes de cada
etapa RK. Por ello no debe documentarse como separación de Strang.

## 10. Cómo extender el código

### 10.1 Añadir una condición inicial

1. Defina constantes y parámetros en `variables.f90`.
2. Añada el nombre en `parse_problem` y `problem_label`.
3. Declare y valide sus claves en `parameters.f90`.
4. Cree el preset en `initialization.f90`.
5. Implemente el dato inicial y las fronteras en `conditions.f90`.
6. Añada al menos una prueba unitaria y un `.par` pequeño reproducible.
7. Documente estado, tolerancias y resultado en `docs/DEVELOPMENT_LOG.md`.

### 10.2 Añadir una métrica

La métrica debe proporcionar, en coordenadas físicas y si procede
logarítmicas:

- `alpha`, `beta^i`, `gamma_ij` y `sqrt(gamma)`;
- inversa espaciotemporal cuando se solicite;
- derivadas espaciales de `g_mu_nu`;
- símbolos de Christoffel;
- derivadas de `ln(alpha)`.

Conecte las rutinas mediante `set_metric_type`. Pruebe determinante, inversa,
simetrías, derivadas, horizonte, transformación de coordenadas y el límite a
una métrica ya conocida.

### 10.3 Añadir un reconstructor o Riemann solver

1. Añada un identificador y parser explícito.
2. Defina el número de celdas fantasma requerido.
3. Mantenga separadas primitivas, conservadas no densitizadas y densitizadas.
4. Añada pruebas de estado uniforme, contactos, choques y escalamiento.
5. Compruebe positividad, causalidad, fallback y convergencia formal.

### 10.4 Camino hacia GRMHD

No basta con aumentar `neq`. Antes de activar campos magnéticos deben fijarse:

- primitivas y conservadas GRMHD;
- tensor de energía--momento magnetizado;
- inducción y control de `div B` mediante GLM o CT;
- velocidades magnetosónicas y solver compatible;
- recuperación de primitivas GRMHD;
- fronteras y topología polar magnéticas;
- límite obligatorio `B -> 0` contra este baseline GRHD.

## 11. Reproducibilidad mínima

Cada resultado que se conserve debe incluir:

- hash del commit y estado del árbol de trabajo;
- copia exacta del archivo `.par`;
- compilador, versión y `FFLAGS`;
- número de hilos OpenMP y plataforma;
- `run_manifest.txt`;
- log de ejecución y código de salida;
- hashes de fuente, ejecutable y datos principales;
- scripts de análisis y versión de dependencias;
- explicación de cualquier reinicio o reparación de celda.

Los directorios de producción, VTK, checkpoints y logs están excluidos de Git.
Conserve datos grandes fuera del repositorio junto con un manifiesto trazable.

## 12. Límites de interpretación científica

- Fondo fijo/Cowling: no hay respuesta de la métrica al fluido.
- Sin autogravedad: se requiere una razón `M_torus/M_BH` físicamente pequeña
  para interpretar directamente el modelo de fluido de prueba.
- Sin MHD: no aparecen MRI, campos magnéticos ni transporte magnético.
- Una celda polar: no resuelve espesor vertical, *warping* ni modos 3D.
- Finn--Evans: proporciona un proxy de frecuencia/fase bajo una misma
  normalización, no amplitud absoluta observable.
- Barrido de espín histórico: cambian masa y densidad iniciales; no aísla el
  efecto causal del espín.
- Una sola semilla o resolución no establece incertidumbre estadística ni
  convergencia.

## 13. Solución de problemas

### El programa imprime el uso y termina

Debe proporcionar un archivo:

```bash
./grhd2 archivo.par
```

### Se rechaza una clave válida para otro problema

Los parámetros de condición inicial se restringen al `problem` seleccionado.
Elimine la clave ajena o cambie el problema.

### La malla radial logarítmica falla

Compruebe `logarithmic_r = true` y `r_min > 0`.

### EF rechaza el espín

La implementación EF representa Schwarzschild; use `spin = 0`. Para Kerr use
`metric = kerr_schild` con `|spin| <= bh_mass`.

### Un toro 3D aparece vacío

Verifique que los centros polares realmente muestreen el cuerpo del toro. Un
número pequeño de celdas en todo `[0,pi]` puede no colocar ninguna celda dentro
de su espesor angular. Inspeccione el VTK de `final_time = 0` antes de lanzar
una campaña larga.

### `GW_signal.dat` parece tener amplitud extraña

Compruebe `ny`, distancia registrada, sustracción de atmósfera, horizonte,
mapeo y normalización. Con `ny=1` la amplitud absoluta no tiene normalización
3D resuelta.

### Una corrida termina sin `N_invalid`

Ese contador sólo existe cuando `ppi_diagnostics = true`. Su ausencia no
demuestra que la corrida sea válida o inválida.

## 14. Documentación relacionada

- `docs/DEVELOPMENT_LOG.md`: bitácora técnica verificable.
- `docs/teoria/teoria_desarrollo_grhd.tex`: teoría, formulación e historia de
  desarrollo en formato LaTeX.
- `CONTEXTO_GRHD_PARA_PRISM.md`: auditoría científica de campañas anteriores.
- `extraccion_gw_gauge_invariant.md`: ruta desde el proxy actual hacia métodos
  perturbativos y extracción relativista.
- `par/README.txt`: referencia compacta de los parámetros por problema.

Para compilar el documento teórico y resolver sus referencias:

```bash
cd docs/teoria
pdflatex -interaction=nonstopmode -halt-on-error teoria_desarrollo_grhd.tex
bibtex teoria_desarrollo_grhd
pdflatex -interaction=nonstopmode -halt-on-error teoria_desarrollo_grhd.tex
pdflatex -interaction=nonstopmode -halt-on-error teoria_desarrollo_grhd.tex
```
