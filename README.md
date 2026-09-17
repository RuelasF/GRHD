# GRHD2

Este README funciona como manual de usuario y desarrollo de GRHD2. Describe cómo compilar, configurar, ejecutar, verificar y extender el código contenido en este repositorio; su fuente son las rutinas Fortran y las pruebas versionadas en `tests/`.

> **Estado físico:** GRHD para un fluido perfecto sobre una métrica analítica fija. No es GRMHD, no evoluciona el espacio-tiempo y no incluye autogravedad del fluido. La salida `GW_signal.dat` es un proxy cuadrupolar Finn--Evans de campo débil, no un *strain* observable ni una extracción gauge-invariant.

> **Estado de entrega, 15 de septiembre de 2026:** la base GRHD queda cerrada como baseline funcional para iniciar GRMHD. Compilan el ejecutable y la suite `make test`; los doce problemas cuentan con corridas de validación, nueve también llegaron a sus tiempos de escala nominal, la matriz de métodos pasó 168/168 y una comparación Fishbone--Moncrief ecuatorial, sagital y 3D con el mismo toro y una órbita terminó 3/3 sin estados inválidos. Esto verifica la implementación y la estabilidad dentro de ese alcance; no equivale a convergencia 3D, validación de PPI o fidelidad astrofísica de una familia de toros.

## 1. Requisitos

El flujo probado utiliza:

- GNU Fortran (`gfortran`);
- GNU Make;
- OpenMP, suministrado por `gfortran`;
- VisIt o ParaView, opcionales, para inspeccionar los VTK.

En Debian/Ubuntu puede instalarse la base con:

```bash
sudo apt update
sudo apt install gfortran make
```

La instalación de paquetes es una operación del sistema y no forma parte del `Makefile`.

## 2. Compilación rápida

Desde la raíz del repositorio:

```bash
make -j
```

El ejecutable resultante es `grhd2`. Las banderas de producción actuales son `-O3`, OpenMP, optimización específica de la máquina y LTO. Para una auditoría de memoria o punto flotante conviene reemplazar temporalmente `FFLAGS` por las banderas de depuración comentadas en el `Makefile`.

### 2.1 OpenACC en NVIDIA P100

La rama `openacc-p100-soa` incluye un backend separado para GPU. Requiere NVIDIA HPC SDK y `nvfortran`; no utiliza OpenMP durante la evolución GPU. El estado hidrodinámico usa el orden `[x,y,z,variable]`, permanece residente en el dispositivo entre etapas RK y todos los estados, flujos, fuentes y conversiones usan `real*8` (FP64). No se habilitan `fastmath`, TF32, FP32 ni precisiones mixtas.

Para Tesla P100 (Pascal GP100, capacidad de cómputo 6.0):

```bash
make openacc-p100
./grhd2 --check ejemplo.par
CUDA_VISIBLE_DEVICES=0 ./grhd2 ejemplo.par
```

Para una prueba funcional en una RTX 4070 de Zotz (Ada, capacidad 8.9):

```bash
make openacc-rtx4070
CUDA_VISIBLE_DEVICES=0 ./grhd2 ejemplo.par
```

La RTX 4070 admite instrucciones FP64 y permite comprobar que los kernels funcionan, pero su hardware está orientado a FP32/IA y no representa el rendimiento FP64 de una P100. Una GTX 1050 Ti también ejecuta CUDA y FP64, pero su capacidad de cómputo es 6.1 y su rendimiento de doble precisión es demasiado bajo para usarla como referencia de desempeño. La validación final del target `openacc-p100` debe hacerse sobre una P100 con una versión de NVIDIA HPC SDK que aún genere código `cc60`.

Para confirmar qué GPU tomó el proceso y observar memoria/ocupación:

```bash
nvidia-smi -L
watch -n 1 nvidia-smi
```

`OMP_NUM_THREADS` no controla esta ruta. Las transferencias dispositivo--host se reservan para VTK, checkpoints, diagnósticos y la perturbación diferida; reducir su frecuencia disminuye el tráfico PCIe.

#### Comparación reproducible CPU--P100 con el jet

La rama incluye un benchmark del jet axisimétrico con la malla nominal `320 x 1 x 800`, WENO5--HLLE, FP64 y `shock_sensor = true`. Ambos ejecutables reciben el mismo archivo de parámetros y evolucionan hasta `t=2`; la prueba mide el tiempo de pared y compara todos los campos del VTK final mediante normas L1, L2, Linf y L2 relativa.

En una máquina con `gfortran`, NVIDIA HPC SDK, controlador NVIDIA y una P100 visible, basta ejecutar:

```bash
git clone --branch openacc-p100-soa --single-branch git@github.com:RuelasF/GRHD.git
cd GRHD
CPU_THREADS=20 GPU_ID=0 make benchmark-jet-p100
```

Los resultados quedan en `.benchmark/jet_cpu_p100_FECHA/`: `system.txt` registra hardware, compiladores y commit; `cpu/` y `gpu/` contienen logs, tiempos y VTK; `comparison.txt` contiene las diferencias numéricas y `timing_summary.txt` la aceleración CPU/GPU. `CPU_THREADS` debe representar la referencia que se desea comparar y no controla el cálculo OpenACC. Para verificar sólo la ruta CPU puede ejecutarse `bash benchmarks/jet/run_cpu_p100.sh --cpu-only`.

Para eliminar objetos, módulos, ejecutable y binarios de prueba:

```bash
make clean
```

## 3. Primer cálculo

El repositorio no versiona configuraciones de producción. Cree localmente un archivo `ejemplo.par`; una configuración mínima que aprovecha el *preset* del problema puede comenzar así:

```text
problem = fishbone_equatorial
metric = kerr_schild
geometry = spheroidal
logarithmic_r = true
spin = 0.9
reconstruction = weno5
riemann_solver = hlle
final_time = 1000.0
```

Compruebe el archivo sin iniciar la simulación:

```bash
./grhd2 --check ejemplo.par
```

Si la validación termina correctamente, ejecute:

```bash
export OMP_NUM_THREADS=8
./grhd2 ejemplo.par
```

`OMP_NUM_THREADS` debe ajustarse al equipo. El programa exige exactamente un archivo `.par`; no permite una ejecución válida sin argumentos.

Si la salida estándar se redirige a un `.log` y se quiere ver el progreso en tiempo real, establezca `GFORTRAN_UNBUFFERED_ALL=y` antes de ejecutar. Esto afecta el búfer de mensajes de GNU Fortran, no los datos ni la evolución:

```bash
export GFORTRAN_UNBUFFERED_ALL=y
./grhd2 ejemplo.par >ejemplo.log 2>&1
```

El ejemplo genera un toro de Fishbone--Moncrief ecuatorial en Kerr--Schild con `a=0.9`, malla `400 x 1 x 200`, WENO5--HLLE y radio logarítmico. La salida se guarda, por defecto, en:

```text
FM_KS_a0.9_data/weno5_hlle/
```

Aunque se evolucionan las direcciones radial y azimutal, una sola celda polar no resuelve la estructura vertical. Esta configuración debe describirse como **modelo ecuatorial 2.5D**, no como simulación 3D.

## 4. Archivos de parámetros

### 4.1 Sintaxis

Cada línea activa tiene la forma:

```text
nombre = valor       # comentario opcional
```

Las claves y valores simbólicos no distinguen mayúsculas de minúsculas. Las rutas conservan su escritura. Una clave desconocida, duplicada, incompatible con el problema o con valor no admisible detiene el programa con un mensaje que indica archivo y línea.

Para dominios `spherical` o `spheroidal`, `theta_min`, `theta_max`, `phi_min` y `phi_max` se escriben como múltiplos de pi:

```text
theta_min = 0
theta_max = 1       # pi
phi_min   = 0
phi_max   = 2       # 2*pi
```

Se admiten cocientes, por ejemplo `theta_min = 1/2`.

Las claves internas `y_min`, `y_max`, `z_min` y `z_max` se leen directamente en unidades de coordenada, sin multiplicarlas por pi. En el jet cilíndrico, `y` representa el ángulo azimutal en radianes; su preset usa `ny=1`, `y_min=0` y `y_max=pi`, de modo que el corte meridional se centra en `phi=pi/2`. Si se modifica desde un `.par`, use valores numéricos en radianes para `y_min/y_max`; `theta_min/theta_max` son coeficientes de pi y se convierten al leer el archivo.

### 4.2 Problemas disponibles

| `problem` | Descripción | Geometría usual | Evidencia disponible al 15-09-2026 |
|---|---|---|---|
| `sod` | Tubo de choque relativista | Minkowski cartesiana 1D | WENO5, HLLE/HLLC, cuatro resoluciones y comparación L1 contra referencia fina; escala nominal completada |
| `strong_shock` | Choque fuerte | Minkowski cartesiana 1D | WENO5, HLLE/HLLC, cuatro resoluciones y comparación L1 contra referencia fina |
| `shu_osher` | Interacción choque--onda | Minkowski cartesiana 1D | WENO5, HLLE/HLLC, cuatro resoluciones y error L1 decreciente |
| `kelvin_helmholtz` | Inestabilidad de cizalla | Minkowski cartesiana 2D | Escala nominal hasta `t=3`; conservación de masa y ausencia de estados inválidos; convergencia física pendiente |
| `axisymmetric_jet` | Jet relativista axisimétrico | Minkowski cilíndrica | Escala nominal hasta `t=120`; orientación VTK corregida y verificada; comparación cuantitativa externa pendiente |
| `convergence` | Onda sinusoidal advectada | Minkowski cartesiana 1D | Orden observado medido para Godunov, TVD, WENO3, MP5 y WENO5 |
| `michel` | Acreción esférica estacionaria | EF esférica | Error contra solución analítica decreciente al refinar; escala nominal hasta `t=500` |
| `dust` | Acreción fría desde reposo en el infinito | EF esférica | Flujo/energía puntuales y error contra solución analítica; escala nominal hasta `t=500` |
| `offaxis_blast` | Explosión fuera del eje | EF esférica | Comparación de resoluciones y escala nominal hasta `t=2000`; dominio abierto, no masa cerrada |
| `fishbone_equatorial` | Toro FM en `r-phi` | EF esférica o KS esferoidal | Escala nominal hasta `t=1000`; toro EF emparejado hasta una órbita; PPI no queda validada por ello |
| `fishbone_sagittal` | Corte FM en `r-theta` | EF esférica | Escala nominal hasta `t=100`; toro EF emparejado hasta una órbita |
| `fishbone_3d` | Toro FM con tres direcciones activas | EF esférica o KS esferoidal | Escala nominal hasta `t=35`; toro EF emparejado hasta una órbita con VTK y diagnósticos auditados; convergencia 3D pendiente |

Cada problema carga un *preset* desde `initialization.f90`; el archivo `.par` sobrescribe después los valores especificados. Es aconsejable declarar de forma explícita todos los parámetros relevantes para una corrida científica.

### 4.3 Métricas y geometrías

Valores admitidos:

```text
metric   = minkowski | eddington_finkelstein | kerr_schild
geometry = cartesian | cylindrical | spherical | spheroidal
```

Compatibilidad actual:

- Minkowski: cartesiana, cilíndrica o esférica, con `bh_mass = 0` y `spin = 0`;
- Eddington--Finkelstein: esférica, Schwarzschild y `spin = 0`;
- Kerr--Schild: esferoidal y `|spin| <= bh_mass`.

Las geometrías esférica y esferoidal requieren `r_min > 0` porque la implementación no incluye una regularización del origen. En coordenada radial logarítmica se usa `x = ln(r)`; esta opción sólo está disponible para geometrías esférica y esferoidal:

```text
logarithmic_r = true
```

La opción esférica plana está implementada de forma explícita y puede seleccionarse con `metric = minkowski`, `geometry = spherical` y `bh_mass = 0`. Incluye tanto la coordenada física `r` como la coordenada lógica `x = ln(r)`; aunque el tensor de Riemann es nulo, los símbolos de Christoffel coordenados y los términos fuente geométricos no lo son.

En la geometría cilíndrica se usa el orden lógico `(x1,x2,x3) = (rho,phi,z)`. El caso de jet axisimétrico mantiene `ny = 1`, por lo que `phi` es la dirección degenerada y `z` es la dirección axial activa; la cara `rho = 0` tiene flujo densitizado nulo y se trata explícitamente sin invocar el solucionador de Riemann sobre una métrica singular.

El corte VTK del preset del jet está centrado en `phi=pi/2` y yace en el plano `Y-Z`. Los VTK generados antes de la corrección de septiembre de 2026 usan `phi=0.5` rad y se ven inclinados; la física axisimétrica no cambia por una rotación del plano. Para comparar imágenes históricas, rote esos VTK `pi/2-0.5` rad alrededor de `Z` o reproyéctelos antes de graficar.

Kerr--Schild utiliza `geometry = spheroidal` porque sus superficies de radio constante son esferoides cuando `spin` es distinto de cero. Por compatibilidad, un archivo o checkpoint antiguo que declare `geometry = spherical` junto con Kerr--Schild se acepta temporalmente, emite una advertencia y se normaliza internamente a `Spheroidal`.

### 4.4 Métodos numéricos

```text
reconstruction = godunov | tvd | weno3 | mp5 | weno5
tvd_limiter    = minmod | superbee | mc
riemann_solver = hlle | hllc
shock_sensor   = true | false
```

La profundidad de celdas fantasma se selecciona automáticamente: una para Godunov, dos para TVD/WENO3 y tres para MP5/WENO5. `HLLD` no es una opción válida: sólo existe un esqueleto que aborta si se invoca internamente.

HLLE es el *baseline* robusto. HLLC está implementado para GRHD mediante un marco ortonormal local de cada cara y vuelve a HLLE ante estados degenerados o no admisibles.

El sensor de choques puede degradar localmente la reconstrucción hacia métodos más robustos. Desactivarlo permite comparar reconstructores puros, pero no es una decisión universalmente más estable.

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
- `untwisted`: representación oblata sin la rotación radial, sólo para visualización.

Esta opción no cambia el estado ni las coordenadas de evolución.

### 4.6 Diagnósticos y perturbaciones

```text
extract_gw = false
extract_mdot = false
ppi_diagnostics = true
diagnostic_stride = 100
extraction_stride = 10
mass_monitor_stride = 1000

perturbation = none | pressure_noise | density_noise | density_mode
perturbation_seed = 3435
perturbation_time = 1000.0
perturbation_amplitude = 0.01
perturbation_mode = 1.0
```

`diagnostic_stride` controla la cadencia de `ppi_modes.dat` y `global_diagnostics.dat`; `extraction_stride` controla `GW_signal.dat` y `m_dot.dat`. `mass_monitor_stride` sólo controla el cálculo e impreso periódico de la masa en la terminal y puede fijarse en `0` para desactivarlo sin alterar los diagnósticos guardados.

`pressure_noise` y `density_noise` usan una semilla reproducible. `density_mode` inyecta un modo azimutal coherente cuyo número es `perturbation_mode`. La amplitud debe satisfacer `0 <= perturbation_amplitude < 1`.

Una perturbación solicitada se aplica al comienzo del primer paso para el cual el tiempo acumulado ya alcanzó `perturbation_time`; el instante efectivo se registra en `perturbation_events.dat`.

### 4.7 Parámetros específicos por problema

Las claves específicas se validan en `parameters.f90`; sus grupos principales son:

| Problema | Claves específicas |
|---|---|
| Sod | `sod_interface`, densidades y presiones izquierda/derecha |
| Choque fuerte | `strong_interface`, densidades y presiones izquierda/derecha |
| Shu--Osher | `shu_interface`, estados izquierdo/derecho, amplitud y número de onda |
| Kelvin--Helmholtz | `khi_half_width`, densidades, velocidades, presión, amplitud, número de onda y `khi_perturbation_width` |
| Jet | estado ambiente, geometría de la boquilla y estado del jet |
| Convergencia | densidad, amplitud, velocidad, presión y número de onda |
| Michel | radio y densidad críticos |
| Polvo | `dust_accretion_constant` |
| Explosión fuera del eje | centro, ancho, fondo y amplitudes |
| Fishbone--Moncrief | `fm_inner_radius`, `fm_pressure_max_radius`, `fm_polytropic_constant` |

Los parámetros de una condición inicial distinta a `problem` son rechazados; esto evita configuraciones híbridas accidentales.

## 5. Reinicio desde checkpoint

Los checkpoints actuales usan el identificador `GRHDCP_V5` y contienen malla, física, métodos, parámetros específicos de la condición inicial, estado conservado, tiempo, paso y estado de la perturbación. V5 registra el estado con el orden `[x,y,z,variable]`; el lector conserva compatibilidad con V2--V4 y transpone automáticamente sus estados `[variable,x,y,z]`. Los parámetros específicos permiten reconstruir exactamente las fronteras de Jet, Michel y Dust después de un reinicio.

Ejemplo mínimo:

```text
restart = true
restart_file = FM_KS_a0.9_data/checkpoint_weno5_hlle_step_000123456.rst

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

En un reinicio sólo pueden cambiarse controles de tiempo, salida, diagnósticos, perturbación y sensor. La malla, métrica, EOS y arquitectura numérica se leen del checkpoint y el lector rechaza intentos de modificarlas.

Antes de reiniciar, valide conjuntamente parámetros y checkpoint:

```bash
./grhd2 --check reinicio.par
```

Los checkpoints son volcados binarios Fortran; deben tratarse como ligados a la versión del formato, precisión, compilador/arquitectura y forma del estado. Conserve siempre el código y manifiesto que los generaron.

## 6. Salidas

La ruta normal es:

```text
output_folder/scheme_name/
```

donde `scheme_name` combina reconstructor y Riemann, por ejemplo `weno5_hlle`.

| Archivo | Contenido |
|---|---|
| `run_manifest.txt` | Configuración efectiva y metadatos físicos |
| `PREFIX_SCHEME_NNNN.vtk` | Estado para VisIt/ParaView |
| `checkpoint_SCHEME_step_NNNNNNNNN.rst` | Reinicio binario identificado de forma única por el paso |
| `m_dot.dat` | `t`, tasa de acreción |
| `ppi_modes.dat` | masa y modos complejos `m=1,...,4` |
| `global_diagnostics.dat` | masas, momento angular, extremos, `v^2` y conteos |
| `perturbation_events.dat` | instante y parámetros de cada inyección |
| `GW_signal.dat` | `t`, `h_plus`, `h_cross` del proxy Finn--Evans |
| `convergence_*.dat`, `michel_convergence_*.dat`, `dust_convergence_*.dat` | normas de error para la onda advectada y las soluciones estacionarias de acreción |

El manifiesto se añade al reiniciar, delimitando cada sesión con `[run]` y `[/run]`.

### 6.1 Interpretación de diagnósticos

La tasa de acreción se integra en la primera celda exterior al horizonte con la velocidad coordenada efectiva `alpha*v^r - beta^r`.

Las fronteras radiales de salida en EF y Kerr--Schild aplican el mismo criterio al impedir entrada: el umbral de flujo coordenado nulo es `v^r = beta^r/alpha`, no `v^r = 0`. Esta distinción evita inyección espuria de masa cuando el *shift* radial es distinto de cero.

Los modos PPI se definen a partir de la masa azimutal:

```text
C_m = integral D exp(-i m phi) d^3x
A_m = |C_m| / C_0
```

Se reportan partes real e imaginaria, amplitud normalizada y fase para `m=1,...,4`.

En `global_diagnostics.dat`, `N_invalid=0` es necesario pero no suficiente para validar una corrida. También deben examinarse conservación, sensibilidad a resolución, floors, fronteras y método.

Con `ny=1`, las masas y amplitudes integradas heredan una cuadratura polar de una sola celda. No deben interpretarse como integrales tridimensionales resueltas.

### 6.2 Campaña PPI `inundaciones`

`inundaciones.sh` automatiza la campaña ecuatorial destinada al congreso: EF con `a=0` y KS con `a=0,0.2,0.4,0.6,0.9`, cada fondo con MP5--HLLE y WENO5--HLLE. Todos los casos usan `800x1x400`, `t_final=5000 M`, ruido de presión de `1%` con semilla 3435 en `t=1000 M`, salidas VTK cada `50 M` y checkpoints cada `250 M`. El sensor de choques se desactiva para comparar los reconstructores sin degradación adaptativa.

El controlador está diseñado para recibir el ejecutable `grhd2` ya compilado junto con el script; no requiere Git, Make ni un compilador en el servidor. Copia el binario como `inundaciones`, registra su SHA-256, genera y valida los doce `.par` con `--check` y sólo entonces inicia las simulaciones. Ejecuta seis casos simultáneos con un presupuesto total de 72 hilos; si el sistema expone menos de 72 pero al menos 60, usa 60. Para transferir los archivos y lanzar una ejecución que sobreviva al cierre de la terminal:

```bash
scp grhd2 inundaciones.sh usuario@servidor:/ruta/campana/
# Ya en /ruta/campana del servidor:
bash inundaciones.sh --detach
```

Si el binario tiene otro nombre o ubicación puede indicarse con `GRHD_EXECUTABLE=/ruta/al/binario`. Al reanudar, el script exige el mismo SHA-256 para impedir que una campaña mezcle ejecutables distintos.

También puede enviarse con `sbatch inundaciones.sh` o `qsub inundaciones.sh` si el clúster usa uno de esos planificadores. El estado general se consulta con:

```bash
tail -f ../inundaciones_ppi_800_t5000/logs/inundaciones.driver.log
```

Una segunda ejecución con la misma `CAMPAIGN_ROOT` omite casos completos y crea una nueva sesión para cada caso incompleto a partir del checkpoint más reciente. Los diagnósticos válidos de las sesiones se consolidan en `results/CASO/combined/`, mientras que los datos originales permanecen separados para conservar procedencia. Con 101 VTK por caso, la campaña completa requiere aproximadamente 35 GB entre campos, checkpoints y archivos auxiliares; debe comprobarse el espacio disponible antes de iniciarla.

Este barrido conserva `r_in`, `r_pmax` y `K`, no la masa o la densidad máxima inicial del toro. Por ello sirve para comparar las familias configuradas y estudiar señales PPI, pero las diferencias entre espines no deben presentarse como efecto causal aislado de `a` sin una campaña adicional que empareje propiedades físicas del disco. `GW_signal.dat` sigue siendo un proxy Finn--Evans con una sola celda polar, no una amplitud tridimensional observable.

El lanzador `inundaciones_weno3.sh` reproduce los mismos seis fondos exclusivamente con WENO3--HLLE en otro directorio, `inundaciones_ppi_weno3_800_t5000`. Está configurado para usar exactamente 27 hilos: ejecuta tres casos simultáneos con 9 hilos cada uno y completa la matriz en dos tandas. Conserva la misma malla, perturbación, tiempos, salidas y reglas de reanudación. Como este lanzador reutiliza el controlador general, al segundo servidor se transfieren el ejecutable y ambos scripts:

```bash
scp grhd2 inundaciones.sh inundaciones_weno3.sh usuario@servidor:/ruta/campana/
# Ya en /ruta/campana del servidor:
bash inundaciones_weno3.sh --detach
```

Su log general queda en `../inundaciones_ppi_weno3_800_t5000/logs/inundaciones.driver.log` y el volumen esperado de resultados es cercano a 18 GB. La campaña histórica mostró que WENO3 podía sembrar modos antes de aplicar la perturbación explícita; por ello esta nueva serie es un control numérico útil, pero el crecimiento anterior a `t=1000 M` debe analizarse y reportarse por separado, no atribuirse a la perturbación programada.

## 7. Visualización y análisis

Las salidas VTK pueden inspeccionarse con VisIt o ParaView. Salvo el controlador reproducible `inundaciones.sh`, los scripts de análisis y las configuraciones exploratorias se mantienen fuera de este repositorio mínimo; cualquier análisis destinado a una publicación debe conservar externamente su versión, dependencias y correspondencia con el manifiesto de la corrida.

## 8. Pruebas disponibles

Ejecute la regresión completa versionada con:

```bash
make test
```

Los objetivos individuales `test-hllc`, `test-gw`, `test-fm`, `test-polar`, `test-cfl`, `test-minkowski-spherical`, `test-cases`, `test-checkpoint` y `test-output` permiten aislar una familia concreta.

Cobertura actual:

- `test-hllc`: estados uniformes, contactos, métrica general, atmósfera, invariancia de escala y salto relativista fuerte;
- `test-gw`: radio del horizonte, transformaciones KS, cinemática cartesiana y tensor de esfuerzo Finn--Evans;
- `test-fm`: propiedades puntuales del dato inicial Fishbone--Moncrief en EF y KS;
- `test-polar`: periodicidad azimutal, cruce transpolar y paridades;
- `test-cfl`: condición CFL aditiva en 1D, 2D y 3D con *shift*;
- `test-minkowski-spherical`: forma analítica, equivalencia con EF para `bh_mass = 0` y curvatura nula en coordenadas radial física y logarítmica;
- `test-cases`: métrica cilíndrica `(rho,phi,z)`, perturbación transversal de Kelvin--Helmholtz, radio físico de OffAxis en malla logarítmica, solución/fronteras de acreción fría y salida radial compatible con *shift*.
- `test-checkpoint`: escritura/lectura V5, conservación del estado y de los parámetros específicos necesarios para reproducir datos iniciales y fronteras.
- `test-acc-reconstruction`: equivalencia de la reconstrucción por caras OpenACC con la ruta de tiras CPU para los cinco métodos, con y sin sensor.
- `test-acc-rhs`: equivalencia del RHS CPU/OpenACC para los cinco reconstructores, HLLE/HLLC y geometría plana/curva.
- `test-output`: construcción segura de nombres VTK largos dentro de los límites públicos de configuración.

Estas pruebas no sustituyen campañas de convergencia ni validación física. Antes de publicar resultados deben registrarse compilador, banderas, commit, archivo `.par`, número de hilos y salida completa de las pruebas.

### 8.1 Cierre funcional GRHD

Las campañas se ejecutaron con una copia aislada del código en Zotz y resultados fuera de Git. La matriz de compatibilidad cubrió 12 problemas, siete opciones efectivas de reconstrucción y dos solucionadores de Riemann: 168/168 corridas cortas completadas. La onda suave permitió medir orden observado; Michel y polvo se compararon contra soluciones estacionarias; Sod, choque fuerte y Shu--Osher se refinaron de `N=200` a `N=1600`; KHI, jet, OffAxis y las tres representaciones FM pasaron pruebas multidimensionales y corridas a escala nominal. La tabla y trazas locales están en `.validation/`, excluido del repositorio público.

La prueba final emparejada usó un toro Fishbone--Moncrief no perturbado sobre Schwarzschild en Eddington--Finkelstein con `M=1`, `r_in=6.25`, `r_pmax=9.25`, `K=0.0015`, `Gamma=4/3`, radio logarítmico `1.2<=r<=40`, WENO5--HLLE, sensor de choques activo y `CFL=0.4`. Las mallas fueron `128x1x128` para el corte ecuatorial, `128x64x1` para el sagital y `128x64x128` para 3D; las tres evolucionaron hasta `t=180 M`, aproximadamente una órbita en `r_pmax`. El sagital resuelve el equilibrio meridional pero no los modos azimutales; el ecuatorial resuelve `r,phi` pero no el espesor vertical; el 3D combina ambas direcciones.

| Caso | Pasos | Cambio de `M_total` | Mayor `v^2` registrado | `N_invalid` máximo | VTK |
|---|---:|---:|---:|---:|---|
| Ecuatorial | 18 137 | -1.329% | 0.9023 | 0 | 7/7 íntegros |
| Sagital | 18 137 | -2.839% | 0.9050 | 0 | 7/7 íntegros |
| 3D | 206 974 | -2.982% | 0.9199 | 0 | 7/7 íntegros |

Los tres procesos salieron con código 0, alcanzaron exactamente `t=180 M` y aprobaron la auditoría de valores finitos, floors, causalidad y diagnósticos. El cambio de `M_total` es compatible con un dominio radial abierto y no debe citarse como error de conservación de un volumen cerrado sin integrar los flujos de frontera y las reparaciones. Los modos `m=1,...,4` del 3D no perturbado permanecieron en el nivel de redondeo (`A_m` máximo inferior a `2.4e-16`), como corresponde a una corrida axisimétrica basal, pero esto no mide crecimiento PPI.

En el sagital y el 3D apareció un máximo transitorio de presión en la primera celda radial, `r≈1.216<r_+=2`, dentro del horizonte. Los VTK finales tienen presión exterior máxima cercana a `0.00202` alrededor de `r≈9.76`, y todos sus campos son finitos; este máximo interior debe excluirse de la interpretación del equilibrio externo y estudiarse por separado si se modifica la frontera o la excisión. La corrida de una órbita acredita estabilidad funcional de esta configuración, no ausencia de problemas en toda combinación de métrica, resolución, EOS, perturbación y tiempo.

### 8.2 Criterio de congelamiento para GRMHD

La implementación GRHD puede usarse como referencia de regresión para GRMHD si se conserva una instantánea identificable de fuente, `FFLAGS`, pruebas y `.par` de aceptación. En esta sesión el árbol de trabajo aún no tiene commit ni tag de cierre; por tanto, no existe todavía un identificador Git inmutable para esta versión. Antes de abrir una rama de implementación magnética conviene decidir cuál será esa instantánea, sin modificar retrospectivamente los resultados actuales. El límite `B->0` deberá reproducir las primitivas, conservadas, flujos, velocidades características y evolución GRHD con tolerancias explícitas.

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

El algoritmo reconstruye tiras por dirección, pero no aplica operadores temporales separados: las contribuciones espaciales se suman antes de cada etapa RK. Por ello no debe documentarse como separación de Strang.

## 10. Cómo extender el código

### 10.1 Añadir una condición inicial

1. Defina constantes y parámetros en `variables.f90`.
2. Añada el nombre en `parse_problem` y `problem_label`.
3. Declare y valide sus claves en `parameters.f90`.
4. Cree el preset en `initialization.f90`.
5. Implemente el dato inicial y las fronteras en `conditions.f90`.
6. Añada al menos una prueba unitaria y un `.par` pequeño reproducible.
7. Documente externamente el estado, las tolerancias y el resultado de verificación.

### 10.2 Añadir una métrica

La métrica debe proporcionar, en coordenadas físicas y si procede logarítmicas:

- `alpha`, `beta^i`, `gamma_ij` y `sqrt(gamma)`;
- inversa espaciotemporal cuando se solicite;
- derivadas espaciales de `g_mu_nu`;
- símbolos de Christoffel;
- derivadas de `ln(alpha)`.

Conecte las rutinas mediante `set_metric_type`. Pruebe determinante, inversa, simetrías, derivadas, horizonte, transformación de coordenadas y el límite a una métrica ya conocida.

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

Los directorios de producción, VTK, checkpoints y logs están excluidos de Git. Conserve datos grandes fuera del repositorio junto con un manifiesto trazable.

## 12. Límites de interpretación científica

- Fondo fijo/Cowling: no hay respuesta de la métrica al fluido.
- Sin autogravedad: se requiere una razón `M_torus/M_BH` físicamente pequeña para interpretar directamente el modelo de fluido de prueba.
- Sin MHD: no aparecen MRI, campos magnéticos ni transporte magnético.
- Una celda polar: no resuelve espesor vertical, *warping* ni modos 3D.
- Finn--Evans: proporciona un proxy de frecuencia/fase bajo una misma normalización, no amplitud absoluta observable.
- Barrido de espín histórico: cambian masa y densidad iniciales; no aísla el efecto causal del espín.
- Una sola semilla o resolución no establece incertidumbre estadística ni convergencia.
- La prueba 3D de una órbita se hizo en Schwarzschild, sin perturbación. No demuestra estabilidad de PPI ni convergencia polar/azimutal/radial, ni sustituye una evolución prolongada del toro Kerr.

## 13. Solución de problemas

### El programa imprime el uso y termina

Debe proporcionar un archivo:

```bash
./grhd2 archivo.par
```

### Se rechaza una clave válida para otro problema

Los parámetros de condición inicial se restringen al `problem` seleccionado. Elimine la clave ajena o cambie el problema.

### La malla radial logarítmica falla

Compruebe `logarithmic_r = true` y `r_min > 0`.

### EF rechaza el espín

La implementación EF representa Schwarzschild; use `spin = 0`. Para Kerr use `metric = kerr_schild` con `|spin| <= bh_mass`.

### Un toro 3D aparece vacío

Verifique que los centros polares realmente muestreen el cuerpo del toro. Un número pequeño de celdas en todo `[0,pi]` puede no colocar ninguna celda dentro de su espesor angular. Inspeccione el VTK de `final_time = 0` antes de lanzar una campaña larga.

### `GW_signal.dat` parece tener amplitud extraña

Compruebe `ny`, distancia registrada, sustracción de atmósfera, horizonte, mapeo y normalización. Con `ny=1` la amplitud absoluta no tiene normalización 3D resuelta.

### Una corrida termina sin `N_invalid`

Ese contador sólo existe cuando `ppi_diagnostics = true`. Su ausencia no demuestra que la corrida sea válida o inválida.
