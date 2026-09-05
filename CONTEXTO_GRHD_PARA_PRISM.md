# Contexto integral del proyecto GRHD para protocolo y futura tesis

**Fecha de corte:** 5 de septiembre de 2026  
**Fuente de verdad inspeccionada:** Zotz, `/home/fruelas/Codigo/GRHD`  
**Repositorio:** `git@github.com:RuelasF/GRHD.git`  
**Rama canónica:** `main`  
**Campañas principales:** `/home/fruelas/Codigo/GRHD_campaigns`

## 1. Propósito de este documento

Este documento concentra el contexto técnico, físico, numérico y científico
necesario para reorganizar el protocolo de investigación y delimitar una tesis
basada en el código GRHD. Distingue explícitamente:

1. lo que el código ya implementa;
2. lo que ha sido verificado numéricamente;
3. lo que sólo cuenta con evidencia exploratoria;
4. lo que todavía no puede afirmarse físicamente;
5. el trabajo mínimo necesario para una tesis defendible y reproducible.

No debe confundirse este proyecto con la línea independiente alojada en
`/home/fruelas/ollinsphere-bib/superaccretion`. Ambos proyectos pueden compartir
motivación relativista, pero utilizan códigos, modelos y objetivos diferentes.
Sólo deben integrarse en una misma tesis si se formula explícitamente un puente
científico entre ellos.

## 2. Resumen ejecutivo

GRHD es un simulador de hidrodinámica relativista general escrito en Fortran.
Emplea volúmenes finitos, reconstrucción de alto orden, solucionadores de
Riemann aproximados, integración SSP-RK3 y paralelismo OpenMP. Evoluciona un
fluido ideal relativista sobre fondos analíticos fijos de Minkowski,
Schwarzschild o Kerr.

El código es actualmente un **prototipo de investigación funcional con
validación parcial**. El núcleo compila, las pruebas disponibles pasan y las
campañas ecuatoriales de toros de Fishbone--Moncrief producen evoluciones
largas, diagnósticos de la inestabilidad de Papaloizou--Pringle (PPI), tasa de
acreción y un proxy cuadrupolar de ondas gravitacionales.

La evidencia más sólida hasta ahora es:

- convergencia correcta de una onda suave advectada para los reconstructores;
- pruebas unitarias del HLLC relativista en marcos cartesianos y métricos
  generales;
- pruebas de las transformaciones geométricas y del estrés Finn--Evans;
- 18 corridas 2.5D completas hasta `T=3000 M` sin celdas inválidas registradas;
- concordancia estrecha entre Eddington--Finkelstein y Kerr--Schild para
  `a=0` antes de la fase no lineal con WENO5 y MP5;
- crecimiento claro de modos no axisimétricos compatibles con la PPI;
- relación de frecuencia compatible con `f_GW aproximadamente 2 f_patron,m=1`.

Las limitaciones más importantes son:

- el espacio-tiempo no evoluciona y el disco no tiene autogravedad;
- la campaña principal tiene una sola celda polar y no es 3D;
- la masa y densidad iniciales del toro cambian mucho al variar el espín;
- sólo se utilizó una semilla de perturbación;
- la salida GW es un proxy de campo débil, no un strain físico calibrado;
- el piloto con cuatro celdas polares no resolvió el toro;
- no todos los doce problemas incluidos tienen validación automatizada;
- antes de la consolidación de septiembre, el código efectivo de agosto no
  estaba representado por el `main` de Git.

La tesis más defendible no debe centrarse todavía en una predicción absoluta de
ondas gravitacionales. El eje recomendado es el desarrollo, verificación y
aplicación del código al estudio controlado de la PPI en toros de
Fishbone--Moncrief. La señal cuadrupolar debe mantenerse como diagnóstico
secundario y como ruta de extensión.

## 3. Estado de control de versiones

### 3.1 Historia encontrada

El antiguo `main` apuntaba al commit `751af93`, fechado el 26 de mayo de 2026.
Ese commit era únicamente una base histórica y no representaba el código que
se ejecutó en agosto.

La rama `integracion/cluster-fm` avanzó hasta `d7a7fc4` el 13 de agosto e
incorporó correcciones Kerr--Schild, limpieza de artefactos, campañas de espín
y comparaciones controladas EF--KS. Después continuó desarrollo sin commit en
el directorio principal de Zotz.

Cronología reconstruida de la fuente efectiva:

- **13 de agosto:** HLLC relativista, pruebas de contacto y campaña de
  convergencia de la onda advectada.
- **23 de agosto:** perturbaciones reproducibles, reinicio configurable,
  diagnósticos PPI y diagnósticos globales.
- **25 de agosto:** protecciones de reconstrucción cerca de la atmósfera,
  discontinuidades fuertes y fronteras relativistas.
- **26 de agosto:** correcciones de la geometría Kerr--Schild, horizonte,
  transformación a coordenadas cartesianas, proxy Finn--Evans, manifiestos y
  pruebas GW.
- **27--28 de agosto:** análisis, auditoría y figuras de campañas; no se
  detectaron cambios posteriores del núcleo Fortran.

No se encontró una fuente GRHD posterior al 26 de agosto en Zotz. Los archivos
Fortran modificados a inicios de septiembre pertenecen a OllinSphere, no a
este repositorio.

### 3.2 Consolidación canónica

La versión promovida a `main` parte de `integracion/cluster-fm` para conservar
la historia de agosto y la eliminación intencional de `.history`, binarios y
artefactos. Sobre ella se integró la fuente efectiva de Zotz, junto con pruebas,
scripts de convergencia y documentación.

Los datos de producción no forman parte de Git. Las campañas conservan sus
propias copias de fuente, manifiestos y hashes SHA-256.

Antes de la consolidación se guardó en Zotz una copia recuperable en:

```text
/home/fruelas/Codigo/GRHD_safety/20260905_main_promotion
```

## 4. Alcance físico exacto

El código resuelve las ecuaciones de **Euler relativistas generales para un
fluido perfecto**. No debe describirse como un código de Navier--Stokes porque
no contiene viscosidad física.

Las variables primitivas son:

```text
rho, p, v^1, v^2, v^3
```

Las variables conservadas densitizadas son:

```text
D, tau, S_1, S_2, S_3
```

La densitización incluye `sqrt(gamma)`. El fluido usa una ecuación de estado
gamma-law/politrópica. Los casos de toros emplean normalmente `Gamma=4/3`.

### 4.1 Aproximación geométrica

La métrica es un fondo analítico fijo:

```text
g_mu_nu(t,x) = gbar_mu_nu(x)
```

Por tanto:

- no se resuelven las ecuaciones de Einstein;
- no existe retroacción del fluido sobre el agujero negro;
- no existe autogravedad del toro;
- no existe movimiento del agujero negro respecto al centro de masa;
- no se generan perturbaciones métricas dinámicas;
- la aproximación corresponde a fluido de prueba o Cowling.

Esta aproximación requiere, para interpretación física, que
`M_toro/M_BH << 1`. La normalización actual de las campañas no satisface todavía
una interpretación directa de esa relación; las masas diagnósticas son
principalmente cantidades internas del modelo 2.5D.

### 4.2 Física no implementada

- Magnetohidrodinámica relativista.
- Solucionador HLLD funcional.
- Campos magnéticos y divergencia de `B`.
- Radiación, enfriamiento o transporte de neutrinos.
- Viscosidad física.
- Autogravedad del fluido.
- Evolución BSSN, CCZ4 u otra formulación de Einstein.
- Extracción de `Psi_4`, RWZ, Teukolsky o CCE.

`is_mhd` y HLLD son marcadores de desarrollo futuro. HLLD se detiene de forma
explícita si se invoca.

## 5. Arquitectura del programa

El núcleo consta de aproximadamente 6,768 líneas distribuidas en diez fuentes
Fortran principales.

| Archivo | Responsabilidad principal |
|---|---|
| `variables.f90` | Parámetros globales, estado, cachés métricas, selectores numéricos y diagnósticos. |
| `metrics.f90` | Métricas 3+1, inversas, determinantes, derivadas, Christoffel y transformaciones geométricas. |
| `conditions.f90` | Condiciones iniciales, fronteras, Michel, Fishbone--Moncrief y perturbaciones. |
| `equations.f90` | Conversión primitiva/conservada, recuperación de primitivas, flujos físicos y fuentes. |
| `fluxes.f90` | HLLE, HLLC, velocidades características y marcador HLLD. |
| `reconstruction.f90` | Godunov, TVD, WENO3, WENO5, MP5 y fallbacks de robustez. |
| `evolution.f90` | Ensamblaje del RHS en las tres direcciones y paso CFL. |
| `output.f90` | VTK, checkpoints, manifiestos, PPI, acreción, proxy GW y convergencia. |
| `initialization.f90` | Selección de caso, parámetros, variables de entorno, malla y caché métrica. |
| `grhd2.f90` | Arranque, reinicio, SSP-RK3, perturbación retardada, diagnósticos y salida. |

El flujo de cálculo es:

```text
configuración
  -> construcción de malla y caché métrica
  -> condiciones iniciales
  -> primitivas a conservadas
  -> cálculo CFL
  -> tres etapas SSP-RK3
       -> reconstrucción por dirección
       -> problema de Riemann en cada cara
       -> divergencia de flujos + fuentes geométricas
       -> recuperación de primitivas
       -> condiciones de frontera
  -> diagnósticos, VTK y checkpoints
```

## 6. Geometrías y métricas

### 6.1 Minkowski

- Coordenadas cartesianas para problemas SRHD.
- Coordenadas cilíndricas para el jet axisimétrico.
- Minkowski esférica aparece como posibilidad futura, pero no está
  implementada.

### 6.2 Eddington--Finkelstein

- Fondo Schwarzschild con `a=0`.
- Coordenadas esféricas penetrantes del horizonte.
- Radio físico o logarítmico.

### 6.3 Kerr--Schild

- Fondo Kerr con `|a| <= M`.
- Coordenadas esféricas/oblatas penetrantes del horizonte.
- Radio físico o logarítmico.
- Incluye componentes espaciales no diagonales y transformaciones cartesianas
  compatibles con el espín.

El radio del horizonte usado es:

```text
r_+ = M + sqrt(M^2-a^2)
```

## 7. Métodos numéricos

### 7.1 Discretización

- Volúmenes finitos centrados en celda.
- Método de líneas.
- Contribuciones de flujo en `x`, `y` y `z` sumadas en el RHS.
- Integrador SSP-RK3 de tercer orden.
- Paso temporal adaptativo a partir del cono de luz coordenado.
- Paralelismo de memoria compartida OpenMP.

### 7.2 Reconstrucción

| ID | Método | Uso previsto |
|---:|---|---|
| 1 | Godunov | Referencia robusta de primer orden. |
| 2 | TVD Minmod/Superbee/MC | Segundo orden y choques. |
| 3 | WENO3 | Menor costo, pero puede sembrar PPI antes de tiempo. |
| 4 | MP5 | Alta resolución y buena preservación preperturbación. |
| 5 | WENO5-Z | Opción principal para regiones suaves y PPI. |

El sensor opcional degrada localmente la reconstrucción:

```text
alto orden -> WENO3 -> TVD-MC -> Godunov
```

Se activa según saltos de presión/densidad, proximidad a la atmósfera y
fronteras geométricas. En la campaña PPI controlada se desactivó para observar
el comportamiento puro de cada reconstructor; otras campañas usan las
protecciones.

### 7.3 Solucionadores de Riemann

- **HLLE:** baseline robusto de las campañas PPI.
- **HLLC:** implementado en una base ortonormal local de la cara. Considera la
  velocidad de la cara inducida por el shift y transforma el flujo al sistema
  coordenado. Ante estados degenerados o no admisibles vuelve a HLLE.
- **HLLD:** no implementado.

### 7.4 Recuperación de primitivas

La inversión conservadas--primitivas resuelve el factor de Lorentz mediante un
Newton protegido por intervalo. Si Newton no converge, utiliza Brent. Incluye:

- pisos locales densitizados;
- limitación causal del momento;
- rechazo de `NaN` o estados no físicos;
- reinicio local a la atmósfera como último recurso.

Esta estrategia aporta robustez, pero el número de celdas reparadas debe
registrarse y discutirse en cualquier análisis de conservación.

## 8. Casos iniciales incluidos

| ID | Caso | Métrica/geometría | Estado de evidencia |
|---:|---|---|---|
| 1 | Sod | Minkowski cartesiana 1D | Implementado; falta reporte automatizado actualizado. |
| 2 | Choque fuerte | Minkowski cartesiana 1D | Implementado; falta comparación de referencia consolidada. |
| 3 | Shu--Osher relativista | Minkowski cartesiana 1D | Implementado; validación cuantitativa pendiente. |
| 4 | Kelvin--Helmholtz | Minkowski cartesiana 2D | Implementado; estudio de convergencia pendiente. |
| 5 | Jet axisimétrico | Minkowski cilíndrica | Implementado; validación sistemática pendiente. |
| 6 | Onda advectada | Minkowski cartesiana 1D | Campaña de convergencia completa. |
| 7 | Acreción de Michel | EF esférica | Solución analítica implementada; convergencia final pendiente. |
| 8 | Acreción de polvo | EF esférica | Implementada; validación cuantitativa pendiente. |
| 9 | Explosión fuera del eje | EF esférica | Implementada; observable GW sólo exploratorio. |
| 10 | Toro FM ecuatorial | EF/KS, `r-phi` | Caso principal; 18 corridas 2.5D completas. |
| 11 | Toro FM sagital | EF, `r-theta` | Implementado; no es el baseline científico actual. |
| 12 | Toro FM 3D | EF/KS | Implementado formalmente; aún no validado ni resuelto adecuadamente. |

## 9. Configuración por defecto

El ejecutable sin variables de entorno selecciona:

```text
case_id              = 10
case                  = FishMoncEqu
geometry              = Spherical
metric                = Kerr-Schild
M_BH                  = 1
a                     = 0.9
grid                  = 400 x 1 x 200
r                     = [1.2, 40]
radial coordinate     = logarithmic
Gamma                 = 4/3
r_in                  = 6.25 M
r_max_density         = 9.25 M
K                     = 0.0015
final_time            = 1000 M
CFL                   = 0.4
save_interval         = 10 M
reconstruction        = WENO5
Riemann solver        = HLLE
shock sensor          = false
PPI diagnostics       = true
delayed perturbation  = false by default
Mdot extraction       = false by default
GW extraction         = false by default
```

Esta configuración no debe confundirse con la usada en las campañas de
`T=3000 M`.

## 10. Control reproducible por variables de entorno

El código permite configurar corridas sin editar las fuentes.

### Arranque y reinicio

```text
GRHD_DO_RESTART
GRHD_RESTART_FILE
```

### Caso, métrica y malla

```text
GRHD_CASE_ID
GRHD_METRIC_TYPE
GRHD_A_SPIN
GRHD_NX
GRHD_NY
GRHD_NZ
GRHD_R_MIN
GRHD_R_MAX
GRHD_USE_LOG_R
```

### Arquitectura numérica

```text
GRHD_RECONSTRUCTION
GRHD_TVD_LIMITER
GRHD_RIEMANN_SOLVER
GRHD_USE_SHOCK_SENSOR
```

### Tiempo y salida

```text
GRHD_FINAL_TIME
GRHD_CFL
GRHD_SAVE_INTERVAL
GRHD_OUTPUT_PREFIX
GRHD_OUTPUT_FOLDER
```

### Diagnósticos y perturbación

```text
GRHD_DO_GW_EXTRACTION
GRHD_DO_MDOT_EXTRACTION
GRHD_DO_PPI_DIAGNOSTICS
GRHD_DIAGNOSTIC_STRIDE
GRHD_APPLY_PERTURBATION
GRHD_PERTURBATION_TYPE
GRHD_PERTURBATION_SEED
GRHD_PERTURBATION_TIME
GRHD_PERTURBATION_AMPLITUDE
GRHD_PERTURBATION_MODE
```

Tipos de perturbación:

```text
0 = ninguna
1 = ruido de presión
2 = ruido de densidad
3 = modo azimutal de densidad
```

## 11. Salidas y diagnósticos

Cada corrida puede producir:

- VTK binario para VisIt/ParaView;
- checkpoints binarios para reinicio;
- `run_manifest.txt` con la configuración;
- `ppi_modes.dat`;
- `global_diagnostics.dat`;
- `perturbation_events.dat`;
- `m_dot.dat`;
- `GW_signal.dat`;
- normas de convergencia para la onda suave o Michel.

### 11.1 Modos PPI

Se calcula la masa azimutal y sus coeficientes complejos para `m=1,...,4`:

```text
A_m = |C_m| / C_0
```

También se almacena la fase del patrón. El modo `m=1` es el observable primario
para la PPI y los modos superiores permiten seguir armónicos y fase no lineal.

### 11.2 Diagnósticos globales

Se almacenan:

- masa total y masa por encima del corte de atmósfera;
- momento angular azimutal;
- mínimos y máximos de densidad y presión;
- máximo de `v^2`;
- número de celdas de atmósfera;
- número de celdas inválidas.

### 11.3 Tasa de acreción

`Mdot` se integra en la primera celda exterior al horizonte usando la velocidad
coordenada efectiva `alpha v^r - beta^r` y la densidad conservada.

## 12. Verificación numérica existente

### 12.1 Pruebas HLLC

La suite comprueba:

- estados uniformes en las tres direcciones cartesianas;
- estados uniformes con métrica espacial general;
- contacto estacionario;
- contactos móviles y tangenciales;
- contacto en una métrica general;
- comportamiento a escala de atmósfera;
- invariancia por cambio de escala;
- salto fuerte de presión relativista.

Todas las pruebas pasan en Zotz.

### 12.2 Onda sinusoidal advectada

Campaña:

```text
AdvectedWave_HLLC_Convergence_20260813
N = 32, 64, 128, 256, 512, 1024
CFL = 0.1
T = 0.7
SSP-RK3 + HLLC
```

Resultado resumido:

| Método | Comportamiento observado |
|---|---|
| Godunov | Orden aproximadamente 1. |
| TVD-MC | Orden aproximadamente 2 en L1; menor en Linf cerca de extremos. |
| WENO3 | Orden aproximadamente 3. |
| MP5 | Cerca de orden 5 en mallas útiles; saturación de redondeo al refinar. |
| WENO5 | Cerca de orden 5 en mallas útiles; saturación de redondeo al refinar. |

En `N=1024`, MP5/WENO5 alcanzan errores del orden de `2e-13`; la caída aparente
del orden en ese extremo está dominada por precisión de máquina.

### 12.3 Pruebas GW auxiliares

- Radio del horizonte y transformaciones geométricas EF/KS.
- Transformación de velocidades eulerianas a tasas coordenadas.
- Coordenadas cartesianas Kerr--Schild.
- Fórmula local de estrés Finn--Evans.

Estas pruebas validan la implementación algebraica del proxy, no convierten el
proxy en una extracción relativista gauge-invariant.

## 13. Campaña PPI principal

Ruta:

```text
/home/fruelas/Codigo/GRHD_campaigns/
PPI_GWFIX_EF_KS_RECON_2D3D_20260826_173840
```

La campaña conserva:

- fuente congelada;
- parche aplicado al caso 3D;
- hashes SHA-256;
- ejecutable y hash;
- manifiesto de casos;
- logs y tiempos;
- checkpoints;
- VTK;
- diagnósticos auxiliares;
- scripts de análisis;
- reportes y figuras.

### 13.1 Matriz completada

Se completaron 18 corridas ecuatoriales:

```text
metric/spin:
  EF a=0
  KS a=0, 0.2, 0.4, 0.6, 0.9

reconstruction:
  WENO3, WENO5, MP5

common:
  HLLE
  grid = 400 x 1 x 200
  r = [1.2, 40]
  linear radial coordinate
  T = 3000 M
  CFL = 0.4
  perturbation = 1% pressure noise
  perturbation time = 1000 M
  seed = 3435
  save interval = 10 M
```

Cada corrida terminó con código cero, produjo 301 VTK y 30 checkpoints. No se
registraron celdas inválidas.

### 13.2 Diferencia respecto al default

La campaña usó radio **lineal**, mientras el default actual usa radio
**logarítmico**. También extendió `T` de 1000 a 3000 y activó perturbación,
`Mdot`, PPI y proxy GW. Esta diferencia debe aparecer en toda tabla metodológica
y pie de figura.

## 14. Hallazgos científicos preliminares

### 14.1 Inicio de la PPI

- WENO3 desarrolla modos no axisimétricos antes de `T=1000`.
- Por tanto, la perturbación explícita no controla completamente el inicio de
  esas corridas.
- WENO5 y MP5 mantienen `A1` cerca del ruido numérico hasta la perturbación.
- WENO5 y MP5 son los candidatos principales para mediciones controladas.
- WENO3 es útil para estudiar sensibilidad numérica, no como único baseline.

### 14.2 Tasas de crecimiento

En los ajustes refinados de WENO5/MP5, los tiempos e-fold del modo `m=1` están
aproximadamente entre `30 M` y `78 M`. Cada valor debe publicarse junto con:

- ventana temporal usada;
- coeficiente `R^2`;
- resolución;
- semilla;
- reconstructor;
- métrica y espín;
- criterio de selección del intervalo exponencial.

No debe reportarse una tasa automática sin inspección visual y análisis de
sensibilidad de la ventana.

### 14.3 Control EF contra KS para `a=0`

EF y KS describen Schwarzschild con foliaciones distintas. Antes de la fase no
lineal:

- WENO5 y MP5 muestran concordancia muy estrecha en masa y acreción;
- diferencias cercanas al ruido de máquina en los modos no deben compararse de
  forma relativa;
- después del crecimiento no lineal, diferencias de redondeo pueden amplificarse
  caóticamente sin implicar que las formulaciones físicas sean incompatibles.

Este control es necesario pero todavía debe complementarse con observables
invariantes o adecuadamente alineados en tiempo retardado.

### 14.4 Frecuencia de patrón y proxy GW

La campaña encontró:

```text
media de f_GW/(2 f_patron,m=1) = 1.01842
intervalo observado            = 0.88215 a 1.14326
```

La coincidencia apoya una relación dinámica entre la estructura `m=1` y el
proxy cuadrupolar. No calibra la amplitud del strain.

### 14.5 Región cercana al horizonte

El mayor `v^2` registrado fue aproximadamente `0.99998222` para `a=0.9`. La
inspección espacial lo localizó dentro del horizonte de Kerr
`r_+=1.43589 M`, no en el cuerpo exterior del toro. Debe seguir vigilándose el
margen de la recuperación de primitivas, pero este valor no invalida por sí
solo la dinámica exterior.

## 15. Confusor principal del barrido de espín

La secuencia actual mantiene `r_in`, `r_max_density` y `K` nominales mientras
el equilibrio cambia con el espín. Como resultado:

```text
M0(a=0.9) / M0(a=0)                 = 8.64969
rho_max,0(a=0.9) / rho_max,0(a=0)   = 6.92319
```

Por ello, cualquier tendencia cruda de:

- tasa de crecimiento;
- amplitud de saturación;
- masa perdida;
- tasa de acreción;
- amplitud Finn--Evans;

mezcla el efecto del espín con cambios en masa, densidad, entalpía y estructura
del toro. La campaña actual es exploratoria y no aísla causalmente el espín.

Para la campaña de tesis debe elegirse y documentarse una familia normalizada,
por ejemplo manteniendo fija una de estas cantidades:

- `rho_max`;
- masa en reposo del toro;
- razón `M_toro/M_BH`;
- geometría/equipotencial y una escala de densidad físicamente declarada.

La elección debe ser físicamente compatible con Cowling y comprobar que las
conclusiones no dependan únicamente de esa normalización.

## 16. Significado exacto de 2.5D

La malla principal `400 x 1 x 200` tiene una sola celda en `theta`, centrada en
el ecuador, y evoluciona `r` y `phi`. Esto permite estudiar modos azimutales,
pero no resuelve la sección vertical del toro.

Consecuencias:

- la morfología es ecuatorial, no tridimensional;
- la masa usa implícitamente el ancho polar completo de una sola celda;
- `M_total` y `M_disk` no son integrales 3D físicas;
- la amplitud Finn--Evans depende de una cuadratura polar no resuelta;
- no se capturan modos verticales, warping ni acoplamientos 3D;
- no se puede demostrar convergencia polar.

En la tesis debe usarse “modelo ecuatorial 2.5D” y evitar llamar a estas
corridas “simulaciones 3D”.

## 17. Estado real del piloto 3D

Se planeó una etapa con `N_theta=4`, pero sólo comenzó el caso
`n4_ef_a0_weno3`.

Hechos:

- malla `400 x 4 x 200`;
- todas las 320,000 celdas fueron clasificadas como atmósfera en `T=0`;
- `M_disk=0` desde el inicio;
- el toro no quedó muestreado por la malla polar;
- la corrida avanzó hasta aproximadamente `T=2734 M`;
- terminó con código 130 porque la campaña recibió una señal externa;
- no fue un colapso numérico del solver;
- la cola 3D se detuvo y los casos restantes no corrieron.

La interpretación correcta es: **no existe todavía una simulación 3D válida de
la PPI con este código**. El primer requisito es escoger `N_theta` a partir del
espesor angular del toro y demostrar que varias celdas cubren la región de
densidad significativa.

## 18. Ondas gravitacionales: implementación actual

`GW_signal.dat` se obtiene mediante una fórmula de estrés Finn--Evans de campo
débil aplicada al fluido GRHD. El cálculo:

- transforma las velocidades eulerianas en tasas coordenadas;
- usa coordenadas cartesianas compatibles con Kerr--Schild;
- resta la atmósfera numérica;
- excluye celdas dentro del horizonte;
- usa el Jacobiano cartesiano plano de las coordenadas KS;
- coloca convencionalmente al observador sobre el eje `+z`;
- usa distancia `R=1000` en unidades geométricas;
- produce `h_plus` y `h_cross` como proxies.

El propio manifiesto lo identifica como:

```text
Finn-Evans weak-field stress proxy
```

### 18.1 Lo que sí puede usarse

- frecuencia dominante;
- fase relativa;
- tiempo de aparición de la señal;
- correlación con modos PPI;
- comparación metodológica bajo una misma normalización.

### 18.2 Lo que no puede afirmarse aún

- strain observable absoluto;
- energía radiada físicamente calibrada;
- independencia de gauge;
- extracción en infinito nulo;
- amplitud tridimensional;
- retroacción de la radiación o autogravedad.

La ruta de extensión está documentada en
`extraccion_gw_gauge_invariant.md`: primero normalización física y resolución
polar, después perturbaciones de Schwarzschild/RWZ con fuente y, sólo tras esa
validación, una posible extensión Teukolsky para Kerr.

## 19. Qué está verificado, validado o sólo explorado

| Afirmación | Nivel actual |
|---|---|
| El código compila en Zotz | Verificado. |
| Las pruebas HLLC pasan | Verificado. |
| Las transformaciones y estrés Finn--Evans pasan sus pruebas | Verificado algebraicamente. |
| La onda suave recupera los órdenes esperados | Verificado. |
| Las 18 corridas 2.5D llegan a `T=3000` | Verificado. |
| No hay celdas inválidas registradas en esas corridas | Verificado por el diagnóstico implementado. |
| WENO3 siembra modos antes que WENO5/MP5 | Evidencia numérica fuerte en una resolución/semilla. |
| Existe crecimiento compatible con PPI | Evidencia exploratoria fuerte. |
| La frecuencia GW proxy sigue aproximadamente `2 f_m1` | Evidencia exploratoria. |
| El espín cambia causalmente la tasa de crecimiento | No demostrado por falta de normalización controlada. |
| La amplitud GW es física | No demostrado; actualmente falso como interpretación estricta. |
| El código está validado en 3D | No. |
| El toro es autogravitante | No. |
| El código es GRMHD | No. |

## 20. Formulación recomendada de la tesis

### 20.1 Título provisional

> Desarrollo, verificación y aplicación de un código de hidrodinámica
> relativista general de volúmenes finitos al estudio de la inestabilidad de
> Papaloizou--Pringle en toros de Fishbone--Moncrief alrededor de agujeros
> negros de Schwarzschild y Kerr.

Una versión más centrada en metodología:

> Verificación de métodos de alto orden para GRHD y su aplicación a modos no
> axisimétricos de toros relativistas en fondos de agujeros negros.

### 20.2 Planteamiento del problema

La PPI puede producir estructuras no axisimétricas, redistribución de momento
angular, acreción variable y emisión cuadrupolar. Su medición numérica depende
de la resolución, reconstrucción, solver, perturbación, geometría de
coordenadas y normalización del toro. Antes de atribuir tendencias físicas al
espín es necesario separar esos efectos numéricos y de condiciones iniciales.

### 20.3 Pregunta central

> ¿Cómo dependen del espín del agujero negro la tasa de crecimiento, frecuencia
> de patrón, saturación y acreción asociadas a la PPI cuando se controlan la
> normalización del toro, la resolución, la semilla y la disipación numérica?

### 20.4 Hipótesis de trabajo

> Para una familia de toros físicamente normalizada y suficientemente resuelta,
> las tasas de crecimiento y frecuencias de la PPI convergerán con la resolución
> y serán consistentes entre WENO5 y MP5. Las diferencias EF--KS para `a=0`
> disminuirán en observables comparables, mientras que las variaciones
> persistentes al cambiar `a` podrán interpretarse como efecto del fondo Kerr y
> no como artefactos de normalización.

La hipótesis debe refinarse después de definir exactamente qué propiedad del
toro se mantiene fija.

### 20.5 Objetivo general

Desarrollar, verificar y aplicar una plataforma GRHD reproducible para
caracterizar la PPI en toros relativistas sobre fondos Schwarzschild/Kerr,
cuantificando incertidumbres numéricas y separando los efectos del espín de los
de la discretización y las condiciones iniciales.

### 20.6 Objetivos específicos

1. Documentar la formulación GRHD, convenciones 3+1, variables y unidades.
2. Consolidar una versión reproducible del código y una suite de regresión.
3. Verificar órdenes de convergencia y soluciones analíticas SRHD/GRHD.
4. Validar la preservación de equilibrio de toros Fishbone--Moncrief.
5. Construir una familia de toros normalizada al variar el espín.
6. Medir crecimiento, fase, frecuencia y saturación de `m=1,...,4`.
7. Cuantificar sensibilidad a resolución, reconstructor, solver y semilla.
8. Verificar equivalencia EF--KS en Schwarzschild.
9. Extender el estudio a resolución polar real y posteriormente a 3D.
10. Relacionar la dinámica PPI con `Mdot` y el proxy cuadrupolar, declarando sus
    límites de interpretación.

## 21. Metodología recomendada

### Fase A: consolidación y regresión

- Mantener `main` como fuente canónica.
- Etiquetar versiones usadas en campañas.
- Ejecutar compilación optimizada y compilación estricta de depuración.
- Automatizar HLLC, GW, onda advectada y smoke tests.
- No guardar VTK, checkpoints, logs o ejecutables en Git.
- Guardar manifiesto, commit, compilador, flags y hash del ejecutable por
  campaña.

### Fase B: verificación hidrodinámica

- Sod, choque fuerte y Shu--Osher frente a referencias reproducibles.
- Convergencia de la onda suave.
- Acreción de Michel frente a la solución analítica.
- Conservación y costo de pisos/reparaciones.
- Comparación HLLE--HLLC en problemas donde el contacto sea relevante.

### Fase C: preservación del equilibrio FM

- Evoluciones sin perturbación.
- Tres resoluciones como mínimo.
- WENO5 y MP5 como métodos principales.
- Medir deriva de masa, densidad máxima, posición del máximo, momento angular
  y amplitudes `A_m` antes de perturbar.
- Comparar radio lineal y logarítmico.
- Comparar EF y KS para `a=0`.

### Fase D: familia controlada de espín

- Definir la normalización física del toro.
- Recalcular la constante politrópica o escala de densidad cuando sea necesario.
- Mantener explícitamente fija la cantidad elegida.
- Usar `a=0,0.2,0.4,0.6,0.9` u otra malla justificada.
- Usar varias semillas, idealmente al menos tres y preferentemente cinco.
- Separar fase lineal, saturación y fase no lineal.
- Incluir incertidumbre por resolución, semilla y método.

### Fase E: resolución polar y 3D

- Medir el espesor angular del toro inicial.
- Elegir `N_theta` por celdas a través de la escala vertical, no por costo
  arbitrario.
- Hacer primero pilotos cortos y verificar `M_disk>0` en `T=0`.
- Realizar convergencia polar.
- Sólo después lanzar campañas largas 3D.

### Fase F: observables

- Correlacionar `A_m`, fase, frecuencia, `Mdot` y proxy GW.
- Para amplitud física, fijar `M_toro/M_BH` y distancia.
- Mantener Finn--Evans como proxy mientras no exista extracción perturbativa.
- Considerar RWZ con fuente en Schwarzschild como extensión de tesis si el
  calendario lo permite.

## 22. Variables de respuesta y análisis estadístico

Variables principales:

- `A_m(t)` para `m=1,...,4`;
- fase y velocidad angular del patrón;
- tasa de crecimiento `sigma_m` y tiempo e-fold;
- tiempo de cruce de umbrales;
- amplitud de saturación;
- masa total y masa del disco;
- momento angular;
- `Mdot(t)`;
- máximos de densidad, presión y `v^2`;
- frecuencia dominante del proxy GW;
- relación `f_GW/(2 f_m1)`;
- número de celdas reparadas/inválidas;
- costo computacional.

Cada ajuste debe reportar:

- ventana temporal;
- método de selección de ventana;
- `R^2` o estadístico equivalente;
- incertidumbre del ajuste;
- dispersión entre semillas;
- diferencia entre resoluciones;
- diferencia entre WENO5 y MP5.

No deben usarse errores relativos de amplitudes cercanas a `1e-16`, porque
carecen de significado físico.

## 23. Criterios mínimos de aceptación

Una corrida científica debe:

1. terminar con código de salida cero;
2. conservar manifiesto y commit de fuente;
3. no contener `NaN`, infinitos o celdas inválidas fuera de regiones excluidas;
4. resolver espacialmente el toro inicial;
5. registrar el momento exacto de perturbación;
6. tener una fase preperturbación suficientemente estable;
7. demostrar que el crecimiento medido no desaparece al refinar;
8. incluir balance de masa compatible con flujo por fronteras y horizonte;
9. separar datos físicos de atmósfera;
10. mantener sin cambios todos los parámetros no incluidos en el barrido.

Para una tendencia con espín:

- la familia inicial debe estar normalizada;
- debe repetirse con más de una semilla;
- debe ser compatible entre al menos dos reconstructores de alto orden;
- debe incluir estimación de error por resolución.

Para un resultado GW:

- frecuencia/fase proxy: aceptable con las advertencias actuales;
- amplitud física: requiere resolución polar, masa física y formalismo de
  extracción justificado.

## 24. Riesgos del proyecto

| Riesgo | Consecuencia | Mitigación |
|---|---|---|
| Fuente no versionada | Resultados irreproducibles | `main`, tags, snapshots y hashes. |
| Normalización cambia con `a` | Falsa tendencia de espín | Familia FM controlada. |
| WENO3 rompe simetría temprano | Semilla no controlada | WENO5/MP5 y estudio de resolución. |
| Una sola semilla | Sin incertidumbre estadística | Ensamble de semillas. |
| `N_theta` insuficiente | Ausencia o deformación del toro | Criterio de celdas por escala vertical. |
| Atmósfera y floors | Pérdida artificial de masa | Presupuesto de masa y conteo de reparaciones. |
| Frontera radial cercana | Reflexiones/contaminación | Prueba de dominio y frontera. |
| Proxy GW sobreinterpretado | Conclusión física inválida | Etiquetado explícito y normalización. |
| Cowling con disco pesado | Modelo inconsistente | Elegir `M_toro/M_BH << 1`. |
| Campañas muy grandes | Costo y organización | Pilotos, manifiestos y criterios de parada. |

## 25. Contribuciones potenciales de la tesis

1. Implementación modular GRHD con métricas penetrantes del horizonte.
2. HLLC relativista mediante marcos ortonormales locales con fallback robusto.
3. Comparación sistemática de WENO3, WENO5 y MP5 para PPI.
4. Control EF--KS para Schwarzschild.
5. Metodología reproducible de perturbaciones y modos azimutales.
6. Separación cuantitativa de errores numéricos y tendencias físicas.
7. Catálogo normalizado de crecimiento y frecuencia PPI frente a espín.
8. Ruta documentada desde proxy cuadrupolar hacia extracción perturbativa.

Las contribuciones 6 y 7 requieren las campañas controladas pendientes; no
deben presentarse todavía como resultados concluidos.

## 26. Estructura sugerida de la tesis

1. **Introducción y motivación astrofísica.**
2. **Hidrodinámica relativista general en formulación 3+1.**
3. **Métodos de volúmenes finitos y arquitectura del código.**
4. **Verificación numérica y soluciones de referencia.**
5. **Toros de Fishbone--Moncrief y PPI.**
6. **Campañas controladas: coordenadas, espín y métodos numéricos.**
7. **Acreción, modos y proxy cuadrupolar.**
8. **Extensión polar/3D y limitaciones de Cowling.**
9. **Conclusiones y trabajo futuro.**

RWZ/Teukolsky debe ser capítulo futuro o adicional, salvo que se implemente y
valide dentro del calendario real de la tesis.

## 27. Prioridades inmediatas

1. Conservar la versión consolidada de `main` y etiquetar el baseline.
2. Actualizar el README y documentar compilación/pruebas.
3. Añadir una compilación estricta reproducible.
4. Consolidar resultados de Sod, choque fuerte, Shu--Osher y Michel.
5. Implementar una campaña de preservación FM sin perturbación.
6. Diseñar la normalización fija al variar el espín.
7. Ejecutar tres resoluciones y varias semillas con WENO5/MP5.
8. Diseñar un piloto polar que contenga realmente el toro.
9. Separar formalmente datos, figuras y fuente.
10. Decidir si la tesis se centra en PPI o si sólo usa PPI como validación para
    otra línea; no mezclar automáticamente con superacreción/OllinSphere.

## 28. Rutas de evidencia en Zotz

Código principal:

```text
/home/fruelas/Codigo/GRHD
```

Campaña PPI/GW principal:

```text
/home/fruelas/Codigo/GRHD_campaigns/
PPI_GWFIX_EF_KS_RECON_2D3D_20260826_173840
```

Reporte refinado:

```text
analysis_2p5d_20260827_v2/analysis_report_v2.md
```

Resumen de strain y acreción:

```text
strain_mdot_2p5d_20260827/README.md
```

Estado de casos:

```text
logs/case_status.tsv
```

Parámetros de la matriz:

```text
provenance/cases.tsv
```

Fuentes y hashes usados:

```text
source/
provenance/
```

Respaldo previo a la promoción de `main`:

```text
/home/fruelas/Codigo/GRHD_safety/20260905_main_promotion
```

## 29. Instrucción lista para Prism

Usar este documento como fuente primaria de contexto para reorganizar el
protocolo. Al proponer texto académico, Prism debe:

1. distinguir verificación numérica, validación física y exploración;
2. llamar al modelo principal “GRHD sobre fondo fijo/Cowling”;
3. llamar a las corridas `400 x 1 x 200` “ecuatorial 2.5D”;
4. llamar a `GW_signal.dat` “proxy Finn--Evans de campo débil”;
5. no presentar amplitudes GW como strain observable;
6. no atribuir tendencias al espín sin corregir la normalización inicial;
7. tratar WENO5 y MP5 como métodos principales y WENO3 como control de
   sensibilidad;
8. incluir resolución y semillas en la metodología;
9. presentar la simulación 3D como trabajo pendiente;
10. mantener separado el proyecto GRHD de OllinSphere/superacreción salvo
    decisión expresa del investigador.

La salida solicitada a Prism debería contener:

- título provisional y variantes;
- resumen y planteamiento del problema;
- preguntas de investigación;
- hipótesis principal y secundarias;
- objetivo general y objetivos específicos;
- marco teórico por secciones;
- metodología con matriz de experimentos;
- criterios de aceptación y análisis de incertidumbre;
- resultados preliminares claramente etiquetados;
- trabajo pendiente;
- cronograma realista;
- índice propuesto de protocolo y tesis;
- lista de afirmaciones que requieren bibliografía o validación adicional.

## 30. Síntesis final para decisiones

El proyecto ya posee suficiente código y evidencia para sustentar una tesis de
desarrollo y verificación GRHD aplicada a la PPI. Todavía no posee una base
suficiente para afirmar de manera definitiva una ley física de dependencia con
el espín ni para predecir amplitudes observables de ondas gravitacionales.

La ruta de menor riesgo científico es:

```text
consolidar código
  -> completar verificación
  -> preservar equilibrio FM
  -> normalizar familia de toros
  -> convergencia + ensamble de semillas
  -> medición PPI frente a espín
  -> resolución polar/3D
  -> observables GW más rigurosos
```

El resultado de agosto es una base exploratoria fuerte. La contribución de
tesis surgirá al convertirla en una campaña controlada, convergente,
reproducible y físicamente interpretable.
