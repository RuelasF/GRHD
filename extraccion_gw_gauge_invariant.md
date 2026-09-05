# Extracción gauge-invariant de ondas gravitacionales en GRHD

## Comparación directa entre Eddington–Finkelstein y Kerr–Schild

La campaña 2.5D de agosto de 2026 ya completó la comparación controlada entre
Eddington–Finkelstein (EF) y Kerr–Schild (KS) para \(a=0\), con WENO3, WENO5 y
MP5, todos con HLLE:

| Métrica/coordenadas | WENO3 | WENO5 |
|---|---:|---:|
| EF, \(a=0\) | completa | completa |
| KS, \(a=0\) | completa | completa |

También se completaron casos KS con \(a=0.2,0.4,0.6,0.9\). EF y KS con
\(a=0\) describen el mismo espacio-tiempo de Schwarzschild, aunque emplean
foliaciones distintas. Los diagnósticos de masa y acreción concuerdan
estrechamente antes de la fase no lineal para WENO5 y MP5. Esta comparación es
un control importante, pero no sustituye una extracción gauge-invariant: los
observables invariantes futuros deberán compararse en un mismo tiempo retardado
y con las transformaciones tensoriales correspondientes.

## Qué significa realmente *gauge-invariant*

Escribimos la métrica como

\[
g_{\mu\nu}=\bar g_{\mu\nu}+h_{\mu\nu},
\]

donde \(\bar g_{\mu\nu}\) es el fondo —Schwarzschild o Kerr— y \(h_{\mu\nu}\) es una perturbación pequeña.

Bajo una transformación infinitesimal de coordenadas

\[
x^\mu\rightarrow x^\mu+\xi^\mu,
\]

la perturbación cambia como

\[
h_{\mu\nu}\rightarrow h_{\mu\nu}
-\bar\nabla_\mu\xi_\nu
-\bar\nabla_\nu\xi_\mu.
\]

Por tanto, una componente como \(h_{rr}\), \(h_{t\phi}\), o incluso una combinación construida directamente en una esfera de radio finito, puede contener:

- radiación física;
- una perturbación puramente coordinada;
- campos cercanos o de Coulomb;
- errores por la elección de la tetrada;
- efectos del radio de extracción finito.

En sentido estricto, la radiación gravitacional se define sin ambigüedad en el infinito nulo futuro \(\mathscr I^+\), donde se puede construir el *strain* de Bondi, el *news* o una \(\Psi_4\) con una tetrada asintótica apropiada. La extracción característica lleva los datos precisamente hasta \(\mathscr I^+\) y evita las ambigüedades de radio finito. [Reisswig et al. (2009)](https://arxiv.org/abs/0912.1285) describen esta construcción mediante *Cauchy-characteristic extraction* (CCE).

## Qué ocurre actualmente en GRHD

En la aproximación de Cowling se impone

\[
g_{\mu\nu}(t,\mathbf{x})=\bar g_{\mu\nu}(\mathbf{x}),
\]

es decir, la métrica nunca responde al fluido. En consecuencia:

- El código no evoluciona \(h_{\mu\nu}\).
- Si calculáramos \(\Psi_4\) solamente a partir de la métrica fija, obtendríamos la del fondo, no la radiación producida por el toro.
- El `gw_strain` actual es una estimación cuadrupolar construida directamente a partir del fluido.
- Su frecuencia puede servir como diagnóstico, pero su amplitud y forma no son una extracción gauge-invariant.

Esto no significa que necesitemos inmediatamente una evolución completa de las ecuaciones de Einstein. Hay una opción intermedia adecuada para el proyecto.

## Perturbaciones con fuente sobre el fondo fijo

Se puede mantener la evolución GRHD en Cowling y, simultáneamente, resolver

\[
G_{\mu\nu}^{(1)}[h]=8\pi T_{\mu\nu},
\]

donde:

- \(G_{\mu\nu}^{(1)}\) es el tensor de Einstein linealizado alrededor del agujero negro;
- \(T_{\mu\nu}\) se calcula con el fluido GRHD;
- \(h_{\mu\nu}\) se evoluciona como una perturbación que no retroactúa sobre el fluido.

Para un fluido perfecto,

\[
T^{\mu\nu}
=
\rho_0 h_{\mathrm{fluid}}u^\mu u^\nu
+p\,g^{\mu\nu},
\]

con

\[
h_{\mathrm{fluid}}=1+\epsilon+\frac{p}{\rho_0}.
\]

Esta aproximación es consistente si

\[
\frac{M_{\mathrm{torus}}}{M_{\mathrm{BH}}}\ll1.
\]

La dinámica del fluido puede ser fuertemente no lineal —la PPI puede saturar y producir choques— mientras la perturbación gravitacional siga siendo lineal, siempre que la masa total del toro sea pequeña comparada con la del agujero negro.

La onda se calcula entonces de manera relativista a partir del \(T_{\mu\nu}\) completo, no mediante un cuadrupolo newtoniano.

### Una ventaja particular del problema

El fondo Schwarzschild/Kerr es vacío:

\[
\bar T_{\mu\nu}=0.
\]

Por eso, el \(T_{\mu\nu}\) de primer orden es gauge-invariant bajo transformaciones infinitesimales del fondo. El problema principal pasa a ser construir variables maestras gauge-invariant a partir de las perturbaciones métricas y de sus fuentes.

## Schwarzschild: formalismo Regge–Wheeler–Zerilli–Moncrief

Para EF y KS con \(a=0\), la recomendación es implementar primero el formalismo Regge–Wheeler–Zerilli (RWZ).

Una perturbación de Schwarzschild se descompone en armónicos esféricos:

\[
h_{\mu\nu}
=
\sum_{\ell,m}
\left(
h_{\mu\nu}^{\ell m,\mathrm{even}}
+
h_{\mu\nu}^{\ell m,\mathrm{odd}}
\right).
\]

Los sectores son:

- Paridad par o polar: perturbaciones tipo Zerilli–Moncrief.
- Paridad impar o axial: perturbaciones tipo Regge–Wheeler/Cunningham–Price–Moncrief.

Para cada \((\ell,m)\), las numerosas componentes de \(h_{\mu\nu}\) se combinan en solamente dos escalares gauge-invariant:

\[
\Psi_{\ell m}^{\mathrm{even}},
\qquad
\Psi_{\ell m}^{\mathrm{odd}}.
\]

Estos satisfacen ecuaciones de onda independientes:

\[
\left[
-\partial_t^2
+\partial_{r_*}^2
-V_\ell^{\mathrm{even/odd}}(r)
\right]
\Psi_{\ell m}^{\mathrm{even/odd}}
=
S_{\ell m}^{\mathrm{even/odd}},
\]

donde

\[
r_*=r+2M\ln\left(\frac{r}{2M}-1\right)
\]

es la coordenada tortuga.

Definiendo

\[
f=1-\frac{2M}{r},
\qquad
\lambda=\frac{(\ell-1)(\ell+2)}{2},
\]

el potencial impar es

\[
V_\ell^{\mathrm{odd}}
=
f\left[
\frac{\ell(\ell+1)}{r^2}
-\frac{6M}{r^3}
\right],
\]

mientras que el potencial par es

\[
V_\ell^{\mathrm{even}}
=
\frac{2f}{
r^3(\lambda r+3M)^2
}
\left[
\lambda^2(\lambda+1)r^3
+3\lambda^2Mr^2
+9\lambda M^2r
+9M^3
\right].
\]

Las fuentes \(S_{\ell m}\) se obtienen proyectando \(T_{\mu\nu}\) sobre armónicos escalares, vectoriales y tensoriales. La formulación de Martel y Poisson proporciona las fuentes completas y covariantes para materia, así como la reconstrucción del *strain* y los flujos radiados. Véase [Martel y Poisson (2005)](https://arxiv.org/abs/gr-qc/0502028).

### Compatibilidad con EF y KS

No es necesario transformar primero toda la evolución a coordenadas Schwarzschild, cuya foliación es singular en el horizonte. Las ecuaciones pueden formularse covariantemente sobre una métrica esféricamente simétrica

\[
ds^2=g_{ab}(x^c)dx^a dx^b+r^2d\Omega^2,
\qquad x^a=(t,r),
\]

y usarse directamente en foliaciones penetrantes del horizonte.

Sarbach y Tiglio derivaron una formulación gauge-invariant para foliaciones arbitrarias de Schwarzschild, incluyendo coordenadas penetrantes como EF y KS. Véase [Sarbach y Tiglio (2001)](https://arxiv.org/abs/gr-qc/0104061).

Esto sería excelente para validar el código:

\[
\Psi_{\ell m}^{\mathrm{EF}}(u)
\approx
\Psi_{\ell m}^{\mathrm{KS},a=0}(u)
\]

después de expresar ambas señales en el mismo tiempo retardado \(u\). Las variables primitivas o componentes de velocidad coordinadas no tienen por qué coincidir componente a componente, pero los escalares maestros y la radiación asintótica sí.

## De las variables maestras al *strain*

En la zona de ondas,

\[
h_+-ih_\times
=
\frac{1}{R}
\sum_{\ell=2}^{\infty}
\sum_{m=-\ell}^{\ell}
\sqrt{
\frac{(\ell+2)!}{(\ell-2)!}
}
\left(
\Psi_{\ell m}^{\mathrm{even}}
+i\Psi_{\ell m}^{\mathrm{odd}}
\right)
{}_{-2}Y_{\ell m}(\iota,\varphi)
+
O(R^{-2}).
\]

Aquí:

- \(R\) es la distancia al observador;
- \((\iota,\varphi)\) es la orientación del observador;
- \({}_{-2}Y_{\ell m}\) son armónicos esféricos de peso de espín \(-2\);
- el signo de \(i\) puede cambiar con la convención usada para \(h_\times\) y los armónicos.

El *strain* será adimensional. En unidades físicas,

\[
h\sim
\frac{GM_{\mathrm{BH}}}{c^2D}
\,H_{\mathrm{code}},
\]

donde \(H_{\mathrm{code}}\) es la combinación adimensional de modos calculada por la simulación.

La frecuencia se convierte mediante

\[
f_{\mathrm{Hz}}
=
f_{\mathrm{code}}
\frac{c^3}{GM_{\mathrm{BH}}},
\]

si \(f_{\mathrm{code}}\) está medida en ciclos por \(M\). Si se utiliza frecuencia angular \(\omega\), se incluye el factor \(2\pi\).

## Cómo se construyen las fuentes

En cada radio se proyecta el tensor de energía-momento:

\[
T^{ab}_{\ell m}(t,r)
=
\int
T^{ab}(t,r,\theta,\phi)
\overline{Y}_{\ell m}(\theta,\phi)
\,d\Omega,
\]

y de manera análoga:

\[
T^{aB}
\longrightarrow
\text{armónicos vectoriales},
\]

\[
T^{AB}
\longrightarrow
\text{armónicos tensoriales pares e impares}.
\]

Esas proyecciones se combinan para obtener

\[
S_{\ell m}^{\mathrm{even}}(t,r),
\qquad
S_{\ell m}^{\mathrm{odd}}(t,r).
\]

Para la PPI convendría comenzar con

\[
2\leq\ell\leq4,
\qquad
-\ell\leq m\leq\ell.
\]

No conviene guardar solamente \((\ell,m)=(2,2)\). El modo hidrodinámico dominante \(m=1\) puede producir señal en varios multipolos y generar armónicos adicionales durante la fase no lineal.

## Limitación importante de la malla 2.5D actual

La corrida ecuatorial tiene una sola celda en \(\theta\). Eso es suficiente para estudiar la dinámica ecuatorial, pero no proporciona la información angular necesaria para evaluar exactamente

\[
\int d\Omega
=
\int_0^\pi\sin\theta\,d\theta
\int_0^{2\pi}d\phi.
\]

Tratar esa celda como si representara todo \([0,\pi]\) no reconstruye la estructura vertical de un toro grueso. Esto explica una parte de la arbitrariedad de la amplitud actual.

Hay tres posibilidades:

1. **Simulación tridimensional o con resolución finita en \(\theta\).**

   Es la opción físicamente correcta y la recomendable para obtener amplitudes publicables.

2. **Modelo de disco infinitamente delgado:**

   \[
   T_{\mu\nu}(r,\theta,\phi)
   =
   \Sigma_{\mu\nu}(r,\phi)
   \frac{\delta(\theta-\pi/2)}{r}.
   \]

   Permite calcular proyecciones analíticas en \(\theta=\pi/2\), pero no representa un toro grueso.

3. **Reconstrucción vertical prescrita:**

   \[
   T_{\mu\nu}(r,\theta,\phi)
   =
   T_{\mu\nu}^{\mathrm{eq}}(r,\phi)F(r,\theta).
   \]

   Podría usarse un perfil vertical de equilibrio Fishbone–Moncrief, pero sería un modelo adicional, no una extracción exacta de la simulación.

Para el congreso, el diagnóstico cuadrupolar puede continuar claramente etiquetado como una estimación. Para una publicación de amplitudes de GW, necesitaríamos al menos la segunda opción bien justificada y, preferentemente, la primera.

## Normalización física de la masa del toro

Antes de hablar de un *strain* observable, necesitamos definir

\[
q=\frac{M_{\mathrm{torus}}}{M_{\mathrm{BH}}}.
\]

Actualmente la densidad del toro está normalizada principalmente para la evolución hidrodinámica. La “masa total” de los diagnósticos no debe interpretarse directamente en unidades de masas del agujero negro.

En una EOS gamma-law, una transformación

\[
\rho\rightarrow C\rho,
\qquad
p\rightarrow Cp
\]

mantiene esencialmente la misma dinámica *test-fluid*, pero cambia linealmente el tensor \(T_{\mu\nu}\) y, por tanto,

\[
\Psi_{\ell m}\rightarrow C\Psi_{\ell m},
\qquad
h\rightarrow Ch.
\]

Deberíamos elegir explícitamente uno o varios valores de \(q\), usando una secuencia de toros ligeros, y comprobar que:

- la dinámica hidrodinámica permanece igual;
- la amplitud de la perturbación métrica escala linealmente con \(q\);
- la energía radiada escala como \(q^2\).

Si el toro es demasiado masivo, la aproximación de Cowling y la perturbación lineal dejan de ser apropiadas. En particular, el movimiento del agujero negro alrededor del centro de masa y la auto-gravedad del toro pueden modificar considerablemente una PPI dominada por \(m=1\).

## Condiciones iniciales y *junk radiation*

No conviene iniciar las ecuaciones RWZ con

\[
\Psi_{\ell m}=0,
\qquad
\partial_t\Psi_{\ell m}=0
\]

si al mismo tiempo existe una fuente no nula desde \(t=0\). Esto generaría una ráfaga artificial porque la perturbación inicial no satisface la respuesta estacionaria al toro.

Las opciones incluyen:

- Resolver primero la ecuación estacionaria para la parte axisimétrica.
- Encender la fuente suavemente.
- Usar como fuente radiativa solamente la parte no axisimétrica:

  \[
  \delta T_{\mu\nu}
  =
  T_{\mu\nu}
  -
  \langle T_{\mu\nu}\rangle_\phi.
  \]

- Iniciar los modos \(m\neq0\) en cero mientras el toro es exactamente axisimétrico.

La última opción encaja especialmente bien con nuestras pruebas: antes de la perturbación, todos los modos \(m\neq0\) deberían permanecer en cero dentro del error numérico.

## Implementación numérica recomendada

No convendría poner toda esta evolución dentro de `output.f90`, porque RWZ introduce variables dinámicas propias. Una arquitectura más limpia sería:

```text
gw_perturbations.f90
├── inicialización de Ψeven/odd
├── proyección de Tμν
├── construcción de Seven/odd
├── evolución de las ecuaciones maestras
└── condiciones de frontera

output.f90
├── escritura de Ψ_lm
├── h+ y h× por orientación
├── flujo de energía
└── flujo de momento angular
```

Las ecuaciones maestras pueden evolucionarse con:

- reducción hiperbólica a primer orden;
- SSPRK3, igual que GRHD;
- diferencias centradas de cuarto orden, porque la solución es suave;
- disipación de Kreiss–Oliger pequeña;
- condición puramente *ingoing* en el horizonte;
- condición *outgoing* en la frontera externa.

La fuente debe evaluarse en las etapas del RK si queremos conservar el orden temporal. Usar solamente los VTK cada \(10M\) sería suficiente para un espectro preliminar, pero no para una extracción precisa.

También convendría emplear:

- varias superficies de extracción;
- tiempo retardado \(u=t-r_*\);
- extrapolación a \(R\rightarrow\infty\), o una coordenada hiperbólica compactificada;
- una frontera externa mucho más lejana que la usada actualmente para el toro.

## Extensión a Kerr: ecuación de Teukolsky

RWZ funciona directamente para Schwarzschild. Para Kerr, el formalismo natural es Teukolsky.

En una tetrada de Newman–Penrose adaptada a las direcciones nulas principales de Kerr, la radiación saliente queda codificada en

\[
\delta\Psi_4
=
-\delta C_{\alpha\beta\gamma\delta}
n^\alpha\bar m^\beta
n^\gamma\bar m^\delta.
\]

La variable de espín \(s=-2\), frecuentemente una versión reescalada de \(\delta\Psi_4\), satisface

\[
\mathcal T_{-2}[\psi]
=
4\pi\,\mathcal S[T_{\mu\nu}].
\]

La fuente involucra proyecciones de materia sobre la tetrada:

\[
T_{nn},
\qquad
T_{n\bar m},
\qquad
T_{\bar m\bar m},
\]

además de operadores diferenciales.

Teukolsky derivó las ecuaciones separables con fuentes para perturbaciones de Kerr. Véase [Teukolsky (1972)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.29.1114). La construcción de fuentes extendidas es considerablemente más difícil que en Schwarzschild; trabajos recientes incluso han demostrado separaciones para una nube de gas ideal alrededor de Kerr. Véase [Spiers (2024)](https://arxiv.org/abs/2402.00604).

Asintóticamente,

\[
\Psi_4
\simeq
\ddot h_+-i\ddot h_\times,
\]

hasta signos y factores dependientes de la convención. Entonces

\[
h_+-ih_\times
=
\int^t dt'
\int^{t'}dt''
\,\Psi_4(t'').
\]

La doble integración directa suele producir deriva secular. En la práctica se utiliza integración en frecuencia:

\[
\tilde h(\omega)
=
-\frac{\tilde\Psi_4(\omega)}{\omega^2},
\]

con un tratamiento controlado de las frecuencias bajas.

Una precisión conceptual: \(\delta\Psi_4\) es gauge-invariant a primer orden sobre Kerr porque la \(\Psi_4\) del fondo se anula en la tetrada principal, pero todavía existe libertad de tetrada. Hay que fijar una tetrada tipo Kinnersley, preferentemente en una versión regular en el horizonte, y llevarla a una tetrada de Bondi en la zona asintótica.

## Métrica dinámica y CCE: el nivel final

Si en el futuro evolucionamos la métrica con BSSN, CCZ4 u otra formulación, podremos construir

\[
\Psi_4
=
(E_{ij}-iB_{ij})\bar m^i\bar m^j,
\]

donde \(E_{ij}\) y \(B_{ij}\) son las partes eléctrica y magnética del tensor de Weyl.

Después:

1. Se extrae \(\Psi_4\) en varias esferas.
2. Se descompone en \({}_{-2}Y_{\ell m}\).
3. Se extrapola a infinito o se utiliza CCE.
4. Se integra para obtener \(h_+\) y \(h_\times\).

La extracción directa de \(\Psi_4\) a radio finito puede conservar efectos de gauge y de campo cercano. CCE es la opción más rigurosa porque propaga los datos hasta \(\mathscr I^+\). Comparaciones directas muestran que CCE puede eliminar diferencias producidas por gauges distintos que permanecen en *waveforms* extrapoladas desde radio finito. Véase [Taylor et al. (2013)](https://arxiv.org/abs/1309.3605).

Simulaciones completamente relativistas de PPI ya han empleado \(\Psi_4\) para identificar la emisión cuasiperiódica asociada con la estructura no axisimétrica. Véase [Kiuchi et al. (2011)](https://arxiv.org/abs/1105.5035). Ese es un buen punto de comparación futuro, aunque su sistema incluye auto-gravedad y métrica dinámica.

## Ruta recomendada para GRHD

1. Mantener por ahora el `gw_strain` cuadrupolar como *proxy* de frecuencia y fase, renombrándolo conceptualmente como `gw_quadrupole_proxy`.

2. Usar la comparación EF contra KS con \(a=0\) ya completada como baseline y
   repetirla con observables gauge-invariant cuando exista el módulo RWZ.

3. Definir una normalización física \(M_{\mathrm{torus}}/M_{\mathrm{BH}}\).

4. Implementar proyecciones de \(T_{\mu\nu}\) y comprobarlas con distribuciones angulares analíticas.

5. Implementar RWZ con fuente para Schwarzschild.

6. Validarlo con:

   - pulso gaussiano sin fuente;
   - *ringdown* y cola de Schwarzschild;
   - convergencia de las variables maestras;
   - ausencia de radiación para una fuente estacionaria y axisimétrica;
   - coincidencia EF–KS con \(a=0\);
   - independencia del radio de extracción;
   - comparación con el cuadrupolo en el límite débil.

7. Aplicarlo primero a una simulación con resolución real en \(\theta\), aunque sea modesta.

8. Implementar Teukolsky con fuente solamente después de dominar RWZ.

## Conclusión

No es indispensable evolucionar una métrica dinámica para obtener una señal lineal gauge-invariant. Podemos mantener GRHD sobre un fondo fijo y resolver RWZ o Teukolsky con \(T_{\mu\nu}\) como fuente.

Esta aproximación no incluiría:

- retroacción de las ondas sobre el fluido;
- auto-gravedad del toro;
- movimiento consistente del agujero negro;
- modificaciones no lineales de la geometría.

Para toros suficientemente ligeros, la aproximación perturbativa puede ser rigurosa y considerablemente más alcanzable para esta tesis que implementar inmediatamente una evolución completa del espacio-tiempo.

## Referencias principales

1. K. Martel y E. Poisson, *Gravitational perturbations of the Schwarzschild spacetime: A practical covariant and gauge-invariant formalism*, Phys. Rev. D 71, 104003 (2005). [arXiv:gr-qc/0502028](https://arxiv.org/abs/gr-qc/0502028).
2. O. Sarbach y M. Tiglio, *Gauge invariant perturbations of Schwarzschild black holes in horizon-penetrating coordinates*, Phys. Rev. D 64, 084016 (2001). [arXiv:gr-qc/0104061](https://arxiv.org/abs/gr-qc/0104061).
3. S. A. Teukolsky, *Rotating Black Holes: Separable Wave Equations for Gravitational and Electromagnetic Perturbations*, Phys. Rev. Lett. 29, 1114 (1972). [DOI](https://doi.org/10.1103/PhysRevLett.29.1114).
4. C. Reisswig, N. T. Bishop, D. Pollney y B. Szilágyi, *Characteristic extraction in numerical relativity: binary black hole merger waveforms at null infinity* (2009). [arXiv:0912.1285](https://arxiv.org/abs/0912.1285).
5. N. W. Taylor et al., *Comparing Gravitational Waveform Extrapolation to Cauchy-Characteristic Extraction in Binary Black Hole Simulations*, Phys. Rev. D 88, 124010 (2013). [arXiv:1309.3605](https://arxiv.org/abs/1309.3605).
6. K. Kiuchi, M. Shibata, P. J. Montero y J. A. Font, *Gravitational waves from the Papaloizou-Pringle instability in black hole-torus systems* (2011). [arXiv:1105.5035](https://arxiv.org/abs/1105.5035).
7. A. Spiers, *Analytically Separating the Source of the Teukolsky Equation* (2024). [arXiv:2402.00604](https://arxiv.org/abs/2402.00604).
