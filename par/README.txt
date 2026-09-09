ARCHIVOS DE PARAMETROS DE GRHD2
================================

Ejecucion:

  ./grhd2 archivo.par
  ./grhd2 --check archivo.par

Cada linea activa tiene la forma nombre = valor. El caracter # inicia un comentario.
Los nombres y valores simbolicos no distinguen mayusculas de minusculas.

COORDENADAS ANGULARES
---------------------

theta_min, theta_max, phi_min y phi_max se escriben en unidades de pi. Se aceptan
decimales y cocientes. Por ejemplo:

  theta_min = 1/2       # pi/2
  theta_max = 1         # pi
  phi_min = 0
  phi_max = 2           # 2*pi

CONDICIONES INICIALES
---------------------

Solo se aceptan los parametros que corresponden al valor de problem.

problem = sod
  sod_interface
  sod_density_left, sod_density_right
  sod_pressure_left, sod_pressure_right

problem = strong_shock
  strong_interface
  strong_density_left, strong_density_right
  strong_pressure_left, strong_pressure_right

problem = shu_osher
  shu_interface
  shu_density_left, shu_pressure_left
  shu_density_right, shu_density_amplitude, shu_wave_number, shu_pressure_right

problem = kelvin_helmholtz
  khi_half_width
  khi_density_inner, khi_density_outer
  khi_velocity_inner, khi_velocity_outer
  khi_pressure, khi_perturbation_amplitude, khi_perturbation_wave_number
  El numero de onda de KHI multiplica a pi dentro del codigo.

problem = axisymmetric_jet
  jet_ambient_density, jet_ambient_pressure, jet_ambient_velocity
  jet_nozzle_radius, jet_nozzle_length
  jet_density, jet_pressure, jet_velocity

problem = convergence
  convergence_density, convergence_amplitude
  convergence_velocity, convergence_pressure, convergence_wave_number
  convergence_wave_number se expresa como coeficiente de pi.

problem = michel
  michel_critical_radius, michel_critical_density

problem = dust
  dust_accretion_constant
  Debe ser negativo para producir una densidad positiva con la convencion actual.

problem = offaxis_blast
  offaxis_radial_center, offaxis_phi_center, offaxis_width
  offaxis_background_density, offaxis_background_pressure
  offaxis_density_amplitude, offaxis_pressure_amplitude
  offaxis_phi_center se expresa en unidades de pi.

problem = fishbone_equatorial, fishbone_sagittal o fishbone_3d
  fm_inner_radius, fm_pressure_max_radius, fm_polytropic_constant
  Los dos radios se expresan en unidades de la masa M del agujero negro.

METODOS NUMERICOS
-----------------

  reconstruction = godunov | tvd | weno3 | mp5 | weno5
  tvd_limiter = minmod | superbee | mc
  riemann_solver = hlle | hllc

Los presets de initialization.f90 proporcionan los valores por omision. El archivo .par
solo necesita indicar problem y los valores que se quieran cambiar.

SALIDA VTK
----------

En geometria esferica, la posicion de los puntos se selecciona con:

  vtk_mapping = physical    # Mapeo cartesiano fisico; conserva la torsion KS.
  vtk_mapping = untwisted   # Mapeo oblato sin la rotacion radial de KS.

El valor por omision es physical. Esta opcion solo cambia la geometria mostrada en el
archivo VTK; no modifica las coordenadas ni el estado usados durante la evolucion.
