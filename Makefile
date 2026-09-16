# ==========================================
# Makefile para Simulador GRHD2
# ==========================================

# Compilador y Banderas de optimización
FC = gfortran
FFLAGS = -O3 -fopenmp -march=native -flto=auto -fno-math-errno -fno-trapping-math -ffree-line-length-none
#FFLAGS = -fopenmp -O0 -g -fbacktrace -Wall -Wextra -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow,underflow,denormal

# Nombre del ejecutable final
TARGET = grhd2
TEST_HLLC_TARGET = tests/test_hllc_fluxes
TEST_GW_GEOMETRY_TARGET = tests/test_gw_geometry
TEST_GW_STRESS_TARGET = tests/test_gw_stress
TEST_FM_TARGET = tests/test_fm_initial
TEST_POLAR_TARGET = tests/test_polar_boundaries
TEST_CFL_TARGET = tests/test_multidimensional_cfl
TEST_MINKOWSKI_SPHERICAL_TARGET = tests/test_minkowski_spherical
TEST_CASE_CORRECTIONS_TARGET = tests/test_case_corrections
TEST_CHECKPOINT_TARGET = tests/test_checkpoint_v4
TEST_OUTPUT_TARGET = tests/test_output_paths

.PHONY: test test-hllc test-gw test-gw-geometry test-gw-stress test-fm test-polar test-cfl test-minkowski-spherical test-cases test-checkpoint test-output clean

# Lista de archivos objeto en el orden estricto de dependencias
OBJS = variables.o \
       parameters.o \
       metrics.o \
       output.o \
       initialization.o \
       conditions.o \
       equations.o \
       fluxes.o \
       reconstruction.o \
       evolution.o \
       grhd2.o

# Regla principal para construir el ejecutable
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS)

test: test-hllc test-gw test-fm test-polar test-cfl test-minkowski-spherical test-cases test-checkpoint test-output

test-hllc: variables.o metrics.o equations.o fluxes.o tests/test_hllc_fluxes.f90
	$(FC) $(FFLAGS) -o $(TEST_HLLC_TARGET) tests/test_hllc_fluxes.f90 variables.o metrics.o equations.o fluxes.o
	./$(TEST_HLLC_TARGET)

test-gw-geometry: variables.o metrics.o conditions.o output.o tests/test_gw_geometry.f90
	$(FC) $(FFLAGS) -o $(TEST_GW_GEOMETRY_TARGET) tests/test_gw_geometry.f90 variables.o metrics.o conditions.o output.o
	./$(TEST_GW_GEOMETRY_TARGET)

test-gw-stress: variables.o metrics.o conditions.o output.o tests/test_gw_stress.f90
	$(FC) $(FFLAGS) -o $(TEST_GW_STRESS_TARGET) tests/test_gw_stress.f90 variables.o metrics.o conditions.o output.o
	./$(TEST_GW_STRESS_TARGET)

test-gw: test-gw-geometry test-gw-stress

test-fm: variables.o metrics.o conditions.o tests/test_fm_initial.f90
	$(FC) $(FFLAGS) -o $(TEST_FM_TARGET) tests/test_fm_initial.f90 variables.o metrics.o conditions.o
	./$(TEST_FM_TARGET)

test-polar: variables.o metrics.o conditions.o tests/test_polar_boundaries.f90
	$(FC) $(FFLAGS) -o $(TEST_POLAR_TARGET) tests/test_polar_boundaries.f90 variables.o metrics.o conditions.o
	./$(TEST_POLAR_TARGET)

test-cfl: variables.o metrics.o equations.o fluxes.o reconstruction.o evolution.o tests/test_multidimensional_cfl.f90
	$(FC) $(FFLAGS) -o $(TEST_CFL_TARGET) tests/test_multidimensional_cfl.f90 variables.o metrics.o equations.o fluxes.o reconstruction.o evolution.o
	./$(TEST_CFL_TARGET)

test-minkowski-spherical: variables.o metrics.o tests/test_minkowski_spherical.f90
	$(FC) $(FFLAGS) -o $(TEST_MINKOWSKI_SPHERICAL_TARGET) tests/test_minkowski_spherical.f90 variables.o metrics.o
	./$(TEST_MINKOWSKI_SPHERICAL_TARGET)

test-cases: variables.o parameters.o metrics.o conditions.o equations.o fluxes.o initialization.o tests/test_case_corrections.f90
	$(FC) $(FFLAGS) -o $(TEST_CASE_CORRECTIONS_TARGET) tests/test_case_corrections.f90 variables.o parameters.o metrics.o conditions.o equations.o fluxes.o initialization.o
	./$(TEST_CASE_CORRECTIONS_TARGET)

test-checkpoint: variables.o metrics.o conditions.o output.o tests/test_checkpoint_v4.f90
	$(FC) $(FFLAGS) -o $(TEST_CHECKPOINT_TARGET) tests/test_checkpoint_v4.f90 variables.o metrics.o conditions.o output.o
	./$(TEST_CHECKPOINT_TARGET)

test-output: variables.o metrics.o conditions.o output.o tests/test_output_paths.f90
	$(FC) $(FFLAGS) -o $(TEST_OUTPUT_TARGET) tests/test_output_paths.f90 variables.o metrics.o conditions.o output.o
	./$(TEST_OUTPUT_TARGET)

# Regla genérica para compilar archivos .f90 a objetos .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# --- DEPENDENCIAS DE LOS MÓDULOS ---
# Esto le dice a Make qué debe compilarse antes que otra cosa
metrics.o: variables.o
parameters.o: variables.o
output.o: variables.o metrics.o conditions.o
initialization.o: variables.o parameters.o metrics.o conditions.o fluxes.o
conditions.o: variables.o metrics.o
equations.o: variables.o metrics.o
fluxes.o: variables.o equations.o metrics.o
reconstruction.o: variables.o
evolution.o: variables.o reconstruction.o fluxes.o equations.o
grhd2.o: variables.o parameters.o metrics.o initialization.o conditions.o equations.o fluxes.o reconstruction.o evolution.o output.o

# Limpiar archivos compilados (comando: make clean)
clean:
	rm -f *.o *.mod $(TARGET) $(TEST_HLLC_TARGET) $(TEST_GW_GEOMETRY_TARGET) $(TEST_GW_STRESS_TARGET) $(TEST_FM_TARGET) $(TEST_POLAR_TARGET) $(TEST_CFL_TARGET) $(TEST_MINKOWSKI_SPHERICAL_TARGET) $(TEST_CASE_CORRECTIONS_TARGET) $(TEST_CHECKPOINT_TARGET) $(TEST_OUTPUT_TARGET)
