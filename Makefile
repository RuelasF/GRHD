# ==========================================
# Makefile para Simulador GRHD2
# ==========================================

# Compilador y Banderas de optimización
FC = gfortran
ARCH_FLAGS ?= -march=native
FFLAGS = -O3 -fopenmp $(ARCH_FLAGS) -flto=auto -fno-math-errno -fno-trapping-math -ffree-line-length-none -cpp
#FFLAGS = -fopenmp -O0 -g -fbacktrace -Wall -Wextra -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow,underflow,denormal

# Nombre del ejecutable final
TARGET = grhd2
.DEFAULT_GOAL := $(TARGET)
TEST_HLLC_TARGET = tests/test_hllc_fluxes
TEST_GW_GEOMETRY_TARGET = tests/test_gw_geometry
TEST_GW_STRESS_TARGET = tests/test_gw_stress
TEST_FM_TARGET = tests/test_fm_initial
TEST_POLAR_TARGET = tests/test_polar_boundaries
TEST_CFL_TARGET = tests/test_multidimensional_cfl
TEST_MINKOWSKI_SPHERICAL_TARGET = tests/test_minkowski_spherical
TEST_CASE_CORRECTIONS_TARGET = tests/test_case_corrections
TEST_CHECKPOINT_TARGET = tests/test_checkpoint_v5
TEST_OUTPUT_TARGET = tests/test_output_paths
TEST_ACC_RECON_TARGET = tests/test_accelerator_reconstruction
TEST_ACC_RHS_TARGET = tests/test_accelerator_rhs

.PHONY: server-broadwell openacc-p100 openacc-rtx4070 test test-hllc test-gw test-gw-geometry test-gw-stress test-fm test-polar test-cfl test-minkowski-spherical test-cases test-checkpoint test-output test-acc-reconstruction test-acc-rhs clean

# Ejecutable transferible para Intel Xeon E5 v4 (Broadwell-EP).
server-broadwell:
	$(MAKE) clean
	$(MAKE) ARCH_FLAGS='-march=broadwell -mtune=broadwell' $(TARGET)

# NVIDIA Tesla P100 = Pascal, compute capability 6.0. No se usa -fast ni
# precisiones mixtas: todos los real*8 permanecen en IEEE FP64.
openacc-p100:
	$(MAKE) clean
	$(MAKE) FC=nvfortran FFLAGS='-O3 -acc=gpu -gpu=cc60,lineinfo -Minfo=accel -Mfree -Mextend -cpp -DUSE_OPENACC' $(TARGET)

# Compilación funcional para Zotz/RTX 4070 (Ada, CC 8.9). Conserva FP64;
# sirve para validar kernels, pero no representa el rendimiento FP64 de P100.
openacc-rtx4070:
	$(MAKE) clean
	$(MAKE) FC=nvfortran FFLAGS='-O3 -acc=gpu -gpu=cc89,lineinfo -Minfo=accel -Mfree -Mextend -cpp -DUSE_OPENACC' $(TARGET)

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
       accelerator.o \
       grhd2.o

# Regla principal para construir el ejecutable
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS)

test: test-hllc test-gw test-fm test-polar test-cfl test-minkowski-spherical test-cases test-checkpoint test-output test-acc-reconstruction test-acc-rhs

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

test-checkpoint: variables.o metrics.o conditions.o output.o tests/test_checkpoint_v5.f90
	$(FC) $(FFLAGS) -o $(TEST_CHECKPOINT_TARGET) tests/test_checkpoint_v5.f90 variables.o metrics.o conditions.o output.o
	./$(TEST_CHECKPOINT_TARGET)

test-output: variables.o metrics.o conditions.o output.o tests/test_output_paths.f90
	$(FC) $(FFLAGS) -o $(TEST_OUTPUT_TARGET) tests/test_output_paths.f90 variables.o metrics.o conditions.o output.o
	./$(TEST_OUTPUT_TARGET)

test-acc-reconstruction: variables.o reconstruction.o tests/test_accelerator_reconstruction.f90
	$(FC) $(FFLAGS) -o $(TEST_ACC_RECON_TARGET) tests/test_accelerator_reconstruction.f90 variables.o reconstruction.o
	./$(TEST_ACC_RECON_TARGET)

test-acc-rhs: variables.o metrics.o equations.o fluxes.o reconstruction.o evolution.o tests/test_accelerator_rhs.f90
	$(FC) $(FFLAGS) -o $(TEST_ACC_RHS_TARGET) tests/test_accelerator_rhs.f90 variables.o metrics.o equations.o fluxes.o reconstruction.o evolution.o
	./$(TEST_ACC_RHS_TARGET)

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
accelerator.o: variables.o
grhd2.o: variables.o parameters.o metrics.o initialization.o conditions.o equations.o fluxes.o reconstruction.o evolution.o accelerator.o output.o

# Limpiar archivos compilados (comando: make clean)
clean:
	rm -f *.o *.mod $(TARGET) $(TEST_HLLC_TARGET) $(TEST_GW_GEOMETRY_TARGET) $(TEST_GW_STRESS_TARGET) $(TEST_FM_TARGET) $(TEST_POLAR_TARGET) $(TEST_CFL_TARGET) $(TEST_MINKOWSKI_SPHERICAL_TARGET) $(TEST_CASE_CORRECTIONS_TARGET) $(TEST_CHECKPOINT_TARGET) $(TEST_OUTPUT_TARGET) $(TEST_ACC_RECON_TARGET) $(TEST_ACC_RHS_TARGET)
