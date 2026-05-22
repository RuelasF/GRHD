# ==========================================
# Makefile para Simulador GRHD2
# ==========================================

# Compilador y Banderas de optimización
FC = gfortran
FFLAGS = -O3 -fopenmp -march=native -flto=auto -fno-math-errno -fno-trapping-math -ffree-line-length-none
#FFLAGS = -fopenmp -O0 -g -fbacktrace -Wall -Wextra -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow,underflow,denormal

# Nombre del ejecutable final
TARGET = grhd2

# Lista de archivos objeto en el orden estricto de dependencias
OBJS = variables.o \
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

# Regla genérica para compilar archivos .f90 a objetos .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# --- DEPENDENCIAS DE LOS MÓDULOS ---
# Esto le dice a Make qué debe compilarse antes que otra cosa
metrics.o: variables.o
output.o: variables.o metrics.o conditions.o
initialization.o: variables.o metrics.o conditions.o fluxes.o
conditions.o: variables.o metrics.o
equations.o: variables.o metrics.o
fluxes.o: variables.o equations.o metrics.o
reconstruction.o: variables.o
evolution.o: variables.o reconstruction.o fluxes.o equations.o
grhd2.o: variables.o metrics.o initialization.o conditions.o equations.o fluxes.o reconstruction.o evolution.o output.o

# Limpiar archivos compilados (comando: make clean)
clean:
	rm -f *.o *.mod $(TARGET)