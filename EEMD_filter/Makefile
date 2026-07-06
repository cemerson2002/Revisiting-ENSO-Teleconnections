##
## ref: https://chatgpt.com/share/67e3f713-63e4-800f-9b2a-c1ed2a3ccc63
## 

FC = gfortran -mcmodel=large -fPIC -o3 -g -fcheck=all -march=native -ffast-math -fopenmp
FFLAGS = -I/usr/lib64/gfortran/modules -lnetcdff
TARGET = main_prog.exe
SRC = netcdf_reader.f90 netcdf_writer_1x1.f90 spatiotemporal_filter.f90 main_prog.f90
OBJ = $(SRC:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) *.mod $(TARGET)
