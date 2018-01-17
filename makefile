CC         = icpc
#CFLAGS     = -pg -fopenmp
CFLAGS     = -O3 -openmp -math -Wall -I/projects/netpub/fftw_openmp/3.3.4/include
LIBS      = -openmp  -lm  -O3 -lfftw3_omp -lfftw3 -lpthread -L/projects/netpub/fftw_openmp/3.3.4/lib
#CFLAGS     = -O3 -openmp -math -Wall -I/home/hchao/Install/fftw3_openmp/include
#LIBS      = -openmp  -lm  -O3 -lfftw3_omp -lfftw3 -lpthread -L/home/hchao/Install/fftw3_openmp/lib


#############################################################################
# nothing should be changed below here

SRCS = time.cpp stress.cpp main.cpp matrix.cpp array_utils.cpp die.cpp  random.cpp grid_utils.cpp \
			 torque.cpp quanterions.cpp fftw_wrappers.cpp initialize.cpp config_utils.cpp io_utils.cpp \
		update_euler.cpp update_positions.cpp forces.cpp integ_utils.cpp read_input.cpp \
			 bonded.cpp calc_unb.cpp 
       
       
			 


OBJS = ${SRCS:.cpp=.o}

.cpp.o:
	${CC} ${CFLAGS} ${DFLAGS} -c  $<

dmft:  ${OBJS}
	$(CC) ${CFLAGS} ${DFLAGS} -o $@ ${OBJS} $(LIBS)

clean:
	rm -f *.o
	rm -f dmft
	rm -f *~

