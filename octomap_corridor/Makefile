CXXFLAGS = -std=c++11 -I/opt/ros/kinetic/include -L/opt/ros/kinetic/lib -O3
CFLAGS = -O3 -x c

CC = $(CXX)

octomap_corridor : octomap_corridor.o \
/opt/ros/kinetic/lib/liboctomap.a /opt/ros/kinetic/lib/liboctomath.a \
cvxgen/ldl.o cvxgen/matrix_support.o cvxgen/solver.o cvxgen/util.o

clean:
	rm -rf cvxgen/*.o
	rm -rf octomap_corridor *.o
