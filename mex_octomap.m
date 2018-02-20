% compiles the optional mex functions that use Octomap. "make octomap" runs this for you.

cd utils
mex CXXFLAGS='$CXXFLAGS -I/opt/ros/kinetic/include/ -std=c++11' read_octomap_bbox_mex.cpp /opt/ros/kinetic/lib/liboctomap.a /opt/ros/kinetic/lib/liboctomath.a
cd ..
