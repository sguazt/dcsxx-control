## DCS USER CONFIGURATIONS
##


### dcs general configurations

# Set here you C/C++ preprocessor flags
USR_CPPFLAGS=
# Set here you C compiler flags
USR_CFLAGS=
# Set here you C++ compiler flags
#USR_CXXFLAGS=
#USR_CXXFLAGS=-I../dcsxx/src
USR_CXXFLAGS="-I../dcsxx/src -I../boost-ublas -I$HOME/tmp/cvs-prj/boost-trunk"
#USR_CXXFLAGS="-I../dcsxx/src -I../boost-ublas -I$HOME/Sys/usr/include/boost"
# Set here you C/C++ linker flags
USR_LDFLAGS=

export USR_CPPFLAGS USR_CFLAGS USR_CXXFLAGS USR_LDFLAGS
