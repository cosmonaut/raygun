FC = /usr/bin/gfortran
PROGRAM = raygun
OBJECTS = raygun.o 
#add -v here to see shenanigans!
VERBOSE = 
ARGS = -fdefault-real-8 -fdefault-double-8

LIB_TAG = d
PKG_CONFIG_ENV = 
MODDIR = -J xmllib

$(PROGRAM): raygun.f90 xml_data_config.o
	echo "BUILDIFYING!"
	$(FC) $(ARGS) $(VERBOSE) raygun.f90 -o $@ $(MODDIR) xml_data_config.o \
	`$(PKG_CONFIG_ENV) pkg-config --cflags --libs plplot$(LIB_TAG)-f95` \
	xmllib/xmlparse.a 

xml_data_config.mod: xml_data_config.o xml_data_config.f90 
	$(FC) -c $(MODDIR) xml_data_config.o

xml_data_config.o: xml_data_config.f90
	$(FC) -c $(MODDIR) xml_data_config.f90

# raygun: raygun.o
# 	$(FC) $(VERBOSE) raygun.o -o raygun

# raygun.o: raygun.f90
# 	$(FC) -c $(VERBOSE) raygun.f90

#%.o: %.f90
#	$(FC) -c $@ $<

clean:
	@if [ -e $(OBJECTS) ]; then \
		rm $(OBJECTS); \
	fi
	@if [ -e *\~ ]; then \
		rm *\~ ; \
	fi
	@if [ -e $(PROGRAM) ]; then \
		rm $(PROGRAM) ; \
	fi
	@if [ -e *.o ]; then \
		rm *.o ; \
	fi
