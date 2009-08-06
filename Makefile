FC = /usr/bin/gfortran
PROGRAM = raygun
OBJECTS = raygun.o
#add -v here to see shenanigans!
VERBOSE = 

LIB_TAG = d
PKG_CONFIG_ENV = 
MODDIR = -J xmllib

$(PROGRAM): raygun.f90 xmllib/xmlparse.a xmllib/xmlparse.mod
	echo "BUILDIFYING!"
	$(FC) $(VERBOSE) raygun.f90 -o $@ $(MODDIR) \
	`$(PKG_CONFIG_ENV) pkg-config --cflags --libs plplot$(LIB_TAG)-f95` \
	xmllib/xmlparse.a

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
