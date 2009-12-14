FC = /usr/bin/gfortran
PROGRAM = raygun
OBJECTS = raygun.o 
#add -v here to see shenanigans!
VERBOSE = 
#ARGS = -Wall -fdefault-real-8 -fdefault-double-8
ARGS = -Wall -fdefault-real-8 -fdefault-double-8

LIB_TAG = d
PKG_CONFIG_ENV = 
MODDIR = -J xmllib


$(PROGRAM): raygun.f90 xml_data_config.o constants_NSWC.o quart.o
	echo "BUILDIFYING!"
	$(FC) $(ARGS) $(VERBOSE) raygun.f90 -o $@ $(MODDIR) constants_NSWC.o quart.o xml_data_config.o \
	`$(PKG_CONFIG_ENV) pkg-config --cflags --libs plplot$(LIB_TAG)-f95` \
	xmllib/xmlparse.a 

xml_data_config.mod: xml_data_config.o xml_data_config.f90 
	$(FC) -c $(MODDIR) xml_data_config.o

xml_data_config.o: xml_data_config.f90
	$(FC) -c $(MODDIR) xml_data_config.f90

quart.mod: quart.o quart.f90 constants_NSWC.mod
	$(FC) -c $(MODDIR) quart.o 

quart.o: quart.f90 constants_NSWC.o constants_NSWC.mod
	$(FC) -c $(MODDIR) quart.f90 

constants_NSWC.mod: constants_NSWC.o constants_NSWC.f90 
	$(FC) -c $(MODDIR) constants_NSWC.o

constants_NSWC.o: constants_NSWC.f90
	$(FC) -c $(MODDIR) constants_NSWC.f90


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
	@if [ -e $(PROGRAM) ]; then \
		rm $(PROGRAM) ; \
	fi
	@if [ -e xmllib/xml_data_config.mod ]; then \
		rm xmllib/xml_data_config.mod; \
	fi
	@if [ -e xmllib/quart.mod ]; then \
		rm xmllib/quart.mod; \
	fi
	@if [ -e xmllib/constants_nswc.mod ]; then \
		rm xmllib/constants_nswc.mod; \
	fi
	-rm -f *.o
	-rm -f *\~
