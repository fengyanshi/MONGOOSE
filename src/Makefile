#
# Makefile for mongoose and bubbles
# 
# Assign some variables
#
F90=gfortran
FFLAGS=
#
# Assign some path variables 
#
OBJ_DIR=./Objects
RES_DIR=../bin
#
# Assign some commands
#
.PHONY: all clean
#
# The object files
#
OBJECTS= bc.o bubble_transf.o ke.o main.o read.o surface.o water.o
#
# This is the main dependency
#
all: ${OBJECTS} makedir mongoose

makedir:
	@if [ ! -d ${OBJ_DIR} ]; then mkdir ${OBJ_DIR} ; fi
	@mv ${OBJECTS} ${OBJ_DIR}
	@if [ ! -d ${RES_DIR} ]; then mkdir ${RES_DIR} ; fi

mongoose: ${OBJECTS}
	@echo 'Objects ===> mongoose'
	@cd ${OBJ_DIR}; \
	${F90} ${FFLAGS} ${OBJECTS} -o mongoose; \
	mv mongoose ../${RES_DIR}; \
	cd ..
	@echo

${OBJECTS}: %.o : %.f
	@${F90} -c ${FFLAGS} $< -o $@

#
# This is the clean dependency
#
clean:
	@echo 'Erasing object files'
	@rm -rf ${OBJ_DIR}

	@echo 'Erasing old files'
	@rm -f *~

	@echo 'Erasing executables'
	@rm -f mongoose

veryclean:  clean
	@echo 'Erasing results'
	@rm -rf ${RES_DIR}
