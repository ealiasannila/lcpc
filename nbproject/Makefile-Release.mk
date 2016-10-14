#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/lib/clipper/cpp/clipper.o \
	${OBJECTDIR}/lib/common/shapes.o \
	${OBJECTDIR}/lib/sweep/advancing_front.o \
	${OBJECTDIR}/lib/sweep/cdt.o \
	${OBJECTDIR}/lib/sweep/sweep.o \
	${OBJECTDIR}/lib/sweep/sweep_context.o \
	${OBJECTDIR}/src/Coords.o \
	${OBJECTDIR}/src/Funnel.o \
	${OBJECTDIR}/src/lcpc.o \
	${OBJECTDIR}/src/lcpfinder.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/lcpc

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/lcpc: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/lcpc ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/lib/clipper/cpp/clipper.o: lib/clipper/cpp/clipper.cpp 
	${MKDIR} -p ${OBJECTDIR}/lib/clipper/cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lib/clipper/cpp/clipper.o lib/clipper/cpp/clipper.cpp

${OBJECTDIR}/lib/common/shapes.o: lib/common/shapes.cc 
	${MKDIR} -p ${OBJECTDIR}/lib/common
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lib/common/shapes.o lib/common/shapes.cc

${OBJECTDIR}/lib/sweep/advancing_front.o: lib/sweep/advancing_front.cc 
	${MKDIR} -p ${OBJECTDIR}/lib/sweep
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lib/sweep/advancing_front.o lib/sweep/advancing_front.cc

${OBJECTDIR}/lib/sweep/cdt.o: lib/sweep/cdt.cc 
	${MKDIR} -p ${OBJECTDIR}/lib/sweep
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lib/sweep/cdt.o lib/sweep/cdt.cc

${OBJECTDIR}/lib/sweep/sweep.o: lib/sweep/sweep.cc 
	${MKDIR} -p ${OBJECTDIR}/lib/sweep
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lib/sweep/sweep.o lib/sweep/sweep.cc

${OBJECTDIR}/lib/sweep/sweep_context.o: lib/sweep/sweep_context.cc 
	${MKDIR} -p ${OBJECTDIR}/lib/sweep
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lib/sweep/sweep_context.o lib/sweep/sweep_context.cc

${OBJECTDIR}/src/Coords.o: src/Coords.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Coords.o src/Coords.cpp

${OBJECTDIR}/src/Funnel.o: src/Funnel.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Funnel.o src/Funnel.cpp

${OBJECTDIR}/src/lcpc.o: src/lcpc.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lcpc.o src/lcpc.cpp

${OBJECTDIR}/src/lcpfinder.o: src/lcpfinder.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lcpfinder.o src/lcpfinder.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/lcpc

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
