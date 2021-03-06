# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/cmake-3.8.2/bin/cmake

# The command to remove a file.
RM = /usr/local/cmake-3.8.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
#CMAKE_SOURCE_DIR = /home/sgt29/dihmm_py
CMAKE_SOURCE_DIR = /home/yk890/diHMM/dihmm_py

# The top-level build directory on which CMake was run.
#CMAKE_BINARY_DIR = /home/sgt29/dihmm_py
CMAKE_BINARY_DIR = /home/yk890/diHMM/dihmm_py

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/local/cmake-3.8.2/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/cmake-3.8.2/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/yk890/diHMM/dihmm_py/CMakeFiles /home/yk890/diHMM/dihmm_py/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/yk890/diHMM/dihmm_py/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named dihmm_ext

# Build rule for target.
dihmm_ext: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 dihmm_ext
.PHONY : dihmm_ext

# fast build rule for target.
dihmm_ext/fast:
	$(MAKE) -f CMakeFiles/dihmm_ext.dir/build.make CMakeFiles/dihmm_ext.dir/build
.PHONY : dihmm_ext/fast

#=============================================================================
# Target rules for targets named dihmm

# Build rule for target.
dihmm: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 dihmm
.PHONY : dihmm

# fast build rule for target.
dihmm/fast:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/build
.PHONY : dihmm/fast

Emissions.o: Emissions.cpp.o

.PHONY : Emissions.o

# target to build an object file
Emissions.cpp.o:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Emissions.cpp.o
.PHONY : Emissions.cpp.o

Emissions.i: Emissions.cpp.i

.PHONY : Emissions.i

# target to preprocess a source file
Emissions.cpp.i:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Emissions.cpp.i
.PHONY : Emissions.cpp.i

Emissions.s: Emissions.cpp.s

.PHONY : Emissions.s

# target to generate assembly for a file
Emissions.cpp.s:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Emissions.cpp.s
.PHONY : Emissions.cpp.s

Emissions_Dec.o: Emissions_Dec.cpp.o

.PHONY : Emissions_Dec.o

# target to build an object file
Emissions_Dec.cpp.o:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Emissions_Dec.cpp.o
.PHONY : Emissions_Dec.cpp.o

Emissions_Dec.i: Emissions_Dec.cpp.i

.PHONY : Emissions_Dec.i

# target to preprocess a source file
Emissions_Dec.cpp.i:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Emissions_Dec.cpp.i
.PHONY : Emissions_Dec.cpp.i

Emissions_Dec.s: Emissions_Dec.cpp.s

.PHONY : Emissions_Dec.s

# target to generate assembly for a file
Emissions_Dec.cpp.s:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Emissions_Dec.cpp.s
.PHONY : Emissions_Dec.cpp.s

Forward_Backward.o: Forward_Backward.cpp.o

.PHONY : Forward_Backward.o

# target to build an object file
Forward_Backward.cpp.o:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Forward_Backward.cpp.o
.PHONY : Forward_Backward.cpp.o

Forward_Backward.i: Forward_Backward.cpp.i

.PHONY : Forward_Backward.i

# target to preprocess a source file
Forward_Backward.cpp.i:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Forward_Backward.cpp.i
.PHONY : Forward_Backward.cpp.i

Forward_Backward.s: Forward_Backward.cpp.s

.PHONY : Forward_Backward.s

# target to generate assembly for a file
Forward_Backward.cpp.s:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Forward_Backward.cpp.s
.PHONY : Forward_Backward.cpp.s

Model.o: Model.cpp.o

.PHONY : Model.o

# target to build an object file
Model.cpp.o:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Model.cpp.o
.PHONY : Model.cpp.o

Model.i: Model.cpp.i

.PHONY : Model.i

# target to preprocess a source file
Model.cpp.i:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Model.cpp.i
.PHONY : Model.cpp.i

Model.s: Model.cpp.s

.PHONY : Model.s

# target to generate assembly for a file
Model.cpp.s:
	$(MAKE) -f CMakeFiles/dihmm.dir/build.make CMakeFiles/dihmm.dir/Model.cpp.s
.PHONY : Model.cpp.s

dihmm_ext.o: dihmm_ext.cpp.o

.PHONY : dihmm_ext.o

# target to build an object file
dihmm_ext.cpp.o:
	$(MAKE) -f CMakeFiles/dihmm_ext.dir/build.make CMakeFiles/dihmm_ext.dir/dihmm_ext.cpp.o
.PHONY : dihmm_ext.cpp.o

dihmm_ext.i: dihmm_ext.cpp.i

.PHONY : dihmm_ext.i

# target to preprocess a source file
dihmm_ext.cpp.i:
	$(MAKE) -f CMakeFiles/dihmm_ext.dir/build.make CMakeFiles/dihmm_ext.dir/dihmm_ext.cpp.i
.PHONY : dihmm_ext.cpp.i

dihmm_ext.s: dihmm_ext.cpp.s

.PHONY : dihmm_ext.s

# target to generate assembly for a file
dihmm_ext.cpp.s:
	$(MAKE) -f CMakeFiles/dihmm_ext.dir/build.make CMakeFiles/dihmm_ext.dir/dihmm_ext.cpp.s
.PHONY : dihmm_ext.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... dihmm_ext"
	@echo "... rebuild_cache"
	@echo "... dihmm"
	@echo "... Emissions.o"
	@echo "... Emissions.i"
	@echo "... Emissions.s"
	@echo "... Emissions_Dec.o"
	@echo "... Emissions_Dec.i"
	@echo "... Emissions_Dec.s"
	@echo "... Forward_Backward.o"
	@echo "... Forward_Backward.i"
	@echo "... Forward_Backward.s"
	@echo "... Model.o"
	@echo "... Model.i"
	@echo "... Model.s"
	@echo "... dihmm_ext.o"
	@echo "... dihmm_ext.i"
	@echo "... dihmm_ext.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

