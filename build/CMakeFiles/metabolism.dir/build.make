# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /qfs/projects/ops/rh7_gpu/rocmapps/views/cmake/3.21.4/bin/cmake

# The command to remove a file.
RM = /qfs/projects/ops/rh7_gpu/rocmapps/views/cmake/3.21.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /people/cann484/metab

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /people/cann484/metab/build

# Include any dependencies generated for this target.
include CMakeFiles/metabolism.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/metabolism.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/metabolism.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/metabolism.dir/flags.make

CMakeFiles/metabolism.dir/src/main.cc.o: CMakeFiles/metabolism.dir/flags.make
CMakeFiles/metabolism.dir/src/main.cc.o: ../src/main.cc
CMakeFiles/metabolism.dir/src/main.cc.o: CMakeFiles/metabolism.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/people/cann484/metab/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/metabolism.dir/src/main.cc.o"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/metabolism.dir/src/main.cc.o -MF CMakeFiles/metabolism.dir/src/main.cc.o.d -o CMakeFiles/metabolism.dir/src/main.cc.o -c /people/cann484/metab/src/main.cc

CMakeFiles/metabolism.dir/src/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/metabolism.dir/src/main.cc.i"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /people/cann484/metab/src/main.cc > CMakeFiles/metabolism.dir/src/main.cc.i

CMakeFiles/metabolism.dir/src/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/metabolism.dir/src/main.cc.s"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /people/cann484/metab/src/main.cc -o CMakeFiles/metabolism.dir/src/main.cc.s

CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.o: CMakeFiles/metabolism.dir/flags.make
CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.o: ../src/ifopt_constraint_classes.cc
CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.o: CMakeFiles/metabolism.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/people/cann484/metab/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.o"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.o -MF CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.o.d -o CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.o -c /people/cann484/metab/src/ifopt_constraint_classes.cc

CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.i"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /people/cann484/metab/src/ifopt_constraint_classes.cc > CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.i

CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.s"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /people/cann484/metab/src/ifopt_constraint_classes.cc -o CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.s

CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.o: CMakeFiles/metabolism.dir/flags.make
CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.o: ../src/ifopt_cost_class.cc
CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.o: CMakeFiles/metabolism.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/people/cann484/metab/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.o"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.o -MF CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.o.d -o CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.o -c /people/cann484/metab/src/ifopt_cost_class.cc

CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.i"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /people/cann484/metab/src/ifopt_cost_class.cc > CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.i

CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.s"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /people/cann484/metab/src/ifopt_cost_class.cc -o CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.s

CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.o: CMakeFiles/metabolism.dir/flags.make
CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.o: ../src/ifopt_variable_class.cc
CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.o: CMakeFiles/metabolism.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/people/cann484/metab/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.o"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.o -MF CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.o.d -o CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.o -c /people/cann484/metab/src/ifopt_variable_class.cc

CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.i"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /people/cann484/metab/src/ifopt_variable_class.cc > CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.i

CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.s"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /people/cann484/metab/src/ifopt_variable_class.cc -o CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.s

CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.o: CMakeFiles/metabolism.dir/flags.make
CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.o: ../src/maximum_entropy_relaxed.cc
CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.o: CMakeFiles/metabolism.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/people/cann484/metab/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.o"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.o -MF CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.o.d -o CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.o -c /people/cann484/metab/src/maximum_entropy_relaxed.cc

CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.i"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /people/cann484/metab/src/maximum_entropy_relaxed.cc > CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.i

CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.s"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /people/cann484/metab/src/maximum_entropy_relaxed.cc -o CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.s

CMakeFiles/metabolism.dir/src/read_data_files.cc.o: CMakeFiles/metabolism.dir/flags.make
CMakeFiles/metabolism.dir/src/read_data_files.cc.o: ../src/read_data_files.cc
CMakeFiles/metabolism.dir/src/read_data_files.cc.o: CMakeFiles/metabolism.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/people/cann484/metab/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/metabolism.dir/src/read_data_files.cc.o"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/metabolism.dir/src/read_data_files.cc.o -MF CMakeFiles/metabolism.dir/src/read_data_files.cc.o.d -o CMakeFiles/metabolism.dir/src/read_data_files.cc.o -c /people/cann484/metab/src/read_data_files.cc

CMakeFiles/metabolism.dir/src/read_data_files.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/metabolism.dir/src/read_data_files.cc.i"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /people/cann484/metab/src/read_data_files.cc > CMakeFiles/metabolism.dir/src/read_data_files.cc.i

CMakeFiles/metabolism.dir/src/read_data_files.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/metabolism.dir/src/read_data_files.cc.s"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /people/cann484/metab/src/read_data_files.cc -o CMakeFiles/metabolism.dir/src/read_data_files.cc.s

CMakeFiles/metabolism.dir/src/metabolism.cc.o: CMakeFiles/metabolism.dir/flags.make
CMakeFiles/metabolism.dir/src/metabolism.cc.o: ../src/metabolism.cc
CMakeFiles/metabolism.dir/src/metabolism.cc.o: CMakeFiles/metabolism.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/people/cann484/metab/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/metabolism.dir/src/metabolism.cc.o"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/metabolism.dir/src/metabolism.cc.o -MF CMakeFiles/metabolism.dir/src/metabolism.cc.o.d -o CMakeFiles/metabolism.dir/src/metabolism.cc.o -c /people/cann484/metab/src/metabolism.cc

CMakeFiles/metabolism.dir/src/metabolism.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/metabolism.dir/src/metabolism.cc.i"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /people/cann484/metab/src/metabolism.cc > CMakeFiles/metabolism.dir/src/metabolism.cc.i

CMakeFiles/metabolism.dir/src/metabolism.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/metabolism.dir/src/metabolism.cc.s"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /people/cann484/metab/src/metabolism.cc -o CMakeFiles/metabolism.dir/src/metabolism.cc.s

CMakeFiles/metabolism.dir/src/find_initial_values.cc.o: CMakeFiles/metabolism.dir/flags.make
CMakeFiles/metabolism.dir/src/find_initial_values.cc.o: ../src/find_initial_values.cc
CMakeFiles/metabolism.dir/src/find_initial_values.cc.o: CMakeFiles/metabolism.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/people/cann484/metab/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/metabolism.dir/src/find_initial_values.cc.o"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/metabolism.dir/src/find_initial_values.cc.o -MF CMakeFiles/metabolism.dir/src/find_initial_values.cc.o.d -o CMakeFiles/metabolism.dir/src/find_initial_values.cc.o -c /people/cann484/metab/src/find_initial_values.cc

CMakeFiles/metabolism.dir/src/find_initial_values.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/metabolism.dir/src/find_initial_values.cc.i"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /people/cann484/metab/src/find_initial_values.cc > CMakeFiles/metabolism.dir/src/find_initial_values.cc.i

CMakeFiles/metabolism.dir/src/find_initial_values.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/metabolism.dir/src/find_initial_values.cc.s"
	/share/apps/gcc/9.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /people/cann484/metab/src/find_initial_values.cc -o CMakeFiles/metabolism.dir/src/find_initial_values.cc.s

# Object files for target metabolism
metabolism_OBJECTS = \
"CMakeFiles/metabolism.dir/src/main.cc.o" \
"CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.o" \
"CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.o" \
"CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.o" \
"CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.o" \
"CMakeFiles/metabolism.dir/src/read_data_files.cc.o" \
"CMakeFiles/metabolism.dir/src/metabolism.cc.o" \
"CMakeFiles/metabolism.dir/src/find_initial_values.cc.o"

# External object files for target metabolism
metabolism_EXTERNAL_OBJECTS =

metabolism: CMakeFiles/metabolism.dir/src/main.cc.o
metabolism: CMakeFiles/metabolism.dir/src/ifopt_constraint_classes.cc.o
metabolism: CMakeFiles/metabolism.dir/src/ifopt_cost_class.cc.o
metabolism: CMakeFiles/metabolism.dir/src/ifopt_variable_class.cc.o
metabolism: CMakeFiles/metabolism.dir/src/maximum_entropy_relaxed.cc.o
metabolism: CMakeFiles/metabolism.dir/src/read_data_files.cc.o
metabolism: CMakeFiles/metabolism.dir/src/metabolism.cc.o
metabolism: CMakeFiles/metabolism.dir/src/find_initial_values.cc.o
metabolism: CMakeFiles/metabolism.dir/build.make
metabolism: ../PREREQS_INSTALL/lib64/libifopt_ipopt.so
metabolism: ../PREREQS_INSTALL/lib64/libifopt_core.so
metabolism: CMakeFiles/metabolism.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/people/cann484/metab/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable metabolism"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/metabolism.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/metabolism.dir/build: metabolism
.PHONY : CMakeFiles/metabolism.dir/build

CMakeFiles/metabolism.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/metabolism.dir/cmake_clean.cmake
.PHONY : CMakeFiles/metabolism.dir/clean

CMakeFiles/metabolism.dir/depend:
	cd /people/cann484/metab/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /people/cann484/metab /people/cann484/metab /people/cann484/metab/build /people/cann484/metab/build /people/cann484/metab/build/CMakeFiles/metabolism.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/metabolism.dir/depend

