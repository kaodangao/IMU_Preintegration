# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rixap/projects/preintegration

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rixap/projects/preintegration/build

# Include any dependencies generated for this target.
include CMakeFiles/preintegration.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/preintegration.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/preintegration.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/preintegration.dir/flags.make

CMakeFiles/preintegration.dir/test_preintegration.cpp.o: CMakeFiles/preintegration.dir/flags.make
CMakeFiles/preintegration.dir/test_preintegration.cpp.o: ../test_preintegration.cpp
CMakeFiles/preintegration.dir/test_preintegration.cpp.o: CMakeFiles/preintegration.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rixap/projects/preintegration/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/preintegration.dir/test_preintegration.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/preintegration.dir/test_preintegration.cpp.o -MF CMakeFiles/preintegration.dir/test_preintegration.cpp.o.d -o CMakeFiles/preintegration.dir/test_preintegration.cpp.o -c /home/rixap/projects/preintegration/test_preintegration.cpp

CMakeFiles/preintegration.dir/test_preintegration.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/preintegration.dir/test_preintegration.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rixap/projects/preintegration/test_preintegration.cpp > CMakeFiles/preintegration.dir/test_preintegration.cpp.i

CMakeFiles/preintegration.dir/test_preintegration.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/preintegration.dir/test_preintegration.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rixap/projects/preintegration/test_preintegration.cpp -o CMakeFiles/preintegration.dir/test_preintegration.cpp.s

# Object files for target preintegration
preintegration_OBJECTS = \
"CMakeFiles/preintegration.dir/test_preintegration.cpp.o"

# External object files for target preintegration
preintegration_EXTERNAL_OBJECTS =

preintegration: CMakeFiles/preintegration.dir/test_preintegration.cpp.o
preintegration: CMakeFiles/preintegration.dir/build.make
preintegration: /usr/local/lib/libglog.so.0.6.0
preintegration: /usr/local/lib/libgflags.so.2.2.2
preintegration: CMakeFiles/preintegration.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rixap/projects/preintegration/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable preintegration"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/preintegration.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/preintegration.dir/build: preintegration
.PHONY : CMakeFiles/preintegration.dir/build

CMakeFiles/preintegration.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/preintegration.dir/cmake_clean.cmake
.PHONY : CMakeFiles/preintegration.dir/clean

CMakeFiles/preintegration.dir/depend:
	cd /home/rixap/projects/preintegration/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rixap/projects/preintegration /home/rixap/projects/preintegration /home/rixap/projects/preintegration/build /home/rixap/projects/preintegration/build /home/rixap/projects/preintegration/build/CMakeFiles/preintegration.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/preintegration.dir/depend

