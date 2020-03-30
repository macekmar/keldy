file(GLOB_RECURSE sources *.cpp *.hpp)
add_library(keldy_c ${sources})

# Link against dependencies
target_link_libraries(keldy_c PUBLIC triqs triqs::boost triqs::blas_lapack GSL::gsl PRIVATE project_warnings triqs_lapack_additions)

# Configure compilation
target_compile_options(keldy_c PUBLIC -fPIC)
target_include_directories(keldy_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
target_compile_definitions(keldy_c PUBLIC
				KELDY_GIT_HASH=${KELDY_GIT_HASH}
				TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
				$<$<CONFIG:Debug>:KELDY_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
			  )

# Install library and headers
install(TARGETS keldy_c EXPORT keldy-targets DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")


# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(keldy_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
	TARGET keldy_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++17
	  --template=gcc
	  --verbose
	  --quiet
          ${sources}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()


# ========= Dynamic Analyzer Checks ==========

option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")

if(ASAN)
  find_package(sanitizer REQUIRED "asan")
  target_link_libraries(keldy_c PUBLIC asan)
  install(TARGETS asan EXPORT keldy-targets)
endif()
if(UBSAN)
  find_package(sanitizer REQUIRED "ubsan")
  target_link_libraries(keldy_c PUBLIC ubsan)
  install(TARGETS ubsan EXPORT keldy-targets)
endif()
