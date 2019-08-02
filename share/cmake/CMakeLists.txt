configure_file(keldy-config.cmake.in keldy-config.cmake @ONLY)
configure_file(keldy-config-version.cmake.in keldy-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/keldy-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/keldy-config-version.cmake
  DESTINATION lib/cmake/keldy
)

install(EXPORT keldy-targets DESTINATION lib/cmake/keldy)
