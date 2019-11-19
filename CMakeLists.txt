cmake_minimum_required(VERSION 3.5)

set(BUILD_TARGET "objlab")

project(${BUILD_TARGET})

find_package(OpenGL REQUIRED)

set(CMAKE_CXX_STANDARD 11)

# [glfw3]
set(GLFW_BUILD_EXAMPLES OFF CACHE INTERNAL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE INTERNAL "" FORCE)
set(GLFW_BUILD_DOCS OFF CACHE INTERNAL "" FORCE)
set(GLFW_INSTALL OFF CACHE INTERNAL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third_party/glfw)

set(OBJLAB_SOURCES
    src/app.cc
    src/main.cc
    )

set(OBJLAB_DEP_SOURCES
    #src/glad.c
    src/trackball.cc
    deps/glad/src/glad.c
    deps/imgui/imgui.cpp
    deps/imgui/imgui_draw.cpp
    deps/imgui/imgui_widgets.cpp
    deps/imgui/imgui_impl_opengl3.cpp
    deps/imgui/imgui_impl_glfw.cpp
    deps/ImGuizmo/ImGuizmo.cpp
    )

add_executable(${BUILD_TARGET} ${OBJLAB_SOURCES} ${OBJLAB_DEP_SOURCES})

# [filesytem]
target_include_directories(${BUILD_TARGET} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party/filesystem/include")

# [imgui]
target_include_directories(${BUILD_TARGET} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/deps/")

# [glad]
add_definitions("-DIMGUI_IMPL_OPENGL_LOADER_GLAD")
target_include_directories(${BUILD_TARGET} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/deps/glad/include")

target_include_directories(${BUILD_TARGET} PRIVATE ${glfw_INCLUDE_DIRS})
target_link_libraries(${BUILD_TARGET} glfw ${OPENGL_LIBRARIES} ${CMAKE_DL_LIBS})

# Increase warning level for clang.
IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set_source_files_properties(${OBJLAB_SOURCES} PROPERTIES COMPILE_FLAGS "-Weverything -Wno-system-headers -Werror -Wno-padded -Wno-c++98-compat-pedantic -Wno-documentation -Wno-documentation-unknown-command -Wno-reserved-id-macro -Wno-unused-variable -Wno-unused-function")
ENDIF ()


# [VisualStudio]
if(WIN32)
  # Set ${BUILD_TARGET} as a startup project for VS IDE
  set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT ${BUILD_TARGET})

  # For easier debugging in VS IDE(cmake 3.8.0 or later required) Set working
  # directory to ${BUILD_TARGET} git repo root.
  if(CMAKE_VERSION VERSION_GREATER 3.8.0)
    set_target_properties(${BUILD_TARGET}
                          PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY
                                     "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
endif()
