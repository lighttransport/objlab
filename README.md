# objlab, Simple Wavefront .obj viewer + some featues

objlab is a simple Wavefront .obj viewer + some features, written in C++11.

## Requirements

* cmake
* OpenGL 2.x

## Build

### Setup

```
$ git submodule update --init --recursive
```

### CMake build

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```


## License

objlab is licensed under MIT license.


### Third party license

* TinyObjLoader : MIT license.
* imgui : MIT license.
* stb_image, stb_image_write, stb_image_resize : Public domain
* glfw3 : zlib/libpng license.
