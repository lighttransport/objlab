rem If you use VS2019, use Visual Stuio 2019's cmake feature(`clone repo` or `open folder`)

rmdir /q /s build
mkdir build

cmake.exe -G "Visual Studio 15 2017" -A x64  -Bbuild .
