代码文件共有 `dsu.hpp`，`graph_generator.hpp`，`naive_pso.hpp`，`hybrid_pso.hpp`，`main.cpp` 共五个，其中朴素粒子群算法的代码在文件 `naive_pso.hpp` 中，混合算法的代码在 `hybrid_pso.hpp` 中。`main.cpp` 用于对两种算法进行测试。

如需运行代码，编译运行 `main.cpp` 即可。

我的编译参数是 `g++ -std=c++20 main.cpp -o main`，但如果编译器不支持 C++20 也无伤大雅，理论上编译器只需支持 C++11 的特性即可完成编译；

![Img](r1.png)

![Img](r11.png)

我的编译器版本是 `TDM-GCC 10.2.0`，操作系统为 Windows10：

![Img](r2.png)

```
Using built-in specs.
COLLECT_GCC=g++
COLLECT_LTO_WRAPPER=C:/mingw/MinGW/bin/../lib/gcc/i686-w64-mingw32/10.2.0/lto-wrapper.exe
Target: i686-w64-mingw32
Configured with: ../gcc-10.2.0/configure --prefix=/mingw32 --with-local-prefix=/mingw32/local --build=i686-w64-mingw32 --host=i686-w64-mingw32 --target=i686-w64-mingw32 --with-native-system-header-dir=/mingw32/i686-w64-mingw32/include --libexecdir=/mingw32/lib --enable-bootstrap --with-arch=i686 --with-tune=generic --enable-languages=c,lto,c++,fortran,ada,objc,obj-c++,jit --enable-shared --enable-static --enable-libatomic --enable-threads=posix --enable-graphite --enable-fully-dynamic-string --enable-libstdcxx-filesystem-ts=yes --enable-libstdcxx-time=yes --disable-libstdcxx-pch --disable-libstdcxx-debug --disable-isl-version-check --enable-lto --enable-libgomp --disable-multilib --enable-checking=release --disable-rpath --disable-win32-registry --disable-nls --disable-werror --disable-symvers --disable-plugin --with-libiconv --with-system-zlib --with-gmp=/mingw32 --with-mpfr=/mingw32 --with-mpc=/mingw32 --with-isl=/mingw32 --with-pkgversion='Rev6, Built by MSYS2 project' --with-bugurl=https://github.com/msys2/MINGW-packages/issues --with-gnu-as --with-gnu-ld --with-boot-ldflags='-pipe -Wl,--dynamicbase,--nxcompat,--no-seh -Wl,--large-address-aware -Wl,--disable-dynamicbase -static-libstdc++ -static-libgcc' 'LDFLAGS_FOR_TARGET=-pipe -Wl,--dynamicbase,--nxcompat,--no-seh -Wl,--large-address-aware' --enable-linker-plugin-flags='LDFLAGS=-static-libstdc++\ -static-libgcc\ -pipe\ -Wl,--dynamicbase,--nxcompat,--no-seh\ -Wl,--large-address-aware\ -Wl,--stack,12582912' --disable-sjlj-exceptions --with-dwarf2
Thread model: posix
Supported LTO compression algorithms: zlib zstd
gcc version 10.2.0 (Rev6, Built by MSYS2 project)

```
