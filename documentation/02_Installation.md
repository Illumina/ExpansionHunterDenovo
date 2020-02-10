# Installation

ExpansionHunter Denovo is designed for Linux and macOS operating systems.
Compiled binaries can be downloaded from the [releases
page](https://github.com/Illumina/ExpansionHunterDenovo/releases).
Alternatively, the program can be built from source following the
instructions below.

## Building from source

Prerequisites:

- GCC or clang compiler supporting C++11 standard
- CMake version 3.10 or above
- Libraries zlib, bzip2, liblzma along with their development files; these
  can be installed on Ubuntu Linux like so:

  ```bash
  sudo apt install zlib1g-dev libbz2-dev liblzma-dev
  ```
  
- Active internet connection (to automatically download Boost libraries)

Once the above prerequisites are satisfied the program can be built as follows:

```bash
cd ExpansionHunterDenovo/
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ../source
make
```

If the build procedure succeeds, the `build` directory will contain the
`ExpansionHunterDenovo` binary file.
