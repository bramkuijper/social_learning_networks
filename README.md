# Social learning using networks in kin-structured populations

## How to run the simulation
Simulations have been written in C++, which means we first need to compile the source code into an executable file and then we can run the executable file to gather results.

### Compiling the simulation code
To compile the simulation code we need to have the building tool `cmake` installed. On platforms like `msys2` this can be installed via the command `pacman -S mingw-w64-x86_64-cmake`. 

Once `cmake` is installed, navigate to the `src/ibm` directory of your local repository and then run
```
cmake - S . -B build/
```
which merely prepares the compilation cycle by telling `cmake` where the source code is (`-S .`, i.e., it is in the current directory (indicated by `.`) and where the executable will be built (`-B biuld/`, i.e., in the `build/` directory.

Next, we then compile the simulation:
```
cmake --build build
```
and the executable can be found in the directory `build`.

## Running the simulation
You can then run the simulation by doing
```
build/NetworksIBM
```


