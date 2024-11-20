This repo contains tools to help visualize the problem of finding Mutually Unbiased Bases (MUBs).

Here there are two versions of most or less the same thing: a C++ version and a JavaScript version. The C++ version is more complete.

### Compiling the C++ version

You need the library Optim installed into the folder $HOME/Libraries/optim/header_only_version. Or use a different folder and change the makefile. You also need SFML and Eigen:
```
sudo apt-get install libsfml-dev libeigen3-dev
```
Then you can use "make" in the C++ folder to compile the executable. 
