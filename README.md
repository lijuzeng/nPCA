# nPCA
neural principal component annlysis
## Compile
- linux: g++
> g++ -o npca -O3 -Wall -march=native main.cpp npca.hpp wyhash.h
## Usage
### Optional parameters
> chmod +x npca
> 
> ./npca

Usage: npca [options] data
- -e: learning rate = 0.001
- -g: max of gr = 4
## Example
> ./npca Frogs_MFCCs.data
