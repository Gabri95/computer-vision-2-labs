#!/bin/bash


./final ../3dframes/ 0 n 1.2 4 -1 -1 0 0
./final ../3dframes/ 0 n 1.2 6 -1 -1 0 0
./final ../3dframes/ 0 n 1.2 8 -1 -1 0 0
./final ../3dframes/ 0 n 1.2 10 -1 -1 0 0


./final ../3dframes/ 1 n 1.2 10 -1 -1 0 0
./final ../3dframes/ 1 n 1.2 40 -1 -1 0 0
./final ../3dframes/ 1 n 1.2 60 -1 -1 0 0
./final ../3dframes/ 1 n 1.2 100 -1 -1 0 0
./final ../3dframes/ 1 n 1.2 130 -1 -1 0 0


#Trying to smooth the meshes

./final ../3dframes/ 0 n 1.2 8 0.1 -1 0 0
./final ../3dframes/ 0 n 1.2 8 0.5 -1 0 0
./final ../3dframes/ 1 n 1.2 60 0.1 -1 0 0
./final ../3dframes/ 1 n 1.2 60 0.5 -1 0 0


#Filling the holes

./final ../3dframes/ 0 n 1.2 8 -1 10 0 0
./final ../3dframes/ 1 n 1.2 60 -1 10 0 0

# Adding textures

./final ../3dframes/ 0 t 1.2 8 -1 10 0 0
./final ../3dframes/ 1 t 1.2 60 -1 10 0 0




