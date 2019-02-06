#!/bin/bash


# -c            :   to also run cmake and build the project again
# no parameters :   to just recompile
# other params  :   if there are other params (at least 1), the project is also run, passing these parameters left as command line argument

CC=/usr/bin/g++-5

rm -f ./final

mkdir -p build

cd ./build

re_make=false

PARAMS=()

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -c)
    re_make=true
    shift
    ;;
    *)    # unknown option
    PARAMS+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

set -- "${PARAMS[@]}" # restore positional parameters

if [ "$re_make" = true ] ; then
    rm -r *
fi

if [ -z "$(ls -A .)" ]; then
   cmake -DCMAKE_CXX_COMPILER:PATH=$CC ../
fi


make -j4

cd ..

if [ "$#" -gt 0 ]; then
    ./final ../3dframes $@ #${PARAMS[@]}
fi


