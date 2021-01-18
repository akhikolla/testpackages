#!/bin/bash
##############################################
# define function to make symbolic link of files from a directory

function make_symbolic_link
{
# loop over header files
for file in $1/*.h
do 
link=./src/$(basename $1)/$(basename $file)
ln -sf ../../$file $link
#echo $file
done

# loop over header files
for file in $1/*.cpp
do 
link=./src/$(basename $1)/$(basename $file)
ln -sf ../../$file $link
#echo $link
done
}

function make_symbolic_link_only_headers
{
# loop over header files
for file in $1/*.h
do 
link=./src/$(basename $1)/$(basename $file)
ln -sf ../../$file $link
#echo $file
done

}
make_symbolic_link ../../../coclust/src/Algorithms
make_symbolic_link ../../../coclust/src/CoClustFacade
make_symbolic_link ../../../coclust/src/InputParameters
make_symbolic_link ../../../coclust/src/Models
make_symbolic_link ../../../coclust/src/Strategy
make_symbolic_link_only_headers ../../../coclust/src/enumerations
make_symbolic_link_only_headers ../../../coclust/src/Initialization
make_symbolic_link_only_headers ../../../coclust/src/typedefs
