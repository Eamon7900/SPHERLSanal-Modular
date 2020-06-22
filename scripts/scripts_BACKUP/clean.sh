#!/bin/bash
#Removes all files that do not end in .py, .sh or .txt
mkdir temp
mv *.py temp
mv *.sh temp
mv *.txt temp
rm -f *
mv temp/* .
rm -d temp