#!/bin/bash

touch curDoc.txt
touch doc.txt

for file in *.py 
do
    $file -h > curDoc.txt
    cat curDoc.txt doc.txt
done

rm -f curDoc.txt