mkdir temp
mv *.py temp
mv *.sh temp
rm -f *
mv temp/* .
rm -d temp