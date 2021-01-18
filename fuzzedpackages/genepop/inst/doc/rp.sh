
sed -i 's/\$\$/\;/g' newFIle.md

sed -i 's/\;.*\;/\\deqn{&}/g' newFIle.md

sed -i 's/\;//g' newFIle.md

sed -i 's/\$.\$/\\eqn{&}/g' newFIle.md

sed -i 's/\$//g' newFIle.md
