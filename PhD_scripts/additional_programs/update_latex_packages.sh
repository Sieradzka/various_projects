#first you have to dowloand index.html/ If the commad belov diesn't work try by wget -r
#if you are in different country check the website: http://ctan.org/mirrors#Europe to change the address of website
wget http://mirror.ox.ac.uk/sites/ctan.org/install/macros/latex/contrib/

#greb name of packages to a file:
sed 's/.*href="\([^"]*\).*/\1/' index.html | cat >name_packages

#dowload package from website:
#while read p; do
#  echo $p
#  wget http://mirror.ox.ac.uk/sites/ctan.org/install/macros/latex/contrib/$p
#done <name_packages

#unzip folders:
#while read p; do
#  echo $p
#  unzip  $p
#done <name_packages

