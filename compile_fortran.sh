#!/bin/sh

yn="Yes"
while [ $yn = "Yes" ]; do
  echo "Add a module?"
  select yn in "Yes" "No"
  do
    case $yn in
    Yes)
      read -p "Enter module name: " mymod
      gfortran -c $mymod.f90
      echo "Add a module?"
      echo "1) Yes"
      echo "2) No"
      ;;
    No)
      break
      ;;
  esac
  done
done
read -p "Enter name of main program: " mymain; gfortran -c $mymain.f90
myexe="execute"
gfortran *.o -o $myexe
./$myexe
rm *.o
rm *.mod
rm $myexe

