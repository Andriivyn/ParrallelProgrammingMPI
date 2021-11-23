#!/bin/bash

usage() {
  echo "Usage:       ${0##*/} [-np N] class [-nobuild]"
  echo "For example: ${0##*/} -np 9 com.mathpar.students.ukma.atsaruk.mpi.HelloWorld"
  echo "Options:"
  echo "  -np N      sets number of processors in MPI World to N. Defaults to 16"
  echo "  -nobuild   runs class without build"
}

exit_on_error() {
  exit_code=$1
  if [ $exit_code -ne 0 ]; then
    echo "Exiting with code $exit_code"
    exit $exit_code
  fi
}

if ([ "$#" -eq 3 ] || ([ $# -eq 4 ] && [ "$4" = "-nobuild" ])) && ([ "$1" != "-np" ] || ! [[ "$2" =~ ^[0-9]+$ ]]); then
  usage
  exit 1
fi

if [ ! -x "$(command -v java)" ]; then
  echo "Java is not found. Please visit https://www.java.com and install it"
  exit 1
fi

if [ ! -x "$(command -v mvn)" ]; then
  echo "Maven is not found. Please visit https://maven.apache.org and install it"
  exit 1
fi

if [ ! -x "$(command -v mpirun)" ]; then
  echo "Open MPI is not found. Please visit https://www.open-mpi.org and install it"
  exit 1
fi

if [ ! -e hostfile ]; then
  echo "hostfile not found. Trying to create default hostfile"
  echo "localhost slots=16" > hostfile
  exit_on_error $?
  echo "hostfile created"
fi

if [ "$#" -eq 1 ] || ([ "$#" -eq 2 ] && [ "$2" = "-nobuild" ]); then
  np=16
  classname=$1
else
  np=$2
  classname=$3
fi


if ([ "$4" != "-nobuild" ]) && ([ "$2" != "-nobuild" ]); then
  mvn compile
  shift 3
else
  shift 4
fi

target=$PWD/target
classes=$target/classes
lib=$target/lib

classpath=$classes

for entry in "$lib"/*.jar
do
  #Skipping mpi jar file as one is provided by mpirun
  if ! [[ ${entry##*/} =~ .*mpi.*  ]]; then
    classpath=$classpath:$entry
  fi
done

mpirun  -X --hostfile hostfile -np $np java -Xmx4g -cp "$classpath" $classname "$@"

