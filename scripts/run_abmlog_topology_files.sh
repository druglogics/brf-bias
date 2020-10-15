#!/bin/bash

# run abmlog@1.6.0 on the scale-free topology .sif files

# Put this script on the root of the `abmlog` repo, as well as
# the `scale_free_gamma2` and `scale_free_gamma2_5` directories
# with the topology files (find them at `data/random_topology_files`)

# gamma = 2.5
topology_files=`ls scale_free_gamma2_5`
files_num=`ls scale_free_gamma2_5 | wc -l`
count=0

echo "Scale-free networks with gamma = 2.5"
for file in ${topology_files}; do
  count=$((count + 1))
  echo Using topology file No. $count/$files_num: $file

  start=`date +%s`
  java -cp target/abmlog-1.6.0-jar-with-dependencies.jar eu.druglogics.abmlog.BooleanModelGenerator --file=scale_free_gamma2_5/$file --attractors=fixpoints --parallel --max-dir-size=10000 > /dev/null 2>&1
  runtime=$(($(date +%s)-$start))
  echo Execution Time: "$(($runtime / 60)) minutes"
done

# gamma = 2
topology_files=`ls scale_free_gamma2`
files_num=`ls scale_free_gamma2 | wc -l`
count=0

echo "Scale-free networks with gamma = 2"
for file in ${topology_files}; do
  count=$((count + 1))
  echo Using topology file No. $count/$files_num: $file

  start=`date +%s`
  java -cp target/abmlog-1.6.0-jar-with-dependencies.jar eu.druglogics.abmlog.BooleanModelGenerator --file=scale_free_gamma2/$file --attractors=fixpoints --parallel --max-dir-size=10000 > /dev/null 2>&1
  runtime=$(($(date +%s)-$start))
  echo Execution Time: "$(($runtime / 60)) minutes"
done

