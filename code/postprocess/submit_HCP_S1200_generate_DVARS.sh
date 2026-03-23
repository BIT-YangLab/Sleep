#!/bin/sh
sub_list=$1
proj_dir=$2
data_dir=$3
code_path=`pwd`
cnt=0
for data_num in `cat ${sub_list}`
do 
{
  if [ $cnt = 20 ] ;then
    csh ${code_path}/HCP_S1200_generate_DVARS.csh ${data_num} ${proj_dir} ${data_dir} 
    cnt=0
  else
    csh ${code_path}/HCP_S1200_generate_DVARS.csh ${data_num} ${proj_dir} ${data_dir} &
  fi
  cnt=$((cnt+1))
}
done
