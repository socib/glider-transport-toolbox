#!/bin/sh

####################################################################################
#                                                                                  #
# This script computes the glider geostrophic transports over a period             #
#  (processing all past missions from a given date)                                #
#                                                                                  #
# For every mission, the script:                                                   #
#  - computes the geostrophic transports for each transect                         #
#  - generates figures of transports                                               #
#                                                                                  #
# Inputs:                                                                          #
#  - L1 glider file (SOCIB thredds)                                                #
#                                                                                  #
# Outputs:                                                                         #
#  - matfiles/figures from transect/transport processing                           #
#  - figures transports over the period                                            #
#                                                                                  #
# Author: Melanie Juza (mjuza@socib.es)                                            #
#                                                                                  #
# Creation date: February 2018 (last modification: June 2022)                      #
#                                                                                  #
####################################################################################

### Parameters
# datatype   : L1
# datamode   : dt
# dataset    : G (glider)
# deployment : "Canales"
# strdate_ini: initital date (yyyymmdd) 
# strdate_end: end date (yyyymmdd)
# section    : IbizaChannel, MallorcaChannel
# smooth     : smooth parameter for transport computation (in km)
# critWIW    : criterion for WIW detection = fixed-range, geometry
# logo       : socib, nologo
# optionplot : 1 (true) or 0 (false) to save plot
# optionfile : 1 (true) or 0 (false) to save file
# optionbgc  : 1 (true) or 0 (false) to process biogeochemical data

####################################################################################

datatype=L1
datamode=dt
dataset=G;
deployment="Canales"
strdate_ini=$(date -d "-3 month" +"%Y%m%d")
strdate_end=$(date +%Y%m%d) 
section=IbizaChannel
smooth=24                   
critWIW=geometry
logo=socib
optionplot=1
optionfile=1
optionbgc=0

### TO BE COMPLETED BY THE CLIENT
main_dir=
glider_dir=${main_dir}/Transports_smooth${smooth}_critWIW_${critWIW}_onlyglider
image_dir=${glider_dir}/timeseries/
tools_dir=${HOME}/glidertoolbox_JS3_D2PTS
log_dir=${tools_dir}/logs

####################################################################################

### Directories
mkdir -p $glider_dir
mkdir -p $image_dir
mkdir -p $log_dir

####################################################################################

### Sorted list (by time) of "good" glider mission files
echo 'Listing the glider mission files from SOCIB thredds...'
filename_glider_mission_list=${glider_dir}/SOCIB_glider_mission_file_list.txt
${tools_dir}/launch/generate_SOCIB_glider_mission_file_list.sh ${datamode} ${datatype} ${deployment} ${strdate_ini} ${strdate_end} ${filename_glider_mission_list}

### Compute glider transports for every mission
echo 'Computing transport for every glider file if not already processed'
while read filename; do

  # Glider mission
  echo $filename 
  # Compute glider transport 
  ${tools_dir}/launch/fake_display.sh "-logfile ${log_dir}/out_matlab_compute_transports_${filename}.log -r addpath(genpath('${tools_dir}/matlab'));compute_transports_glider('${filename}','${datamode}','${dataset}','${section}',${smooth},'${critWIW}','${glider_dir}','${tools_dir}',${optionfile},${optionplot},${optionbgc},'${logo}','off');exit;"

  echo $filename
done < ${filename_glider_mission_list}

####################################################################################

### List of processed mission by dates
echo 'Listing the processed glider files' 
cd ${glider_dir}/matfiles_glider_${datamode}/
missionList=''
dates=`ls -d canales* | cut -c 8-14`

for mydate in $dates; do
  yyyy=`echo $mydate | cut -c 4-7`
  month=`echo $mydate | cut -c 1-3`
  mm=`mmm=$month;  LC_ALL=C date -d "1 $mmm $yyyy" +%m`
  echo $yyyy-$mm >> listdates.txt
done 

sort listdates.txt >> sortedlistdates.txt

while read filename; do
  yyyy=`echo $filename | cut -c 1-4`
  month=`echo $filename | cut -c 6-7`
  mmm=`mm=$month;  LC_ALL=C date -d ${yyyy}/${mm}/01 +%b`
  missionList="$missionList 'canales${mmm}${yyyy}'"
done < sortedlistdates.txt
export missionList=`echo $missionList | sed -e "s/ /;/g"`
\rm listdates.txt sortedlistdates.txt

### Plot time series
echo 'Plotting the transport time series'
${tools_dir}/launch/fake_display.sh "-r addpath(genpath('${tools_dir}/matlab'));plot_glideronly_transports({${missionList}},'${datamode}','${section}','${logo}','${glider_dir}','${image_dir}');exit;"
${tools_dir}/launch/fake_display.sh "-r addpath(genpath('${tools_dir}/matlab'));plot_glideronly_transports_bar({${missionList}},'${datamode}','${section}','${logo}','${glider_dir}','${image_dir}');exit;"

