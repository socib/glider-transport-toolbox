
###############################################################################
#                                                                             #
#  Generate glider data file list from catalog files read on SOCIB thredds    #
#                                                                             #
#   Arguments:                                                                #
#     - data mode: real time (rt) / delayed time (dt)                         #
#     - data level: L1/L2                                                     #
#     - output: text file containing the list of netcdf glider data files     #
#                                                                             #
#  Author: Baptiste Mourre (bmourre@socib.es)                                 #
#  Creation: 08-Feb-2016                                                      #
#                                                                             #
#  Latest modification: 30-May-2017 (Melanie Juza)                            #
#                       11-Jul-2017 (Melanie Juza) - deployment specificity   #
#                                                    "good" files only        #
#                       10-Jul-2018 (Melanie Juza) - refine mission date      #
#                                                                             #
###############################################################################

# ===== Glider data mode (rt=Real Time; dt=Delayed Time)   =====
datamode=$1

# ===== Glider data type (L1, L2)   =====
datatype=$2

# ===== Deployement name (Canales,all)
deployment=$3

# ===== First and last date (YYYYMMDD)
strdate_ini=$4
strdate_end=$5

# ===== Output glider data file list   =====
outfile=$6

# ===== Glider directory on thredds  =====
thredds_dir=http://thredds.socib.es/thredds/catalog/auv/glider
thredds=http://thredds.socib.es/thredds/dodsC/auv/glider

# ===== Remove former glider mission file list if present   =====
\rm ${outfile}

# ===== Get upper level catalog.html to determine the glider names  =====
wget ${thredds_dir}/catalog.html

mv catalog.html catalog_allgliders_list.html

sed -n '/Folder/p' catalog_allgliders_list.html > tmp1.txt
nbgliders=`wc -l < tmp1.txt`


# ===== Loop over the glider names  =====
kg=2   # begin line 2 since first line is the general glider directory

while test ${kg} -le ${nbgliders}
  do

   line=$(head -n ${kg} tmp1.txt | tail -1)
   line=${line%%/catalog*}
   glider_name=${line##*href=\'}

   # ===== Loop over the years  =====
   first_year=`echo $strdate_ini | cut -c 1-4 `
   last_year=`echo $strdate_end | cut -c 1-4 `

   myyear=${first_year}

   while test ${myyear} -le ${last_year}
   do

    wget ${thredds_dir}/${glider_name}/${datatype}/${myyear}/catalog.html

    mv catalog.html catalog_${glider_name}_${myyear}.html

    ext=_${datamode}.nc
    sed -n 's/'"$ext"'/&/p' catalog_${glider_name}_${myyear}.html > tmp2.txt

    nbfiles=`wc -l < tmp2.txt`

    # ===== Loop over the number of files  =====
    kf=1
    while test ${kf} -le ${nbfiles}
    do

     line=$(head -n ${kf} tmp2.txt | tail -1)
     line=${line%%_${datamode}.nc*}
     line=${line##*/dep}
     gliderfile=dep${line}_${datamode}.nc

     mymission="$deployment $myyear"
     touch filename.txt
     ncdump -h ${thredds}/${glider_name}/${datatype}/${myyear}/${gliderfile} | grep "${mymission}" >> filename.txt
     tmp=`echo ${gliderfile} | sed -e "s/${datatype}_/ /g" -e "s/_data/ /g"`; set $tmp; 
     mydate=$2; mydate=$(date -d "$mydate" +%Y%m%d)
     if [ -s filename.txt ] || [ $deployment = "all" ]; then
      touch filename2.txt
       ncdump -h ${thredds}/${glider_name}/${datatype}/${myyear}/${gliderfile} | grep -e 'butterfly' -e 'configuration error' -e 'aborted because' -e 'leak' -e 'Premature recovery' -e 'Cancelled' >> filename2.txt
       if [ ! -s filename2.txt ] && [ ${mydate} -ge ${strdate_ini} ] && [ ${mydate} -lt ${strdate_end} ] ; then
        echo ${thredds}/${glider_name}/${datatype}/${myyear}/${gliderfile} >> out.txt
       fi
       \rm filename2.txt
     fi
     \rm filename.txt 

     kf=`expr ${kf} + 1`

   done

   \rm catalog_${glider_name}_${myyear}.html
   \rm tmp2.txt

   myyear=`expr ${myyear} + 1`

  done

  kg=`expr ${kg} + 1`

done

cat out.txt | sort -t'_' -k5 >> ${outfile}

\rm catalog_allgliders_list.html
\rm tmp1.txt out.txt
                                                                                                                                                                                                   
