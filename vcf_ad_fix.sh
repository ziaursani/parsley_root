#!/bin/bash

# vcf file to read in on the command line, e.g. vcf_ad_fix.sh foo.vcf
infile=${1}

# loop over every line in the vcf
cat ${infile} | while read -r line; do
  if [[ ! $line =~ ^#+ ]] ; then
    # get the info section
    attrinfo=$(echo ${line} | awk '{print $8}')
    # split on ; into attrs array
    IFS=';' read -ra attrs <<< "${attrinfo}"
    # get the AO field (6th element in array) value
    alt=$(echo ${attrs[5]} | sed 's/AO=//g')

    # get the fields section
    fieldinfo=$(echo ${line} | awk '{print $10}')
    #split on : into fields array
    IFS=':' read -ra fields <<< "${fieldinfo}"
    # get the reference count
    ref=${fields[1]}
    # replace the split AD field in the array with the new ref,alt value
    fields[2]="$ref,$alt"

    # rejoin the array with :
    newfieldinfo=$(IFS=":"; echo "${fields[*]}")
    # replace the old split field with the new ref:alt one (we have to use a different sed delimiter (@) because some vcf fields contain / which break it
    newline=$(echo "${line}" | sed -e "s@$fieldinfo@$newfieldinfo@g")

    # print out new reformatted line
    echo "$newline"
  else
    echo "$line"
  fi
done
