#!/bin/bash

# process argument
input_file="$1"

# get only the file name devoid of all arguments and extensions
file_name=$(basename $input_file)
file_name="${file_name%.*}"

# get output file name
output_file="../compare-images/$file_name.png"

# convert single dot file with multiple digraphs into pdfs
dot -Tpdf $input_file | csplit --quiet --elide-empty-files -b "%1d.pdf" --prefix=../compare-images/temp-$file_name - "/%%EOF/+1" "{*}"

# convert all the pdfs to pngs
mogrify -flatten -density 500 -resize 800 -quality 100 -format png ../compare-images/temp-$file_name*.pdf

# stitch all the pngs together
montage ../compare-images/temp-$file_name*.png -tile "2x1" -geometry +5+0 ${output_file}

# remove the temporary files
rm -f ../compare-images/temp-$file_name*
