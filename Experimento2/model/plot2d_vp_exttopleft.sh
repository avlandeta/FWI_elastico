#!/bin/bash

# Define the binary file
binfile=${1}
# binfile="marmousi_II_marine.vp"
echo "Binary file: $binfile"

# Extract the base name and extension of the binary file
fbase=$(basename "$binfile")
fext="${fbase##*.}"
fbase="${fbase%.*}"
echo "Base name: $fbase"
echo "Extension: $fext"

# redefine the base name for include the parameter
fbase="${fbase}_vp"

# Set dimensions and ranges
n1=362
n2=1740

d1=10.0
d2=10.0

x1min=0.0
x1max=3610.0

x2min=0.0
x2max=17390.0

o1=0.0
o2=0.0

label="Vp(m/s)"
vmin=1500.0
vmax=4800
vbias=$vmin

# Convert the binary file to RSF format using sfbin2rsf
sfbin2rsf bfile=${binfile} n1=${n1} n2=${n2} d1=${d1} d2=${d2} o1=${o1} o2={o2} | \
sfput o1=${o1} o2=${o2} d1=${d1} d2=${d2} label1="Depth" label2="Distance" unit1="m" unit2="m" > ${fbase}.rsf

# Process the data using the specified commands directly on the RSF file
sfwindow < ${fbase}.rsf min1=${x1min} max1=${x1max} \
                       min2=${x2min} max2=${x2max} | \
sfbyte bar=bvel2d.rsf point1=0.5= point2=0.75 obarnum=${vmin} allpos=y \
       bias=${vbias} minval=${vmin} maxval=${vmax}  barlabel=${label}  > ${fbase}_byte.rsf

# Calculate the aspect ratio
# aspect_ratio=$(echo "$n2 / $n1" | bc -l)
# scrt=${aspect_ratio} ;
# echo $scrt
# # Set the desired width (adjust as necessary)
# width=15
# scwd=${width}
# echo $scwd
# # Calculate the height based on the aspect ratio
# height=$(echo "$width / $aspect_ratio" | bc -l)
# scht=${height}
# echo $scht
# # Set the desired height (adjust as necessary)
# scht=9
# echo $scht
# # Calculate the width based on the aspect ratio
# scwd=$(echo " $scht * $aspect_ratio" | bc)
# echo $scwd

scrt=0.5; 
scwd=14;
scht=7;
# Create a visualization while keeping the aspect ratio
sfgrey < ${fbase}_byte.rsf bar=bvel2d.rsf point1=0.5 point2=0.75 \
       minval=${vmin} maxval=${vmax} bias=${vbias} barlabel=${label} allpos=y color=j scalebar=y \
       screenht=$scht screenwd=$scwd screeratio=$scrt title="" > ${fbase}.vpl

# Convert the visualization file to JPEG format
vpconvert ${fbase}.vpl format=pdf 
vpconvert ${fbase}.vpl format=jpeg 
#pdfcrop ${bname}.pdf
#mv ${bname}-crop.pdf ${bname}.pdf
#rm tmp.bin
#echo "Created: ${outfig}.pdf"
