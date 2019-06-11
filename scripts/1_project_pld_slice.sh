#!/bin/bash -l

cat 0_parsed_output.txt | while read S
do
    s=( $S )
    structure=${s[0]}
    axis=${s[1]}
    node1=${s[2]}
    node2=${s[3]}
    x1=${s[4]}
    y1=${s[5]}
    z1=${s[6]}
    x2=${s[7]}
    y2=${s[8]}
    z2=${s[9]}
    dx=${s[10]}
    dy=${s[11]}
    dz=${s[12]}
    rad=${s[13]}
    a=${s[14]}
    b=${s[15]}
    c=${s[16]}
    alpha=${s[17]}
    beta=${s[18]}
    gamma=${s[19]}
    echo $structure >> 1_verbose.txt
    echo "arguments:  ./ ${structure} $axis 8.0 ./UFF.rad $x1 $y1 $z1 $x2 $y2 $z2 $dx $dy $dz $a $b $c $alpha $beta $gamma "
    echo "" >> 1_verbose.txt
    ../projection_pld_slice_PBC.py ./ ${structure} $axis 8.0 ./UFF.rad $x1 $y1 $z1 $x2 $y2 $z2 $dx $dy $dz $a $b $c $alpha $beta $gamma >> 1_verbose.txt
    echo "" >> 1_verbose.txt
done

