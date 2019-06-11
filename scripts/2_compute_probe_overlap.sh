# ellipse radii
r1=3.5
r2=5.0

echo "name overlap angle xshift yshift" > 2_intersect_area.out
for i in ./output/*cut8.0.xyz
do
    echo $i
    ../compute_probe_overlap_2D_minimize.py $i 0.0 0.0 $r1 $r2 >> 2_intersect_area.out
done
