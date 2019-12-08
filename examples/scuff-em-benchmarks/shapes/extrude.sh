#!/bin/sh

sed -e "s/RADIUS/30./g" -e "s/HEIGHT/30./g" -e "s/ELEMSIZ/4./g" cylinder.geo.template >tmp.geo \
        && gmsh tmp.geo -o cylinder_r30_h30_fine.msh -2 

sed -e "s/RADIUS/30./g" -e "s/HEIGHT/30./g" -e "s/ELEMSIZ/7./g" cylinder.geo.template >tmp.geo \
        && gmsh tmp.geo -o cylinder_r30_h30.msh -2 

sed -e "s/RADIUS/30./g" -e "s/HEIGHT/30./g" -e "s/ELEMSIZ/15./g" cylinder.geo.template >tmp.geo \
        && gmsh tmp.geo -o cylinder_r30_h30_rough.msh -2 

sed -e "s/RADIUS/30./g" -e "s/HEIGHT/30./g" -e "s/ELEMSIZ/25./g" cylinder.geo.template >tmp.geo \
        && gmsh tmp.geo -o cylinder_r30_h30_veryrough.msh -2 

sed -e "s/RADIUS/100./g" -e "s/HEIGHT/50./g" -e "s/ELEMSIZ/13.3/g" cylinder.geo.template >tmp.geo \
        && gmsh tmp.geo -o cylinder_r100_h50_fine.msh -2 

sed -e "s/RADIUS/100./g" -e "s/HEIGHT/50./g" -e "s/ELEMSIZ/23.3/g" cylinder.geo.template >tmp.geo \
        && gmsh tmp.geo -o cylinder_r100_h50.msh -2 

sed -e "s/RADIUS/100./g" -e "s/HEIGHT/50./g" -e "s/ELEMSIZ/50./g" cylinder.geo.template >tmp.geo \
        && gmsh tmp.geo -o cylinder_r100_h50_rough.msh -2 

sed -e "s/RADIUS/100./g" -e "s/HEIGHT/50./g" -e "s/ELEMSIZ/83.3/g" cylinder.geo.template >tmp.geo \
        && gmsh tmp.geo -o cylinder_r100_h50_veryrough.msh -2 

rm tmp.geo

