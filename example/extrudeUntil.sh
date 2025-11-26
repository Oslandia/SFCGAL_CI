#!/bin/sh
# Définition des polygones
RECTANGLE='POLYGON((0 0,10 0,10 6,0 6,0 0))'
L_SHAPE='POLYGON((3 0,3 6,6 6,6 3,10 3,10 0,3 0))'
BIG_T_SHAPE='POLYGON((0 0,12 0,12 3,9 3,9 9,3 9,3 3,0 3,0 0))'
T_SHAPE='POLYGON((0 0,9 0,9 3,6 3,6 6,3 6,3 3,0 3,0 0))'

# Types de toits et leurs paramètres
# format: "operation:paramètres"
ROOFS="
generate_flat_roof:height=2.0
generate_hipped_roof:height=3.0
generate_gable_roof:slope_angle=30
generate_skillion_roof:slope_angle=15
"

# Itération sur les polygones
for poly_name in RECTANGLE L_SHAPE BIG_T_SHAPE T_SHAPE; do
    poly=$(eval echo \$$poly_name)
    echo "=============================="
    echo "Traitement du polygone : $poly_name"
    echo "WKT : $poly"

    # Export du footprint
    echo "Export du footprint : ${poly_name}_footprint.obj"
    # echo "./build/sfcgalop/sfcgalop -a \"$poly\" -f obj > ${poly_name}_footprint.obj"
    ./build/sfcgalop/sfcgalop -a "$poly" -f obj > "${poly_name}_footprint.obj"

    # Itération sur les types de toits
    echo "$ROOFS" | while IFS=: read roof_type roof_param; do
        if [ -z "$roof_type" ]; then continue; fi
        echo "------------------------------"
        #echo "Génération du toit : $roof_type avec $roof_param pour $poly_name"

        # Génération du toit et translation
        echo "Appel à sfcgalop pour générer le toit avec translation dz=5"
        #echo "./build/sfcgalop/sfcgalop -a \"$poly\" \"$roof_type\" \"$roof_param\" | ./build/sfcgalop/sfcgalop translate \"dz=5\" --"
        roof=$(./build/sfcgalop/sfcgalop -a "$poly" "$roof_type" "$roof_param" | \
               ./build/sfcgalop/sfcgalop translate "dz=5" --)

        # Export du toit
        echo "Export du toit : ${poly_name}_${roof_type}_roof.obj"
        #echo "./build/sfcgalop/sfcgalop -a \"$poly\" \"$roof_type\" \"$roof_param\" | ./build/sfcgalop/sfcgalop translate \"dz=5\" -f obj > ${poly_name}_${roof_type}_roof.obj"
        ./build/sfcgalop/sfcgalop -a "$poly" "$roof_type" "$roof_param" | \
            ./build/sfcgalop/sfcgalop translate "dz=5" -f obj > "${poly_name}_${roof_type}_roof.obj"

        # Génération du solide avec extrudeUntil
        echo "Génération du solide : ${poly_name}_${roof_type}_solid.obj"
        #echo "./build/sfcgalop/sfcgalop -a \"$poly\" -b \"$roof\" extrudeUntil -f obj > ${poly_name}_${roof_type}_solid.obj"
        ./build/sfcgalop/sfcgalop -a "$poly" -b "$roof" extrudeUntil -f obj > "${poly_name}_${roof_type}_solid.obj"
    done
done

