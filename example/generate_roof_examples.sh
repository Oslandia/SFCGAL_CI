#!/usr/bin/env bash

# Script to generate OBJ files for all roof types with different geometries
# Usage: ./generate_roof_examples.sh

set -e

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== SFCGAL Roof Generation Examples ===${NC}"

# Ensure sfcgalop exists
SFCGALOP="../build/sfcgalop/sfcgalop"
if [[ ! -x "$SFCGALOP" ]]; then
    echo -e "${RED}Error: $SFCGALOP not found or not executable${NC}"
    echo "Make sure to build the project first: cmake --build build"
    exit 1
fi

# Create output directory
OUTPUT_DIR="roof_examples"
mkdir -p "$OUTPUT_DIR"

# Define test geometries
RECTANGLE='POLYGON((0 0,10 0,10 6,0 6,0 0))'
L_SHAPE='POLYGON((3 0, 3 6, 6 6, 6 3, 10 3, 10 0, 3 0))'
BIG_T_SHAPE='POLYGON((0 0,12 0,12 3,9 3,9 9,3 9,3 3,0 3,0 0))'
T_SHAPE='POLYGON((0 0, 9 0, 9 3, 6 3, 6 6, 3 6, 3 3, 0 3, 0 0))'

# Function to generate OBJ file
generate_obj() {
    local geometry="$1"
    local operation="$2"
    local params="$3"
    local filename="$4"
    local description="$5"

    echo -e "${YELLOW}Generating: ${description}${NC}"

    if $SFCGALOP -a "$geometry" "$operation" "$params" -f obj > "$OUTPUT_DIR/$filename" 2>/dev/null; then
        echo -e "  ✓ Created: $OUTPUT_DIR/$filename"
    else
        echo -e "${RED}  ✗ Failed to create: $filename${NC}"
        return 1
    fi
}

echo -e "\n${GREEN}=== RECTANGLE GEOMETRY ===${NC}"

# Rectangle - Gable Roofs
generate_obj "$RECTANGLE" "generate_gable_roof" "slope_angle=30" \
    "rect_gable_slope30.obj" \
    "Rectangle gable roof, 30° slope, roof only"

generate_obj "$RECTANGLE" "generate_gable_roof" "slope_angle=30,add_vertical_faces=1" \
    "rect_gable_slope30_faces.obj" \
    "Rectangle gable roof, 30° slope, with vertical faces"

generate_obj "$RECTANGLE" "generate_gable_roof" "slope_angle=30,building_height=3" \
    "rect_gable_slope30_building3.obj" \
    "Rectangle gable roof, 30° slope, with 3m building"

generate_obj "$RECTANGLE" "generate_gable_roof" "slope_angle=30,add_vertical_faces=1,building_height=3" \
    "rect_gable_slope30_faces_building3.obj" \
    "Rectangle gable roof, 30° slope, with faces and 3m building"

generate_obj "$RECTANGLE" "generate_gable_roof" "slope_angle=45,building_height=5" \
    "rect_gable_slope45_building5.obj" \
    "Rectangle gable roof, 45° slope, with 5m building"

# Rectangle - Flat Roofs
generate_obj "$RECTANGLE" "generate_flat_roof" "height=2" \
    "rect_flat_height2.obj" \
    "Rectangle flat roof, 2m height, roof only"

generate_obj "$RECTANGLE" "generate_flat_roof" "height=0,building_height=4" \
    "rect_flat_height0_building4.obj" \
    "Rectangle flat roof, 0m height, with 4m building"

# Rectangle - Skillion Roofs
generate_obj "$RECTANGLE" "generate_skillion_roof" "slope_angle=15" \
    "rect_skillion_slope15.obj" \
    "Rectangle skillion roof, 15° slope, roof only"

generate_obj "$RECTANGLE" "generate_skillion_roof" "slope_angle=15,add_vertical_faces=1" \
    "rect_skillion_slope15_faces.obj" \
    "Rectangle skillion roof, 15° slope, with vertical faces"

generate_obj "$RECTANGLE" "generate_skillion_roof" "slope_angle=15,building_height=3" \
    "rect_skillion_slope15_building3.obj" \
    "Rectangle skillion roof, 15° slope, with 3m building"

generate_obj "$RECTANGLE" "generate_skillion_roof" "slope_angle=25,building_height=4" \
    "rect_skillion_slope25_building4.obj" \
    "Rectangle skillion roof, 25° slope, with 4m building"

generate_obj "$RECTANGLE" "generate_skillion_roof" "slope_angle=20,add_vertical_faces=1,building_height=3" \
    "rect_skillion_slope20_faces_building3.obj" \
    "Rectangle skillion roof, 20° slope, with faces and 3m building"

# Rectangle - Hipped Roofs
generate_obj "$RECTANGLE" "generate_hipped_roof" "height=3" \
    "rect_hipped_height3.obj" \
    "Rectangle hipped roof, 3m height, roof only"

generate_obj "$RECTANGLE" "generate_hipped_roof" "height=2.5,building_height=4" \
    "rect_hipped_height2.5_building4.obj" \
    "Rectangle hipped roof, 2.5m height, with 4m building"

generate_obj "$RECTANGLE" "generate_hipped_roof" "height=4,building_height=6" \
    "rect_hipped_height4_building6.obj" \
    "Rectangle hipped roof, 4m height, with 6m building"

echo -e "\n${GREEN}=== L-SHAPE GEOMETRY ===${NC}"

# L-Shape - Gable Roofs
generate_obj "$L_SHAPE" "generate_gable_roof" "slope_angle=30" \
    "lshape_gable_slope30.obj" \
    "L-shape gable roof, 30° slope, roof only"

generate_obj "$L_SHAPE" "generate_gable_roof" "slope_angle=30,add_vertical_faces=1" \
    "lshape_gable_slope30_faces.obj" \
    "L-shape gable roof, 30° slope, with vertical faces"

generate_obj "$L_SHAPE" "generate_gable_roof" "slope_angle=35,building_height=4" \
    "lshape_gable_slope35_building4.obj" \
    "L-shape gable roof, 35° slope, with 4m building"

generate_obj "$L_SHAPE" "generate_gable_roof" "slope_angle=40,add_vertical_faces=1,building_height=3.5" \
    "lshape_gable_slope40_faces_building3.5.obj" \
    "L-shape gable roof, 40° slope, with faces and 3.5m building"

# L-Shape - Skillion Roofs
generate_obj "$L_SHAPE" "generate_skillion_roof" "slope_angle=18" \
    "lshape_skillion_slope18.obj" \
    "L-shape skillion roof, 18° slope, roof only"

generate_obj "$L_SHAPE" "generate_skillion_roof" "slope_angle=18,add_vertical_faces=1" \
    "lshape_skillion_slope18_faces.obj" \
    "L-shape skillion roof, 18° slope, with vertical faces"

generate_obj "$L_SHAPE" "generate_skillion_roof" "slope_angle=22,building_height=3.5" \
    "lshape_skillion_slope22_building3.5.obj" \
    "L-shape skillion roof, 22° slope, with 3.5m building"

generate_obj "$L_SHAPE" "generate_skillion_roof" "slope_angle=25,add_vertical_faces=1,building_height=4" \
    "lshape_skillion_slope25_faces_building4.obj" \
    "L-shape skillion roof, 25° slope, with faces and 4m building"

# L-Shape - Flat Roofs
generate_obj "$L_SHAPE" "generate_flat_roof" "height=1.8" \
    "lshape_flat_height1.8.obj" \
    "L-shape flat roof, 1.8m height, roof only"

generate_obj "$L_SHAPE" "generate_flat_roof" "height=2.2,building_height=5" \
    "lshape_flat_height2.2_building5.obj" \
    "L-shape flat roof, 2.2m height, with 5m building"

# L-Shape - Hipped Roofs
generate_obj "$L_SHAPE" "generate_hipped_roof" "height=2.8" \
    "lshape_hipped_height2.8.obj" \
    "L-shape hipped roof, 2.8m height, roof only"

generate_obj "$L_SHAPE" "generate_hipped_roof" "height=3.5,building_height=4.5" \
    "lshape_hipped_height3.5_building4.5.obj" \
    "L-shape hipped roof, 3.5m height, with 4.5m building"

echo -e "\n${GREEN}=== T-SHAPE GEOMETRY ===${NC}"

# T-Shape - Gable Roofs
generate_obj "$T_SHAPE" "generate_gable_roof" "slope_angle=30" \
    "tshape_gable_slope30.obj" \
    "T-shape gable roof, 30° slope, roof only"

generate_obj "$T_SHAPE" "generate_gable_roof" "slope_angle=30,add_vertical_faces=1" \
    "tshape_gable_slope30_faces.obj" \
    "T-shape gable roof, 30° slope, with vertical faces"

generate_obj "$T_SHAPE" "generate_gable_roof" "slope_angle=25,building_height=4" \
    "tshape_gable_slope25_building4.obj" \
    "T-shape gable roof, 25° slope, with 4m building"

generate_obj "$T_SHAPE" "generate_gable_roof" "slope_angle=35,add_vertical_faces=1,building_height=5" \
    "tshape_gable_slope35_faces_building5.obj" \
    "T-shape gable roof, 35° slope, with faces and 5m building"

# T-Shape - Skillion Roofs
generate_obj "$T_SHAPE" "generate_skillion_roof" "slope_angle=12" \
    "tshape_skillion_slope12.obj" \
    "T-shape skillion roof, 12° slope, roof only"

generate_obj "$T_SHAPE" "generate_skillion_roof" "slope_angle=12,add_vertical_faces=1" \
    "tshape_skillion_slope12_faces.obj" \
    "T-shape skillion roof, 12° slope, with vertical faces"

generate_obj "$T_SHAPE" "generate_skillion_roof" "slope_angle=20,building_height=5" \
    "tshape_skillion_slope20_building5.obj" \
    "T-shape skillion roof, 20° slope, with 5m building"

generate_obj "$T_SHAPE" "generate_skillion_roof" "slope_angle=18,add_vertical_faces=1,building_height=4" \
    "tshape_skillion_slope18_faces_building4.obj" \
    "T-shape skillion roof, 18° slope, with faces and 4m building"

# T-Shape - Flat Roofs
generate_obj "$T_SHAPE" "generate_flat_roof" "height=2" \
    "tshape_flat_height2.obj" \
    "T-shape flat roof, 2m height, roof only"

generate_obj "$T_SHAPE" "generate_flat_roof" "height=1.5,building_height=6" \
    "tshape_flat_height1.5_building6.obj" \
    "T-shape flat roof, 1.5m height, with 6m building"

# T-Shape - Hipped Roofs
generate_obj "$T_SHAPE" "generate_hipped_roof" "height=3.2" \
    "tshape_hipped_height3.2.obj" \
    "T-shape hipped roof, 3.2m height, roof only"

generate_obj "$T_SHAPE" "generate_hipped_roof" "height=2.8,building_height=5.5" \
    "tshape_hipped_height2.8_building5.5.obj" \
    "T-shape hipped roof, 2.8m height, with 5.5m building"

echo -e "\n${GREEN}=== BIG T-SHAPE GEOMETRY ===${NC}"

# Big T-Shape - Gable Roofs
generate_obj "$BIG_T_SHAPE" "generate_gable_roof" "slope_angle=30" \
    "big_tshape_gable_slope30.obj" \
    "Big T-shape gable roof, 30° slope, roof only"

generate_obj "$BIG_T_SHAPE" "generate_gable_roof" "slope_angle=30,add_vertical_faces=1" \
    "big_tshape_gable_slope30_faces.obj" \
    "Big T-shape gable roof, 30° slope, with vertical faces"

generate_obj "$BIG_T_SHAPE" "generate_gable_roof" "slope_angle=25,building_height=4" \
    "big_tshape_gable_slope25_building4.obj" \
    "Big T-shape gable roof, 25° slope, with 4m building"

generate_obj "$BIG_T_SHAPE" "generate_gable_roof" "slope_angle=35,add_vertical_faces=1,building_height=5" \
    "big_tshape_gable_slope35_faces_building5.obj" \
    "Big T-shape gable roof, 35° slope, with faces and 5m building"

# Big T-Shape - Skillion Roofs
generate_obj "$BIG_T_SHAPE" "generate_skillion_roof" "slope_angle=12" \
    "big_tshape_skillion_slope12.obj" \
    "Big T-shape skillion roof, 12° slope, roof only"

generate_obj "$BIG_T_SHAPE" "generate_skillion_roof" "slope_angle=12,add_vertical_faces=1" \
    "big_tshape_skillion_slope12_faces.obj" \
    "Big T-shape skillion roof, 12° slope, with vertical faces"

generate_obj "$BIG_T_SHAPE" "generate_skillion_roof" "slope_angle=20,building_height=5" \
    "big_tshape_skillion_slope20_building5.obj" \
    "Big T-shape skillion roof, 20° slope, with 5m building"

generate_obj "$BIG_T_SHAPE" "generate_skillion_roof" "slope_angle=18,add_vertical_faces=1,building_height=4" \
    "big_tshape_skillion_slope18_faces_building4.obj" \
    "Big T-shape skillion roof, 18° slope, with faces and 4m building"

# Big T-Shape - Flat Roofs
generate_obj "$BIG_T_SHAPE" "generate_flat_roof" "height=2" \
    "big_tshape_flat_height2.obj" \
    "Big T-shape flat roof, 2m height, roof only"

generate_obj "$BIG_T_SHAPE" "generate_flat_roof" "height=1.5,building_height=6" \
    "big_tshape_flat_height1.5_building6.obj" \
    "Big T-shape flat roof, 1.5m height, with 6m building"

# Big T-Shape - Hipped Roofs
generate_obj "$BIG_T_SHAPE" "generate_hipped_roof" "height=3.2" \
    "big_tshape_hipped_height3.2.obj" \
    "Big T-shape hipped roof, 3.2m height, roof only"

generate_obj "$BIG_T_SHAPE" "generate_hipped_roof" "height=2.8,building_height=5.5" \
    "big_tshape_hipped_height2.8_building5.5.obj" \
    "Big T-shape hipped roof, 2.8m height, with 5.5m building"


# Additional interesting combinations
echo -e "\n${GREEN}=== ADDITIONAL VARIATIONS ===${NC}"

# Extreme slope angles
generate_obj "$RECTANGLE" "generate_gable_roof" "slope_angle=15,building_height=4" \
    "rect_gable_slope15_building4.obj" \
    "Rectangle gable roof, low 15° slope, with 4m building"

generate_obj "$RECTANGLE" "generate_gable_roof" "slope_angle=60,building_height=3" \
    "rect_gable_slope60_building3.obj" \
    "Rectangle gable roof, steep 60° slope, with 3m building"

# High buildings
generate_obj "$L_SHAPE" "generate_hipped_roof" "height=5,building_height=8" \
    "lshape_hipped_height5_building8.obj" \
    "L-shape hipped roof, 5m height, with 8m building"

generate_obj "$T_SHAPE" "generate_flat_roof" "height=3,building_height=10" \
    "tshape_flat_height3_building10.obj" \
    "T-shape flat roof, 3m height, with 10m building"

# Different ridge edges for skillion roofs
generate_obj "$RECTANGLE" "generate_skillion_roof" "slope_angle=20,ridge_edge=1,building_height=3" \
    "rect_skillion_ridge1_slope20_building3.obj" \
    "Rectangle skillion roof, ridge on edge 1, 20° slope, with 3m building"

generate_obj "$RECTANGLE" "generate_skillion_roof" "slope_angle=20,ridge_edge=2,add_vertical_faces=1" \
    "rect_skillion_ridge2_slope20_faces.obj" \
    "Rectangle skillion roof, ridge on edge 2, 20° slope, with vertical faces"

generate_obj "$L_SHAPE" "generate_skillion_roof" "slope_angle=18,ridge_edge=1,add_vertical_faces=1,building_height=2" \
    "lshape_skillion_ridge1_slope18_faces_building2.obj" \
    "L-shape skillion roof, ridge on edge 1, 18° slope, with faces and 2m building"

echo -e "\n${GREEN}=== SUMMARY ===${NC}"
echo "Generated OBJ files in directory: $OUTPUT_DIR/"
echo -e "Total files created: $(find "$OUTPUT_DIR" -name "*.obj" | wc -l)"
echo
echo -e "${YELLOW}You can now open these OBJ files in 3D viewers like:${NC}"
echo "  - Blender"
echo "  - MeshLab"
echo "  - Online viewers (e.g., https://3dviewer.net/)"
echo "  - ParaView"
echo
echo -e "${YELLOW}Example geometries:${NC}"
echo "  - Rectangle: 10m × 6m rectangular footprint"
echo "  - L-Shape: Complex L-shaped building footprint"
echo "  - T-Shape: T-shaped building footprint"
echo
echo -e "${GREEN}Generation complete!${NC}"
