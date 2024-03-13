reinitialize
bg_color white
load ../data/internal/alphafold/monomers/HEY1-1/ranked_0.pdb, HEY1_1
fetch 4H10
align HEY1_1 and resi 50-107, 4H10 and chain A
hide cartoon, 4H10 and chain A
color grey70, 4H10 and chain B
color forest, 4H10 and chain C
color forest, 4H10 and chain D
color deepblue, HEY1_1
color red, HEY1_1 and resi 84-87
remove solvent
set depth_cue, 0
set cartoon_loop_radius, 0.3
reset
set antialias,2
set hash_max, 300
set ray_shadows, 0
set ray_trace_fog, 0
set ray_opaque_background, off
# run these commands after positioning structure
#ray 1600,2000
#png ../figures/HEY1-1_with-DNA-and-dimer.png, dpi=300
#rotate y, 90
#ray 2000,2000
#png ../figures/HEY1-1_with-DNA-and-dimer_rotate-y-90.png, dpi=300