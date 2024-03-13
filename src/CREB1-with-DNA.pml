reinitialize
bg_color white
load ../data/internal/alphafold/monomers/CREB1-1/ranked_0.pdb, CREB1_1
fetch 1DH3
align CREB1_1 and resi 281-340, 1DH3 and chain A
hide cartoon, 1DH3 and chain A
color grey70, 1DH3 and chain C
color forest, 1DH3 and chain B
color forest, 1DH3 and chain D
color deepblue, CREB1_1
color red, CREB1_1 and resi 88-101
remove solvent
remove element mg
set depth_cue, 0
reset
set cartoon_loop_radius, 0.3
reset
set antialias,2
set hash_max, 300
set ray_shadows, 0
set ray_trace_fog, 0
set ray_opaque_background, off
# run these commands after positioning structure
#ray 1600,2000
#png ../figures/CREB1-1_with-DNA-and-dimer.png, dpi=300
#rotate y, 90
#ray 1600,2000
#png ../figures/CREB1-1_with-DNA-and-dimer_rotate-y-90.png, dpi=300
