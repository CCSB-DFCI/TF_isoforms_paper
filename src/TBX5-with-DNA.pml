reinitialize
bg_color white
load ../data/internal/alphafold/monomers/TBX5-1/ranked_0.pdb, TBX5_1
fetch 5FLV
align TBX5_1 and resi 56-238, 5FLV and chain A and resi 1004-1238
hide cartoon, 5FLV
show cartoon, 5FLV and chain B
show cartoon, 5FLV and chain C
color forest, 5FLV and chain B
color forest, 5FLV and chain C
color deepblue, TBX5_1
color 0xd95f02, TBX5_1 and resi 1-50
color 0xe7298a, TBX5_1 and resi 328-518
remove solvent
remove organic
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
#png ../figures/TBX5-1_with-DNA.png, dpi=300
#rotate y, 90
#ray 2000,2000
#png ../figures/TBX5-1_with-DNA_rotate-y-90.png, dpi=300