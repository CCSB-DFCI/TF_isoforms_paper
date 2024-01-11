reinitialize
bg_color white
load ../data/internal/alphafold/monomers/THRA-2/ranked_0.pdb, THRA_2
load ../data/internal/alphafold/monomers/THRA-1/ranked_0.pdb, THRA_1
load ../data/internal/alphafold/monomers/THRB-2/ranked_0.pdb, THRB_2
load ../data/internal/alphafold/monomers/THRB-1/ranked_0.pdb, THRB_1

align THRA_1, THRA_2
align THRB_2, THRA_2
align THRB_1, THRB_2

set grid_mode, 1
set grid_slot, 1, THRA_2
set grid_slot, 2, THRA_1
set grid_slot, 3, THRB_2
set grid_slot, 4, THRB_1
refresh
viewport 2000, 1000
zoom

color grey, THRA_2
color grey, THRA_1
color grey, THRB_2
color grey, THRB_1

# color regions that are different between isoforms or paralogs
color 0xd95f02, THRA_2 AND resi 1-41
color 0xd95f02, THRA_1 AND resi 1-41
color red, THRA_1 AND resi 370-451
color 0xe7298a, THRB_2 AND resi 1-95
color 0x66a61e, THRB_1 AND resi 1-110


set depth_cue, 0
set cartoon_loop_radius, 0.3
reset
set antialias,2
set hash_max, 300
set ray_shadows, 0
set ray_trace_fog, 0
set ray_opaque_background, off

# run these commands after positioning structure
#ray 2000,2000
#png ../figures/alphafold_THR_paralogs-vs-isoforms.png, dpi=300