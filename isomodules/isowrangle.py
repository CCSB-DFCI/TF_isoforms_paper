# functions that wrangle (reshape, regroup, restructure) iso-objects


import isoclass
import isogroup

def get_same_prot_seq_grps_of_orfs_from_gene_(gen_obj):
    """Given a gen_obj dict., find groups of ORFs with the same protein seq.
       and return as same-seq-grp objects.
    """
    same_seq_orfs = gen_obj.same_seq_orfs
    grps = []  # grps of same-protein-seq orf_objs
    for orfs in same_seq_orfs:
        grp = isogroup.SameSequenceGencodeGroup(orfs)
        grps.append(grp)
    return grps
