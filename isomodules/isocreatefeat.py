#!/usr/bin/env python
# title           :isocreate.py
# description     :Functions to create features mapped to isoforms.
# author          :Gloria Sheynkman
# date            :June 17th, 2019
# version         :1
# python_version  :2.6.6
# ==============================================================================

from collections import defaultdict
import isofeature

# *****************************************************************************
# create and map features to iso_obj
def create_and_map_domain(orf, domain_info, with_evals=True):
    """Given an orf and info. on a domain mapped to its AAs, create a
       hierarchical feat_obj mapping.
       Input:
        orf - orf_obj
        domain_info - [pfam_acc, pfam_name, cat (DBD, other), 1-base_start,
                       1-base_end]
      Return (optional):
       featf (domain object)
    """
    if with_evals:
        pfam, name, cat, eval, start, end = domain_info
    else:
        pfam, name, cat, start, end = domain_info
        eval = 'NA'
    featf = isofeature.FeatureFull(orf, 'dom', cat, eval=eval, acc=pfam, desc=name)
    domain_ress = orf.res_chain[int(start) - 1:int(end)]
    for res in domain_ress:
        featr = isofeature.FeatureResidue(featf, res)
        featf.chain.append(featr)
    featr_groupings = find_groups_of_featr_in_subblocks(featf.chain)
    for acc in sorted(featr_groupings):
        cds_ord, cds = acc
        featrs = featr_groupings[acc]
        featsb = isofeature.FeatureSubblock(featf, cds, featrs)
        # update upwards/downwards references
        featf.subblocks.append(featsb)
        for featr in featrs:
            featr.featsb = featsb
    return featf

def find_groups_of_featr_in_subblocks(featrs):
    """Given a list of featr (repr. domain mapped to orf), group the featrs
       into those mapping to the same cds.

       Return:
        {(cds_ordinal, cds_obj) -> [featr, featr, featr]}
       cds_ordinal allows for sorting of keys so featsb can be populated
       into featf in upstream-to-downstream order
    """
    sb_grps = defaultdict(list)
    for featr in featrs:
        cds_ord = featr.res.cds.ord
        cds = featr.res.cds
        acc = (cds_ord, cds)
        sb_grps[acc].append(featr)
    return sb_grps
