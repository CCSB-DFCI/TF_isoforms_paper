#!/bin/bash
rsync -r --progress --include='ranked_0*' --include='*/' --exclude='*' ll309@transfer.rc.hms.harvard.edu:/n/data2/dfci/genetics/vidal/alphafold/TF_isoforms/monomers/ data/alphafold
