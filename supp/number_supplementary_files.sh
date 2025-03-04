#!/bin/bash

zip Data_S1.zip SuppTable_CloneList.txt SuppTable_CloningPrimers.txt SuppTable_SinglesCloneList.txt SuppTable_DBD_definition.txt
zip Data_S2.zip SuppTable_eY1HResults.txt SuppTable_DNABaits.txt SuppTable_PDI_validation.txt
cp SuppTable_M1HResults.txt Data_S3.tsv
zip Data_S4.zip SuppTable_PairwiseY2HResults.txt SuppTable_N2HResults.txt
zip Data_S5.zip SuppTable_CREB1-PBM.txt SuppTable_TBX5-PBM.txt
cp SuppTable_Paralogs.txt Data_S6.tsv
cp SuppTable_Condensates.txt Data_S7.tsv
cp SuppTable_NegRegs.txt Data_S8.tsv
zip Data_S9.zip SuppTable_TCGASamps.txt SuppTable_TCGAResults.txt
cp SuppTable_JoungResults.txt Data_S10.tsv
