rename table horfeome_cdna_mapping_hg38 to horfeome_cdna_mapping_hg38_orig ;
create table horfeome_cdna_mapping_hg38 like horfeome_cdna_mapping_hg38_orig;



alter table orf_utr_hg38_blat_best  add orf_id int(8) , add gi int(12);

update orf_utr_hg38_blat_best set orf_id = substr(qName , 8, locate ('|', qName) -8 ), gi = substr(qName, locate('gi=', qName)+3);
alter table orf_utr_hg38_blat_best add index (orf_id), add index (gi) ;

insert into horfeome_cdna_mapping_hg38 (orf_id, gi, chrom, strand, ref_start, ref_end, cdna_start, cdna_end) select distinct a.orf_id , a.gi , a.chrom , a.strand, ref_start, ref_end, cdna_start, cdna_end from horfeome_cdna_mapping_hg38_orig a, orf_utr_hg38_blat_best b, orf_utr_hg38_blat_exon c  where a.orf_id = b.orf_id and a.gi = b.gi and a.chrom = tName  and a.ref_start >= c.tStart-10 and ref_end <= c.tEnd+10 and a.strand = b.strand and a.strand = '+' and b.LINE_ID = c.LINE_ID and (a.ref_start  = c.tStart or a.ref_end = c.tEnd ) ;


insert into horfeome_cdna_mapping_hg38 (orf_id, gi, chrom, strand, ref_start, ref_end, cdna_start, cdna_end)   select distinct a.orf_id , a.gi , a.chrom , a.strand, ref_start, ref_end, cdna_start, cdna_end from horfeome_cdna_mapping_hg38_orig a, orf_utr_hg38_blat_best b, orf_utr_hg38_blat_exon c  where a.orf_id = b.orf_id and a.gi = b.gi and  a.chrom = tName  and a.ref_start >= c.tEnd-10  and a.ref_end<= c.tStart+10 and a.strand = b.strand and a.strand = '-' and b.LINE_ID = c.LINE_ID and (a.ref_start  = c.tEnd or a.ref_end = c.tStart) ;
