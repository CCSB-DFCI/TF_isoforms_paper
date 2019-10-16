package MY_SeqAlign;
use strict;
use DBI;
#use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO::blast;
#use Bio::Tools::BPbl2seq;
use Bio::Seq;
use Bio::SeqIO;
use POSIX qw(ceil floor);
use Switch;
use lib '/home/yun/lib/';
use MY_DButility;
use MY_BioUtility;
use Cwd;
use constant MAX_GAP    => 10000000000000;
use constant MIN_GAP    => 0;
use constant MIN_HIT_OVERLAP => 10;
use constant MIN_BL_SCORE => 100;

### define constant

my $MAX_WELL_PER_PRIMER_PLATE = 96; # ??? shoud leave 2 well empty for control
my $MAX_WELL_PER_CONSO_PLATE = 94; #leave g12, h12 to be control well
my $TOTAL_WELL_PER_PLATE=96;


#*************************************************************************

##########################################################################
# BEGIN of Bioinformatics functions
##########################################################################

##########################################################################
# CROSS MATCH
##########################################################################
sub run_vector_trim {
  my ($sql, $query, $rtn);
  my ($pkg, $tables, $rtinfo) = @_;
  
  my $result_dir = defined $rtinfo->{RESULT_DIR}?$rtinfo->{RESULT_DIR}:"";
  if(!-e $result_dir) {
    system("mkdir -p $result_dir");
  }
  my $tmp_dir = defined $rtinfo->{TMP_DIR}?$rtinfo->{TMP_DIR}:"";
  if(!-e $tmp_dir) {
    system("mkdir -p $tmp_dir");
  }

# get vector fasta
  $sql = "select SEQ_NAME, DIR, VECTOR_SEQ from $tables->{VECTOR}{TABLENAME} a ";     
  $query = $tables->{VECTOR}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  my $vector_fasta = defined($rtinfo->{VECTOR_FASTA}) ? $tmp_dir . $rtinfo->{VECTOR_FASTA}:$tmp_dir . "vector.fa";
  my %vector = ();
  open OUT, ">$vector_fasta";
  while (my ($SEQ_NAME, $DIR, $SEQ) = $query -> fetchrow()) {
    print OUT ">$SEQ_NAME\n";
    print OUT "$SEQ\n";
    $vector{$SEQ_NAME}{SEQ} = $SEQ;
    $vector{$SEQ_NAME}{DIR} = $DIR;
  }
  close OUT;
  $query -> finish;

# get clone seq fasta:
  $sql = "select TRACE_ID, TRACE_DIR, TRACE_SEQ, TRACE_QUAL from $tables->{TRACE_SEQ}{TABLENAME} " . (defined($tables->{TRACE_SEQ}{FILTER})? "where $tables->{TRACE_SEQ}{FILTER} ":"") . "order by TRACE_ID";
  $query = $tables->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();

  my $trace_fasta = defined($rtinfo->{TRACE_FASTA}) ? $tmp_dir .$rtinfo->{TRACE_FASTA}:$tmp_dir . "trace.fa";
  my $trace_qual = defined($rtinfo->{TRACE_FASTA}) ? $tmp_dir .$rtinfo->{TRACE_FASTA} . ".qual" : $tmp_dir . "trace.fa.qual";
  open OUT, ">$trace_fasta";
  open QUAL, ">$trace_qual";
  while(my ($TRACE_ID, $TRACE_DIR, $TRACE_SEQ, $TRACE_QUAL) = $query -> fetchrow()) {
    print OUT ">$TRACE_ID" . "_" . $TRACE_DIR . "\n";
    print OUT "$TRACE_SEQ\n";    
    print QUAL ">$TRACE_ID" . "_" . $TRACE_DIR . "\n";
    print QUAL "$TRACE_QUAL\n";    
  }
  $query -> finish;  
  close OUT;
  close QUAL;

# run cross match
  my $screen_out = $tmp_dir . "crossmatch.out";
  system("cross_match $trace_fasta $vector_fasta -masklevel 101 -screen > $screen_out");

# parse screen file:
#  $sql = "update $tables->{TRACE_SEQ}{TABLENAME} set $tables->{TRACE_SEQ}{VTRIM_START_FIELD}=? where TRACE_ID=? and ($tables->{TRACE_SEQ}{VTRIM_START_FIELD} is null or $tables->{TRACE_SEQ}{VTRIM_START_FIELD}<?)";
#  my $qu_start = $tables->{TRACE_SEQ}{DBH} -> prepare($sql);

#  $sql = "update $tables->{TRACE_SEQ}{TABLENAME} set $tables->{TRACE_SEQ}{VTRIM_END_FIELD}=? where TRACE_ID=? and ($tables->{TRACE_SEQ}{VTRIM_END_FIELD} is null or $tables->{TRACE_SEQ}{VTRIM_END_FIELD}>?)";
#  my $qu_end = $tables->{TRACE_SEQ}{DBH} -> prepare($sql);

  $sql = "update $tables->{TRACE_SEQ}{TABLENAME} set " . (defined($tables->{TRACE_SEQ}{VTRIM_START_FIELD})?$tables->{TRACE_SEQ}{VTRIM_START_FIELD}:"VTRIM_START") . "=? where TRACE_ID=?";
  my $qu_start = $tables->{TRACE_SEQ}{DBH} -> prepare($sql);

  $sql = "update $tables->{TRACE_SEQ}{TABLENAME} set " . (defined($tables->{TRACE_SEQ}{VTRIM_END_FIELD})?$tables->{TRACE_SEQ}{VTRIM_END_FIELD}:"VTRIM_END") . "=? where TRACE_ID=?";
  my $qu_end = $tables->{TRACE_SEQ}{DBH} -> prepare($sql);

# cross match output Interpretation:
# Example:
# 440  2.38 1.39 0.79  hh44a1.s1       33   536 (    0)  C 00311     ( 3084)  8277   7771  *
# item[0] : 440, smith-waterman score of the match (complexity-adjusted, by default).
# item[1] :  2.38, %substitutions in matching region
# item[2] :  1.39, %deletions (in 1st seq rel to 2d) in matching region
# item[3] :  0.79,  %insertions (in 1st seq rel to 2d) in matching region 
# item[4] :  hh44a1.s1, id of 1st sequence
# item[5] :  33, starting position of match in 1st sequence
# item[6] :  536, ending position of match in 1st sequence
# item[7] :  (0), no. of bases in 1st sequence past the ending position of match
#         (so 0 means that the match extended all the way to the end of 
#           the 1st sequence)
# item[8] : "C 00311",  match is with the Complement of sequence 00311
#           (3084) : there are 3084 bases in (complement of) 2d sequence prior to 
#        beginning of the match
# item[9]:  8277 = starting position of match in 2d sequence (using top-strand 
# item[10]:         numbering)
# item[11]:  7771 =  ending position of match in 2d sequence
#  * indicates that there is a higher-scoring match whose domain partly
#           includes the domain of this match.


  my ($for_5_count, $for_3_count) = (0,0);
  open IN, "<$screen_out";
  my %best_trim =();
  my $vector_string = join '|', keys %vector;
  while (<IN>) {
    if(/$vector_string/) {
        chomp;
        my ($score) = $_ =~ /^\s+(\d+)\s/;
        $_ =~ s/^(\s+\d+\s)//;
        my @items = split /\s\s+/;
        $items[3] =~ s/(\(\d+\))//; 
        my ($id, $dir) = $items[1]=~ /(\d+)_(FOR|REV)/;
#        next if($dir eq 'FOR' && $items[4]=~ /C $vector_string/);
#        next if($dir eq 'REV' && !$items[4]=~ /C $vector_string/); 
     
        if(!defined $best_trim{TRACE}) {
            $best_trim{TRACE} = $items[1];
            $best_trim{VECTOR} = $items[4];
            $best_trim{SCORE} = $score; 
            $best_trim{START} = $items[2];
            $best_trim{END} = $items[3];
            $best_trim{ID} = $id;
            $best_trim{DIR} = $dir;
        }
        elsif($best_trim{TRACE} eq $items[1] && $best_trim{VECTOR} eq $items[4]) {
          if($best_trim{SCORE}<$score) { # get the best trim so far
            $best_trim{SCORE} = $score; # update best trim score
            $best_trim{START} = $items[2];
            $best_trim{END} = $items[3];
          }
        }

        else { # before move on, update the current best trim first:
          foreach my $SEQ_NAME (keys %vector) {
              if ($best_trim{DIR} eq $vector{$SEQ_NAME}{DIR} && $best_trim{VECTOR} eq $SEQ_NAME){
		$qu_start->execute($best_trim{END}+1, $best_trim{ID}); 
	        $for_5_count++;
              }
              elsif ($best_trim{DIR} ne $vector{$SEQ_NAME}{DIR} && $best_trim{VECTOR} eq 'C ' . $SEQ_NAME) {#Forward sequence matches reverse tail anti-sense and determines insert end
		$qu_end->execute($best_trim{START}+1, $best_trim{ID});
                $for_3_count++;
              }
         } # end of parsing vector

# reset best_trim to the current value      
            %best_trim = (); # reset first
#=pod
            $best_trim{TRACE} = $items[1];
            $best_trim{VECTOR} = $items[4];
            $best_trim{ID} = $id;
            $best_trim{DIR} = $dir;
            $best_trim{SCORE} = $score;
            $best_trim{START} = $items[2];
            $best_trim{END} = $items[3];
#=cut
        } # end of update database and move on
     } # end of finding matched line
  } # end of read in file
  close IN;

# update the last trim result:
          foreach my $SEQ_NAME (keys %vector) {
              if ($best_trim{DIR} eq $vector{$SEQ_NAME}{DIR} && $best_trim{VECTOR} eq $SEQ_NAME){
		$qu_start->execute($best_trim{END}+1, $best_trim{ID}); 
	        $for_5_count++;
              }
              elsif ($best_trim{DIR} ne $vector{$SEQ_NAME}{DIR} && $best_trim{VECTOR} eq 'C ' . $SEQ_NAME) {#Forward sequence matches reverse tail anti-sense and determines insert end
		$qu_end->execute($best_trim{START}+1, $best_trim{ID});
                $for_3_count++;
              }
        } # end of parsing vector
  
  print "start $for_5_count\nend $for_3_count\n";
#  system("rm $trace_fasta $vector_fasta $screen_out");

} # end of run_vector_trim

##########################################################################
# END of Bioinformatics functions
##########################################################################

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# BEGIN of Helper function for load data DB
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



##########################################################################
# create necessary result tables
# xxx_TRACE_PLATEMAP,
# xxx_CHERRYMAP
# xxx_TRACE_SEQ
##########################################################################
sub create_result_tables {
  my ($pkg, $dbh, $dbname, $trace_plates, $result_tables, $fields)=@_;
 
  my ($sql, $query, $rtn);
  my $exist;
  my $filter = "TRACE_PLA in ( ". join(',', @{$trace_plates}) . ")";

  my $nonstd_platemap_type_table = 'yun_template_table.G_NONSTD_PLATEMAP_TYPE';

# remove all the records in xxx_TRACE_PLATEMAP assoc. w/ @{$trace_plates}
  $exist = MY_DButility->check_table_existence($dbh, $dbname, $result_tables->{PLATE_MAP_TABLE});
  if($exist == 0) { # not exists
    $sql = "create table $result_tables->{PLATE_MAP_TABLE} ( "
      . "PROJ_ABBR char(3), "
      . "TRACE_PLA int(6), "  # ref. to tbPlate(PLATE_ID)
      . "TEMPL_PLA_NAME varchar(50), " # ref. to tbPlate(PLATE_NAME)
      . "SEQ_TYPE enum('PCR', '3RACE', '5RACE', 'RTPOOL', 'ISOCOL', 'PCRPOOL'), "
      . "IS_STD_MAP enum('NO', 'YES') default 'YES', "
      . "NONSTD_MAP_TYPE int(2), "
      . "NOTE varchar(50), "
      . "CREATE_TIME timestamp default CURRENT_TIMESTAMP, "
      . "PRIMARY KEY(TRACE_PLA), "
      . "FOREIGN KEY(NONSTD_MAP_TYPE) references $nonstd_platemap_type_table (ID) on update cascade on delete restrict "
      . ");";
    $rtn = $dbh->do($sql); 
  }

# make sure CHERRYMAP always there for use
  $exist = MY_DButility->check_table_existence($dbh, $dbname, $result_tables->{CHERRYMAP_TABLE});
  if($exist == 0) { # not exists
    printf "warning - no CHERRYMAP table is ready, please create one first. \n";
    return -1;
  }

# make sure PLATE_DETAIL is there for use
  $exist = MY_DButility->check_table_existence($dbh, $dbname, $result_tables->{PLATE_DETAIL_TABLE});
  if($exist == 0) { # not exists
    printf "warning - no PLATE_DETAIL table is ready, please create one first. \n";
    return -1;
  }

  $exist = MY_DButility->check_column_existence($dbh, $dbname, $result_tables->{CHERRYMAP_TABLE}, $fields->{RUN_TRACE_PLA});
  if($exist==0 ) { # column not exists
    $sql = "alter table " . $result_tables->{CHERRYMAP_TABLE} . " add ( "
      . $fields->{RUN_TRACE_PLA} . " int(6), " #trace plate id assigned after blast in trace_upload
      . $fields->{RUN_TRACE_POS} . " char(3), "
#      . $fields->{RUN_TRACE_POSID} . " int(3), "
      . $fields->{RUN_BL2_RESULT} . " enum('BOTH_BAD','FOR_GOOD','REV_GOOD','BOTH_GOOD') "
      . ");";
    $rtn = $dbh->do($sql); 
  } # only add when not added yet

  $exist = MY_DButility->check_table_existence($dbh, $dbname, $result_tables->{TRACE_TABLE});
  if($exist == 0) { # not exists
    $sql = "create table $result_tables->{TRACE_TABLE} ( "
      . "TRACE_ID bigint(20), "
      . "TRACE_PLA int(6), "
      . "TRACE_POS char(3), "
      . "TRACE_DIR enum('FOR', 'REV'), "
      . "TRACE_SEQ text, "
      . "TRACE_QUAL blob, "
      . "TRACE_AVG_QUAL float, "
      . "BL2_SCORE int(5), "
      . "INSERT_START int(6), "
      . "INSERT_END int(6), "
      . "PRIMARY KEY(TRACE_ID) "
      . ");";
    $rtn = $dbh->do($sql); 
  }

}

##########################################################################
# clean up the database for a clean process - has error, and redo later
# reset the following tables:
# xxx_TRACE_PLATEMAP,
# xxx_CHERRYMAP
# xxx_TRACE_SEQ
##########################################################################
sub reset_result_tables {
  my ($pkg, $dbh,$trace_plates, $result_tables, $fields)=@_;
 
  my ($sql, $query, $rtn);
#  my $filter = "$fields->{RUN_TRACE_PLA} in ( ". join(',', @{$trace_plates}) . ")";

# remove all the records in xxx_TRACE_PLATEMAP assoc. w/ @{$trace_plates}
  $sql = "delete from $result_tables->{PLATE_MAP_TABLE} "
       . " where TRACE_PLA in (" . join(',', @{$trace_plates}) . ")" 
       . ";";
  $rtn = $dbh->do($sql); 

# remove all the results for the seq/primer map
  $sql = "update $result_tables->{CHERRYMAP_TABLE} set "
       . "$fields->{RUN_TRACE_PLA}=null, "
       . "$fields->{RUN_TRACE_POS}=null, "
       . "$fields->{RUN_BL2_RESULT}=null "
       . "where $fields->{RUN_TRACE_PLA} in (" . join(',', @{$trace_plates}) . ")" 
       . ";";
  $rtn = $dbh->do($sql); 

# remove all the results for the trace seq info
  $sql = "delete from $result_tables->{TRACE_TABLE} "
       . " where TRACE_PLA in (" . join(',', @{$trace_plates}) . ")" 
         . ";";
  $rtn = $dbh->do($sql); 
  
  return 0;
}



####################################################################
# create the mapping between trace sequence plates and primer plates
####################################################################
sub load_plate_map {
    my ($pkg, $dbh, $dbh_trace, $proj_info, $trace_plates, $result_tables, $fields, $note)=@_; #input trace seq plate ids, i.e. (9144..9151);
    my ($sql, $query, $rtn);

    my $seq_type = $proj_info->{SEQ_TYPE};
    my $proj_name = $proj_info->{PROJ_ABBR};

    $sql = "insert IGNORE into $result_tables->{PLATE_MAP_TABLE} \
             (PROJ_ABBR, TRACE_PLA, TEMPL_PLA_NAME, SEQ_TYPE)\
              values(?, ?, ?, ?);";
     my $query_i = $dbh->prepare($sql);

     $sql = "select PLATE_ID, PLATE_NAME from tbPlate where plate_id in ("
          . join(",", @{$trace_plates})
          . ")";
     $query=$dbh_trace->prepare($sql);
     $query->execute();

     $sql = "update $result_tables->{PLATE_MAP_TABLE} set note= ? where TRACE_PLA=?";
     my $query_add_note = $dbh->prepare($sql);

     my $count_plates=0;
     while (my ($trace_plate, $templ_pla_name)=$query->fetchrow()) {

=pod
     $sql = "select distinct PLATE_NAME from " 
          . $result_tables->{CHERRYMAP_TABLE} . " a inner join "
          . $result_tables->{PLATE_DETAIL_TABLE} . " b "
          . "on a." . $fields->{RUN_TEMPL_PLA} . "=b.PLATE_ID "
          . "or a." . $fields->{RUN_TEMPL_PLA} . "=b.PLATE_NAME "
          . "where b.PLATE_ID=? or b.PLATE_NAME=?";
     my $query_templ_pl = $dbh->prepare($sql);

# check conso_plate, before insert into the map; - for there could
# be mistake between actual conso plate in cherry pick map and the
# one manually labeled
    	$query_templ_pl->execute($trace_plate, $templ_pla_name);
        my $templ_pla = $query_templ_pl->fetchrow();
        if ( !defined $templ_pla ) {
	  print "warning -- no plate with name = $templ_pla_name, double check please ...\n"; # probably not yet added to cherry map file
#	  next;
        }
=cut
        $query_i->execute($proj_name, $trace_plate, $templ_pla_name, $seq_type);
        if ( defined $note) {
	  $query_add_note->execute($note, $trace_plate);
        }
        $count_plates++;
   }
   return $count_plates;
}


############################################################################
# flag the rotated plates and map the correct positions
############################################################################

sub flag_nonstd_map {
    my ($pkg, $dbh, $trace_plate, $nonstd_type, $result_tables) = @_;
    my ($sql, $query, $rtn);
    my $nstd_map_type_table = 'yun_template_table.G_NONSTD_PLATEMAP_TYPE';

    $sql = "select id from $nstd_map_type_table where type=?";
    my $query_id = $dbh->prepare($sql);
    $query_id->execute($nonstd_type);
    my $id = $query_id->fetchrow();
    if( !defined $id)  { # if it is a new type
        my $detail;
        switch(uc($nonstd_type)) {
                case "ROTATE" {
		   $detail = "plate rotation happened in either cloning or sequencing stage";  }
                else {
                   print "no such type supported, please check and redo.\n";
                   return;
                }
        }
        $sql = "insert into $nstd_map_type_table (type, detail_desc) values('"
              . $nonstd_type . "', '". $detail . "');";
        $query = $dbh->prepare($sql);
        $query->execute();

#fetch newly inserted type id:
        $query_id->execute($nonstd_type);
        $id = $query_id->fetchrow();
    }
    $sql = "update $result_tables->{PLATE_MAP_TABLE} set is_std_map=?, non_std_map_type=? where trace_pla=?";
    my $query_u = $dbh->prepare($sql);
    $query_u->execute('NO', $id, $trace_plate);
}


############################################################################
# flag the rotated plates and map the correct positions
############################################################################
sub fix_rotate_plate {
    my ($pkg, $dbh, $plates_to_rotate, $result_tables, $fields) = @_;
    my @rows = ('A'..'H');
    my @cols = (1..12);
    my ($sql, $query, $rtn);

    $sql = "update $result_tables->{CHERRYMAP_TABLE} a "
         . "inner join $result_tables->{PLATE_DETAIL_TABLE} b on a.$fields->{RUN_TEMPL_PLA}=b.PLATE_ID "
         . "inner join $result_tables->{PLATE_MAP_TABLE} c on b.PLATE_NAME=c.TEMPL_PLA_NAME "
         . "set "
         . "$fields->{RUN_TRACE_PLA}= c.TRACE_PLA, "
         . "$fields->{RUN_TRACE_POS}= ? "
         . "where c.TRACE_PLA=? "
         . "and a.$fields->{RUN_TEMPL_POS}=? "
         . ";";
    my $query_u = $dbh->prepare($sql);

    my ($trace_pos, $templ_pos);
    my $count_pos_rotated = 0;
    foreach my $trace_pla (@{$plates_to_rotate}) {
      foreach my $r  (0..scalar(@rows)-1) {
        foreach my $c (0..scalar(@cols)-1) {
          $trace_pos = sprintf("%s%02s", uc($rows[$r]), $cols[$c]);
          $templ_pos = sprintf("%s%02s",uc($rows[-$r-1]), $cols[-$c-1]);
          $query_u->execute($trace_pos, $trace_pla, $templ_pos);
          $count_pos_rotated++;
        }
      }
    }
    return $count_pos_rotated;
}


#########################################################################
# after check non-std seq plate map, do this primer/seq map for later use
# load seq and primer map, to get the ref. seq
#########################################################################
sub load_seq_primer_map {
    my  ($pkg, $dbh, $trace_plates, $result_tables, $fields) = @_;
    my ($sql, $query, $rtn);

    my $nstd_map_type_table = 'yun_template_table.G_NONSTD_PLATEMAP_TYPE';

# add std plate/pos
    $sql = "update $result_tables->{CHERRYMAP_TABLE} a, "
         . "$result_tables->{PLATE_MAP_TABLE} b, "
         . "$result_tables->{PLATE_DETAIL_TABLE} c set "
         . "a.$fields->{RUN_TRACE_PLA}=b.TRACE_PLA, "
         . "a.$fields->{RUN_TRACE_POS}=(select POS from yun_template_table.G_POS_ID_MAP where POS_ID=a.$fields->{RUN_TEMPL_POSID}) "
         . "where b.TRACE_PLA in (" . join(",", @{$trace_plates}) . ") "
#         . "and c.PLATE_ID=a.$fields->{RUN_TEMPL_PLA} "
#         . "and c.PLATE_NAME=b.TEMPL_PLA_NAME "
         . "and (c.PLATE_NAME=a.$fields->{RUN_TEMPL_PLA} or c.PLATE_ID=a.$fields->{RUN_TEMPL_PLA}) " # eiterh use plate name or plate id to map cherry table and plate detail table
         . "and c.PLATE_NAME=b.TEMPL_PLA_NAME "
         . "and is_std_map='YES'"
         . ";";
     $query = $dbh->prepare($sql);
     $query->execute();
     $query->finish;

# add rotated plate/pos
    my @rotated_plates=();
    $sql = "select TRACE_PLA from "
         . "$result_tables->{PLATE_MAP_TABLE} a "
         . "inner join $nstd_map_type_table b "
         . "on a.NONSTD_MAP_TYPE=b.ID "
         . "where a.TRACE_PLA in (" . join(',', @{$trace_plates}) . ") "
         . "and a.is_std_map='NO' "
         . "and b.TYPE='ROTATE' "
         . ";";
    my $query_plates = $dbh->prepare($sql);
    $query_plates->execute();
    while (my ($trace_pla) = $query_plates->fetchrow()) {
      push @rotated_plates, $trace_pla;
    }

    if(scalar(@rotated_plates)>0) {
#     MY_SeqAlign->fix_rotate_plate($dbh, \@rotated_plates, $result_tables, $fields);
# the rotated plate, satisfies the following rule:
# original_pos_id+rotated_pos_id=97, => rotated_pos_id=97-original_pos_id
# the code following is based upon the above observation:
# add rotated plate/pos
     
       $sql = "update $result_tables->{CHERRYMAP_TABLE} a, "
         . "$result_tables->{PLATE_MAP_TABLE} b, "
         . "$result_tables->{PLATE_DETAIL_TABLE} c set "
         . "a.$fields->{RUN_TRACE_PLA}=b.TRACE_PLA, "
         . "a.$fields->{RUN_TRACE_POS}=(select POS from yun_template_table.G_POS_ID_MAP where POS_ID=97-a.$fields->{RUN_TEMPL_POSID}) "
         . "where b.TRACE_PLA in (" . join(",", @rotated_plates) . ") "
#         . "and c.PLATE_ID=a.$fields->{RUN_TEMPL_PLA} "
#         . "and c.PLATE_NAME=b.TEMPL_PLA_NAME "
         . "and (c.PLATE_NAME=a.$fields->{RUN_TEMPL_PLA} or c.PLATE_ID=a.$fields->{RUN_TEMPL_PLA}) " # eiterh use plate name or plate id to map cherry table and plate detail table
         . "and c.PLATE_NAME=b.TEMPL_PLA_NAME "
         . ";";
       $query = $dbh->prepare($sql);
       $query->execute();
       $query->finish;
     }


}


##########################################################################
# load traces from bali2 to the local database
##########################################################################
sub load_trace_seq {
    my ($pkg, $dbh, $dbh_trace, $trace_plates, $result_tables, $fields) = @_; # get all plate ids, in which seq's to be aligned
    my ($sql, $query, $rtn);

    $sql = "select TRACE_ID, PLATE_ID as TRACE_PLA, substr(trace_name, -7,3) as TRACE_POS, FASTA, QUAL_DATA, AVG_QUAL_VAL as TRACE_AVG_QUAL \
            from tbTrace where PLATE_ID in (" . join(',',@{$trace_plates}) . ")";
    $query = $dbh_trace->prepare($sql);
    $query->execute();
    $sql = "insert IGNORE into $result_tables->{TRACE_TABLE} "
         . "(TRACE_ID, TRACE_PLA, TRACE_POS, TRACE_DIR, TRACE_SEQ, TRACE_QUAL, TRACE_AVG_QUAL) "
         . "values(" . "?,"x6 . "?)"
         . ";";
    my $query_i = $dbh->prepare($sql);

    $sql = "select count(*) from $result_tables->{CHERRYMAP_TABLE} where $fields->{RUN_TRACE_PLA}=? and $fields->{RUN_TRACE_POS}=?";
    my $query_valid = $dbh->prepare($sql);

    my $count_trace_seq=0;
    while(my ($trace_id, $trace_pla, $trace_pos, $fasta, $qual_data, $avg_qual_val)=$query->fetchrow()) {
       $fasta =~ s/\n//g;
       $qual_data =~ s/\n//g;
       my ($dir, $pos, $trace_seq) = $fasta =~ /\w+\d+.+(FOR|REV)_(\w\d+)\.scf(.+)/i;
       my ($qual) = $qual_data =~ />.+SCF(.+)$/;
       $query_valid->execute($trace_pla, $trace_pos);
       my $count = $query_valid->fetchrow();
       if($count>0) {
	 $query_i->execute($trace_id, $trace_pla, $trace_pos, $dir, $trace_seq, $qual, $avg_qual_val);
         $count_trace_seq++;
       }
    }
    return $count_trace_seq;
}


##########################################################################
# load traces from bali2 to the local database with no mapping info
##########################################################################
sub load_trace_seq_nomap {
    my ($pkg, $t, $trace_plates) = @_; # get all plate ids, in which seq's to be aligned
    my ($sql, $query, $rtn);

    $sql = "select TRACE_ID, TRACE_NAME, b.plate_name as trace_plate_name,  a.PLATE_ID as TRACE_PLATE, substr(trace_name, -7,3) as TRACE_POS, FASTA, QUAL_DATA, AVG_QUAL_VAL as TRACE_AVG_QUAL, IDENTITY as TRACE_IDENTITY from $t->{TRACE_ORIG}{TABLENAME} a inner join tbPlate b on a.plate_id=b.plate_id where a.PLATE_ID in (" . join(',',@{$trace_plates}) . ")";# . " and trace_id=1289057";
    $query = $t->{TRACE_ORIG}{DBH}->prepare($sql);
    $rtn = $query->execute();
    
    $sql = "insert IGNORE into $t->{TRACE_SEQ}{TABLENAME} "
         . "(TRACE_ID, TRACE_NAME, TRACE_PLATE_NAME, TRACE_PLATE, TRACE_POS, TRACE_DIR, TRACE_SEQ, TRACE_QUAL, TRACE_AVG_QUAL, TRACE_IDENTITY) "
         . "values(" . "?,"x9 . "?)"
         . ";";
    my $query_i = $t->{TRACE_SEQ}{DBH}->prepare($sql);

    my $count_trace_seq=0;
    while(my ($trace_id, $trace_name, $trace_plate_name, $trace_pla, $trace_pos, $fasta, $qual_data, $avg_qual_val, $trace_identity)=$query->fetchrow()) {
       $fasta =~ s/\n//g;
       $fasta =~ s/FWD|Forward/FOR/i;
       $fasta =~ s/Reverse/REV/i;
       $qual_data =~ s/\n//g;
#       my ($dir, $pos, $trace_seq) = $fasta =~ /\w+\d+.+(FOR|REV)_(\w\d+)\.scf(.+)/i;
#       my ($dir, $pos, $trace_seq) = $fasta =~ /M13.*(FOR|REV).*_(\w\d+)\.scf(.+)/i;
       my ($dir, $pos, $trace_seq) = $fasta =~ />.+(FOR|REV)_(\w\d+)\.scf(.+)/i;
       my ($qual) = $qual_data =~ />.+SCF(.+)$/;
       $rtn = $query_i->execute($trace_id, $trace_name, $trace_plate_name, $trace_pla, $trace_pos, $dir, $trace_seq, $qual, $avg_qual_val, $trace_identity);
         $count_trace_seq++;
    }
    return $count_trace_seq;
} # load_trace_seq_nomap



##########################################################################
# clean up seq_align database for a clean process 
##########################################################################
sub reset_result_tables_seq_align {
  my ($pkg, $dbh, $dbh_align, $trace_plates, $result_tables, $fields) = @_;
  my ($sql, $query, $rtn);
  my $filter = "PLATE_ID in ( ". join(',', @{$trace_plates}) . ")";

# clean OST_SEQ
  $sql = "delete a.* from $result_tables->{OST_SEQ_TABLE} a "
         . "where  $filter "
         . ";";
  $rtn = $dbh_align->do($sql); 
 
  $sql = "select TEMPL_PLA_NAME from $result_tables->{PLATE_MAP_TABLE} where TRACE_PLA in (" . join(',', @{$trace_plates}) . ")";
  $query = $dbh->prepare($sql);
  $query->execute();
  my @trace_plate_names=();
  while (my $trace_pla_name = $query->fetchrow()) {
     push @trace_plate_names, $trace_pla_name;
  }
  
# clean ALIGN_SEQ table
  $sql = "delete a.* from $result_tables->{ALIGN_SEQ_TABLE} a "
       . "inner join $result_tables->{SEQ_MAP_TABLE} b on a.TRACE_ID=b.TRACE_ID "
       . "inner join $result_tables->{EXP_SEQ_TABLE} c on b.EXP_SEQ_ID=c.SEQ_ID "
       . "where c.EPLATE_NAME in ('" . join("','", @trace_plate_names) . "')"
       . ";";
  $rtn = $dbh_align->do($sql); 

# clean SEQ_MAP table
  $sql = "delete b.* from $result_tables->{SEQ_MAP_TABLE} b  "
       . "inner join $result_tables->{EXP_SEQ_TABLE} c on b.EXP_SEQ_ID=c.SEQ_ID "
       . "where c.EPLATE_NAME in ('" . join("','", @trace_plate_names) . "')"
       . ";";
  $rtn = $dbh_align->do($sql); 

# clean EXP_SEQ table
  $sql = "delete c.* from $result_tables->{EXP_SEQ_TABLE} c "
       . "where c.EPLATE_NAME in ('" . join("','", @trace_plate_names) . "')"
       . ";";
  $rtn = $dbh_align->do($sql); 

}



##########################################################################
# load refseq into seq_align db
##########################################################################
sub load_refseq_map_seq_align {
  my ($pkg, $dbh, $dbh_align, $trace_plates, $result_tables, $fields, $proj_info) = @_;
  my ($sql, $query, $rtn);
  my $filter = "PLATE_ID in ( ". join(',', @{$trace_plates}) . ")";

# get refseq
  $sql = "select a.ORF_ID, a.$fields->{RUN_TRACE_PLA}, a.$fields->{RUN_TRACE_POS}, b.CDS_SEQ, c.PLATE_NAME "
      .  "from $result_tables->{CHERRYMAP_TABLE} a "
      .  "inner join $result_tables->{REFSEQ_TABLE} b using(ORF_ID) "
	."inner join $result_tables->{PLATE_DETAIL_TABLE} c on a.$fields->{RUN_TEMPL_PLA}=c.PLATE_ID "
      .  "where a.$fields->{RUN_TRACE_PLA} in (" . join(',',@{$trace_plates}) . ")"
      . ";";
  $query = $dbh -> prepare($sql);

   $sql = "insert IGNORE into EXP_SEQ (EX_ID, EPLATE_NAME, PLATE_POS, NOTE, ORGANISM,SEQ_TYPE, EXP_SEQ) values (" . "?,"x6 . "?) ";
   my $query_add_expseq = $dbh_align->prepare($sql);

   $sql = "insert IGNORE into SEQ_MAP (EXP_SEQ_ID, TRACE_ID, NOTE) values (" . "?,"x2 . "?) ";
   my $query_add_seqmap = $dbh_align->prepare($sql);

   $sql = "select last_insert_id();";
   my $query_lastid = $dbh_align->prepare($sql);

   $sql = "select trace_id from $result_tables->{TRACE_TABLE} where TRACE_PLA=? and TRACE_POS=?";
   my $query_trace = $dbh->prepare($sql);

   $query->execute();
   while (my ($orf_id, $trace_pla, $trace_pos, $refseq, $templ_pla_name) = $query->fetchrow) {
     $query_add_expseq->execute($proj_info->{PROJ_ABBR} . $orf_id, 
                                $templ_pla_name, 
                                $trace_pos, 
                                "",
                                $proj_info->{ORGANISM}, 
                                $proj_info->{SEQ_TYPE}, 
                                $refseq); 
     $query_lastid->execute();
     my $expseq_id = $query_lastid->fetchrow();
     $query_trace->execute($trace_pla, $trace_pos);
     while(my $trace_id=$query_trace->fetchrow()) {
        $query_add_seqmap->execute($expseq_id, $trace_id, "");           
     }
      
   } # end of all traces loop
}





#*************************************************************************
#* END of subroutines
#*************************************************************************



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$ BEGIN of entry functions
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

##########################################################################
# map the seq. and load map to db so as to get ready for later alignment
##########################################################################
sub preprocess_seq {
  my ($pkg, $dbh, $dbname, $dbh_trace, $trace_plates, $result_tables, $fields, $nonstd_plates, $proj_info) = @_;
  my $seq_type = $proj_info->{SEQ_TYPE};

# clean up thoroughly
  MY_SeqAlign->reset_result_tables($dbh, $trace_plates, $result_tables, $fields);

# load all the prep. data and trace seqs into human orfdb
  MY_SeqAlign->load_plate_map($dbh, $dbh_trace, $proj_info, $trace_plates, $result_tables, $fields);

  # below is not routine, only when plate rotation or etc. happens
  if(defined $nonstd_plates ) {
    foreach my $nonstd_plate (@{$nonstd_plates}) {
      flag_nonstd_map($nonstd_plate->{TRACE_PLA}, $nonstd_plate->{NONSTD});
      if(uc($nonstd_plate->{NONSTD}) eq 'ROTATE') {
	MY_SeqAlign->fix_rotate_plate($dbh, $nonstd_plate->{PLATE_ID});
      }
    }
  }
  MY_SeqAlign->load_seq_primer_map($dbh, $trace_plates, $result_tables, $fields);
  MY_SeqAlign->load_trace_seq($dbh, $dbh_trace, $trace_plates, $result_tables, $fields);
}



####################################################
# dealing with the experiment result
# INPUT:
#   @trace_plates = (10805..10814);
#   %result_tables = (
#     REFSEQ_TABLE => 'horfeome7.cDNA',
#     TRACE_TABLE => 'isoform.SANGER_TRACE_SEQ',
#     CHERRYMAP_TABLE => 'isoform.ISOFORM_CHERRYMAP',
#   );
#
#   %fields = (
#     RUN_TRACE_PLA =>'RUN1_TRACE_PLA',
#     RUN_TRACE_POS =>'RUN1_TRACE_POS',
#   );
#
####################################################
sub run_bl2seq {
    my ($pkg, $dbh, $trace_plates, $result_tables, $fields) = @_;

    my $count_hi_score_seq = 0;
    my $count_total_seq = 0;
    my $HSP_THRESHOLD = 60; # originally it set to be 100, and the result in DB are based on 100 setting. 
    my $refseq_table = $result_tables->{REFSEQ_TABLE};
    my ($sql, $query, $rtn);

    $sql = "select distinct a.TRACE_ID, b.ORF_ID, a.TRACE_PLA, a.TRACE_POS, TRACE_SEQ, CDS_SEQ \
            from $result_tables->{TRACE_TABLE} a inner join $result_tables->{CHERRYMAP_TABLE} b \
            on a.TRACE_PLA=b." . $fields->{RUN_TRACE_PLA} . " and a.TRACE_POS=b." . $fields->{RUN_TRACE_POS} . " \  
            inner join $refseq_table c using(ORF_ID) \
            where a.TRACE_PLA in (" . join(',',@{$trace_plates}) . ") "
#            . "and a.TRACE_ID=992750"
#	    . "and a.TRACE_POS = 'A11'" 
            . ";";
    $query = $dbh->prepare($sql);
    $query->execute();
    
    $sql = "update $result_tables->{TRACE_TABLE} set BL2_SCORE=?, INSERT_START=?, INSERT_END=? where TRACE_ID=?";
    my $update_q = $dbh->prepare($sql);
    my $factory = Bio::Tools::Run::StandAloneBlast->new(
                'program' => 'blastn',
                'F'=>'F', #must set filter to be false to avoid missing the longest alignment b/c low complexity filter 
                 );
    while(my ($trace_id, $orf_id, $trace_pla, $trace_pos, $trace_seq, $cds_seq)=$query->fetchrow_array()){
      #print "$trace_id\n";
      $count_total_seq++;
      my $traceobj= Bio::Seq->new( -display_id => $trace_id,
                                   -seq => $trace_seq);
      my $cdsobj=Bio::Seq->new( -display_id => $orf_id,
                                -seq => $cds_seq);
      $factory->io->_io_cleanup();
      my $report = $factory->bl2seq($traceobj, $cdsobj);
      my $hi_score = 0;
      my $insert_start=-1;
      my $insert_end=-1;
      if ($report) {
        while (my $res=$report->next_result){
          while (my $hit= $res->next_hit){
            while(my $hsp=$hit->next_hsp) {
              next if ($hsp->length('hit')<$HSP_THRESHOLD);
              #print "score ", $hsp->score," Percent ", $hsp->frac_identical('hit'), " Match ",$hsp->length('hit')," Length ",  $hsp->length('total'),"\n";
              #print "$trace_id\t", $hsp->score, "\n";
              if($hsp->score > $hi_score) {
                 $hi_score = $hsp->score;
		 $insert_start = $hsp->hit->start;
		 $insert_end = $hsp->hit->end;
	      }
            }
          }
        }
        if ($hi_score != 0){
#          $update_q->execute($hi_score, $insert_start, $insert_end, $trace_id);
          $update_q->execute($hi_score, $insert_start-1, $insert_end-1, $trace_id); # the bl2 result seems 1 base shift
          $count_hi_score_seq++;
        }
      }
    }
printf "total seq. # : $count_total_seq, high score seq # : $count_hi_score_seq\n";   

    return $count_hi_score_seq;

}

##########################################################################
# general bl2seq code when the source doesn't match with above standard format
##########################################################################
sub  calc_phred_score {
    my ($pkg, $tables) = @_;
    my ($sql, $query, $rtn);
 
    if(!defined($tables->{TRACE_SEQ}{TRACE_ID_FIELD})) {
        $tables->{TRACE_SEQ}{TRACE_ID_FIELD} = 'TRACE_ID';
    }
    if(!defined($tables->{TRACE_SEQ}{TRACE_PLATE_FIELD})) {
        $tables->{TRACE_SEQ}{TRACE_PLATE_FIELD} = 'TRACE_PLATE';
    }
    if(!defined($tables->{TRACE_SEQ}{TRACE_POS_FIELD})) {
        $tables->{TRACE_SEQ}{TRACE_POS_FIELD} = 'TRACE_POS';
    }
    if(!defined($tables->{TRACE_SEQ}{INSERT_START_FIELD})) {
        $tables->{TRACE_SEQ}{INSERT_START_FIELD} = 'VTRIM_START';
    }
    if(!defined($tables->{TRACE_SEQ}{INSERT_END_FIELD})) {
        $tables->{TRACE_SEQ}{INSERT_END_FIELD} = 'VTRIM_END';
    }
    if(!defined($tables->{TRACE_SEQ}{TRACE_QUAL_FIELD})) {
        $tables->{TRACE_SEQ}{TRACE_QUAL_FIELD} = 'TRACE_QUAL'; 
    }
   
    $sql = "select $tables->{TRACE_SEQ}{TRACE_ID_FIELD}, $tables->{TRACE_SEQ}{TRACE_PLATE_FIELD},$tables->{TRACE_SEQ}{TRACE_POS_FIELD}, $tables->{TRACE_SEQ}{INSERT_START_FIELD},$tables->{TRACE_SEQ}{INSERT_END_FIELD}, $tables->{TRACE_SEQ}{TRACE_QUAL_FIELD} from $tables->{TRACE_SEQ}{TABLENAME} ";
    $sql .= " where $tables->{FILTER}{PHRED} " if defined($tables->{FILTER}{PHRED});
    $sql .= " order by $tables->{ORDER_BY}{PHRED} " if defined($tables->{ORDER_BY}{PHRED});
    $query = $tables->{TRACE_SEQ}{DBH} -> prepare($sql);
    $sql = "update $tables->{TRACE_SEQ}{TABLENAME} set $tables->{TRACE_SEQ}{AVG_PHRED_FIELD} = ? where $tables->{TRACE_SEQ}{TRACE_ID_FIELD} = ?";
    my $qu_phred = $tables->{TRACE_SEQ}{DBH} -> prepare($sql);
    $rtn = $query -> execute();
    while (my ($TRACE_ID, $TRACE_PLATE, $TRACE_POS, $INSERT_START, $INSERT_END, $TRACE_QUAL) = $query->fetchrow()) {
      my @qual_array = split('\s', $TRACE_QUAL);
      my $qual_total = 0;
      $INSERT_START = defined($INSERT_START)?$INSERT_START:1;
      $INSERT_END = defined($INSERT_END)?$INSERT_END:scalar(@qual_array);
      for(my $i=$INSERT_START-1; $i<$INSERT_END; $i++) {
        $qual_total += $qual_array[$i];
      } #end of foreach
      $qu_phred -> execute($qual_total*1.0/($INSERT_END-$INSERT_START+1), $TRACE_ID);
    } # end of while  
    $query -> finish;
    $qu_phred -> finish;
} # end of calc_phred_sore();

##########################################################################
# general bl2seq code when the source doesn't match with above standard format
##########################################################################
sub general_bl2seq {
    my ($pkg, $dbh, $tables, $bl2seq_type) = @_;
    my $HSP_THRESHOLD = 60; # originally it set to be 100, and the result in DB are based on 100 setting. 
    my ($sql, $query, $rtn);
    
    my $trace_table = $tables->{TRACE_SEQ};
    my $refseq_table = $tables->{REFSEQ};
    my $cherrymap_table = $tables->{CHERRYMAP};

   if(!defined($tables->{TRACE_SEQ}{BL2_SCORE_FIELD})) {
        $tables->{TRACE_SEQ}{BL2_SCORE_FIELD} = 'BL2_SCORE';
    }
    if(!defined($trace_table->{INSERT_START_FIELD})) {
        $trace_table->{INSERT_START_FIELD} = 'INSERT_START';
    }
    if(!defined($trace_table->{INSERT_END_FIELD})) {
        $trace_table->{INSERT_END_FIELD} = 'INSERT_END';
    }
    if(!defined($trace_table->{ALIGN_LEN_FIELD})) {
        $trace_table->{ALIGN_LEN_FIELD} = 'ALIGN_LEN';
    }
    if(!defined($trace_table->{IDENTITY_FIELD})) {
        $trace_table->{IDENTITY_FIELD} = 'TRACE_IDENTITY';
    }
    if(!defined($trace_table->{TRACE_ID_FIELD})) {
        $trace_table->{TRACE_ID_FIELD} = 'TRACE_ID';
    }
    if(!defined($trace_table->{SCORE_FIELD})) {
        $trace_table->{SCORE_FIELD} = 'SCORE';
    }
    if(!defined($trace_table->{BITS_FIELD})) {
        $trace_table->{BITS_FIELD} = 'BITS';
    }
    if(!defined($trace_table->{numHSP_FIELD})) {
	$trace_table->{numHSP_FIELD} = 'numHSP';
    }
    if(!defined($trace_table->{qLenAlign_FIELD})) {
        $trace_table->{qLenAlign_FIELD} = 'qLenAlign';
    }
    if(!defined($trace_table->{sLenAlign_FIELD})) {
        $trace_table->{sLenAlign_FIELD} = 'sLenAlign';
    }
    if(!defined($trace_table->{qStart_FIELD})) {
        $trace_table->{qStart_FIELD} = 'qStart';
    }
    if(!defined($trace_table->{qEnd_FIELD})) {
        $trace_table->{qEnd_FIELD} = 'qEnd';
    }
    if(!defined($trace_table->{sStart_FIELD})) {
        $trace_table->{sStart_FIELD} = 'sStart';
    }
    if(!defined($trace_table->{sEnd_FIELD})) {
        $trace_table->{sEnd_FIELD} = 'sEnd';
    }
    if(!defined($trace_table->{qNumUnaligned_FIELD})) {
        $trace_table->{qNumUnaligned_FIELD} = 'qNumUnaligned';
    }
    if(!defined($trace_table->{sNumUnaligned_FIELD})) {
        $trace_table->{sNumUnaligned_FIELD} = 'sNumUnaligned';
    }
    if(!defined($trace_table->{qSeqIndIdentical_FIELD})) {
        $trace_table->{qSeqIndIdentical_FIELD} = 'qSeqIndIdentical';
    }
    if(!defined($trace_table->{sSeqIndIdentical_FIELD})) {
         $trace_table->{sSeqIndIdentical_FIELD} = 'sSeqIndIdentical';
    }
    if(!defined($trace_table->{qStrand_FIELD})) {
         $trace_table->{qStrand_FIELD} = 'qStrand';
    }
    if(!defined($trace_table->{sStrand_FIELD})) {
         $trace_table->{sStrand_FIELD} = 'sStrand';
    }
    if(!defined($trace_table->{RANK_FIELD})) {
	$trace_table->{RANK_FIELD} = 'RANK';
    }
    if(!defined($trace_table->{qSeqLenIdentical_FIELD})) {
         $trace_table->{qSeqLenIdentical_FIELD} = 'qSeqLenIdentical';
    }
    if(!defined($trace_table->{sSeqLenIdentical_FIELD})) {
         $trace_table->{sSeqLenIdentical_FIELD} = 'sSeqLenIdentical';
    }
    if(!defined($trace_table->{qSeqGap_FIELD})) {
	$trace_table->{qSeqGap_FIELD} = 'qSeqGap';
    }
    if(!defined($trace_table->{sSeqGap_FIELD})) {
	$trace_table->{sSeqGap_FIELD} = 'sSeqGap';
    }
    if(!defined($tables->{REFSEQ}{REFSEQ_ID_FIELD})) {
        $tables->{REFSEQ}{REFSEQ_ID_FIELD} = 'ORF_ID';
    }
    if(!defined($tables->{REFSEQ}{REFSEQ_FIELD})) {
        $tables->{REFSEQ}{REFSEQ_FIELD} = 'CDS_SEQ';
    }

    $sql = "update $trace_table->{TABLENAME} set $trace_table->{BL2_SCORE_FIELD}=?, $trace_table->{INSERT_START_FIELD}=$trace_table->{SEQ_OFFSET_FIELD}-1+?, $trace_table->{INSERT_END_FIELD}=$trace_table->{SEQ_OFFSET_FIELD}-1+? , $trace_table->{ALIGN_LEN_FIELD}=?, $trace_table->{IDENTITY_FIELD}=? where $trace_table->{TRACE_ID_FIELD}=?";
    my $update_q = $dbh->prepare($sql);

    my $update_q_blast;

#### to check on Monday below########

    if($bl2seq_type == 1) {
      $sql = "update $trace_table->{TABLENAME} set $trace_table->{SCORE_FIELD}=?, $trace_table->{BITS_FIELD}=?, $trace_table->{numHSP_FIELD}=?, $trace_table->{qLenAlign_FIELD}=?, $trace_table->{sLenAlign_FIELD}=?, $trace_table->{qStart_FIELD}=?, $trace_table->{qEnd_FIELD}=?, $trace_table->{sStart_FIELD}=?, $trace_table->{sEnd_FIELD}=?, $trace_table->{qNumUnaligned_FIELD}=?, $trace_table->{sNumUnaligned_FIELD}=?, $trace_table->{qSeqIndIdentical_FIELD}=?,$trace_table->{sSeqIndIdentical_FIELD}=?, $trace_table->{qStrand_FIELD}=?, $trace_table->{sStrand_FIELD}=?, $trace_table->{RANK_FIELD}=?, $trace_table->{qSeqLenIdentical_FIELD}=?, $trace_table->{sSeqLenIdentical_FIELD}=?, $trace_table->{qSeqGap_FIELD}=?, $trace_table->{sSeqGap_FIELD}=? where $trace_table->{TRACE_ID_FIELD}=?";
      $update_q_blast = $dbh->prepare($sql);
    };

    my $factory = Bio::Tools::Run::StandAloneBlast->new(
                'program' => 'blastn',
                'F'=>'F', #must set filter to be false to avoid missing the longest alignment b/c low complexity filter 
                 'outfile' => 'bl2seq20141117.out',
                 );

=pod
#    $sql = "select distinct a.$trace_table->{TRACE_ID_FIELD}, a.$trace_table->{TRACE_PLATE_FIELD}, a.$trace_table->{TRACE_POS_FIELD}, $trace_table->{TRACE_SEQ_FIELD}, c.$refseq_table->{REFSEQ_ID_FIELD}, c.$refseq_table->{REFSEQ_FIELD} from (select * from $trace_table->{TABLENAME} ";
#    $sql .= " where " . $trace_table->{FILTER} if defined $trace_table->{FILTER};
#    $sql .= ") a inner join $cherrymap_table->{TABLENAME} b on $tables->{JOIN_CONDITION}{TRACE_CHERRY} inner join $refseq_table->{DBNAME}.$refseq_table->{TABLENAME} c on $tables->{JOIN_CONDITION}{CHERRY_REFSEQ} ";
#    $sql .= " where $tables->{FILTER}{BL2SEQ} " if defined($tables->{FILTER}{BL2SEQ});
#    $sql .= " order by $tables->{ORDER_BY}{BL2SEQ}" if defined($tables->{ORDER_BY}{BL2SEQ});
=cut
     $sql = "select distinct a." . $trace_table->{TRACE_ID_FIELD} . ", a." . $trace_table->{TRACE_PLATE_FIELD} . ", a." . $trace_table->{TRACE_POS_FIELD} .", " . $trace_table->{TRACE_SEQ_FIELD} . ", c." . $refseq_table->{REFSEQ_ID_FIELD} . ", c." . $refseq_table->{REFSEQ_FIELD} . " from $trace_table->{TABLENAME} ";
    $sql .= " a inner join $cherrymap_table->{TABLENAME} b on $tables->{JOIN_CONDITION}{TRACE_CHERRY} inner join $refseq_table->{DBNAME}.$refseq_table->{TABLENAME} c on $tables->{JOIN_CONDITION}{CHERRY_REFSEQ} ";
    $sql .= " where 1";
    $sql .= " and " . $tables->{FILTER}{BL2SEQ} if defined($tables->{FILTER}{BL2SEQ});
    $sql .= " and " . $trace_table->{FILTER} if defined $trace_table->{FILTER};
    $sql .= " order by $tables->{ORDER_BY}{BL2SEQ}" if defined($tables->{ORDER_BY}{BL2SEQ});
    $query = $dbh->prepare($sql);
    $rtn = $query -> execute();
    my $count_hi_score_seq = 0;
    my $count_total_seq = 0;
    while (my ($trace_id, $trace_plate, $trace_pos, $trace_seq, $refseq_id, $refseq) = $query->fetchrow()) {

print "now checking TRACE_ID:$trace_id\n";

      $count_total_seq++;
      if ($trace_seq eq "") { # skip empty/unqualified trace
         next;
      }
#print "$trace_id\t$trace_plate\t$trace_pos\n";

      my $traceobj= Bio::Seq->new( -display_id => $trace_id,
                                   -seq => $trace_seq);
      my $cdsobj=Bio::Seq->new( -display_id => $refseq_id,
                                -seq => $refseq);
      $factory->io->_io_cleanup();
      my $report = $factory->bl2seq($traceobj, $cdsobj);
      my $hi_score = 0;
      my $insert_start=-1;
      my $insert_end=-1;
      my $aln_length=-1;
      my $percent_identity = -1;
      my %col2insert = ();

      if ($report) {
        while (my $res=$report->next_result){
          my ($previous_qStart, $previous_qEnd) = (0,0);
          while (my $hit= $res->next_hit){
# first check if need to record bl2 details:
if($bl2seq_type == 1) {
            last if($hit->rank>2 || $hit->score<MIN_BL_SCORE || ($previous_qStart<=$hit->end('query') && $previous_qStart>=$hit->start('query') && $hit->end('query')-$previous_qStart>MIN_HIT_OVERLAP) || ($hit->start('query')>=$previous_qStart && $hit->start('query')<=$previous_qEnd && $previous_qEnd-$hit->start('query')>MIN_HIT_OVERLAP)); # discard rank 3 and above hits, and those score are too low, and hit overlaps with higher rank
            $count_hi_score_seq ++;
            $previous_qStart = $hit->start('query');
            $previous_qEnd = $hit->end('query');      
            $col2insert{score} = $hit->score;
            $col2insert{bits} = $hit->bits;
            $col2insert{numHSP} = $hit->num_hsps;
            $col2insert{qLenAlign} = $hit->length_aln('query');
            $col2insert{sLenAlign} = $hit->length_aln('sbjct');
            $col2insert{qStart} = $hit->start('query');
            $col2insert{qEnd} = $hit->end('query');
            $col2insert{sStart} = $hit->start('sbjct');
            $col2insert{sEnd} = $hit->end('sbjct');
            $col2insert{qNumUnaligned} = $hit->num_unaligned_query;
            $col2insert{sNumUnaligned} = $hit->num_unaligned_hit;
            $col2insert{qSeqIndIdentical} = join('::', $hit->seq_inds('query', 'identical', 1));
            $col2insert{sSeqIndIdentical} = join('::', $hit->seq_inds('sbjct', 'identical', 1));
            $col2insert{qStrand} = join('::', $hit->strand('query'));
            $col2insert{sStrand} = join('::', $hit->strand('sbjct'));
            $col2insert{rank} = join ('::', $hit->rank);

# filling gap info:

            my $qLenGap =  getAlignmentGap($col2insert{qSeqIndIdentical});
	    if (defined $qLenGap->{SEQLEN}) {
	      $col2insert{qSeqLenIdentical} = join(',',@{$qLenGap->{SEQLEN}});  
	    } else {
	      $col2insert{qSeqLenIdentical} = 0; 
	    }

	    if (defined $qLenGap->{SEQGAP}) {
	      $col2insert{qSeqGap} = join(',',@{$qLenGap->{SEQGAP}});  
	    } else {
	      $col2insert{qSeqGap} = 0;
	    }

	    my $sLenGap =  getAlignmentGap($col2insert{sSeqIndIdentical});
	    if (defined $sLenGap->{SEQLEN}) {
	      $col2insert{sSeqLenIdentical} = join(',',@{$sLenGap->{SEQLEN}});  
	    } else {
	      $col2insert{sSeqLenIdentical} = 0; 
	    }

	    if (defined $sLenGap->{SEQGAP}) {
	      $col2insert{sSeqGap} = join(',',@{$sLenGap->{SEQGAP}});  
	    } else {
	      $col2insert{sSeqGap} = 0;
	    }

} # end of $bl2seq_type

            while(my $hsp=$hit->next_hsp) {
              next if ($hsp->length('hit')<$HSP_THRESHOLD);
              #print "score ", $hsp->score," Percent ", $hsp->frac_identical('hit'), " Match ",$hsp->length('hit')," Length ",  $hsp->length('total'),"\n";
              #print "$trace_id\t", $hsp->score, "\n";
              if($hsp->score > $hi_score) {
                 $hi_score = $hsp->score;
#		 $insert_start = $hsp->query->start-1; # 1-base shift
#		 $insert_end = $hsp->query->end-1;
		 $insert_start = $hsp->query->start; # 1-base shift
		 $insert_end = $hsp->query->end;
                 $aln_length = $hsp->length('total');
                 $percent_identity = $hsp->percent_identity();
                 if($bl2seq_type==1) {
# fix the bug caused by bl2seq tiled_hsps():
                   if($col2insert{sSeqGap} == 0) { # get seq ind from first hsp:
                     $col2insert{qSeqIndIdentical} = join('::', $hsp->seq_inds('query', 'identical', 1));
                     $col2insert{sSeqIndIdentical} = join('::', $hsp->seq_inds('sbjct', 'identical', 1));
                     my $qLenGap =  getAlignmentGap($col2insert{qSeqIndIdentical});
	             if (defined $qLenGap->{SEQLEN}) {
	                $col2insert{qSeqLenIdentical} = join(',',@{$qLenGap->{SEQLEN}});  
	             } else {
	               $col2insert{qSeqLenIdentical} = 0; 
	             }
 
	             if (defined $qLenGap->{SEQGAP}) {
	               $col2insert{qSeqGap} = join(',',@{$qLenGap->{SEQGAP}});  
	             } else {
	               $col2insert{qSeqGap} = 0;
	             }

	             my $sLenGap =  getAlignmentGap($col2insert{sSeqIndIdentical});
	             if (defined $sLenGap->{SEQLEN}) {
	               $col2insert{sSeqLenIdentical} = join(',',@{$sLenGap->{SEQLEN}});  
	             } else {
	               $col2insert{sSeqLenIdentical} = 0; 
	             }

	             if (defined $sLenGap->{SEQGAP}) {
	               $col2insert{sSeqGap} = join(',',@{$sLenGap->{SEQGAP}});  
	             } else {
	               $col2insert{sSeqGap} = 0;
	             }
		 } # end of if $bl2seq_type==1
	      } # end if-4
            } # end while-3
          } # end while-2
        } # end while-1
        if ($hi_score != 0){
          $rtn = $update_q->execute($hi_score, $insert_start, $insert_end, $aln_length, $percent_identity, $trace_id);
          if($bl2seq_type==1) {
              $rtn = $update_q_blast -> execute($col2insert{score}, 
                                        $col2insert{bits},
                                        $col2insert{numHSP},
					$col2insert{qLenAlign},
					$col2insert{sLenAlign},
                                        $col2insert{qStart},
                                        $col2insert{qEnd},
                                        $col2insert{sStart},
                                        $col2insert{sEnd},
                                        $col2insert{qNumUnaligned},
                                        $col2insert{sNumUnaligned},
                                        $col2insert{qSeqIndIdentical},
                                        $col2insert{sSeqIndIdentical},
                                        $col2insert{qStrand},
                                        $col2insert{sStrand},
					$col2insert{rank},
                                        $col2insert{qSeqLenIdentical},
                                        $col2insert{sSeqLenIdentical},
					$col2insert{qSeqGap},
                                        $col2insert{sSeqGap},
                                        $trace_id);	
	     }
          }  
          $count_hi_score_seq ++;
        } # end if-2
      } # end if-1
    } # end of while-0
    $query -> finish;
    $update_q -> finish;
    if ($bl2seq_type==1) {
      $update_q_blast -> finish;
    }
    return $count_hi_score_seq;  


} # end of general_bl2seq



##########################################################################
# general bl2seq code recording  brief blast result
##########################################################################
sub general_bl2seq_brief_blast {
    my ($pkg, $dbh, $tables) = @_;
#    my $HSP_THRESHOLD = 60; # originally it set to be 100, and the result in DB are based on 100 setting. 
    my ($sql, $query, $rtn);
    
    my $trace_table = $tables->{TRACE_SEQ};
    my $refseq_table = $tables->{REFSEQ};
    my $cherrymap_table = $tables->{CHERRYMAP};

    $sql = "update $trace_table->{TABLENAME} set $trace_table->{SCORE_FIELD}=?, $trace_table->{BITS_FIELD}=?, $trace_table->{numHSP_FIELD}=?, $trace_table->{qLenAlign_FIELD}=?, $trace_table->{sLenAlign_FIELD}=?, $trace_table->{qStart_FIELD}=?, $trace_table->{qEnd_FIELD}=?, $trace_table->{sStart_FIELD}=?, $trace_table->{sEnd_FIELD}=?, $trace_table->{qNumUnaligned_FIELD}=?, $trace_table->{sNumUnaligned_FIELD}=?, $trace_table->{qSeqIndIdentical_FIELD}=?,$trace_table->{sSeqIndIdentical_FIELD}=?, $trace_table->{qStrand_FIELD}=?, $trace_table->{sStrand_FIELD}=?, $trace_table->{RANK_FIELD}=? where $trace_table->{TRACE_ID_FIELD}=?";
    my $update_q = $dbh->prepare($sql);

    my $factory = Bio::Tools::Run::StandAloneBlast->new(
                'program' => 'blastn',
                'F'=>'F', #must set filter to be false to avoid missing the longest alignment b/c low complexity filter 
                 );

    $sql = "select distinct a.$trace_table->{TRACE_ID_FIELD}, a.$trace_table->{TRACE_PLATE_FIELD}, a.$trace_table->{TRACE_POS_FIELD}, $trace_table->{TRACE_SEQ_FIELD}, c.$refseq_table->{REFSEQ_ID_FIELD}, c.$refseq_table->{REFSEQ_FIELD} from (select * from $trace_table->{TABLENAME} ";
    $sql .= " where " . $trace_table->{FILTER} if defined $trace_table->{FILTER};
    $sql .= ") a inner join $cherrymap_table->{TABLENAME} b on $tables->{JOIN_CONDITION}{TRACE_CHERRY} inner join $refseq_table->{DBNAME}.$refseq_table->{TABLENAME} c on $tables->{JOIN_CONDITION}{CHERRY_REFSEQ} ";
    $sql .= " where $tables->{FILTER}{BL2SEQ} " if defined($tables->{FILTER}{BL2SEQ});
    $sql .= " order by $tables->{ORDER_BY}{BL2SEQ}" if defined($tables->{ORDER_BY}{BL2SEQ});
    $query = $dbh->prepare($sql);
    $rtn = $query -> execute();
    my $count_hi_score_seq = 0;
    my $count_total_seq = 0;
    while (my ($trace_id, $trace_plate, $trace_pos, $trace_seq, $refseq_id, $refseq) = $query->fetchrow()) {
      $count_total_seq++;
      if ($trace_seq eq "") { # skip empty/unqualified trace
         next;
      }
#print "$trace_id\t$trace_plate\t$trace_pos\n";

      my $traceobj= Bio::Seq->new( -display_id => $trace_id,
                                   -seq => $trace_seq);
      my $cdsobj=Bio::Seq->new( -display_id => $refseq_id,
                                -seq => $refseq);
      $factory->io->_io_cleanup();
      my $report = $factory->bl2seq($traceobj, $cdsobj);
#      my $hi_score = 0;
#      my $insert_start=-1;
#      my $insert_end=-1;
#      my $aln_length=-1;
#      my $percent_identity = -1;
      my %col2insert = ();

      if ($report) {
        while (my $res=$report->next_result){
          my ($previous_qStart, $previous_qEnd) = (0,0);
          while (my $hit= $res->next_hit){
            last if($hit->rank>2 || $hit->score<MIN_BL_SCORE || ($previous_qStart<=$hit->end('query') && $previous_qStart>=$hit->start('query') && $hit->end('query')-$previous_qStart>MIN_HIT_OVERLAP) || ($hit->start('query')>=$previous_qStart && $hit->start('query')<=$previous_qEnd && $previous_qEnd-$hit->start('query')>MIN_HIT_OVERLAP)); # discard rank 3 and above hits, and those score are too low, and hit overlaps with higher rank
            $count_hi_score_seq ++;
            $previous_qStart = $hit->start('query');
            $previous_qEnd = $hit->end('query');      
            $col2insert{score} = $hit->score;
            $col2insert{bits} = $hit->bits;
            $col2insert{numHSP} = $hit->num_hsps;
            $col2insert{qLenAlign} = $hit->length_aln('query');
            $col2insert{sLenAlign} = $hit->length_aln('sbjct');
            $col2insert{qStart} = $hit->start('query');
            $col2insert{qEnd} = $hit->end('query');
            $col2insert{sStart} = $hit->start('sbjct');
            $col2insert{sEnd} = $hit->end('sbjct');
            $col2insert{qNumUnaligned} = $hit->num_unaligned_query;
            $col2insert{sNumUnaligned} = $hit->num_unaligned_hit;
            $col2insert{qSeqIndIdentical} = join('::', $hit->seq_inds('query', 'identical', 1));
            $col2insert{sSeqIndIdentical} = join('::', $hit->seq_inds('sbjct', 'identical', 1));
            $col2insert{qStrand} = join('::', $hit->strand('query'));
            $col2insert{sStrand} = join('::', $hit->strand('sbjct'));
            $col2insert{rank} = join ('::', $hit->rank);

            $rtn = $update_q -> execute($col2insert{score}, 
                                        $col2insert{bits},
                                        $col2insert{numHSP},
					$col2insert{qLenAlign},
					$col2insert{sLenAlign},
                                        $col2insert{qStart},
                                        $col2insert{qEnd},
                                        $col2insert{sStart},
                                        $col2insert{sEnd},
                                        $col2insert{qNumUnaligned},
                                        $col2insert{sNumUnaligned},
                                        $col2insert{qSeqIndIdentical},
                                        $col2insert{sSeqIndIdentical},
                                        $col2insert{qStrand},
                                        $col2insert{sStrand},
                                        $col2insert{rank},
                                        $trace_id);	
          } #end of hit-while
        } #end of result-while
      } # end of if-report
    } # end of while-0

    $query -> finish;
    $update_q -> finish;

    return $count_hi_score_seq;  

} # end of general_bl2seq_brief_blasatResult



##########################################################################
# get bl2 result
##########################################################################
sub analyze_bl2_result {
  my ($pkg, $dbh, $trace_plates, $result_tables, $fields, $avg_qual_cutoff) = @_;
  my $GOOD_RESULT_SCORE = 100;
  my $filter = "TRACE_PLA in ( ". join(',', @{$trace_plates}) . ")";
  my $filter_qual = "TRACE_AVG_QUAL>=$avg_qual_cutoff";
  my $filter_bl2 = "BL2_SCORE>=$GOOD_RESULT_SCORE";

  my %result_list=();
  my ($sql, $query, $rtn);
  my $bl2_result_for_table = "("
   . "select TRACE_PLA, TRACE_POS, BL2_SCORE, TRACE_AVG_QUAL from $result_tables->{TRACE_TABLE} "
   . "where TRACE_DIR='FOR' and $filter"
   . ")";
  my $bl2_result_rev_table = "("
   . "select TRACE_PLA, TRACE_POS, BL2_SCORE, TRACE_AVG_QUAL from $result_tables->{TRACE_TABLE} "
   . "where TRACE_DIR='REV' and $filter"
   . ")";

# first find total seq number of all assoc. plates:
  $sql = "select count(distinct TRACE_PLA, TRACE_POS) from ($bl2_result_for_table union $bl2_result_rev_table) a";
   $query = $dbh->prepare($sql); 
   $query -> execute();
   $result_list{TOTAL_TRIED_SEQ} = $query->fetchrow();
   $query->finish;
   printf "TOTAL TRIED SEQ in either direction passed quality cutoff($avg_qual_cutoff) : $result_list{TOTAL_TRIED_SEQ}\n"; 

# first find total seq number that pass the quality cutoff:
  $sql = "select count(distinct TRACE_PLA, TRACE_POS) from ($bl2_result_for_table union $bl2_result_rev_table) a"
       . " where $filter_qual";
   $query = $dbh->prepare($sql); 
   $query -> execute();
   $result_list{TOTAL_VALID_SEQ_EITHER} = $query->fetchrow();
   $query->finish;
   printf "TOTAL VALID SEQ in either direction passed quality cutoff($avg_qual_cutoff) : $result_list{TOTAL_VALID_SEQ_EITHER}\n"; 

# then find total seq number that pass the quality cutoff for forward:
  $sql = "select count(distinct TRACE_PLA, TRACE_POS) from $bl2_result_for_table a "
       . " where $filter_qual";
   $query = $dbh->prepare($sql); 
   $query -> execute();
   $result_list{TOTAL_VALID_SEQ_FOR} = $query->fetchrow();
   $query->finish;
   printf "TOTAL VALID SEQ in forward direction passed quality cutoff($avg_qual_cutoff) : $result_list{TOTAL_VALID_SEQ_FOR}\n"; 

# then find total seq number that pass the quality cutoff for reverse:
  $sql = "select count(distinct TRACE_PLA, TRACE_POS) from $bl2_result_rev_table a "
       . " where $filter_qual";
   $query = $dbh->prepare($sql); 
   $query -> execute();
   $result_list{TOTAL_VALID_SEQ_REV} = $query->fetchrow();
   $query->finish;
   printf "TOTAL VALID SEQ in reverse direction passed quality cutoff($avg_qual_cutoff) : $result_list{TOTAL_VALID_SEQ_REV}\n"; 

# then find total seq number that pass the quality cutoff for both directions:
  $sql = "select count(distinct a.TRACE_PLA, a.TRACE_POS) from $bl2_result_for_table a inner join $bl2_result_rev_table b using(TRACE_PLA, TRACE_POS)"
       . " where a.$filter_qual and b.$filter_qual";
   $query = $dbh->prepare($sql); 
   $query -> execute();
   $result_list{TOTAL_VALID_SEQ_BOTH}= $query->fetchrow();
   $query->finish;
   printf "TOTAL VALID SEQ in both direction passed quality cutoff($avg_qual_cutoff) : $result_list{TOTAL_VALID_SEQ_BOTH}\n"; 

# find BOTH GOOD pass bl2 check
  $sql = "update $result_tables->{CHERRYMAP_TABLE} a set a.$fields->{RUN_BL2_RESULT}='BOTH_GOOD' where a.$fields->{RUN_TRACE_PLA} in (" .join(',', @{$trace_plates}) . ") and exists (select * from $bl2_result_for_table b where a.$fields->{RUN_TRACE_PLA}=b.TRACE_PLA and a.$fields->{RUN_TRACE_POS}=b.TRACE_POS and b.$filter_bl2 and b.$filter_qual) and exists (select * from $bl2_result_rev_table c where a.$fields->{RUN_TRACE_PLA}=c.TRACE_PLA and a.$fields->{RUN_TRACE_POS}=c.TRACE_POS  and c.$filter_bl2 and c.$filter_qual) "; 
  $rtn = $dbh->do($sql); 
  $result_list{GOOD_SEQ_BOTH} = $rtn;
  printf("GOOD SEQUENCES : %d, PERCENTAGE\%: %.2f\n", $rtn,  100*$rtn/$result_list{TOTAL_VALID_SEQ_EITHER}); 
#  printf("BOTH GOOD SEQUENCES : %d\n", $rtn); 
 
# find FOR GOOD out of those pass qual
  $sql = "update $result_tables->{CHERRYMAP_TABLE} a set a.$fields->{RUN_BL2_RESULT}='FOR_GOOD' where a.$fields->{RUN_TRACE_PLA} in (" .join(',', @{$trace_plates}) . ") and exists (select * from $bl2_result_for_table b where a.$fields->{RUN_TRACE_PLA}=b.TRACE_PLA and a.$fields->{RUN_TRACE_POS}=b.TRACE_POS and b.$filter_bl2 and b.$filter_qual) and not exists (select * from $bl2_result_rev_table c where a.$fields->{RUN_TRACE_PLA}=c.TRACE_PLA and a.$fields->{RUN_TRACE_POS}=c.TRACE_POS  and c.$filter_bl2 and c.$filter_qual)";
  $rtn = $dbh->do($sql); 
  $result_list{GOOD_SEQ_FOR} = $rtn;
  printf("ONLY FORWARD GOOD SEQUENCES : %d, PERCENTAGE\%: %.2f\n", $rtn,  100*$rtn/$result_list{TOTAL_VALID_SEQ_EITHER}); 
#  printf("ONLY FORWARD GOOD SEQUENCES : %d\n", $rtn); 

# find REV GOOD 
  $sql = "update $result_tables->{CHERRYMAP_TABLE} a set a.$fields->{RUN_BL2_RESULT}='REV_GOOD' where a.$fields->{RUN_TRACE_PLA} in (" .join(',', @{$trace_plates}) . ") and not exists (select * from $bl2_result_for_table b where a.$fields->{RUN_TRACE_PLA}=b.TRACE_PLA and a.$fields->{RUN_TRACE_POS}=b.TRACE_POS and b.$filter_bl2 and b.$filter_qual) and exists (select * from $bl2_result_rev_table c where a.$fields->{RUN_TRACE_PLA}=c.TRACE_PLA and a.$fields->{RUN_TRACE_POS}=c.TRACE_POS  and c.$filter_bl2 and c.$filter_qual)";
  $rtn = $dbh->do($sql); 
  $result_list{GOOD_SEQ_REV} = $rtn;
  printf("ONLY REVERSE GOOD SEQUENCES : %d, PERCENTAGE\%: %.2f\n", $rtn,  100*$rtn/$result_list{TOTAL_VALID_SEQ_EITHER}); 
#  printf("ONLY REVERSE GOOD SEQUENCES : %d\n", $rtn); 

# find BOTH BAD 
  $sql = "update $result_tables->{CHERRYMAP_TABLE} a set a.$fields->{RUN_BL2_RESULT}='BOTH_BAD' where a.$fields->{RUN_TRACE_PLA} in (" .join(',', @{$trace_plates}) . ") and not exists (select * from $bl2_result_for_table b where a.$fields->{RUN_TRACE_PLA}=b.TRACE_PLA and a.$fields->{RUN_TRACE_POS}=b.TRACE_POS and b.$filter_bl2 and b.$filter_qual) and not exists (select * from $bl2_result_rev_table c where a.$fields->{RUN_TRACE_PLA}=c.TRACE_PLA and a.$fields->{RUN_TRACE_POS}=c.TRACE_POS  and c.$filter_bl2 and c.$filter_qual)";
  $rtn = $dbh->do($sql); 
  $result_list{BAD_SEQ} = $rtn;
  $result_list{BAD_SEQ_PASS_QUAL} = $result_list{BAD_SEQ}-$result_list{TOTAL_TRIED_SEQ}+$result_list{TOTAL_VALID_SEQ_EITHER};
  printf("BOTH BAD SEQUENCES of all: %d, BOTH BAD SEQ, but at least one end pass quality check: %d, PERCENTAGE\%: %.2f\n", $rtn, $result_list{BAD_SEQ_PASS_QUAL}, 100*$result_list{BAD_SEQ_PASS_QUAL}/$result_list{TOTAL_VALID_SEQ_EITHER}); 
#  printf("BOTH BAD SEQUENCES : %d\n", $rtn); 

  return \%result_list; 
}

##########################################################################
# run seq align using Changyu's seq align utility
##########################################################################
sub run_seq_align {
  my ($pkg, $dbh, $dbh_align, $trace_plates, $result_tables, $fields, $proj_info) = @_;

# clean up previous result
  MY_SeqAlign->reset_result_tables_seq_align($dbh, $dbh_align, $trace_plates, $result_tables, $fields);
 
  MY_SeqAlign->load_refseq_map_seq_align($dbh, $dbh_align, $trace_plates, $result_tables, $fields, $proj_info);

  # run seq_align
  system("~/bin/perl/seq_align.pl -p " . join(",", @{$trace_plates}));

  printf "DONE with running seq_align.\n";

}

################################################################
# SUBROUTINE: general_qualtrim
#
# $tables contains TRACE_SEQ-related table info 
#  
#######################################################

sub general_qualtrim {
  my ($pkg, $tables, $rtinfo) = @_;
  my ($sql, $query, $rtn);

  if(!defined($tables->{TRACE_SEQ}{TRACE_ID_FIELD})) {
      $tables->{TRACE_SEQ}{TRACE_ID_FIELD} = 'TRACE_ID';
  }
  if(!defined($tables->{TRACE_SEQ}{TRACE_QUAL_FIELD})) {
      $tables->{TRACE_SEQ}{TRACE_QUAL_FIELD} = 'TRACE_QUAL';
  }
  if(!defined($tables->{TRACE_SEQ}{TRACE_START_FIELD})) {
      $tables->{TRACE_SEQ}{TRACE_START_FIELD} = 'VTRIM_START';
  }
  if(!defined($tables->{TRACE_SEQ}{HIGH_QUAL_START_FIELD})) {
      $tables->{TRACE_SEQ}{HIGH_QUAL_START_FIELD} = 'HQ_START';
  }
  if(!defined($tables->{TRACE_SEQ}{HIGH_QUAL_END_FIELD})) {
      $tables->{TRACE_SEQ}{HIGH_QUAL_END_FIELD} = 'HQ_END';
  }


  $sql = "select $tables->{TRACE_SEQ}{TRACE_ID_FIELD}, $tables->{TRACE_SEQ}{TRACE_QUAL_FIELD}, $tables->{TRACE_SEQ}{TRACE_START_FIELD} as TRACE_START from $tables->{TRACE_SEQ}{TABLENAME} " . (defined($tables->{FILTER}{QUAL_TRIM})? "where $tables->{FILTER}{QUAL_TRIM} " : "")
#       . " and TRACE_ID=1192988"
        . ";";
  $query = $tables->{TRACE_SEQ}{DBH}->prepare($sql);
  
  $sql = "update $tables->{TRACE_SEQ}{TABLENAME} set $tables->{TRACE_SEQ}{HIGH_QUAL_START_FIELD}=?, $tables->{TRACE_SEQ}{HIGH_QUAL_END_FIELD}=? where $tables->{TRACE_SEQ}{TRACE_ID_FIELD}=?";
  my $query_u = $tables->{TRACE_SEQ}{DBH}->prepare($sql);

  $query->execute();
  while (my ($id, $qual, $vectrim_start) = $query->fetchrow_array()) {
     my ($high_qual_block_start, $high_qual_block_end) = (0,0);
     my ($high_qual_start, $high_qual_end) = (-1,-1); # set flag
#     my ($TRIM_WINSIZE, $AVG_QUAL, $IGNORE_LEN) = (20, 30, 50);
     my $TRIM_WINSIZE = defined($rtinfo->{TRIM_WINSIZE})?$rtinfo->{TRIM_WINSIZE}:20;
     my $AVG_QUAL = defined($rtinfo->{AVG_QUAL})?$rtinfo->{AVG_QUAL}:30; 
     my $IGNORE_LEN = defined($rtinfo->{IGNORE_LEN})?$rtinfo->{IGNORE_LEN}:50;

     my @qual_array = split /\s/, $qual;

#Search for a block starting from insert_start, that is 20 nt long with at least avg score of 30
     my $search_start = 0;
     my $low_qual_nt = 0;
     my $MAX_LOW_QUAL_GAP =  defined($rtinfo->{MAX_LOW_QUAL_GAP})?$rtinfo->{MAX_LOW_QUAL_GAP}:10; # maxium seq gap of low quality allowed, the bigger, the looser, default=5
     my $AVG_QUAL_LOW = $AVG_QUAL-10; # qual threshold set for low qual seq gap

     (defined $vectrim_start) and ($search_start=$vectrim_start-1);
     for (my $i=$search_start;$i<scalar(@qual_array)-$TRIM_WINSIZE+1;$i++) {
        my $sub_array_qual=0;
        for (my $j=$i;$j<=$i+$TRIM_WINSIZE-1;$j++) {
           $sub_array_qual += $qual_array[$j];
        }
        my $avg_sub_array_qual = $sub_array_qual / $TRIM_WINSIZE;
        if ($avg_sub_array_qual >= $AVG_QUAL) {
          $low_qual_nt=0; # reset low qual seq counter
          if ($high_qual_block_start>0) {
	     next;
          }
          else {
	     $high_qual_block_start = $i;
          }
        }
        elsif ($i-$search_start<=$IGNORE_LEN)  { # add robust by ignoring the begining of the sequence
          $high_qual_block_start=0; # discard the block, restart search
          $low_qual_nt=0; # reset low qual seq counter
          next;
        }
#=pod
        elsif ($low_qual_nt<$MAX_LOW_QUAL_GAP && $avg_sub_array_qual>=$AVG_QUAL_LOW)  { # add more robust by allowing certain gap of low quality data
          $low_qual_nt++;
          next;
        }
#=cut
        else  { 
          if ($high_qual_block_start>0) {
	     $high_qual_block_end = $i+$TRIM_WINSIZE-1;
	     last;
          }
          else {
	     next;
          }
        }
    } # end of outside for-loop

    if (($high_qual_block_start>0) && ($high_qual_block_end==0)) {
        $high_qual_block_end = @qual_array - 1;
    }

  #$high_qual_start = $high_qual_block_start;
  #$high_qual_end = $high_qual_block_end;
  #($high_qual_start==0) and print "$id\n";

#chop off from both end, individual nts with scores below 10
    my $LEAST_QUAL = 10;
    if ($high_qual_block_start>0) {
      for (my $k=$high_qual_block_start; $k<$high_qual_block_end+1; $k++) {
         $high_qual_start = $k+1; #1-based
         ($qual_array[$k]>= $LEAST_QUAL) and last;
      }
    
      for (my $l=$high_qual_block_end;$l>=$high_qual_block_start;$l--){
         $high_qual_end = $l+1; #1-based
         ($qual_array[$l] >= $LEAST_QUAL) and last;
      }
    }

  #print "$id\t$high_qual_block_start\t$high_qual_block_end\t$high_qual_start\t$high_qual_end\n";

    $query_u->execute($high_qual_start, $high_qual_end, $id);

  } # end of while
  $query->finish;
  $query_u->finish;
} # general_qualtrim

################################################################
# SUBROUTINE: qualtrim algorithm 2: using conventional algorithm
# INPUT PARAM:
#  my %result_tables = (
#     PLATE_MAP_TABLE => 'SACOME_TRACE_PLATEMAP',
#     CHERRYMAP_TABLE => 'SACOME_TRACE_CHERRYMAP', # to distinguish with the original existing SACOME_CHERRYMAP which records all the plate/pos mapping info for each step
#     PLATE_DETAIL_TABLE => 'SACOME_PLATE_DETAIL',
#     TRACE_TABLE => 'SACOME_TRACE_SEQ',
#     REFSEQ_TABLE => 'horfeome7.cDNA',
#     CONTIG_TABLE=>'SACOME_TRACE_CONTIG',
#     SINGLET_TABLE=>'SACOME_TRACE_SINGLET',
#  );
#
#  my %fields = (
#     RUN_TRACE_PLA => 'TRACE_PLA',
#     RUN_TRACE_POS => 'TRACE_POS',
#     RUN_TRACE_POSID => 'TRACE_POS_ID',
#     RUN_BL2_RESULT => 'BL2_RESULT',
#     DstPlate => 'DstPlate',
#     DstPosID => 'DstPosID',
#  );
#
#  my %trace_plates = (
#     10818=>'SACHs_ORF000001', #mapping the trace id and its templ name. 
#     10819=>'SACHs_ORF000001', 
#     10820=>'SACHs_ORF000001', 
#     10821=>'SACHs_ORF000001',
#  );
#
#  my %proj_info = (
#     PROJ_ABBR => 'SAC',
#     SEQ_TYPE => 'RTPOOL',
#     ORGANISM => 'Human',
#     NOTE => 'WildType',
#  );
#  
#######################################################

sub qualtrim_trace_seq2 {
  my ($pkg, $dbh, $trace_plates, $result_tables, $fields) = @_;
  my ($sql, $query, $rtn);

  $sql = "select TRACE_ID, TRACE_QUAL, VECTRIM_START from $result_tables->{TRACE_TABLE}{TABLENAME} where TRACE_PLA in (" . join (',', @{$trace_plates}). ")"
#       . " and TRACE_ID=1072577"
        . ";";
  $query = $dbh->prepare($sql);
  
  $sql = "update $result_tables->{TRACE_TABLE}{TABLENAME} set HIGH_QUAL_START=?, HIGH_QUAL_END=? where TRACE_ID=?";
  my $query_u = $dbh->prepare($sql);

  $query->execute();
  while (my ($id, $qual, $vectrim_start) = $query->fetchrow_array()) {
     my ($high_qual_block_start, $high_qual_block_end) = (0,0);
     my ($high_qual_start, $high_qual_end) = (-1,-1); # set flag
     my ($TRIM_WINSIZE, $AVG_QUAL) = (20, 30);
     my @qual_array = split /\s/, $qual;

#Search for a block starting from insert_start, that is 20 nt long with at least avg score of 30
     my $search_start = 0;
     (defined $vectrim_start) and ($search_start=$vectrim_start-1);
     for (my $i=$search_start;$i<scalar(@qual_array)-$TRIM_WINSIZE+1;$i++) {
        my $sub_array_qual=0;
        for (my $j=$i;$j<=$i+$TRIM_WINSIZE-1;$j++) {
           $sub_array_qual += $qual_array[$j];
        }
        my $avg_sub_array_qual = $sub_array_qual / $TRIM_WINSIZE;
        if ($avg_sub_array_qual >= $AVG_QUAL) {
          if ($high_qual_block_start>0) {
	     next;
          }
          else {
	     $high_qual_block_start = $i;
          }
        }
        else {
          if ($high_qual_block_start>0) {
	     $high_qual_block_end = $i+$TRIM_WINSIZE-1;
	     last;
          }
          else {
	     next;
          }
        }
    } # end of outside for-loop

    if (($high_qual_block_start>0) && ($high_qual_block_end==0)) {
        $high_qual_block_end = @qual_array - 1;
    }

  #$high_qual_start = $high_qual_block_start;
  #$high_qual_end = $high_qual_block_end;
  #($high_qual_start==0) and print "$id\n";

#chop off from both end, individual nts with scores below 10
    my $LEAST_QUAL = 10;
    if ($high_qual_block_start>0) {
      for (my $k=$high_qual_block_start; $k<$high_qual_block_end+1; $k++) {
         $high_qual_start = $k+1;
         ($qual_array[$k]>= $LEAST_QUAL) and last;
      }
    
      for (my $l=$high_qual_block_end;$l>=$high_qual_block_start;$l--){
         $high_qual_end = $l+1;
         ($qual_array[$l] >= $LEAST_QUAL) and last;
      }
    }

  #print "$id\t$high_qual_block_start\t$high_qual_block_end\t$high_qual_start\t$high_qual_end\n";

    $query_u->execute($high_qual_start, $high_qual_end, $id);

  } # end of while
  $query->finish;
  $query_u->finish;
} # qualtrim_trace_seq2


###############################################################
# sub: run_phrap
# desc: run phrap for single seq. after vector and quality trim 
#       has been done
###############################################################
sub run_phrap {
  my ($pkg, $dbh, $dbname, $trace_plates, $tables, $fields, $phrap_dir) = @_;
  my ($sql, $query, $rtn);
  my $seq_fileprefix;
  my $seq_fasta;
  my $seq_qual; 
  my $phrap_out;
  my $phrap_err;

  $sql = "select distinct a.ORF_ID, b.TRACE_PLATE, b.TRACE_POS from $tables->{CHERRYMAP}{TABLENAME} a inner join $tables->{TRACE_SEQ}{TABLENAME} b on a.dstPlate=b.TRACE_PLATE and a.dstPos=b.TRACE_POS where b.TRACE_PLATE in (" . join(',', @{$trace_plates}) . ") ";
#  $sql .= " and b.TRACE_POS='F02' "; # debug
  $sql .= " order by a.ORF_ID, b.TRACE_PLATE, b.TRACE_POS";
  $query = $dbh->prepare($sql);
  $query -> execute();

# use vector trim instead of quality trim, and only consider those succeeded in vector trim:
#  $sql = "select b.TRACE_ID, b.TRACE_DIR, b.HIGH_QUAL_START, b.HIGH_QUAL_END, b.TRACE_SEQ, b.TRACE_QUAL from $tables->{TRACE_SEQ}{TABLENAME} b where b.TRACE_PLATE=? and b.TRACE_POS=?  and b.HIGH_QUAL_START>-1 order by b.TRACE_DIR";

#  $sql = "select b.TRACE_ID, b.TRACE_DIR, substr(b.TRACE_SEQ, b.HIGH_QUAL_START, b.HIGH_QUAL_END-b.HIGH_QUAL_START+1) as TRIM_TRACE_SEQ, substring_index(substring_index(b.TRACE_QUAL,' ', b.HIGH_QUAL_END+1), ' ',-b.HIGH_QUAL_END+b.HIGH_QUAL_START+1) as TRIM_TRACE_QUAL from $tables->{TRACE_SEQ}{TABLENAME} b where b.TRACE_PLATE=? and b.TRACE_POS=?  and b.HIGH_QUAL_START>-1 order by b.TRACE_DIR";

# for trim start/end is null, set the start from very beginning, and end to be the end of the trace:
#  $sql = "select b.TRACE_ID, b.TRACE_DIR, if(b.VTRIM_START is null, 1, b.VTRIM_START), if(b.VTRIM_END is null, length(b.TRACE_SEQ), b.VTRIM_END), b.TRACE_SEQ, b.TRACE_QUAL from $tables->{TRACE_SEQ}{TABLENAME} b where b.TRACE_PLATE=? and b.TRACE_POS=? and b.VTRIM_START is not null order by b.TRACE_DIR";

  $sql = "select b.TRACE_ID, b.TRACE_DIR, if(b.VTRIM_START is null, 1, b.VTRIM_START), if(b.VTRIM_END is null, length(b.TRACE_SEQ), b.VTRIM_END), b.TRACE_SEQ, b.TRACE_QUAL from $tables->{TRACE_SEQ}{TABLENAME} b where b.TRACE_PLATE=? and b.TRACE_POS=? order by b.TRACE_DIR"; # remove trim start restriction

#  $sql = "select b.TRACE_ID, b.TRACE_DIR, 1 as VTRIM_START, length(b.TRACE_SEQ) as VTRIM_END, b.TRACE_SEQ, b.TRACE_QUAL from $tables->{TRACE_SEQ}{TABLENAME} b where b.TRACE_PLATE=? and b.TRACE_POS=? order by b.TRACE_DIR";

  my $query_trace = $dbh->prepare($sql);
  
  if (!-e $phrap_dir) {
    system("mkdir -p $phrap_dir");
  }

# change to the result directory:
  my $script_dir = Cwd::getcwd();
  chdir($phrap_dir);
  
  my $seq_name;
  my $trim_seq;
  my @trim_qual;
  while(my ($ORF_ID, $TRACE_PLA, $TRACE_POS) = $query->fetchrow()) {
    $seq_fileprefix = 'O' . $ORF_ID . "." . $TRACE_PLA . $TRACE_POS;
    $seq_fasta = $seq_fileprefix . '.fa';
    $seq_qual =  $seq_fileprefix . '.fa.qual';
    $phrap_out =  $seq_fileprefix . '.phrap.out';
    $phrap_err =  $seq_fileprefix . '.phrap.err';

    open OUT, ">$seq_fasta";
    open OUT_Q, ">$seq_qual";
    $query_trace->execute($TRACE_PLA, $TRACE_POS);
    while(my ($TRACE_ID, $TRACE_DIR, $TRIM_START, $TRIM_END, $TRACE_SEQ, $TRACE_QUAL) = $query_trace->fetchrow()) {
      $seq_name = 'O' . $ORF_ID . "." . 'T' . $TRACE_ID . "." . $TRACE_PLA . $TRACE_POS . $TRACE_DIR;
      print OUT ">" . $seq_name . " DIRECTION: $TRACE_DIR CHEM: unknown  DYE: unknown TEMPLATE: HORF$ORF_ID" . "_" . $TRACE_PLA . $TRACE_POS ."\n";
      $trim_seq = substr($TRACE_SEQ, $TRIM_START-1, $TRIM_END-$TRIM_START+1); # notice all position are 1-based in DB
      print "seq len : " . length($trim_seq) . "\n";
      print OUT "$trim_seq\n";
      print OUT_Q ">" . $seq_name . " DIRECTION: $TRACE_DIR CHEM: unknown  DYE: unknown TEMPLATE: HORF$ORF_ID" . "_" . $TRACE_PLA . $TRACE_POS . "\n";
      my @qual_array = split /\s/, $TRACE_QUAL;
      @trim_qual = splice(@qual_array, $TRIM_START-1, $TRIM_END-$TRIM_START+1);
      print "qual len : " . scalar(@trim_qual) . "\n";
      print OUT_Q join(" ", @trim_qual);
      print OUT_Q "\n";
    }
    close OUT;
    close OUT_Q;

# run phrap:
    system("phrap -new_ace $seq_fasta" . " >$phrap_out >& $phrap_err");  

# load phrap result to db:
    MY_SeqAlign->loadPhrapResult($dbh, $dbname, $seq_fasta, $tables, $phrap_dir);
  }
  $query->finish;
  $query_trace->finish;

# restore working directory:
  chdir($script_dir);
    
} # run_phrap

###############################################################
# sub: general_tcoffee
# desc: prepare fasta, and then run t-coffee
#       
###############################################################
sub general_tcoffee {
  my ($pkg, $tables, $rtinfo)=@_;
  my ($sql, $query, $rtn);

# change to the result directory:
  my $script_dir = Cwd::getcwd();
  my $tcoffee_dir = $rtinfo->{RESULT_DIR} . "tcoffee/";
  if(!-e $tcoffee_dir) {
    system("mkdir -p $tcoffee_dir");
  }
  chdir($tcoffee_dir);
  
# prepare fasta file for t-coffee:
  
  $sql = "select distinct GROUP_ID from $tables->{SEQ}{TABLENAME} a order by GROUP_ID";
  $query = $tables->{SEQ}{DBH} -> prepare($sql);
  $sql = "select SEQ_NAME, SEQ from $tables->{SEQ}{TABLENAME} a where GROUP_ID=? order by SEQ_NAME";  
  my $query_seq = $tables->{SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  while (my ($GROUP_ID) = $query->fetchrow()) {
    open OUT, ">$GROUP_ID" . ".fa";
    $rtn = $query_seq -> execute($GROUP_ID);
    while (my($SEQ_NAME, $SEQ) = $query_seq->fetchrow()) {
      print OUT ">$SEQ_NAME\n";
      print OUT "$SEQ\n";
    }
    close OUT;
  }
  $query_seq -> finish;
  $query -> finish;

# run t-coffee:
  
  my $filelist = "temp_tcoffee_input.list";
  system("find -name '*.fa' > $filelist");
  open LIST, $filelist;
  while (<LIST>) {
    chomp;
    system("t_coffee $_>t_coffee.out >& t_coffee.err");
  }
  close LIST;

# remove temporary file:
  system("rm temp_tcoffee_input.list");

# chg back to the script dir:
  chdir($script_dir);

} # end of tcoffee_general


###############################################################
# sub: run_tcoffee
# desc: prepare fasta, and then run t-coffee
#       
###############################################################
sub run_tcoffee {
  my ($pkg, $dbh, $trace_plates, $result_tables, $fields, $tcoffee_dir, $inc_trace)=@_;
  my ($sql, $query, $rtn);

  if (!-e $tcoffee_dir) {
    system("mkdir -p $tcoffee_dir");
  }

# change to the result directory:
  my $script_dir = Cwd::getcwd();
  chdir($tcoffee_dir);
  
  my $seq_name;

# prepare fasta file for t-coffee:

# get all reference seq:
  $sql = "select distinct a.$result_tables->{JOIN}{CHERRY_REFSEQ}, b.$result_tables->{REFSEQ}{FIELD_LIST}{CDS_SEQ} from $result_tables->{CHERRYMAP}{TABLENAME} a inner join $result_tables->{REFSEQ}{TABLENAME} b using($result_tables->{JOIN}{CHERRY_REFSEQ}) where a.dstPlate in (" . join(',', @{$trace_plates}) . ")";
  $query = $dbh->prepare($sql);
  
# get all contig:
  $sql = "select a.CONTIG_NAME, a.CONTIG_SEQ from $result_tables->{CONTIG_TABLE}{TABLENAME} a where substring(substring_index(a.CONTIG_NAME, '.', 1), 2)=? ORDER BY a.CONTIG_NAME";
  my $query_contig = $dbh->prepare($sql);

# get all singlets:
  $sql = "select substring_index(a.SINGLET_NAME, ' ', 1) as SINGLET_NAME , a.SINGLET_SEQ from $result_tables->{SINGLET_TABLE}{TABLENAME} a where substring(substring_index(a.SINGLET_NAME, '.', 1), 2)=? ORDER BY a.SINGLET_NAME";
  my $query_singlet = $dbh->prepare($sql);

# get all traces:
  $sql = "select a.TRACE_ID, a.TRACE_PLA, a.TRACE_POS, a.TRACE_DIR, substr(a.TRACE_SEQ, HIGH_QUAL_START, HIGH_QUAL_END-HIGH_QUAL_START+1) as TRIMMED_TRACE_SEQ from $result_tables->{TRACE_TABLE}{TABLENAME} a where exists (select * from $result_tables->{CHERRYMAP}{TABLENAME} b where b.$result_tables->{JOIN}{CHERRY_REFSEQ}=? and b.dstPlate=a.TRACE_PLA and b.dstPos=a.TRACE_POS) ORDER BY a.TRACE_DIR";
  my $query_trace = $dbh->prepare($sql);

  $query -> execute();
  while (my ($orf_id, $orf_seq) = $query->fetchrow) {
     my $seq_fasta = 'ORF' . $orf_id . '.fa';
     open OUT, ">" . $seq_fasta;
     $seq_name = "ORF" . $orf_id;
     print OUT ">" . $seq_name . "\n";
     print OUT "$orf_seq\n";

# write traces:
     if($inc_trace==1) { 
        $query_trace->execute($orf_id);
        while(my ($TRACE_ID, $TRACE_PLA, $TRACE_POS, $TRACE_DIR, $TRIMMED_TRACE_SEQ) = $query_trace->fetchrow()) {
          $seq_name = 'O' . $orf_id . "." . $TRACE_PLA . $TRACE_POS . $TRACE_DIR . "T" . $TRACE_ID . "." . ".trace";
          print OUT ">" . $seq_name . "\n";
          print OUT "$TRIMMED_TRACE_SEQ\n";
        }
     }

# write contig if available:
     $query_contig->execute($orf_id);
     while(my ($CONTIG_NAME, $CONTIG_SEQ) = $query_contig->fetchrow()) {
      $seq_name = $CONTIG_NAME;
      print OUT ">" . $seq_name . "\n";
      print OUT "$CONTIG_SEQ\n";
     }

# write singlet if available:
     $query_singlet->execute($orf_id);
     while(my ($SINGLET_NAME, $SINGLET_SEQ) = $query_singlet->fetchrow()) {
      $seq_name = $SINGLET_NAME . ".singlet";
      print OUT ">" . $seq_name . "\n";
      print OUT "$SINGLET_SEQ\n";
     }

     close OUT;

  }

  $query->finish;
  $query_trace->finish;
  $query_contig->finish;
  $query_singlet->finish;

# run t-coffee:
  
  my $filelist = "temp_tcoffee_input.list";
  system("ls *.fa > $filelist");
  open LIST, $filelist;
  while (<LIST>) {
    chomp;
    system("t_coffee $_>t_coffee.out >& t_coffee.err");
  }
  close LIST;

# remove temporary file:
  system("rm temp_tcoffee_input.list");

# chg back to the script dir:
  chdir($script_dir);

} # run_tcoffee


###############################################################
# sub: run_tcoffee_batch
# desc: run tcoffee against the batch list of files, assuming the 
#       input fasta file already prepared
#       
###############################################################
sub run_tcoffee_batch {
  my ($pkg, $tcoffee_dir)=@_;
  my ($sql, $query, $rtn);

  if (!-e $tcoffee_dir) {
    print "ERROR: $tcoffee_dir does not exist. \n";
    return;
  }

# change to the result directory:
  my $script_dir = Cwd::getcwd();
  chdir($tcoffee_dir);

# run tcoffee agains each fasta:
  my @fasta_list = <*.fa>;
  foreach my $fasta (@fasta_list) {
    system("t_coffee $fasta >$fasta" . ".out >& $fasta" . ".err");
  }

# chg back to the script dir:
  chdir($script_dir);

} # run_tcoffee


###############################################################
# INPUT: 
#	load contig to database after running phrap
#          SACTEPM_PL2,10818,SACHs_ORF000001, count=67
#          SACTEPM_PL4,10819,SACHs_ORF000001
#          SACTEPM_PL5,10820,SACHs_ORF000001
#          SACTEPM_PL8,10821,SACHs_ORF000001
# OUTPUT:  
#       MODIFY xxx_TRACE_CHERRYMAP, add ASM_TYPE enum('SINGLET','CONTIG');
#       ADD xxx_TRACE_CONTIG, to store contigs for assembled seq.
#       ADD xxx_TRACE_SINGLET, to store contigs for assembled seq.
###############################################################
sub loadPhrapResult {
  my ($pkg, $dbh, $dbname, $baseName, $result_tables, $phrap_dir) = @_;
  my ($sql, $query, $rtn);

# create table if not exists
  my $exist = MY_DButility->check_table_existence($dbh, $dbname, $result_tables->{CONTIG_TABLE}{TABLENAME});
  if($exist == 0) { # not exists
    $sql = "create table $result_tables->{CONTIG_TABLE}{TABLENAME} ( "
      . "CONTIG_NAME varchar(100), " # store identifier name for later mapping
      . "CONTIG_SEQ mediumtext, "
      . "PRIMARY KEY(CONTIG_NAME) "
      . ");";
    $rtn = $dbh->do($sql); 
  }

  $exist = MY_DButility->check_table_existence($dbh, $dbname, $result_tables->{SINGLET_TABLE}{TABLENAME});
  if($exist == 0) { # not exists
    $sql = "create table $result_tables->{SINGLET_TABLE}{TABLENAME} ( "
      . "SINGLET_NAME varchar(100), " # store identifier name for later mapping
      . "SINGLET_SEQ mediumtext, "
      . "PRIMARY KEY(SINGLET_NAME) "
      . ");";
    $rtn = $dbh->do($sql); 
  }

  my %tableInfo = (
     TABLE_NAME => $result_tables->{CONTIG_TABLE}{TABLENAME},
#     CREATE_NEW => 1,
     COL_LIST => [
	{COL_NO => 0,
         COL_NAME => 'CONTIG_NAME',
         COL_TYPE => 'varchar(100)',
         COL_MODIFIER => 'PRIMARY KEY',
        },
	{COL_NO => 1,
         COL_NAME => 'CONTIG_SEQ',
         COL_TYPE => 'mediumtext',
        }
     ],
  );

#  my $contig_fasta = $phrap_dir . $baseName . ".contigs";
  my $contig_fasta = $baseName . ".contigs"; # already in phrap dir
  if (-e $contig_fasta) {
    MY_BioUtility->load_table_from_fasta($dbh, $dbname, '\n',\%tableInfo, $contig_fasta);
  }
   %tableInfo = (
     TABLE_NAME => $result_tables->{SINGLET_TABLE}{TABLENAME},
#     CREATE_NEW => 1,
     COL_LIST => [
	{COL_NO => 0,
         COL_NAME => 'SINGLET_NAME',
         COL_TYPE => 'varchar(100)',
         COL_MODIFIER => 'PRIMARY KEY',
        },
	{COL_NO => 1,
         COL_NAME => 'SINGLET_SEQ',
         COL_TYPE => 'mediumtext',
        }
     ],
  );

#  my $singlet_fasta = $phrap_dir . $baseName . ".singlets";
  my $singlet_fasta = $baseName . ".singlets";
  if (-e $singlet_fasta) {
    MY_BioUtility->load_table_from_fasta($dbh, $dbname, '\n',\%tableInfo, $singlet_fasta);
  }

} # end of loadPhrapResult


#############################################
# psl sample:
#
#psLayout version 3
#
#match	mismatch 	repmatch 	NCount	Qgap	Qgap	Tgap	Tgap	strand	Qcount        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
#    	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#1509	0	0	0	0	0	5	29109	+	527|51013@B08|GI:51271	1509	0	1509	chr9	140273252	34210912	34241530	6	34,125,924,183,102,141,	0,34,159,1083,1266,1368,	34210912,34224213,34231182,34239776,34240655,34241389,
#						
#
#############################################

sub load_blat_result_old {
  my ($sql, $query, $rtn);
  my ($pkg, $psl, $dbh, $dbname, $tablename_prefix) = @_;
  my %tables = (
    BLAT => {
      DBH => $dbh,
      DBNAME => $dbname, 
      TABLENAME => $tablename_prefix . "_BLAT_BEST",
      COL_LIST => ordered_hash_ref(
        LINE_ID => 'int(8)',
        tName => 'varchar(50)',
        qName => 'varchar(50)',
        strand => "enum('+','-')",
        match_nts => 'int(8)',
        mismatch_nts => 'int(8)',
        tSize => 'int(12)',
        tStart => 'int(12)',
        tEnd => 'int(12)',
        qSize => 'int(12)',
        qStart => 'int(12)',
        qEnd => 'int(12)', 
        blockCount => 'int(3)',
        qNumInsert => 'int(3)',
        qBaseInsert => 'int(3)',
        repMatches => 'int(3)',
        nCount => 'int(3)',
      ),
      PRIMARY_KEY => 'LINE_ID',
    },

    EXON => {
      DBH => $dbh,
      DBNAME => $dbname, 
      TABLENAME => $tablename_prefix . "_BLAT_EXON",
      COL_LIST => ordered_hash_ref(
        LINE_ID => 'int(8)',
        EXON_NR => 'tinyint(2)',
        qStart => 'int(12)',
        qEnd => 'int(12)',
        tStart => 'int(12)',
        tEnd => 'int(12)'
      ),
      PRIMARY_KEY => 'LINE_ID, EXON_NR',
    },
  );

  my %PSL_COL_LIST = (
    matches => { #1509
      COL_NO => 0,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Number of bases that match that aren't repeats",
    },
    misMatches => { #0
      COL_NO => 1, 
      COL_TYPE => 'int unsigned',
      COL_DESC => "Number of bases that don't match",
    },
    repMatches => { #0
      COL_NO => 2,
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of bases that match but are part of repeats",
    },
    nCount => { #0
      COL_NO => 3, 
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of 'N' bases",
    },
    qNumInsert => { #0
      COL_NO => 4, 
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of inserts in query",
    },
    qBaseInsert => { #0
      COL_NO => 5,
      COL_TYPE => 'int unsigned' ,
      COL_TYPE => "Number of bases inserted in query",
    },
    tNumInsert => { #5
      COL_NO => 6,
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of inserts in target",
    },
    tBaseInsert => { #29109
      COL_NO => 7, 
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of bases inserted in target",
    },
    strand => { # +
      COL_NO => 8,
      COL_TYPE => 'char(2)' ,
      COL_DESC => "+ or - for query strand, optionally followed by + or - for target strand",
    },
    qName => { #527|51013@B08|GI:51271
      COL_NO => 9,
      COL_TYEP => 'varchar(255)',
      COL_DESC => "Query sequence name",
    },
    qSize => { #1509
      COL_NO => 10,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Query sequence size",
    },
    qStart => {#0
      COL_NO => 11,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Alignment start position in query",
    },
    qEnd => { #1509
      COL_NO => 12,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Alignment end position in query",
    },
    tName => { #chr9
      COL_NO => 13, 
      COL_TYPE => 'varchar(255)' ,
      COL_DESC => "Target sequence name",
    },
    tSize => {#	140273252
      COL_NO => 14,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Target sequence size",
    },
    tStart => {# 34210912
      COL_NO => 15,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Alignment start position in target",
    },
    tEnd => { #34241530
      COL_NO => 16,
      COL_TYPE => 'int unsigned',
      COL_DESC => "alignment end position in target",
    },
    blockCount => { # 6
      COL_NO => 17, 
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of blocks in alignment",
    },
    blockSizes => { #34,125,924,183,102,141,
       COL_NO => 18,
       COL_TYPE => 'longblob' ,
       COL_DESC => "Size of each block in a comma separated list",
    },
    qStarts => { #0,34,159,1083,1266,1368,
       COL_NO => 19,
       COL_TYPE => 'longblob' ,
       COL_DESC => "Start of each block in query in a comma separated list",
    },
    tStarts => { #34210912,34224213,34231182,34239776,34240655,34241389,
      COL_NO => 20, 
      COL_TYPE => 'longblob' ,
      COL_DESC => "Start of each block in target in a comma separated list",
    },

  );

  MY_DButility -> prepare_tables(\%tables);
  
  $sql = "insert ignore into $tables{BLAT}{TABLENAME} (LINE_ID, tName, qName, strand, match_nts, mismatch_nts, tSize, tStart, tEnd, qSize, qStart, qEnd, blockCount, qNumInsert, qBaseInsert, repMatches, nCount) values (" . "?, "x16 . "?)";
  my $qi_psl = $tables{BLAT}{DBH} -> prepare($sql);
  $sql = "insert ignore into $tables{EXON}{TABLENAME} (LINE_ID, EXON_NR, qStart, qEnd, tStart, tEnd) values (" . "?, "x5 . "?)";
  my $qi_exon = $tables{EXON}{DBH} -> prepare($sql);
  
  open IN, "<$psl";

  my @current_items = ();
  my $current_line_id = 0;
  my @items = ();
  my $score = 0;
#  my $current_line_id;
  my $line;
  my ($qName, $current_qName) = ('', '');

# skip the first 5 line
  for (1..5) { 
    $line = <IN>;
  }

# parse psl lines:

  while ($line = <IN>) {
    next if $line=~/chr.+_.+/; # skip those chromosome like 'chr9_random';
    chomp $line;
    @items = split /\t/, $line;

    ($qName) = $items[$PSL_COL_LIST{qName}{COL_NO}];
    if ($qName ne $current_qName) {#read a new id
      if ($current_qName ne '') {#if current values are not empty
        ($current_qName) = $current_items[$PSL_COL_LIST{qName}{COL_NO}];
        $qi_psl->execute(
                          $current_line_id,
                          $current_items[$PSL_COL_LIST{tName}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qName}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{strand}{COL_NO}],
                          $current_items[$PSL_COL_LIST{matches}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{misMatches}{COL_NO}],
                          $current_items[$PSL_COL_LIST{tSize}{COL_NO}],
                          $current_items[$PSL_COL_LIST{tStart}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{tEnd}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qSize}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qStart}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qEnd}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{blockCount}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qNumInsert}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qBaseInsert}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{repMatches}{COL_NO}],
                          $current_items[$PSL_COL_LIST{nCount}{COL_NO}]
                         );
      
        my @sizes = split /,/, $current_items[$PSL_COL_LIST{blockSizes}{COL_NO}];
        my @query_starts = split /,/, $current_items[$PSL_COL_LIST{qStarts}{COL_NO}];
        my @hit_starts = split /,/, $current_items[$PSL_COL_LIST{tStarts}{COL_NO}];

 #######IMPORTANT!!!!!!!!!!!!!!####################
 #Fix - strand hits so that they can be processed in the same way as the + Stand hits
 
       if ($current_items[$PSL_COL_LIST{strand}{COL_NO}] eq '-') {
	  for (my $i=0; $i<scalar(@sizes);$i++) {
	    $query_starts[$i] = $current_items[$PSL_COL_LIST{qSize}{COL_NO}] - $query_starts[$i] - $sizes[$i]; #The query starts in - sequences are off 1 bit
	  }
	
	  @query_starts = reverse @query_starts;
	  @hit_starts = reverse @hit_starts;
          @sizes = reverse @sizes;
        }

        for (my $i=0;$i<@sizes;$i++) {
	  my $exon_nr = $i+1;
	  my $q_start = $query_starts[$i];
	  my $h_start = $hit_starts[$i];
	  my $q_end = $q_start + $sizes[$i] -1;
	  my $h_end = $h_start + $sizes[$i] - 1;

	  if ($current_items[$PSL_COL_LIST{strand}{COL_NO}] eq '-'){ #for - strand, swap start and end in hit
	    my $tmp = $h_start;
	    $h_start = $h_end;
	    $h_end = $tmp;
	  }

	  $qi_exon->execute($current_line_id, $exon_nr, $q_start+1, $q_end+1, $h_start+1, $h_end+1); # set coordinates to be 1-based
	  #print "$exon_nr\t $q_start\t $q_end\t $hit_starts[$i]\t $h_end\t $current_qName\n";
        }

        $current_qName = $qName;
        @current_items = @items;
        $score = $items[$PSL_COL_LIST{matches}{COL_NO}]-$items[$PSL_COL_LIST{misMatches}{COL_NO}];
        $current_line_id = $.;
      }
      else { # first record
        @current_items = @items;
        $current_qName = $qName;
        $current_line_id = $.;
      }
    } # end of check new orf alignment
    elsif ($items[$PSL_COL_LIST{matches}{COL_NO}]-$items[$PSL_COL_LIST{misMatches}{COL_NO}] > $score) {#replace current items when the score is higher, so that the current values alway represent the highest score hit of the query
      $score = $items[$PSL_COL_LIST{matches}{COL_NO}]-$items[$PSL_COL_LIST{misMatches}{COL_NO}];
      @current_items = @items;
      $current_line_id = $.;
    }
  } # end of parse psl file

  close IN;

#Now load the last hit

  ($current_qName) = $current_items[$PSL_COL_LIST{qName}{COL_NO}];
   $qi_psl->execute(
                    $current_line_id,
                    $current_items[$PSL_COL_LIST{tName}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qName}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{strand}{COL_NO}],
                    $current_items[$PSL_COL_LIST{matches}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{misMatches}{COL_NO}],
                    $current_items[$PSL_COL_LIST{tSize}{COL_NO}],
                    $current_items[$PSL_COL_LIST{tStart}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{tEnd}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qSize}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qStart}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qEnd}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{blockCount}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qNumInsert}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qBaseInsert}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{repMatches}{COL_NO}],
                    $current_items[$PSL_COL_LIST{nCount}{COL_NO}]
                   );

  my @sizes = split /,/, $current_items[$PSL_COL_LIST{blockSizes}{COL_NO}];
  my @query_starts = split /,/, $current_items[$PSL_COL_LIST{qStarts}{COL_NO}];
  my @hit_starts = split /,/, $current_items[$PSL_COL_LIST{tStarts}{COL_NO}];

#######IMPORTANT!!!!!!!!!!!!!!####################
#Fix - strand hits so that they can be processed in the same way as the + Stand hits
  if ($current_items[$PSL_COL_LIST{strand}{COL_NO}] eq '-') {
    for (my $i=0; $i<scalar(@sizes);$i++) {
      $query_starts[$i] = $current_items[$PSL_COL_LIST{qSize}{COL_NO}] - $query_starts[$i] - $sizes[$i]; #The query starts in - sequences are off 1 bit
    }
	
    @query_starts = reverse @query_starts;
    @hit_starts = reverse @hit_starts;
    @sizes = reverse @sizes;
  } # end of reverse '-' strand blocks

  for (my $i=0;$i<scalar(@sizes);$i++) {
    my $exon_nr = $i+1;
    my $q_start = $query_starts[$i];
    my $h_start = $hit_starts[$i];
    my $q_end = $q_start + $sizes[$i] -1;
    my $h_end = $h_start + $sizes[$i] - 1;

    if ($current_items[$PSL_COL_LIST{strand}{COL_NO}] eq '-'){ #for - strand, swap start and end in hit
      my $tmp = $h_start;
      $h_start = $h_end;
      $h_end = $tmp;
    }

    $qi_exon->execute($current_line_id, $exon_nr, $q_start+1, $q_end+1, $h_start+1, $h_end+1);  

  } # end of insert each exon 

  $qi_psl -> finish;
  $qi_exon -> finish;
  
} # end of load_blat_result_old


#############################################
# psl sample:
#
#psLayout version 3
#
#match	mismatch 	repmatch 	NCount	Qgap	Qgap	Tgap	Tgap	strand	Qcount        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
#    	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#1509	0	0	0	0	0	5	29109	+	527|51013@B08|GI:51271	1509	0	1509	chr9	140273252	34210912	34241530	6	34,125,924,183,102,141,	0,34,159,1083,1266,1368,	34210912,34224213,34231182,34239776,34240655,34241389,
#						
#
#############################################

sub load_blat_result {
  my ($sql, $query, $rtn);
  my ($pkg, $tables, $rtInfo) = @_;

  my %PSL_COL_LIST = (
    matches => { #1509
      COL_NO => 0,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Number of bases that match that aren't repeats",
    },
    misMatches => { #0
      COL_NO => 1, 
      COL_TYPE => 'int unsigned',
      COL_DESC => "Number of bases that don't match",
    },
    repMatches => { #0
      COL_NO => 2,
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of bases that match but are part of repeats",
    },
    nCount => { #0
      COL_NO => 3, 
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of 'N' bases",
    },
    qNumInsert => { #0
      COL_NO => 4, 
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of inserts in query",
    },
    qBaseInsert => { #0
      COL_NO => 5,
      COL_TYPE => 'int unsigned' ,
      COL_TYPE => "Number of bases inserted in query",
    },
    tNumInsert => { #5
      COL_NO => 6,
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of inserts in target",
    },
    tBaseInsert => { #29109
      COL_NO => 7, 
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of bases inserted in target",
    },
    strand => { # +
      COL_NO => 8,
      COL_TYPE => 'char(2)' ,
      COL_DESC => "+ or - for query strand, optionally followed by + or - for target strand",
    },
    qName => { #527|51013@B08|GI:51271
      COL_NO => 9,
      COL_TYEP => 'varchar(255)',
      COL_DESC => "Query sequence name",
    },
    qSize => { #1509
      COL_NO => 10,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Query sequence size",
    },
    qStart => {#0
      COL_NO => 11,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Alignment start position in query",
    },
    qEnd => { #1509
      COL_NO => 12,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Alignment end position in query",
    },
    tName => { #chr9
      COL_NO => 13, 
      COL_TYPE => 'varchar(255)' ,
      COL_DESC => "Target sequence name",
    },
    tSize => {#	140273252
      COL_NO => 14,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Target sequence size",
    },
    tStart => {# 34210912
      COL_NO => 15,
      COL_TYPE => 'int unsigned',
      COL_DESC => "Alignment start position in target",
    },
    tEnd => { #34241530
      COL_NO => 16,
      COL_TYPE => 'int unsigned',
      COL_DESC => "alignment end position in target",
    },
    blockCount => { # 6
      COL_NO => 17, 
      COL_TYPE => 'int unsigned' ,
      COL_DESC => "Number of blocks in alignment",
    },
    blockSizes => { #34,125,924,183,102,141,
       COL_NO => 18,
       COL_TYPE => 'longblob' ,
       COL_DESC => "Size of each block in a comma separated list",
    },
    qStarts => { #0,34,159,1083,1266,1368,
       COL_NO => 19,
       COL_TYPE => 'longblob' ,
       COL_DESC => "Start of each block in query in a comma separated list",
    },
    tStarts => { #34210912,34224213,34231182,34239776,34240655,34241389,
      COL_NO => 20, 
      COL_TYPE => 'longblob' ,
      COL_DESC => "Start of each block in target in a comma separated list",
    },

  );

  $sql = "insert ignore into $tables->{BLAT}{TABLENAME} (LINE_ID, tName, qName, strand, match_nts, mismatch_nts, tSize, tStart, tEnd, qSize, qStart, qEnd, blockCount, qNumInsert, qBaseInsert, repMatches, nCount) values (" . "?, "x16 . "?)";
  my $qi_psl = $tables->{BLAT}{DBH} -> prepare($sql);

  $sql = "insert ignore into $tables->{BLAT_EXON}{TABLENAME} (LINE_ID, EXON_NR, qStart, qEnd, tStart, tEnd) values (" . "?, "x5 . "?)";
  my $qi_exon = $tables->{BLAT_EXON}{DBH} -> prepare($sql);

  open IN, "<",  $rtInfo->{RESULT_DIR} . $rtInfo->{BLAT_RESULT} or die $!;

  my @current_items = ();
  my $current_line_id = 0;
  my @items = ();
  my $score = 0;
#  my $current_line_id;
  my $line;
  my ($qName, $current_qName) = ('', '');

# skip the first 5 line
  for (1..5) { 
    $line = <IN>;
  }

# parse psl lines:

  while ($line = <IN>) {
    next if $line=~/chr.+_.+/; # skip those chromosome like 'chr9_random';
    chomp $line;
    @items = split /\t/, $line;

    ($qName) = $items[$PSL_COL_LIST{qName}{COL_NO}];
    if ($qName ne $current_qName) {#read a new id
      if ($current_qName ne '') {#if current values are not empty
        ($current_qName) = $current_items[$PSL_COL_LIST{qName}{COL_NO}];
        $qi_psl->execute(
                          $current_line_id,
                          $current_items[$PSL_COL_LIST{tName}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qName}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{strand}{COL_NO}],
                          $current_items[$PSL_COL_LIST{matches}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{misMatches}{COL_NO}],
                          $current_items[$PSL_COL_LIST{tSize}{COL_NO}],
                          $current_items[$PSL_COL_LIST{tStart}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{tEnd}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qSize}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qStart}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qEnd}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{blockCount}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qNumInsert}{COL_NO}],
                          $current_items[$PSL_COL_LIST{qBaseInsert}{COL_NO}], 
                          $current_items[$PSL_COL_LIST{repMatches}{COL_NO}],
                          $current_items[$PSL_COL_LIST{nCount}{COL_NO}]
                         );
        my @sizes = split /,/, $current_items[$PSL_COL_LIST{blockSizes}{COL_NO}];
        my @query_starts = split /,/, $current_items[$PSL_COL_LIST{qStarts}{COL_NO}];
        my @hit_starts = split /,/, $current_items[$PSL_COL_LIST{tStarts}{COL_NO}];

 #######IMPORTANT!!!!!!!!!!!!!!####################
 #Fix - strand hits so that they can be processed in the same way as the + Stand hits
 
       if ($current_items[$PSL_COL_LIST{strand}{COL_NO}] eq '-') {
	  for (my $i=0; $i<scalar(@sizes);$i++) {
	    $query_starts[$i] = $current_items[$PSL_COL_LIST{qSize}{COL_NO}] - $query_starts[$i] - $sizes[$i]; #The query starts in - sequences are off 1 bit
	  }
	
	  @query_starts = reverse @query_starts;
	  @hit_starts = reverse @hit_starts;
          @sizes = reverse @sizes;
        }

        for (my $i=0;$i<@sizes;$i++) {
	  my $exon_nr = $i+1;
	  my $q_start = $query_starts[$i];
	  my $h_start = $hit_starts[$i];
	  my $q_end = $q_start + $sizes[$i] -1;
	  my $h_end = $h_start + $sizes[$i] - 1;

	  if ($current_items[$PSL_COL_LIST{strand}{COL_NO}] eq '-'){ #for - strand, swap start and end in hit
	    my $tmp = $h_start;
	    $h_start = $h_end;
	    $h_end = $tmp;
	  }

	  $qi_exon->execute($current_line_id, $exon_nr, $q_start+1, $q_end+1, $h_start+1, $h_end+1); # set coordinates to be 1-based
	 # print "$exon_nr\t $q_start\t $q_end\t $hit_starts[$i]\t $h_end\t $current_qName\n";
        }

        $current_qName = $qName;
        @current_items = @items;
        $score = $items[$PSL_COL_LIST{matches}{COL_NO}]-$items[$PSL_COL_LIST{misMatches}{COL_NO}];
        $current_line_id = $.;
      }
      else { # first record
        @current_items = @items;
        $current_qName = $qName;
        $current_line_id = $.;
      }
    } # end of check new orf alignment
    elsif ($items[$PSL_COL_LIST{matches}{COL_NO}]-$items[$PSL_COL_LIST{misMatches}{COL_NO}] > $score) {#replace current items when the score is higher, so that the current values alway represent the highest score hit of the query
      $score = $items[$PSL_COL_LIST{matches}{COL_NO}]-$items[$PSL_COL_LIST{misMatches}{COL_NO}];
      @current_items = @items;
      $current_line_id = $.;
    }
  } # end of parse psl file

  close IN;

#Now load the last hit

  ($current_qName) = $current_items[$PSL_COL_LIST{qName}{COL_NO}];
   $qi_psl->execute(
                    $current_line_id,
                    $current_items[$PSL_COL_LIST{tName}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qName}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{strand}{COL_NO}],
                    $current_items[$PSL_COL_LIST{matches}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{misMatches}{COL_NO}],
                    $current_items[$PSL_COL_LIST{tSize}{COL_NO}],
                    $current_items[$PSL_COL_LIST{tStart}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{tEnd}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qSize}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qStart}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qEnd}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{blockCount}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qNumInsert}{COL_NO}],
                    $current_items[$PSL_COL_LIST{qBaseInsert}{COL_NO}], 
                    $current_items[$PSL_COL_LIST{repMatches}{COL_NO}],
                    $current_items[$PSL_COL_LIST{nCount}{COL_NO}]
                   );

  my @sizes = split /,/, $current_items[$PSL_COL_LIST{blockSizes}{COL_NO}];
  my @query_starts = split /,/, $current_items[$PSL_COL_LIST{qStarts}{COL_NO}];
  my @hit_starts = split /,/, $current_items[$PSL_COL_LIST{tStarts}{COL_NO}];

#######IMPORTANT!!!!!!!!!!!!!!####################
#Fix - strand hits so that they can be processed in the same way as the + Stand hits
  if ($current_items[$PSL_COL_LIST{strand}{COL_NO}] eq '-') {
    for (my $i=0; $i<scalar(@sizes);$i++) {
      $query_starts[$i] = $current_items[$PSL_COL_LIST{qSize}{COL_NO}] - $query_starts[$i] - $sizes[$i]; #The query starts in - sequences are off 1 bit
    }
	
    @query_starts = reverse @query_starts;
    @hit_starts = reverse @hit_starts;
    @sizes = reverse @sizes;
  } # end of reverse '-' strand blocks

  for (my $i=0;$i<scalar(@sizes);$i++) {
    my $exon_nr = $i+1;
    my $q_start = $query_starts[$i];
    my $h_start = $hit_starts[$i];
    my $q_end = $q_start + $sizes[$i] -1;
    my $h_end = $h_start + $sizes[$i] - 1;

    if ($current_items[$PSL_COL_LIST{strand}{COL_NO}] eq '-'){ #for - strand, swap start and end in hit
      my $tmp = $h_start;
      $h_start = $h_end;
      $h_end = $tmp;
    }

    $qi_exon->execute($current_line_id, $exon_nr, $q_start+1, $q_end+1, $h_start+1, $h_end+1);  

  } # end of insert each exon 

  $qi_psl -> finish;
  $qi_exon -> finish;
  
} # end of load_blat_result


############################################################
# run blat
############################################################
sub run_blat {
  my ($pkg, $t, $rtInfo) = @_;
  my ($sql, $query, $rtn);
  
  my $result_dir = (defined $rtInfo->{RESULT_DIR}) ? $rtInfo->{RESULT_DIR}:"";
  my $subject_dir = (defined $rtInfo->{SUBJECT_DIR}) ? $rtInfo->{SUBJECT_DIR}:"";
  my $bin_dir = (defined $rtInfo->{BIN_DIR}) ? $rtInfo->{BIN_DIR}:"";
  my $subject_fa = $rtInfo -> {SUBJECT_FASTA};
  my $query_fa = defined($rtInfo->{QUERY_FASTA}) ? $rtInfo->{QUERY_FASTA} : "query.fa";

# create temp fasta dir:
  my $fasta_dir = $result_dir . "fasta/";
  if(!-e $fasta_dir) {
    system("mkdir -p $fasta_dir");
  } 

# create blat dir:
  my $blat_dir = $result_dir . "blat/";
  if(!-e $blat_dir) {
    system("mkdir -p $blat_dir");
  } 

  my $script_dir = Cwd::getcwd();
  chdir($fasta_dir);

my $DO_QUERY_FASTA = 1;
if($DO_QUERY_FASTA==1) {
# get query fasta:
  $sql = "select $t->{QUERY}{QUERY_NAME_FIELD},$t->{QUERY}{QUERY_SEQ_FIELD} from $t->{QUERY}{TABLENAME} a ";
  $sql .= "where $t->{QUERY}{FILTER}" if defined($t->{QUERY}{FILTER});
  $query = $t->{QUERY}{DBH} -> prepare($sql);
  $rtn = $query -> execute();  
  open OUT, ">$query_fa";
  while (my ($query_name, $query_seq) = $query->fetchrow()) {
    print OUT ">$query_name";
    print OUT "\n";
    print OUT "$query_seq";
    print OUT "\n";
  } # end of while
  $query -> finish;
  close OUT;
}

# link subject db:
my $LINK_SUBJECT = 1;
if($LINK_SUBJECT == 1) {
  my $filelist = 'temp.list';
  system("find $subject_dir -name '" . $subject_fa . "*'>$filelist");
  open IN, "<$filelist";
  while(<IN>) {
    chomp;
    my ($fa) = $_ =~ /\/([^\/]+?)$/;
    system("ln -s $_ $fa"); 
  } # end while
  close IN;
  system("rm $filelist");
} # end of creating symlink


  MY_BioUtilitySeq -> enlist_blat_table($t, 'BLAT', $t->{QUERY}{DBH}, $t->{QUERY}{DBNAME}, $rtInfo->{RESULT_TABLE_PREFIX} . '_BLAT_BEST', 1, 1);
  MY_BioUtilitySeq -> enlist_blat_table($t, 'BLAT_EXON', $t->{QUERY}{DBH}, $t->{QUERY}{DBNAME}, $rtInfo->{RESULT_TABLE_PREFIX} . '_BLAT_EXON', 1, 1);
  MY_DButility -> prepare_tables($t);

# prepare blat:

my $DO_BLAT = 1;
if($DO_BLAT==1) {
  my $blat_result = "../blat/" . $query_fa . ".psl";
  system("blat -t=dna -q=dna -tileSize=8 -out=psl $subject_fa $query_fa $blat_result");
} # end of DO_BLAT

  chdir($script_dir);
  
# load blat result to db:
  my %blat_rtInfo = (
    RESULT_DIR => $blat_dir,
    BLAT_RESULT => $query_fa . ".psl", 
  );

  MY_SeqAlign -> load_blat_result($t, \%blat_rtInfo); 
  
} # end of run_blat


###############################
# New routine for load traces
# created on 20110517
###############################
sub load_raw_trace {
  my ($sql, $query, $rtn);
  my ($pkg, $plate_ids, $t) = @_;
  my @trace_plate_ids = @{$plate_ids};

# first delete all that is loaded:
  $sql = "delete from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(",", @trace_plate_ids) . ")";
  $rtn = $t->{TRACE_SEQ}{DBH} -> do($sql);
   
  $sql = "select trace_id, plate_id, plate_pos, fasta, qual_data, avg_qual_val, identity from $t->{TRACE_SEQ_ORIG}{TABLENAME} where plate_id in (" . join(",", @trace_plate_ids) 
   . ") "
#   . " and trace_id = 1192940 " # debug
#   . " and plate_id = 11946 and plate_pos = 'C02' " # debug, not reliable since plate_pos not correctly parsed
   . " order by plate_id, plate_pos";
  $query = $t->{TRACE_SEQ_ORIG}{DBH} -> prepare($sql);
  $sql = "insert into $t->{TRACE_SEQ}{TABLENAME} (TRACE_ID, TRACE_PLATE, TRACE_POS, TRACE_DIR, TRACE_SEQ, TRACE_QUAL, TRACE_AVG_QUAL, TRACE_IDENTITY) values (" . "?,"x7 . "?)";
  my $qi_trace = $t->{TRACE_SEQ}{DBH} -> prepare($sql);

  $rtn = $query -> execute();
  while (my ($trace_id, $plate_id, $plate_pos, $fasta, $qual, $avg_qual_val, $trace_identity) = $query -> fetchrow()) {
    $fasta =~ s/\n//g;
    $qual =~ s/\n//g;

#parse plate/pos and dir thru name: >NS_RC4_WTPlate1_24451_H03.scf 
 
    my ($trace_pos, $trace_seq) = $fasta =~ />.+_([A-H]\d+)\.scf(.*)\s*$/;
    my ($trace_qual) = $qual =~ />.+SCF(.*)\s*$/g;
    my ($trace_dir) = $fasta =~ />.+(FOR|FWD|REV|AD|DB|TERM|24668|24460)_[A-H]\d+\.scf.*\s*$/;
    if($trace_dir eq 'AD' || $trace_dir eq 'DB' || $trace_dir eq 'FOR' || $trace_dir eq 'FWD') {
      $trace_dir = 'FOR';
    }
    elsif ($trace_dir eq 'TERM' || $trace_dir eq 'REV') {
      $trace_dir = 'REV';
    }
    elsif ($trace_dir eq '24460') { # TERM
      $trace_dir = 'REV';
    }
    elsif ($trace_dir eq '24668') { # AD
      $trace_dir = 'FOR';
    }    
    else {
      $trace_dir = 'FOR'; #default dir='FOR' if not specified
    }
    $rtn = $qi_trace -> execute($trace_id, $plate_id, $trace_pos, $trace_dir, $trace_seq, $trace_qual, $avg_qual_val, $trace_identity);

  } # end of load-traces

  $query -> finish;
  $qi_trace -> finish;

} # end of load_raw_trace


#*************************************************************************
#* END of subroutines
#*************************************************************************

###################################################
# General QC pipeline
# %t must contain the definition of the following tables:
#    CHERRYMAP,
#    PLATE_DETAIL,
#    TRACE_SEQ_ORIG,
#    TRACE_SEQ,
#    REFSEQ,
#    
###################################################
sub QC_helper_new {
  my ($pkg, $t, $qc_rtInfo) = @_;
  my @trace_plate_ids = @{$qc_rtInfo->{TRACE_PLATE_IDS}};
  my @trimseq_ids = defined($qc_rtInfo->{TRIM_SEQ_IDS}) ? @{$qc_rtInfo->{TRIM_SEQ_IDS}} : (3);
  my ($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time);
  my $curdate = sprintf("%04d%02d%02d", $Year+1900, $Month+1, $Day);
  my $result_dir = defined($qc_rtInfo->{RESULT_DIR}) ? $qc_rtInfo->{RESULT_DIR} : "result/QC/" . $qc_rtInfo->{SOURCE} . "_" . $curdate . "/";

  if(!-e $result_dir) {
    system("mkdir -p $result_dir");
  }

  my ($sql, $query, $rtn);
  
# MY_DButility -> prepare_tables($t); # assume tables are all prepared before hand.

# STEP 1: LOAD TRACE SEQ:

  my $LOAD_TRACE = defined($qc_rtInfo->{LOAD_TRACE})? $qc_rtInfo->{LOAD_TRACE}:0;
  
  if ($LOAD_TRACE==1) {

# first delete all that is loaded:
    $sql = "delete from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in ('" . join("','", @trace_plate_ids) . "')"; # treat trace plate as string
    $rtn = $t->{TRACE_SEQ}{DBH} -> do($sql);
    $sql = "select trace_id, trace_name, plate_name, plate_id, plate_pos, fasta, qual_data, avg_qual_val, identity from $t->{TRACE_SEQ_ORIG}{TABLENAME} inner join tbPlate b using(plate_id) where plate_id in ('" . join("','", @trace_plate_ids) 
      . "') "
#      . " and trace_id = 1192940 " # debug
#      . " and plate_id = 11946 and plate_pos = 'C02' " # debug, not reliable since plate_pos not correctly parsed
      . " order by plate_id, plate_pos";
    $query = $t->{TRACE_SEQ_ORIG}{DBH} -> prepare($sql);
    $sql = "insert into $t->{TRACE_SEQ}{TABLENAME} (TRACE_ID, TRACE_NAME, TRACE_PLATE_NAME,  TRACE_PLATE, TRACE_POS, TRACE_DIR, TRACE_SEQ, TRACE_QUAL, TRACE_AVG_QUAL, TRACE_IDENTITY) values (" . "?,"x9 . "?)";
    my $qi_trace = $t->{TRACE_SEQ}{DBH} -> prepare($sql);

    $rtn = $query -> execute();
    while (my ($trace_id, $trace_name, $trace_plate_name, $plate_id, $plate_pos, $fasta, $qual, $avg_qual_val, $trace_identity) = $query -> fetchrow()) {
      $fasta =~ s/\n//g;
      $qual =~ s/\n//g;

#parse plate/pos and dir thru name: >NS_RC4_WTPlate1_24451_H03.scf 
 
      my ($trace_pos, $trace_seq) = $fasta =~ />.+_([A-H]\d+)\.scf(.*)\s*$/;
      my ($trace_qual) = $qual =~ />.+SCF(.*)\s*$/g;
      my ($trace_dir) = $fasta =~ />.+(FOR|FWD|FORWARD|REV|AD|DB|TERM|24773|ADHprom|Mdm2chkR).*_[A-H]\d+\.scf.*\s*$/i;
      $trace_dir = uc($trace_dir);
      if($trace_dir eq 'AD' || $trace_dir eq 'DB' || $trace_dir eq 'FOR' || $trace_dir eq 'FWD' || $trace_dir eq 'ADHPROM') {
        $trace_dir = 'FOR';
      }
      elsif ($trace_dir eq 'TERM' || $trace_dir eq 'REV'|| $trace_dir eq uc('Mdm2chkR')) {
        $trace_dir = 'REV';
      }
      elsif ($trace_dir eq '24755') {
        $trace_dir = 'FOR';
      }
      elsif ($trace_dir eq '24756') {
        $trace_dir = 'REV';
      }
      else {
        $trace_dir = 'FOR'; #default dir='FOR' if not specified
      }
      $rtn = $qi_trace -> execute($trace_id, $trace_name, $trace_plate_name, $plate_id, $trace_pos, $trace_dir, $trace_seq, $trace_qual, $avg_qual_val, $trace_identity);

  } # end of load-traces

  $query -> finish;
  $qi_trace -> finish;

# update trace plate name: 
  

} # end of load traces

# step2: trim vector sequence:
# do gateway tail trim 
#  my $VECTOR_TRIM = 0;
  my $VECTOR_TRIM = defined($qc_rtInfo->{VECTOR_TRIM})? $qc_rtInfo->{VECTOR_TRIM}:0;

if($VECTOR_TRIM==1) { 
=pod
  my %vtrim_tables = (
    VECTOR => {
      DBH => $t->{TRACE_SEQ}{DBH},
      DBNAME => $t->{TRACE_SEQ}{DBNAME},
      TABLENAME => "(select SEQ_NAME, DIR, right(SEQ, 70) as VECTOR_SEQ from bio_utility.g_special_seq where SEQ_ID in (" . join(",", @trimseq_ids) . "))", #only FOR-M13-G used
    }, 

    TRACE_SEQ => {
      DBH => $t->{TRACE_SEQ}{DBH},
      DBNAME => $t->{TRACE_SEQ}{DBNAME},
      TABLENAME => $t->{TRACE_SEQ}{TABLENAME},
      FILTER => "TRACE_PLATE in ('" . join("','",@trace_plate_ids) . "')", 
 
      VTRIM_START_FIELD => 'VTRIM_START',
      VTRIM_END_FIELD => 'VTRIM_END',
    }, 
  
  );
=cut
    my @trimseq_ids = defined($qc_rtInfo->{TRIM_SEQ_IDS}) ? @{$qc_rtInfo->{TRIM_SEQ_IDS}} : (3);
    my %vtrim_tables = ();
    if(defined($t->{VECTOR})) {
       $vtrim_tables{VECTOR} = $t->{VECTOR};
    } 
    else {
      $vtrim_tables{VECTOR} = {
	  DBH => $t->{TRACE_SEQ}{DBH},
	  DBNAME => $t->{TRACE_SEQ}{DBNAME},
	  TABLENAME => "(select SEQ_NAME, DIR, right(SEQ, 70) as VECTOR_SEQ from bio_utility.g_special_seq where SEQ_ID in (" . join(",", @trimseq_ids) . "))", #only FOR-M13-G used
      }; 
    }

      $vtrim_tables{TRACE_SEQ} = {
	  DBH => $t->{TRACE_SEQ}{DBH},
	  DBNAME => $t->{TRACE_SEQ}{DBNAME},
	  TABLENAME => $t->{TRACE_SEQ}{TABLENAME},
	  FILTER => "TRACE_PLATE in ('" . join("','",@trace_plate_ids) . "')",
	  VTRIM_START_FIELD => 'VTRIM_START',
	  VTRIM_END_FIELD => 'VTRIM_END',
      }; 

  my %vtrim_rtinfo = (
    RESULT_DIR => $result_dir,
    TMP_DIR => $result_dir . 'VTRIM/',
    VECTOR_FASTA => 'gal4vector.fa',
    TRACE_FASTA => "qc_" . $qc_rtInfo->{SOURCE} .".fa", 
  );

  MY_SeqAlign -> run_vector_trim(\%vtrim_tables, \%vtrim_rtinfo);

} # end of vector trimming

# NOW DO QUAL TRIM
  my $QUAL_TRIM = defined($qc_rtInfo->{QUAL_TRIM})? $qc_rtInfo->{QUAL_TRIM}:0;

if ($QUAL_TRIM==1) {
  my %qtrim_tables = (
    TRACE_SEQ => $t->{TRACE_SEQ},
  );
  $qtrim_tables{TRACE_SEQ}{TRACE_ID_FIELD} = 'TRACE_ID';
  $qtrim_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = 'TRACE_SEQ';
  $qtrim_tables{TRACE_SEQ}{TRACE_QUAL_FIELD} = 'TRACE_QUAL';
  $qtrim_tables{TRACE_SEQ}{TRACE_START_FIELD} = 'VTRIM_START'; 
  $qtrim_tables{TRACE_SEQ}{HIGH_QUAL_START_FIELD} = 'HQ_START';
  $qtrim_tables{TRACE_SEQ}{HIGH_QUAL_END_FIELD} = 'HQ_END';
  $qtrim_tables{FILTER}{QUAL_TRIM} = "TRACE_PLATE in ('" . join("','", @trace_plate_ids) . "')"; #chg to string type
#  $qtrim_tables{FILTER}{QUAL_TRIM} .= " and TRACE_ID=1153342"; # debug
#  $qtrim_tables{FILTER}{QUAL_TRIM} .= " and TRACE_PLATE=11687 and TRACE_POS='A06'"; # debug

# quality trim:
  MY_SeqAlign -> general_qualtrim(\%qtrim_tables, $qc_rtInfo);

} # end of qual trimming

# run bl2seq:

  my $BL2SEQ = defined($qc_rtInfo->{BL2SEQ})? $qc_rtInfo->{BL2SEQ}:0;
  my $BL2SEQ_W_BRIEF_BLAST = defined($qc_rtInfo->{BL2SEQ_W_BRIEF_BLAST})? $qc_rtInfo->{BL2SEQ_W_BRIEF_BLAST}:0;
  my $BL2SEQ_W_FULL_BLAST = defined($qc_rtInfo->{BL2SEQ_W_FULL_BLAST})? $qc_rtInfo->{BL2SEQ_W_FULL_BLAST}:0;

if($BL2SEQ == 1) {
  my $HQ_SEQ_CUTOFF = defined($qc_rtInfo->{HQ_SEQ_CUTOFF})?$qc_rtInfo->{HQ_SEQ_CUTOFF}:550; #set the cutoff for High Quality sequence by the intrinsic seq. error tolerance 

  my %bl2_tables = (
    TRACE_SEQ => $t->{TRACE_SEQ},
    REFSEQ => $t->{REFSEQ},
    CHERRYMAP => $t->{CHERRYMAP},    
  );

  $bl2_tables{REFSEQ}{REFSEQ_ID_FIELD} = defined $t->{REFSEQ}{REFSEQ_ID_FIELD}?$t->{REFSEQ}{REFSEQ_ID_FIELD}:"ORF_ID";
  $bl2_tables{REFSEQ}{REFSEQ_FIELD} = defined $t->{REFSEQ}{REFSEQ_FIELD} ? $t->{REFSEQ}{REFSEQ_FIELD}:'REFSEQ'; # need to use alias if not exactly same;

  $bl2_tables{TRACE_SEQ}{TRACE_ID_FIELD} = 'TRACE_ID';
  $bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD} = 'TRACE_PLATE';
  $bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD} = 'TRACE_POS';
#  $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = "if(HQ_START is null, TRACE_SEQ, substr(TRACE_SEQ, HQ_START, if(HQ_END<$HQ_SEQ_CUTOFF, HQ_END,$HQ_SEQ_CUTOFF)-HQ_START+1))"; # if no QUAL trim is found, use whole seq. 
  if($qc_rtInfo->{USE_RAW_TRACE} == 1) {
    $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = 'TRACE_SEQ';
  }
  else {
    if($HQ_SEQ_CUTOFF != 0) {
      $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = "substr(TRACE_SEQ, HQ_START, if(HQ_END<$HQ_SEQ_CUTOFF, HQ_END,$HQ_SEQ_CUTOFF)-HQ_START+1)";
    }
    else {
      $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = "substr(TRACE_SEQ, HQ_START, HQ_END-HQ_START+1)";
    }
  }
  $bl2_tables{TRACE_SEQ}{BL2_SCORE_FIELD} = 'BL2_SCORE';
  $bl2_tables{TRACE_SEQ}{INSERT_START_FIELD} = 'INSERT_START';
  $bl2_tables{TRACE_SEQ}{INSERT_END_FIELD} = 'INSERT_END';
  $bl2_tables{TRACE_SEQ}{ALIGN_LEN_FIELD} = 'ALIGN_LEN';
  $bl2_tables{TRACE_SEQ}{IDENTITY_FIELD} = 'IDENTITY';
  $bl2_tables{TRACE_SEQ}{SEQ_OFFSET_FIELD} = ($qc_rtInfo->{USE_RAW_TRACE}==1)?1:"HQ_START"; 

  if ($BL2SEQ_W_BRIEF_BLAST==1) {
     $bl2_tables{TRACE_SEQ}{SCORE_FIELD}= 'SCORE';
     $bl2_tables{TRACE_SEQ}{BITS_FIELD}= 'BITS';
     $bl2_tables{TRACE_SEQ}{numHSP_FIELD}= 'numHSP';
     $bl2_tables{TRACE_SEQ}{qLenAlign_FIELD}= 'qLenAlign';
     $bl2_tables{TRACE_SEQ}{sLenAlign_FIELD}= 'sLenAlign';
     $bl2_tables{TRACE_SEQ}{qStart_FIELD}= 'qStart';
     $bl2_tables{TRACE_SEQ}{qEnd_FIELD}= 'qEnd';
     $bl2_tables{TRACE_SEQ}{sStart_FIELD}= 'sStart';
     $bl2_tables{TRACE_SEQ}{sEnd_FIELD}= 'sEnd';
     $bl2_tables{TRACE_SEQ}{qNumUnaligned_FIELD}= 'qNumUnaligned';
     $bl2_tables{TRACE_SEQ}{sNumUnaligned_FIELD}= 'sNumUnaligned';
     $bl2_tables{TRACE_SEQ}{qSeqIndIdentical_FIELD}= 'qSeqIndIdentical';
     $bl2_tables{TRACE_SEQ}{sSeqIndIdentical_FIELD}= 'sSeqIndIdentical';
     $bl2_tables{TRACE_SEQ}{qStrand_FIELD}= 'qStrand';
     $bl2_tables{TRACE_SEQ}{sStrand_FIELD}= 'sStrand';
     $bl2_tables{TRACE_SEQ}{RANK_FIELD}= 'RANK';
     $bl2_tables{TRACE_SEQ}{qSeqLenIdentical_FIELD} = 'qSeqLenIdentical';
     $bl2_tables{TRACE_SEQ}{sSeqLenIdentical_FIELD} = 'sSeqLenIdentical';
     $bl2_tables{TRACE_SEQ}{qSeqGap_FIELD} = 'qSeqGap';
     $bl2_tables{TRACE_SEQ}{sSeqGap_FIELD} = 'sSeqGap';
  }
  
# put the filter in the single table:
  $bl2_tables{TRACE_SEQ}{FILTER} = "$bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD} in ('" . join("','",@trace_plate_ids) . "') ";
  $bl2_tables{TRACE_SEQ}{FILTER} .= " and HQ_START>-1" if ($qc_rtInfo->{USE_RAW_TRACE} != 1) ; # only fetch those passing quality trim

#  $bl2_tables{TRACE_SEQ}{FILTER} .= " and $bl2_tables{TRACE_SEQ}{TRACE_ID_FIELD}=1196305"; # debug
#  $bl2_tables{TRACE_SEQ}{FILTER} .= "  $bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}=205003 and $bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD}='A11'"; # debug
#    $bl2_tables{TRACE_SEQ}{FILTER} .= " $bl2_tables{TRACE_SEQ}{TRACE_ID_FIELD} not in (1793, 521, 582, 767, 1438)"; # debug
  
  $bl2_tables{CHERRYMAP}{SRC_PLATE_FIELD} = 'srcPlate';
  $bl2_tables{CHERRYMAP}{SRC_POS_FIELD} = 'srcPos';
  $bl2_tables{CHERRYMAP}{DST_PLATE_FIELD} = 'dstPlate';
  $bl2_tables{CHERRYMAP}{DST_POS_FIELD} = 'dstPos';
  $bl2_tables{CHERRYMAP}{REFSEQ_ID_FIELD} = 'ORF_ID';
  
  $bl2_tables{JOIN_CONDITION}{TRACE_CHERRY} = "a.$bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}=b.$bl2_tables{CHERRYMAP}{DST_PLATE_FIELD} and a.$bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD}=b.$bl2_tables{CHERRYMAP}{DST_POS_FIELD}";

  $bl2_tables{JOIN_CONDITION}{CHERRY_REFSEQ} = "b.$bl2_tables{CHERRYMAP}{REFSEQ_ID_FIELD}=c.$bl2_tables{REFSEQ}{REFSEQ_ID_FIELD}";

  $bl2_tables{FILTER}{BL2SEQ} = " bl2_score is null"; # use to put joined condition
 
  $bl2_tables{ORDER_BY}{BL2SEQ} = "$bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}, $bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD}";

  my $bl2seq_type = 0;
  if($BL2SEQ_W_BRIEF_BLAST==1) {
      $bl2seq_type = 1;
  }
  MY_SeqAlign -> general_bl2seq($bl2_tables{TRACE_SEQ}{DBH}, \%bl2_tables, $bl2seq_type);

# catch those didn't pass quality/vector trim by using raw trace:

  if($qc_rtInfo->{USE_RAW_TRACE} != 1) {
    $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = 'TRACE_SEQ';
    $bl2_tables{TRACE_SEQ}{SEQ_OFFSET_FIELD} = 1;
# put the filter in the single table:
    $bl2_tables{TRACE_SEQ}{FILTER} = "$bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD} in ('" . join("','",@trace_plate_ids) . "') ";
    $bl2_tables{TRACE_SEQ}{FILTER} .= " and $bl2_tables{TRACE_SEQ}{BL2_SCORE_FIELD} is null"; # only fetch those didn't align first time 
#    $bl2_tables{TRACE_SEQ}{FILTER} .= " and $bl2_tables{TRACE_SEQ}{TRACE_ID_FIELD} not in (1793)"; # debug
#  $bl2_tables{TRACE_SEQ}{FILTER} .= " and $bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}=11687 and $bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD}='A06'"; # debug

#rerun bl2seq using raw trace seq    
     MY_SeqAlign -> general_bl2seq($bl2_tables{TRACE_SEQ}{DBH}, \%bl2_tables, $bl2seq_type);
  }
} # end bl2seq


} # end QC_helper_new



###################################################
# General QC pipeline
# %t must contain the definition of the following tables:
#    CHERRYMAP,
#    PLATE_DETAIL,
#    TRACE_SEQ_ORIG,
#    TRACE_SEQ,
#    REFSEQ,
#    
###################################################
sub QC_helper {
  my ($pkg, $qc_rtInfo, $source, $plate_ids, $t) = @_;
  my @trace_plate_ids = @{$plate_ids};
  my @trimseq_ids = defined($qc_rtInfo->{TRIM_SEQ_IDS}) ? @{$qc_rtInfo->{TRIM_SEQ_IDS}} : (3);
  my ($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time);
  my $curdate = sprintf("%04d%02d%02d", $Year+1900, $Month+1, $Day);
  my $result_dir = defined($qc_rtInfo->{RESULT_DIR}) ? $qc_rtInfo->{RESULT_DIR} : "result/QC/" . $source . "_" . $curdate . "/";

  if(!-e $result_dir) {
    system("mkdir -p $result_dir");
  }

  my ($sql, $query, $rtn);
  
# MY_DButility -> prepare_tables($t);

# LOAD TRACE SEQ:

#  my $LOAD_TRACE = 1;
  my $LOAD_TRACE = defined($qc_rtInfo->{LOAD_TRACE})? $qc_rtInfo->{LOAD_TRACE}:0;
  
if ($LOAD_TRACE==1) {
# first delete all that is loaded:
  $sql = "delete from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in ('" . join("','", @trace_plate_ids) . "')";
  $rtn = $t->{TRACE_SEQ}{DBH} -> do($sql);
   
  $sql = "select trace_id, plate_id, plate_pos, fasta, qual_data, avg_qual_val, identity from $t->{TRACE_SEQ_ORIG}{TABLENAME} where plate_id in ('" . join("','", @trace_plate_ids) 
   . "') "
#   . " and trace_id = 1289057 " # debug
#   . " and plate_id = 11946 and plate_pos = 'C02' " # debug, not reliable since plate_pos not correctly parsed
   . " order by plate_id, plate_pos";
  $query = $t->{TRACE_SEQ_ORIG}{DBH} -> prepare($sql);
  $sql = "insert into $t->{TRACE_SEQ}{TABLENAME} (TRACE_ID, TRACE_PLATE, TRACE_POS, TRACE_DIR, TRACE_SEQ, TRACE_QUAL, TRACE_AVG_QUAL, TRACE_IDENTITY) values (" . "?,"x7 . "?)";
  my $qi_trace = $t->{TRACE_SEQ}{DBH} -> prepare($sql);

  $rtn = $query -> execute();
  while (my ($trace_id, $plate_id, $plate_pos, $fasta, $qual, $avg_qual_val, $trace_identity) = $query -> fetchrow()) {
    $fasta =~ s/\n//g;
    $qual =~ s/\n//g;

#parse plate/pos and dir thru name: >NS_RC4_WTPlate1_24451_H03.scf 
 
    my ($trace_pos, $trace_seq) = $fasta =~ />.+_([A-H]\d+)\.scf(.*)\s*$/;
    my ($trace_qual) = $qual =~ />.+SCF(.*)\s*$/g;
    my ($trace_dir) = $fasta =~ />.+(FOR|FWD|FORWARD|REV|AD|DB|TERM|24773|ADHprom|Mdm2chkR|REVERSE).*_[A-H]\d+\.scf.*\s*$/i;
    $trace_dir = uc($trace_dir);
    if($trace_dir eq 'AD' || $trace_dir eq 'DB' || $trace_dir eq 'FOR' || $trace_dir eq 'FWD' || $trace_dir eq 'ADHPROM' || $trace_dir eq 'FORWARD'  ) {
      $trace_dir = 'FOR';
    }
    elsif ($trace_dir eq 'TERM' || $trace_dir eq 'REV'|| $trace_dir eq uc('Mdm2chkR') || $trace_dir eq 'REVERSE') {
      $trace_dir = 'REV';
    }
=pod
    elsif ($trace_dir eq '24755') {
      $trace_dir = 'FOR';
    }
    elsif ($trace_dir eq '24756') {
      $trace_dir = 'REV';
    }
    else {
      $trace_dir = 'FOR'; #default dir='FOR' if not specified
    }
=cut
    $rtn = $qi_trace -> execute($trace_id, $plate_id, $trace_pos, $trace_dir, $trace_seq, $trace_qual, $avg_qual_val, $trace_identity);

  } # end of load-traces

  $query -> finish;
  $qi_trace -> finish;

} # end of load traces

# do gateway tail trim 
#  my $VECTOR_TRIM = 0;
  my $VECTOR_TRIM = defined($qc_rtInfo->{VECTOR_TRIM})? $qc_rtInfo->{VECTOR_TRIM}:0;

if($VECTOR_TRIM==1) {  
=pod
  my %vtrim_tables = (
    VECTOR => {
      DBH => $t->{TRACE_SEQ}{DBH},
      DBNAME => $t->{TRACE_SEQ}{DBNAME},
      TABLENAME => "(select SEQ_NAME, DIR, right(SEQ, 70) as VECTOR_SEQ from bio_utility.g_special_seq where SEQ_ID in (" . join(",", @trimseq_ids) . "))", #only FOR-M13-G used
    }, 

    TRACE_SEQ => {
      DBH => $t->{TRACE_SEQ}{DBH},
      DBNAME => $t->{TRACE_SEQ}{DBNAME},
      TABLENAME => $t->{TRACE_SEQ}{TABLENAME},
      FILTER => "TRACE_PLATE in ('" . join("','",@trace_plate_ids) . "')", 
 
      VTRIM_START_FIELD => 'VTRIM_START',
      VTRIM_END_FIELD => 'VTRIM_END',
    }, 
  
  );
=cut
# .  " and trace_id=1349998 ", " and TRACE_PLATE=13269 and TRACE_POS='H10'", . " and TRACE_DIR = 'REV'",
    my %vtrim_tables = ();
    if(defined($t->{VECTOR})) {
       $vtrim_tables{VECTOR} = $t->{VECTOR};
    } 
    else {
      $vtrim_tables{VECTOR} = {
	  DBH => $t->{TRACE_SEQ}{DBH},
	  DBNAME => $t->{TRACE_SEQ}{DBNAME},
	  TABLENAME => "(select SEQ_NAME, DIR, right(SEQ, 70) as VECTOR_SEQ from bio_utility.g_special_seq where SEQ_ID in (" . join(",", @trimseq_ids) . "))", #only FOR-M13-G used
      }; 
    }

      $vtrim_tables{TRACE_SEQ} = {
	  DBH => $t->{TRACE_SEQ}{DBH},
	  DBNAME => $t->{TRACE_SEQ}{DBNAME},
	  TABLENAME => $t->{TRACE_SEQ}{TABLENAME},
	  FILTER => "TRACE_PLATE in ('" . join("','",@trace_plate_ids) . "')", #.  " and trace_id=1349998 ", #" and TRACE_PLATE=13269 and TRACE_POS='H10'", # . " and TRACE_DIR = 'REV'",
	  VTRIM_START_FIELD => 'VTRIM_START',
	  VTRIM_END_FIELD => 'VTRIM_END',
      }; 

  my %vtrim_rtinfo = (
    RESULT_DIR => $result_dir,
    TMP_DIR => $result_dir . 'VTRIM/',
    VECTOR_FASTA => 'gal4vector.fa',
    TRACE_FASTA => "qc_" . $source .".fa", 
  );

  MY_SeqAlign -> run_vector_trim(\%vtrim_tables, \%vtrim_rtinfo);

} # end of vector trimming


# NOW DO QUAL TRIM
#  my $QUAL_TRIM = 1;
  my $QUAL_TRIM = defined($qc_rtInfo->{QUAL_TRIM})? $qc_rtInfo->{QUAL_TRIM}:0;

if ($QUAL_TRIM==1) {
  my %qtrim_tables = (
    TRACE_SEQ => $t->{TRACE_SEQ},
  );
  $qtrim_tables{TRACE_SEQ}{TRACE_ID_FIELD} = 'TRACE_ID';
  $qtrim_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = 'TRACE_SEQ';
  $qtrim_tables{TRACE_SEQ}{TRACE_QUAL_FIELD} = 'TRACE_QUAL';
  $qtrim_tables{TRACE_SEQ}{TRACE_START_FIELD} = 'VTRIM_START'; 
  $qtrim_tables{TRACE_SEQ}{HIGH_QUAL_START_FIELD} = 'HQ_START';
  $qtrim_tables{TRACE_SEQ}{HIGH_QUAL_END_FIELD} = 'HQ_END';
  $qtrim_tables{FILTER}{QUAL_TRIM} = "TRACE_PLATE in ('" . join("','", @trace_plate_ids) . "')"; #chg to string type
#  $qtrim_tables{FILTER}{QUAL_TRIM} .= " and TRACE_ID=1153342"; # debug
#  $qtrim_tables{FILTER}{QUAL_TRIM} .= " and TRACE_PLATE=11687 and TRACE_POS='A06'"; # debug

# quality trim:
  MY_SeqAlign -> general_qualtrim(\%qtrim_tables, $qc_rtInfo);

} # end of qual trimming

# run bl2seq:
#  my $BL2SEQ = 1;
  my $BL2SEQ = defined($qc_rtInfo->{BL2SEQ})? $qc_rtInfo->{BL2SEQ}:0;
  my $BL2SEQ_W_BRIEF_BLAST = defined($qc_rtInfo->{BL2SEQ_W_BRIEF_BLAST})? $qc_rtInfo->{BL2SEQ_W_BRIEF_BLAST}:0;
  my $BL2SEQ_W_FULL_BLAST = defined($qc_rtInfo->{BL2SEQ_W_FULL_BLAST})? $qc_rtInfo->{BL2SEQ_W_FULL_BLAST}:0;

if($BL2SEQ == 1) {
# || $BL2SEQ_W_BRIEF_BLAST == 1 || $BL2SEQ_W_FULL_BLAST == 1) {
  my $HQ_SEQ_CUTOFF = defined($qc_rtInfo->{HQ_SEQ_CUTOFF})?$qc_rtInfo->{HQ_SEQ_CUTOFF}:550; #set the cutoff for High Quality sequence by the intrinsic seq. error tolerance 

  my %bl2_tables = (
    TRACE_SEQ => $t->{TRACE_SEQ},
    REFSEQ => $t->{REFSEQ},
    CHERRYMAP => $t->{CHERRYMAP},    
  );

  $bl2_tables{REFSEQ}{REFSEQ_ID_FIELD} = defined $t->{REFSEQ}{REFSEQ_ID_FIELD}?$t->{REFSEQ}{REFSEQ_ID_FIELD}:"ORF_ID";
  $bl2_tables{REFSEQ}{REFSEQ_FIELD} = defined $t->{REFSEQ}{REFSEQ_FIELD} ? $t->{REFSEQ}{REFSEQ_FIELD}:'REFSEQ'; # need to use alias if not exactly same;

  $bl2_tables{TRACE_SEQ}{TRACE_ID_FIELD} = 'TRACE_ID';
  $bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD} = 'TRACE_PLATE';
  $bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD} = 'TRACE_POS';
#  $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = "if(HQ_START is null, TRACE_SEQ, substr(TRACE_SEQ, HQ_START, if(HQ_END<$HQ_SEQ_CUTOFF, HQ_END,$HQ_SEQ_CUTOFF)-HQ_START+1))"; # if no QUAL trim is found, use whole seq. 
  if($qc_rtInfo->{USE_RAW_TRACE} == 1) {
    $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = 'TRACE_SEQ';
  }
  else {
    if($HQ_SEQ_CUTOFF != 0) {
      $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = "substr(TRACE_SEQ, HQ_START, if(HQ_END<$HQ_SEQ_CUTOFF, HQ_END,$HQ_SEQ_CUTOFF)-HQ_START+1)";
    }
    else {
      $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = "substr(TRACE_SEQ, HQ_START, HQ_END-HQ_START+1)";
    }
  }
  $bl2_tables{TRACE_SEQ}{BL2_SCORE_FIELD} = 'BL2_SCORE';
  $bl2_tables{TRACE_SEQ}{INSERT_START_FIELD} = 'INSERT_START';
  $bl2_tables{TRACE_SEQ}{INSERT_END_FIELD} = 'INSERT_END';
  $bl2_tables{TRACE_SEQ}{ALIGN_LEN_FIELD} = 'ALIGN_LEN';
  $bl2_tables{TRACE_SEQ}{IDENTITY_FIELD} = 'IDENTITY';
  $bl2_tables{TRACE_SEQ}{SEQ_OFFSET_FIELD} = ($qc_rtInfo->{USE_RAW_TRACE}==1)?1:"HQ_START"; 

  if ($BL2SEQ_W_BRIEF_BLAST==1) {
     $bl2_tables{TRACE_SEQ}{SCORE_FIELD}= 'SCORE';
     $bl2_tables{TRACE_SEQ}{BITS_FIELD}= 'BITS';
     $bl2_tables{TRACE_SEQ}{numHSP_FIELD}= 'numHSP';
     $bl2_tables{TRACE_SEQ}{qLenAlign_FIELD}= 'qLenAlign';
     $bl2_tables{TRACE_SEQ}{sLenAlign_FIELD}= 'sLenAlign';
     $bl2_tables{TRACE_SEQ}{qStart_FIELD}= 'qStart';
     $bl2_tables{TRACE_SEQ}{qEnd_FIELD}= 'qEnd';
     $bl2_tables{TRACE_SEQ}{sStart_FIELD}= 'sStart';
     $bl2_tables{TRACE_SEQ}{sEnd_FIELD}= 'sEnd';
     $bl2_tables{TRACE_SEQ}{qNumUnaligned_FIELD}= 'qNumUnaligned';
     $bl2_tables{TRACE_SEQ}{sNumUnaligned_FIELD}= 'sNumUnaligned';
     $bl2_tables{TRACE_SEQ}{qSeqIndIdentical_FIELD}= 'qSeqIndIdentical';
     $bl2_tables{TRACE_SEQ}{sSeqIndIdentical_FIELD}= 'sSeqIndIdentical';
     $bl2_tables{TRACE_SEQ}{qStrand_FIELD}= 'qStrand';
     $bl2_tables{TRACE_SEQ}{sStrand_FIELD}= 'sStrand';
     $bl2_tables{TRACE_SEQ}{RANK_FIELD}= 'RANK';
     $bl2_tables{TRACE_SEQ}{qSeqLenIdentical_FIELD} = 'qSeqLenIdentical';
     $bl2_tables{TRACE_SEQ}{sSeqLenIdentical_FIELD} = 'sSeqLenIdentical';
     $bl2_tables{TRACE_SEQ}{qSeqGap_FIELD} = 'qSeqGap';
     $bl2_tables{TRACE_SEQ}{sSeqGap_FIELD} = 'sSeqGap';
  }
  
# put the filter in the single table:
#  $bl2_tables{TRACE_SEQ}{FILTER} = "$bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD} in ('" . join("','",@trace_plate_ids) . "') ";
  $bl2_tables{TRACE_SEQ}{FILTER} .= " and HQ_START>-1" if ($qc_rtInfo->{USE_RAW_TRACE} != 1) ; # only fetch those passing quality trim

#  $bl2_tables{TRACE_SEQ}{FILTER} .= " and $bl2_tables{TRACE_SEQ}{TRACE_ID_FIELD}=1196305"; # debug
#  $bl2_tables{TRACE_SEQ}{FILTER} .= " and $bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}=205003 and $bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD}='A11'"; # debug
#    $bl2_tables{TRACE_SEQ}{FILTER} .= " $bl2_tables{TRACE_SEQ}{TRACE_ID_FIELD} not in (1793, 521, 582, 767, 1438)"; # debug
  
  $bl2_tables{CHERRYMAP}{SRC_PLATE_FIELD} = 'srcPlate';
  $bl2_tables{CHERRYMAP}{SRC_POS_FIELD} = 'srcPos';
  $bl2_tables{CHERRYMAP}{DST_PLATE_FIELD} = 'dstPlate';
  $bl2_tables{CHERRYMAP}{DST_POS_FIELD} = 'dstPos';
  $bl2_tables{CHERRYMAP}{REFSEQ_ID_FIELD} = 'ORF_ID';
  
  $bl2_tables{JOIN_CONDITION}{TRACE_CHERRY} = "a.$bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}=b.$bl2_tables{CHERRYMAP}{DST_PLATE_FIELD} and a.$bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD}=b.$bl2_tables{CHERRYMAP}{DST_POS_FIELD}";

  $bl2_tables{JOIN_CONDITION}{CHERRY_REFSEQ} = "b.$bl2_tables{CHERRYMAP}{REFSEQ_ID_FIELD}=c.$bl2_tables{REFSEQ}{REFSEQ_ID_FIELD}";

  $bl2_tables{FILTER}{BL2SEQ} = " bl2_score is null"; # use to put joined condition
 
  $bl2_tables{ORDER_BY}{BL2SEQ} = "$bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}, $bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD}";

  my $bl2seq_type = 0;
  if($BL2SEQ_W_BRIEF_BLAST==1) {
      $bl2seq_type = 1;
  }
  MY_SeqAlign -> general_bl2seq($bl2_tables{TRACE_SEQ}{DBH}, \%bl2_tables, $bl2seq_type);

# catch those didn't pass quality/vector trim by using raw trace:

  if($qc_rtInfo->{USE_RAW_TRACE} != 1) {
    $bl2_tables{TRACE_SEQ}{TRACE_SEQ_FIELD} = 'TRACE_SEQ';
    $bl2_tables{TRACE_SEQ}{SEQ_OFFSET_FIELD} = 1;
# put the filter in the single table:
    $bl2_tables{TRACE_SEQ}{FILTER} = "$bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD} in ('" . join("','",@trace_plate_ids) . "') ";
    $bl2_tables{TRACE_SEQ}{FILTER} .= " and $bl2_tables{TRACE_SEQ}{BL2_SCORE_FIELD} is null"; # only fetch those didn't align first time 
#    $bl2_tables{TRACE_SEQ}{FILTER} .= " and $bl2_tables{TRACE_SEQ}{TRACE_ID_FIELD} not in (1793)"; # debug
#  $bl2_tables{TRACE_SEQ}{FILTER} .= " and $bl2_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}=11687 and $bl2_tables{TRACE_SEQ}{TRACE_POS_FIELD}='A06'"; # debug

#rerun bl2seq using raw trace seq    
     MY_SeqAlign -> general_bl2seq($bl2_tables{TRACE_SEQ}{DBH}, \%bl2_tables, $bl2seq_type);
  }
} # end bl2seq

# calc phred score:
#  my $CALC_PHRED = 1;
  my $CALC_PHRED = defined($qc_rtInfo->{CALC_PHRED})? $qc_rtInfo->{CALC_PHRED}:0;

if ($CALC_PHRED == 1) {
  my %phred_tables = (
    TRACE_SEQ => $t->{TRACE_SEQ},
  );

  $phred_tables{TRACE_SEQ}{TRACE_ID_FIELD} = 'TRACE_ID';
  $phred_tables{TRACE_SEQ}{TRACE_PLATE_FIELD} = 'TRACE_PLATE';
  $phred_tables{TRACE_SEQ}{TRACE_POS_FIELD} = 'TRACE_POS';
  $phred_tables{TRACE_SEQ}{TRACE_QUAL_FIELD} = 'TRACE_QUAL';
  $phred_tables{TRACE_SEQ}{AVG_PHRED_FIELD} = 'AVG_PHRED';
  $phred_tables{TRACE_SEQ}{INSERT_START_FIELD} = 'INSERT_START';
  $phred_tables{TRACE_SEQ}{INSERT_END_FIELD} = 'INSERT_END';

  $phred_tables{FILTER}{PHRED} = "$phred_tables{TRACE_SEQ}{TRACE_PLATE_FIELD} in ('" . join("','",@trace_plate_ids) . "') and $phred_tables{TRACE_SEQ}{INSERT_START_FIELD} is not null"; # only fetch those passing quality trim

#  $phred_tables{FILTER}{PHRED} .= " and $phred_tables{TRACE_SEQ}{TRACE_ID_FIELD}=1153342"; # debug
#  $phred_tables{FILTER}{PHRED} .= " and $phred_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}=11687 and $phred_tables{TRACE_SEQ}{TRACE_POS_FIELD}='A06'"; # debug

  $phred_tables{ORDER_BY}{PHRED} = "$phred_tables{TRACE_SEQ}{TRACE_PLATE_FIELD}, $phred_tables{TRACE_SEQ}{TRACE_POS_FIELD}";

  MY_SeqAlign -> calc_phred_score(\%phred_tables);

} # end of calc phred score

# run t-coffee:
#  my $TCOFFEE = 0;
  my $TCOFFEE = defined($qc_rtInfo->{TCOFFEE})? $qc_rtInfo->{TCOFFEE}:0;

if ($TCOFFEE==1) {
  my $seq_table = "("
    . "select distinct concat(y.dstPlate, y.dstPos,'_', y.CUSTOM_ID, '_ORF', y.ORF_ID) as GROUP_ID, concat(y.ORF_ID, '_REF') as SEQ_NAME, x.$t->{REFSEQ}{REFSEQ_FIELD} as SEQ from $t->{REFSEQ}{DBNAME}.$t->{REFSEQ}{TABLENAME} x inner join $t->{CHERRYMAP}{TABLENAME} y on x.$t->{REFSEQ}{REFSEQ_ID_FIELD}=y.ORF_ID where y.DstPlate in ('" . join("','",@trace_plate_ids) . "')" # ref seq
    . " union "
    . " select distinct concat(y.dstPlate, y.dstPos,'_', y.CUSTOM_ID, '_ORF', y.ORF_ID) as GROUP_ID, concat(y.ORF_ID, '_t', x.TRACE_ID) as SEQ_NAME, x.TRACE_SEQ as SEQ from $t->{TRACE_SEQ}{TABLENAME} x inner join $t->{CHERRYMAP}{TABLENAME} y on x.TRACE_PLATE=y.DstPlate and x.TRACE_POS=y.DstPos where y.DstPlate in ('" . join("','",@trace_plate_ids) . "')" # raw trace
    . " union "
    . " select distinct concat(y.dstPlate, y.dstPos,'_', y.CUSTOM_ID, '_ORF', y.ORF_ID) as GROUP_ID, concat(y.ORF_ID, '_t', x.TRACE_ID, '_TRIM') as SEQ_NAME, substr(x.TRACE_SEQ, x.INSERT_START, x.INSERT_END-x.INSERT_START+1) as SEQ from $t->{TRACE_SEQ}{TABLENAME} x inner join $t->{CHERRYMAP}{TABLENAME} y on x.TRACE_PLATE=y.DstPlate and x.TRACE_POS=y.DstPos where y.DstPlate in ('" . join("','",@trace_plate_ids) . "') and x.INSERT_START is not null" # trim seq
    .")";
  my %tcoffee_tables = (
    SEQ => {
      DBH => $t->{TRACE_SEQ}{DBH},
      TABLENAME => $seq_table,
    },
  );
  my %tcoffee_rtinfo = (
    RESULT_DIR => $result_dir,
  );

  MY_SeqAlign -> general_tcoffee(\%tcoffee_tables, \%tcoffee_rtinfo);

} # end of tcoffee


# analyze result:

#  my $ANALYZE_RESULT = 1;
  my $ANALYZE_RESULT = defined($qc_rtInfo->{ANALYZE_RESULT})? $qc_rtInfo->{ANALYZE_RESULT}:0;

if($ANALYZE_RESULT == 1) {
# good seq criteria: IDENTITY=100, ALIGN_LEN>60
  my $file_prefix = $source; #use $source as the file prefix
  my %cutoff = (
     AVG_QUAL => 30,
     BL2_SCORE => defined($qc_rtInfo->{CUTOFF_BL2_SCORE})? $qc_rtInfo->{CUTOFF_BL2_SCORE}:60,
     IDENTITY => 100,
     ALIGN_LEN => defined($qc_rtInfo->{CUTOFF_ALIGN_LEN})? $qc_rtInfo->{CUTOFF_ALIGN_LEN}:60,
  );
  my %summary = ();
  my %overall = ();

# 0. get overall result:
  $sql = "select count(distinct dstPlate, dstPos) from $t->{CHERRYMAP}{TABLENAME} a where dstPlate in (" . join (',', @trace_plate_ids) .")";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{TOTAL_PICK} = $query->fetchrow();
  $query -> finish;

  $sql = "select count(distinct TRACE_ID) from $t->{TRACE_SEQ}{TABLENAME} a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and TRACE_PLATE in (" . join(',', @trace_plate_ids) .")";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{PASS_QUAL} = $query -> fetchrow();
  $query -> finish;

# 2b.  sequences pass quality trim: only by each POS, no consider TRACE_DIR. 
  $sql = "select count(distinct TRACE_PLATE, TRACE_POS) from $t->{TRACE_SEQ}{TABLENAME} a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and TRACE_PLATE in (" . join(',', @trace_plate_ids) .")";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{PASS_QUAL_WELL} = $query -> fetchrow();
  $query -> finish;

# 3. sequences pass bl2 score criteria for those pass quality criteria:
  $sql = "select count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE}";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{PASS_BL2} = $query -> fetchrow();
  $query -> finish;

# 3b. sequences pass bl2 score criteria for those pass quality criteria (combine FOR/REV):
  $sql = "select count(distinct TRACE_PLATE, TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE}";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{PASS_BL2_WELL} = $query -> fetchrow();
  $query -> finish;

# 4. perfect alignment (no mutation) count, start from ATG:
  $sql = "select count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY>=$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3) = 'ATG' and ALIGN_LEN>=$cutoff{ALIGN_LEN}";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{WILD_TYPE} = $query -> fetchrow();
  $query -> finish;

# 4b. perfect alignment (no mutation) count, start from ATG (combined FOR/REV):
  $sql = "select count(distinct TRACE_PLATE, TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY>=$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3) = 'ATG' and ALIGN_LEN>=$cutoff{ALIGN_LEN}";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{WILD_TYPE_WELL} = $query -> fetchrow();
  $query -> finish;

# 5a.grab all those have ATG, not perfect match:
  $sql = "select count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)='ATG'";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{FAIL_IDENTITY} = $query -> fetchrow();
  $query -> finish;

# 5a'.grab all those have ATG, not perfect match (combine FOR/REV):
  $sql = "select count(distinct TRACE_PLATE, TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)='ATG'";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{FAIL_IDENTITY_WELL} = $query -> fetchrow();
  $query -> finish;

# 5b.grab all those not perfect match, no start from ATG:
  $sql = "select count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY>=$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG'";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{NO_ATG} = $query -> fetchrow();
  $query -> finish;

# 5b'.grab all those not perfect match, no start from ATG (combine FOR/REV):
  $sql = "select count(distinct TRACE_PLATE, TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY>=$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG'";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{NO_ATG_WELL} = $query -> fetchrow();
  $query -> finish;

# 5c.grab all those not perfect match, no start from ATG:
  $sql = "select count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG'";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{FAIL_IDENTITY_NOATG} = $query -> fetchrow();
  $query -> finish;

# 5c'.grab all those not perfect match, no start from ATG (combine FOR/REV):
  $sql = "select count(distinct TRACE_PLATE, TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG'";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $overall{FAIL_IDENTITY_NOATG_WELL} = $query -> fetchrow();
  $query -> finish;

# print out overall:
if($qc_rtInfo->{BOTH_DIR} == 1) {
  open OUT, ">$result_dir$file_prefix" . "_overall.csv";
  print OUT "Below is the overall QC result for plates $trace_plate_ids[0] to $trace_plate_ids[scalar(@trace_plate_ids)]:\n\n";
  print OUT join "\t", "TOTAL_PICK", "PASS_QUAL(combineF/R)", "PASS_BL2(Cmb_F/R)", "WILD_TYPE(Cmb_F/R)", "IDENTITY<100(Cmb_F/R)", "NO_ATG(Cmb_F/R)", "IDENTITY<100_and_NOATG(Cmb_F/R)", "%PASS_QUAL(Cmb_F/R)", "%PASS_BL2 vs TOTAL_PICK(Cmb_F/R)", "%PASS_BL2 vs PASS_QUAL(Cmb_F/R)", "%WILD_TYPE(Cmb_F/R)", "%IDENTITY<100(Cmb_F/R)", "%NO_ATG(Cmb_F/R)", "%IDENTITY<100_and_NOATG(Cmb_F/R)";
  print OUT "\n"; 
  print OUT join "\t", 
      defined($overall{TOTAL_PICK})?$overall{TOTAL_PICK}:0,
      sprintf("%d(%d)", defined($overall{PASS_QUAL})? $overall{PASS_QUAL}: 0, defined($overall{PASS_QUAL_WELL})? $overall{PASS_QUAL_WELL}: 0),
      sprintf("%d(%d)", defined($overall{PASS_BL2})? $overall{PASS_BL2}: 0, defined($overall{PASS_BL2_WELL})? $overall{PASS_BL2_WELL}: 0),
      sprintf("%d(%d)", defined($overall{WILD_TYPE})? $overall{WILD_TYPE}: 0, defined($overall{WILD_TYPE_WELL})? $overall{WILD_TYPE_WELL}: 0),
      sprintf("%d(%d)", defined($overall{FAIL_IDENTITY})? $overall{FAIL_IDENTITY}: 0, defined($overall{FAIL_IDENTITY_WELL})? $overall{FAIL_IDENTITY_WELL}: 0),
      sprintf("%d(%d)", defined($overall{NO_ATG})? $overall{NO_ATG}: 0, defined($overall{NO_ATG_WELL})? $overall{NO_ATG_WELL}: 0),
      sprintf("%d(%d)", defined($overall{FAIL_IDENTITY_NOATG})? $overall{FAIL_IDENTITY_NOATG}: 0, defined($overall{FAIL_IDENTITY_NOATG_WELL})? $overall{FAIL_IDENTITY_NOATG_WELL}: 0),
      sprintf("%4.2f(%4.2f)", defined($overall{PASS_QUAL}) && defined($overall{TOTAL_PICK}) ? $overall{PASS_QUAL}*100/$overall{TOTAL_PICK} : 0, defined($overall{PASS_QUAL_WELL}) && defined($overall{TOTAL_PICK}) ? $overall{PASS_QUAL_WELL}*100/$overall{TOTAL_PICK} : 0),
      sprintf("%4.2f(%4.2f)", defined($overall{PASS_BL2}) && defined($overall{TOTAL_PICK}) ? $overall{PASS_BL2}*100/$overall{TOTAL_PICK} : 0, defined($overall{PASS_BL2_WELL}) && defined($overall{TOTAL_PICK}) ? $overall{PASS_BL2_WELL}*100/$overall{TOTAL_PICK} : 0),
      sprintf("%4.2f(%4.2f)", defined($overall{PASS_BL2}) && defined($overall{PASS_QUAL}) ? $overall{PASS_BL2}*100/$overall{PASS_QUAL} : 0, defined($overall{PASS_BL2_WELL}) && defined($overall{PASS_QUAL_WELL}) ? $overall{PASS_BL2_WELL}*100/$overall{PASS_QUAL_WELL} : 0),
      sprintf("%4.2f(%4.2f)", defined($overall{WILD_TYPE}) && defined($overall{PASS_QUAL}) ? $overall{WILD_TYPE}*100/$overall{PASS_QUAL} : 0, defined($overall{WILD_TYPE_WELL}) && defined($overall{PASS_QUAL_WELL}) ? $overall{WILD_TYPE_WELL}*100/$overall{PASS_QUAL_WELL} : 0),
      sprintf("%4.2f(%4.2f)", defined($overall{FAIL_IDENTITY}) && defined($overall{PASS_QUAL}) ? $overall{FAIL_IDENTITY}*100/$overall{PASS_QUAL} : 0, defined($overall{FAIL_IDENTITY_WELL}) && defined($overall{PASS_QUAL_WELL}) ? $overall{FAIL_IDENTITY_WELL}*100/$overall{PASS_QUAL_WELL} : 0),
      sprintf("%4.2f(%4.2f)", defined($overall{NO_ATG}) && defined($overall{PASS_QUAL}) ? $overall{NO_ATG}*100/$overall{PASS_QUAL} : 0, defined($overall{NO_ATG_WELL}) && defined($overall{PASS_QUAL_WELL}) ? $overall{NO_ATG_WELL}*100/$overall{PASS_QUAL_WELL} : 0),
      sprintf("%4.2f(%4.2f)", defined($overall{FAIL_IDENTITY_NOATG}) && defined($overall{PASS_QUAL}) ? $overall{FAIL_IDENTITY_NOATG}*100/$overall{PASS_QUAL} : 0, defined($overall{FAIL_IDENTITY_NOATG_WELL}) && defined($overall{PASS_QUAL_WELL}) ? $overall{FAIL_IDENTITY_NOATG_WELL}*100/$overall{PASS_QUAL_WELL} : 0);
  print OUT "\n";
  close OUT;
} # end of both dir

else { # single dir. trace
  open OUT, ">$result_dir$file_prefix" . "_overall.csv";
  print OUT join "\t", "TOTAL_PICK", "PASS_QUAL", "PASS_BL2", "WILD_TYPE", "IDENTITY<100", "NO_ATG", "IDENTITY<100_and_NOATG", "%PASS_QUAL", "%PASS_BL2 vs TOTAL_PICK", "%PASS_BL2 vs PASS_QUAL", "%WILD_TYPE", "%IDENTITY<100", "%NO_ATG", "%IDENTITY<100_and_NOATG";
  print OUT "\n"; 
  print OUT join "\t", 
      defined($overall{TOTAL_PICK})?$overall{TOTAL_PICK}:0,
      sprintf("%d", defined($overall{PASS_QUAL})? $overall{PASS_QUAL}: 0),
      sprintf("%d", defined($overall{PASS_BL2})? $overall{PASS_BL2}: 0),
      sprintf("%d", defined($overall{WILD_TYPE})? $overall{WILD_TYPE}: 0),
      sprintf("%d", defined($overall{FAIL_IDENTITY})? $overall{FAIL_IDENTITY}: 0),
      sprintf("%d", defined($overall{NO_ATG})? $overall{NO_ATG}: 0),
      sprintf("%d", defined($overall{FAIL_IDENTITY_NOATG})? $overall{FAIL_IDENTITY_NOATG}: 0),
      sprintf("%4.2f", defined($overall{PASS_QUAL}) && defined($overall{TOTAL_PICK}) ? $overall{PASS_QUAL}*100/$overall{TOTAL_PICK} : 0),
      sprintf("%4.2f", defined($overall{PASS_BL2}) && defined($overall{TOTAL_PICK}) ? $overall{PASS_BL2}*100/$overall{TOTAL_PICK} : 0),
      sprintf("%4.2f", defined($overall{PASS_BL2}) && defined($overall{PASS_QUAL}) ? $overall{PASS_BL2}*100/$overall{PASS_QUAL} : 0),
      sprintf("%4.2f", defined($overall{WILD_TYPE}) && defined($overall{PASS_QUAL}) ? $overall{WILD_TYPE}*100/$overall{PASS_QUAL} : 0),
      sprintf("%4.2f", defined($overall{FAIL_IDENTITY}) && defined($overall{PASS_QUAL}) ? $overall{FAIL_IDENTITY}*100/$overall{PASS_QUAL} : 0),
      sprintf("%4.2f", defined($overall{NO_ATG}) && defined($overall{PASS_QUAL}) ? $overall{NO_ATG}*100/$overall{PASS_QUAL} : 0),
      sprintf("%4.2f", defined($overall{FAIL_IDENTITY_NOATG}) && defined($overall{PASS_QUAL}) ? $overall{FAIL_IDENTITY_NOATG}*100/$overall{PASS_QUAL} : 0);
  print OUT "\n";
  close OUT;
}
  
# get stat. data by plates:
# 1. all sequence picked:
  $sql = "select dstPlate, count(distinct dstPos) from $t->{CHERRYMAP}{TABLENAME} a where dstPlate in (" . join (',', @trace_plate_ids) .") group by dstPlate order by dstPlate";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_total_picked.tab";
  print OUT join "\t", "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{TOTAL_PICK} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 2. sequences pass quality trim:
  $sql = "select TRACE_PLATE, count(distinct TRACE_ID) from $t->{TRACE_SEQ}{TABLENAME} a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and TRACE_PLATE in (" . join(',', @trace_plate_ids) .") group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_qual_count.tab";
  print OUT join "\t", "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{PASS_QUAL} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 2b.  sequences pass quality trim: only by each POS, no consider TRACE_DIR. 
  $sql = "select TRACE_PLATE, count(distinct TRACE_POS) from $t->{TRACE_SEQ}{TABLENAME} a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and TRACE_PLATE in (" . join(',', @trace_plate_ids) .") group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_qual_count_well.tab";
  print OUT join "\t", "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{PASS_QUAL_WELL} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 3. sequences pass bl2 score criteria for those pass quality criteria:
  $sql = "select TRACE_PLATE, count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_score60_count.tab";
  print OUT "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{PASS_BL2} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 3b. sequences pass bl2 score criteria for those pass quality criteria (combine FOR/REV):
  $sql = "select TRACE_PLATE, count(distinct TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_score60_count_well.tab";
  print OUT join "\t", "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{PASS_BL2_WELL} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 4. perfect alignment (no mutation) count, start from ATG:
  $sql = "select TRACE_PLATE, count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY>=$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3) = 'ATG' and ALIGN_LEN>=$cutoff{ALIGN_LEN} group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_wt_count.tab";
  print OUT join "\t", "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{WILD_TYPE} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 4b. perfect alignment (no mutation) count, start from ATG (combined FOR/REV):
  $sql = "select TRACE_PLATE, count(distinct TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY>=$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3) = 'ATG' and ALIGN_LEN>=$cutoff{ALIGN_LEN} group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_wt_count_well.tab";
  print OUT join "\t", "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{WILD_TYPE_WELL} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 5a.grab all those have ATG, not perfect match:
  $sql = "select TRACE_PLATE, count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)='ATG' group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_idLess100.tab";
  print OUT "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{FAIL_IDENTITY} = $COUNT;    
  }
  $query -> finish;
  close OUT;

# 5a'.grab all those have ATG, not perfect match (combine FOR/REV):
  $sql = "select TRACE_PLATE, count(distinct TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)='ATG' group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_idLess100_well.tab";
  print OUT join "\t", "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{FAIL_IDENTITY_WELL} = $COUNT;    
  }
  $query -> finish;
  close OUT;

# 5b.grab all those not perfect match, no start from ATG:
  $sql = "select TRACE_PLATE, count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY>=$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG' group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_noATG.tab";
  print OUT join "\t", "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{NO_ATG} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 5b'.grab all those not perfect match, no start from ATG (combine FOR/REV):
  $sql = "select TRACE_PLATE, count(distinct TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY>=$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG' group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_noATG_well.tab";
  print OUT "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{NO_ATG_WELL} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 5c.grab all those not perfect match, no start from ATG:
  $sql = "select TRACE_PLATE, count(distinct TRACE_ID) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG' group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_noATG_less100.tab";
  print OUT "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{FAIL_IDENTITY_NOATG} = $COUNT;
  }
  $query -> finish;
  close OUT;

# 5c'.grab all those not perfect match, no start from ATG (combine FOR/REV):
  $sql = "select TRACE_PLATE, count(distinct TRACE_POS) from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG' group by TRACE_PLATE  order by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_noATG_less100.tab";
  print OUT "TRACE_PLATE, COUNT";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $COUNT) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $COUNT;
    print OUT "\n";
    $summary{$TRACE_PLATE}{FAIL_IDENTITY_NOATG_WELL} = $COUNT;
  }
  $query -> finish;
  close OUT;

# now start to get details:

# 5d.grab all those have ATG, but not perfect match in detail:

# put these fields needs to be output in $qc_rtInfo->{DISP_CUSTOM_FIELDS}, such as dnaSeqCode, etc

  if ($qc_rtInfo->{BOTH_DIR} == 1) {
    push @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}, "TRACE_DIR";
  }

  $sql = "select distinct " . (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join(",", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) . "," : "") . " a.AVG_PHRED,b.CUSTOM_ID as GENE_ID, b.ORF_ID, a.TRACE_PLATE, a.TRACE_POS, a.TRACE_SEQ, substr(a.TRACE_SEQ, a.INSERT_START, a.INSERT_END-a.INSERT_START+1) as TRACE_TRIM, a.BL2_SCORE, a.IDENTITY, a.TRACE_AVG_QUAL, length($t->{REFSEQ}{REFSEQ_FIELD}) as ORF_SIZE, c.$t->{REFSEQ}{REFSEQ_FIELD} from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS inner join $t->{REFSEQ}{TABLENAME} c on c.$t->{REFSEQ}{REFSEQ_ID_FIELD} = b.ORF_ID where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)='ATG' order by TRACE_PLATE, TRACE_POS";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_idLess100_detail.tab";
  print OUT join "\t", (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join(",", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) : "") , "AVG_PHRED", "GENE_ID", "ORF_ID", "TRACE_PLATE", "TRACE_POS", "TRACE_SEQ", "TRACE_TRIM", "BL2_SCORE", "IDENTITY", "TRACE_AVG_QUAL", "ORF_SIZE", $t->{REFSEQ}{REFSEQ_FIELD};
  print  OUT "\n";
  while (my (@cols) = $query -> fetchrow()) {
    print OUT join "\t", @cols;
    print OUT "\n";
  }

  $query -> finish;
  close OUT;

# 5e.grab all those perfect match, but no ATG in detail:
  $sql = "select distinct " . (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join(",", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) . "," : "") . " a.AVG_PHRED, b.CUSTOM_ID as GENE_ID, b.ORF_ID, a.TRACE_PLATE, a.TRACE_POS, a.TRACE_SEQ, substr(a.TRACE_SEQ, a.INSERT_START, a.INSERT_END-a.INSERT_START+1) as TRACE_TRIM, a.BL2_SCORE, a.IDENTITY, a.TRACE_AVG_QUAL, length($t->{REFSEQ}{REFSEQ_FIELD}) as ORF_SIZE, c.$t->{REFSEQ}{REFSEQ_FIELD} from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS inner join $t->{REFSEQ}{TABLENAME} c on c.$t->{REFSEQ}{REFSEQ_ID_FIELD} = b.ORF_ID where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY>=$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG' order by TRACE_PLATE, TRACE_POS";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_noATG_detail.tab";
  print OUT join "\t", (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join(",", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) : "") , "AVG_PHRED", "GENE_ID", "ORF_ID", "TRACE_PLATE", "TRACE_POS", "TRACE_SEQ", "TRACE_TRIM", "BL2_SCORE", "IDENTITY", "TRACE_AVG_QUAL", "ORF_SIZE", $t->{REFSEQ}{REFSEQ_FIELD};
  print  OUT "\n";
  while (my (@cols) = $query -> fetchrow()) {
    print OUT join "\t", @cols;
    print OUT "\n";
  }

  $query -> finish;
  close OUT;

# 5f.grab all those not perfect match, no ATG in detail:
  $sql = "select distinct " . (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join(",", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) . "," : "") . " a.AVG_PHRED, b.CUSTOM_ID as GENE_ID, b.ORF_ID, a.TRACE_PLATE, a.TRACE_POS, a.TRACE_SEQ, substr(a.TRACE_SEQ, a.INSERT_START, a.INSERT_END-a.INSERT_START+1) as TRACE_TRIM, a.BL2_SCORE, a.IDENTITY, a.TRACE_AVG_QUAL, length(c.$t->{REFSEQ}{REFSEQ_FIELD}) as ORF_SIZE, c.$t->{REFSEQ}{REFSEQ_FIELD} from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS inner join $t->{REFSEQ}{TABLENAME} c on c.$t->{REFSEQ}{REFSEQ_ID_FIELD} = b.ORF_ID where TRACE_AVG_QUAL>$cutoff{AVG_QUAL} and BL2_SCORE>$cutoff{BL2_SCORE} and IDENTITY<$cutoff{IDENTITY} and substr(TRACE_SEQ, INSERT_START, 3)!='ATG' order by TRACE_PLATE, TRACE_POS";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_noATG_less100_detail.tab";
  print OUT join "\t", (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join(",", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) : "") , "AVG_PHRED", "GENE_ID", "ORF_ID", "TRACE_PLATE", "TRACE_POS", "TRACE_SEQ", "TRACE_TRIM", "BL2_SCORE", "IDENTITY", "TRACE_AVG_QUAL", "ORF_SIZE", $t->{REFSEQ}{REFSEQ_FIELD};
  print  OUT "\n";
  while (my (@cols) = $query -> fetchrow()) {
    print OUT join "\t", @cols;
    print OUT "\n";
  }
  $query -> finish;
  close OUT;

# 6. grab all those not aligned at all in detail:
  $sql = "select distinct " . (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join(",", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) . "," : "") . " a.AVG_PHRED, b.CUSTOM_ID as GENE_ID, b.ORF_ID, a.TRACE_PLATE, a.TRACE_POS, a.TRACE_SEQ, substr(a.TRACE_SEQ, a.INSERT_START, a.INSERT_END-a.INSERT_START+1) as TRACE_TRIM, a.BL2_SCORE, a.IDENTITY, a.TRACE_AVG_QUAL, length(c.$t->{REFSEQ}{REFSEQ_FIELD}) as ORF_SIZE, ";
  if($BL2SEQ_W_BRIEF_BLAST == 1) {
    $sql .= "a.qSeqIndIdentical, a.sSeqIndIdentical, a.qSeqLenIdentical, a.sSeqLenIdentical, a.qSeqGap, a.sSeqGap, ";
  }
  $sql .= "c.$t->{REFSEQ}{REFSEQ_FIELD} from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS inner join $t->{REFSEQ}{TABLENAME} c on c.$t->{REFSEQ}{REFSEQ_ID_FIELD} = b.ORF_ID where TRACE_AVG_QUAL>$cutoff{AVG_QUAL}  and IDENTITY is null  order by TRACE_PLATE, TRACE_POS";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_no_hit_detail.tab";
  print OUT join "\t", (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join("\t", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) : "") , "AVG_PHRED", "GENE_ID", "ORF_ID", "TRACE_PLATE", "TRACE_POS", "TRACE_SEQ", "TRACE_TRIM", "BL2_SCORE", "IDENTITY", "TRACE_AVG_QUAL", "ORF_SIZE"; 
  if($BL2SEQ_W_BRIEF_BLAST == 1) {
    print OUT "\t";
    print OUT join "\t", "qSeqIndIdentical", "sSeqIndIdentical", "qSeqLenIdentical", "sSeqLenIdentical", "qSeqGap", "sSeqGap";
  }
  print OUT "\t";
  print OUT join "\t", $t->{REFSEQ}{REFSEQ_FIELD};
  print  OUT "\n";
  while (my (@cols) = $query -> fetchrow()) {
    print OUT join "\t", @cols;
    print OUT "\n";
  }
  $query -> finish;
  close OUT;

#8. grab all sequences:
  $sql = "select distinct " . (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join(",", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) . "," : "") . " a.AVG_PHRED, b.CUSTOM_ID as GENE_ID, b.ORF_ID, a.TRACE_PLATE, a.TRACE_POS, a.TRACE_SEQ, substr(a.TRACE_SEQ, a.INSERT_START, a.INSERT_END-a.INSERT_START+1) as TRACE_TRIM, a.BL2_SCORE, a.IDENTITY, a.TRACE_AVG_QUAL, ";
  if($BL2SEQ_W_BRIEF_BLAST == 1) {
    $sql .= "a.qSeqIndIdentical, a.sSeqIndIdentical, a.qSeqLenIdentical, a.sSeqLenIdentical, a.qSeqGap, a.sSeqGap, ";
  }
  $sql .= "length(c.$t->{REFSEQ}{REFSEQ_FIELD}) as ORF_SIZE, c.$t->{REFSEQ}{REFSEQ_FIELD} from (select * from $t->{TRACE_SEQ}{TABLENAME} where TRACE_PLATE in (" . join(',', @trace_plate_ids) .")) a inner join $t->{CHERRYMAP}{TABLENAME} b on b.dstPlate=a.TRACE_PLATE and b.dstPos=a.TRACE_POS inner join $t->{REFSEQ}{TABLENAME} c on c.$t->{REFSEQ}{REFSEQ_ID_FIELD} = b.ORF_ID order by TRACE_PLATE, TRACE_POS";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_all_detail.tab";
  print OUT join "\t", (defined($qc_rtInfo->{DISP_CUSTOM_FIELDS})? join("\t", @{$qc_rtInfo->{DISP_CUSTOM_FIELDS}}) : "") , "AVG_PHRED", "GENE_ID", "ORF_ID", "TRACE_PLATE", "TRACE_POS", "TRACE_SEQ", "TRACE_TRIM", "BL2_SCORE", "IDENTITY", "TRACE_AVG_QUAL", "ORF_SIZE";
  if($BL2SEQ_W_BRIEF_BLAST == 1) {
    print OUT "\t";
    print OUT join "\t", "qSeqIndIdentical", "sSeqIndIdentical", "qSeqLenIdentical", "sSeqLenIdentical", "qSeqGap", "sSeqGap";
  }
  print OUT "\t";
  print OUT join "\t", $t->{REFSEQ}{REFSEQ_FIELD};
  print  OUT "\n";
  while (my (@cols) = $query -> fetchrow()) {
    print OUT join "\t", @cols;
    print OUT "\n";
  }
  $query -> finish;
  close OUT;

#9. calc phred score:
  $sql = "select distinct TRACE_PLATE, sum(AVG_PHRED*(INSERT_END-INSERT_START+1))/sum(INSERT_END-INSERT_START+1) as AVG_PHRED from $t->{TRACE_SEQ}{TABLENAME} a where TRACE_PLATE in (" . join(",", @trace_plate_ids) .") group by TRACE_PLATE";
  $query = $t->{TRACE_SEQ}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_avg_phred.tab";
  print OUT join "\t", "TRACE_PLATE_ID", "AVG_PHRED";
  print  OUT "\n";
  while (my ($TRACE_PLATE, $AVG_PHRED) = $query -> fetchrow()) {
    print OUT join "\t", $TRACE_PLATE, $AVG_PHRED;
    print OUT "\n";
    $summary{$TRACE_PLATE}{AVG_PHRED} = $AVG_PHRED;
  }
  $query -> finish;
  close OUT;

# print out summary:
if($qc_rtInfo->{BOTH_DIR} == 1) {
  open OUT, ">$result_dir$file_prefix" . "_summary.tab";
  print OUT join "\t", "PLATE_ID", "AVG_PHRED", "TOTAL_PICK", "PASS_QUAL(combineF/R)", "PASS_BL2(Cmb_F/R)", "WILD_TYPE(Cmb_F/R)", "IDENTITY<100(Cmb_F/R)", "NO_ATG(Cmb_F/R)", "IDENTITY<100_and_NOATG(Cmb_F/R)", "%PASS_QUAL(Cmb_F/R)", "%PASS_BL2 vs TOTAL_PICK(Cmb_F/R)", "%PASS_BL2 vs PASS_QUAL(Cmb_F/R)", "%WILD_TYPE(Cmb_F/R)", "%IDENTITY<100(Cmb_F/R)", "%NO_ATG(Cmb_F/R)", "%IDENTITY<100_and_NOATG(Cmb_F/R)";
  print OUT "\n"; 
  foreach my $trace_plate_id (sort keys %summary) {
    print OUT join "\t", 
      $trace_plate_id, 
      defined($summary{$trace_plate_id}{AVG_PHRED})?$summary{$trace_plate_id}{AVG_PHRED}:"n/a",
      defined($summary{$trace_plate_id}{TOTAL_PICK})?$summary{$trace_plate_id}{TOTAL_PICK}:0,
      sprintf("%d(%d)", defined($summary{$trace_plate_id}{PASS_QUAL})? $summary{$trace_plate_id}{PASS_QUAL}: 0, defined($summary{$trace_plate_id}{PASS_QUAL_WELL})? $summary{$trace_plate_id}{PASS_QUAL_WELL}: 0),
      sprintf("%d(%d)", defined($summary{$trace_plate_id}{PASS_BL2})? $summary{$trace_plate_id}{PASS_BL2}: 0, defined($summary{$trace_plate_id}{PASS_BL2_WELL})? $summary{$trace_plate_id}{PASS_BL2_WELL}: 0),
      sprintf("%d(%d)", defined($summary{$trace_plate_id}{WILD_TYPE})? $summary{$trace_plate_id}{WILD_TYPE}: 0, defined($summary{$trace_plate_id}{WILD_TYPE_WELL})? $summary{$trace_plate_id}{WILD_TYPE_WELL}: 0),
      sprintf("%d(%d)", defined($summary{$trace_plate_id}{FAIL_IDENTITY})? $summary{$trace_plate_id}{FAIL_IDENTITY}: 0, defined($summary{$trace_plate_id}{FAIL_IDENTITY_WELL})? $summary{$trace_plate_id}{FAIL_IDENTITY_WELL}: 0),
      sprintf("%d(%d)", defined($summary{$trace_plate_id}{NO_ATG})? $summary{$trace_plate_id}{NO_ATG}: 0, defined($summary{$trace_plate_id}{NO_ATG_WELL})? $summary{$trace_plate_id}{NO_ATG_WELL}: 0),
      sprintf("%d(%d)", defined($summary{$trace_plate_id}{FAIL_IDENTITY_NOATG})? $summary{$trace_plate_id}{FAIL_IDENTITY_NOATG}: 0, defined($summary{$trace_plate_id}{FAIL_IDENTITY_NOATG_WELL})? $summary{$trace_plate_id}{FAIL_IDENTITY_NOATG_WELL}: 0),
      sprintf("%4.2f(%4.2f)", defined($summary{$trace_plate_id}{PASS_QUAL}) && defined($summary{$trace_plate_id}{TOTAL_PICK}) ? $summary{$trace_plate_id}{PASS_QUAL}*100/$summary{$trace_plate_id}{TOTAL_PICK} : 0, defined($summary{$trace_plate_id}{PASS_QUAL_WELL}) && defined($summary{$trace_plate_id}{TOTAL_PICK}) ? $summary{$trace_plate_id}{PASS_QUAL_WELL}*100/$summary{$trace_plate_id}{TOTAL_PICK} : 0),
      sprintf("%4.2f(%4.2f)", defined($summary{$trace_plate_id}{PASS_BL2}) && defined($summary{$trace_plate_id}{TOTAL_PICK}) ? $summary{$trace_plate_id}{PASS_BL2}*100/$summary{$trace_plate_id}{TOTAL_PICK} : 0, defined($summary{$trace_plate_id}{PASS_BL2_WELL}) && defined($summary{$trace_plate_id}{TOTAL_PICK}) ? $summary{$trace_plate_id}{PASS_BL2_WELL}*100/$summary{$trace_plate_id}{TOTAL_PICK} : 0),
      sprintf("%4.2f(%4.2f)", defined($summary{$trace_plate_id}{PASS_BL2}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{PASS_BL2}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0, defined($summary{$trace_plate_id}{PASS_BL2_WELL}) && defined($summary{$trace_plate_id}{PASS_QUAL_WELL}) ? $summary{$trace_plate_id}{PASS_BL2_WELL}*100/$summary{$trace_plate_id}{PASS_QUAL_WELL} : 0),
      sprintf("%4.2f(%4.2f)", defined($summary{$trace_plate_id}{WILD_TYPE}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{WILD_TYPE}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0, defined($summary{$trace_plate_id}{WILD_TYPE_WELL}) && defined($summary{$trace_plate_id}{PASS_QUAL_WELL}) ? $summary{$trace_plate_id}{WILD_TYPE_WELL}*100/$summary{$trace_plate_id}{PASS_QUAL_WELL} : 0),
      sprintf("%4.2f(%4.2f)", defined($summary{$trace_plate_id}{FAIL_IDENTITY}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{FAIL_IDENTITY}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0, defined($summary{$trace_plate_id}{FAIL_IDENTITY_WELL}) && defined($summary{$trace_plate_id}{PASS_QUAL_WELL}) ? $summary{$trace_plate_id}{FAIL_IDENTITY_WELL}*100/$summary{$trace_plate_id}{PASS_QUAL_WELL} : 0),
      sprintf("%4.2f(%4.2f)", defined($summary{$trace_plate_id}{NO_ATG}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{NO_ATG}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0, defined($summary{$trace_plate_id}{NO_ATG_WELL}) && defined($summary{$trace_plate_id}{PASS_QUAL_WELL}) ? $summary{$trace_plate_id}{NO_ATG_WELL}*100/$summary{$trace_plate_id}{PASS_QUAL_WELL} : 0),
      sprintf("%4.2f(%4.2f)", defined($summary{$trace_plate_id}{FAIL_IDENTITY_NOATG}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{FAIL_IDENTITY_NOATG}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0, defined($summary{$trace_plate_id}{FAIL_IDENTITY_NOATG_WELL}) && defined($summary{$trace_plate_id}{PASS_QUAL_WELL}) ? $summary{$trace_plate_id}{FAIL_IDENTITY_NOATG_WELL}*100/$summary{$trace_plate_id}{PASS_QUAL_WELL} : 0);
    print OUT "\n";
  }
  close OUT;
} # end of both dir

else { # single dir. trace
  open OUT, ">$result_dir$file_prefix" . "_summary.tab";
  print OUT join "\t", "PLATE_ID", "AVG_PHRED", "TOTAL_PICK", "PASS_QUAL", "PASS_BL2", "WILD_TYPE", "IDENTITY<100", "NO_ATG", "IDENTITY<100_and_NOATG", "%PASS_QUAL", "%PASS_BL2_vs_TOTAL_PICK", "%PASS_BL2_vs_PASS_QUAL", "%WILD_TYPE", "%IDENTITY<100", "%NO_ATG", "%IDENTITY<100_and_NOATG";
  print OUT "\n"; 
  foreach my $trace_plate_id (sort keys %summary) {
    print OUT join "\t", 
      $trace_plate_id, 
      defined($summary{$trace_plate_id}{AVG_PHRED})?$summary{$trace_plate_id}{AVG_PHRED}:"n/a",
      defined($summary{$trace_plate_id}{TOTAL_PICK})?$summary{$trace_plate_id}{TOTAL_PICK}:0,
      sprintf("%d", defined($summary{$trace_plate_id}{PASS_QUAL})? $summary{$trace_plate_id}{PASS_QUAL}: 0),
      sprintf("%d", defined($summary{$trace_plate_id}{PASS_BL2})? $summary{$trace_plate_id}{PASS_BL2}: 0),
      sprintf("%d", defined($summary{$trace_plate_id}{WILD_TYPE})? $summary{$trace_plate_id}{WILD_TYPE}: 0),
      sprintf("%d", defined($summary{$trace_plate_id}{FAIL_IDENTITY})? $summary{$trace_plate_id}{FAIL_IDENTITY}: 0),
      sprintf("%d", defined($summary{$trace_plate_id}{NO_ATG})? $summary{$trace_plate_id}{NO_ATG}: 0),
      sprintf("%d", defined($summary{$trace_plate_id}{FAIL_IDENTITY_NOATG})? $summary{$trace_plate_id}{FAIL_IDENTITY_NOATG}: 0),
      sprintf("%4.2f", defined($summary{$trace_plate_id}{PASS_QUAL}) && defined($summary{$trace_plate_id}{TOTAL_PICK}) ? $summary{$trace_plate_id}{PASS_QUAL}*100/$summary{$trace_plate_id}{TOTAL_PICK} : 0),
      sprintf("%4.2f", defined($summary{$trace_plate_id}{PASS_BL2}) && defined($summary{$trace_plate_id}{TOTAL_PICK}) ? $summary{$trace_plate_id}{PASS_BL2}*100/$summary{$trace_plate_id}{TOTAL_PICK} : 0),
      sprintf("%4.2f", defined($summary{$trace_plate_id}{PASS_BL2}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{PASS_BL2}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0),
      sprintf("%4.2f", defined($summary{$trace_plate_id}{WILD_TYPE}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{WILD_TYPE}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0),
      sprintf("%4.2f", defined($summary{$trace_plate_id}{FAIL_IDENTITY}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{FAIL_IDENTITY}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0),
      sprintf("%4.2f", defined($summary{$trace_plate_id}{NO_ATG}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{NO_ATG}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0),
      sprintf("%4.2f", defined($summary{$trace_plate_id}{FAIL_IDENTITY_NOATG}) && defined($summary{$trace_plate_id}{PASS_QUAL}) ? $summary{$trace_plate_id}{FAIL_IDENTITY_NOATG}*100/$summary{$trace_plate_id}{PASS_QUAL} : 0);
    print OUT "\n";
        
  }
  close OUT;
}

# print out plate detail info:
  $sql = "select PLATE_ID, PLATE_NAME, WELL_COUNT, MIN_ORF_SIZE, MAX_ORF_SIZE from $t->{PLATE_DETAIL}{TABLENAME} where plate_id in (" . join(",", @trace_plate_ids) . ") order by PLATE_ID";
  $query = $t->{PLATE_DETAIL}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  open OUT, ">$result_dir$file_prefix" . "_plate_detail.tab";
  print OUT join "\t", "PLATE_ID", "PLATE_NAME", "WELL_COUNT", "MIN_ORF_SIZE", "MAX_ORF_SIZE";
  print  OUT "\n";
  while (my ($PLATE_ID, $PLATE_NAME, $WELL_COUNT, $MIN_ORF_SIZE, $MAX_ORF_SIZE) = $query -> fetchrow()) {
    print OUT join "\t", $PLATE_ID, $PLATE_NAME, $WELL_COUNT, $MIN_ORF_SIZE, $MAX_ORF_SIZE;
    print OUT "\n";
  }
  $query -> finish;
  close OUT;

} # end of ANALYZE RESULT

} # end QC_helper

#*#################################################################
# SUBROUTINE: load same formatted alignment result file
# FUNCTION:   
#################################################################
sub load_sam_report {
  my ($pkg, $dbh, $dbname, $sam_prefix, $sam_result) = @_;
  my ($sql, $query, $rtn);

  my %t = (
    SAM_HEAD => {
      DBH => $dbh,
      DBNAME => $dbname,
      TABLENAME => $sam_prefix . "_SAMHEAD",
      COL_LIST => ordered_hash_ref(
        line_id => 'bigint(20) auto_increment',
        head_type => 'char(2)',
        head_def => 'mediumblob', 
      ),
      PRIMARY_KEY => 'line_id',
      CREATE_IF_NOT_EXIST => 1,
      NEED_BACKUP => 1,
    },

    SAM_ALN => {
      DBH => $dbh,
      DBNAME => $dbname,
      TABLENAME => $sam_prefix . "_SAMALN",
      COL_LIST => ordered_hash_ref(
        line_id => 'bigint(20) unsigned auto_increment',
        QNAME => 'varchar(255)',
        FLAG => ' smallint unsigned', 
        RNAME => 'varchar(255)',
        POS => 'bigint(29) unsigned',
        MAPQ => 'tinyint(3) unsigned',
        CIGAR => 'mediumblob',
        RNEXT => 'mediumblob',
        PNEXt => 'bigint(29) unsigned',
        TLEN => 'bigint(30)',
        SEQ => 'mediumblob', 
        QUAL => 'mediumblob',
        OPT => 'mediumblob', 
        numOpts => 'tinyint(2)',
      ),
      PRIMARY_KEY => 'line_id',
      CREATE_IF_NOT_EXIST => 1,
      NEED_BACKUP => 1,
    },     
  );
  
  MY_DButility -> prepare_tables(\%t);
  $sql = "insert into $t{SAM_HEAD}{TABLENAME} (head_type, head_def) values(?,?)";
  my $q_head = $t{SAM_HEAD}{DBH} -> prepare($sql);

  $sql = "insert into $t{SAM_ALN}{TABLENAME} (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, OPT, numOpts) values(?" . ",?"x12 . ")";
  my $q_aln = $t{SAM_HEAD}{DBH} -> prepare($sql);
  
  open IN, "<$sam_result";
  while(my $line = <IN>) {
    chomp $line;
# process header if present:
    if($line=~/^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/) {
	my ($type) = $line =~ /^\@([A-Za-z][A-Za-z]).+$/;
	my ($def) = $line =~ /^\@[A-Za-z][A-Za-z]\t(.+)\s*$/;
        $rtn = $q_head -> execute($type, $def);
    }
    elsif($line =~ /^\@CO\t.*/) {
        my ($type) = $line =~ /^\@(CO)\t.+$/;
	my ($def) =~ $line =~ /^\@CO\t(.+)\s*$/;
        $rtn = $q_head -> execute($type, $def);
    }
    else { # alignment line
	my @items = split("\t", $line);
        my @opt = splice(@items, 11);
        $rtn = $q_aln -> execute(@items, join("\t", @opt), scalar(@opt));
    }  
  }
  close IN;
  $q_head -> finish;
  $q_aln -> finish;

} #end load_sam_report

#*#################################################################
# SUBROUTINE: loadSamReport
# FUNCTION:  load as tab-limited file
#################################################################
sub load_sam_report {
  my ($pkg, $dbh, $dbname, $sam_table, $sam_result) = @_;
  my ($sql, $query, $rtn);
  my %t = ();

  if(MY_DButility->check_table_existence($dbh, $dbname, $sam_table)==1) { # if blast table already exists, print prompt and return;
      print "WARNING - $sam_table has existed already in $dbname, please remove it first if need to reroad.\n";
      return -1;
  }

  MY_BioUtilitySeq ->enlist_sam_table(\%t, 'SAM', $dbh, $dbname, $sam_table, 1, 1);
  MY_DButility -> prepare_tables(\%t);

  my $sql = "load data local infile '$sam_result' into table $sam_table  FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '\"' ESCAPED BY '\\' LINES TERMINATED BY '\n' (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, OP1, OP2, OP3, OP4, OP5, OP6, OP7, OP8, OP9, OP10, OP11, OP12, OP13, OP14, OP15, OP16)";
  $rtn = $dbh -> do($sql);

    $sql = "alter table $sam_table add index(QNAME)";
    $rtn = $dbh -> do($sql);

    $sql = "alter table $sam_table add index(RNAME)";
    $rtn = $dbh -> do($sql);

# remove empty columns:
  foreach my $c (16..1) {
    $sql = "select distinct OP" . $c . " from $sam_table";
    $query = $dbh->prepare($sql);
    $rtn = $query->execute();
    if ($rtn == 0) {
      $sql = "alter table $sam_table drop OP" . $c;
    }
  } #end of op. column check
 
} #end of load_sam_report

#*#################################################################
# SUBROUTINE: load agrep result for well index alignment
#################################################################
sub load_widx_report {
  my ($pkg, $dbh, $dbname, $widx_table, $widx_result) = @_;
  my ($sql, $query, $rtn);
  my %t = ();

  if(MY_DButility->check_table_existence($dbh, $dbname, $widx_table)==1) { # if table already exists, print prompt and return;
      print "WARNING - $widx_table has existed already in $dbname, please remove it first if need to reroad.\n";
      return -1;
  }

  MY_BioUtilitySeq ->enlist_widx_table(\%t, 'WIDX', $dbh, $dbname, $widx_table, 1, 1);
  MY_DButility -> prepare_tables(\%t);

  $sql = "load data local infile '$widx_result' into table $widx_table  FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '\"' ESCAPED BY '\' LINES TERMINATED BY '\n' (wellID, readID, readSeq)";
  $rtn = $dbh -> do($sql);

# note readID needs to be calculated:
  $sql = "update $widx_table set readID=(readID+2)/4";
  $rtn = $dbh -> do($sql);

  $sql = "alter table $widx_table add index(readID)";
  $rtn = $dbh -> do($sql);

  $sql = "alter table $widx_table add index(wellID)";
  $rtn = $dbh -> do($sql);

} #end of load_widx_report

#*#################################################################
# SUBROUTINE: load SOAP output result 
#################################################################
sub load_soap_report {
  my ($pkg, $dbh, $dbname, $soap_table, $soap_result) = @_;
  my ($sql, $query, $rtn);
  my %t = ();
 
  my $current_dir = Cwd::getcwd();
  if(MY_DButility->check_table_existence($dbh, $dbname, $soap_table)==1) { # if soap table already exists, print prompt and return;
      print "WARNING - $soap_table has existed already in $dbname, please remove it first if need to reroad.\n";
      return -1;
  }

  MY_BioUtilitySeq ->enlist_soap_table(\%t, 'SOAP', $dbh, $dbname, $soap_table, 1, 1);
  MY_DButility -> prepare_tables(\%t);

  my $sql = "load data local infile '$soap_result' into table $soap_table  FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '\"' ESCAPED BY '\' LINES TERMINATED BY '\n' (readID, readSeq, readQual, numHit, endType, readLen, strand, refName, refStart, matchType, mInfo1, mInfo2)";
  $rtn = $dbh -> do($sql);

  $sql = "update $soap_table set readID=readID+1";
  $rtn = $dbh -> do($sql);

  $sql = "alter table $soap_table add index(readID)";
  $rtn = $dbh -> do($sql);

  $sql = "alter table $soap_table add index(refName)";
#    $rtn = $dbh -> do($sql);

} #end of load_soap_report


#*#################################################################
# SUBROUTINE: load SOAP output result 
#################################################################
sub load_soap_report_bak {
  my ($pkg, $dbh, $dbname, $soap_table, $soap_result, $ext) = @_;
  my ($sql, $query, $rtn);
  my %t = ();

  if(MY_DButility->check_table_existence($dbh, $dbname, $soap_table)==1) { # if soap table already exists, print prompt and return;
      print "WARNING - $soap_table has existed already in $dbname, please remove it first if need to reroad.\n";
      return -1;
  }

  if(defined $ext && $ext==1) {
     MY_BioUtilitySeq ->enlist_soap_widx_table(\%t, 'SOAP', $dbh, $dbname, $soap_table, 1, 1);
  }
  else {
     MY_BioUtilitySeq ->enlist_soap_table(\%t, 'SOAP', $dbh, $dbname, $soap_table, 1, 1);
  } 
  MY_DButility -> prepare_tables(\%t);

  my $sql = "load data local infile '$soap_result' into table $soap_table  FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '\"' ESCAPED BY '\' LINES TERMINATED BY '\n' (readID, readSeq, readQual, numHit, endType, readLen, strand, refName, refStart, matchType, mInfo1, mInfo2)";
  $rtn = $dbh -> do($sql);

    $sql = "alter table $soap_table add index(readID)";
    $rtn = $dbh -> do($sql);

    $sql = "alter table $soap_table add index(refName)";
#    $rtn = $dbh -> do($sql);

} #end of load_soap_report_bak

#*#################################################################
# SUBROUTINE: loadBlastReport2
# FUNCTION:   improved version of load blast result, 2nd way, just 
#             load as tab-limited file
#################################################################
sub load_blast_report2 {
  my ($pkg, $dbh, $dbname, $blast_table, $blast_result, $best_hit) = @_;
  my ($sql, $query, $rtn);
  my %t = ();

  if(MY_DButility->check_table_existence($dbh, $dbname, $blast_table)==1) { # if blast table already exists, print prompt and return;
      print "WARNING - $blast_table has existed already in $dbname, please remove it first if need to reroad.\n";
      return -1;
  }

  MY_BioUtilitySeq ->enlist_blast_table(\%t, 'BLAST', $dbh, $dbname, $blast_table, 1, 1);
  MY_DButility -> prepare_tables(\%t);

  my $sql = "load data local infile '$blast_result' into table $blast_table  FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '\"' ESCAPED BY '\' LINES TERMINATED BY '\n' (QUERY, SUBJECT, PERCENT_IDENTITY, ALIGN_LEN, MISMATCH, GAP, Q_START, Q_END, S_START, S_END, E_VALUE, BIT_SCORE)";
  $rtn = $dbh -> do($sql);

    $sql = "alter table $blast_table add index(QUERY)";
    $rtn = $dbh -> do($sql);

    $sql = "alter table $blast_table add index(SUBJECT)";
    $rtn = $dbh -> do($sql);

# get best hit (only do it when required):
  if($best_hit == 1) {
    my $tmp_table = "tmp" . $blast_table . "_maxscore";
    $sql = "create table $tmp_table select QUERY, max(bit_score) as MAX_SCORE from $blast_table group by QUERY";
    $rtn = $dbh -> do($sql);
    $sql = "alter table $tmp_table add index(QUERY)";
    $rtn = $dbh -> do($sql);
 
    $sql = "create table $blast_table" . "_BEST select distinct a.* from $blast_table a inner join $tmp_table b using(QUERY) where a.BIT_SCORE=b.MAX_SCORE";
    $rtn = $dbh -> do($sql);
    
    $sql = "drop table $tmp_table";
    $rtn = $dbh -> do($sql);

  }
} #end of load_blast_report2

#*#################################################################
# SUBROUTINE: loadBlastReport
# FUNCTION:   load blast result, 2nd way, just load as tab-limited 
#             report file
#################################################################
sub load_blast_report {
  my ($pkg, $dbh, $dbname, $blast_table, $blast_result) = @_;
  my %col2insert = ();
  my ($sql, $query, $rtn);

  my %tableInfo = (
	TABLE_NAME=>$blast_table,    
  );
  my @colList = (
      {COL_NAME=>'QUERY',
       COL_TYPE=>'varchar(100)',
       COL_NO=>0
      },
      {COL_NAME=>'SUBJECT',
       COL_TYPE=>'varchar(100)',
       COL_NO=>1,
      },
      {COL_NAME=>'PERCENT_IDENTITY',
       COL_TYPE=>'decimal(5,2)',
       COL_NO=>2,
      },
      {COL_NAME=>'ALIGN_LEN',
       COL_TYPE=>'int(5)',
       COL_NO=>3,
      },
      {COL_NAME=>'MISMATCH',
       COL_TYPE=>'int(3)',
       COL_NO=>4,
      },
      {COL_NAME=>'GAP',
       COL_TYPE=>'int(3)',
       COL_NO=>5,
      },
      {COL_NAME=>'Q_START',
       COL_TYPE=>'int(6)',
       COL_NO=>6,
      },
      {COL_NAME=>'Q_END',
       COL_TYPE=>'int(6)',
       COL_NO=>7,
      },
      {COL_NAME=>'S_START',
       COL_TYPE=>'int(6)',
       COL_NO=>8,
      },
      {COL_NAME=>'S_END',
       COL_TYPE=>'int(6)',
       COL_NO=>9,
      },
      {COL_NAME=>'E_VALUE',
       COL_TYPE=>'double',
       COL_NO=>10,
      },
      {COL_NAME=>'BIT_SCORE',
       COL_TYPE=>'int(6)',
       COL_NO=>11,
      },     
  );

  @{$tableInfo{COL_LIST}}=@colList;
  my $delimiter = '\t';

  MY_BioUtilitySeq->load_table_from_linefile($dbh, $dbname, $delimiter, \%tableInfo, $blast_result);

# add Q_ORF_ID, S_ORF_ID, Q_ORF_SIZE, S_ORF_SIZE:
my $AUGMENT_TABLE= 0;
if($AUGMENT_TABLE == 1) {
  $sql = "alter table $dbname.$blast_table add (Q_ORF_ID int(8), S_ORF_ID int(8), Q_ORF_SIZE int(5),S_ORF_SIZE int(5))";
  $rtn = $dbh -> do($sql);

  $sql = "alter table $dbname.$blast_table add index (Q_ORF_ID)";
  $rtn = $dbh -> do($sql);

  $sql = "alter table $dbname.$blast_table add index (S_ORF_ID)";
  $rtn = $dbh -> do($sql);

# extract  Q_ORF_ID, S_ORF_ID, Q_ORF_SIZE, S_ORF_SIZE:
  $sql = "update $dbname.$blast_table set Q_ORF_ID=substring_index(substring_index(QUERY,'|',1), '=',-1), S_ORF_ID=substring_index(substring_index(SUBJECT,'|',1), '=',-1), Q_ORF_SIZE=substring_index(substring_index(QUERY,'|',2), '=',-1), S_ORF_SIZE=substring_index(substring_index(SUBJECT,'|',2), '=',-1)";
  $rtn = $dbh -> do($sql);

# get best hit:
  $sql = "create table $blast_table" . "_BEST select * from $blast_table group by QUERY, substring_index(QUERY, '|', 1)=substring_index(SUBJECT, '|', 1) order by QUERY, substring_index(QUERY, '|', 1)=substring_index(SUBJECT, '|', 1) desc, BIT_SCORE desc"; # find the best hit in both match and unmatch categories
  $rtn = $dbh -> do($sql);
} # end of augment table
  
} # end of loadBlastReport

##########################################################
# LOAD SPLIGN result
# ADDED DATE: 2011/02/28
##########################################################
sub load_splign_result {
  my ($pkg, $t, $rtInfo) = @_;
  my ($sql, $query, $rtn);

  my $result_dir = (defined $rtInfo -> {RESULT_DIR}) ? $rtInfo->{RESULT_DIR}:"";
  my $splign_result = $rtInfo -> {SPLIGN_RESULT};

my $LOAD_SPLIGN = 0;
if($LOAD_SPLIGN == 1) {
  $sql = "insert into $t->{SPLIGN}{TABLENAME} (COMP_DIR, COMP_ID, QUERY, SUBJECT, IDENTITY, ALIGN_LEN, Q_START, Q_END, S_START, S_END, ACCEPTOR_SITE, TYPE, DONOR_SITE, ALIGN_TRANSCRIPT) values (?" . ",?"x12 . ",?)";
  $query = $t->{SPLIGN}{DBH} -> prepare($sql);

  open IN, "<$result_dir" . $splign_result;
  while(my $line = <IN>) {
    last if $line =~ /^# END$/;
    chomp $line;
    my @items = split("\t", $line);
    my $comp_dir = substr($items[0],0,1);
    my $comp_id = int(substr($items[0],1));
    my ($acceptor_site, $type, $donor_site) = split("<|>", $items[9]);
    $rtn = $query -> execute($comp_dir, $comp_id, $items[1], $items[2], $items[3], $items[4], $items[5], $items[6], $items[7], $items[8], uc($acceptor_site), uc($type), uc($donor_site), $items[10]);     
  } #end of while
  close IN;
  $query -> finish;

# add index:
  $sql = "alter table $t->{SPLIGN}{TABLENAME} add index (COMP_ID, COMP_DIR)";
  $rtn = $t->{SPLIGN}{DBH} -> do($sql);

}

my $DO_CLEANUP = 0;
if($DO_CLEANUP==1) {
# clean up the result and load into SPLIGN_CLEAN table:
  $sql = "insert into $t->{SPLIGN_CLEAN}{TABLENAME} (COMP_DIR, COMP_ID, QUERY, SUBJECT, EXON_NR, IDENTITY, ALIGN_LEN, Q_START, Q_END, S_START, S_END, ACCEPTOR_SITE, TYPE, DONOR_SITE, ALIGN_TRANSCRIPT) values (?" . ",?"x13 . ",?)";
  my $qi_clean = $t->{SPLIGN}{DBH} -> prepare($sql);
  $sql = "select QUERY, COMP_DIR, COMP_ID, max(AVG_IDENTITY) from (select QUERY, COMP_DIR, COMP_ID, sum(IDENTITY*ALIGN_LEN)/sum(ALIGN_LEN) as AVG_IDENTITY from  $t->{SPLIGN}{TABLENAME} ";
#  $sql .= " where QUERY='2070_3254_0'";
  $sql .= " group by QUERY, COMP_DIR,COMP_ID order by QUERY, AVG_IDENTITY DESC) a group by QUERY";
  $query = $t->{SPLIGN}{DBH} -> prepare($sql);
  $rtn = $query -> execute();
  $sql = "select SUBJECT, IDENTITY, ALIGN_LEN, Q_START, Q_END, S_START, S_END, ACCEPTOR_SITE, TYPE, DONOR_SITE, ALIGN_TRANSCRIPT from   $t->{SPLIGN}{TABLENAME} where QUERY = ? and COMP_DIR = ? and COMP_ID = ? order by Q_START";
  my $q_exon = $t->{SPLIGN}{DBH} -> prepare($sql);

  while(my ($QUERY, $COMP_DIR, $COMP_ID, $AVG_IDENTITY) = $query -> fetchrow()) {
    $rtn = $q_exon -> execute($QUERY, $COMP_DIR, $COMP_ID);    
    my $exon_nr = 1;
    my $gap_nr = -1; # using negtive order number to identify gap
    while (my ($SUBJECT, $IDENTITY, $ALIGN_LEN, $Q_START, $Q_END, $S_START, $S_END, $ACCEPTOR_SITE, $TYPE, $DONOR_SITE, $ALIGN_TRANSCRIPT) = $q_exon->fetchrow()) {
      if($TYPE ne 'EXON') {
        $qi_clean -> execute($COMP_DIR, $COMP_ID, $QUERY, $SUBJECT, $gap_nr--, $IDENTITY, $ALIGN_LEN, $Q_START, $Q_END, $S_START, $S_END, $ACCEPTOR_SITE, $TYPE, $DONOR_SITE, $ALIGN_TRANSCRIPT);
      }
      else {
        $qi_clean -> execute($COMP_DIR, $COMP_ID, $QUERY, $SUBJECT, $exon_nr++, $IDENTITY, $ALIGN_LEN, $Q_START, $Q_END, $S_START, $S_END, $ACCEPTOR_SITE, $TYPE, $DONOR_SITE, $ALIGN_TRANSCRIPT);
      }
    } # end inner while exon
  } # end while

  $qi_clean -> finish;
  $query -> finish;
  $q_exon -> finish;

# add index:
  $sql = "alter table $t->{SPLIGN_CLEAN}{TABLENAME} add index (COMP_ID, COMP_DIR)";
  $rtn = $t->{SPLIGN}{DBH} -> do($sql);
} # end of DO_CLEANUP

# distinguish good/bad hit:
  my $bad_comp_id = "tmp_comp_id";
  $sql = "create table $bad_comp_id select distinct COMP_ID from $t->{SPLIGN_CLEAN}{TABLENAME} where TYPE != 'EXON'";
  $rtn = $t->{SPLIGN}{DBH} -> do($sql); 
  $sql = "alter table $bad_comp_id add index (COMP_ID)";
  $rtn = $t->{SPLIGN}{DBH} -> do($sql);
  
  $sql = "create table $t->{SPLIGN_CLEAN}{TABLENAME}" . "_GOOD like $t->{SPLIGN_CLEAN}{TABLENAME}";
  $rtn =  $t->{SPLIGN_CLEAN}{DBH} -> do($sql);
  $sql = "create table $t->{SPLIGN_CLEAN}{TABLENAME}" . "_BAD like $t->{SPLIGN_CLEAN}{TABLENAME}";
  $rtn =  $t->{SPLIGN_CLEAN}{DBH} -> do($sql);
  
#  $sql = "create table $t->{SPLIGN_CLEAN}{TABLENAME}" . "_GOOD select * from $t->{SPLIGN_CLEAN}{TABLENAME} where COMP_ID not in (select COMP_ID from $bad_comp_id)";
  $sql = "insert into $t->{SPLIGN_CLEAN}{TABLENAME}" . "_GOOD select * from $t->{SPLIGN_CLEAN}{TABLENAME} where COMP_ID not in (select COMP_ID from $bad_comp_id)";
  $rtn =  $t->{SPLIGN_CLEAN}{DBH} -> do($sql);

#  $sql = "alter table $t->{SPLIGN_CLEAN}{TABLENAME}" . "_GOOD add index (COMP_ID, COMP_DIR)";
#  $rtn = $t->{SPLIGN}{DBH} -> do($sql);

  $sql = "insert into $t->{SPLIGN_CLEAN}{TABLENAME}" . "_BAD select * from $t->{SPLIGN_CLEAN}{TABLENAME} where COMP_ID  in (select COMP_ID from $bad_comp_id)";
  $rtn =  $t->{SPLIGN_CLEAN}{DBH} -> do($sql);
#  $sql = "alter table $t->{SPLIGN_CLEAN}{TABLENAME}" . "_BAD add index (COMP_ID, COMP_DIR)";
#  $rtn = $t->{SPLIGN}{DBH} -> do($sql);

  $sql = "drop table $bad_comp_id";
  $rtn = $t->{SPLIGN}{DBH} -> do($sql);
} # end load_splign_result 

##########################################################
# RUN SPLIGN UTILITY
# input: input fasta, 
#        input database
# output: splign result
# rtInfo: result_dir
#         input_dir
#         bin_dir # where splign utility located
##########################################################
sub run_splign {
  my ($pkg, $t, $rtInfo) = @_;
  my ($sql, $query, $rtn);
  
  my $result_dir = (defined $rtInfo->{RESULT_DIR}) ? $rtInfo->{RESULT_DIR}:"";
  my $subject_dir = (defined $rtInfo->{SUBJECT_DIR}) ? $rtInfo->{SUBJECT_DIR}:"";
  my $bin_dir = (defined $rtInfo->{BIN_DIR}) ? $rtInfo->{BIN_DIR}:"";
  my $subject_fa = $rtInfo -> {SUBJECT_FASTA};

# create temp fasta dir:
  my $fasta_dir = $result_dir . "fasta/";
  if(!-e $fasta_dir) {
    system("mkdir -p $fasta_dir");
  } 
  my $query_fa = defined($rtInfo->{QUERY_FASTA}) ? $rtInfo->{QUERY_FASTA} : "query.fa";

# create splign dir:
  my $splign_dir = $result_dir . "splign/";
  if(!-e $splign_dir) {
    system("mkdir -p $splign_dir");
  } 
  my $splign_result =  $query_fa . ".splign";

  my $script_dir = Cwd::getcwd();
  chdir($fasta_dir);

# link subject db:
my $LINK_SUBJECT = 0;
if($LINK_SUBJECT == 1) {
  my $filelist = 'temp.list';
  system("find $subject_dir -name '" . $subject_fa . "*'>$filelist");
  open IN, "<$filelist";
  while(<IN>) {
    chomp;
    my ($fa) = $_ =~ /\/([^\/]+?)$/;
    system("ln -s $_ $fa"); 
  } # end while
  close IN;
  system("rm $filelist");
} # end of creating symlink

my $DO_QUERY_FASTA = 0;
if($DO_QUERY_FASTA==1) {

# get query fasta:
  $sql = "select $t->{QUERY}{QUERY_NAME_FIELD},$t->{QUERY}{QUERY_SEQ_FIELD} from $t->{QUERY}{TABLENAME} a ";
  $sql .= "where $t->{QUERY}{FILTER}" if defined($t->{QUERY}{FILTER});
  $query = $t->{QUERY}{DBH} -> prepare($sql);
  $rtn = $query -> execute();  
  open OUT, ">$query_fa";
  while (my ($query_name, $query_seq) = $query->fetchrow()) {
    print OUT ">$query_name";
    print OUT "\n";
    print OUT "$query_seq";
    print OUT "\n";
  } # end of while
  $query -> finish;
  close OUT;
  system("formatdb -pF -oT -i $query_fa");
}

my $DO_SPLIGN = 0;
if($DO_SPLIGN==1) {
#   system("splign -mklds .");
#   system("compart -qdb $query_fa -sdb $subject_fa > $query_fa" . ".compart"); 
#   system("splign -ldsdir .  -comps $query_fa" . ".compart > ../splign/$splign_result");

   system($bin_dir ."splign64 -mklds .");
   system($bin_dir ."compart64 -qdb $query_fa -sdb $subject_fa > $query_fa" . ".compart"); 
   system($bin_dir ."splign64 -ldsdir .  -comps $query_fa" . ".compart > ../splign/$splign_result");

} # end of DO_SPLIGN

  chdir($script_dir);
  
# load splign result to db:
  my %splign_rtInfo = (
    RESULT_DIR => $splign_dir,
    SPLIGN_RESULT => $splign_result,
  );

# prepare splign:
  MY_BioUtilitySeq -> enlist_splign_table($t, 'SPLIGN', $t->{QUERY}{DBH}, $t->{QUERY}{DBNAME}, $rtInfo->{RESULT_TABLE_PREFIX} . "_SPLIGN_ORIG", 1, 1);
  MY_BioUtilitySeq -> enlist_splign_clean_table($t, 'SPLIGN_CLEAN', $t->{QUERY}{DBH}, $t->{QUERY}{DBNAME}, $rtInfo->{RESULT_TABLE_PREFIX} . "_SPLIGN_CLEAN", 1, 1);
  
# prepare tables:
  MY_DButility -> prepare_tables($t);

  MY_SeqAlign -> load_splign_result($t, \%splign_rtInfo); 

} #end run_splign


##########################################################
#  BEGIN of utility function
##########################################################
sub ordered_hash_ref {
    tie my %hash, 'Tie::IxHash', @_;
    return \%hash;
}

#########################################################
# customized my_tiled_hsps()
#########################################################
sub my_tiled_hsps {
  my ($hit) = @_;
  my @hsp_list_f=();
  my @hsp_list_r=();
  my @seq_ind_f = ();
  my @seq_ind_r = ();

  while(my $hsp = $hit->next_hsp()) {
    if($hsp->strand('query') == 1 && $hsp->strand('hit') == 1) { 
       push @hsp_list_f, (
          Q_STRAND => $hsp->strand('query'), 
          S_STRAND => $hsp->strand('hit'),
          Q_START => $hsp->{QUERY_START},
          Q_END => $hsp->{QUERY_END},
          S_START => $hsp->{HIT_START},
          S_END => $hsp->{HIT_END},
       );
    }
    elsif ($hsp->strand('query') == 0 && $hsp->strand('hit') == 0) { 
       push @hsp_list_r, (
          Q_STRAND => $hsp->strand('query'), 
          S_STRAND => $hsp->strand('hit'),
          Q_START => $hsp->{QUERY_START},
          Q_END => $hsp->{QUERY_END},
          S_START => $hsp->{HIT_START},
          S_END => $hsp->{HIT_END},
          S_INDS => $hsp->seq_inds('HIT'),
          Q_INDS => $hsp->seq_inds('QUERY'),
       );
    }
  }
  
  my %seq_inds = ();
  $seq_inds{HIT} = 0;
  $seq_inds{QUERY} = 0;

  foreach my $hsp_info (@hsp_list_f) {
#    $seq_inds{HIT} = 
  }


  
  
} # end my_tiled_hsps

##########################################################
# getAlignmentGap for blast result
##########################################################
sub getAlignmentGap {
  my ($seqInd) = @_;
  my %retval = ();
  
  my $curSeqLen = 0;
  my $curGap = 0;
  my $curPos = 0;
  my $curSign = undef;
  my $min_gap = MAX_GAP; 
  my $max_gap = MIN_GAP;

  my @items = split(/(-|::)/, $seqInd);
 

  foreach my $item (@items) {
    if($item =~ /\d+/) { # number
      if($curSign eq "-") {
	$curSeqLen += $item-$curPos;
        $curPos = $item;
      }
      elsif ($curSign eq "::") {
        $curGap = $item-$curPos-1;
        if($curGap<0) {
	  printf "getAlignmentGap: something weird happened at $item\n";
	  $curGap = 0;
          return;
        } 
        elsif ($curGap==0) { # actual continuous, no gap
	  $curSeqLen +=1;
          $curPos = $item;
        }
        else { # gap>0
	  push @{$retval{SEQLEN}}, $curSeqLen;
          push @{$retval{SEQGAP}}, $curGap; 
          $curPos = $item;
	  $curSeqLen=1; #reset segment
          if ($curGap>$max_gap) {
	    $max_gap = $curGap;
          }
	  if ($curGap<$min_gap) {
	    $min_gap = $curGap;
          }
        }
      }
      else { # undef
	$curPos = $item;
        $curSeqLen = 1;	#the len must be at least 1
      } 
    }
    elsif ($item =~ /-/) { # alignment
      $curSign = $item;
    }
    elsif ($item =~ /::/) { #gap
      $curSign = "::";
    }
    else {
      printf "unexpected character in seqInd text : $item\n";
    }
  } 

# record last segment
  if($curSign ne "-") {
    printf "getAlignmentGap: end with single number, probably truncated.\n";
  }
  else {
    push @{$retval{SEQLEN}}, $curSeqLen;
    $retval{MIN_GAP} = $min_gap;
    $retval{MAX_GAP} = $max_gap;
  }

  return \%retval;

} # end of getAlignmentGap

##########################################################
#  END of utility function
##########################################################
