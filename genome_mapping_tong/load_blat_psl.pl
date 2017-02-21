#!/usr/local/bin/perl

use DBI;
use strict;
use Switch;
use Getopt::Std;
use Cwd;
use lib '/home/tong/lib/';
use MY_DButility;
use MY_BioUtilitySeq;
use MY_SeqAlign;
use MY_PrimerDesign;
use MY_PrimerOrder;
use MY_ProtUtility;



my $user = "tong";
my $passwd = "70n9";
my $host = 'paros';
my $db = "horfeome_annotation";

my $dbh = DBI->connect("dbi:mysql:$db;host=$host", $user, $passwd, {RaiseError => 1});

my %t;
 

MY_BioUtilitySeq -> enlist_blat_table(\%t, 'BLAT', $dbh, $db, 'orf_utr_hg38_blat_best', 0, 1);
MY_BioUtilitySeq -> enlist_blat_exon_table(\%t, 'BLAT_EXON', $dbh, $db, 'orf_utr_hg38_blat_exon', 0, 1);
my %rtInfo = (
    RESULT_DIR => '/home/tong/horfeome/annotation/blat/',
    BLAT_RESULT => 'all_seq.psl',
  );

MY_DButility -> prepare_tables(\%t);

MY_SeqAlign -> load_blat_result(\%t, \%rtInfo);


$dbh ->disconnect ;
