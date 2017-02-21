package MY_DButility;
use strict;
use DBI;
use POSIX qw(ceil floor);

1;
########################################################
# this script is to put together all the templates for
# tables used in various projects;
########################################################

# create on 2009-06-01
sub connect_db {
  my ($pkg, $db_server, $db_to_work_on, $user, $passwd) = @_;
  my $db_path = "dbi:mysql:$db_to_work_on;host=$db_server";
  return  DBI->connect($db_path, $user, $passwd, {RaiseError => 1});
}

# create on 2009-06-01
sub create_table_on_db {
  my ($pkg,$db, $template_name, $fields_hash, $primary_key) = @_;
  my $table_desc = "";
  while (my ($key, $value) = each %{$fields_hash}) {
          $table_desc .= "," if($table_desc ne "");
          $table_desc .= "$key $value\n";
  }
  my $sql = "CREATE TABLE if not exists $template_name (\n"
           . $table_desc;
  $sql .= ", PRIMARY KEY (" . $primary_key . ")" if defined ($primary_key);
  $sql .= "\n);";
  my $query = $db->prepare($sql);
  $query->execute();
}

# create on 2009-06-01
sub copy_table_struct_from {
  my ($pkg, $db, $src_dbname, $template_name) = @_;
  my $sql = "CREATE TABLE if not exists $template_name like $src_dbname.$template_name;";
  my $query = $db->prepare($sql);
  $query->execute();
}

# create on 2009-06-03
sub copy_table_with_data_from {
  my ($pkg, $db, $src_dbname, $template_name) = @_;
  my $sql = "CREATE TABLE if not exists $template_name like $src_dbname.$template_name;";
  my $query = $db->prepare($sql);
  $query->execute();
  $sql = "INSERT INTO $template_name SELECT * FROM $src_dbname.$template_name;";
  $query = $db->prepare($sql);
  $query->execute();
}

# create on 2009-06-01
sub rename_table {
  my ($pkg,$db, $orig_template_name, $new_template_name) = @_;
  my $sql = "ALTER TABLE $orig_template_name rename to $new_template_name;";
  my $query = $db->prepare($sql);
  $query->execute();
}

# create on 2009-06-03
sub add_column_to_table {
  my ($pkg,$db, $template_name, $col_name, $type, $modifier) = @_;
  my $sql = "ALTER TABLE $template_name ADD $col_name $type \n";
  while( my ($k, $v) = each %{$modifier}) {
    $sql .= " $k $v \n";
  }
  my $query = $db->prepare($sql);
  $query->execute();
}

# create on 2009-06-03
sub insert_row {
  my ($pkg,$db, $template_name, $fields) = @_;
  my $sql = "INSERT IGNORE INTO $template_name (" .
            join(",", keys %{$fields}) . ") VALUES (" . 
            join(",", values %{$fields}) . ")";
  my $query = $db->prepare($sql);
  my $rtn = $query->execute();
  $query->finish();
  return $rtn;
}

# create on 2010-06-08
# called after insertion to get the auto_increment id
sub last_insert_id {
  my ($pkg,$db) = @_;
  my $sql = "SELECT last_insert_id()";
  my $query = $db->prepare($sql);
  $query -> execute();
  my $last_id = $query->fetchrow();
  $query->finish;
  return $last_id;
}

sub delete_table {
  my ($pkg,$db, $template_name) = @_;
  my $sql = "DROP TABLE $template_name";
  my $query = $db->prepare($sql);
  $query->execute();
}

# create on 2009-12-09
sub update_row {
  my ($pkg,$db, $template_name, $fields, $filter) = @_;
  my @update_columns = ();
  foreach my $k (keys %{$fields}) {
    push @update_columns, join('=', $k, $fields->{$k});
  }
  my $sql = "UPDATE $template_name set " . join(",", @update_columns);
  $sql .= " where $filter " if defined $filter && $filter ne "";
  my $query = $db->prepare($sql);
  $query->execute();
}


###########################################################################
# Helper : get_filter_str
###########################################################################
sub get_filter_str {
     my ($pkg, $filter_list) = @_;

     my @filter_items = ();
     foreach my $col_name (keys %{$filter_list}) {
       push @filter_items, "$col_name = $filter_list->{$col_name}";
     }
     my $filter_str = join (" and ", @filter_items);

     return $filter_str;
}

###############################################################################
# sub check_db_existence - check if the database exists
# INPUT: $dbh - db handler
#        $dbname - the database name for the table to check existence
# OUTPUT: 0 not exists 1 exists
###############################################################################
sub check_db_existence {
   my ($pkg, $dbh, $dbname, $table_name) = @_;
   my ($sql, $query, $count);

   $sql = "SELECT count(*) FROM information_schema.schemata WHERE schema_name = '"
        . $dbname
        . "'" #remove possible prefix in table name
        . ";";
   $query = $dbh->prepare($sql);
   $query->execute();
   $count=$query->fetchrow();
   return $count; # 1: exist 0:not exist
}

###############################################################################
# sub check_table_existence - check if the table to be created already exists 
#        in database
# INPUT: $dbh - db handler
#        $dbname - the database name for the table to check existence
#        $table_name - the table needs to check existence, 
#            can be long name if it is not directly assoc. w/ db handler
# OUTPUT: 0 not exists 1 exists
###############################################################################
sub check_table_existence {
   my ($pkg, $dbh, $dbname, $table_name) = @_;
   my ($sql, $query, $count);

   $sql = "SELECT count(*) FROM information_schema.tables WHERE table_schema = '"
        . $dbname
        . "' AND table_name = substring_index('$table_name', '.', -1)" #remove possible prefix in table name
        . ";";
   $query = $dbh->prepare($sql);
   $query->execute();
   $count=$query->fetchrow();
   return $count; # 1: exist 0:not exist
}

###############################################################################
# sub check_column_existence - check if the column in the table to be created 
#        already exists in database
# INPUT: $dbh - db handler
#        $dbname - the database name for the table to check existence
#        $table_name - the table needs to check existence, 
#            can be long name if it is not directly assoc. w/ db handler
#        $column_name
# OUTPUT: 0 not exists 1 exists
###############################################################################
sub check_column_existence {
   my ($pkg, $dbh, $dbname, $table_name, $column_name) = @_;
   my ($sql, $query, $count);

   $sql = "SELECT count(*) FROM information_schema.columns WHERE table_schema ='" 
        . $dbname 
        . "' AND table_name = substring_index('$table_name', '.', -1) " #remove possible prefix in table name
        . " AND column_name = substring_index('$column_name', '.', -1) " #remove possible prefix in column name
        . ";";
   $query = $dbh->prepare($sql);
   $query->execute();
   $count=$query->fetchrow();
   return $count; # 1: exist 0:not exist
}

########################################################
# Helper : Assign column value
########################################################

sub assign_col_value {
  my ($pkg, $col_type, $col_value) = @_;
 
  if($col_type =~ /int|float|double|decimal|numeric/) {
    return $col_value;
  }  
  elsif ($col_type =~ /char|text|enum/) {
    return "'$col_value'";
  }
  elsif ($col_type =~ /datetime/) {
    return "str_to_date('$col_value', '%a %b %e %H:%i:%s %Y')";
  }
}

#######################################################
# batch check table existence/augment columns
#######################################################
sub prepare_tables {
  my ($pkg, $tables) = @_;
  
  my ($sql, $query, $rtn);
  my ($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time);
  my $curdate = sprintf("%04d%02d%02d", $Year+1900, $Month+1, $Day);

  foreach my $t (keys %{$tables}) {
# first check if table already exists:
    if (!MY_DButility->check_table_existence($tables->{$t}{DBH}, $tables->{$t}{DBNAME},$tables->{$t}{TABLENAME}) ) { # table not exists
      if( $tables->{$t}{CREATE_IF_NOT_EXIST}==1) { # need to creat new table
        $rtn = MY_DButility->create_table_on_db($tables->{$t}{DBH}, $tables->{$t}{DBNAME} . "." . $tables->{$t}{TABLENAME}, $tables->{$t}{COL_LIST}, $tables->{$t}{PRIMARY_KEY});
      }
      else {
        print "Warning -- Table $tables->{$t}{TABLENAME} does not exist.\n";
      }
    } # table not existence
    else { # table already exists
# check if need to back up, only back up with the first time of the day:
      my $backup_table = join("_", $tables->{$t}{TABLENAME},'bak', $curdate);
      if($tables->{$t}{NEED_BACKUP} == 1 &&
        !MY_DButility->check_table_existence($tables->{$t}{DBH}, $tables->{$t}{DBNAME},$backup_table)
      ) {
        $sql = "create table $tables->{$t}{DBNAME}.$backup_table select * from $tables->{$t}{DBNAME}.$tables->{$t}{TABLENAME}";
        $rtn = $tables->{$t}{DBH} -> do($sql); 
      }
     
# augment column if any
      foreach my $c (keys %{$tables->{$t}{COL_LIST}}) {
        if( ! MY_DButility->check_column_existence($tables->{$t}{DBH}, $tables->{$t}{DBNAME}, $tables->{$t}{TABLENAME}, $c)) {
           $rtn = MY_DButility->add_column_to_table($tables->{$t}{DBH}, $tables->{$t}{DBNAME} . "." . $tables->{$t}{TABLENAME}, $c, $tables->{$t}{COL_LIST}{$c}, undef);
        } # end of adding column
      } # end of $c
    } # table already exists
  } # end of $t
} # end of prepare_tables

#######################################################
# batch analyze data tables
#######################################################
sub analyze_tables {
  my ($pkg, $tables, $outfile) = @_;
  
  my ($sql, $query, $rtn);

  open OUT, ">$outfile";
  print OUT "THIS IS the analysis report for the following tables.\n";
  print OUT "="x80 . "\n\n";

  foreach my $t (keys %{$tables}) {
    if (!MY_DButility->check_table_existence($tables->{$t}{DBH}, $tables->{$t}{DBNAME},$tables->{$t}{TABLENAME})) { # table not exists
      print OUT "Warning - TABLE $t doesn't exist\n";
    }
    else { # table already exists

# get total count for the entire table:

      $sql = "select count(*) from $tables->{$t}{TABLENAME}";
      $query = $tables->{$t}{DBH} -> prepare($sql);
      $query -> execute();
      my $count = $query -> fetchrow();
      print OUT "total records in TABLE $tables->{$t}{TABLENAME}\t$count\n";
      $query -> finish;

      foreach my $c (keys %{$tables->{$t}{COL_LIST}}) {
        if( !MY_DButility->check_column_existence($tables->{$t}{DBH}, $tables->{$t}{DBNAME}, $tables->{$t}{TABLENAME}, $c)) {
          print OUT "Warning - Column $c is not in TABLE $tables->{$t}{TABLENAME}.\n";
        } # end of adding column
        else {
          $sql = "select count(distinct $c) as count_$c";
          if($tables->{$t}{COL_LIST}{$c} =~ /int|float|double|decimal|numeric|datetime/) {
             $sql  .= ", max($c) as max_$c, min($c) as min_$c";
          }  
          elsif ($tables->{$t}{COL_LIST}{$c} =~ /char|text|enum/) {
             $sql  .= ", max(length($c)) as maxlen_$c, min(length($c)) as minlen_$c";
          }
          else {
            print OUT "Warning - Column $c has unrecognized data type\n";
            next;
          }
          $sql .= " from $tables->{$t}{TABLENAME}";
          $query = $tables->{$t}{DBH} -> prepare($sql);
          $rtn = $query -> execute();
          my ($count_c, $max_c, $min_c) = $query->fetchrow();
          print OUT "Total Unique Values for Col $c\t$count_c\t";
          if($tables->{$t}{COL_LIST}{$c} =~ /int|float|double|decimal|numeric|datetime|enum/) {
              print OUT "max value\t$max_c\tmin_value\t$min_c";
          }  
          elsif($tables->{$t}{COL_LIST}{$c} =~ /char|text|enum/) {
              print OUT "max length\t$max_c\tmin_length\t$min_c";
          }
          print OUT "\n";
          $query -> finish;
           
          if($count_c<20) { # if total values is not many, type the stat.
            print OUT "Categories by $c :\n";
            $sql = "select $c, count(*) from $tables->{$t}{TABLENAME} group by $c order by $c";
            $query = $tables->{$t}{DBH} -> prepare($sql);
            $rtn = $query -> execute();
            while (my ($value_c, $count_by_c) = $query->fetchrow) {
              print OUT "$value_c\t$count_by_c\n";
            }
            $query -> finish;
          }           
        }
      } # end of $c
    }
  } # end of $t

  $query -> finish;
  close OUT;

} # end of analyze_tables


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$ BEGIN of Date Time functions
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#####################################################################
# getdate
#####################################################################
sub getdate {
  my ($pkg) = @_;

  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
  return sprintf("%d%02d%02d", $year+1900, $mon+1, $mday+1);  
}

#********************************************************************
#* END of Date Time functions
#*********************************************************************


#######################################################################
# utility : kill process 
#######################################################################
sub kill_hanging_process {
  my ($pkg, $dbh, $my_user, $my_command) = @_;
  my ($sql, $query, $rtn);

  $query = $dbh->prepare("SHOW FULL PROCESSLIST");
  $query->execute();
  my $query_kill = $dbh->prepare("KILL ? ");

  while (my ($process_id, $user, $host, $db, $command, $time, $state, $info)= $query->fetchrow()) {
    printf join ',',  $process_id, $user, $host, $db, $command, $time, $state, $info, "\n";
    if ($user eq $my_user && uc($command) eq uc($my_command) ) {
#      $sql="KILL $process_id";
      $query_kill->execute($process_id);    
    }
  }
  $query->finish();
  $query_kill->finish();
}
