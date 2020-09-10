#!/usr/bin/perl
use strict; 
use Bio::SearchIO;
 
my $fh;
my $fasta   = "/Users/donatcl1/Downloads/fasta-35.4.11/bin/fasta35";
my $library = $ARGV[1];
my $query   = $ARGV[0];
my $options = "-d 1 -b 1 -H -Q";
my $command = "$fasta $options $query $library";
 
open $fh,"$command |" || die("cannot run fasta cmd of $command: $!\n");
 
my $searchio  = Bio::SearchIO->new(-format => 'fasta', -fh => $fh);

          print "Query",   "\t",
            " Hit",        "\t",
            " Length",     "\t",
            " Length q",   "\t",
            " Length_h",   "\t",
            " Percent_id", "\t", 
            " Start q",    "\t",
            " Start h",    "\t",
            " Stop q",    "\t",
            " Stop h","\t",
            " SGGGSGR","\n";


while( my $result = $searchio->next_result ) {
  ## $result is a Bio::Search::Result::ResultI compliant object
  while( my $hit = $result->next_hit ) {
    ## $hit is a Bio::Search::Hit::HitI compliant object
    while( my $hsp = $hit->next_hsp ) {
      ## $hsp is a Bio::Search::HSP::HSPI compliant object
      if ( $hsp->percent_identity >= 50 ) {
      	      
      	      print $result->query_name,"\t",
      	      $hit->name,"\t",
      	      $hsp->length('total'),"\t",
      	      $hsp->length('query'),"\t",
      	      $hsp->length('hit'),"\t",
      	      $hsp->percent_identity,"\t", 
      	      $hsp->start("query"),"\t",
      	      $hsp->start("hit"),"\t",
      	      $hsp->end("query"),"\t",
      	      $hsp->end("hit"),
      	      "\t";
      	      my $seq_q=$hsp->query_string;
      	      $seq_q=~ s/^\s+//;
      	      $seq_q=~ s/-+//;
 #     	      print $seq_q,"\n";
      	      my $ind_end=index $seq_q,"SGGGSGR";
      	      print $ind_end,"\n";
      }
    }  
  }
}

