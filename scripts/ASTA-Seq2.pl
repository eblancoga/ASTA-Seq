#!/usr/bin/perl -w                                
#                                                 #
# Script ASTA-Seq2.pl                             #
#                                                 #
#                                                 #
# Input : List of reads with gene,SNP,CG (COMPACT)# 
# Output: totals per gene and SNP/CG (FINAL)      #
#                                                 #
#                                                 #
# by Enrique Blanco (2024)                        #
###################################################

use strict;
use Getopt::Std;
use Term::ANSIColor;

#DEFINEs
my $TRUE = 1;
my $FALSE = 0;
my $PSEUDOCOUNT = 0.001;
my $NOTFOUND = -1;
my $NOSEQ ="*";
my $NOSPECIES="X";
my $LINES_INFO=100000;
my $TOKEN = "_SEPARATOR_";


## Step 1. Reading arguments
my %opt;

(getopts('v',\%opt)) or print_error("parser: Problems reading options\n");

print_mess("ASTA-Seq2.pl by Enrique Blanco (CRG 2024)\n\n");
print_mess("Stage 1.  Reading options");

my ($compact_file);
my $n_files;

$n_files = $#ARGV+1;
($n_files == 1) or print_error("USAGE: One parameter (compact file of results) is required but $n_files are provided!");
($compact_file) = @ARGV;

print_ok();
##


## Step 2 Summarizing the COMPACT lines of results per gene
my $line;
my $n_lines;
my $n_notfound;
my @record;
my ($readid,$gene,$SNP,$CG,$species);
my $key;
my (%voidCG,%voidSNP,%xSNP,%GENE);
my $key2;
my %totalOK;


print_mess("Stage 2.  Loading the COMPACT file values ($compact_file)\n");

open(FILE,$compact_file) or print_error("FILE $compact_file can not be opened");
$n_lines=0;
$n_notfound = 0;
while($line=<FILE>)
{
    ###print ">>> $line";

    if ($line =~ /FOUND/)
    {
	###print "NOT FOUND";
	$n_notfound++;
    }
    else
    {
	@record = split(/\t+/,$line);

	$readid = $record[0];
	$gene = $record[1];
	$SNP = $record[7];
	$CG = $record[11];
	$species = $record[14];
	chomp($species);
	
	# no CG flanking sequence was found
	if ($CG eq "*")
	{
	    if (!defined($voidCG{$gene}))
	    {
		$voidCG{$gene} = 1;
	    }
	    else
	    {
		$voidCG{$gene}++;
	    }
	}
	else
	{
	    # no species flanking sequence was found or SNP REF did not fit well
	    if ($SNP eq "*")
	    {
		if (!defined($voidSNP{$gene}))
		{
		    $voidSNP{$gene} = 1;
		}
		else
		{
		    $voidSNP{$gene}++;
		}
	    }
	    else
	    {
		if ($species eq "X")
		{
		    if (!defined($xSNP{$gene}))
		    {
			$xSNP{$gene} = 1;
		    }
		    else
		    {
			$xSNP{$gene}++;
		    }
		}
		else
		{
		    ###print join("\t",$gene,$CG,$species),"\n";	    
		    $key = $gene.$TOKEN.$CG;
		    if (!defined($GENE{$key}))
		    {
			$GENE{$key}{$species} = 1;
		    }
		    else
		    {
			$GENE{$key}{$species}++;
		    }
		    
		    #register total reads ok per gene
		    if (!defined($totalOK{$gene}))
		    {
			$totalOK{$gene} = 1;
		    }
		    else
		    {
			$totalOK{$gene}++;
		    }
		}
	    }
	}
    }    

    if ($n_lines % $LINES_INFO==0)
    {
	print_mess("...$n_lines\n");
    }
    
    $n_lines++;    

}
close(FILE);

# FINAL VALUES (CG per gene divided into MUS and CAS)
foreach $key (sort keys(%GENE))
{
    @record = split($TOKEN,$key);
    $gene = $record[0]; 
    $CG = $record[1];
    print $gene,"\t",$CG,"\t";
    
    for $key2 (sort keys %{$GENE{$key}})
    {
	print $key2,"\t",$GENE{$key}{$key2},"\t",($GENE{$key}{$key2}/$totalOK{$gene})*100,"\t";
    }
    print "\n";
}
# wrong lines per gene
foreach $key (sort keys(%voidSNP))
{
    print "noSNP: ",$key,"\t",$voidSNP{$key},"\n";
}
foreach $key (sort keys(%voidCG))
{
    print "noCG: ",$key,"\t",$voidCG{$key},"\n";
}
foreach $key (sort keys(%xSNP))
{
    print "SNP REF no fit: ",$key,"\t",$xSNP{$key},"\n";
}

print_mess("\n");

print_mess("$n_lines reads from $compact_file have been acquired\n");
print_mess("$n_notfound reads are NOT FOUND (no gene can be associated)");
print_ok();


## Step X. Finishing successful program execution

print_mess("Successful termination:");
print_ok();
exit(0);
##


############ Subroutines

sub print_mess
{
        my @mess = @_;

        print STDERR color("bold green"),"%%%% @mess" if (exists($opt{v}));
	print STDERR color("reset");
}

sub print_error
{
        my @mess = @_;

        print STDERR color("bold green"),"%%%% @mess\n";
	print STDERR color("reset");
	exit();
}

sub print_ok
{
    if (exists($opt{v}))
    {
	print STDERR color("bold green"), "\t\t[OK]\n\n";
	print STDERR color("reset");
    }
}

