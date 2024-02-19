#!/usr/bin/perl -w                                
#                                                  #
# Script ASTA-Seq1                                 #
#                                                  #
#                                                  #
# Input : Configuration file,FASTQ file            # 
# Output: Read assignment of gene,SNP and CG       #
#                                                  #
#                                                  #
# by Enrique Blanco (2024)                         #
####################################################

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


## Step 1. Reading arguments
my %opt;

(getopts('v',\%opt)) or print_error("parser: Problems reading options\n");

print_mess("ASTA-Seq1.pl by Enrique Blanco (CRG 2024)\n\n");
print_mess("Stage 1.  Reading options");

my ($config_file,$reads_file);
my $n_files;

$n_files = $#ARGV+1;
($n_files == 2) or print_error("USAGE: Two parameters (configuration file and FastQ file) are required but $n_files are provided!");
($config_file,$reads_file) = @ARGV;

print_ok();
##


## Step 2. Loading the information from the configuration file
my $line;
my $n_lines;
my @record;
my ($key,$primer1,$snp,$cg,$ref,$alt);
my @KEY;
my @PRIMER1;
my @SNP;
my @CG;
my @REF;
my @ALT;
my $i;
my $j;


print_mess("Stage 2.  Loading the configuration file values\n");

open(FILE,$config_file) or print_error("FILE $config_file can not be opened");
$n_lines=0;
while($line=<FILE>)
{
    print_mess($line);

    next if ($line =~ /\#/);

    @record = split(/\s+/,$line);

    $key = $record[0];
    $primer1 = $record[1];
    $snp = $record[2];
    $cg = $record[3];
    $ref = $record[4];
    $alt = $record[5];

    push(@KEY,$key);
    push(@PRIMER1,$primer1);

    # SNP: multiple bases separated by comma
    @record = split(",",$snp);
    for($i=0;$i<scalar(@record);$i++)
    {
	$SNP[$n_lines][$i] = $record[$i];
    }
    
    # CG: multiple bases separated by comma
    @record = split(",",$cg);
    for($i=0;$i<scalar(@record);$i++)
    {
	$CG[$n_lines][$i] = $record[$i];
    }

    # REF: multiple bases separated by comma
    @record = split(",",$ref);
    for($i=0;$i<scalar(@record);$i++)
    {
	$REF[$n_lines][$i] = $record[$i];
    }

    # ALT (CAS) SNP: single occurrence 
    push(@ALT,$alt);
    
    $n_lines++;
}
close(FILE);

print_mess("\n");
print_mess("KEYS: ",join("\t",@KEY),"\n");

print_mess("$n_lines elements from $config_file have been acquired");
print_ok();

# DEBUGGING LINES
#for($i=0; $i<$n_lines; $i++)
#{
#    print $i,"\n";
#    for($j=0; $j<scalar(@{$SNP[$i]}); $j++)
#    {
#	print $SNP[$i][$j],"\n";
#    }
#    print "\n";
#}

#for($i=0; $i<$n_lines; $i++)
#{
#    print $i,"\n";
#    for($j=0; $j<scalar(@{$CG[$i]}); $j++)
#    {
#	print $CG[$i][$j],"\n";
#    }
#    print "\n";
#}

#for($i=0; $i<$n_lines; $i++)
#{
#    print $i,"\n";
#    for($j=0; $j<scalar(@{$REF[$i]}); $j++)
#    {
#	print $REF[$i][$j],"\n";
#    }
#    print "\n";
#}


## Step 3. Identify the source of each read
my ($line1,$line2,$line3,$line4);
my $n_reads;
my $readID;
my ($pos,$posS,$posC);
my $start;
my ($SNP,$CG);
my $read;
my $species;
my ($foundref,$foundsnp,$foundcg);
my ($jSNP,$jCG);


print_mess("Stage 3.  Reading the FASTQ rawdata file\n");

open(FILE,$reads_file) or print_error("FILE $reads_file can not be opened");
$n_reads=0;

# first read
$line1 = <FILE>;
$line2 = <FILE>;
$line3 = <FILE>;
$line4 = <FILE>;
$species = $NOSPECIES;

# get the read ID
@record = split(/\s+/,$line1);
$readID = $record[0];
$read =  $line2;


# visit all candidates to identify which genes from the list of X are associated to this primer (copies with specific CGs)
$i=0;
while($i<$n_lines)
{
    # search the position of the primer inside the sequence
    $pos = index($line2,$PRIMER1[$i]);

    # positive hit in the search
    if ($pos!=$NOTFOUND)
    {
	# DEBUGGING LINES
	#print_mess("***** ".$readID."\t".$read."\n");
	#print_mess("***** gene ".$i." | pos ".$pos."\n");

	# SNP SEARCH
	$foundsnp = $FALSE;
	$j = 0;
	while($j<scalar(@{$SNP[$i]}) && !$foundsnp)
	{
	    # search the sequence upstream of the SNP
	    $posS = index($line2,$SNP[$i][$j]);
	    
	    if ($posS != $NOTFOUND)
	    {
		$foundsnp= $TRUE;
	    }
	    else
	    {
		$j++;
	    }
	}

	if ($foundsnp == $TRUE)
	{
	    #extract the SNP after the flanking sequence
	    $SNP = substr($line2,$posS+length($SNP[$i][$j]),1);
	    # save this SNP j
	    $jSNP = $j;
	    
	    #check for coincidence with ALT
	    if ($SNP eq $ALT[$i])
	    {
		$species = "CAS";
	    }
	    else
	    {
		# check the SNP is one of the available REF ones
		$foundref= $FALSE;
		$j=0;
		while($j<scalar(@{$REF[$i]}) && !$foundref)
		{
		    if ($SNP eq $REF[$i][$j])
		    {
			$foundref= $TRUE;
		    }
		    else
		    {
			$j++;
		    }
		}
		
		if ($foundref==$TRUE)
		{
		    $species = "MUS";
		}
	    }
	}
	else
	{
	    $SNP = $NOSEQ;
	}
	
	# CG SEARCH
	$foundcg = $FALSE;
	$j = 0;
	while($j<scalar(@{$CG[$i]}) && !$foundcg)
	{
	    # search the sequence upstream of the CG
	    $posC = index($line2,$CG[$i][$j]);

	    # if CG flanking sequence is right at the end, no CG can be extracted (results => not found) 
	    if ($posC != $NOTFOUND && ($posC+length($CG[$i][$j])<300))
	    {
		$foundcg= $TRUE;
	    }
	    else
	    {
		$j++;
	    }
	}
	
	if ($foundcg == $TRUE)
	{
	    #extract the CG after the flanking sequence
	    $CG = substr($line2,$posC+length($CG[$i][$j]),2);
	    # save this CG j
	    $jCG = $j;
	}
	else
	{
	    $CG = $NOSEQ;
	}
	
	# final output
	print join("\t",$readID,
		   $KEY[$i],
		   $PRIMER1[$i],
		   $pos+1,
		   ($foundsnp == $TRUE)? $SNP[$i][$jSNP]:$NOSEQ,
		   ($foundsnp == $TRUE)? $posS+1:$NOTFOUND,
		   ($foundsnp == $TRUE)? $posS+length($SNP[$i][$jSNP]):$NOTFOUND,
		   ($foundsnp == $TRUE)? $SNP:$NOSEQ,
		   ($foundcg == $TRUE)? $CG[$i][$jCG]:$NOSEQ,
		   ($foundcg == $TRUE)? $posC+1:$NOTFOUND,
		   ($foundcg == $TRUE)? $posC+length($CG[$i][$jCG]):$NOTFOUND,
		   ($foundcg == $TRUE)? $CG:$NOSEQ,
		   join(",",@{$REF[$i]}),
		   $ALT[$i],             
		   $species),"\n";
	
	# (debug) print the sequence of the read together with the output line
	if (exists($opt{v}))
	{
	    print $read,"\n";	
	}
    }
    else
    {
	# skip this candidate gene line
    }

    # next candidate gene
    $i++;
}


# first read processed
$n_reads++;

# analysis of the rest of reads
while($line1=<FILE>)
{
    $line2 = <FILE>;
    $line3 = <FILE>;
    $line4 = <FILE>;
    $species = $NOSPECIES;
    
    # get the read ID
    @record = split(/\s+/,$line1);
    $readID = $record[0];
    $read =  $line2;
    
    # visit all candidates to identify which genes from the list of X are associated to this primer (copies with specific CGs)
    $i=0;
    while($i<$n_lines)
    {
	# search the position of the primer inside the sequence
	$pos = index($line2,$PRIMER1[$i]);
	
	# positive hit in the search
	if ($pos!=$NOTFOUND)
	{
	    # DEBUGGING LINES
	    #print_mess("** $n_reads *** ".$readID."\n".$read);
	    #print_mess("***** FOUND gene ".$KEY[$i]." | pos ".$pos."\n");

	    # SNP SEARCH
	    $foundsnp = $FALSE;
	    $j = 0;
	    while($j<scalar(@{$SNP[$i]}) && !$foundsnp)
	    {
		# search the sequence upstream of the SNP
		$posS = index($line2,$SNP[$i][$j]);
		
		if ($posS != $NOTFOUND)
		{
		    $foundsnp= $TRUE;
		}
		else
		{
		    $j++;
		}
	    }
	    
	    if ($foundsnp == $TRUE)
	    {
		#extract the SNP after the flanking sequence
		$SNP = substr($line2,$posS+length($SNP[$i][$j]),1);
		# save this SNP j
		$jSNP = $j;
		
		#check for coincidence with ALT
		if ($SNP eq $ALT[$i])                                                                                                               
		{
		    $species = "CAS";
		}
		else
		{
		    # check the SNP is one of the available ones
		    $foundref = $FALSE;
		    $j=0;
		    while($j<scalar(@{$REF[$i]}) && !$foundref)
		    {
			if ($SNP eq $REF[$i][$j])
			{
			    $foundref= $TRUE;
			}
			else
			{
			    $j++;
			}
		    }
		    
		    if ($foundref==$TRUE)
		    {
			$species = "MUS";
		    }
		}  
	    }
	    else
	    {
		$SNP = $NOSEQ;
	    }
	    
	    # CG SEARCH
	    $foundcg = $FALSE;
	    $j = 0;
	    while($j<scalar(@{$CG[$i]}) && !$foundcg)
	    {
		# search the sequence upstream of the CG
		$posC = index($line2,$CG[$i][$j]);

		# if CG flanking sequence is right at the end, no CG can be extracted (results => not found) 
		if ($posC != $NOTFOUND && ($posC+length($CG[$i][$j])<300))
		{
		    $foundcg = $TRUE;
		}
		else
		{
		    $j++;
		}
	    }
	    
	    if ($foundcg == $TRUE)
	    {
		#extract the CG after the flanking sequence
		$CG = substr($line2,$posC+length($CG[$i][$j]),2);
		# save this CG j
		$jCG = $j;
	    }
	    else
	    {
		$CG = $NOSEQ;
	    }
	    
	    # final output
	    print join("\t",$readID,
		       $KEY[$i],
		       $PRIMER1[$i],
		       $pos+1,
		       ($foundsnp == $TRUE)? $SNP[$i][$jSNP]:$NOSEQ,
		       ($foundsnp == $TRUE)? $posS+1:$NOTFOUND,
		       ($foundsnp == $TRUE)? $posS+length($SNP[$i][$jSNP]):$NOTFOUND,
		       ($foundsnp == $TRUE)? $SNP:$NOSEQ,
		       ($foundcg == $TRUE)? $CG[$i][$jCG]:$NOSEQ,
		       ($foundcg == $TRUE)? $posC+1:$NOTFOUND,
		       ($foundcg == $TRUE)? $posC+length($CG[$i][$jCG]):$NOTFOUND,
		       ($foundcg == $TRUE)? $CG:$NOSEQ,
		       join(",",@{$REF[$i]}),
		       $ALT[$i],             
		       $species),"\n";
	    
	    # (debug) print the sequence of the read
	    if (exists($opt{v}))
	    {
		print $read,"\n";   
	    }
	}
	else
	{
	    # skip this candidate gene line
	}
	
	# next candidate gene
	$i++;
    }
   
    $n_reads++;
}
    
print_mess("$n_reads elements from $reads_file have been acquired");
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

