#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);

my $mock =0;


sub fileExists{
  my ($exeFile) = @_;

  if (!( -e $exeFile)) {
    die "Executable ".$exeFile." does not exist\n";
  }
}


my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $minq=30;
my $minl=35;
my $subsam=10000000;
my $length=10;


my $rohan    = $pathdir."/rohan";
my $bam2prof = $pathdir."/../bam2prof/bam2prof";

my $all=1;
my $double=0;
my $single=0;
my $map="0";


print STDERR "detecting ROHan...";
fileExists($rohan);
print STDERR  "..success\n";
$rohan = abs_path($rohan);

print STDERR  "detecting bam2prof...";
fileExists($bam2prof);
print STDERR  "..success\n";
$bam2prof = abs_path($bam2prof);

sub findcommand{
  my ($cmdtorun) = @_;
  print STDERR "trying to find executable ". $cmdtorun."..";
  my $command;

  if($mock != 1){
    my $cmdtodetect = "bash -i -c \"type ".$cmdtorun."\"";
    my $out = `$cmdtodetect`;

    if( $out =~ /^(\S+)\s+is\s+(\S+)$/){
      my @ar = split(" ",$out);
      $command=$ar[$#ar]."\n";
      print STDERR ".. found: ".$command."\n";
    }else{
      if( $out =~ /^(\S+)\s+is aliased to\s+(\S+)$/){
	my @ar = split(" ",$out);
	$command=$ar[$#ar]."\n";
	$command =~ s/^`//g;
	$command =~ s/'$//g;

	print STDERR ".. alias found: ".$command."\n";
      }else{
	print STDERR "cannot find executable ".$cmdtorun."\n";
	die;
      }
    }

  }else{
    print STDERR ".. running as mock: ".$cmdtorun."\n";
    $command = $cmdtorun;
  }
  return $command;

}

sub runcmd{
  my ($cmdtorun) = @_;

  print STDERR "running cmd ". $cmdtorun."\n";
  if($mock != 1){
    my @argstorun = ( "bash", "-i","-c", $cmdtorun );

    if(system(@argstorun) != 0){
      die "system  cmd $cmdtorun failed: $?"
    }else{
    }
  }

}

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
 This script is a wrapper to estimate the damage patterns with the
 potentially polymorphic positions masked.

\n\n usage:\t".$0." <options> [bam file] \n\n".

" Options:\n".
  "\t-o\t\t\t\t\tOutput prefix (default: [bam file])\n".
  "\t--mock\t\t\t\t\tDo nothing, just print the commands that will be run\n".

 "\n".
  " Damage estimate\n".
  " ===================\n".

 		"\t\t--minq\t\t\tRequire the base to have at least this quality to be considered (Default: ".$minq.")\n".
		"\t\t--minl\t\t\tRequire the read to have at least this mapping quality to be considered (Default: ".$minl.")\n".#TODO
		"\t\t--subsam\t\t\tSubsample this number of reads to consider (Default: ".$subsam.")\n".#TODO

		"\t\t--length\t[length]\t\tDo not consider bases beyond this length  (Default: ".$length." )\n".
		"\t\t--map\t[bed file]\t\tUse these mappability tracked in bed format  (Default: ".$length." )\n".

		  "\n".
                " \tSubstitutions reported: specify either one of the 3 possible options:

		  --all			Report all potential substitutions                      (Default: used)
		  --single		Use the deamination profile of a single strand library  (Default: not on/not used)
		  --double		Use the deamination profile of a double strand library  (Default: not on/not used)".

  "\n";
  exit;
}

my $help;
my $bamfile = $ARGV[$#ARGV];
my $outputprefix;

usage() if ( @ARGV < 1 or
             ! GetOptions('help|?' => \$help,'o=s' => \$outputprefix,'mock' => \$mock,'minq=i' => \$minq,'minl=i' => \$minl,'subsam=i' => \$subsam,'length=i' => \$length,'all' => \$all,'single' => \$single,'double' => \$double,'map=s' => \$map) or defined $help );


if( !(defined $outputprefix) ){
  $outputprefix = substr($bamfile,0,-4)."";
}else{
}
print STDERR  "Using output prefix ... ".$outputprefix."\n";

if($double != 1){
  $all=0;
}

if($single != 1){
  $all=0;
}

#samtools
my $cmdsamtools    = "samtools";#todo fix detection


#print STDERR  "typing to find samtools ";

$cmdsamtools=findcommand($cmdsamtools);





my $bam2profcmd1 = $cmdsamtools." -h view ".$bamfile." ";
if($map ne "0"){
  $bam2profcmd1 .= " -L ".$map." ";
}

$bam2profcmd1 .= "| head -n ".$subsam." | ".$cmdsamtools." -bS view /dev/stdin | ".$bam2prof." -minq ".$minq." -minl ".$minl." -length ".$length." ";
if($all){
  #do nothing
}
if($single){
  $bam2profcmd1 = $bam2profcmd1." -single ";
}
if($double){
  $bam2profcmd1 = $bam2profcmd1." -double ";
}
$bam2profcmd1 = $bam2profcmd1." -5p  ".$outputprefix."_1.5p.prof ";
$bam2profcmd1 = $bam2profcmd1." -3p  ".$outputprefix."_1.3p.prof ";

$bam2profcmd1 = $bam2profcmd1." /dev/stdin ";


runcmd($bam2profcmd1);



my $rohancommand = $cmdsamtools." -h view ".$bamfile." ";
die $outputprefix;


