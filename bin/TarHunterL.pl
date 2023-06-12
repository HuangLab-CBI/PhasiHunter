#!/usr/bin/perl

#########################################################################
#  Minor edits by yuhanfei@g.harvard.edu
#########################################################################

use strict;
use warnings;
use 5.010;
use Getopt::Long qw/:config no_ignore_case/;
use Cwd 'abs_path';
use File::Temp qw/tempfile tempdir/;
use File::Basename;
use IO::File;



##########################  Global variables  ###########################


my $TarHunterL_version = "1.0";
my $fasta_version      = "fasta36";

# %target_IDs_to_annotations hash added to retain the full annotations from the input target fasta file
my ($temp_dir, $qmir_file, $targ_file, $output_file, $score_cutoff,
    $mimics, $tab_format, $help, $mir_seq_href, %target_IDs_to_annotations);
#    $mimics, $tab_format, $help, $mir_seq_href);

my $mispair_cutoff = 100;         #set loose cutoff
my $seed_mispair_cutoff = 100;    #set loose cutoff
my $mimics_arm_misp_cutoff = 0;   #default: no misp in eTM arms
my $threads = 1;                  #default: 1 CPU
my $fasta_output_ID = 0;

GetOptions (
    "q|qmir=s"         => \$qmir_file,
    "b|targ=s"         => \$targ_file,
    "M|total_misp=i"   => \$mispair_cutoff,
    "m|seed_misp=i"    => \$seed_mispair_cutoff,
    "f|score=f"        => \$score_cutoff,
    "I|mimics"         => \$mimics,
    "i|mimics_str=i"   => \$mimics_arm_misp_cutoff,
    "T|threads=i"      => \$threads,
    "t|tab"            => \$tab_format,
    "o|output=s"       => \$output_file,
    "h|help"           => \$help
);

$score_cutoff = 100 if $mimics;            #eTM: set loose score cutoff
$score_cutoff = 4 if (! $mimics and ! $score_cutoff); #default score: 4



############################  Main program  #############################


command_line();
check_arguments();
create_temp_folder();
check_tools();
read_mir();
# this subroutine was added to put in the full annotations from the input target fasta file into the output file
get_target_annotations();
target_search();

STDERR->say ("TarHunterL analysis is completed.\n");



#############################  Subroutines  #############################


sub command_line {
    STDERR->say ("\n", "=="x40);
    STDERR->say ("TarHunterL(version $TarHunterL_version)    Plant miRNA target searching tool");
    STDERR->say ("=="x40);

    my $time = scalar localtime;
    STDERR->say ("$time");

    my @cmdline = ();
    my $buffer;
    if (open my $cmd_fh, "<:raw", "/proc/$$/cmdline") {
        read($cmd_fh, $buffer, 8388608);
        close $cmd_fh;
        @cmdline = split /\0/s, $buffer;
    }
    STDERR->say ("Command Line: ", join(" ", @cmdline), "\n");
}



sub check_arguments {
    unless ($qmir_file) {
        STDERR->say ("FATAL: query miR file does not exist.");
        usage();
    }
    unless ($targ_file) {
        STDERR->say ("FATAL: target file is not provided.");
        usage();
    }
    unless ($output_file) {
        STDERR->say ("FATAL: output file is not provided.");
        usage();
    }
    if ($help) { usage(); }
    if (-e $output_file) {
        STDERR->say ("Output file already exists. TarHunterL exits.");
        exit;
    }
}



sub check_tools {
    my $FASTA_command = qx(which $fasta_version 2>>$temp_dir/msg.txt);
    unless ($FASTA_command) {
        STDERR->say ("FATAL: $fasta_version not found.");
        usage();
    }
}



sub create_temp_folder {
    #get TarHunterL directory path
    my $program_path = abs_path($0);
    my $dir = dirname($program_path);

    #create temporary folder
    system "mkdir $dir/temp" unless -e "$dir/temp";
    $temp_dir = tempdir('TarHunterL_job_XXXXXXXXX', DIR => "$dir/temp");
    STDERR->say ("Creating temporary folder ($temp_dir).");
}



sub usage {
my $usage = << "USAGE";

Usage:
    perl $0 -q <mir_file> -b <targ_file> -o <out_file> [Options]

Required arguments:
    -q (--qmir):         query miRNA file
    -b (--targ):         target file
    -o (--output):       output file

Options:

    -M (--total_misp):   max. total mispairs               [Default: off]
    -m (--seed_misp):    max. seed mispairs                [Default: off]
    -f (--score):        score cutoff                      [Default:  4 ]

    -I (--mimics):       eTM search                        [Default: off]
    -i (--mimics_str):   eTM stringency
                         (0: strict, 1: relaxed)           [Default:  0 ]

    -T (--threads):      FASTA threads                     [Default:  1 ]
    -t (--tab):          tabular format output             [Default: off]
    -h (--help):         help information

Dependencies:
    $fasta_version

USAGE

    print $usage;
    exit;
}



sub read_mir {
    STDERR->say ( "Reading miRNA file..." );
    $mir_seq_href = read_seq_file($qmir_file);
    if (! %{$mir_seq_href}) {
        STDERR->say ( "No miRNAs. TarHunterL exits." );
        exit;
    }
}



sub read_seq_file {
    my ($seqfile) = @_;
    my %seqhash;

    local $/ = '>';
    my $infile = IO::File->new($seqfile, 'r') or die $!;
    while (<$infile>) {
        s/\r//g;
        chomp;
        my ($id, $seq) = /(\S+)(?:.*?)\n(.+)/s or next;
        ($seq = uc $seq) =~ s/\n//g;
        $seqhash{$id} = $seq;
    }
    return \%seqhash;
}



sub target_search {
    if (defined $tab_format) {
        open my $out, '>>', $output_file or die $!;
        #print $out "targ_ID\ttarg_seq\tmir_ID\tmir_seq\tTotal_misp\tScore\t",
        #            "Seed_misp\tCleavage\tStart_pos\tSlice_pos\n";
        #print $out "miRNA_name\ttarget_name\tExpectation\tUPE\$\tmiRNA_start\tmiRNA_end\t",
        print $out "miRNA_name\ttarget_name\tExpectation\tUPE\$\tmiRNA_start\tmiRNA_end\t",
                    "target_start\ttarget_end\tmiRNA_seq\ttarget_seq\t",
                    "cleavage\tannotation\tmultiplicity\tseed_mismatch\ttotal_mismatch\tmiRNA_aligned_fragment\ttarget_aligned_fragment\tsymb\n";

        close $out;
    }

    for my $query_mir (sort keys %{$mir_seq_href}) {
        my $mir_tmpfile = IO::File->new("$temp_dir/mir.fa.tmp", 'w') or die $!;
        $mir_tmpfile->say (">$query_mir\n$mir_seq_href->{$query_mir}");
        $mir_tmpfile->close;

        run_FASTA($query_mir, "$temp_dir/mir.fa.tmp", $targ_file);
    }
}


# new subroutine for adding full annotation from the target FASTA into the output file from this script
sub get_target_annotations {
    open my $targfile_fh, "<", $targ_file, or die $!;
    
    while (my $line = <$targfile_fh>){
        chomp $line; 

        if ($line =~/^>(\S+)\s+(.+\S+.*)$/){
            # this captures lines with gene IDs, then space, then annotations
            my $seq_id = $1; 
            my $seq_annotation = $2;
            #print "$line\n$seq_id\n$seq_annotation\n";
            $target_IDs_to_annotations{$seq_id}= $seq_annotation;
        }elsif($line =~/^>(.*)$/){
            # this will capture lines with no annotations
            my $seq_id = $1;
            my $seq_annotation = "N/A";
            #print "$line\n$seq_id\n$seq_annotation\n";
            $target_IDs_to_annotations{$seq_id} = $seq_annotation;
        }else{
            # This will capture sequence lines 
        }
    }   
 
    close $targfile_fh;
}

sub run_FASTA {
    my ($mirName, $mirfile, $tarfile) = @_;

    my $mirAbbr = substr($mirName, 0, 5);
    my %fasta_hash;
    my $find_first = 0;
    my ($geneName, $rna_pair_beg_pos, $rna_pair_end_pos, $rna_start_pos, $rna_end_pos,
        $mir_end_pos, $mir_rc, $aln, $aln_copy, $tar, $target_annotation, $aln_line, $tar_line, $left_pos);
        #$mir_end_pos, $mir_rc, $aln, $aln_copy, $tar, $aln_line, $tar_line, $left_pos);
    my (@mir, @tar, @aln);

    STDERR->print ("Running FASTA for $mirName ...  ");

    open my $run_fasta, "$fasta_version -A -n -Q -i -U -T $threads -E 100000 $mirfile $tarfile 1 |" or die $!;
    #open my $output_fasta36, ">", "fasta36_output.txt" or die $!;
    while (<$run_fasta>) {
        #print $output_fasta36 $_;
        $fasta_output_ID ++;
        chomp;
        unless ($find_first) {
            next unless /^>>(\S+)/;
            #next unless /^>>(\S+)\s?.?\s(.*\S)\s+\(.+\)$/;
            #next unless /^>>(\S+)\s*(.*\S*)\s+\(.*\)$/;
            $geneName = $1;
            #$target_annotation=$2; 
            $find_first = 1;

        } else {
            if (/^>>(\S+)/) {
            #if (/^>>(\S+)\s?.?\s(.*\S)\s+\(.+\)$/) {
            #if (/^>>(\S+)\s*(.*\S*)\s+\(.*\)$/) {
                $geneName = $1;
                #$target_annotation =$2;
            } elsif (/Smith-Waterman.+\d+-\d+:(\d+)-\d+/) {
                $rna_pair_beg_pos = $1;

            } elsif (/^($mirAbbr\-\s+)(\S+)/) {
                $mir_rc = $2;

                $left_pos = length($1);
                $aln_line = <$run_fasta>;  #aln symbol
                if ($aln_line =~ /^(\s+)(\S.*\S)(\s*)$/) {
                    $aln = $2;
                    $aln_copy = $aln =~ s/://gr;
                }

                $tar_line = <$run_fasta>;
                $tar = substr($tar_line, $left_pos, length($mir_rc));
                $tar =~ s/\s/-/g;    #some targs have unpairings at 5' and 3' end

                #primary filtering
                next if (defined $aln_copy and length($aln_copy) > $mispair_cutoff);

                #FASTA_parser($mirName, $geneName, $mir_rc, $tar, $rna_pair_beg_pos, $target_annotation);
                FASTA_parser($mirName, $geneName, $mir_rc, $tar, $rna_pair_beg_pos);
            }else{
                print "Aberrant line encountered: $_\n" if /^>>/;
            }
        }
    }
    close $run_fasta;
    #close $output_fasta36;
    STDERR->say ("done.");
}



sub FASTA_parser {
    #my ($mirID, $geneID, $mir_rc, $tar, $rna_start_pos, $target_ann) = @_;
    my ($mirID, $geneID, $mir_rc, $tar, $rna_start_pos) = @_;
    my $tar_copy = $tar =~ s/-//gr;
    my $tar_len  = length($tar_copy);
    return -1 if $tar_len < 10; #remove some low-quality results

    (my $mir = reverse uc $mir_rc) =~ tr/AUGC/UACG/;
    my $mir_copy = $mir =~ s/-//gr;
    my $mir_len  = length($mir_copy);
    $tar = reverse uc $tar;
    my @mir = split //, $mir; #5'-3'
    my @tar = split //, $tar; #3'-5'
    my $aln_symb;
    my $cleavage = 'Yes';

    my ($count, $count2, $score, $g3g5, $total_misp, $seed_misp, $p8_misp, $middle_misp, $end_misp)
        = (0,0,0,0,0,0,0,0,0);

    my %pairing = (
    #class1:pairing, 2:G:U, 3:mismatch, 4:bulge miR strand, 5:bulge tar strand
    'AU' => 'class1',    'UA' => 'class1',    'GC' => 'class1',    'CG' => 'class1',
    'UG' => 'class2',    'GU' => 'class2',    'AG' => 'class3',    'AC' => 'class3',
    'UC' => 'class3',    'GA' => 'class3',    'CA' => 'class3',    'CU' => 'class3',
    'AA' => 'class3',    'UU' => 'class3',    'GG' => 'class3',    'CC' => 'class3',
    'A-' => 'class4',    'U-' => 'class4',    'G-' => 'class4',    'C-' => 'class4',
    '-A' => 'class5',    '-U' => 'class5',    '-G' => 'class5',    '-C' => 'class5',
    );

    open my $out, '>>', $output_file or die $!;
    for my $i (0..$#mir) {
        my $mir_nt = $mir[$i];
        my $tar_nt = $tar[$i];
        return -1 unless exists $pairing{"$mir_nt$tar_nt"};
        if ($pairing{"$mir_nt$tar_nt"} eq 'class1') {
            $count ++;
            $count2 ++ if $count <= 10;
            $aln_symb .= '|';

        } elsif ($pairing{"$mir_nt$tar_nt"} eq 'class2') {
            $count ++;
            $count2 ++ if $count <= 10;
            $total_misp ++;
            return -1 if $total_misp > $mispair_cutoff;
            $aln_symb .= ':';

            #Score
            if ($count >= 2 and $count <= 12) {
                $score += 1;
            } else {
                $score += 0.5;
            }
            return -1 if $score > $score_cutoff;
            $g3g5 = 1 if $count >= 3 and $count <= 5;

            #seed mispairs, cleavage
            if ($count >= 2 and $count <= 7) {
                $seed_misp ++;
                return -1 if $seed_misp > $seed_mispair_cutoff;
            } elsif ($count == 8 and defined $mimics) {
                $p8_misp ++;
            } elsif ($count >= 9 and $count <= 11) {
                $cleavage = 'No';
                $middle_misp ++ if defined $mimics;
            } elsif ($count >= 12 and $count <= 19 and defined $mimics) {
                $end_misp ++;

                #mimics arm (positions 2-8, 12-19) misp <= $mimics_arm_misp_cutoff
                return -1 if $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff;
            }

        } elsif ($pairing{"$mir_nt$tar_nt"} eq 'class3') {
            $count ++;
            $count2 ++ if $count <= 10;
            $total_misp ++;
            return -1 if $total_misp > $mispair_cutoff;
            $aln_symb .= '.';

            if ($count >= 2 and $count <= 12) {
                $score += 2;
            } else {
                $score += 1;
            }
            return -1 if $score > $score_cutoff;
            $g3g5 = 1 if $count >= 3 and $count <= 5;

            if ($count >= 2 and $count <= 7) {
                $seed_misp ++;
                return -1 if $seed_misp > $seed_mispair_cutoff;
            } elsif ($count == 8 and defined $mimics) {
                $p8_misp ++;
            } elsif ($count >= 9 and $count <= 11) {
                $cleavage = 'No';
                $middle_misp ++ if defined $mimics;
            } elsif ($count >= 12 and $count <= 19 and defined $mimics) {
                $end_misp ++;
                return -1 if $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff;
            }

        } elsif ($pairing{"$mir_nt$tar_nt"} eq 'class4') {
            $count ++;
            $total_misp ++;
            return -1 if $total_misp > $mispair_cutoff;
            $aln_symb .= 'o';

            if ($count >= 2 and $count <= 12) {
                $score += 2;
            } else {
                $score += 1;
            }
            return -1 if $score > $score_cutoff;
            $g3g5 = 1 if $count >= 3 and $count <= 5;

            if ($count >= 2 and $count <= 7) {
                $seed_misp ++;
                return -1 if $seed_misp > $seed_mispair_cutoff;
            } elsif ($count == 8 and defined $mimics) {
                $p8_misp ++;
            } elsif ($count >= 9 and $count <= 11) {
                $cleavage = 'No';
                $middle_misp ++ if defined $mimics;
            } elsif ($count >= 12 and $count <= 19 and defined $mimics) {
                $end_misp ++;
                return -1 if $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff;
            }

        } elsif ($pairing{"$mir_nt$tar_nt"} eq 'class5') {
            $count2 ++ if $count < 10;
            $total_misp ++;
            return -1 if $total_misp > $mispair_cutoff;
            $aln_symb .= 'o';

            if ($count >= 2 and $count <= 11) {
                $score += 2;
            } else {
                $score += 1;
            }
            return -1 if $score > $score_cutoff;
            $g3g5 = 1 if $count >= 3 and $count <= 4;

            if ($count >= 2 and $count <= 7) {
                $seed_misp ++;
                return -1 if $seed_misp > $seed_mispair_cutoff;
            } elsif ($count == 8 and defined $mimics) {
                $p8_misp ++;
            } elsif ($count >= 9 and $count <= 10) {
                $cleavage = 'No';
                $middle_misp ++ if defined $mimics;
            } elsif ($count >= 11 and $count <= 19 and defined $mimics) {# pos 12-19
                $end_misp ++;
                return -1 if $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff;
            }
        }
    }

    $score ++ if $g3g5;
    return -1 if $score > $score_cutoff;

    #mimics positions 9-11 misp <= 5
    if ( defined $mimics and
        ($cleavage eq 'Yes' or $middle_misp > 5 or
        $seed_misp + $p8_misp + $end_misp > $mimics_arm_misp_cutoff) ) {
        return -1;
    }

    $tar      = reverse $tar;  #change to 5'-3'
    my $mir_r = reverse $mir;  #change to 3'-5'
    $aln_symb = reverse $aln_symb;

    my $slice_pos = $rna_start_pos + $tar_len - $count2;
    my $upe="-1.0";
    my $miRNA_start="1";
    my $miRNA_end=$mir_len;
    my $rna_end_pos=$rna_start_pos + $tar_len - 1;
    #my $annotation=$geneAnnotation;
    my $annotation="-";
    my $multiplicity="1";
    #output detailed info
    if (defined $tab_format) {
        #print $out "$geneID\t$tar\t$mirID\t$mir\t$total_misp\t",
        #           "$score\t$seed_misp\t$cleavage\t$rna_start_pos\t$slice_pos\n";


        #$mir remove gap - $mir_copy
        #$tar remove gap - $tar_copy
        #$annotation from -b parameters

        # pull out the annotation from the hash with gene ID -> annotation mapping
        my $target_ann = $target_IDs_to_annotations{$geneID}; 
       
        print $out "$mirID\t$geneID\t$score\t$upe\t$miRNA_start\t$miRNA_end\t",
                    "$rna_start_pos\t$rna_end_pos\t$mir_copy\t$tar_copy\t",
                    "$cleavage\t$target_ann\t$multiplicity\t$seed_misp\t$total_misp\t$mir\t$tar\t$aln_symb\n";

    } else {
        print  $out "\#F$fasta_output_ID", "=="x40, "\n";
        printf $out "%-35s\t%s\n", $geneID, "5' $tar 3'";
        printf $out "%-35s\t%s\n", "", "   $aln_symb";
        printf $out "%-35s\t%s\n\n", $mirID, "3' $mir_r 5'";
        printf $out "%-35s\t%s\n", "Total mispairs:", $total_misp;
        printf $out "%-35s\t%s\n", "Score:", $score;
        printf $out "%-35s\t%s\n", "Seed mispairs:", $seed_misp;
        printf $out "%-35s\t%s\n", "Cleavage:", $cleavage;
        printf $out "%-35s\t%s\n", "Start position:", $rna_start_pos;
        printf $out "%-35s\t%s\n\n", "Slice position:", $slice_pos;
    }
    close $out;
}

__END__
