#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $version = "0.0.23";

# version update
# 0.0.23
#   -added option to require a minimum level of evidence for sample level genotyping
#
# 0.0.22
#   -updated version to be consistent with dinumt package
#
# 0.0.16
#   -fixed program with too many supporting reads crashing GL calculation
#   -updated output to include more original INFO tag information
#
# 0.0.15
#   -added option to turn off EM iterations
#
# 0.0.14
#   -added MT mapping check for clipped sequences at potential breakpoints
#
# 0.0.13
#   -fixed error with likelihood calculation
#   -added support for by_chr_dir genotyping
#
# 0.0.12
#   -added per-read mapQ errors with likelihood calculation
#
# 0.0.11
#   -updated supporting read determination
#
# 0.0.10
#   -integrated RP and SR liklihoods, when available
#
# 0.0.9
#   -incorporated EM algorithm for implementing allele frequency priors
#
# 0.0.8
#   -updated support counting around breakpoints
#
# 0.0.7
#   -modified scoring to implement bayesian scoring of evidence
#   -incorporated insert length cutoff for alternative allele determination
#
# 0.0.6
#   -fixed bug where altSR was not being initialized correctly
#   -removed max_read_cov filter in report()
#
# 0.0.5
#   -added --chr option for by-chromosome analysis
#
# 0.0.4
#   -bug fixes
#
# 0.0.3
#   -modified scoring function to better assess clipping in breakpoint determination and scoring
#   -modified rp scoring to not consider reads clipped at breakpoint as reference
#
# 0.0.2
#   -modified input parameters to expect VCF format
#
# 0.0.1
#   -initial version

my %opts = ();

$opts{len_cluster_include} = 600;
$opts{len_cluster_link}    = 800;
$opts{min_reads_cluster}   = 1;
$opts{min_clipped_seq}     = 5;
$opts{clipped_flank}       = 50;
$opts{max_num_clipped}     = 5;
$opts{include_mask}        = 0;
$opts{min_evidence}        = 3;
$opts{min_map_qual}        = 10;
$opts{filter_qual}         = 13;
$opts{filter_depth}        = 5;
$opts{max_read_cov}        = 200;
$opts{use_priors}          = 0;
$opts{mask_filename}       = "/home2/remills/projects/numts/numtS.bed";
$opts{info_filename}       = "/scratch/remills_flux/remills/sampleInfo.txt";
$opts{reference}           = "/nfs/remills-scratch/reference/hs37d5/hs37d5.fa";
$opts{mt_filename}         = "/nfs/remills-scratch/reference/GRCh37/Sequence/Chromosomes/MT.fa";
$opts{samtools}            = "/home2/remills/bin/samtools";
$opts{exonerate}           = "/home2/remills/apps/exonerate-2.2.0/bin/exonerate";
$opts{dir_tmp}             = "/tmp";
$opts{prefix}              = "numt";
$opts{len_mt}              = 16596;                                                                #eventually should be read in by BAM header
$opts{ploidy}              = 2;

my $optResult = GetOptions(
    "input_filename=s"      => \$opts{input_filename},
    "output_filename=s"     => \$opts{output_filename},
    "mask_filename=s"       => \$opts{mask_filename},
    "info_filename=s"       => \$opts{info_filename},
    "mt_filename=s"         => \$opts{mt_filename},
    "dir_tmp=s"             => \$opts{dir_tmp},
    "chr=s"                 => \$opts{chr},
    "include_mask"          => \$opts{include_mask},
    "len_cluster_include=i" => \$opts{len_cluster_include},
    "len_cluster_link=i"    => \$opts{len_cluster_link},
    "min_reads_cluster=i"   => \$opts{min_reads_cluster},
    "min_evidence=i"        => \$opts{min_evidence},
    "min_clipped_seq=i"     => \$opts{min_clipped_seq},
    "max_num_clipped=i"     => \$opts{max_num_clipped},
    "min_map_qual=i"        => \$opts{min_map_qual},
    "max_read_cov=i"        => \$opts{max_read_cov},
    "read_groups"           => \$opts{read_groups},
    "breakpoint"            => \$opts{breakpoint},
    "reference=s"           => \$opts{reference},
    "samtools=s"            => \$opts{samtools},
    "exonerate=s"           => \$opts{exonerate},
    "use_priors"            => \$opts{use_priors},
    "by_chr_dir"            => \$opts{by_chr_dir},
    "prefix=s"              => \$opts{prefix},
    "ucsc"                  => \$opts{ucsc},
    "help"                  => \$opts{help},
    "verbose"               => \$opts{verbose}
);

checkOptions( $optResult, \%opts, $version );

my $seq_num  = 0;
my %seq_hash = ();

my %sorted_hash = ();

my $i            = 1;
my %sample_hash  = ();
my %infile_hash  = ();
my %group_hash   = ();
my %outfile_hash = ();
my %mask_hash    = ();
my %data_hash    = ();

getInput( \%infile_hash, \%mask_hash, \%sample_hash );
getData( \%infile_hash, \%mask_hash, \%sample_hash, \%data_hash );
if ($opts{breakpoint}) {
    assessBreaks( \%infile_hash, \%sample_hash, \%data_hash );
}
refineData( \%infile_hash, \%mask_hash, \%sample_hash, \%data_hash );
scoreData( \%data_hash, \%sample_hash );
report( \%infile_hash, \%data_hash );

################################################################################################################
sub mappedClippedSeq {
    my ($clipSeq) = @_;

    if ( $clipSeq =~ /N/ ) { return 0; }

    my $random    = int( rand(100000) );
    my @timeData  = localtime(time);
    my $timestamp = join( '', @timeData );
    my $fileClp   = "$opts{dir_tmp}/clp$timestamp$random.fa";
    open( CLP, ">$fileClp" );
    print CLP ">clp\n$clipSeq\n";
    close CLP;

    my $results = `$opts{exonerate} --model affine:local --exhaustive yes --percent 70 --showvulgar no --showalignment no --showcigar yes $fileClp $opts{mt_filename} 2> /dev/null`;
    my ( $query_id, $query_start, $query_end, $query_strand, $target_id, $target_start, $target_end, $target_strand, $score, $cigar ) = $results =~ /cigar\:\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|-)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|-)\s+(\d+)\s+(.*?)\n/;
    `rm $fileClp`;
    if ($cigar) { return 1; }
    else { return 0; } 
}

sub assessBreaks {
    my ( $infile_hash, $sample_hash, $data_hash ) = @_;
    print "Entering assessBreaks()\n" if $opts{verbose};
    foreach my $var ( keys %$data_hash ) {
        print "-variant $var\n" if $opts{verbose};
        my %clippedPos = ();
        my %clippedSeq = ();
        my $sumClipped = 0;
        my $winStart   = 1e10;
        my $winEnd     = 0;
        my $maxFirst   = 0;
        my $maxSecond  = 0;

        foreach my $sample ( keys %{ $data_hash{$var} } ) {
            foreach my $cPos ( keys %{ $data_hash{$var}{$sample}{clipPos} } ) {
                next if $cPos == -1;
                if ( defined( $data_hash{$var}{$sample}{clipPos}{$cPos}{1} ) ) {
                    if ( mappedClippedSeq( $data_hash{$var}{$sample}{clipSeq}{$cPos} ) ) {
                        $clippedPos{$cPos} += $data_hash{$var}{$sample}{clipPos}{$cPos}{1};
                        $sumClipped += $data_hash{$var}{$sample}{clipPos}{$cPos}{1};

                        #if (!defined($clippedSeq{$cPos} || length($data_hash{$var}{$sample}{clipSeq}{$cPos}) > $clippedSeq{$cPos})) { $clippedSeq{$cPos} = $data_hash{$var}{$sample}{clipSeq}{$cPos}; }
                    }
                }
            }
        }
        my $winLen = scalar keys %clippedPos;
        my @sorted = sort { $clippedPos{$b} <=> $clippedPos{$a} } keys %clippedPos;

        if ( $winLen == 0 ) { next; }
        elsif ( $winLen == 1 ) {
            $infile_hash{$var}{leftBkpt}  = $sorted[0];
            $infile_hash{$var}{rightBkpt} = $sorted[0] + 1;
        }
        elsif ( $winLen == 2 ) {
            $infile_hash{$var}{leftBkpt}  = $sorted[0];
            $infile_hash{$var}{rightBkpt} = $sorted[0] + 1;
            if ( $sorted[1] > $sorted[0] ) {
                $infile_hash{$var}{rightBkpt} = $sorted[1];
            }
            else {
                $infile_hash{$var}{leftBkpt} = $sorted[1];
                $infile_hash{$var}{rightBkpt}--;
            }
        }
        else {
            my @sorted      = sort { $clippedPos{$b} <=> $clippedPos{$a} } keys %clippedPos;
            my $meanClipped = $sumClipped / $winLen;
            my $sumSquares  = 0;
            print "\tmean clipped: $meanClipped\n" if $opts{verbose};

            #for ( my $p = $winStart ; $p <= $winEnd ; $p++ ) {
            #    if ( !defined( $clippedPos{$p} ) ) { $clippedPos{$p} = 0; }
            #    $sumSquares += ( $clippedPos{$p} - $meanClipped )**2;
            #}
            foreach my $cPos ( keys %clippedPos ) {
                $sumSquares += ( $clippedPos{$cPos} - $meanClipped )**2;
            }
            my $sdClipped = sqrt( $sumSquares / ( $winLen - 1 ) );
            print "\tsd clipped: $sdClipped\n"                 if $opts{verbose};
            print "\t1st clipped: $clippedPos{ $sorted[0] }\n" if defined( $sorted[0] ) && $opts{verbose};
            print "\t2nd clipped: $clippedPos{ $sorted[1] }\n" if defined( $sorted[1] ) && $opts{verbose};

            if ( defined( $sorted[0] ) && $clippedPos{ $sorted[0] } > $meanClipped + 1 * $sdClipped ) {
                print "\t*updating to 1st\n" if $opts{verbose};
                $infile_hash{$var}{leftBkpt}  = $sorted[0];
                $infile_hash{$var}{rightBkpt} = $sorted[0] + 1;
            }
            if ( defined( $sorted[1] ) && $clippedPos{ $sorted[1] } > $meanClipped + 1 * $sdClipped ) {
                print "\t*updating to 1st and 2nd\n" if $opts{verbose};
                if ( $sorted[1] > $sorted[0] ) {
                    $infile_hash{$var}{rightBkpt} = $sorted[1];
                }
                else {
                    $infile_hash{$var}{leftBkpt} = $sorted[1];
                    $infile_hash{$var}{rightBkpt}--;
                }
            }
        }
        print "-variant $var completed\n" if $opts{verbose};
    }
    print "Exiting assessBreaks()\n\n" if $opts{verbose};
}

sub scoreData {
    my ( $data_hash, $sample_hash ) = @_;
    print "Entering scoreData()\n" if $opts{verbose};

    foreach my $var ( keys %$data_hash ) {
        print "Variant: $var\n" if $opts{verbose};

        my %priors   = ();
        my %genoFreq = ();
        foreach my $geno ( 0 .. $opts{ploidy} ) {

            #start with uniform priors
            $priors{$geno} = 1 / ( $opts{ploidy} + 1 );
            $genoFreq{$geno}{old} = 0;
            $genoFreq{$geno}{new} = 0;
        }

        my $numIteration = 0;
        my $sumGenoFreq  = 0;
        while ( $numIteration < 10 ) {
            $numIteration++;
            foreach my $sample ( keys %{$sample_hash} ) {

                print "\tSample: $sample\n" if $opts{verbose};
                
                #zero out alternative supporting evidence if below threshold
                if ($$data_hash{$var}{$sample}{numAltRP} + $$data_hash{$var}{$sample}{numAltSR} < $opts{min_evidence}) {
                    $$data_hash{$var}{$sample}{numAltRP} = 0;
                    $$data_hash{$var}{$sample}{numAltSR} = 0;
                    $$data_hash{$var}{$sample}{qualAltRP} = ();
                    $$data_hash{$var}{$sample}{qualAltSR} = ();
                }
                my $numRefRP = $$data_hash{$var}{$sample}{numRefRP};
                my $numAltRP = $$data_hash{$var}{$sample}{numAltRP};
                my $numRefSR = $$data_hash{$var}{$sample}{numRefSR};
                my $numAltSR = $$data_hash{$var}{$sample}{numAltSR};

                foreach my $geno ( 0 .. $opts{ploidy} ) {
                    $$data_hash{$var}{$sample}{pl}{$geno}  = 0;
                    $$data_hash{$var}{$sample}{gl}{$geno}  = 0;
                    $$data_hash{$var}{$sample}{gl0}{$geno} = 0;
                }
                $$data_hash{$var}{$sample}{gq} = 0;
                $$data_hash{$var}{$sample}{gt} = "./.";
                $$data_hash{$var}{$sample}{ft} = "LowQual";

                if ( $numAltRP + $numRefRP + $numAltSR + $numRefSR == 0 ) { next; }
                if ( $$data_hash{$var}{$sample}{avgQ} <= 0 || $$data_hash{$var}{$sample}{avgQ} > 1 ) { $$data_hash{$var}{$sample}{avgQ} = 0.999999; }
                $$data_hash{$var}{$sample}{avgQ} = 0.00001;
                foreach my $g ( 0 .. $opts{ploidy} ) {
                    my $geno = $opts{ploidy} - $g;    #need to reverse as calculation is reference allele based
                    if ( $numAltRP + $numRefRP > 0 &&  1 / $opts{ploidy} ** ($numAltRP + $numRefRP) > 0) {
                        $$data_hash{$var}{$sample}{gl0}{$geno} += calcGl( $opts{ploidy}, $g, $numAltRP + $numRefRP, $numRefRP, $$data_hash{$var}{$sample}{qualRefRP}, $$data_hash{$var}{$sample}{qualAltRP} );
                    }
                    if ( $numAltSR + $numRefSR > 0 && $opts{breakpoint} &&  1 / $opts{ploidy} ** ($numAltSR + $numRefSR) > 0) {
                        $$data_hash{$var}{$sample}{gl0}{$geno} += calcGl( $opts{ploidy}, $g, $numAltSR + $numRefSR, $numRefSR, $$data_hash{$var}{$sample}{qualRefSR}, $$data_hash{$var}{$sample}{qualAltSR} );
                    }
                    print "\tgl0 returned from calcGL() for geno $geno: $$data_hash{$var}{$sample}{gl0}{$geno}\n" if $opts{verbose};
                    if ( $$data_hash{$var}{$sample}{gl0}{$geno} < -255 ) { $$data_hash{$var}{$sample}{gl0}{$geno} = -255; }    #capped
                    $$data_hash{$var}{$sample}{gl}{$geno} = log10( $priors{$geno} ) + $$data_hash{$var}{$sample}{gl0}{$geno};  #adjust with population inferred prior
                }

                my @sortedGeno = sort { $$data_hash{$var}{$sample}{gl}{$b} <=> $$data_hash{$var}{$sample}{gl}{$a} } keys %{ $$data_hash{$var}{$sample}{gl} };

                #calculate PL from GL
                foreach my $geno ( 0 .. $opts{ploidy} ) {
                    print "\tcalculating pl from geno ($$data_hash{$var}{$sample}{gl}{$geno})\n" if $opts{verbose};
                    $$data_hash{$var}{$sample}{pl}{$geno} = int( -10 * $$data_hash{$var}{$sample}{gl}{$geno} );
                    if ( $$data_hash{$var}{$sample}{pl}{$geno} > 255 ) { $$data_hash{$var}{$sample}{pl}{$geno} = 255; }
                    print "\t...$$data_hash{$var}{$sample}{pl}{$geno}\n" if $opts{verbose};
                }

                #normalize PL to most likely genotype
                foreach my $geno ( 0 .. $opts{ploidy} ) {
                    $$data_hash{$var}{$sample}{pl}{$geno} -= $$data_hash{$var}{$sample}{pl}{ $sortedGeno[0] };
                }

                #determine genotype quality
                $$data_hash{$var}{$sample}{gq} = int( 10 * ( $$data_hash{$var}{$sample}{gl}{ $sortedGeno[0] } - $$data_hash{$var}{$sample}{gl}{ $sortedGeno[1] } ) );
                print "\t...$$data_hash{$var}{$sample}{gq}\n" if $opts{verbose};

                my $gt = "0/0";
                if    ( $sortedGeno[0] == 1 ) { $gt = "0/1"; }
                elsif ( $sortedGeno[0] == 2 ) { $gt = "1/1"; }
                $$data_hash{$var}{$sample}{gt} = $gt;
                if   ( $numAltRP + $numRefRP + $numAltSR + $numRefSR < $opts{filter_depth} || $$data_hash{$var}{$sample}{gq} < $opts{filter_qual} ) { $$data_hash{$var}{$sample}{ft} = "LowQual"; }
                else                                                                                                                                { $$data_hash{$var}{$sample}{ft} = "PASS"; }

                $genoFreq{ $sortedGeno[0] }{new}++;
                $sumGenoFreq++;
            }

            my $isConverged = 1;
            if ( $sumGenoFreq == 0 ) { warn "no genotypes passing filters\n"; last; }
            foreach my $geno ( 0 .. $opts{ploidy} ) {
                if ( $genoFreq{$geno}{new} == 0 ) {
                    $priors{$geno} = 1 / $sumGenoFreq;
                }
                else {
                    $priors{$geno} = $genoFreq{$geno}{new} / $sumGenoFreq;
                }
                if ( $genoFreq{$geno}{old} != $genoFreq{$geno}{new} ) { $isConverged = 0; }
                $genoFreq{$geno}{old} = $genoFreq{$geno}{new};
                $genoFreq{$geno}{new} = 0;
            }
            if ($isConverged || !$opts{use_priors}) { last; }
        }
    }
    print "Exiting scoreData()\n\n" if $opts{verbose};
}

sub getDate {
    my ( $second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings ) = localtime();
    my $year = 1900 + $yearOffset;
    $month++;
    my $fmonth = sprintf( "%.2d", $month );
    my $fday   = sprintf( "%.2d", $dayOfMonth );
    return "$year$fmonth$fday";
}

sub report {
    my ( $infile_hash, $data_hash ) = @_;
    print "Entering report()\n" if $opts{verbose};

    #open output file
    if ( defined( $opts{output_filename} ) ) {
        open( foutname1, ">$opts{output_filename}" ) or die("error opening file $opts{output_filename}\n");
    }
    else {
        open( foutname1, ">&", \*STDOUT ) or die;
    }

    my $filedate = getDate();
    print foutname1 <<HEADER;
##fileformat=VCFv4.1
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS:MT,Description="Nuclear Mitochondrial Insertion">
##FILTER=<ID=LowQual,Description="No PASS calls in any sample after merging">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNL,Number=.,Type=Float,Description="Copy number likelihoods">
##FORMAT=<ID=CNL0,Number=.,Type=Float,Description="Copy number likelihoods with no frequency prior">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=GL0,Number=G,Type=Float,Description="Genotype likelihoods with no frequency prior">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=MSTART,Number=1,Type=Flag,Description="Mitochondrial start coordinate of inserted sequence">
##INFO=<ID=MEND,Number=1,Type=Flag,Description="Mitochondrial end coordinate of inserted sequence">
##INFO=<ID=MLEN,Number=1,Type=Flag,Description="Estimated length of mitochondrial insert">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">
##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample(s) in which site was originally discovered">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr15_KI270905v1_alt,length=5161414>
##contig=<ID=chr6_GL000256v2_alt,length=4929269>
##contig=<ID=chr6_GL000254v2_alt,length=4827813>
##contig=<ID=chr6_GL000251v2_alt,length=4795265>
##contig=<ID=chr6_GL000253v2_alt,length=4677643>
##contig=<ID=chr6_GL000250v2_alt,length=4672374>
##contig=<ID=chr6_GL000255v2_alt,length=4606388>
##contig=<ID=chr6_GL000252v2_alt,length=4604811>
##contig=<ID=chr17_KI270857v1_alt,length=2877074>
##contig=<ID=chr16_KI270853v1_alt,length=2659700>
##contig=<ID=chr16_KI270728v1_random,length=1872759>
##contig=<ID=chr17_GL000258v2_alt,length=1821992>
##contig=<ID=chr5_GL339449v2_alt,length=1612928>
##contig=<ID=chr14_KI270847v1_alt,length=1511111>
##contig=<ID=chr17_KI270908v1_alt,length=1423190>
##contig=<ID=chr14_KI270846v1_alt,length=1351393>
##contig=<ID=chr5_KI270897v1_alt,length=1144418>
##contig=<ID=chr7_KI270803v1_alt,length=1111570>
##contig=<ID=chr19_GL949749v2_alt,length=1091841>
##contig=<ID=chr19_KI270938v1_alt,length=1066800>
##contig=<ID=chr19_GL949750v2_alt,length=1066390>
##contig=<ID=chr19_GL949748v2_alt,length=1064304>
##contig=<ID=chr19_GL949751v2_alt,length=1002683>
##contig=<ID=chr19_GL949746v1_alt,length=987716>
##contig=<ID=chr19_GL949752v1_alt,length=987100>
##contig=<ID=chr8_KI270821v1_alt,length=985506>
##contig=<ID=chr1_KI270763v1_alt,length=911658>
##contig=<ID=chr6_KI270801v1_alt,length=870480>
##contig=<ID=chr19_GL949753v2_alt,length=796479>
##contig=<ID=chr19_GL949747v2_alt,length=729520>
##contig=<ID=chr8_KI270822v1_alt,length=624492>
##contig=<ID=chr4_GL000257v2_alt,length=586476>
##contig=<ID=chr12_KI270904v1_alt,length=572349>
##contig=<ID=chr4_KI270925v1_alt,length=555799>
##contig=<ID=chr15_KI270852v1_alt,length=478999>
##contig=<ID=chr15_KI270727v1_random,length=448248>
##contig=<ID=chr9_KI270823v1_alt,length=439082>
##contig=<ID=chr15_KI270850v1_alt,length=430880>
##contig=<ID=chr1_KI270759v1_alt,length=425601>
##contig=<ID=chr12_GL877876v1_alt,length=408271>
##contig=<ID=chrUn_KI270442v1,length=392061>
##contig=<ID=chr17_KI270862v1_alt,length=391357>
##contig=<ID=chr15_GL383555v2_alt,length=388773>
##contig=<ID=chr19_GL383573v1_alt,length=385657>
##contig=<ID=chr4_KI270896v1_alt,length=378547>
##contig=<ID=chr4_GL383528v1_alt,length=376187>
##contig=<ID=chr17_GL383563v3_alt,length=375691>
##contig=<ID=chr8_KI270810v1_alt,length=374415>
##contig=<ID=chr1_GL383520v2_alt,length=366580>
##contig=<ID=chr1_KI270762v1_alt,length=354444>
##contig=<ID=chr15_KI270848v1_alt,length=327382>
##contig=<ID=chr17_KI270909v1_alt,length=325800>
##contig=<ID=chr14_KI270844v1_alt,length=322166>
##contig=<ID=chr8_KI270900v1_alt,length=318687>
##contig=<ID=chr10_GL383546v1_alt,length=309802>
##contig=<ID=chr13_KI270838v1_alt,length=306913>
##contig=<ID=chr8_KI270816v1_alt,length=305841>
##contig=<ID=chr22_KI270879v1_alt,length=304135>
##contig=<ID=chr8_KI270813v1_alt,length=300230>
##contig=<ID=chr11_KI270831v1_alt,length=296895>
##contig=<ID=chr15_GL383554v1_alt,length=296527>
##contig=<ID=chr8_KI270811v1_alt,length=292436>
##contig=<ID=chr18_GL383567v1_alt,length=289831>
##contig=<ID=chrX_KI270880v1_alt,length=284869>
##contig=<ID=chr8_KI270812v1_alt,length=282736>
##contig=<ID=chr19_KI270921v1_alt,length=282224>
##contig=<ID=chr17_KI270729v1_random,length=280839>
##contig=<ID=chr17_JH159146v1_alt,length=278131>
##contig=<ID=chrX_KI270913v1_alt,length=274009>
##contig=<ID=chr6_KI270798v1_alt,length=271782>
##contig=<ID=chr7_KI270808v1_alt,length=271455>
##contig=<ID=chr22_KI270876v1_alt,length=263666>
##contig=<ID=chr15_KI270851v1_alt,length=263054>
##contig=<ID=chr22_KI270875v1_alt,length=259914>
##contig=<ID=chr1_KI270766v1_alt,length=256271>
##contig=<ID=chr19_KI270882v1_alt,length=248807>
##contig=<ID=chr3_KI270778v1_alt,length=248252>
##contig=<ID=chr15_KI270849v1_alt,length=244917>
##contig=<ID=chr4_KI270786v1_alt,length=244096>
##contig=<ID=chr12_KI270835v1_alt,length=238139>
##contig=<ID=chr17_KI270858v1_alt,length=235827>
##contig=<ID=chr19_KI270867v1_alt,length=233762>
##contig=<ID=chr16_KI270855v1_alt,length=232857>
##contig=<ID=chr8_KI270926v1_alt,length=229282>
##contig=<ID=chr5_GL949742v1_alt,length=226852>
##contig=<ID=chr3_KI270780v1_alt,length=224108>
##contig=<ID=chr17_GL383565v1_alt,length=223995>
##contig=<ID=chr2_KI270774v1_alt,length=223625>
##contig=<ID=chr4_KI270790v1_alt,length=220246>
##contig=<ID=chr11_KI270927v1_alt,length=218612>
##contig=<ID=chr19_KI270932v1_alt,length=215732>
##contig=<ID=chr11_KI270903v1_alt,length=214625>
##contig=<ID=chr2_KI270894v1_alt,length=214158>
##contig=<ID=chr14_GL000225v1_random,length=211173>
##contig=<ID=chrUn_KI270743v1,length=210658>
##contig=<ID=chr11_KI270832v1_alt,length=210133>
##contig=<ID=chr7_KI270805v1_alt,length=209988>
##contig=<ID=chr4_GL000008v2_random,length=209709>
##contig=<ID=chr7_KI270809v1_alt,length=209586>
##contig=<ID=chr19_KI270887v1_alt,length=209512>
##contig=<ID=chr4_KI270789v1_alt,length=205944>
##contig=<ID=chr3_KI270779v1_alt,length=205312>
##contig=<ID=chr19_KI270914v1_alt,length=205194>
##contig=<ID=chr19_KI270886v1_alt,length=204239>
##contig=<ID=chr11_KI270829v1_alt,length=204059>
##contig=<ID=chr14_GL000009v2_random,length=201709>
##contig=<ID=chr21_GL383579v2_alt,length=201197>
##contig=<ID=chr11_JH159136v1_alt,length=200998>
##contig=<ID=chr19_KI270930v1_alt,length=200773>
##contig=<ID=chrUn_KI270747v1,length=198735>
##contig=<ID=chr18_GL383571v1_alt,length=198278>
##contig=<ID=chr19_KI270920v1_alt,length=198005>
##contig=<ID=chr6_KI270797v1_alt,length=197536>
##contig=<ID=chr3_KI270935v1_alt,length=197351>
##contig=<ID=chr17_KI270861v1_alt,length=196688>
##contig=<ID=chr15_KI270906v1_alt,length=196384>
##contig=<ID=chr5_KI270791v1_alt,length=195710>
##contig=<ID=chr14_KI270722v1_random,length=194050>
##contig=<ID=chr16_GL383556v1_alt,length=192462>
##contig=<ID=chr13_KI270840v1_alt,length=191684>
##contig=<ID=chr14_GL000194v1_random,length=191469>
##contig=<ID=chr11_JH159137v1_alt,length=191409>
##contig=<ID=chr19_KI270917v1_alt,length=190932>
##contig=<ID=chr7_KI270899v1_alt,length=190869>
##contig=<ID=chr19_KI270923v1_alt,length=189352>
##contig=<ID=chr10_KI270825v1_alt,length=188315>
##contig=<ID=chr19_GL383576v1_alt,length=188024>
##contig=<ID=chr19_KI270922v1_alt,length=187935>
##contig=<ID=chrUn_KI270742v1,length=186739>
##contig=<ID=chr22_KI270878v1_alt,length=186262>
##contig=<ID=chr19_KI270929v1_alt,length=186203>
##contig=<ID=chr11_KI270826v1_alt,length=186169>
##contig=<ID=chr6_KB021644v2_alt,length=185823>
##contig=<ID=chr17_GL000205v2_random,length=185591>
##contig=<ID=chr1_KI270765v1_alt,length=185285>
##contig=<ID=chr19_KI270916v1_alt,length=184516>
##contig=<ID=chr19_KI270890v1_alt,length=184499>
##contig=<ID=chr3_KI270784v1_alt,length=184404>
##contig=<ID=chr12_GL383551v1_alt,length=184319>
##contig=<ID=chr20_KI270870v1_alt,length=183433>
##contig=<ID=chrUn_GL000195v1,length=182896>
##contig=<ID=chr1_GL383518v1_alt,length=182439>
##contig=<ID=chr22_KI270736v1_random,length=181920>
##contig=<ID=chr10_KI270824v1_alt,length=181496>
##contig=<ID=chr14_KI270845v1_alt,length=180703>
##contig=<ID=chr3_GL383526v1_alt,length=180671>
##contig=<ID=chr13_KI270839v1_alt,length=180306>
##contig=<ID=chr22_KI270733v1_random,length=179772>
##contig=<ID=chrUn_GL000224v1,length=179693>
##contig=<ID=chr10_GL383545v1_alt,length=179254>
##contig=<ID=chrUn_GL000219v1,length=179198>
##contig=<ID=chr5_KI270792v1_alt,length=179043>
##contig=<ID=chr17_KI270860v1_alt,length=178921>
##contig=<ID=chr19_GL000209v2_alt,length=177381>
##contig=<ID=chr11_KI270830v1_alt,length=177092>
##contig=<ID=chr9_KI270719v1_random,length=176845>
##contig=<ID=chrUn_GL000216v2,length=176608>
##contig=<ID=chr22_KI270928v1_alt,length=176103>
##contig=<ID=chr1_KI270712v1_random,length=176043>
##contig=<ID=chr6_KI270800v1_alt,length=175808>
##contig=<ID=chr1_KI270706v1_random,length=175055>
##contig=<ID=chr2_KI270776v1_alt,length=174166>
##contig=<ID=chr18_KI270912v1_alt,length=174061>
##contig=<ID=chr3_KI270777v1_alt,length=173649>
##contig=<ID=chr5_GL383531v1_alt,length=173459>
##contig=<ID=chr3_JH636055v2_alt,length=173151>
##contig=<ID=chr14_KI270725v1_random,length=172810>
##contig=<ID=chr5_KI270796v1_alt,length=172708>
##contig=<ID=chrEBV,length=171823>
##contig=<ID=chr9_GL383541v1_alt,length=171286>
##contig=<ID=chr19_KI270885v1_alt,length=171027>
##contig=<ID=chr19_KI270919v1_alt,length=170701>
##contig=<ID=chr19_KI270889v1_alt,length=170698>
##contig=<ID=chr19_KI270891v1_alt,length=170680>
##contig=<ID=chr19_KI270915v1_alt,length=170665>
##contig=<ID=chr19_KI270933v1_alt,length=170537>
##contig=<ID=chr19_KI270883v1_alt,length=170399>
##contig=<ID=chr19_GL383575v2_alt,length=170222>
##contig=<ID=chr19_KI270931v1_alt,length=170148>
##contig=<ID=chr12_GL383550v2_alt,length=169178>
##contig=<ID=chr13_KI270841v1_alt,length=169134>
##contig=<ID=chrUn_KI270744v1,length=168472>
##contig=<ID=chr18_KI270863v1_alt,length=167999>
##contig=<ID=chr18_GL383569v1_alt,length=167950>
##contig=<ID=chr12_GL877875v1_alt,length=167313>
##contig=<ID=chr21_KI270874v1_alt,length=166743>
##contig=<ID=chr3_KI270924v1_alt,length=166540>
##contig=<ID=chr1_KI270761v1_alt,length=165834>
##contig=<ID=chr3_KI270937v1_alt,length=165607>
##contig=<ID=chr22_KI270734v1_random,length=165050>
##contig=<ID=chr18_GL383570v1_alt,length=164789>
##contig=<ID=chr5_KI270794v1_alt,length=164558>
##contig=<ID=chr4_GL383527v1_alt,length=164536>
##contig=<ID=chrUn_GL000213v1,length=164239>
##contig=<ID=chr3_KI270936v1_alt,length=164170>
##contig=<ID=chr3_KI270934v1_alt,length=163458>
##contig=<ID=chr9_GL383539v1_alt,length=162988>
##contig=<ID=chr3_KI270895v1_alt,length=162896>
##contig=<ID=chr22_GL383582v2_alt,length=162811>
##contig=<ID=chr3_KI270782v1_alt,length=162429>
##contig=<ID=chr1_KI270892v1_alt,length=162212>
##contig=<ID=chrUn_GL000220v1,length=161802>
##contig=<ID=chr2_KI270767v1_alt,length=161578>
##contig=<ID=chr2_KI270715v1_random,length=161471>
##contig=<ID=chr2_KI270893v1_alt,length=161218>
##contig=<ID=chrUn_GL000218v1,length=161147>
##contig=<ID=chr18_GL383572v1_alt,length=159547>
##contig=<ID=chr8_KI270817v1_alt,length=158983>
##contig=<ID=chr4_KI270788v1_alt,length=158965>
##contig=<ID=chrUn_KI270749v1,length=158759>
##contig=<ID=chr7_KI270806v1_alt,length=158166>
##contig=<ID=chr7_KI270804v1_alt,length=157952>
##contig=<ID=chr18_KI270911v1_alt,length=157710>
##contig=<ID=chrUn_KI270741v1,length=157432>
##contig=<ID=chr17_KI270910v1_alt,length=157099>
##contig=<ID=chr19_KI270884v1_alt,length=157053>
##contig=<ID=chr19_GL383574v1_alt,length=155864>
##contig=<ID=chr19_KI270888v1_alt,length=155532>
##contig=<ID=chr3_GL000221v1_random,length=155397>
##contig=<ID=chr11_GL383547v1_alt,length=154407>
##contig=<ID=chr2_KI270716v1_random,length=153799>
##contig=<ID=chr12_GL383553v2_alt,length=152874>
##contig=<ID=chr6_KI270799v1_alt,length=152148>
##contig=<ID=chr22_KI270731v1_random,length=150754>
##contig=<ID=chrUn_KI270751v1,length=150742>
##contig=<ID=chrUn_KI270750v1,length=148850>
##contig=<ID=chr8_KI270818v1_alt,length=145606>
##contig=<ID=chrX_KI270881v1_alt,length=144206>
##contig=<ID=chr21_KI270873v1_alt,length=143900>
##contig=<ID=chr2_GL383521v1_alt,length=143390>
##contig=<ID=chr8_KI270814v1_alt,length=141812>
##contig=<ID=chr12_GL383552v1_alt,length=138655>
##contig=<ID=chrUn_KI270519v1,length=138126>
##contig=<ID=chr2_KI270775v1_alt,length=138019>
##contig=<ID=chr17_KI270907v1_alt,length=137721>
##contig=<ID=chrUn_GL000214v1,length=137718>
##contig=<ID=chr8_KI270901v1_alt,length=136959>
##contig=<ID=chr2_KI270770v1_alt,length=136240>
##contig=<ID=chr16_KI270854v1_alt,length=134193>
##contig=<ID=chr8_KI270819v1_alt,length=133535>
##contig=<ID=chr17_GL383564v2_alt,length=133151>
##contig=<ID=chr2_KI270772v1_alt,length=133041>
##contig=<ID=chr8_KI270815v1_alt,length=132244>
##contig=<ID=chr5_KI270795v1_alt,length=131892>
##contig=<ID=chr5_KI270898v1_alt,length=130957>
##contig=<ID=chr20_GL383577v2_alt,length=128386>
##contig=<ID=chr1_KI270708v1_random,length=127682>
##contig=<ID=chr7_KI270807v1_alt,length=126434>
##contig=<ID=chr5_KI270793v1_alt,length=126136>
##contig=<ID=chr6_GL383533v1_alt,length=124736>
##contig=<ID=chr2_GL383522v1_alt,length=123821>
##contig=<ID=chr19_KI270918v1_alt,length=123111>
##contig=<ID=chr12_GL383549v1_alt,length=120804>
##contig=<ID=chr2_KI270769v1_alt,length=120616>
##contig=<ID=chr4_KI270785v1_alt,length=119912>
##contig=<ID=chr12_KI270834v1_alt,length=119498>
##contig=<ID=chr7_GL383534v2_alt,length=119183>
##contig=<ID=chr20_KI270869v1_alt,length=118774>
##contig=<ID=chr21_GL383581v2_alt,length=116689>
##contig=<ID=chr3_KI270781v1_alt,length=113034>
##contig=<ID=chr17_KI270730v1_random,length=112551>
##contig=<ID=chrUn_KI270438v1,length=112505>
##contig=<ID=chr4_KI270787v1_alt,length=111943>
##contig=<ID=chr18_KI270864v1_alt,length=111737>
##contig=<ID=chr2_KI270771v1_alt,length=110395>
##contig=<ID=chr1_GL383519v1_alt,length=110268>
##contig=<ID=chr2_KI270768v1_alt,length=110099>
##contig=<ID=chr1_KI270760v1_alt,length=109528>
##contig=<ID=chr3_KI270783v1_alt,length=109187>
##contig=<ID=chr17_KI270859v1_alt,length=108763>
##contig=<ID=chr11_KI270902v1_alt,length=106711>
##contig=<ID=chr18_GL383568v1_alt,length=104552>
##contig=<ID=chr22_KI270737v1_random,length=103838>
##contig=<ID=chr13_KI270843v1_alt,length=103832>
##contig=<ID=chr22_KI270877v1_alt,length=101331>
##contig=<ID=chr5_GL383530v1_alt,length=101241>
##contig=<ID=chr11_KI270721v1_random,length=100316>
##contig=<ID=chr22_KI270738v1_random,length=99375>
##contig=<ID=chr22_GL383583v2_alt,length=96924>
##contig=<ID=chr2_GL582966v2_alt,length=96131>
##contig=<ID=chrUn_KI270748v1,length=93321>
##contig=<ID=chrUn_KI270435v1,length=92983>
##contig=<ID=chr5_GL000208v1_random,length=92689>
##contig=<ID=chrUn_KI270538v1,length=91309>
##contig=<ID=chr17_GL383566v1_alt,length=90219>
##contig=<ID=chr16_GL383557v1_alt,length=89672>
##contig=<ID=chr17_JH159148v1_alt,length=88070>
##contig=<ID=chr5_GL383532v1_alt,length=82728>
##contig=<ID=chr21_KI270872v1_alt,length=82692>
##contig=<ID=chrUn_KI270756v1,length=79590>
##contig=<ID=chr6_KI270758v1_alt,length=76752>
##contig=<ID=chr12_KI270833v1_alt,length=76061>
##contig=<ID=chr6_KI270802v1_alt,length=75005>
##contig=<ID=chr21_GL383580v2_alt,length=74653>
##contig=<ID=chr22_KB663609v1_alt,length=74013>
##contig=<ID=chr22_KI270739v1_random,length=73985>
##contig=<ID=chr9_GL383540v1_alt,length=71551>
##contig=<ID=chrUn_KI270757v1,length=71251>
##contig=<ID=chr2_KI270773v1_alt,length=70887>
##contig=<ID=chr17_JH159147v1_alt,length=70345>
##contig=<ID=chr11_KI270827v1_alt,length=67707>
##contig=<ID=chr1_KI270709v1_random,length=66860>
##contig=<ID=chrUn_KI270746v1,length=66486>
##contig=<ID=chr16_KI270856v1_alt,length=63982>
##contig=<ID=chr21_GL383578v2_alt,length=63917>
##contig=<ID=chrUn_KN707963v1_decoy,length=62955>
##contig=<ID=chrUn_KI270753v1,length=62944>
##contig=<ID=chr19_KI270868v1_alt,length=61734>
##contig=<ID=chr9_GL383542v1_alt,length=60032>
##contig=<ID=chr20_KI270871v1_alt,length=58661>
##contig=<ID=chr12_KI270836v1_alt,length=56134>
##contig=<ID=chr19_KI270865v1_alt,length=52969>
##contig=<ID=chr1_KI270764v1_alt,length=50258>
##contig=<ID=chrUn_KI270589v1,length=44474>
##contig=<ID=chr14_KI270726v1_random,length=43739>
##contig=<ID=chr19_KI270866v1_alt,length=43156>
##contig=<ID=chr22_KI270735v1_random,length=42811>
##contig=<ID=chr1_KI270711v1_random,length=42210>
##contig=<ID=chrUn_KI270745v1,length=41891>
##contig=<ID=chr1_KI270714v1_random,length=41717>
##contig=<ID=chr22_KI270732v1_random,length=41543>
##contig=<ID=chr1_KI270713v1_random,length=40745>
##contig=<ID=chrUn_KI270754v1,length=40191>
##contig=<ID=chr1_KI270710v1_random,length=40176>
##contig=<ID=chr12_KI270837v1_alt,length=40090>
##contig=<ID=chr9_KI270717v1_random,length=40062>
##contig=<ID=chr14_KI270724v1_random,length=39555>
##contig=<ID=chr9_KI270720v1_random,length=39050>
##contig=<ID=chr14_KI270723v1_random,length=38115>
##contig=<ID=chr9_KI270718v1_random,length=38054>
##contig=<ID=chrUn_KI270317v1,length=37690>
##contig=<ID=chr13_KI270842v1_alt,length=37287>
##contig=<ID=chrY_KI270740v1_random,length=37240>
##contig=<ID=chrUn_KI270755v1,length=36723>
##contig=<ID=chr8_KI270820v1_alt,length=36640>
##contig=<ID=chrUn_KN707896v1_decoy,length=33125>
##contig=<ID=chr1_KI270707v1_random,length=32032>
##contig=<ID=chrUn_KI270579v1,length=31033>
##contig=<ID=chrUn_KI270752v1,length=27745>
##contig=<ID=chrUn_JTFH01000001v1_decoy,length=25139>
##contig=<ID=chrUn_KI270512v1,length=22689>
##contig=<ID=chrUn_KI270322v1,length=21476>
##contig=<ID=chrUn_KN707908v1_decoy,length=19837>
##contig=<ID=chrUn_KN707673v1_decoy,length=19544>
##contig=<ID=chrUn_JTFH01001735v1_decoy,length=19311>
##contig=<ID=chrUn_JTFH01000002v1_decoy,length=18532>
##contig=<ID=chrUn_JTFH01001257v1_decoy,length=17929>
##contig=<ID=chrUn_KN707798v1_decoy,length=17599>
##contig=<ID=chrM,length=16569>
##contig=<ID=chrUn_KN707862v1_decoy,length=16513>
##contig=<ID=HLA-DRB1*07:01:01:02,length=16120>
##contig=<ID=HLA-DRB1*07:01:01:01,length=16110>
##contig=<ID=HLA-DRB1*09:21,length=16039>
##contig=<ID=chrUn_KN707967v1_decoy,length=15368>
##contig=<ID=HLA-DRB1*04:03:01,length=15246>
##contig=<ID=chrUn_JTFH01000003v1_decoy,length=15240>
##contig=<ID=chrUn_GL000226v1,length=15008>
##contig=<ID=chrUn_KN707728v1_decoy,length=14625>
##contig=<ID=HLA-DRB1*13:02:01,length=13941>
##contig=<ID=HLA-DRB1*14:54:01,length=13936>
##contig=<ID=HLA-DRB1*13:01:01,length=13935>
##contig=<ID=HLA-DRB1*14:05:01,length=13933>
##contig=<ID=HLA-DRB1*11:01:02,length=13931>
##contig=<ID=HLA-DRB1*11:01:01,length=13921>
##contig=<ID=HLA-DRB1*11:04:01,length=13919>
##contig=<ID=HLA-DRB1*03:01:01:01,length=13908>
##contig=<ID=chrUn_JTFH01000004v1_decoy,length=13739>
##contig=<ID=HLA-DRB1*08:03:02,length=13562>
##contig=<ID=HLA-DRB1*10:01:01,length=13501>
##contig=<ID=HLA-DRB1*03:01:01:02,length=13426>
##contig=<ID=HLA-DRB1*12:01:01,length=13404>
##contig=<ID=chrUn_KN707971v1_decoy,length=12943>
##contig=<ID=chrUn_KN707779v1_decoy,length=12444>
##contig=<ID=chrUn_KI270311v1,length=12399>
##contig=<ID=chrUn_JTFH01001736v1_decoy,length=11713>
##contig=<ID=HLA-DRB1*15:01:01:02,length=11571>
##contig=<ID=HLA-DRB1*15:03:01:02,length=11569>
##contig=<ID=HLA-DRB1*15:03:01:01,length=11567>
##contig=<ID=chrUn_KN707925v1_decoy,length=11555>
##contig=<ID=chrUn_KN707957v1_decoy,length=11499>
##contig=<ID=chrUn_JTFH01000005v1_decoy,length=11297>
##contig=<ID=chrUn_JTFH01001737v1_decoy,length=11263>
##contig=<ID=HLA-DRB1*12:17,length=11260>
##contig=<ID=HLA-DRB1*01:02:01,length=11229>
##contig=<ID=HLA-DRB1*15:01:01:01,length=11080>
##contig=<ID=HLA-DRB1*15:01:01:04,length=11056>
##contig=<ID=HLA-DRB1*15:01:01:03,length=11056>
##contig=<ID=HLA-DRB1*16:02:01,length=11005>
##contig=<ID=HLA-DRB1*01:01:01,length=10741>
##contig=<ID=chrUn_KN707675v1_decoy,length=10556>
##contig=<ID=chrUn_KN707854v1_decoy,length=10451>
##contig=<ID=HLA-DRB1*15:02:01,length=10313>
##contig=<ID=chrUn_KN707777v1_decoy,length=10311>
##contig=<ID=chrUn_KN707740v1_decoy,length=10282>
##contig=<ID=chrUn_JTFH01000006v1_decoy,length=10074>
##contig=<ID=chrUn_JTFH01000007v1_decoy,length=9891>
##contig=<ID=chrUn_JTFH01001738v1_decoy,length=9779>
##contig=<ID=chrUn_JTFH01000008v1_decoy,length=9774>
##contig=<ID=chrUn_JTFH01001258v1_decoy,length=9749>
##contig=<ID=chrUn_JTFH01000009v1_decoy,length=9727>
##contig=<ID=chrUn_KN707734v1_decoy,length=9652>
##contig=<ID=chrUn_KN707853v1_decoy,length=9597>
##contig=<ID=chrUn_KN707968v1_decoy,length=9572>
##contig=<ID=chrUn_JTFH01001739v1_decoy,length=9568>
##contig=<ID=chrUn_JTFH01000010v1_decoy,length=9358>
##contig=<ID=chrUn_JTFH01001740v1_decoy,length=9344>
##contig=<ID=chrUn_JTFH01001741v1_decoy,length=9188>
##contig=<ID=chrUn_JTFH01001742v1_decoy,length=9100>
##contig=<ID=chrUn_KN707676v1_decoy,length=9066>
##contig=<ID=chrUn_KN707841v1_decoy,length=8992>
##contig=<ID=chrUn_JTFH01000011v1_decoy,length=8920>
##contig=<ID=chrUn_KN707772v1_decoy,length=8915>
##contig=<ID=chrUn_JTFH01001743v1_decoy,length=8771>
##contig=<ID=chrUn_JTFH01001744v1_decoy,length=8690>
##contig=<ID=chrUn_KN707664v1_decoy,length=8683>
##contig=<ID=chrUn_KN707869v1_decoy,length=8632>
##contig=<ID=chrUn_JTFH01001745v1_decoy,length=8566>
##contig=<ID=chrUn_JTFH01000012v1_decoy,length=8479>
##contig=<ID=chrUn_KN707660v1_decoy,length=8410>
##contig=<ID=chrUn_JTFH01000962v1_decoy,length=8358>
##contig=<ID=chrUn_KI270366v1,length=8320>
##contig=<ID=chrUn_JTFH01000013v1_decoy,length=8312>
##contig=<ID=chrUn_KN707820v1_decoy,length=8307>
##contig=<ID=chrUn_JTFH01000014v1_decoy,length=8261>
##contig=<ID=chrUn_JTFH01000015v1_decoy,length=8131>
##contig=<ID=chrUn_KI270511v1,length=8127>
##contig=<ID=chrUn_JTFH01001746v1_decoy,length=8058>
##contig=<ID=chrUn_JTFH01001259v1_decoy,length=8053>
##contig=<ID=chrUn_JTFH01000016v1_decoy,length=8051>
##contig=<ID=chrUn_KN707893v1_decoy,length=8049>
##contig=<ID=chrUn_KI270448v1,length=7992>
##contig=<ID=chrUn_JTFH01000963v1_decoy,length=7932>
##contig=<ID=chrUn_JTFH01000017v1_decoy,length=7832>
##contig=<ID=chrUn_JTFH01001260v1_decoy,length=7826>
##contig=<ID=chrUn_KN707863v1_decoy,length=7821>
##contig=<ID=chrUn_KN707907v1_decoy,length=7776>
##contig=<ID=chrUn_JTFH01001261v1_decoy,length=7768>
##contig=<ID=chrUn_JTFH01001747v1_decoy,length=7759>
##contig=<ID=chrUn_JTFH01000018v1_decoy,length=7710>
##contig=<ID=chrUn_JTFH01000019v1_decoy,length=7702>
##contig=<ID=chrUn_KI270521v1,length=7642>
##contig=<ID=chrUn_KN707873v1_decoy,length=7591>
##contig=<ID=chrUn_JTFH01001748v1_decoy,length=7585>
##contig=<ID=chrUn_KN707733v1_decoy,length=7503>
##contig=<ID=HLA-DQB1*02:01:01,length=7480>
##contig=<ID=chrUn_JTFH01000020v1_decoy,length=7479>
##contig=<ID=chrUn_JTFH01001749v1_decoy,length=7471>
##contig=<ID=HLA-DQB1*02:02:01,length=7471>
##contig=<ID=chrUn_JTFH01001750v1_decoy,length=7461>
##contig=<ID=chrUn_KN707937v1_decoy,length=7445>
##contig=<ID=chrUn_KN707729v1_decoy,length=7431>
##contig=<ID=chrUn_JTFH01000021v1_decoy,length=7368>
##contig=<ID=chrUn_JTFH01001751v1_decoy,length=7342>
##contig=<ID=chrUn_KN707677v1_decoy,length=7267>
##contig=<ID=HLA-DQB1*03:01:01:01,length=7231>
##contig=<ID=HLA-DQB1*03:01:01:03,length=7231>
##contig=<ID=HLA-DQB1*03:01:01:02,length=7230>
##contig=<ID=chrUn_JTFH01001752v1_decoy,length=7223>
##contig=<ID=chrUn_KN707794v1_decoy,length=7219>
##contig=<ID=chrUn_JTFH01000022v1_decoy,length=7162>
##contig=<ID=chrUn_KN707816v1_decoy,length=7150>
##contig=<ID=HLA-DQB1*03:03:02:01,length=7126>
##contig=<ID=HLA-DQB1*03:02:01,length=7126>
##contig=<ID=HLA-DQB1*03:03:02:02,length=7126>
##contig=<ID=HLA-DQB1*06:01:01,length=7111>
##contig=<ID=HLA-DQB1*06:03:01,length=7103>
##contig=<ID=HLA-DQB1*06:02:01,length=7102>
##contig=<ID=HLA-DQB1*06:09:01,length=7102>
##contig=<ID=HLA-DQB1*05:01:01:02,length=7090>
##contig=<ID=HLA-DQB1*05:01:01:01,length=7090>
##contig=<ID=HLA-DQB1*05:03:01:01,length=7089>
##contig=<ID=HLA-DQB1*05:03:01:02,length=7089>
##contig=<ID=chrUn_JTFH01000023v1_decoy,length=7065>
##contig=<ID=chrUn_JTFH01001753v1_decoy,length=7064>
##contig=<ID=chrUn_KI270581v1,length=7046>
##contig=<ID=chrUn_JTFH01000024v1_decoy,length=7019>
##contig=<ID=chrUn_KN707749v1_decoy,length=7010>
##contig=<ID=chrUn_JTFH01000025v1_decoy,length=6997>
##contig=<ID=chrUn_JTFH01000026v1_decoy,length=6994>
##contig=<ID=chrUn_JTFH01000027v1_decoy,length=6979>
##contig=<ID=HLA-DQB1*03:05:01,length=6934>
##contig=<ID=chrUn_JTFH01001754v1_decoy,length=6916>
##contig=<ID=chrUn_JTFH01001755v1_decoy,length=6897>
##contig=<ID=chrUn_JTFH01001756v1_decoy,length=6880>
##contig=<ID=chrUn_JTFH01001757v1_decoy,length=6857>
##contig=<ID=chrUn_JTFH01000964v1_decoy,length=6846>
##contig=<ID=chrUn_JTFH01001758v1_decoy,length=6840>
##contig=<ID=HLA-DQB1*03:03:02:03,length=6800>
##contig=<ID=chrUn_JTFH01000028v1_decoy,length=6797>
##contig=<ID=chrUn_JTFH01001759v1_decoy,length=6728>
##contig=<ID=chrUn_KN707701v1_decoy,length=6722>
##contig=<ID=chrUn_KN707842v1_decoy,length=6698>
##contig=<ID=chrUn_JTFH01001760v1_decoy,length=6688>
##contig=<ID=HLA-DQA1*05:05:01:02,length=6597>
##contig=<ID=HLA-DQA1*05:05:01:01,length=6593>
##contig=<ID=HLA-DQA1*05:11,length=6589>
##contig=<ID=chrUn_KN707903v1_decoy,length=6581>
##contig=<ID=chrUn_JTFH01001761v1_decoy,length=6553>
##contig=<ID=HLA-DQA1*05:01:01:02,length=6529>
##contig=<ID=chrUn_JTFH01000029v1_decoy,length=6525>
##contig=<ID=chrUn_KI270582v1,length=6504>
##contig=<ID=HLA-DQA1*01:03:01:02,length=6492>
##contig=<ID=HLA-DQA1*01:02:01:04,length=6492>
##contig=<ID=HLA-DQA1*01:01:02,length=6489>
##contig=<ID=chrUn_KN707966v1_decoy,length=6486>
##contig=<ID=HLA-DQA1*01:04:01:02,length=6485>
##contig=<ID=HLA-DQA1*01:05:01,length=6485>
##contig=<ID=HLA-DQA1*01:02:01:02,length=6485>
##contig=<ID=HLA-DQA1*01:03:01:01,length=6485>
##contig=<ID=HLA-DQA1*01:02:01:03,length=6485>
##contig=<ID=HLA-DQA1*01:02:01:01,length=6484>
##contig=<ID=HLA-DQA1*01:04:01:01,length=6484>
##contig=<ID=chrUn_KN707702v1_decoy,length=6466>
##contig=<ID=HLA-DQA1*03:02,length=6437>
##contig=<ID=HLA-DQA1*03:01:01,length=6437>
##contig=<ID=HLA-DQA1*03:03:01,length=6437>
##contig=<ID=chrUn_KN707969v1_decoy,length=6413>
##contig=<ID=HLA-DQA1*02:01,length=6403>
##contig=<ID=chrUn_JTFH01001762v1_decoy,length=6396>
##contig=<ID=HLA-DQA1*05:05:01:03,length=6393>
##contig=<ID=chrUn_KI270515v1,length=6361>
##contig=<ID=chrUn_JTFH01001763v1_decoy,length=6345>
##contig=<ID=chrUn_JTFH01001764v1_decoy,length=6295>
##contig=<ID=chrUn_KN707964v1_decoy,length=6282>
##contig=<ID=chrUn_JTFH01001765v1_decoy,length=6266>
##contig=<ID=chrUn_KN707806v1_decoy,length=6252>
##contig=<ID=chrUn_JTFH01000030v1_decoy,length=6246>
##contig=<ID=chrUn_KN707975v1_decoy,length=6211>
##contig=<ID=HLA-DQA1*04:02,length=6210>
##contig=<ID=chrUn_JTFH01001766v1_decoy,length=6173>
##contig=<ID=chrUn_JTFH01001767v1_decoy,length=6171>
##contig=<ID=chrUn_KI270588v1,length=6158>
##contig=<ID=HLA-DQA1*05:03,length=6121>
##contig=<ID=chrUn_JTFH01001768v1_decoy,length=6120>
##contig=<ID=chrUn_JTFH01001769v1_decoy,length=6105>
##contig=<ID=chrUn_JTFH01001770v1_decoy,length=6099>
##contig=<ID=chrUn_KN707720v1_decoy,length=6028>
##contig=<ID=chrUn_KN707765v1_decoy,length=6020>
##contig=<ID=chrUn_KN707775v1_decoy,length=5997>
##contig=<ID=chrUn_KN707857v1_decoy,length=5985>
##contig=<ID=HLA-DQA1*01:07,length=5959>
##contig=<ID=chrUn_JTFH01000031v1_decoy,length=5926>
##contig=<ID=HLA-DQA1*01:11,length=5926>
##contig=<ID=chrUn_JTFH01000032v1_decoy,length=5914>
##contig=<ID=chrUn_JTFH01000033v1_decoy,length=5898>
##contig=<ID=chrUn_JTFH01001771v1_decoy,length=5893>
##contig=<ID=chrUn_JTFH01000034v1_decoy,length=5879>
##contig=<ID=HLA-DQA1*06:01:01,length=5878>
##contig=<ID=HLA-DQA1*04:01:02:01,length=5853>
##contig=<ID=chrUn_JTFH01000035v1_decoy,length=5834>
##contig=<ID=chrUn_JTFH01001772v1_decoy,length=5829>
##contig=<ID=chrUn_KN707689v1_decoy,length=5823>
##contig=<ID=chrUn_KN707627v1_decoy,length=5821>
##contig=<ID=HLA-DQA1*05:01:01:01,length=5806>
##contig=<ID=chrUn_KI270591v1,length=5796>
##contig=<ID=chrUn_JTFH01001773v1_decoy,length=5793>
##contig=<ID=HLA-DQA1*01:10,length=5790>
##contig=<ID=chrUn_JTFH01001774v1_decoy,length=5776>
##contig=<ID=chrUn_KN707730v1_decoy,length=5776>
##contig=<ID=chrUn_JTFH01001775v1_decoy,length=5759>
##contig=<ID=chrUn_JTFH01000036v1_decoy,length=5743>
##contig=<ID=chrUn_JTFH01001776v1_decoy,length=5716>
##contig=<ID=chrUn_JTFH01001777v1_decoy,length=5708>
##contig=<ID=chrUn_KN707710v1_decoy,length=5701>
##contig=<ID=chrUn_KN707780v1_decoy,length=5691>
##contig=<ID=chrUn_JTFH01001262v1_decoy,length=5691>
##contig=<ID=chrUn_KI270522v1,length=5674>
##contig=<ID=HLA-DQA1*04:01:02:02,length=5666>
##contig=<ID=chrUn_KN707626v1_decoy,length=5623>
##contig=<ID=chrUn_KN707741v1_decoy,length=5601>
##contig=<ID=chrUn_KN707939v1_decoy,length=5600>
##contig=<ID=chrUn_JTFH01001778v1_decoy,length=5590>
##contig=<ID=chrUn_JTFH01000037v1_decoy,length=5577>
##contig=<ID=chrUn_JTFH01001779v1_decoy,length=5566>
##contig=<ID=chrUn_JTFH01001780v1_decoy,length=5558>
##contig=<ID=chrUn_KN707661v1_decoy,length=5534>
##contig=<ID=chrUn_KN707805v1_decoy,length=5446>
##contig=<ID=chrUn_JTFH01001263v1_decoy,length=5444>
##contig=<ID=chrUn_JTFH01001781v1_decoy,length=5418>
##contig=<ID=chrUn_JTFH01000038v1_decoy,length=5413>
##contig=<ID=chrUn_KN707815v1_decoy,length=5410>
##contig=<ID=chrUn_KN707786v1_decoy,length=5385>
##contig=<ID=chrUn_JTFH01001782v1_decoy,length=5375>
##contig=<ID=chrUn_KN707637v1_decoy,length=5354>
##contig=<ID=chrUn_KI270507v1,length=5353>
##contig=<ID=chrUn_KN707755v1_decoy,length=5343>
##contig=<ID=chrUn_KN707803v1_decoy,length=5336>
##contig=<ID=chrUn_KN707943v1_decoy,length=5325>
##contig=<ID=chrUn_JTFH01001783v1_decoy,length=5300>
##contig=<ID=chrUn_KN707782v1_decoy,length=5277>
##contig=<ID=chrUn_JTFH01001784v1_decoy,length=5255>
##contig=<ID=chrUn_JTFH01000039v1_decoy,length=5250>
##contig=<ID=chrUn_JTFH01000040v1_decoy,length=5246>
##contig=<ID=chrUn_KN707874v1_decoy,length=5202>
##contig=<ID=chrUn_JTFH01001785v1_decoy,length=5157>
##contig=<ID=chrUn_JTFH01001786v1_decoy,length=5130>
##contig=<ID=chrUn_JTFH01000041v1_decoy,length=5118>
##contig=<ID=chrUn_KN707799v1_decoy,length=5095>
##contig=<ID=chrUn_KN707746v1_decoy,length=5083>
##contig=<ID=chrUn_JTFH01001264v1_decoy,length=5077>
##contig=<ID=chrUn_JTFH01000042v1_decoy,length=5058>
##contig=<ID=chrUn_JTFH01001265v1_decoy,length=4990>
##contig=<ID=chrUn_JTFH01001787v1_decoy,length=4978>
##contig=<ID=chrUn_KN707718v1_decoy,length=4969>
##contig=<ID=chrUn_JTFH01000043v1_decoy,length=4959>
##contig=<ID=chrUn_JTFH01001788v1_decoy,length=4957>
##contig=<ID=chrUn_JTFH01001789v1_decoy,length=4947>
##contig=<ID=chrUn_KN707773v1_decoy,length=4902>
##contig=<ID=chrUn_JTFH01001790v1_decoy,length=4897>
##contig=<ID=chrUn_KN707691v1_decoy,length=4885>
##contig=<ID=chrUn_KN707843v1_decoy,length=4880>
##contig=<ID=chrUn_KN707756v1_decoy,length=4877>
##contig=<ID=chrUn_JTFH01001791v1_decoy,length=4867>
##contig=<ID=chrUn_JTFH01000044v1_decoy,length=4853>
##contig=<ID=chrUn_KN707919v1_decoy,length=4850>
##contig=<ID=chrUn_JTFH01001792v1_decoy,length=4845>
##contig=<ID=chrUn_JTFH01000045v1_decoy,length=4828>
##contig=<ID=chrUn_KN707731v1_decoy,length=4820>
##contig=<ID=chrUn_JTFH01000046v1_decoy,length=4819>
##contig=<ID=chrUn_KN707692v1_decoy,length=4813>
##contig=<ID=chrUn_JTFH01000047v1_decoy,length=4809>
##contig=<ID=chrUn_KN707830v1_decoy,length=4807>
##contig=<ID=chrUn_KN707814v1_decoy,length=4764>
##contig=<ID=chrUn_KN707742v1_decoy,length=4758>
##contig=<ID=chrUn_KN707871v1_decoy,length=4728>
##contig=<ID=chrUn_KN707750v1_decoy,length=4718>
##contig=<ID=chrUn_JTFH01000048v1_decoy,length=4710>
##contig=<ID=chrUn_KI270590v1,length=4685>
##contig=<ID=chrUn_JTFH01000049v1_decoy,length=4680>
##contig=<ID=chrUn_JTFH01001793v1_decoy,length=4678>
##contig=<ID=chrUn_KN707793v1_decoy,length=4667>
##contig=<ID=chrUn_KN707972v1_decoy,length=4650>
##contig=<ID=chrUn_JTFH01000050v1_decoy,length=4645>
##contig=<ID=chrUn_JTFH01001794v1_decoy,length=4641>
##contig=<ID=chrUn_KN707705v1_decoy,length=4632>
##contig=<ID=chrUn_KN707807v1_decoy,length=4616>
##contig=<ID=chrUn_KN707924v1_decoy,length=4596>
##contig=<ID=chrUn_KN707851v1_decoy,length=4593>
##contig=<ID=chrUn_JTFH01001795v1_decoy,length=4592>
##contig=<ID=chrUn_JTFH01000965v1_decoy,length=4591>
##contig=<ID=chrUn_KN707884v1_decoy,length=4568>
##contig=<ID=chrUn_KN707810v1_decoy,length=4563>
##contig=<ID=chrUn_JTFH01001266v1_decoy,length=4545>
##contig=<ID=chrUn_JTFH01001267v1_decoy,length=4544>
##contig=<ID=chrUn_JTFH01001796v1_decoy,length=4543>
##contig=<ID=chrUn_KN707881v1_decoy,length=4543>
##contig=<ID=chrUn_JTFH01001797v1_decoy,length=4532>
##contig=<ID=chrUn_JTFH01000051v1_decoy,length=4514>
##contig=<ID=chrUn_KI270584v1,length=4513>
##contig=<ID=chrUn_JTFH01001798v1_decoy,length=4503>
##contig=<ID=chrUn_JTFH01001799v1_decoy,length=4495>
##contig=<ID=chrUn_JTFH01001800v1_decoy,length=4444>
##contig=<ID=chrUn_JTFH01000052v1_decoy,length=4439>
##contig=<ID=chrUn_JTFH01000053v1_decoy,length=4416>
##contig=<ID=chrUn_KI270320v1,length=4416>
##contig=<ID=chrUn_JTFH01001801v1_decoy,length=4414>
##contig=<ID=chrUn_JTFH01001802v1_decoy,length=4409>
##contig=<ID=chrUn_JTFH01000054v1_decoy,length=4409>
##contig=<ID=chrUn_JTFH01000055v1_decoy,length=4392>
##contig=<ID=chrUn_KN707804v1_decoy,length=4383>
##contig=<ID=chrUn_KN707681v1_decoy,length=4379>
##contig=<ID=chrUn_KN707783v1_decoy,length=4373>
##contig=<ID=chrUn_KN707738v1_decoy,length=4365>
##contig=<ID=chrUn_JTFH01000056v1_decoy,length=4359>
##contig=<ID=chrUn_KN707879v1_decoy,length=4346>
##contig=<ID=chrUn_KN707707v1_decoy,length=4339>
##contig=<ID=chrUn_JTFH01000057v1_decoy,length=4319>
##contig=<ID=chrUn_JTFH01001803v1_decoy,length=4302>
##contig=<ID=chrUn_JTFH01001804v1_decoy,length=4300>
##contig=<ID=chrUn_JTFH01000058v1_decoy,length=4290>
##contig=<ID=chrUn_KN707895v1_decoy,length=4284>
##contig=<ID=chrUn_KN707739v1_decoy,length=4284>
##contig=<ID=chrUn_JTFH01001805v1_decoy,length=4277>
##contig=<ID=chrUn_KN707688v1_decoy,length=4248>
##contig=<ID=chrUn_KN707797v1_decoy,length=4243>
##contig=<ID=chrUn_JTFH01000059v1_decoy,length=4242>
##contig=<ID=chrUn_KN707875v1_decoy,length=4241>
##contig=<ID=chrUn_JTFH01000060v1_decoy,length=4228>
##contig=<ID=chrUn_KN707706v1_decoy,length=4225>
##contig=<ID=chrUn_JTFH01000061v1_decoy,length=4222>
##contig=<ID=chrUn_JTFH01000062v1_decoy,length=4216>
##contig=<ID=chrUn_KI270382v1,length=4215>
##contig=<ID=chrUn_JTFH01000063v1_decoy,length=4210>
##contig=<ID=chrUn_KN707682v1_decoy,length=4208>
##contig=<ID=chrUn_JTFH01000064v1_decoy,length=4206>
##contig=<ID=chrUn_JTFH01001268v1_decoy,length=4202>
##contig=<ID=chrUn_JTFH01001269v1_decoy,length=4195>
##contig=<ID=chrUn_JTFH01001806v1_decoy,length=4173>
##contig=<ID=chrUn_JTFH01001807v1_decoy,length=4169>
##contig=<ID=chrUn_KN707711v1_decoy,length=4154>
##contig=<ID=chrUn_JTFH01001808v1_decoy,length=4136>
##contig=<ID=chrUn_KN707876v1_decoy,length=4131>
##contig=<ID=chrUn_KN707866v1_decoy,length=4109>
##contig=<ID=chrUn_KN707970v1_decoy,length=4104>
##contig=<ID=chrUn_JTFH01000065v1_decoy,length=4102>
##contig=<ID=chrUn_JTFH01000066v1_decoy,length=4101>
##contig=<ID=chrUn_JTFH01001809v1_decoy,length=4101>
##contig=<ID=chrUn_JTFH01001810v1_decoy,length=4089>
##contig=<ID=chrUn_JTFH01000067v1_decoy,length=4083>
##contig=<ID=chrUn_KN707683v1_decoy,length=4068>
##contig=<ID=chrUn_KI270468v1,length=4055>
##contig=<ID=chrUn_KN707990v1_decoy,length=4048>
##contig=<ID=chrUn_JTFH01000966v1_decoy,length=4041>
##contig=<ID=chrUn_KN707936v1_decoy,length=4031>
##contig=<ID=chrUn_KN707744v1_decoy,length=4024>
##contig=<ID=chrUn_JTFH01001811v1_decoy,length=4015>
##contig=<ID=chrUn_JTFH01001812v1_decoy,length=4000>
##contig=<ID=chrUn_JTFH01001813v1_decoy,length=3973>
##contig=<ID=chrUn_KN707823v1_decoy,length=3970>
##contig=<ID=chrUn_KN707904v1_decoy,length=3968>
##contig=<ID=chrUn_JTFH01000068v1_decoy,length=3967>
##contig=<ID=chrUn_JTFH01000069v1_decoy,length=3955>
##contig=<ID=chrUn_JTFH01000070v1_decoy,length=3945>
##contig=<ID=chrUn_KN707685v1_decoy,length=3938>
##contig=<ID=chrUn_JTFH01000071v1_decoy,length=3930>
##contig=<ID=chrUn_JTFH01000072v1_decoy,length=3929>
##contig=<ID=chrUn_JTFH01000073v1_decoy,length=3924>
##contig=<ID=chrUn_KI270467v1,length=3920>
##contig=<ID=chrUn_JTFH01000074v1_decoy,length=3919>
##contig=<ID=chrUn_KN707764v1_decoy,length=3912>
##contig=<ID=chrUn_JTFH01000075v1_decoy,length=3908>
##contig=<ID=chrUn_JTFH01000076v1_decoy,length=3892>
##contig=<ID=chrUn_JTFH01000077v1_decoy,length=3890>
##contig=<ID=chrUn_KN707986v1_decoy,length=3869>
##contig=<ID=chrUn_JTFH01000078v1_decoy,length=3859>
##contig=<ID=chrUn_JTFH01000079v1_decoy,length=3846>
##contig=<ID=chrUn_KN707812v1_decoy,length=3845>
##contig=<ID=chrUn_JTFH01000967v1_decoy,length=3841>
##contig=<ID=chrUn_KN707898v1_decoy,length=3840>
##contig=<ID=chrUn_KN707965v1_decoy,length=3836>
##contig=<ID=chrUn_JTFH01000080v1_decoy,length=3835>
##contig=<ID=chrUn_JTFH01000081v1_decoy,length=3830>
##contig=<ID=chrUn_JTFH01000082v1_decoy,length=3828>
##contig=<ID=chrUn_KN707696v1_decoy,length=3828>
##contig=<ID=chrUn_JTFH01000083v1_decoy,length=3825>
##contig=<ID=chrUn_JTFH01000084v1_decoy,length=3821>
##contig=<ID=chrUn_JTFH01000085v1_decoy,length=3809>
##contig=<ID=chrUn_JTFH01001270v1_decoy,length=3807>
##contig=<ID=chrUn_JTFH01000086v1_decoy,length=3801>
##contig=<ID=chrUn_JTFH01000087v1_decoy,length=3799>
##contig=<ID=chrUn_KN707623v1_decoy,length=3784>
##contig=<ID=chrUn_JTFH01000968v1_decoy,length=3754>
##contig=<ID=chrUn_JTFH01000969v1_decoy,length=3743>
##contig=<ID=chrUn_JTFH01001271v1_decoy,length=3741>
##contig=<ID=chrUn_JTFH01000088v1_decoy,length=3737>
##contig=<ID=chrUn_JTFH01001814v1_decoy,length=3732>
##contig=<ID=chrUn_KN707690v1_decoy,length=3715>
##contig=<ID=chrUn_JTFH01001815v1_decoy,length=3709>
##contig=<ID=chrUn_KN707700v1_decoy,length=3703>
##contig=<ID=chrUn_JTFH01000970v1_decoy,length=3702>
##contig=<ID=chrUn_JTFH01000089v1_decoy,length=3701>
##contig=<ID=chrUn_JTFH01001272v1_decoy,length=3699>
##contig=<ID=chrUn_JTFH01000090v1_decoy,length=3698>
##contig=<ID=chrUn_JTFH01000091v1_decoy,length=3692>
##contig=<ID=chrUn_JTFH01001816v1_decoy,length=3686>
##contig=<ID=chrUn_JTFH01000092v1_decoy,length=3686>
##contig=<ID=chrUn_KN707787v1_decoy,length=3678>
##contig=<ID=chrUn_JTFH01000093v1_decoy,length=3677>
##contig=<ID=chrUn_JTFH01001817v1_decoy,length=3676>
##contig=<ID=chrUn_JTFH01001818v1_decoy,length=3673>
##contig=<ID=chrUn_JTFH01001819v1_decoy,length=3672>
##contig=<ID=chrUn_KN707809v1_decoy,length=3667>
##contig=<ID=chrUn_JTFH01000094v1_decoy,length=3664>
##contig=<ID=chrUn_KN707768v1_decoy,length=3656>
##contig=<ID=chrUn_KN707654v1_decoy,length=3644>
##contig=<ID=chrUn_JTFH01001273v1_decoy,length=3640>
##contig=<ID=chrUn_JTFH01001820v1_decoy,length=3633>
##contig=<ID=chrUn_JTFH01001821v1_decoy,length=3633>
##contig=<ID=chrUn_JTFH01000971v1_decoy,length=3625>
##contig=<ID=chrUn_JTFH01001822v1_decoy,length=3613>
##contig=<ID=chrUn_JTFH01000095v1_decoy,length=3613>
##contig=<ID=chrUn_JTFH01000096v1_decoy,length=3611>
##contig=<ID=chrUn_JTFH01000097v1_decoy,length=3606>
##contig=<ID=chrUn_JTFH01001823v1_decoy,length=3605>
##contig=<ID=chrUn_JTFH01001824v1_decoy,length=3592>
##contig=<ID=chrUn_JTFH01001825v1_decoy,length=3586>
##contig=<ID=chrUn_JTFH01001826v1_decoy,length=3584>
##contig=<ID=chrUn_KN707920v1_decoy,length=3584>
##contig=<ID=chrUn_JTFH01000098v1_decoy,length=3584>
##contig=<ID=chrUn_JTFH01000099v1_decoy,length=3581>
##contig=<ID=chrUn_JTFH01001827v1_decoy,length=3577>
##contig=<ID=chrUn_KN707748v1_decoy,length=3553>
##contig=<ID=chrUn_KN707751v1_decoy,length=3546>
##contig=<ID=chrUn_JTFH01000100v1_decoy,length=3543>
##contig=<ID=chrUn_JTFH01001828v1_decoy,length=3537>
##contig=<ID=chrUn_KN707887v1_decoy,length=3534>
##contig=<ID=chrUn_KN707845v1_decoy,length=3532>
##contig=<ID=chrUn_JTFH01001274v1_decoy,length=3531>
##contig=<ID=chrUn_KI270362v1,length=3530>
##contig=<ID=chrUn_JTFH01000972v1_decoy,length=3529>
##contig=<ID=chrUn_JTFH01000101v1_decoy,length=3528>
##contig=<ID=chrUn_JTFH01000102v1_decoy,length=3527>
##contig=<ID=HLA-A*74:02:01:02,length=3518>
##contig=<ID=HLA-A*32:01:01,length=3518>
##contig=<ID=HLA-A*29:02:01:02,length=3518>
##contig=<ID=HLA-A*29:01:01:01,length=3518>
##contig=<ID=HLA-A*31:01:02,length=3518>
##contig=<ID=HLA-A*33:03:01,length=3518>
##contig=<ID=HLA-A*29:02:01:01,length=3518>
##contig=<ID=HLA-A*33:01:01,length=3518>
##contig=<ID=HLA-A*02:48,length=3517>
##contig=<ID=HLA-A*02:251,length=3517>
##contig=<ID=HLA-A*26:01:01,length=3517>
##contig=<ID=HLA-A*68:01:01:02,length=3517>
##contig=<ID=HLA-A*34:01:01,length=3517>
##contig=<ID=HLA-A*02:06:01,length=3517>
##contig=<ID=HLA-A*68:02:01:01,length=3517>
##contig=<ID=HLA-A*02:10,length=3517>
##contig=<ID=HLA-A*66:01:01,length=3517>
##contig=<ID=HLA-A*02:05:01,length=3517>
##contig=<ID=HLA-A*02:07:01,length=3517>
##contig=<ID=HLA-A*02:03:01,length=3517>
##contig=<ID=HLA-A*02:01:01:01,length=3517>
##contig=<ID=HLA-A*68:01:02:01,length=3517>
##contig=<ID=HLA-A*02:32N,length=3517>
##contig=<ID=HLA-A*02:01:01:04,length=3516>
##contig=<ID=chrUn_JTFH01001829v1_decoy,length=3510>
##contig=<ID=chrUn_JTFH01001830v1_decoy,length=3509>
##contig=<ID=chrUn_JTFH01000973v1_decoy,length=3508>
##contig=<ID=HLA-A*68:02:01:02,length=3506>
##contig=<ID=HLA-A*01:03,length=3503>
##contig=<ID=HLA-A*30:04:01,length=3503>
##contig=<ID=HLA-A*01:01:01:01,length=3503>
##contig=<ID=HLA-A*24:11N,length=3503>
##contig=<ID=HLA-A*30:01:01,length=3503>
##contig=<ID=HLA-A*11:01:01,length=3503>
##contig=<ID=HLA-A*11:01:18,length=3503>
##contig=<ID=HLA-A*11:02:01,length=3503>
##contig=<ID=HLA-A*03:01:01:01,length=3502>
##contig=<ID=HLA-A*24:20,length=3502>
##contig=<ID=HLA-A*03:02:01,length=3502>
##contig=<ID=HLA-A*24:02:01:02L,length=3502>
##contig=<ID=HLA-A*24:09N,length=3502>
##contig=<ID=HLA-A*24:10:01,length=3502>
##contig=<ID=HLA-A*24:07:01,length=3502>
##contig=<ID=HLA-A*23:01:01,length=3502>
##contig=<ID=HLA-A*24:02:01:01,length=3502>
##contig=<ID=HLA-A*24:08,length=3502>
##contig=<ID=HLA-A*24:03:01,length=3502>
##contig=<ID=HLA-A*11:69N,length=3500>
##contig=<ID=chrUn_JTFH01000103v1_decoy,length=3496>
##contig=<ID=chrUn_JTFH01000104v1_decoy,length=3493>
##contig=<ID=chrUn_JTFH01001831v1_decoy,length=3488>
##contig=<ID=chrUn_JTFH01000105v1_decoy,length=3484>
##contig=<ID=chrUn_JTFH01001832v1_decoy,length=3473>
##contig=<ID=chrUn_KN707796v1_decoy,length=3473>
##contig=<ID=chrUn_JTFH01001275v1_decoy,length=3455>
##contig=<ID=chrUn_KN707762v1_decoy,length=3450>
##contig=<ID=chrUn_JTFH01001833v1_decoy,length=3445>
##contig=<ID=chrUn_JTFH01000106v1_decoy,length=3435>
##contig=<ID=chrUn_JTFH01001834v1_decoy,length=3427>
##contig=<ID=HLA-A*24:86N,length=3415>
##contig=<ID=chrUn_JTFH01001276v1_decoy,length=3411>
##contig=<ID=chrUn_KN707961v1_decoy,length=3410>
##contig=<ID=HLA-A*03:11N,length=3404>
##contig=<ID=chrUn_JTFH01001835v1_decoy,length=3395>
##contig=<ID=chrUn_JTFH01000107v1_decoy,length=3391>
##contig=<ID=HLA-A*32:06,length=3389>
##contig=<ID=HLA-A*33:07,length=3389>
##contig=<ID=HLA-A*43:01,length=3388>
##contig=<ID=HLA-A*02:95,length=3388>
##contig=<ID=HLA-A*68:01:02:02,length=3388>
##contig=<ID=HLA-A*02:65,length=3387>
##contig=<ID=chrUn_JTFH01001277v1_decoy,length=3387>
##contig=<ID=chrUn_JTFH01000108v1_decoy,length=3374>
##contig=<ID=HLA-A*01:01:38L,length=3374>
##contig=<ID=HLA-A*01:11N,length=3374>
##contig=<ID=HLA-A*01:02,length=3374>
##contig=<ID=HLA-A*30:02:01:02,length=3374>
##contig=<ID=HLA-A*11:05,length=3373>
##contig=<ID=HLA-A*03:01:01:02N,length=3373>
##contig=<ID=HLA-A*02:89,length=3371>
##contig=<ID=chrUn_JTFH01000109v1_decoy,length=3371>
##contig=<ID=HLA-A*02:77,length=3371>
##contig=<ID=HLA-C*17:01:01:01,length=3368>
##contig=<ID=HLA-C*17:01:01:02,length=3368>
##contig=<ID=HLA-C*17:01:01:03,length=3368>
##contig=<ID=chrUn_JTFH01001836v1_decoy,length=3367>
##contig=<ID=HLA-A*11:50Q,length=3362>
##contig=<ID=chrUn_JTFH01000110v1_decoy,length=3361>
##contig=<ID=chrUn_JTFH01000974v1_decoy,length=3359>
##contig=<ID=chrUn_JTFH01001278v1_decoy,length=3358>
##contig=<ID=HLA-A*24:02:10,length=3356>
##contig=<ID=HLA-C*07:66,length=3354>
##contig=<ID=HLA-C*07:02:01:05,length=3354>
##contig=<ID=HLA-C*07:06,length=3354>
##contig=<ID=HLA-C*07:02:01:01,length=3354>
##contig=<ID=HLA-C*07:02:01:03,length=3354>
##contig=<ID=HLA-C*07:04:01,length=3354>
##contig=<ID=HLA-C*07:385,length=3354>
##contig=<ID=HLA-C*07:67,length=3354>
##contig=<ID=HLA-C*07:391,length=3354>
##contig=<ID=HLA-C*07:01:19,length=3354>
##contig=<ID=HLA-C*07:02:06,length=3354>
##contig=<ID=HLA-C*07:01:45,length=3354>
##contig=<ID=HLA-C*07:56:02,length=3354>
##contig=<ID=HLA-C*07:392,length=3354>
##contig=<ID=HLA-C*07:02:64,length=3354>
##contig=<ID=HLA-C*07:01:01:01,length=3354>
##contig=<ID=HLA-C*07:18,length=3353>
##contig=<ID=HLA-C*07:02:01:04,length=3353>
##contig=<ID=HLA-C*07:01:02,length=3352>
##contig=<ID=chrUn_JTFH01000111v1_decoy,length=3351>
##contig=<ID=HLA-C*08:01:01,length=3349>
##contig=<ID=HLA-C*04:01:01:02,length=3349>
##contig=<ID=HLA-C*01:02:01,length=3349>
##contig=<ID=HLA-C*16:01:01,length=3349>
##contig=<ID=HLA-C*08:20,length=3349>
##contig=<ID=HLA-C*14:02:01,length=3349>
##contig=<ID=HLA-C*08:21,length=3349>
##contig=<ID=HLA-C*15:96Q,length=3349>
##contig=<ID=HLA-C*15:17,length=3349>
##contig=<ID=HLA-C*08:22,length=3349>
##contig=<ID=HLA-C*15:05:02,length=3349>
##contig=<ID=HLA-C*12:03:01:01,length=3349>
##contig=<ID=HLA-C*06:02:01:01,length=3349>
##contig=<ID=HLA-C*05:01:01:01,length=3349>
##contig=<ID=HLA-C*14:03,length=3349>
##contig=<ID=HLA-C*06:02:01:03,length=3349>
##contig=<ID=HLA-C*05:01:01:02,length=3349>
##contig=<ID=HLA-C*15:02:01,length=3349>
##contig=<ID=HLA-C*12:99,length=3349>
##contig=<ID=HLA-C*06:23,length=3349>
##contig=<ID=HLA-C*06:24,length=3349>
##contig=<ID=HLA-C*15:05:01,length=3349>
##contig=<ID=HLA-C*08:02:01:02,length=3349>
##contig=<ID=HLA-C*08:02:01:01,length=3349>
##contig=<ID=HLA-C*04:01:01:03,length=3349>
##contig=<ID=HLA-C*16:04:01,length=3349>
##contig=<ID=HLA-C*06:02:01:02,length=3349>
##contig=<ID=HLA-C*01:08,length=3349>
##contig=<ID=HLA-C*07:384,length=3349>
##contig=<ID=HLA-C*12:02:02,length=3349>
##contig=<ID=HLA-C*04:03:01,length=3349>
##contig=<ID=HLA-C*04:06,length=3349>
##contig=<ID=HLA-C*04:177,length=3349>
##contig=<ID=HLA-C*01:30,length=3349>
##contig=<ID=HLA-C*12:19,length=3349>
##contig=<ID=HLA-C*01:02:29,length=3349>
##contig=<ID=HLA-C*08:03:01,length=3349>
##contig=<ID=HLA-C*04:01:01:01,length=3349>
##contig=<ID=HLA-C*08:27,length=3349>
##contig=<ID=HLA-C*01:03,length=3349>
##contig=<ID=HLA-C*03:02:02:01,length=3348>
##contig=<ID=HLA-C*03:04:01:02,length=3348>
##contig=<ID=HLA-C*03:02:02:03,length=3348>
##contig=<ID=HLA-C*03:03:01,length=3348>
##contig=<ID=HLA-C*03:261,length=3348>
##contig=<ID=HLA-C*12:03:01:02,length=3348>
##contig=<ID=HLA-C*03:04:01:01,length=3348>
##contig=<ID=HLA-C*02:85,length=3347>
##contig=<ID=HLA-C*02:02:02:01,length=3347>
##contig=<ID=HLA-C*02:02:02:02,length=3347>
##contig=<ID=HLA-C*02:86,length=3347>
##contig=<ID=HLA-C*18:01,length=3346>
##contig=<ID=chrUn_KN707901v1_decoy,length=3344>
##contig=<ID=HLA-C*07:04:02,length=3343>
##contig=<ID=HLA-B*49:32,length=3340>
##contig=<ID=chrUn_JTFH01000112v1_decoy,length=3340>
##contig=<ID=HLA-B*49:01:01,length=3340>
##contig=<ID=HLA-B*50:01:01,length=3340>
##contig=<ID=HLA-B*45:04,length=3339>
##contig=<ID=HLA-B*45:01:01,length=3338>
##contig=<ID=HLA-B*57:29,length=3337>
##contig=<ID=chrUn_JTFH01001837v1_decoy,length=3337>
##contig=<ID=HLA-B*57:01:01,length=3337>
##contig=<ID=HLA-B*15:83,length=3337>
##contig=<ID=HLA-B*15:07:01,length=3336>
##contig=<ID=HLA-B*15:01:01:01,length=3336>
##contig=<ID=HLA-B*15:18:01,length=3336>
##contig=<ID=HLA-B*15:32:01,length=3336>
##contig=<ID=HLA-B*46:01:01,length=3336>
##contig=<ID=HLA-B*15:58,length=3336>
##contig=<ID=HLA-B*15:77,length=3336>
##contig=<ID=HLA-B*15:11:01,length=3336>
##contig=<ID=HLA-B*58:01:01,length=3336>
##contig=<ID=HLA-B*15:02:01,length=3335>
##contig=<ID=HLA-B*15:25:01,length=3335>
##contig=<ID=HLA-C*07:32N,length=3334>
##contig=<ID=HLA-C*01:02:30,length=3333>
##contig=<ID=HLA-B*59:01:01:01,length=3333>
##contig=<ID=HLA-B*55:02:01,length=3333>
##contig=<ID=HLA-B*15:42,length=3333>
##contig=<ID=HLA-B*54:01:01,length=3332>
##contig=<ID=HLA-B*55:01:03,length=3332>
##contig=<ID=HLA-B*55:01:01,length=3332>
##contig=<ID=HLA-B*59:01:01:02,length=3332>
##contig=<ID=HLA-B*55:24,length=3332>
##contig=<ID=HLA-B*55:12,length=3332>
##contig=<ID=HLA-C*04:01:62,length=3329>
##contig=<ID=HLA-C*03:41:02,length=3328>
##contig=<ID=HLA-B*35:01:01:02,length=3327>
##contig=<ID=HLA-B*52:01:01:01,length=3327>
##contig=<ID=HLA-B*52:01:01:02,length=3327>
##contig=<ID=HLA-B*78:01:01,length=3327>
##contig=<ID=HLA-B*52:01:01:03,length=3327>
##contig=<ID=HLA-B*51:07:01,length=3327>
##contig=<ID=HLA-B*52:01:02,length=3327>
##contig=<ID=HLA-B*53:01:01,length=3327>
##contig=<ID=HLA-B*51:01:01,length=3327>
##contig=<ID=HLA-B*35:14:02,length=3327>
##contig=<ID=HLA-B*51:02:01,length=3327>
##contig=<ID=HLA-B*35:01:01:01,length=3327>
##contig=<ID=HLA-B*35:41,length=3327>
##contig=<ID=HLA-B*35:02:01,length=3327>
##contig=<ID=HLA-B*40:06:01:01,length=3325>
##contig=<ID=HLA-B*27:04:01,length=3325>
##contig=<ID=HLA-B*27:05:02,length=3325>
##contig=<ID=HLA-B*27:32,length=3325>
##contig=<ID=HLA-B*27:06,length=3325>
##contig=<ID=HLA-B*27:131,length=3325>
##contig=<ID=HLA-B*13:01:01,length=3324>
##contig=<ID=HLA-B*13:02:01,length=3324>
##contig=<ID=HLA-B*37:01:01,length=3324>
##contig=<ID=HLA-B*13:08,length=3324>
##contig=<ID=chrUn_JTFH01001838v1_decoy,length=3324>
##contig=<ID=chrUn_KN707774v1_decoy,length=3324>
##contig=<ID=HLA-B*18:26,length=3323>
##contig=<ID=HLA-B*48:01:01,length=3323>
##contig=<ID=HLA-B*40:01:02,length=3323>
##contig=<ID=HLA-B*13:15,length=3323>
##contig=<ID=HLA-B*13:02:03,length=3323>
##contig=<ID=HLA-B*18:01:01:02,length=3323>
##contig=<ID=HLA-B*44:02:01:01,length=3323>
##contig=<ID=HLA-B*44:02:17,length=3323>
##contig=<ID=HLA-B*18:03,length=3323>
##contig=<ID=HLA-B*44:23N,length=3323>
##contig=<ID=HLA-B*44:03:01,length=3323>
##contig=<ID=HLA-B*73:01,length=3323>
##contig=<ID=HLA-B*18:01:01:01,length=3323>
##contig=<ID=HLA-B*07:02:01,length=3323>
##contig=<ID=HLA-B*48:08,length=3323>
##contig=<ID=HLA-B*07:50,length=3323>
##contig=<ID=HLA-B*44:46,length=3323>
##contig=<ID=HLA-B*08:33,length=3322>
##contig=<ID=HLA-C*05:09:01,length=3322>
##contig=<ID=HLA-B*42:01:01,length=3322>
##contig=<ID=HLA-B*41:01:01,length=3322>
##contig=<ID=HLA-B*08:20,length=3322>
##contig=<ID=HLA-B*08:19N,length=3322>
##contig=<ID=HLA-B*08:01:01,length=3322>
##contig=<ID=HLA-B*41:02:01,length=3322>
##contig=<ID=HLA-C*03:20N,length=3321>
##contig=<ID=HLA-B*27:05:18,length=3321>
##contig=<ID=HLA-C*02:11,length=3320>
##contig=<ID=chrUn_JTFH01000113v1_decoy,length=3320>
##contig=<ID=chrUn_JTFH01000975v1_decoy,length=3320>
##contig=<ID=HLA-B*44:09,length=3317>
##contig=<ID=chrUn_JTFH01001839v1_decoy,length=3315>
##contig=<ID=chrUn_JTFH01001840v1_decoy,length=3313>
##contig=<ID=HLA-B*38:01:01,length=3312>
##contig=<ID=HLA-B*39:01:01:03,length=3312>
##contig=<ID=HLA-B*38:02:01,length=3312>
##contig=<ID=HLA-B*39:01:21,length=3312>
##contig=<ID=HLA-B*14:02:01,length=3312>
##contig=<ID=HLA-B*14:01:01,length=3312>
##contig=<ID=HLA-B*67:01:01,length=3312>
##contig=<ID=HLA-A*29:46,length=3310>
##contig=<ID=HLA-A*02:81,length=3309>
##contig=<ID=HLA-B*67:02,length=3307>
##contig=<ID=HLA-B*57:11,length=3306>
##contig=<ID=HLA-A*02:53N,length=3305>
##contig=<ID=HLA-B*40:10:01,length=3304>
##contig=<ID=HLA-A*29:01:01:02N,length=3303>
##contig=<ID=HLA-B*40:06:01:02,length=3299>
##contig=<ID=chrUn_KN707840v1_decoy,length=3298>
##contig=<ID=HLA-A*01:01:01:02N,length=3291>
##contig=<ID=HLA-A*02:01:01:02L,length=3287>
##contig=<ID=chrUn_JTFH01001279v1_decoy,length=3285>
##contig=<ID=HLA-B*57:06,length=3284>
##contig=<ID=HLA-B*15:108,length=3283>
##contig=<ID=chrUn_JTFH01001841v1_decoy,length=3283>
##contig=<ID=HLA-B*40:72:01,length=3283>
##contig=<ID=chrUn_JTFH01000114v1_decoy,length=3282>
##contig=<ID=chrUn_JTFH01000115v1_decoy,length=3278>
##contig=<ID=chrUn_KN707795v1_decoy,length=3277>
##contig=<ID=chrUn_KN707916v1_decoy,length=3275>
##contig=<ID=HLA-B*53:11,length=3274>
##contig=<ID=chrUn_KN707709v1_decoy,length=3273>
##contig=<ID=chrUn_JTFH01001280v1_decoy,length=3273>
##contig=<ID=HLA-B*07:44,length=3270>
##contig=<ID=chrUn_KN707719v1_decoy,length=3270>
##contig=<ID=chrUn_KN707956v1_decoy,length=3270>
##contig=<ID=HLA-B*07:41,length=3266>
##contig=<ID=HLA-A*80:01:01:01,length=3263>
##contig=<ID=chrUn_JTFH01001281v1_decoy,length=3262>
##contig=<ID=chrUn_JTFH01000116v1_decoy,length=3260>
##contig=<ID=chrUn_JTFH01001282v1_decoy,length=3259>
##contig=<ID=chrUn_JTFH01000117v1_decoy,length=3258>
##contig=<ID=HLA-B*40:02:01,length=3258>
##contig=<ID=HLA-B*40:79,length=3257>
##contig=<ID=HLA-B*14:07N,length=3255>
##contig=<ID=HLA-B*39:13:02,length=3255>
##contig=<ID=HLA-B*39:34,length=3254>
##contig=<ID=chrUn_KI270517v1,length=3253>
##contig=<ID=chrUn_JTFH01000118v1_decoy,length=3253>
##contig=<ID=chrUn_JTFH01001842v1_decoy,length=3250>
##contig=<ID=chrUn_JTFH01001843v1_decoy,length=3247>
##contig=<ID=chrUn_JTFH01000119v1_decoy,length=3247>
##contig=<ID=HLA-A*24:02:03Q,length=3247>
##contig=<ID=HLA-A*11:60,length=3241>
##contig=<ID=chrUn_KN707791v1_decoy,length=3240>
##contig=<ID=HLA-B*44:04,length=3239>
##contig=<ID=HLA-B*07:33:01,length=3239>
##contig=<ID=HLA-C*04:161,length=3237>
##contig=<ID=HLA-A*68:18N,length=3237>
##contig=<ID=chrUn_JTFH01001844v1_decoy,length=3237>
##contig=<ID=chrUn_JTFH01001845v1_decoy,length=3235>
##contig=<ID=HLA-A*11:77,length=3233>
##contig=<ID=chrUn_JTFH01000976v1_decoy,length=3231>
##contig=<ID=chrUn_JTFH01000120v1_decoy,length=3230>
##contig=<ID=HLA-A*11:74,length=3227>
##contig=<ID=chrUn_JTFH01000121v1_decoy,length=3224>
##contig=<ID=chrUn_KN707784v1_decoy,length=3224>
##contig=<ID=chrUn_JTFH01001283v1_decoy,length=3222>
##contig=<ID=HLA-C*07:19,length=3222>
##contig=<ID=chrUn_JTFH01000977v1_decoy,length=3220>
##contig=<ID=chrUn_KN707855v1_decoy,length=3219>
##contig=<ID=HLA-A*02:43N,length=3218>
##contig=<ID=HLA-A*02:533,length=3217>
##contig=<ID=HLA-A*26:15,length=3217>
##contig=<ID=chrUn_JTFH01000122v1_decoy,length=3216>
##contig=<ID=chrUn_KN707836v1_decoy,length=3213>
##contig=<ID=chrUn_JTFH01000123v1_decoy,length=3212>
##contig=<ID=chrUn_JTFH01000978v1_decoy,length=3212>
##contig=<ID=chrUn_JTFH01001846v1_decoy,length=3200>
##contig=<ID=HLA-A*68:71,length=3198>
##contig=<ID=HLA-C*17:03,length=3197>
##contig=<ID=chrUn_JTFH01001847v1_decoy,length=3195>
##contig=<ID=HLA-C*07:01:27,length=3195>
##contig=<ID=chrUn_JTFH01000124v1_decoy,length=3194>
##contig=<ID=chrUn_JTFH01000979v1_decoy,length=3192>
##contig=<ID=chrUn_JTFH01000125v1_decoy,length=3189>
##contig=<ID=HLA-A*11:75,length=3184>
##contig=<ID=HLA-C*07:386,length=3183>
##contig=<ID=HLA-C*08:112,length=3178>
##contig=<ID=chrUn_JTFH01000126v1_decoy,length=3177>
##contig=<ID=HLA-A*24:152,length=3176>
##contig=<ID=chrUn_KN707771v1_decoy,length=3176>
##contig=<ID=chrUn_JTFH01000127v1_decoy,length=3176>
##contig=<ID=chrUn_JTFH01001848v1_decoy,length=3175>
##contig=<ID=chrUn_JTFH01000128v1_decoy,length=3173>
##contig=<ID=chrUn_JTFH01000129v1_decoy,length=3170>
##contig=<ID=chrUn_JTFH01000130v1_decoy,length=3166>
##contig=<ID=HLA-B*42:08,length=3165>
##contig=<ID=chrUn_JTFH01000131v1_decoy,length=3163>
##contig=<ID=chrUn_JTFH01001849v1_decoy,length=3158>
##contig=<ID=HLA-B*39:01:03,length=3155>
##contig=<ID=HLA-B*39:01:16,length=3155>
##contig=<ID=HLA-B*39:01:01:01,length=3155>
##contig=<ID=HLA-B*39:01:01:02L,length=3153>
##contig=<ID=HLA-B*44:02:01:02S,length=3152>
##contig=<ID=HLA-B*44:02:01:03,length=3152>
##contig=<ID=HLA-A*02:03:03,length=3148>
##contig=<ID=HLA-A*02:265,length=3148>
##contig=<ID=chrUn_KN707860v1_decoy,length=3148>
##contig=<ID=chrUn_JTFH01001850v1_decoy,length=3143>
##contig=<ID=chrUn_JTFH01000132v1_decoy,length=3143>
##contig=<ID=HLA-A*03:36N,length=3142>
##contig=<ID=HLA-A*26:50,length=3141>
##contig=<ID=chrUn_JTFH01001851v1_decoy,length=3139>
##contig=<ID=chrUn_JTFH01001852v1_decoy,length=3138>
##contig=<ID=chrUn_JTFH01000133v1_decoy,length=3137>
##contig=<ID=HLA-A*01:04N,length=3136>
##contig=<ID=chrUn_JTFH01001853v1_decoy,length=3136>
##contig=<ID=HLA-A*68:17,length=3134>
##contig=<ID=chrUn_KN707974v1_decoy,length=3134>
##contig=<ID=chrUn_JTFH01001855v1_decoy,length=3132>
##contig=<ID=chrUn_JTFH01001854v1_decoy,length=3132>
##contig=<ID=chrUn_JTFH01001284v1_decoy,length=3127>
##contig=<ID=chrUn_KN707614v1_decoy,length=3122>
##contig=<ID=HLA-A*68:08:01,length=3120>
##contig=<ID=HLA-A*68:22,length=3119>
##contig=<ID=chrUn_KN707695v1_decoy,length=3119>
##contig=<ID=HLA-A*02:455,length=3118>
##contig=<ID=chrUn_JTFH01000134v1_decoy,length=3116>
##contig=<ID=HLA-A*24:215,length=3116>
##contig=<ID=chrUn_KN707608v1_decoy,length=3112>
##contig=<ID=HLA-A*02:60:01,length=3112>
##contig=<ID=chrUn_KN707616v1_decoy,length=3111>
##contig=<ID=chrUn_JTFH01001285v1_decoy,length=3110>
##contig=<ID=HLA-A*02:51,length=3109>
##contig=<ID=HLA-A*02:68,length=3109>
##contig=<ID=chrUn_JTFH01000135v1_decoy,length=3106>
##contig=<ID=HLA-A*01:09,length=3105>
##contig=<ID=HLA-A*01:20,length=3105>
##contig=<ID=HLA-A*02:376,length=3104>
##contig=<ID=chrUn_JTFH01001286v1_decoy,length=3104>
##contig=<ID=HLA-A*23:09,length=3104>
##contig=<ID=HLA-A*02:279,length=3103>
##contig=<ID=HLA-A*02:269,length=3101>
##contig=<ID=HLA-C*14:21N,length=3099>
##contig=<ID=HLA-C*07:149,length=3098>
##contig=<ID=HLA-C*08:36N,length=3097>
##contig=<ID=HLA-A*34:02:01,length=3096>
##contig=<ID=HLA-A*03:21N,length=3095>
##contig=<ID=chrUn_JTFH01001856v1_decoy,length=3095>
##contig=<ID=HLA-A*01:14,length=3095>
##contig=<ID=chrUn_JTFH01001857v1_decoy,length=3094>
##contig=<ID=HLA-A*03:01:01:03,length=3094>
##contig=<ID=HLA-C*07:01:01:02,length=3093>
##contig=<ID=chrUn_JTFH01000136v1_decoy,length=3093>
##contig=<ID=chrUn_JTFH01001858v1_decoy,length=3093>
##contig=<ID=chrUn_JTFH01000980v1_decoy,length=3092>
##contig=<ID=HLA-A*26:11N,length=3091>
##contig=<ID=HLA-A*31:14N,length=3090>
##contig=<ID=chrUn_JTFH01000981v1_decoy,length=3087>
##contig=<ID=chrUn_KN707899v1_decoy,length=3087>
##contig=<ID=HLA-C*04:71,length=3086>
##contig=<ID=HLA-C*08:62,length=3086>
##contig=<ID=HLA-C*04:128,length=3086>
##contig=<ID=HLA-A*02:266,length=3084>
##contig=<ID=chrUn_KN707973v1_decoy,length=3080>
##contig=<ID=chrUn_JTFH01000137v1_decoy,length=3079>
##contig=<ID=HLA-A*31:46,length=3075>
##contig=<ID=HLA-A*24:02:01:03,length=3075>
##contig=<ID=HLA-A*66:17,length=3075>
##contig=<ID=HLA-C*07:02:01:02,length=3074>
##contig=<ID=HLA-A*11:25,length=3073>
##contig=<ID=chrUn_JTFH01001287v1_decoy,length=3071>
##contig=<ID=HLA-A*68:113,length=3070>
##contig=<ID=HLA-C*03:219,length=3070>
##contig=<ID=HLA-C*07:26,length=3069>
##contig=<ID=HLA-C*12:08,length=3066>
##contig=<ID=HLA-C*15:16,length=3066>
##contig=<ID=HLA-C*03:13:01,length=3065>
##contig=<ID=HLA-C*02:87,length=3064>
##contig=<ID=chrUn_JTFH01001288v1_decoy,length=3063>
##contig=<ID=chrUn_JTFH01001289v1_decoy,length=3059>
##contig=<ID=chrUn_JTFH01001859v1_decoy,length=3059>
##contig=<ID=HLA-C*05:08,length=3059>
##contig=<ID=HLA-C*12:13,length=3058>
##contig=<ID=HLA-C*04:70,length=3058>
##contig=<ID=HLA-C*01:02:11,length=3057>
##contig=<ID=HLA-A*80:01:01:02,length=3055>
##contig=<ID=HLA-A*02:57,length=3054>
##contig=<ID=chrUn_JTFH01000138v1_decoy,length=3053>
##contig=<ID=HLA-B*15:04:01,length=3052>
##contig=<ID=chrUn_JTFH01000139v1_decoy,length=3051>
##contig=<ID=HLA-B*15:17:01:01,length=3051>
##contig=<ID=HLA-B*15:17:01:02,length=3051>
##contig=<ID=HLA-B*82:02:01,length=3050>
##contig=<ID=chrUn_JTFH01000982v1_decoy,length=3048>
##contig=<ID=chrUn_KN707715v1_decoy,length=3044>
##contig=<ID=HLA-B*51:01:02,length=3043>
##contig=<ID=HLA-A*24:61,length=3043>
##contig=<ID=HLA-B*44:138Q,length=3043>
##contig=<ID=HLA-B*35:241,length=3042>
##contig=<ID=chrUn_KN707922v1_decoy,length=3041>
##contig=<ID=HLA-B*47:01:01:01,length=3041>
##contig=<ID=chrUn_KI270593v1,length=3041>
##contig=<ID=HLA-B*47:01:01:02,length=3041>
##contig=<ID=chrUn_KN707825v1_decoy,length=3040>
##contig=<ID=HLA-B*44:49,length=3039>
##contig=<ID=HLA-B*08:08N,length=3035>
##contig=<ID=chrUn_KN707757v1_decoy,length=3034>
##contig=<ID=HLA-C*03:100,length=3034>
##contig=<ID=chrUn_KN707607v1_decoy,length=3033>
##contig=<ID=HLA-C*02:16:02,length=3029>
##contig=<ID=HLA-B*39:10:01,length=3027>
##contig=<ID=HLA-B*15:01:01:03,length=3026>
##contig=<ID=HLA-A*02:01:01:03,length=3023>
##contig=<ID=chrUn_KN707808v1_decoy,length=3021>
##contig=<ID=HLA-A*23:38N,length=3020>
##contig=<ID=HLA-C*08:41,length=3019>
##contig=<ID=chrUn_JTFH01000140v1_decoy,length=3015>
##contig=<ID=chrUn_JTFH01000141v1_decoy,length=3012>
##contig=<ID=HLA-C*04:01:01:04,length=3012>
##contig=<ID=chrUn_JTFH01000142v1_decoy,length=3009>
##contig=<ID=chrUn_JTFH01000983v1_decoy,length=3005>
##contig=<ID=HLA-B*58:31N,length=3004>
##contig=<ID=chrUn_JTFH01000984v1_decoy,length=3004>
##contig=<ID=HLA-A*02:264,length=3002>
##contig=<ID=HLA-C*08:01:03,length=2998>
##contig=<ID=chrUn_JTFH01000143v1_decoy,length=2997>
##contig=<ID=chrUn_JTFH01000144v1_decoy,length=2997>
##contig=<ID=HLA-C*03:46,length=2997>
##contig=<ID=HLA-C*04:09N,length=2991>
##contig=<ID=chrUn_JTFH01001290v1_decoy,length=2990>
##contig=<ID=HLA-C*06:46N,length=2987>
##contig=<ID=chrUn_JTFH01001291v1_decoy,length=2986>
##contig=<ID=HLA-A*01:16N,length=2985>
##contig=<ID=chrUn_JTFH01001860v1_decoy,length=2985>
##contig=<ID=chrUn_KI270528v1,length=2983>
##contig=<ID=chrUn_JTFH01000145v1_decoy,length=2983>
##contig=<ID=HLA-B*55:48,length=2980>
##contig=<ID=HLA-B*18:17N,length=2979>
##contig=<ID=chrUn_JTFH01000146v1_decoy,length=2979>
##contig=<ID=HLA-A*02:259,length=2978>
##contig=<ID=HLA-C*08:40,length=2978>
##contig=<ID=HLA-C*14:23,length=2976>
##contig=<ID=chrUn_JTFH01001861v1_decoy,length=2975>
##contig=<ID=chrUn_KN707980v1_decoy,length=2973>
##contig=<ID=HLA-B*18:94N,length=2970>
##contig=<ID=chrUn_KI270587v1,length=2969>
##contig=<ID=HLA-C*01:40,length=2968>
##contig=<ID=chrUn_JTFH01001862v1_decoy,length=2967>
##contig=<ID=HLA-B*07:156,length=2967>
##contig=<ID=chrUn_JTFH01000147v1_decoy,length=2967>
##contig=<ID=chrUn_JTFH01000148v1_decoy,length=2967>
##contig=<ID=HLA-C*03:04:04,length=2966>
##contig=<ID=chrUn_JTFH01000149v1_decoy,length=2966>
##contig=<ID=HLA-B*51:42,length=2962>
##contig=<ID=chrUn_JTFH01001863v1_decoy,length=2961>
##contig=<ID=chrUn_KN707988v1_decoy,length=2960>
##contig=<ID=chrUn_KN707628v1_decoy,length=2960>
##contig=<ID=HLA-B*08:134,length=2959>
##contig=<ID=chrUn_KN707858v1_decoy,length=2959>
##contig=<ID=chrUn_JTFH01000985v1_decoy,length=2959>
##contig=<ID=chrUn_JTFH01001864v1_decoy,length=2955>
##contig=<ID=chrUn_JTFH01000150v1_decoy,length=2954>
##contig=<ID=chrUn_JTFH01000151v1_decoy,length=2952>
##contig=<ID=HLA-C*05:93,length=2946>
##contig=<ID=chrUn_KN707642v1_decoy,length=2943>
##contig=<ID=chrUn_KN707684v1_decoy,length=2940>
##contig=<ID=HLA-C*07:49,length=2935>
##contig=<ID=chrUn_JTFH01001865v1_decoy,length=2935>
##contig=<ID=chrUn_JTFH01000986v1_decoy,length=2934>
##contig=<ID=chrUn_JTFH01000152v1_decoy,length=2934>
##contig=<ID=chrUn_JTFH01000987v1_decoy,length=2933>
##contig=<ID=HLA-C*02:69,length=2933>
##contig=<ID=chrUn_JTFH01001866v1_decoy,length=2933>
##contig=<ID=HLA-C*04:01:01:05,length=2931>
##contig=<ID=HLA-A*68:01:01:01,length=2930>
##contig=<ID=chrUn_KN707985v1_decoy,length=2929>
##contig=<ID=chrUn_JTFH01001292v1_decoy,length=2928>
##contig=<ID=chrUn_KN707838v1_decoy,length=2926>
##contig=<ID=chrUn_JTFH01001293v1_decoy,length=2922>
##contig=<ID=chrUn_KN707714v1_decoy,length=2922>
##contig=<ID=chrUn_KN707902v1_decoy,length=2921>
##contig=<ID=HLA-B*13:02:09,length=2919>
##contig=<ID=HLA-A*74:01,length=2918>
##contig=<ID=HLA-A*31:04,length=2918>
##contig=<ID=HLA-A*74:02:01:01,length=2918>
##contig=<ID=chrUn_JTFH01000153v1_decoy,length=2918>
##contig=<ID=HLA-A*31:01:23,length=2918>
##contig=<ID=HLA-A*68:03:01,length=2917>
##contig=<ID=HLA-A*69:01,length=2917>
##contig=<ID=HLA-A*02:02:01,length=2917>
##contig=<ID=HLA-A*25:01:01,length=2917>
##contig=<ID=HLA-A*68:02:02,length=2916>
##contig=<ID=chrUn_KN707827v1_decoy,length=2913>
##contig=<ID=HLA-A*68:02:01:03,length=2909>
##contig=<ID=chrUn_JTFH01001867v1_decoy,length=2909>
##contig=<ID=chrUn_JTFH01001868v1_decoy,length=2904>
##contig=<ID=HLA-A*30:89,length=2903>
##contig=<ID=HLA-A*11:110,length=2903>
##contig=<ID=HLA-A*30:02:01:01,length=2903>
##contig=<ID=HLA-C*07:02:05,length=2903>
##contig=<ID=HLA-A*36:01,length=2903>
##contig=<ID=HLA-C*07:30,length=2903>
##contig=<ID=HLA-B*15:66,length=2902>
##contig=<ID=chrUn_KN707801v1_decoy,length=2901>
##contig=<ID=chrUn_KN707693v1_decoy,length=2899>
##contig=<ID=HLA-C*03:02:02:02,length=2896>
##contig=<ID=HLA-C*08:24,length=2895>
##contig=<ID=HLA-C*08:04:01,length=2895>
##contig=<ID=chrUn_JTFH01000154v1_decoy,length=2895>
##contig=<ID=HLA-C*01:21,length=2895>
##contig=<ID=HLA-C*16:02:01,length=2895>
##contig=<ID=HLA-C*01:06,length=2895>
##contig=<ID=HLA-C*01:14,length=2895>
##contig=<ID=HLA-C*15:13,length=2895>
##contig=<ID=HLA-C*12:22,length=2895>
##contig=<ID=HLA-C*03:06,length=2894>
##contig=<ID=HLA-C*03:05,length=2894>
##contig=<ID=HLA-C*03:02:01,length=2894>
##contig=<ID=HLA-C*03:61,length=2894>
##contig=<ID=HLA-C*03:40:01,length=2894>
##contig=<ID=HLA-C*02:10,length=2893>
##contig=<ID=chrUn_JTFH01001869v1_decoy,length=2892>
##contig=<ID=HLA-B*46:01:05,length=2891>
##contig=<ID=chrUn_KN707716v1_decoy,length=2888>
##contig=<ID=chrUn_JTFH01000155v1_decoy,length=2887>
##contig=<ID=chrUn_JTFH01001870v1_decoy,length=2886>
##contig=<ID=chrUn_JTFH01001871v1_decoy,length=2885>
##contig=<ID=chrUn_KN707722v1_decoy,length=2884>
##contig=<ID=chrUn_JTFH01000156v1_decoy,length=2879>
##contig=<ID=chrUn_JTFH01001872v1_decoy,length=2878>
##contig=<ID=chrUn_JTFH01000157v1_decoy,length=2878>
##contig=<ID=HLA-B*15:220,length=2878>
##contig=<ID=HLA-C*03:04:02,length=2877>
##contig=<ID=chrUn_JTFH01001294v1_decoy,length=2875>
##contig=<ID=chrUn_JTFH01001873v1_decoy,length=2875>
##contig=<ID=chrUn_KN707752v1_decoy,length=2873>
##contig=<ID=HLA-B*44:02:27,length=2872>
##contig=<ID=chrUn_JTFH01000158v1_decoy,length=2872>
##contig=<ID=chrUn_KN707704v1_decoy,length=2871>
##contig=<ID=chrUn_JTFH01000159v1_decoy,length=2868>
##contig=<ID=chrUn_JTFH01000160v1_decoy,length=2866>
##contig=<ID=chrUn_JTFH01000161v1_decoy,length=2865>
##contig=<ID=chrUn_JTFH01000162v1_decoy,length=2864>
##contig=<ID=chrUn_JTFH01001874v1_decoy,length=2861>
##contig=<ID=chrUn_JTFH01000163v1_decoy,length=2859>
##contig=<ID=chrUn_JTFH01001295v1_decoy,length=2859>
##contig=<ID=chrUn_KN707643v1_decoy,length=2857>
##contig=<ID=chrUn_JTFH01001875v1_decoy,length=2856>
##contig=<ID=chrUn_KI270364v1,length=2855>
##contig=<ID=chrUn_JTFH01000164v1_decoy,length=2854>
##contig=<ID=chrUn_JTFH01001296v1_decoy,length=2850>
##contig=<ID=chrUn_KN707674v1_decoy,length=2848>
##contig=<ID=chrUn_JTFH01001876v1_decoy,length=2838>
##contig=<ID=chrUn_JTFH01000165v1_decoy,length=2830>
##contig=<ID=chrUn_JTFH01000166v1_decoy,length=2828>
##contig=<ID=chrUn_JTFH01000988v1_decoy,length=2827>
##contig=<ID=chrUn_KN707758v1_decoy,length=2826>
##contig=<ID=chrUn_JTFH01000167v1_decoy,length=2824>
##contig=<ID=chrUn_JTFH01000168v1_decoy,length=2819>
##contig=<ID=HLA-B*54:18,length=2813>
##contig=<ID=chrUn_JTFH01001297v1_decoy,length=2813>
##contig=<ID=chrUn_JTFH01000169v1_decoy,length=2813>
##contig=<ID=chrUn_JTFH01000170v1_decoy,length=2809>
##contig=<ID=HLA-B*35:01:22,length=2806>
##contig=<ID=chrUn_KI270371v1,length=2805>
##contig=<ID=HLA-B*44:26,length=2804>
##contig=<ID=chrUn_JTFH01000171v1_decoy,length=2802>
##contig=<ID=chrUn_JTFH01001877v1_decoy,length=2801>
##contig=<ID=HLA-B*40:150,length=2800>
##contig=<ID=chrUn_JTFH01001878v1_decoy,length=2797>
##contig=<ID=chrUn_KN707699v1_decoy,length=2795>
##contig=<ID=chrUn_JTFH01000989v1_decoy,length=2794>
##contig=<ID=chrUn_KN707672v1_decoy,length=2794>
##contig=<ID=chrUn_JTFH01000172v1_decoy,length=2791>
##contig=<ID=chrUn_KN707849v1_decoy,length=2790>
##contig=<ID=chrUn_JTFH01001879v1_decoy,length=2788>
##contig=<ID=chrUn_KN707671v1_decoy,length=2786>
##contig=<ID=chrUn_JTFH01001298v1_decoy,length=2785>
##contig=<ID=chrUn_KN707655v1_decoy,length=2785>
##contig=<ID=chrUn_JTFH01000173v1_decoy,length=2783>
##contig=<ID=chrUn_JTFH01000174v1_decoy,length=2778>
##contig=<ID=chrUn_JTFH01000175v1_decoy,length=2777>
##contig=<ID=chrUn_JTFH01001880v1_decoy,length=2773>
##contig=<ID=chrUn_KN707882v1_decoy,length=2772>
##contig=<ID=chrUn_JTFH01000176v1_decoy,length=2770>
##contig=<ID=chrUn_JTFH01000177v1_decoy,length=2769>
##contig=<ID=chrUn_JTFH01000178v1_decoy,length=2767>
##contig=<ID=HLA-B*39:14,length=2765>
##contig=<ID=chrUn_JTFH01000179v1_decoy,length=2763>
##contig=<ID=chrUn_JTFH01001881v1_decoy,length=2755>
##contig=<ID=chrUn_JTFH01001882v1_decoy,length=2754>
##contig=<ID=chrUn_JTFH01000990v1_decoy,length=2749>
##contig=<ID=chrUn_JTFH01000991v1_decoy,length=2745>
##contig=<ID=chrUn_JTFH01000180v1_decoy,length=2745>
##contig=<ID=chrUn_JTFH01001883v1_decoy,length=2743>
##contig=<ID=chrUn_JTFH01000181v1_decoy,length=2742>
##contig=<ID=HLA-B*38:14,length=2738>
##contig=<ID=chrUn_JTFH01001299v1_decoy,length=2736>
##contig=<ID=chrUn_JTFH01000182v1_decoy,length=2736>
##contig=<ID=chrUn_JTFH01000992v1_decoy,length=2733>
##contig=<ID=chrUn_JTFH01000183v1_decoy,length=2729>
##contig=<ID=chrUn_JTFH01000184v1_decoy,length=2726>
##contig=<ID=chrUn_JTFH01001884v1_decoy,length=2725>
##contig=<ID=chrUn_JTFH01001885v1_decoy,length=2722>
##contig=<ID=chrUn_JTFH01000185v1_decoy,length=2719>
##contig=<ID=chrUn_KN707781v1_decoy,length=2717>
##contig=<ID=chrUn_JTFH01000186v1_decoy,length=2715>
##contig=<ID=chrUn_JTFH01000187v1_decoy,length=2708>
##contig=<ID=chrUn_JTFH01000188v1_decoy,length=2704>
##contig=<ID=chrUn_KI270333v1,length=2699>
##contig=<ID=chrUn_KN707886v1_decoy,length=2699>
##contig=<ID=chrUn_JTFH01000993v1_decoy,length=2698>
##contig=<ID=chrUn_JTFH01000189v1_decoy,length=2692>
##contig=<ID=chrUn_JTFH01000190v1_decoy,length=2691>
##contig=<ID=HLA-B*35:05:01,length=2690>
##contig=<ID=chrUn_JTFH01000191v1_decoy,length=2690>
##contig=<ID=HLA-B*15:03:01,length=2689>
##contig=<ID=HLA-B*57:03:01,length=2689>
##contig=<ID=HLA-B*15:10:01,length=2689>
##contig=<ID=HLA-B*15:27:01,length=2689>
##contig=<ID=HLA-B*13:25,length=2689>
##contig=<ID=HLA-B*35:03:01,length=2689>
##contig=<ID=HLA-B*35:08:01,length=2689>
##contig=<ID=HLA-B*56:01:01,length=2688>
##contig=<ID=HLA-B*15:13:01,length=2688>
##contig=<ID=chrUn_JTFH01001300v1_decoy,length=2688>
##contig=<ID=HLA-B*56:03,length=2688>
##contig=<ID=HLA-B*56:04,length=2688>
##contig=<ID=HLA-B*15:16:01,length=2688>
##contig=<ID=HLA-B*37:01:05,length=2687>
##contig=<ID=chrUn_JTFH01000192v1_decoy,length=2687>
##contig=<ID=HLA-B*18:02,length=2686>
##contig=<ID=chrUn_JTFH01001886v1_decoy,length=2682>
##contig=<ID=HLA-B*40:03,length=2677>
##contig=<ID=HLA-B*27:07:01,length=2677>
##contig=<ID=HLA-B*27:25,length=2677>
##contig=<ID=chrUn_JTFH01000193v1_decoy,length=2677>
##contig=<ID=HLA-B*40:40,length=2677>
##contig=<ID=HLA-B*27:24,length=2677>
##contig=<ID=HLA-B*07:06,length=2676>
##contig=<ID=HLA-B*44:03:02,length=2676>
##contig=<ID=HLA-B*08:79,length=2676>
##contig=<ID=HLA-B*81:01,length=2676>
##contig=<ID=HLA-B*48:04,length=2676>
##contig=<ID=HLA-B*44:56N,length=2676>
##contig=<ID=HLA-B*48:03:01,length=2676>
##contig=<ID=HLA-B*07:05:01,length=2676>
##contig=<ID=HLA-B*44:150,length=2676>
##contig=<ID=HLA-B*40:01:01,length=2676>
##contig=<ID=HLA-B*39:05:01,length=2675>
##contig=<ID=HLA-B*39:38Q,length=2675>
##contig=<ID=HLA-B*42:02,length=2675>
##contig=<ID=HLA-B*08:132,length=2675>
##contig=<ID=HLA-B*67:01:02,length=2675>
##contig=<ID=HLA-B*39:06:02,length=2674>
##contig=<ID=chrUn_KN707763v1_decoy,length=2674>
##contig=<ID=chrUn_KN707665v1_decoy,length=2670>
##contig=<ID=chrUn_JTFH01001887v1_decoy,length=2669>
##contig=<ID=chrUn_JTFH01000195v1_decoy,length=2668>
##contig=<ID=chrUn_JTFH01000194v1_decoy,length=2668>
##contig=<ID=chrUn_KN707802v1_decoy,length=2666>
##contig=<ID=chrUn_JTFH01000994v1_decoy,length=2665>
##contig=<ID=chrUn_JTFH01001888v1_decoy,length=2663>
##contig=<ID=chrUn_JTFH01000196v1_decoy,length=2663>
##contig=<ID=chrUn_JTFH01001301v1_decoy,length=2658>
##contig=<ID=chrUn_KI270374v1,length=2656>
##contig=<ID=chrUn_JTFH01000197v1_decoy,length=2655>
##contig=<ID=chrUn_JTFH01001889v1_decoy,length=2652>
##contig=<ID=chrUn_KN707979v1_decoy,length=2648>
##contig=<ID=chrUn_JTFH01001890v1_decoy,length=2647>
##contig=<ID=chrUn_KI270411v1,length=2646>
##contig=<ID=chrUn_JTFH01000198v1_decoy,length=2644>
##contig=<ID=chrUn_JTFH01001302v1_decoy,length=2643>
##contig=<ID=chrUn_JTFH01000199v1_decoy,length=2642>
##contig=<ID=chrUn_JTFH01001891v1_decoy,length=2635>
##contig=<ID=chrUn_JTFH01000995v1_decoy,length=2634>
##contig=<ID=chrUn_JTFH01001892v1_decoy,length=2633>
##contig=<ID=chrUn_JTFH01000200v1_decoy,length=2632>
##contig=<ID=chrUn_JTFH01000201v1_decoy,length=2632>
##contig=<ID=chrUn_KN707785v1_decoy,length=2631>
##contig=<ID=chrUn_KN707727v1_decoy,length=2630>
##contig=<ID=chrUn_JTFH01001893v1_decoy,length=2629>
##contig=<ID=chrUn_JTFH01000202v1_decoy,length=2628>
##contig=<ID=chrUn_KN707918v1_decoy,length=2623>
##contig=<ID=chrUn_JTFH01000203v1_decoy,length=2623>
##contig=<ID=chrUn_JTFH01000204v1_decoy,length=2622>
##contig=<ID=chrUn_JTFH01000205v1_decoy,length=2619>
##contig=<ID=chrUn_JTFH01001303v1_decoy,length=2618>
##contig=<ID=chrUn_KN707776v1_decoy,length=2618>
##contig=<ID=chrUn_JTFH01001894v1_decoy,length=2612>
##contig=<ID=chrUn_KN707983v1_decoy,length=2606>
##contig=<ID=chrUn_JTFH01001304v1_decoy,length=2605>
##contig=<ID=chrUn_JTFH01000206v1_decoy,length=2605>
##contig=<ID=chrUn_KN707659v1_decoy,length=2605>
##contig=<ID=chrUn_JTFH01000207v1_decoy,length=2603>
##contig=<ID=chrUn_JTFH01000208v1_decoy,length=2601>
##contig=<ID=chrUn_JTFH01001895v1_decoy,length=2599>
##contig=<ID=chrUn_JTFH01000209v1_decoy,length=2598>
##contig=<ID=chrUn_JTFH01000210v1_decoy,length=2597>
##contig=<ID=chrUn_JTFH01000211v1_decoy,length=2596>
##contig=<ID=chrUn_JTFH01000212v1_decoy,length=2594>
##contig=<ID=chrUn_JTFH01000213v1_decoy,length=2586>
##contig=<ID=chrUn_JTFH01000214v1_decoy,length=2585>
##contig=<ID=chrUn_JTFH01000215v1_decoy,length=2583>
##contig=<ID=chrUn_JTFH01001305v1_decoy,length=2583>
##contig=<ID=chrUn_JTFH01000216v1_decoy,length=2578>
##contig=<ID=chrUn_KN707822v1_decoy,length=2575>
##contig=<ID=chrUn_JTFH01000217v1_decoy,length=2569>
##contig=<ID=chrUn_JTFH01000218v1_decoy,length=2569>
##contig=<ID=chrUn_JTFH01001896v1_decoy,length=2566>
##contig=<ID=chrUn_KN707725v1_decoy,length=2565>
##contig=<ID=chrUn_KN707921v1_decoy,length=2561>
##contig=<ID=chrUn_JTFH01001897v1_decoy,length=2556>
##contig=<ID=chrUn_KN707767v1_decoy,length=2552>
##contig=<ID=chrUn_JTFH01001898v1_decoy,length=2551>
##contig=<ID=chrUn_JTFH01001899v1_decoy,length=2551>
##contig=<ID=chrUn_JTFH01000219v1_decoy,length=2551>
##contig=<ID=chrUn_JTFH01000221v1_decoy,length=2548>
##contig=<ID=chrUn_JTFH01000220v1_decoy,length=2548>
##contig=<ID=chrUn_JTFH01000222v1_decoy,length=2546>
##contig=<ID=chrUn_KN707617v1_decoy,length=2545>
##contig=<ID=chrUn_JTFH01000223v1_decoy,length=2545>
##contig=<ID=chrUn_JTFH01001900v1_decoy,length=2538>
##contig=<ID=chrUn_JTFH01001901v1_decoy,length=2538>
##contig=<ID=chrUn_JTFH01001306v1_decoy,length=2534>
##contig=<ID=chrUn_JTFH01000224v1_decoy,length=2534>
##contig=<ID=chrUn_JTFH01000225v1_decoy,length=2533>
##contig=<ID=chrUn_KN707892v1_decoy,length=2530>
##contig=<ID=chrUn_KN707915v1_decoy,length=2527>
##contig=<ID=chrUn_JTFH01001902v1_decoy,length=2525>
##contig=<ID=chrUn_JTFH01000227v1_decoy,length=2522>
##contig=<ID=chrUn_JTFH01000226v1_decoy,length=2522>
##contig=<ID=chrUn_KN707981v1_decoy,length=2520>
##contig=<ID=chrUn_JTFH01000228v1_decoy,length=2515>
##contig=<ID=chrUn_JTFH01000229v1_decoy,length=2513>
##contig=<ID=chrUn_JTFH01001307v1_decoy,length=2512>
##contig=<ID=chrUn_JTFH01000230v1_decoy,length=2507>
##contig=<ID=chrUn_JTFH01000231v1_decoy,length=2504>
##contig=<ID=chrUn_JTFH01001308v1_decoy,length=2500>
##contig=<ID=chrUn_JTFH01001903v1_decoy,length=2498>
##contig=<ID=chrUn_JTFH01000232v1_decoy,length=2497>
##contig=<ID=chrUn_JTFH01001904v1_decoy,length=2496>
##contig=<ID=chrUn_JTFH01000996v1_decoy,length=2492>
##contig=<ID=chrUn_KN707883v1_decoy,length=2490>
##contig=<ID=chrUn_KI270414v1,length=2489>
##contig=<ID=chrUn_JTFH01000997v1_decoy,length=2489>
##contig=<ID=chrUn_JTFH01001905v1_decoy,length=2483>
##contig=<ID=chrUn_JTFH01001309v1_decoy,length=2481>
##contig=<ID=chrUn_JTFH01001310v1_decoy,length=2478>
##contig=<ID=chrUn_JTFH01001906v1_decoy,length=2475>
##contig=<ID=chrUn_KN707958v1_decoy,length=2474>
##contig=<ID=chrUn_JTFH01001311v1_decoy,length=2473>
##contig=<ID=chrUn_JTFH01000233v1_decoy,length=2471>
##contig=<ID=chrUn_JTFH01001907v1_decoy,length=2469>
##contig=<ID=chrUn_JTFH01000998v1_decoy,length=2468>
##contig=<ID=chrUn_JTFH01001312v1_decoy,length=2467>
##contig=<ID=chrUn_KN707736v1_decoy,length=2467>
##contig=<ID=chrUn_JTFH01000234v1_decoy,length=2465>
##contig=<ID=chrUn_JTFH01000235v1_decoy,length=2464>
##contig=<ID=chrUn_KN707946v1_decoy,length=2463>
##contig=<ID=chrUn_KN707678v1_decoy,length=2462>
##contig=<ID=chrUn_JTFH01000236v1_decoy,length=2459>
##contig=<ID=chrUn_JTFH01000237v1_decoy,length=2457>
##contig=<ID=chrUn_JTFH01001908v1_decoy,length=2455>
##contig=<ID=chrUn_JTFH01000238v1_decoy,length=2450>
##contig=<ID=chrUn_JTFH01001909v1_decoy,length=2444>
##contig=<ID=chrUn_JTFH01001313v1_decoy,length=2442>
##contig=<ID=chrUn_KN707778v1_decoy,length=2440>
##contig=<ID=chrUn_JTFH01001910v1_decoy,length=2437>
##contig=<ID=chrUn_JTFH01000239v1_decoy,length=2435>
##contig=<ID=chrUn_JTFH01001911v1_decoy,length=2435>
##contig=<ID=chrUn_JTFH01000240v1_decoy,length=2434>
##contig=<ID=chrUn_JTFH01000241v1_decoy,length=2432>
##contig=<ID=chrUn_JTFH01001314v1_decoy,length=2430>
##contig=<ID=chrUn_JTFH01001912v1_decoy,length=2427>
##contig=<ID=chrUn_JTFH01000242v1_decoy,length=2427>
##contig=<ID=chrUn_KN707888v1_decoy,length=2424>
##contig=<ID=chrUn_JTFH01000243v1_decoy,length=2421>
##contig=<ID=chrUn_KN707666v1_decoy,length=2420>
##contig=<ID=chrUn_JTFH01000244v1_decoy,length=2420>
##contig=<ID=chrUn_JTFH01001913v1_decoy,length=2419>
##contig=<ID=chrUn_JTFH01001315v1_decoy,length=2417>
##contig=<ID=chrUn_KI270510v1,length=2415>
##contig=<ID=chrUn_JTFH01000245v1_decoy,length=2414>
##contig=<ID=chrUn_JTFH01000999v1_decoy,length=2414>
##contig=<ID=chrUn_JTFH01001914v1_decoy,length=2413>
##contig=<ID=chrUn_JTFH01001915v1_decoy,length=2412>
##contig=<ID=chrUn_JTFH01001316v1_decoy,length=2408>
##contig=<ID=chrUn_JTFH01000246v1_decoy,length=2404>
##contig=<ID=chrUn_JTFH01000247v1_decoy,length=2403>
##contig=<ID=chrUn_JTFH01000248v1_decoy,length=2402>
##contig=<ID=chrUn_JTFH01001916v1_decoy,length=2400>
##contig=<ID=chrUn_JTFH01001917v1_decoy,length=2399>
##contig=<ID=chrUn_JTFH01000249v1_decoy,length=2397>
##contig=<ID=chrUn_JTFH01001918v1_decoy,length=2396>
##contig=<ID=chrUn_JTFH01000250v1_decoy,length=2395>
##contig=<ID=chrUn_JTFH01001000v1_decoy,length=2395>
##contig=<ID=chrUn_JTFH01001317v1_decoy,length=2395>
##contig=<ID=chrUn_JTFH01000251v1_decoy,length=2394>
##contig=<ID=chrUn_JTFH01001919v1_decoy,length=2393>
##contig=<ID=chrUn_KN707828v1_decoy,length=2389>
##contig=<ID=chrUn_JTFH01000252v1_decoy,length=2388>
##contig=<ID=chrUn_KI270390v1,length=2387>
##contig=<ID=chrUn_JTFH01001920v1_decoy,length=2386>
##contig=<ID=chrUn_JTFH01001921v1_decoy,length=2384>
##contig=<ID=chrUn_JTFH01001922v1_decoy,length=2382>
##contig=<ID=chrUn_JTFH01001923v1_decoy,length=2382>
##contig=<ID=chrUn_JTFH01000253v1_decoy,length=2382>
##contig=<ID=chrUn_JTFH01000254v1_decoy,length=2381>
##contig=<ID=chrUn_JTFH01000255v1_decoy,length=2380>
##contig=<ID=chrUn_KI270375v1,length=2378>
##contig=<ID=chrUn_JTFH01000256v1_decoy,length=2368>
##contig=<ID=chrUn_JTFH01001924v1_decoy,length=2367>
##contig=<ID=chrUn_JTFH01001925v1_decoy,length=2366>
##contig=<ID=chrUn_JTFH01000257v1_decoy,length=2364>
##contig=<ID=chrUn_JTFH01000258v1_decoy,length=2363>
##contig=<ID=chrUn_JTFH01001926v1_decoy,length=2362>
##contig=<ID=chrUn_JTFH01001927v1_decoy,length=2361>
##contig=<ID=chrUn_JTFH01001001v1_decoy,length=2356>
##contig=<ID=chrUn_JTFH01001928v1_decoy,length=2353>
##contig=<ID=chrUn_JTFH01001318v1_decoy,length=2352>
##contig=<ID=chrUn_JTFH01001929v1_decoy,length=2349>
##contig=<ID=chrUn_JTFH01000259v1_decoy,length=2348>
##contig=<ID=chrUn_JTFH01001930v1_decoy,length=2348>
##contig=<ID=chrUn_JTFH01001931v1_decoy,length=2340>
##contig=<ID=chrUn_KN707905v1_decoy,length=2339>
##contig=<ID=chrUn_JTFH01001002v1_decoy,length=2339>
##contig=<ID=chrUn_JTFH01000260v1_decoy,length=2339>
##contig=<ID=chrUn_JTFH01001932v1_decoy,length=2339>
##contig=<ID=chrUn_JTFH01001319v1_decoy,length=2337>
##contig=<ID=chrUn_JTFH01001933v1_decoy,length=2336>
##contig=<ID=chrUn_JTFH01000261v1_decoy,length=2335>
##contig=<ID=chrUn_JTFH01001934v1_decoy,length=2333>
##contig=<ID=chrUn_JTFH01000262v1_decoy,length=2332>
##contig=<ID=chrUn_JTFH01000263v1_decoy,length=2331>
##contig=<ID=chrUn_JTFH01001935v1_decoy,length=2330>
##contig=<ID=chrUn_JTFH01000264v1_decoy,length=2330>
##contig=<ID=chrUn_JTFH01001936v1_decoy,length=2327>
##contig=<ID=chrUn_JTFH01000265v1_decoy,length=2323>
##contig=<ID=chrUn_JTFH01001320v1_decoy,length=2322>
##contig=<ID=chrUn_KI270420v1,length=2321>
##contig=<ID=chrUn_KN707761v1_decoy,length=2319>
##contig=<ID=chrUn_JTFH01000266v1_decoy,length=2319>
##contig=<ID=chrUn_JTFH01001937v1_decoy,length=2318>
##contig=<ID=chrUn_KI270509v1,length=2318>
##contig=<ID=chrUn_KN707982v1_decoy,length=2318>
##contig=<ID=chrUn_KN707630v1_decoy,length=2315>
##contig=<ID=chrUn_JTFH01000267v1_decoy,length=2314>
##contig=<ID=chrUn_JTFH01001003v1_decoy,length=2310>
##contig=<ID=chrUn_JTFH01000268v1_decoy,length=2308>
##contig=<ID=chrUn_JTFH01001321v1_decoy,length=2307>
##contig=<ID=chrUn_JTFH01000269v1_decoy,length=2306>
##contig=<ID=chrUn_JTFH01001322v1_decoy,length=2306>
##contig=<ID=chrUn_KN707708v1_decoy,length=2305>
##contig=<ID=chrUn_KN707766v1_decoy,length=2303>
##contig=<ID=chrUn_KN707856v1_decoy,length=2300>
##contig=<ID=chrUn_KN707846v1_decoy,length=2297>
##contig=<ID=chrUn_JTFH01000270v1_decoy,length=2296>
##contig=<ID=chrUn_KN707618v1_decoy,length=2295>
##contig=<ID=chrUn_KN707959v1_decoy,length=2294>
##contig=<ID=chrUn_JTFH01001938v1_decoy,length=2293>
##contig=<ID=chrUn_JTFH01001323v1_decoy,length=2292>
##contig=<ID=chrUn_JTFH01001939v1_decoy,length=2292>
##contig=<ID=chrUn_KN707953v1_decoy,length=2289>
##contig=<ID=chrUn_JTFH01001004v1_decoy,length=2288>
##contig=<ID=chrUn_JTFH01000271v1_decoy,length=2287>
##contig=<ID=chrUn_JTFH01001940v1_decoy,length=2287>
##contig=<ID=chrUn_JTFH01001005v1_decoy,length=2285>
##contig=<ID=chrUn_JTFH01000272v1_decoy,length=2279>
##contig=<ID=chrUn_JTFH01000273v1_decoy,length=2276>
##contig=<ID=chrUn_KN707821v1_decoy,length=2276>
##contig=<ID=chrUn_KI270315v1,length=2276>
##contig=<ID=chrUn_JTFH01001941v1_decoy,length=2274>
##contig=<ID=chrUn_KI270302v1,length=2274>
##contig=<ID=chrUn_JTFH01001942v1_decoy,length=2274>
##contig=<ID=chrUn_JTFH01000274v1_decoy,length=2273>
##contig=<ID=chrUn_KN707877v1_decoy,length=2272>
##contig=<ID=chrUn_KN707813v1_decoy,length=2272>
##contig=<ID=chrUn_JTFH01001324v1_decoy,length=2271>
##contig=<ID=chrUn_JTFH01001006v1_decoy,length=2269>
##contig=<ID=chrUn_JTFH01001943v1_decoy,length=2267>
##contig=<ID=chrUn_JTFH01001325v1_decoy,length=2265>
##contig=<ID=chrUn_JTFH01000275v1_decoy,length=2262>
##contig=<ID=chrUn_JTFH01001944v1_decoy,length=2260>
##contig=<ID=chrUn_JTFH01001326v1_decoy,length=2260>
##contig=<ID=chrUn_JTFH01001945v1_decoy,length=2257>
##contig=<ID=chrUn_JTFH01000276v1_decoy,length=2254>
##contig=<ID=chrUn_JTFH01001007v1_decoy,length=2253>
##contig=<ID=chrUn_JTFH01000277v1_decoy,length=2252>
##contig=<ID=chrUn_JTFH01000278v1_decoy,length=2245>
##contig=<ID=chrUn_KN707754v1_decoy,length=2243>
##contig=<ID=chrUn_KN707861v1_decoy,length=2242>
##contig=<ID=chrUn_JTFH01001327v1_decoy,length=2240>
##contig=<ID=chrUn_JTFH01001946v1_decoy,length=2240>
##contig=<ID=chrUn_JTFH01001947v1_decoy,length=2239>
##contig=<ID=chrUn_JTFH01000279v1_decoy,length=2239>
##contig=<ID=chrUn_JTFH01001328v1_decoy,length=2238>
##contig=<ID=chrUn_KN707800v1_decoy,length=2237>
##contig=<ID=chrUn_KN707703v1_decoy,length=2235>
##contig=<ID=chrUn_JTFH01001948v1_decoy,length=2232>
##contig=<ID=chrUn_JTFH01001949v1_decoy,length=2230>
##contig=<ID=chrUn_JTFH01001950v1_decoy,length=2230>
##contig=<ID=chrUn_JTFH01001329v1_decoy,length=2228>
##contig=<ID=chrUn_JTFH01000280v1_decoy,length=2223>
##contig=<ID=chrUn_JTFH01001951v1_decoy,length=2222>
##contig=<ID=chrUn_JTFH01000281v1_decoy,length=2220>
##contig=<ID=chrUn_JTFH01000282v1_decoy,length=2218>
##contig=<ID=chrUn_JTFH01001952v1_decoy,length=2216>
##contig=<ID=chrUn_JTFH01000283v1_decoy,length=2215>
##contig=<ID=chrUn_JTFH01001330v1_decoy,length=2215>
##contig=<ID=chrUn_JTFH01001953v1_decoy,length=2214>
##contig=<ID=chrUn_JTFH01000284v1_decoy,length=2213>
##contig=<ID=chrUn_JTFH01001954v1_decoy,length=2210>
##contig=<ID=chrUn_JTFH01001331v1_decoy,length=2205>
##contig=<ID=chrUn_KN707984v1_decoy,length=2205>
##contig=<ID=chrUn_JTFH01001955v1_decoy,length=2203>
##contig=<ID=chrUn_KN707955v1_decoy,length=2203>
##contig=<ID=chrUn_JTFH01000285v1_decoy,length=2203>
##contig=<ID=chrUn_JTFH01001008v1_decoy,length=2203>
##contig=<ID=chrUn_KN707831v1_decoy,length=2201>
##contig=<ID=chrUn_KN707606v1_decoy,length=2200>
##contig=<ID=chrUn_JTFH01000286v1_decoy,length=2200>
##contig=<ID=chrUn_JTFH01001956v1_decoy,length=2197>
##contig=<ID=chrUn_JTFH01000287v1_decoy,length=2197>
##contig=<ID=chrUn_JTFH01001958v1_decoy,length=2196>
##contig=<ID=chrUn_JTFH01001957v1_decoy,length=2196>
##contig=<ID=chrUn_JTFH01000288v1_decoy,length=2194>
##contig=<ID=chrUn_KN707991v1_decoy,length=2193>
##contig=<ID=chrUn_JTFH01001333v1_decoy,length=2191>
##contig=<ID=chrUn_JTFH01001332v1_decoy,length=2191>
##contig=<ID=chrUn_JTFH01001334v1_decoy,length=2190>
##contig=<ID=chrUn_KN707667v1_decoy,length=2189>
##contig=<ID=chrUn_KN707638v1_decoy,length=2189>
##contig=<ID=chrUn_KI270518v1,length=2186>
##contig=<ID=chrUn_JTFH01001335v1_decoy,length=2184>
##contig=<ID=chrUn_JTFH01000289v1_decoy,length=2183>
##contig=<ID=chrUn_JTFH01001959v1_decoy,length=2179>
##contig=<ID=chrUn_JTFH01000290v1_decoy,length=2179>
##contig=<ID=chrUn_JTFH01001961v1_decoy,length=2178>
##contig=<ID=chrUn_JTFH01001960v1_decoy,length=2178>
##contig=<ID=chrUn_JTFH01000293v1_decoy,length=2177>
##contig=<ID=chrUn_JTFH01000291v1_decoy,length=2177>
##contig=<ID=chrUn_JTFH01000292v1_decoy,length=2177>
##contig=<ID=chrUn_JTFH01001009v1_decoy,length=2176>
##contig=<ID=chrUn_KN707662v1_decoy,length=2173>
##contig=<ID=chrUn_JTFH01001962v1_decoy,length=2172>
##contig=<ID=chrUn_JTFH01001963v1_decoy,length=2170>
##contig=<ID=chrUn_JTFH01000294v1_decoy,length=2168>
##contig=<ID=chrUn_KI270530v1,length=2168>
##contig=<ID=chrUn_JTFH01001964v1_decoy,length=2167>
##contig=<ID=chrUn_JTFH01001965v1_decoy,length=2167>
##contig=<ID=chrUn_JTFH01001336v1_decoy,length=2166>
##contig=<ID=chrUn_JTFH01001337v1_decoy,length=2165>
##contig=<ID=chrUn_KN707872v1_decoy,length=2165>
##contig=<ID=chrUn_KI270304v1,length=2165>
##contig=<ID=chrUn_JTFH01001338v1_decoy,length=2162>
##contig=<ID=chrUn_JTFH01000295v1_decoy,length=2160>
##contig=<ID=chrUn_JTFH01001010v1_decoy,length=2159>
##contig=<ID=chrUn_KN707864v1_decoy,length=2159>
##contig=<ID=chrUn_JTFH01001966v1_decoy,length=2157>
##contig=<ID=chrUn_JTFH01001011v1_decoy,length=2155>
##contig=<ID=chrUn_JTFH01000296v1_decoy,length=2155>
##contig=<ID=chrUn_JTFH01001967v1_decoy,length=2153>
##contig=<ID=chrUn_JTFH01001968v1_decoy,length=2151>
##contig=<ID=chrUn_JTFH01001012v1_decoy,length=2149>
##contig=<ID=chrUn_KN707726v1_decoy,length=2149>
##contig=<ID=chrUn_JTFH01001969v1_decoy,length=2147>
##contig=<ID=chrUn_JTFH01001339v1_decoy,length=2146>
##contig=<ID=chrUn_JTFH01001970v1_decoy,length=2145>
##contig=<ID=chrUn_KI270418v1,length=2145>
##contig=<ID=chrUn_KN707753v1_decoy,length=2144>
##contig=<ID=chrUn_JTFH01000297v1_decoy,length=2144>
##contig=<ID=chrUn_JTFH01000298v1_decoy,length=2143>
##contig=<ID=chrUn_JTFH01001971v1_decoy,length=2142>
##contig=<ID=chrUn_JTFH01001972v1_decoy,length=2142>
##contig=<ID=chrUn_KI270424v1,length=2140>
##contig=<ID=chrUn_KN707897v1_decoy,length=2137>
##contig=<ID=chrUn_JTFH01001973v1_decoy,length=2136>
##contig=<ID=chrUn_JTFH01000299v1_decoy,length=2136>
##contig=<ID=chrUn_JTFH01000300v1_decoy,length=2134>
##contig=<ID=chrUn_JTFH01001974v1_decoy,length=2130>
##contig=<ID=chrUn_JTFH01001013v1_decoy,length=2129>
##contig=<ID=chrUn_JTFH01000301v1_decoy,length=2129>
##contig=<ID=chrUn_JTFH01001975v1_decoy,length=2128>
##contig=<ID=chrUn_JTFH01000302v1_decoy,length=2128>
##contig=<ID=chrUn_JTFH01001977v1_decoy,length=2126>
##contig=<ID=chrUn_JTFH01001976v1_decoy,length=2126>
##contig=<ID=chrUn_JTFH01000304v1_decoy,length=2125>
##contig=<ID=chrUn_JTFH01000303v1_decoy,length=2125>
##contig=<ID=chrUn_JTFH01000305v1_decoy,length=2122>
##contig=<ID=chrUn_JTFH01001978v1_decoy,length=2119>
##contig=<ID=chrUn_JTFH01001014v1_decoy,length=2116>
##contig=<ID=chrUn_JTFH01001340v1_decoy,length=2116>
##contig=<ID=chrUn_KN707865v1_decoy,length=2114>
##contig=<ID=chrUn_JTFH01001015v1_decoy,length=2113>
##contig=<ID=chrUn_JTFH01001341v1_decoy,length=2112>
##contig=<ID=chrUn_JTFH01000306v1_decoy,length=2111>
##contig=<ID=chrUn_JTFH01001342v1_decoy,length=2108>
##contig=<ID=chrUn_JTFH01001979v1_decoy,length=2107>
##contig=<ID=chrUn_JTFH01001343v1_decoy,length=2106>
##contig=<ID=chrUn_JTFH01000307v1_decoy,length=2106>
##contig=<ID=chrUn_JTFH01001344v1_decoy,length=2106>
##contig=<ID=chrUn_JTFH01001345v1_decoy,length=2106>
##contig=<ID=chrUn_JTFH01001016v1_decoy,length=2098>
##contig=<ID=chrUn_JTFH01001346v1_decoy,length=2097>
##contig=<ID=chrUn_JTFH01000308v1_decoy,length=2094>
##contig=<ID=chrUn_KN707668v1_decoy,length=2093>
##contig=<ID=chrUn_JTFH01000309v1_decoy,length=2093>
##contig=<ID=chrUn_JTFH01001980v1_decoy,length=2091>
##contig=<ID=chrUn_JTFH01000310v1_decoy,length=2088>
##contig=<ID=chrUn_JTFH01001981v1_decoy,length=2087>
##contig=<ID=chrUn_JTFH01000311v1_decoy,length=2086>
##contig=<ID=chrUn_JTFH01000312v1_decoy,length=2086>
##contig=<ID=chrUn_JTFH01001982v1_decoy,length=2086>
##contig=<ID=chrUn_KN707878v1_decoy,length=2085>
##contig=<ID=chrUn_JTFH01000313v1_decoy,length=2084>
##contig=<ID=chrUn_JTFH01001983v1_decoy,length=2083>
##contig=<ID=chrUn_JTFH01001347v1_decoy,length=2081>
##contig=<ID=chrUn_JTFH01000314v1_decoy,length=2080>
##contig=<ID=chrUn_JTFH01000315v1_decoy,length=2079>
##contig=<ID=chrUn_JTFH01000316v1_decoy,length=2076>
##contig=<ID=chrUn_KN707747v1_decoy,length=2075>
##contig=<ID=chrUn_JTFH01001984v1_decoy,length=2075>
##contig=<ID=chrUn_JTFH01001985v1_decoy,length=2075>
##contig=<ID=chrUn_JTFH01001986v1_decoy,length=2072>
##contig=<ID=chrUn_KN707686v1_decoy,length=2072>
##contig=<ID=chrUn_JTFH01000317v1_decoy,length=2071>
##contig=<ID=chrUn_KN707826v1_decoy,length=2070>
##contig=<ID=chrUn_JTFH01001987v1_decoy,length=2068>
##contig=<ID=chrUn_JTFH01001988v1_decoy,length=2067>
##contig=<ID=chrUn_JTFH01000318v1_decoy,length=2066>
##contig=<ID=chrUn_JTFH01001017v1_decoy,length=2066>
##contig=<ID=chrUn_JTFH01001018v1_decoy,length=2066>
##contig=<ID=chrUn_JTFH01000319v1_decoy,length=2061>
##contig=<ID=chrUn_JTFH01001019v1_decoy,length=2059>
##contig=<ID=chrUn_JTFH01001348v1_decoy,length=2058>
##contig=<ID=chrUn_JTFH01001349v1_decoy,length=2055>
##contig=<ID=chrUn_KN707914v1_decoy,length=2055>
##contig=<ID=chrUn_JTFH01000320v1_decoy,length=2055>
##contig=<ID=chrUn_JTFH01001989v1_decoy,length=2055>
##contig=<ID=chrUn_JTFH01001350v1_decoy,length=2054>
##contig=<ID=chrUn_JTFH01000321v1_decoy,length=2053>
##contig=<ID=chrUn_JTFH01001990v1_decoy,length=2051>
##contig=<ID=chrUn_JTFH01001991v1_decoy,length=2050>
##contig=<ID=chrUn_JTFH01001020v1_decoy,length=2047>
##contig=<ID=chrUn_KN707620v1_decoy,length=2046>
##contig=<ID=chrUn_KI270417v1,length=2043>
##contig=<ID=chrUn_KN707944v1_decoy,length=2043>
##contig=<ID=chrUn_KN707900v1_decoy,length=2041>
##contig=<ID=chrUn_JTFH01000322v1_decoy,length=2040>
##contig=<ID=chrUn_JTFH01001021v1_decoy,length=2040>
##contig=<ID=chrUn_KN707933v1_decoy,length=2038>
##contig=<ID=chrUn_JTFH01001351v1_decoy,length=2037>
##contig=<ID=chrUn_JTFH01000323v1_decoy,length=2036>
##contig=<ID=chrUn_JTFH01000324v1_decoy,length=2035>
##contig=<ID=chrUn_JTFH01000325v1_decoy,length=2034>
##contig=<ID=chrUn_JTFH01001992v1_decoy,length=2033>
##contig=<ID=chrUn_JTFH01001352v1_decoy,length=2032>
##contig=<ID=chrUn_JTFH01001353v1_decoy,length=2032>
##contig=<ID=chrUn_JTFH01000326v1_decoy,length=2032>
##contig=<ID=chrUn_JTFH01001022v1_decoy,length=2030>
##contig=<ID=chrUn_JTFH01000327v1_decoy,length=2029>
##contig=<ID=chrUn_JTFH01000328v1_decoy,length=2025>
##contig=<ID=chrUn_JTFH01001023v1_decoy,length=2024>
##contig=<ID=chrUn_JTFH01001993v1_decoy,length=2024>
##contig=<ID=chrUn_JTFH01000329v1_decoy,length=2021>
##contig=<ID=chrUn_JTFH01001354v1_decoy,length=2020>
##contig=<ID=chrUn_JTFH01001355v1_decoy,length=2018>
##contig=<ID=chrUn_JTFH01000330v1_decoy,length=2018>
##contig=<ID=chrUn_JTFH01001994v1_decoy,length=2016>
##contig=<ID=chrUn_JTFH01000331v1_decoy,length=2015>
##contig=<ID=chrUn_JTFH01001356v1_decoy,length=2014>
##contig=<ID=chrUn_KN707651v1_decoy,length=2013>
##contig=<ID=chrUn_JTFH01001995v1_decoy,length=2011>
##contig=<ID=chrUn_JTFH01000332v1_decoy,length=2009>
##contig=<ID=chrUn_JTFH01001996v1_decoy,length=2009>
##contig=<ID=chrUn_JTFH01000333v1_decoy,length=2007>
##contig=<ID=chrUn_KN707850v1_decoy,length=2006>
##contig=<ID=chrUn_JTFH01000334v1_decoy,length=2005>
##contig=<ID=chrUn_JTFH01001997v1_decoy,length=2003>
##contig=<ID=chrUn_JTFH01000335v1_decoy,length=2003>
##contig=<ID=chrUn_JTFH01000336v1_decoy,length=2001>
##contig=<ID=chrUn_JTFH01000337v1_decoy,length=2001>
##contig=<ID=chrUn_JTFH01001358v1_decoy,length=2001>
##contig=<ID=chrUn_JTFH01001357v1_decoy,length=2001>
##contig=<ID=chrUn_JTFH01001998v1_decoy,length=2001>
##contig=<ID=chrUn_JTFH01001024v1_decoy,length=2001>
##contig=<ID=chrUn_JTFH01000338v1_decoy,length=2000>
##contig=<ID=chrUn_JTFH01000339v1_decoy,length=1996>
##contig=<ID=chrUn_JTFH01000340v1_decoy,length=1992>
##contig=<ID=chrUn_JTFH01001025v1_decoy,length=1992>
##contig=<ID=chrUn_JTFH01001359v1_decoy,length=1991>
##contig=<ID=chrUn_JTFH01001360v1_decoy,length=1990>
##contig=<ID=chrUn_JTFH01000341v1_decoy,length=1985>
##contig=<ID=chrUn_JTFH01001361v1_decoy,length=1983>
##contig=<ID=chrUn_KN707647v1_decoy,length=1982>
##contig=<ID=chrUn_JTFH01000342v1_decoy,length=1981>
##contig=<ID=chrUn_JTFH01001362v1_decoy,length=1981>
##contig=<ID=chrUn_JTFH01001026v1_decoy,length=1981>
##contig=<ID=chrUn_JTFH01001363v1_decoy,length=1981>
##contig=<ID=chrUn_JTFH01001027v1_decoy,length=1979>
##contig=<ID=chrUn_JTFH01001364v1_decoy,length=1979>
##contig=<ID=chrUn_JTFH01000343v1_decoy,length=1977>
##contig=<ID=chrUn_JTFH01000344v1_decoy,length=1971>
##contig=<ID=chrUn_JTFH01000345v1_decoy,length=1968>
##contig=<ID=chrUn_JTFH01001365v1_decoy,length=1963>
##contig=<ID=chrUn_JTFH01000346v1_decoy,length=1962>
##contig=<ID=chrUn_JTFH01000347v1_decoy,length=1961>
##contig=<ID=chrUn_KN707833v1_decoy,length=1961>
##contig=<ID=chrUn_JTFH01000349v1_decoy,length=1960>
##contig=<ID=chrUn_JTFH01000348v1_decoy,length=1960>
##contig=<ID=chrUn_JTFH01001028v1_decoy,length=1957>
##contig=<ID=chrUn_JTFH01000350v1_decoy,length=1954>
##contig=<ID=chrUn_JTFH01001029v1_decoy,length=1953>
##contig=<ID=chrUn_JTFH01000351v1_decoy,length=1952>
##contig=<ID=chrUn_KI270508v1,length=1951>
##contig=<ID=chrUn_JTFH01000352v1_decoy,length=1947>
##contig=<ID=chrUn_KN707631v1_decoy,length=1945>
##contig=<ID=chrUn_JTFH01000353v1_decoy,length=1944>
##contig=<ID=chrUn_JTFH01001030v1_decoy,length=1944>
##contig=<ID=chrUn_JTFH01000354v1_decoy,length=1943>
##contig=<ID=chrUn_KI270303v1,length=1942>
##contig=<ID=chrUn_JTFH01000355v1_decoy,length=1941>
##contig=<ID=chrUn_JTFH01000356v1_decoy,length=1937>
##contig=<ID=chrUn_JTFH01001031v1_decoy,length=1936>
##contig=<ID=chrUn_KN707615v1_decoy,length=1934>
##contig=<ID=chrUn_JTFH01000357v1_decoy,length=1934>
##contig=<ID=chrUn_JTFH01001366v1_decoy,length=1932>
##contig=<ID=chrUn_KN707835v1_decoy,length=1932>
##contig=<ID=chrUn_JTFH01001032v1_decoy,length=1932>
##contig=<ID=chrUn_KI270381v1,length=1930>
##contig=<ID=chrUn_JTFH01001367v1_decoy,length=1929>
##contig=<ID=chrUn_JTFH01000358v1_decoy,length=1929>
##contig=<ID=chrUn_JTFH01000359v1_decoy,length=1924>
##contig=<ID=chrUn_JTFH01000360v1_decoy,length=1924>
##contig=<ID=chrUn_JTFH01000361v1_decoy,length=1923>
##contig=<ID=chrUn_JTFH01000362v1_decoy,length=1921>
##contig=<ID=chrUn_JTFH01000363v1_decoy,length=1918>
##contig=<ID=chrUn_JTFH01000364v1_decoy,length=1915>
##contig=<ID=chrUn_JTFH01000365v1_decoy,length=1915>
##contig=<ID=chrUn_KN707792v1_decoy,length=1915>
##contig=<ID=chrUn_JTFH01000366v1_decoy,length=1914>
##contig=<ID=chrUn_JTFH01000367v1_decoy,length=1912>
##contig=<ID=chrUn_JTFH01000368v1_decoy,length=1910>
##contig=<ID=chrUn_KN707698v1_decoy,length=1908>
##contig=<ID=chrUn_JTFH01000369v1_decoy,length=1907>
##contig=<ID=chrUn_JTFH01000370v1_decoy,length=1904>
##contig=<ID=chrUn_KI270529v1,length=1899>
##contig=<ID=chrUn_JTFH01000371v1_decoy,length=1897>
##contig=<ID=chrUn_KN707911v1_decoy,length=1893>
##contig=<ID=chrUn_JTFH01000372v1_decoy,length=1891>
##contig=<ID=chrUn_JTFH01000373v1_decoy,length=1890>
##contig=<ID=chrUn_KN707653v1_decoy,length=1890>
##contig=<ID=chrUn_JTFH01000374v1_decoy,length=1888>
##contig=<ID=chrUn_JTFH01000375v1_decoy,length=1888>
##contig=<ID=chrUn_JTFH01000376v1_decoy,length=1885>
##contig=<ID=chrUn_KI270425v1,length=1884>
##contig=<ID=chrUn_JTFH01001033v1_decoy,length=1882>
##contig=<ID=chrUn_KN707940v1_decoy,length=1882>
##contig=<ID=chrUn_JTFH01000377v1_decoy,length=1881>
##contig=<ID=chrUn_JTFH01001368v1_decoy,length=1881>
##contig=<ID=chrUn_KI270396v1,length=1880>
##contig=<ID=chrUn_JTFH01000378v1_decoy,length=1879>
##contig=<ID=chrUn_JTFH01001034v1_decoy,length=1878>
##contig=<ID=chrUn_JTFH01000379v1_decoy,length=1877>
##contig=<ID=chrUn_JTFH01000381v1_decoy,length=1876>
##contig=<ID=chrUn_JTFH01000380v1_decoy,length=1876>
##contig=<ID=chrUn_JTFH01001369v1_decoy,length=1874>
##contig=<ID=chrUn_JTFH01000382v1_decoy,length=1874>
##contig=<ID=chrUn_JTFH01000383v1_decoy,length=1872>
##contig=<ID=chrUn_JTFH01001035v1_decoy,length=1870>
##contig=<ID=chrUn_JTFH01000384v1_decoy,length=1869>
##contig=<ID=chrUn_JTFH01000385v1_decoy,length=1866>
##contig=<ID=chrUn_JTFH01000386v1_decoy,length=1865>
##contig=<ID=chrUn_JTFH01000388v1_decoy,length=1865>
##contig=<ID=chrUn_JTFH01000387v1_decoy,length=1865>
##contig=<ID=chrUn_JTFH01000390v1_decoy,length=1862>
##contig=<ID=chrUn_JTFH01000389v1_decoy,length=1862>
##contig=<ID=chrUn_JTFH01000391v1_decoy,length=1859>
##contig=<ID=chrUn_JTFH01000393v1_decoy,length=1856>
##contig=<ID=chrUn_JTFH01000392v1_decoy,length=1856>
##contig=<ID=chrUn_JTFH01000394v1_decoy,length=1854>
##contig=<ID=chrUn_JTFH01000395v1_decoy,length=1850>
##contig=<ID=chrUn_JTFH01001371v1_decoy,length=1849>
##contig=<ID=chrUn_JTFH01000396v1_decoy,length=1849>
##contig=<ID=chrUn_JTFH01000397v1_decoy,length=1849>
##contig=<ID=chrUn_JTFH01001370v1_decoy,length=1849>
##contig=<ID=chrUn_KN707629v1_decoy,length=1848>
##contig=<ID=chrUn_JTFH01000398v1_decoy,length=1847>
##contig=<ID=chrUn_JTFH01000399v1_decoy,length=1839>
##contig=<ID=chrUn_KN707829v1_decoy,length=1835>
##contig=<ID=chrUn_JTFH01000400v1_decoy,length=1834>
##contig=<ID=chrUn_JTFH01001372v1_decoy,length=1833>
##contig=<ID=chrUn_JTFH01001373v1_decoy,length=1832>
##contig=<ID=chrUn_KN707640v1_decoy,length=1831>
##contig=<ID=chrUn_KN707992v1_decoy,length=1830>
##contig=<ID=chrUn_JTFH01001374v1_decoy,length=1826>
##contig=<ID=chrUn_JTFH01000401v1_decoy,length=1821>
##contig=<ID=chrUn_JTFH01001036v1_decoy,length=1821>
##contig=<ID=chrUn_JTFH01000402v1_decoy,length=1815>
##contig=<ID=chrUn_JTFH01001375v1_decoy,length=1814>
##contig=<ID=chrUn_JTFH01001376v1_decoy,length=1814>
##contig=<ID=chrUn_JTFH01001037v1_decoy,length=1813>
##contig=<ID=chrUn_JTFH01000403v1_decoy,length=1811>
##contig=<ID=chrUn_JTFH01001038v1_decoy,length=1809>
##contig=<ID=chrUn_JTFH01000404v1_decoy,length=1808>
##contig=<ID=chrUn_JTFH01000405v1_decoy,length=1808>
##contig=<ID=chrUn_JTFH01000407v1_decoy,length=1807>
##contig=<ID=chrUn_JTFH01000406v1_decoy,length=1807>
##contig=<ID=chrUn_JTFH01001039v1_decoy,length=1804>
##contig=<ID=chrUn_KI270363v1,length=1803>
##contig=<ID=chrUn_JTFH01000408v1_decoy,length=1802>
##contig=<ID=chrUn_JTFH01000409v1_decoy,length=1801>
##contig=<ID=chrUn_JTFH01000410v1_decoy,length=1800>
##contig=<ID=chrUn_JTFH01001040v1_decoy,length=1797>
##contig=<ID=chrUn_JTFH01000411v1_decoy,length=1795>
##contig=<ID=chrUn_JTFH01000412v1_decoy,length=1794>
##contig=<ID=chrUn_JTFH01000413v1_decoy,length=1792>
##contig=<ID=chrUn_JTFH01001377v1_decoy,length=1791>
##contig=<ID=chrUn_JTFH01001041v1_decoy,length=1791>
##contig=<ID=chrUn_JTFH01001378v1_decoy,length=1789>
##contig=<ID=chrUn_KI270386v1,length=1788>
##contig=<ID=chrUn_JTFH01000414v1_decoy,length=1788>
##contig=<ID=chrUn_JTFH01000415v1_decoy,length=1786>
##contig=<ID=chrUn_JTFH01001379v1_decoy,length=1786>
##contig=<ID=chrUn_JTFH01000417v1_decoy,length=1782>
##contig=<ID=chrUn_JTFH01000416v1_decoy,length=1782>
##contig=<ID=chrUn_JTFH01000419v1_decoy,length=1781>
##contig=<ID=chrUn_JTFH01000418v1_decoy,length=1781>
##contig=<ID=chrUn_JTFH01001042v1_decoy,length=1781>
##contig=<ID=chrUn_JTFH01000420v1_decoy,length=1779>
##contig=<ID=chrUn_JTFH01001380v1_decoy,length=1778>
##contig=<ID=chrUn_JTFH01000421v1_decoy,length=1777>
##contig=<ID=chrUn_JTFH01001381v1_decoy,length=1776>
##contig=<ID=chrUn_KN707885v1_decoy,length=1776>
##contig=<ID=chrUn_KN707649v1_decoy,length=1775>
##contig=<ID=chrUn_KI270465v1,length=1774>
##contig=<ID=chrUn_KN707679v1_decoy,length=1774>
##contig=<ID=chrUn_KN707938v1_decoy,length=1770>
##contig=<ID=chrUn_KN707844v1_decoy,length=1766>
##contig=<ID=chrUn_JTFH01001043v1_decoy,length=1766>
##contig=<ID=chrUn_JTFH01000422v1_decoy,length=1764>
##contig=<ID=chrUn_JTFH01001044v1_decoy,length=1764>
##contig=<ID=chrUn_JTFH01000423v1_decoy,length=1762>
##contig=<ID=chrUn_KN707817v1_decoy,length=1762>
##contig=<ID=chrUn_JTFH01001382v1_decoy,length=1762>
##contig=<ID=chrUn_JTFH01001383v1_decoy,length=1758>
##contig=<ID=chrUn_JTFH01001384v1_decoy,length=1757>
##contig=<ID=chrUn_JTFH01000424v1_decoy,length=1755>
##contig=<ID=chrUn_JTFH01001385v1_decoy,length=1754>
##contig=<ID=chrUn_JTFH01001386v1_decoy,length=1752>
##contig=<ID=chrUn_JTFH01001387v1_decoy,length=1751>
##contig=<ID=chrUn_KI270383v1,length=1750>
##contig=<ID=chrUn_JTFH01000425v1_decoy,length=1749>
##contig=<ID=chrUn_JTFH01001388v1_decoy,length=1749>
##contig=<ID=chrUn_KN707889v1_decoy,length=1747>
##contig=<ID=chrUn_JTFH01000426v1_decoy,length=1747>
##contig=<ID=chrUn_JTFH01000427v1_decoy,length=1746>
##contig=<ID=chrUn_JTFH01000428v1_decoy,length=1745>
##contig=<ID=chrUn_JTFH01000429v1_decoy,length=1744>
##contig=<ID=chrUn_JTFH01001045v1_decoy,length=1743>
##contig=<ID=chrUn_JTFH01000430v1_decoy,length=1742>
##contig=<ID=chrUn_KN707717v1_decoy,length=1742>
##contig=<ID=chrUn_JTFH01001046v1_decoy,length=1741>
##contig=<ID=chrUn_JTFH01000431v1_decoy,length=1740>
##contig=<ID=chrUn_JTFH01000432v1_decoy,length=1740>
##contig=<ID=chrUn_JTFH01001389v1_decoy,length=1738>
##contig=<ID=chrUn_KN707909v1_decoy,length=1737>
##contig=<ID=chrUn_JTFH01000433v1_decoy,length=1736>
##contig=<ID=chrUn_JTFH01000434v1_decoy,length=1735>
##contig=<ID=chrUn_KN707646v1_decoy,length=1735>
##contig=<ID=chrUn_JTFH01000435v1_decoy,length=1732>
##contig=<ID=chrUn_JTFH01000436v1_decoy,length=1732>
##contig=<ID=chrUn_JTFH01000437v1_decoy,length=1730>
##contig=<ID=chrUn_JTFH01001390v1_decoy,length=1729>
##contig=<ID=chrUn_JTFH01000438v1_decoy,length=1727>
##contig=<ID=chrUn_JTFH01001391v1_decoy,length=1726>
##contig=<ID=chrUn_KN707636v1_decoy,length=1725>
##contig=<ID=chrUn_JTFH01000439v1_decoy,length=1722>
##contig=<ID=chrUn_JTFH01000440v1_decoy,length=1718>
##contig=<ID=chrUn_JTFH01001392v1_decoy,length=1716>
##contig=<ID=chrUn_JTFH01000441v1_decoy,length=1716>
##contig=<ID=chrUn_JTFH01001393v1_decoy,length=1712>
##contig=<ID=chrUn_JTFH01001394v1_decoy,length=1711>
##contig=<ID=chrUn_JTFH01000442v1_decoy,length=1710>
##contig=<ID=chrUn_JTFH01001047v1_decoy,length=1709>
##contig=<ID=chrUn_JTFH01000443v1_decoy,length=1708>
##contig=<ID=chrUn_JTFH01000444v1_decoy,length=1707>
##contig=<ID=chrUn_JTFH01001048v1_decoy,length=1706>
##contig=<ID=chrUn_JTFH01000445v1_decoy,length=1706>
##contig=<ID=chrUn_JTFH01000446v1_decoy,length=1705>
##contig=<ID=chrUn_JTFH01000447v1_decoy,length=1704>
##contig=<ID=chrUn_JTFH01001395v1_decoy,length=1703>
##contig=<ID=chrUn_JTFH01001396v1_decoy,length=1702>
##contig=<ID=chrUn_JTFH01001049v1_decoy,length=1701>
##contig=<ID=chrUn_JTFH01000448v1_decoy,length=1699>
##contig=<ID=chrUn_JTFH01001397v1_decoy,length=1699>
##contig=<ID=chrUn_JTFH01000449v1_decoy,length=1698>
##contig=<ID=chrUn_JTFH01000450v1_decoy,length=1697>
##contig=<ID=chrUn_JTFH01000451v1_decoy,length=1697>
##contig=<ID=chrUn_JTFH01000452v1_decoy,length=1695>
##contig=<ID=chrUn_JTFH01000453v1_decoy,length=1695>
##contig=<ID=chrUn_JTFH01000454v1_decoy,length=1693>
##contig=<ID=chrUn_JTFH01001050v1_decoy,length=1689>
##contig=<ID=chrUn_JTFH01000455v1_decoy,length=1687>
##contig=<ID=chrUn_JTFH01001398v1_decoy,length=1686>
##contig=<ID=chrUn_JTFH01000456v1_decoy,length=1686>
##contig=<ID=chrUn_JTFH01001399v1_decoy,length=1684>
##contig=<ID=chrUn_JTFH01001400v1_decoy,length=1680>
##contig=<ID=chrUn_JTFH01000457v1_decoy,length=1680>
##contig=<ID=chrUn_JTFH01000459v1_decoy,length=1679>
##contig=<ID=chrUn_JTFH01000458v1_decoy,length=1679>
##contig=<ID=chrUn_JTFH01001401v1_decoy,length=1678>
##contig=<ID=chrUn_JTFH01001402v1_decoy,length=1678>
##contig=<ID=chrUn_JTFH01000460v1_decoy,length=1678>
##contig=<ID=chrUn_JTFH01001403v1_decoy,length=1677>
##contig=<ID=chrUn_JTFH01001404v1_decoy,length=1676>
##contig=<ID=chrUn_JTFH01000461v1_decoy,length=1674>
##contig=<ID=chrUn_JTFH01000462v1_decoy,length=1674>
##contig=<ID=chrUn_JTFH01001405v1_decoy,length=1672>
##contig=<ID=chrUn_JTFH01000463v1_decoy,length=1671>
##contig=<ID=chrUn_JTFH01001406v1_decoy,length=1669>
##contig=<ID=chrUn_JTFH01000464v1_decoy,length=1669>
##contig=<ID=chrUn_JTFH01001407v1_decoy,length=1668>
##contig=<ID=chrUn_JTFH01000465v1_decoy,length=1665>
##contig=<ID=chrUn_JTFH01000466v1_decoy,length=1663>
##contig=<ID=chrUn_JTFH01001408v1_decoy,length=1663>
##contig=<ID=chrUn_JTFH01001410v1_decoy,length=1660>
##contig=<ID=chrUn_JTFH01001409v1_decoy,length=1660>
##contig=<ID=chrUn_KI270384v1,length=1658>
##contig=<ID=chrUn_JTFH01001411v1_decoy,length=1658>
##contig=<ID=chrUn_JTFH01000467v1_decoy,length=1657>
##contig=<ID=chrUn_JTFH01001413v1_decoy,length=1656>
##contig=<ID=chrUn_JTFH01001412v1_decoy,length=1656>
##contig=<ID=chrUn_JTFH01000468v1_decoy,length=1653>
##contig=<ID=chrUn_JTFH01001414v1_decoy,length=1652>
##contig=<ID=chrUn_KI270330v1,length=1652>
##contig=<ID=chrUn_JTFH01000469v1_decoy,length=1652>
##contig=<ID=chrUn_JTFH01000470v1_decoy,length=1650>
##contig=<ID=chrUn_KI270372v1,length=1650>
##contig=<ID=chrUn_JTFH01000472v1_decoy,length=1649>
##contig=<ID=chrUn_JTFH01000471v1_decoy,length=1649>
##contig=<ID=chrUn_KN707954v1_decoy,length=1648>
##contig=<ID=chrUn_KN707641v1_decoy,length=1647>
##contig=<ID=chrUn_JTFH01001415v1_decoy,length=1647>
##contig=<ID=chrUn_JTFH01001051v1_decoy,length=1646>
##contig=<ID=chrUn_JTFH01001416v1_decoy,length=1645>
##contig=<ID=chrUn_KN707609v1_decoy,length=1642>
##contig=<ID=chrUn_JTFH01001052v1_decoy,length=1641>
##contig=<ID=chrUn_JTFH01001417v1_decoy,length=1641>
##contig=<ID=chrUn_JTFH01000473v1_decoy,length=1640>
##contig=<ID=chrUn_JTFH01001053v1_decoy,length=1639>
##contig=<ID=chrUn_JTFH01000474v1_decoy,length=1638>
##contig=<ID=chrUn_JTFH01001418v1_decoy,length=1638>
##contig=<ID=chrUn_JTFH01001054v1_decoy,length=1636>
##contig=<ID=chrUn_JTFH01000475v1_decoy,length=1636>
##contig=<ID=chrUn_JTFH01001419v1_decoy,length=1633>
##contig=<ID=chrUn_JTFH01000476v1_decoy,length=1632>
##contig=<ID=chrUn_JTFH01001055v1_decoy,length=1632>
##contig=<ID=chrUn_JTFH01000477v1_decoy,length=1631>
##contig=<ID=chrUn_JTFH01000478v1_decoy,length=1630>
##contig=<ID=chrUn_JTFH01001056v1_decoy,length=1629>
##contig=<ID=chrUn_JTFH01000479v1_decoy,length=1627>
##contig=<ID=chrUn_JTFH01001420v1_decoy,length=1626>
##contig=<ID=chrUn_JTFH01000480v1_decoy,length=1624>
##contig=<ID=chrUn_KN707743v1_decoy,length=1624>
##contig=<ID=chrUn_JTFH01001057v1_decoy,length=1623>
##contig=<ID=chrUn_JTFH01001058v1_decoy,length=1622>
##contig=<ID=chrUn_JTFH01001059v1_decoy,length=1622>
##contig=<ID=chrUn_JTFH01001060v1_decoy,length=1619>
##contig=<ID=chrUn_KN707613v1_decoy,length=1619>
##contig=<ID=chrUn_JTFH01000481v1_decoy,length=1617>
##contig=<ID=chrUn_JTFH01000482v1_decoy,length=1616>
##contig=<ID=chrUn_JTFH01000483v1_decoy,length=1615>
##contig=<ID=chrUn_JTFH01001421v1_decoy,length=1614>
##contig=<ID=chrUn_JTFH01001422v1_decoy,length=1612>
##contig=<ID=chrUn_JTFH01000484v1_decoy,length=1611>
##contig=<ID=chrUn_JTFH01000485v1_decoy,length=1611>
##contig=<ID=chrUn_JTFH01000486v1_decoy,length=1606>
##contig=<ID=chrUn_JTFH01001061v1_decoy,length=1606>
##contig=<ID=chrUn_JTFH01000488v1_decoy,length=1605>
##contig=<ID=chrUn_JTFH01001423v1_decoy,length=1605>
##contig=<ID=chrUn_JTFH01000487v1_decoy,length=1605>
##contig=<ID=chrUn_JTFH01001424v1_decoy,length=1603>
##contig=<ID=chrUn_JTFH01000489v1_decoy,length=1600>
##contig=<ID=chrUn_KI270548v1,length=1599>
##contig=<ID=chrUn_JTFH01001425v1_decoy,length=1599>
##contig=<ID=chrUn_JTFH01000491v1_decoy,length=1598>
##contig=<ID=chrUn_JTFH01000490v1_decoy,length=1598>
##contig=<ID=chrUn_JTFH01000492v1_decoy,length=1597>
##contig=<ID=chrUn_JTFH01000493v1_decoy,length=1596>
##contig=<ID=chrUn_JTFH01000494v1_decoy,length=1595>
##contig=<ID=chrUn_JTFH01001062v1_decoy,length=1593>
##contig=<ID=chrUn_JTFH01001063v1_decoy,length=1592>
##contig=<ID=chrUn_JTFH01000495v1_decoy,length=1592>
##contig=<ID=chrUn_KN707769v1_decoy,length=1591>
##contig=<ID=chrUn_JTFH01001426v1_decoy,length=1589>
##contig=<ID=chrUn_JTFH01000496v1_decoy,length=1589>
##contig=<ID=chrUn_JTFH01001427v1_decoy,length=1588>
##contig=<ID=chrUn_JTFH01001428v1_decoy,length=1585>
##contig=<ID=chrUn_JTFH01000497v1_decoy,length=1585>
##contig=<ID=chrUn_JTFH01001430v1_decoy,length=1584>
##contig=<ID=chrUn_JTFH01001429v1_decoy,length=1584>
##contig=<ID=chrUn_JTFH01001431v1_decoy,length=1580>
##contig=<ID=chrUn_KN707852v1_decoy,length=1579>
##contig=<ID=chrUn_JTFH01000498v1_decoy,length=1579>
##contig=<ID=chrUn_JTFH01000499v1_decoy,length=1578>
##contig=<ID=chrUn_JTFH01000502v1_decoy,length=1577>
##contig=<ID=chrUn_JTFH01000500v1_decoy,length=1577>
##contig=<ID=chrUn_JTFH01000501v1_decoy,length=1577>
##contig=<ID=chrUn_JTFH01000503v1_decoy,length=1576>
##contig=<ID=chrUn_JTFH01000504v1_decoy,length=1575>
##contig=<ID=chrUn_JTFH01000505v1_decoy,length=1574>
##contig=<ID=chrUn_JTFH01001432v1_decoy,length=1572>
##contig=<ID=chrUn_JTFH01000506v1_decoy,length=1572>
##contig=<ID=chrUn_JTFH01000507v1_decoy,length=1571>
##contig=<ID=chrUn_JTFH01001433v1_decoy,length=1570>
##contig=<ID=chrUn_JTFH01001434v1_decoy,length=1569>
##contig=<ID=chrUn_JTFH01001435v1_decoy,length=1568>
##contig=<ID=chrUn_JTFH01001436v1_decoy,length=1567>
##contig=<ID=chrUn_JTFH01001437v1_decoy,length=1565>
##contig=<ID=chrUn_KN707648v1_decoy,length=1564>
##contig=<ID=chrUn_JTFH01000508v1_decoy,length=1563>
##contig=<ID=chrUn_JTFH01000509v1_decoy,length=1561>
##contig=<ID=chrUn_JTFH01000510v1_decoy,length=1561>
##contig=<ID=chrUn_JTFH01000511v1_decoy,length=1560>
##contig=<ID=chrUn_JTFH01000512v1_decoy,length=1560>
##contig=<ID=chrUn_JTFH01001438v1_decoy,length=1559>
##contig=<ID=chrUn_JTFH01001439v1_decoy,length=1559>
##contig=<ID=chrUn_JTFH01001064v1_decoy,length=1558>
##contig=<ID=chrUn_JTFH01001440v1_decoy,length=1556>
##contig=<ID=chrUn_JTFH01000513v1_decoy,length=1554>
##contig=<ID=chrUn_JTFH01001441v1_decoy,length=1554>
##contig=<ID=chrUn_KI270580v1,length=1553>
##contig=<ID=chrUn_JTFH01000514v1_decoy,length=1552>
##contig=<ID=chrUn_KN707619v1_decoy,length=1551>
##contig=<ID=chrUn_JTFH01001442v1_decoy,length=1549>
##contig=<ID=chrUn_JTFH01000515v1_decoy,length=1548>
##contig=<ID=chrUn_JTFH01000516v1_decoy,length=1546>
##contig=<ID=chrUn_JTFH01001065v1_decoy,length=1545>
##contig=<ID=chrUn_KN707867v1_decoy,length=1544>
##contig=<ID=chrUn_JTFH01001443v1_decoy,length=1542>
##contig=<ID=chrUn_JTFH01001066v1_decoy,length=1542>
##contig=<ID=chrUn_JTFH01001444v1_decoy,length=1541>
##contig=<ID=chrUn_JTFH01000517v1_decoy,length=1541>
##contig=<ID=chrUn_KN707650v1_decoy,length=1540>
##contig=<ID=chrUn_JTFH01001067v1_decoy,length=1540>
##contig=<ID=chrUn_JTFH01001445v1_decoy,length=1538>
##contig=<ID=chrUn_JTFH01001446v1_decoy,length=1537>
##contig=<ID=chrUn_KI270387v1,length=1537>
##contig=<ID=chrUn_JTFH01000518v1_decoy,length=1536>
##contig=<ID=chrUn_KN707622v1_decoy,length=1535>
##contig=<ID=chrUn_JTFH01001447v1_decoy,length=1535>
##contig=<ID=chrUn_JTFH01000519v1_decoy,length=1533>
##contig=<ID=chrUn_JTFH01000520v1_decoy,length=1532>
##contig=<ID=chrUn_JTFH01000521v1_decoy,length=1532>
##contig=<ID=chrUn_JTFH01001448v1_decoy,length=1530>
##contig=<ID=chrUn_JTFH01000522v1_decoy,length=1530>
##contig=<ID=chrUn_JTFH01001068v1_decoy,length=1529>
##contig=<ID=chrUn_JTFH01001449v1_decoy,length=1528>
##contig=<ID=chrUn_JTFH01000523v1_decoy,length=1527>
##contig=<ID=chrUn_KN707913v1_decoy,length=1527>
##contig=<ID=chrUn_JTFH01000524v1_decoy,length=1526>
##contig=<ID=chrUn_JTFH01000525v1_decoy,length=1524>
##contig=<ID=chrUn_JTFH01000526v1_decoy,length=1523>
##contig=<ID=chrUn_JTFH01000527v1_decoy,length=1523>
##contig=<ID=chrUn_KN707962v1_decoy,length=1523>
##contig=<ID=chrUn_JTFH01000529v1_decoy,length=1522>
##contig=<ID=chrUn_JTFH01001450v1_decoy,length=1522>
##contig=<ID=chrUn_JTFH01000528v1_decoy,length=1522>
##contig=<ID=chrUn_JTFH01000530v1_decoy,length=1519>
##contig=<ID=chrUn_JTFH01001069v1_decoy,length=1518>
##contig=<ID=chrUn_JTFH01001070v1_decoy,length=1515>
##contig=<ID=chrUn_JTFH01001451v1_decoy,length=1514>
##contig=<ID=chrUn_JTFH01001071v1_decoy,length=1513>
##contig=<ID=chrUn_JTFH01000531v1_decoy,length=1513>
##contig=<ID=chrUn_JTFH01001452v1_decoy,length=1509>
##contig=<ID=chrUn_JTFH01000533v1_decoy,length=1508>
##contig=<ID=chrUn_JTFH01000532v1_decoy,length=1508>
##contig=<ID=chrUn_JTFH01001453v1_decoy,length=1507>
##contig=<ID=chrUn_JTFH01001072v1_decoy,length=1507>
##contig=<ID=chrUn_JTFH01000534v1_decoy,length=1505>
##contig=<ID=chrUn_JTFH01001073v1_decoy,length=1504>
##contig=<ID=chrUn_JTFH01000535v1_decoy,length=1503>
##contig=<ID=chrUn_JTFH01001454v1_decoy,length=1500>
##contig=<ID=chrUn_JTFH01001456v1_decoy,length=1499>
##contig=<ID=chrUn_JTFH01001455v1_decoy,length=1499>
##contig=<ID=chrUn_JTFH01001074v1_decoy,length=1499>
##contig=<ID=chrUn_JTFH01001457v1_decoy,length=1497>
##contig=<ID=chrUn_JTFH01000536v1_decoy,length=1496>
##contig=<ID=chrUn_JTFH01001458v1_decoy,length=1496>
##contig=<ID=chrUn_JTFH01001075v1_decoy,length=1495>
##contig=<ID=chrUn_JTFH01001076v1_decoy,length=1495>
##contig=<ID=chrUn_JTFH01001077v1_decoy,length=1492>
##contig=<ID=chrUn_JTFH01001078v1_decoy,length=1492>
##contig=<ID=chrUn_JTFH01000537v1_decoy,length=1491>
##contig=<ID=chrUn_JTFH01000538v1_decoy,length=1490>
##contig=<ID=chrUn_JTFH01000539v1_decoy,length=1490>
##contig=<ID=chrUn_JTFH01001079v1_decoy,length=1489>
##contig=<ID=chrUn_JTFH01001459v1_decoy,length=1488>
##contig=<ID=chrUn_JTFH01000540v1_decoy,length=1487>
##contig=<ID=chrUn_JTFH01000541v1_decoy,length=1486>
##contig=<ID=chrUn_JTFH01001460v1_decoy,length=1486>
##contig=<ID=chrUn_JTFH01000542v1_decoy,length=1485>
##contig=<ID=chrUn_JTFH01001461v1_decoy,length=1485>
##contig=<ID=chrUn_JTFH01001080v1_decoy,length=1485>
##contig=<ID=chrUn_KI270391v1,length=1484>
##contig=<ID=chrUn_JTFH01000543v1_decoy,length=1484>
##contig=<ID=chrUn_JTFH01001081v1_decoy,length=1483>
##contig=<ID=chrUn_JTFH01000544v1_decoy,length=1483>
##contig=<ID=chrUn_JTFH01001462v1_decoy,length=1481>
##contig=<ID=chrUn_JTFH01001463v1_decoy,length=1479>
##contig=<ID=chrUn_JTFH01000546v1_decoy,length=1479>
##contig=<ID=chrUn_JTFH01000545v1_decoy,length=1479>
##contig=<ID=chrUn_JTFH01000547v1_decoy,length=1476>
##contig=<ID=chrUn_JTFH01000548v1_decoy,length=1475>
##contig=<ID=chrUn_JTFH01001082v1_decoy,length=1473>
##contig=<ID=chrUn_KI270305v1,length=1472>
##contig=<ID=chrUn_JTFH01000549v1_decoy,length=1472>
##contig=<ID=chrUn_JTFH01001465v1_decoy,length=1472>
##contig=<ID=chrUn_JTFH01001464v1_decoy,length=1472>
##contig=<ID=chrUn_JTFH01001466v1_decoy,length=1470>
##contig=<ID=chrUn_JTFH01001083v1_decoy,length=1470>
##contig=<ID=chrUn_JTFH01000550v1_decoy,length=1469>
##contig=<ID=chrUn_JTFH01000551v1_decoy,length=1468>
##contig=<ID=chrUn_JTFH01000552v1_decoy,length=1467>
##contig=<ID=chrUn_JTFH01001467v1_decoy,length=1466>
##contig=<ID=chrUn_JTFH01001468v1_decoy,length=1465>
##contig=<ID=chrUn_JTFH01000553v1_decoy,length=1465>
##contig=<ID=chrUn_JTFH01000554v1_decoy,length=1464>
##contig=<ID=chrUn_JTFH01000555v1_decoy,length=1463>
##contig=<ID=chrUn_JTFH01001084v1_decoy,length=1463>
##contig=<ID=chrUn_JTFH01000556v1_decoy,length=1463>
##contig=<ID=chrUn_JTFH01001469v1_decoy,length=1461>
##contig=<ID=chrUn_JTFH01001085v1_decoy,length=1460>
##contig=<ID=chrUn_JTFH01000557v1_decoy,length=1459>
##contig=<ID=chrUn_JTFH01000558v1_decoy,length=1459>
##contig=<ID=chrUn_JTFH01000559v1_decoy,length=1458>
##contig=<ID=chrUn_JTFH01001086v1_decoy,length=1458>
##contig=<ID=chrUn_JTFH01000560v1_decoy,length=1458>
##contig=<ID=chrUn_JTFH01001470v1_decoy,length=1458>
##contig=<ID=chrUn_JTFH01001471v1_decoy,length=1457>
##contig=<ID=chrUn_JTFH01001087v1_decoy,length=1456>
##contig=<ID=chrUn_JTFH01000561v1_decoy,length=1454>
##contig=<ID=chrUn_KN707724v1_decoy,length=1454>
##contig=<ID=chrUn_JTFH01001088v1_decoy,length=1453>
##contig=<ID=chrUn_KI270373v1,length=1451>
##contig=<ID=chrUn_JTFH01000562v1_decoy,length=1449>
##contig=<ID=chrUn_JTFH01000563v1_decoy,length=1449>
##contig=<ID=chrUn_JTFH01000564v1_decoy,length=1448>
##contig=<ID=chrUn_JTFH01001472v1_decoy,length=1448>
##contig=<ID=chrUn_JTFH01001473v1_decoy,length=1447>
##contig=<ID=chrUn_JTFH01000565v1_decoy,length=1446>
##contig=<ID=chrUn_KI270422v1,length=1445>
##contig=<ID=chrUn_JTFH01001474v1_decoy,length=1444>
##contig=<ID=chrUn_KI270316v1,length=1444>
##contig=<ID=chrUn_KN707789v1_decoy,length=1443>
##contig=<ID=chrUn_JTFH01001475v1_decoy,length=1443>
##contig=<ID=chrUn_JTFH01001089v1_decoy,length=1443>
##contig=<ID=chrUn_JTFH01001476v1_decoy,length=1443>
##contig=<ID=chrUn_JTFH01000566v1_decoy,length=1442>
##contig=<ID=chrUn_JTFH01000567v1_decoy,length=1441>
##contig=<ID=chrUn_JTFH01001090v1_decoy,length=1441>
##contig=<ID=chrUn_JTFH01000568v1_decoy,length=1440>
##contig=<ID=chrUn_JTFH01000569v1_decoy,length=1439>
##contig=<ID=chrUn_JTFH01001477v1_decoy,length=1438>
##contig=<ID=chrUn_JTFH01000570v1_decoy,length=1437>
##contig=<ID=chrUn_JTFH01000571v1_decoy,length=1436>
##contig=<ID=chrUn_KN707948v1_decoy,length=1432>
##contig=<ID=chrUn_JTFH01001478v1_decoy,length=1432>
##contig=<ID=chrUn_JTFH01001479v1_decoy,length=1430>
##contig=<ID=chrUn_JTFH01001480v1_decoy,length=1430>
##contig=<ID=chrUn_JTFH01001483v1_decoy,length=1429>
##contig=<ID=chrUn_JTFH01001482v1_decoy,length=1429>
##contig=<ID=chrUn_JTFH01000572v1_decoy,length=1429>
##contig=<ID=chrUn_JTFH01001481v1_decoy,length=1429>
##contig=<ID=chrUn_JTFH01000573v1_decoy,length=1429>
##contig=<ID=chrUn_KI270338v1,length=1428>
##contig=<ID=chrUn_KI270340v1,length=1428>
##contig=<ID=chrUn_JTFH01000574v1_decoy,length=1427>
##contig=<ID=chrUn_JTFH01000575v1_decoy,length=1426>
##contig=<ID=chrUn_JTFH01001484v1_decoy,length=1426>
##contig=<ID=chrUn_JTFH01001485v1_decoy,length=1426>
##contig=<ID=chrUn_JTFH01001091v1_decoy,length=1426>
##contig=<ID=chrUn_JTFH01000576v1_decoy,length=1425>
##contig=<ID=chrUn_JTFH01001092v1_decoy,length=1425>
##contig=<ID=chrUn_JTFH01000578v1_decoy,length=1424>
##contig=<ID=chrUn_KN707632v1_decoy,length=1424>
##contig=<ID=chrUn_JTFH01000577v1_decoy,length=1424>
##contig=<ID=chrUn_JTFH01000581v1_decoy,length=1423>
##contig=<ID=chrUn_JTFH01000579v1_decoy,length=1423>
##contig=<ID=chrUn_JTFH01000580v1_decoy,length=1423>
##contig=<ID=chrUn_JTFH01001486v1_decoy,length=1420>
##contig=<ID=chrUn_JTFH01001093v1_decoy,length=1418>
##contig=<ID=chrUn_JTFH01001487v1_decoy,length=1416>
##contig=<ID=chrUn_JTFH01001488v1_decoy,length=1416>
##contig=<ID=chrUn_JTFH01001489v1_decoy,length=1415>
##contig=<ID=chrUn_JTFH01001490v1_decoy,length=1415>
##contig=<ID=chrUn_KN707635v1_decoy,length=1414>
##contig=<ID=chrUn_JTFH01000582v1_decoy,length=1414>
##contig=<ID=chrUn_JTFH01000583v1_decoy,length=1414>
##contig=<ID=chrUn_JTFH01001491v1_decoy,length=1414>
##contig=<ID=chrUn_JTFH01001492v1_decoy,length=1413>
##contig=<ID=chrUn_JTFH01001095v1_decoy,length=1413>
##contig=<ID=chrUn_JTFH01001094v1_decoy,length=1413>
##contig=<ID=chrUn_JTFH01000584v1_decoy,length=1413>
##contig=<ID=chrUn_JTFH01000585v1_decoy,length=1413>
##contig=<ID=chrUn_JTFH01001096v1_decoy,length=1412>
##contig=<ID=chrUn_KN707788v1_decoy,length=1412>
##contig=<ID=chrUn_JTFH01000586v1_decoy,length=1410>
##contig=<ID=chrUn_JTFH01001493v1_decoy,length=1410>
##contig=<ID=chrUn_JTFH01000588v1_decoy,length=1409>
##contig=<ID=chrUn_KN707923v1_decoy,length=1409>
##contig=<ID=chrUn_JTFH01000587v1_decoy,length=1409>
##contig=<ID=chrUn_JTFH01001097v1_decoy,length=1407>
##contig=<ID=chrUn_JTFH01001098v1_decoy,length=1406>
##contig=<ID=chrUn_JTFH01000589v1_decoy,length=1406>
##contig=<ID=chrUn_JTFH01000590v1_decoy,length=1405>
##contig=<ID=chrUn_JTFH01001494v1_decoy,length=1405>
##contig=<ID=chrUn_JTFH01000591v1_decoy,length=1405>
##contig=<ID=chrUn_JTFH01000593v1_decoy,length=1404>
##contig=<ID=chrUn_JTFH01000592v1_decoy,length=1404>
##contig=<ID=chrUn_JTFH01000596v1_decoy,length=1402>
##contig=<ID=chrUn_JTFH01000594v1_decoy,length=1402>
##contig=<ID=chrUn_JTFH01000595v1_decoy,length=1402>
##contig=<ID=chrUn_JTFH01000597v1_decoy,length=1402>
##contig=<ID=chrUn_JTFH01001495v1_decoy,length=1402>
##contig=<ID=chrUn_KI270583v1,length=1400>
##contig=<ID=chrUn_JTFH01000598v1_decoy,length=1400>
##contig=<ID=chrUn_JTFH01001496v1_decoy,length=1398>
##contig=<ID=chrUn_JTFH01000599v1_decoy,length=1398>
##contig=<ID=chrUn_JTFH01001497v1_decoy,length=1397>
##contig=<ID=chrUn_JTFH01000600v1_decoy,length=1396>
##contig=<ID=chrUn_JTFH01001099v1_decoy,length=1396>
##contig=<ID=chrUn_JTFH01000601v1_decoy,length=1395>
##contig=<ID=chrUn_JTFH01001498v1_decoy,length=1395>
##contig=<ID=chrUn_JTFH01000602v1_decoy,length=1394>
##contig=<ID=chrUn_KN707610v1_decoy,length=1393>
##contig=<ID=chrUn_JTFH01000603v1_decoy,length=1393>
##contig=<ID=chrUn_JTFH01001499v1_decoy,length=1392>
##contig=<ID=chrUn_JTFH01000604v1_decoy,length=1391>
##contig=<ID=chrUn_JTFH01001100v1_decoy,length=1390>
##contig=<ID=chrUn_JTFH01000606v1_decoy,length=1389>
##contig=<ID=chrUn_JTFH01000605v1_decoy,length=1389>
##contig=<ID=chrUn_JTFH01001500v1_decoy,length=1388>
##contig=<ID=chrUn_JTFH01000607v1_decoy,length=1388>
##contig=<ID=chrUn_JTFH01000608v1_decoy,length=1387>
##contig=<ID=chrUn_JTFH01001501v1_decoy,length=1386>
##contig=<ID=chrUn_JTFH01000609v1_decoy,length=1384>
##contig=<ID=chrUn_KN707952v1_decoy,length=1383>
##contig=<ID=chrUn_JTFH01001502v1_decoy,length=1382>
##contig=<ID=chrUn_JTFH01001101v1_decoy,length=1382>
##contig=<ID=chrUn_JTFH01000611v1_decoy,length=1381>
##contig=<ID=chrUn_JTFH01001503v1_decoy,length=1381>
##contig=<ID=chrUn_JTFH01000610v1_decoy,length=1381>
##contig=<ID=chrUn_JTFH01001504v1_decoy,length=1379>
##contig=<ID=chrUn_JTFH01000612v1_decoy,length=1379>
##contig=<ID=chrUn_JTFH01000613v1_decoy,length=1377>
##contig=<ID=chrUn_JTFH01000615v1_decoy,length=1376>
##contig=<ID=chrUn_JTFH01000614v1_decoy,length=1376>
##contig=<ID=chrUn_JTFH01001102v1_decoy,length=1376>
##contig=<ID=chrUn_JTFH01001505v1_decoy,length=1376>
##contig=<ID=chrUn_JTFH01000616v1_decoy,length=1375>
##contig=<ID=chrUn_JTFH01001103v1_decoy,length=1375>
##contig=<ID=chrUn_JTFH01001506v1_decoy,length=1374>
##contig=<ID=chrUn_JTFH01000617v1_decoy,length=1374>
##contig=<ID=chrUn_JTFH01001507v1_decoy,length=1374>
##contig=<ID=chrUn_JTFH01001509v1_decoy,length=1373>
##contig=<ID=chrUn_JTFH01001508v1_decoy,length=1373>
##contig=<ID=chrUn_JTFH01001510v1_decoy,length=1372>
##contig=<ID=chrUn_JTFH01000618v1_decoy,length=1372>
##contig=<ID=chrUn_JTFH01001104v1_decoy,length=1371>
##contig=<ID=chrUn_JTFH01000619v1_decoy,length=1371>
##contig=<ID=chrUn_JTFH01000621v1_decoy,length=1370>
##contig=<ID=chrUn_JTFH01001511v1_decoy,length=1370>
##contig=<ID=chrUn_JTFH01000620v1_decoy,length=1370>
##contig=<ID=chrUn_KI270334v1,length=1368>
##contig=<ID=chrUn_JTFH01001105v1_decoy,length=1367>
##contig=<ID=chrUn_JTFH01001512v1_decoy,length=1367>
##contig=<ID=chrUn_KN707894v1_decoy,length=1366>
##contig=<ID=chrUn_JTFH01000622v1_decoy,length=1366>
##contig=<ID=chrUn_JTFH01001513v1_decoy,length=1365>
##contig=<ID=chrUn_JTFH01001514v1_decoy,length=1364>
##contig=<ID=chrUn_JTFH01001106v1_decoy,length=1364>
##contig=<ID=chrUn_JTFH01000623v1_decoy,length=1363>
##contig=<ID=chrUn_JTFH01001516v1_decoy,length=1361>
##contig=<ID=chrUn_KI270429v1,length=1361>
##contig=<ID=chrUn_JTFH01001515v1_decoy,length=1361>
##contig=<ID=chrUn_JTFH01000624v1_decoy,length=1360>
##contig=<ID=chrUn_JTFH01000625v1_decoy,length=1356>
##contig=<ID=chrUn_JTFH01001107v1_decoy,length=1356>
##contig=<ID=chrUn_JTFH01000626v1_decoy,length=1355>
##contig=<ID=chrUn_JTFH01001518v1_decoy,length=1355>
##contig=<ID=chrUn_JTFH01001517v1_decoy,length=1355>
##contig=<ID=chrUn_JTFH01000627v1_decoy,length=1355>
##contig=<ID=chrUn_JTFH01001108v1_decoy,length=1355>
##contig=<ID=chrUn_JTFH01001519v1_decoy,length=1354>
##contig=<ID=chrUn_JTFH01001520v1_decoy,length=1353>
##contig=<ID=chrUn_KN707824v1_decoy,length=1352>
##contig=<ID=chrUn_JTFH01001109v1_decoy,length=1352>
##contig=<ID=chrUn_JTFH01000628v1_decoy,length=1352>
##contig=<ID=chrUn_JTFH01001110v1_decoy,length=1350>
##contig=<ID=chrUn_JTFH01001521v1_decoy,length=1349>
##contig=<ID=chrUn_JTFH01001111v1_decoy,length=1346>
##contig=<ID=chrUn_JTFH01000629v1_decoy,length=1345>
##contig=<ID=chrUn_JTFH01001522v1_decoy,length=1345>
##contig=<ID=chrUn_JTFH01001112v1_decoy,length=1345>
##contig=<ID=chrUn_JTFH01000630v1_decoy,length=1344>
##contig=<ID=chrUn_JTFH01000631v1_decoy,length=1344>
##contig=<ID=chrUn_JTFH01001523v1_decoy,length=1344>
##contig=<ID=chrUn_JTFH01001524v1_decoy,length=1343>
##contig=<ID=chrUn_JTFH01000632v1_decoy,length=1342>
##contig=<ID=chrUn_JTFH01000633v1_decoy,length=1342>
##contig=<ID=chrUn_JTFH01001113v1_decoy,length=1340>
##contig=<ID=chrUn_KN707859v1_decoy,length=1340>
##contig=<ID=chrUn_JTFH01001526v1_decoy,length=1338>
##contig=<ID=chrUn_JTFH01001525v1_decoy,length=1338>
##contig=<ID=chrUn_JTFH01001527v1_decoy,length=1338>
##contig=<ID=chrUn_JTFH01000634v1_decoy,length=1336>
##contig=<ID=chrUn_JTFH01001528v1_decoy,length=1336>
##contig=<ID=chrUn_JTFH01000635v1_decoy,length=1334>
##contig=<ID=chrUn_JTFH01000636v1_decoy,length=1334>
##contig=<ID=chrUn_JTFH01001530v1_decoy,length=1333>
##contig=<ID=chrUn_JTFH01001529v1_decoy,length=1333>
##contig=<ID=chrUn_JTFH01000637v1_decoy,length=1333>
##contig=<ID=chrUn_JTFH01000638v1_decoy,length=1332>
##contig=<ID=chrUn_JTFH01001531v1_decoy,length=1332>
##contig=<ID=chrUn_KN707819v1_decoy,length=1331>
##contig=<ID=chrUn_JTFH01001114v1_decoy,length=1330>
##contig=<ID=chrUn_JTFH01001115v1_decoy,length=1329>
##contig=<ID=chrUn_KN707624v1_decoy,length=1329>
##contig=<ID=chrUn_JTFH01000639v1_decoy,length=1328>
##contig=<ID=chrUn_JTFH01000641v1_decoy,length=1328>
##contig=<ID=chrUn_JTFH01000640v1_decoy,length=1328>
##contig=<ID=chrUn_JTFH01000642v1_decoy,length=1327>
##contig=<ID=chrUn_JTFH01000643v1_decoy,length=1325>
##contig=<ID=chrUn_JTFH01001116v1_decoy,length=1324>
##contig=<ID=chrUn_JTFH01001532v1_decoy,length=1324>
##contig=<ID=chrUn_JTFH01001533v1_decoy,length=1323>
##contig=<ID=chrUn_JTFH01001534v1_decoy,length=1323>
##contig=<ID=chrUn_JTFH01000644v1_decoy,length=1322>
##contig=<ID=chrUn_JTFH01001536v1_decoy,length=1320>
##contig=<ID=chrUn_JTFH01001535v1_decoy,length=1320>
##contig=<ID=chrUn_JTFH01000645v1_decoy,length=1320>
##contig=<ID=chrUn_JTFH01000646v1_decoy,length=1319>
##contig=<ID=chrUn_JTFH01000647v1_decoy,length=1318>
##contig=<ID=chrUn_JTFH01001537v1_decoy,length=1317>
##contig=<ID=chrUn_JTFH01001117v1_decoy,length=1316>
##contig=<ID=chrUn_JTFH01001538v1_decoy,length=1316>
##contig=<ID=chrUn_JTFH01000648v1_decoy,length=1315>
##contig=<ID=chrUn_JTFH01000649v1_decoy,length=1314>
##contig=<ID=chrUn_JTFH01000651v1_decoy,length=1313>
##contig=<ID=chrUn_JTFH01000650v1_decoy,length=1313>
##contig=<ID=chrUn_JTFH01000652v1_decoy,length=1312>
##contig=<ID=chrUn_KN707935v1_decoy,length=1312>
##contig=<ID=chrUn_JTFH01000653v1_decoy,length=1310>
##contig=<ID=chrUn_JTFH01000654v1_decoy,length=1309>
##contig=<ID=chrUn_JTFH01000655v1_decoy,length=1309>
##contig=<ID=chrUn_KI270393v1,length=1308>
##contig=<ID=chrUn_KN707713v1_decoy,length=1308>
##contig=<ID=chrUn_JTFH01001118v1_decoy,length=1307>
##contig=<ID=chrUn_JTFH01000657v1_decoy,length=1307>
##contig=<ID=chrUn_JTFH01000656v1_decoy,length=1307>
##contig=<ID=chrUn_JTFH01000658v1_decoy,length=1305>
##contig=<ID=chrUn_JTFH01001539v1_decoy,length=1304>
##contig=<ID=chrUn_JTFH01001540v1_decoy,length=1304>
##contig=<ID=chrUn_JTFH01001120v1_decoy,length=1304>
##contig=<ID=chrUn_JTFH01001119v1_decoy,length=1304>
##contig=<ID=chrUn_JTFH01000659v1_decoy,length=1304>
##contig=<ID=chrUn_JTFH01000660v1_decoy,length=1303>
##contig=<ID=chrUn_JTFH01001121v1_decoy,length=1303>
##contig=<ID=chrUn_JTFH01001541v1_decoy,length=1303>
##contig=<ID=chrUn_JTFH01001542v1_decoy,length=1302>
##contig=<ID=chrUn_JTFH01000662v1_decoy,length=1302>
##contig=<ID=chrUn_JTFH01000661v1_decoy,length=1302>
##contig=<ID=chrUn_JTFH01001543v1_decoy,length=1301>
##contig=<ID=chrUn_JTFH01001122v1_decoy,length=1301>
##contig=<ID=chrUn_JTFH01000663v1_decoy,length=1301>
##contig=<ID=chrUn_JTFH01000664v1_decoy,length=1301>
##contig=<ID=chrUn_KN707942v1_decoy,length=1300>
##contig=<ID=chrUn_JTFH01000665v1_decoy,length=1300>
##contig=<ID=chrUn_JTFH01001544v1_decoy,length=1300>
##contig=<ID=chrUn_KI270516v1,length=1300>
##contig=<ID=chrUn_JTFH01001123v1_decoy,length=1300>
##contig=<ID=chrUn_JTFH01000666v1_decoy,length=1299>
##contig=<ID=chrUn_JTFH01001545v1_decoy,length=1298>
##contig=<ID=chrUn_KI270389v1,length=1298>
##contig=<ID=chrUn_JTFH01000667v1_decoy,length=1297>
##contig=<ID=chrUn_JTFH01001124v1_decoy,length=1297>
##contig=<ID=chrUn_KN707680v1_decoy,length=1297>
##contig=<ID=chrUn_JTFH01001546v1_decoy,length=1297>
##contig=<ID=chrUn_JTFH01001125v1_decoy,length=1296>
##contig=<ID=chrUn_JTFH01001547v1_decoy,length=1295>
##contig=<ID=chrUn_JTFH01000668v1_decoy,length=1295>
##contig=<ID=chrUn_JTFH01000669v1_decoy,length=1294>
##contig=<ID=chrUn_KN707639v1_decoy,length=1294>
##contig=<ID=chrUn_JTFH01000670v1_decoy,length=1293>
##contig=<ID=chrUn_JTFH01000672v1_decoy,length=1291>
##contig=<ID=chrUn_JTFH01000671v1_decoy,length=1291>
##contig=<ID=chrUn_JTFH01001126v1_decoy,length=1290>
##contig=<ID=chrUn_JTFH01000673v1_decoy,length=1289>
##contig=<ID=chrUn_JTFH01000674v1_decoy,length=1288>
##contig=<ID=chrUn_JTFH01000675v1_decoy,length=1288>
##contig=<ID=chrUn_JTFH01000676v1_decoy,length=1287>
##contig=<ID=chrUn_JTFH01000678v1_decoy,length=1287>
##contig=<ID=chrUn_JTFH01000677v1_decoy,length=1287>
##contig=<ID=chrUn_JTFH01000679v1_decoy,length=1286>
##contig=<ID=chrUn_JTFH01001548v1_decoy,length=1284>
##contig=<ID=chrUn_JTFH01001127v1_decoy,length=1284>
##contig=<ID=chrUn_JTFH01000680v1_decoy,length=1283>
##contig=<ID=chrUn_JTFH01001550v1_decoy,length=1283>
##contig=<ID=chrUn_JTFH01001549v1_decoy,length=1283>
##contig=<ID=chrUn_JTFH01001128v1_decoy,length=1282>
##contig=<ID=chrUn_JTFH01001129v1_decoy,length=1281>
##contig=<ID=chrUn_KN707912v1_decoy,length=1281>
##contig=<ID=chrUn_JTFH01000681v1_decoy,length=1281>
##contig=<ID=chrUn_JTFH01001130v1_decoy,length=1280>
##contig=<ID=chrUn_JTFH01001131v1_decoy,length=1279>
##contig=<ID=chrUn_JTFH01001551v1_decoy,length=1279>
##contig=<ID=chrUn_JTFH01001552v1_decoy,length=1278>
##contig=<ID=chrUn_JTFH01000682v1_decoy,length=1277>
##contig=<ID=chrUn_KN707745v1_decoy,length=1276>
##contig=<ID=chrUn_JTFH01000683v1_decoy,length=1274>
##contig=<ID=chrUn_KN707633v1_decoy,length=1274>
##contig=<ID=chrUn_JTFH01001132v1_decoy,length=1272>
##contig=<ID=chrUn_JTFH01001554v1_decoy,length=1271>
##contig=<ID=chrUn_JTFH01001553v1_decoy,length=1271>
##contig=<ID=chrUn_KN707737v1_decoy,length=1270>
##contig=<ID=chrUn_JTFH01000684v1_decoy,length=1270>
##contig=<ID=chrUn_JTFH01001555v1_decoy,length=1268>
##contig=<ID=chrUn_JTFH01001134v1_decoy,length=1267>
##contig=<ID=chrUn_JTFH01001133v1_decoy,length=1267>
##contig=<ID=chrUn_JTFH01000685v1_decoy,length=1267>
##contig=<ID=chrUn_JTFH01001135v1_decoy,length=1266>
##contig=<ID=chrUn_KN707926v1_decoy,length=1266>
##contig=<ID=chrUn_JTFH01000686v1_decoy,length=1266>
##contig=<ID=chrUn_KN707832v1_decoy,length=1265>
##contig=<ID=chrUn_KN707917v1_decoy,length=1265>
##contig=<ID=chrUn_JTFH01001556v1_decoy,length=1264>
##contig=<ID=chrUn_JTFH01001138v1_decoy,length=1264>
##contig=<ID=chrUn_JTFH01001137v1_decoy,length=1264>
##contig=<ID=chrUn_JTFH01001136v1_decoy,length=1264>
##contig=<ID=chrUn_JTFH01001139v1_decoy,length=1263>
##contig=<ID=chrUn_JTFH01001557v1_decoy,length=1263>
##contig=<ID=chrUn_JTFH01001558v1_decoy,length=1262>
##contig=<ID=chrUn_JTFH01001559v1_decoy,length=1261>
##contig=<ID=chrUn_JTFH01001560v1_decoy,length=1260>
##contig=<ID=chrUn_JTFH01000687v1_decoy,length=1260>
##contig=<ID=chrUn_JTFH01001562v1_decoy,length=1259>
##contig=<ID=chrUn_JTFH01001561v1_decoy,length=1259>
##contig=<ID=chrUn_JTFH01000688v1_decoy,length=1259>
##contig=<ID=chrUn_JTFH01000689v1_decoy,length=1258>
##contig=<ID=chrUn_JTFH01000691v1_decoy,length=1258>
##contig=<ID=chrUn_JTFH01000690v1_decoy,length=1258>
##contig=<ID=chrUn_JTFH01001563v1_decoy,length=1258>
##contig=<ID=chrUn_JTFH01001564v1_decoy,length=1256>
##contig=<ID=chrUn_JTFH01000692v1_decoy,length=1256>
##contig=<ID=chrUn_JTFH01000693v1_decoy,length=1255>
##contig=<ID=chrUn_JTFH01000694v1_decoy,length=1254>
##contig=<ID=chrUn_JTFH01000695v1_decoy,length=1254>
##contig=<ID=chrUn_JTFH01001565v1_decoy,length=1253>
##contig=<ID=chrUn_JTFH01000696v1_decoy,length=1253>
##contig=<ID=chrUn_JTFH01000697v1_decoy,length=1250>
##contig=<ID=chrUn_JTFH01000698v1_decoy,length=1249>
##contig=<ID=chrUn_JTFH01001140v1_decoy,length=1249>
##contig=<ID=chrUn_JTFH01000700v1_decoy,length=1248>
##contig=<ID=chrUn_JTFH01001567v1_decoy,length=1248>
##contig=<ID=chrUn_JTFH01001566v1_decoy,length=1248>
##contig=<ID=chrUn_JTFH01000699v1_decoy,length=1248>
##contig=<ID=chrUn_JTFH01000701v1_decoy,length=1247>
##contig=<ID=chrUn_JTFH01001569v1_decoy,length=1246>
##contig=<ID=chrUn_JTFH01001568v1_decoy,length=1246>
##contig=<ID=chrUn_JTFH01001570v1_decoy,length=1244>
##contig=<ID=chrUn_KN707906v1_decoy,length=1243>
##contig=<ID=chrUn_KN707712v1_decoy,length=1243>
##contig=<ID=chrUn_JTFH01000703v1_decoy,length=1242>
##contig=<ID=chrUn_JTFH01000702v1_decoy,length=1242>
##contig=<ID=chrUn_JTFH01000704v1_decoy,length=1241>
##contig=<ID=chrUn_JTFH01000705v1_decoy,length=1241>
##contig=<ID=chrUn_JTFH01000706v1_decoy,length=1241>
##contig=<ID=chrUn_JTFH01001141v1_decoy,length=1240>
##contig=<ID=chrUn_JTFH01000707v1_decoy,length=1239>
##contig=<ID=chrUn_JTFH01001142v1_decoy,length=1239>
##contig=<ID=chrUn_JTFH01000708v1_decoy,length=1238>
##contig=<ID=chrUn_JTFH01001572v1_decoy,length=1238>
##contig=<ID=chrUn_JTFH01001571v1_decoy,length=1238>
##contig=<ID=chrUn_KN707960v1_decoy,length=1238>
##contig=<ID=chrUn_KN707625v1_decoy,length=1238>
##contig=<ID=chrUn_JTFH01000709v1_decoy,length=1237>
##contig=<ID=chrUn_JTFH01000710v1_decoy,length=1236>
##contig=<ID=chrUn_JTFH01001573v1_decoy,length=1236>
##contig=<ID=chrUn_JTFH01001144v1_decoy,length=1235>
##contig=<ID=chrUn_JTFH01000711v1_decoy,length=1235>
##contig=<ID=chrUn_JTFH01001143v1_decoy,length=1235>
##contig=<ID=chrUn_JTFH01000713v1_decoy,length=1234>
##contig=<ID=chrUn_JTFH01000714v1_decoy,length=1234>
##contig=<ID=chrUn_JTFH01001574v1_decoy,length=1234>
##contig=<ID=chrUn_KN707847v1_decoy,length=1234>
##contig=<ID=chrUn_JTFH01000712v1_decoy,length=1234>
##contig=<ID=chrUn_JTFH01001575v1_decoy,length=1234>
##contig=<ID=chrUn_KI270466v1,length=1233>
##contig=<ID=chrUn_JTFH01000715v1_decoy,length=1233>
##contig=<ID=chrUn_JTFH01001145v1_decoy,length=1233>
##contig=<ID=chrUn_JTFH01001146v1_decoy,length=1232>
##contig=<ID=chrUn_JTFH01000716v1_decoy,length=1232>
##contig=<ID=chrUn_JTFH01000717v1_decoy,length=1232>
##contig=<ID=chrUn_JTFH01001576v1_decoy,length=1231>
##contig=<ID=chrUn_JTFH01001577v1_decoy,length=1231>
##contig=<ID=chrUn_JTFH01000718v1_decoy,length=1231>
##contig=<ID=chrUn_JTFH01001579v1_decoy,length=1230>
##contig=<ID=chrUn_JTFH01001147v1_decoy,length=1230>
##contig=<ID=chrUn_JTFH01000719v1_decoy,length=1230>
##contig=<ID=chrUn_JTFH01001578v1_decoy,length=1230>
##contig=<ID=chrUn_JTFH01000720v1_decoy,length=1228>
##contig=<ID=chrUn_KN707694v1_decoy,length=1228>
##contig=<ID=chrUn_JTFH01001580v1_decoy,length=1228>
##contig=<ID=chrUn_JTFH01001581v1_decoy,length=1227>
##contig=<ID=chrUn_JTFH01000722v1_decoy,length=1227>
##contig=<ID=chrUn_KN707732v1_decoy,length=1227>
##contig=<ID=chrUn_JTFH01000721v1_decoy,length=1227>
##contig=<ID=chrUn_JTFH01000723v1_decoy,length=1226>
##contig=<ID=chrUn_KN707929v1_decoy,length=1226>
##contig=<ID=chrUn_JTFH01001148v1_decoy,length=1226>
##contig=<ID=chrUn_JTFH01000724v1_decoy,length=1224>
##contig=<ID=chrUn_JTFH01000725v1_decoy,length=1224>
##contig=<ID=chrUn_JTFH01001149v1_decoy,length=1223>
##contig=<ID=chrUn_JTFH01001583v1_decoy,length=1222>
##contig=<ID=chrUn_JTFH01001582v1_decoy,length=1222>
##contig=<ID=chrUn_KN707621v1_decoy,length=1222>
##contig=<ID=chrUn_JTFH01001584v1_decoy,length=1221>
##contig=<ID=chrUn_JTFH01001585v1_decoy,length=1221>
##contig=<ID=chrUn_KN707759v1_decoy,length=1221>
##contig=<ID=chrUn_JTFH01001586v1_decoy,length=1220>
##contig=<ID=chrUn_JTFH01000727v1_decoy,length=1220>
##contig=<ID=chrUn_JTFH01000726v1_decoy,length=1220>
##contig=<ID=chrUn_JTFH01000728v1_decoy,length=1219>
##contig=<ID=chrUn_JTFH01001587v1_decoy,length=1218>
##contig=<ID=chrUn_JTFH01001588v1_decoy,length=1218>
##contig=<ID=chrUn_JTFH01000729v1_decoy,length=1217>
##contig=<ID=chrUn_JTFH01000730v1_decoy,length=1216>
##contig=<ID=chrUn_JTFH01001589v1_decoy,length=1216>
##contig=<ID=chrUn_JTFH01001590v1_decoy,length=1216>
##contig=<ID=chrUn_KI270388v1,length=1216>
##contig=<ID=chrUn_JTFH01000731v1_decoy,length=1215>
##contig=<ID=chrUn_JTFH01000732v1_decoy,length=1214>
##contig=<ID=chrUn_JTFH01000733v1_decoy,length=1214>
##contig=<ID=chrUn_JTFH01000734v1_decoy,length=1214>
##contig=<ID=chrUn_JTFH01001150v1_decoy,length=1214>
##contig=<ID=chrUn_JTFH01001151v1_decoy,length=1213>
##contig=<ID=chrUn_JTFH01000735v1_decoy,length=1213>
##contig=<ID=chrUn_JTFH01001591v1_decoy,length=1212>
##contig=<ID=chrUn_JTFH01000736v1_decoy,length=1212>
##contig=<ID=chrUn_JTFH01001152v1_decoy,length=1211>
##contig=<ID=chrUn_JTFH01001592v1_decoy,length=1210>
##contig=<ID=chrUn_JTFH01000737v1_decoy,length=1209>
##contig=<ID=chrUn_KN707770v1_decoy,length=1209>
##contig=<ID=chrUn_JTFH01001153v1_decoy,length=1209>
##contig=<ID=chrUn_JTFH01001593v1_decoy,length=1209>
##contig=<ID=chrUn_JTFH01000738v1_decoy,length=1208>
##contig=<ID=chrUn_KN707880v1_decoy,length=1208>
##contig=<ID=chrUn_JTFH01001594v1_decoy,length=1208>
##contig=<ID=HLA-B*15:01:01:02N,length=1208>
##contig=<ID=chrUn_JTFH01001595v1_decoy,length=1208>
##contig=<ID=chrUn_JTFH01000739v1_decoy,length=1207>
##contig=<ID=chrUn_JTFH01000740v1_decoy,length=1207>
##contig=<ID=chrUn_KN707818v1_decoy,length=1207>
##contig=<ID=chrUn_JTFH01000741v1_decoy,length=1207>
##contig=<ID=chrUn_JTFH01000743v1_decoy,length=1206>
##contig=<ID=chrUn_JTFH01000742v1_decoy,length=1206>
##contig=<ID=chrUn_JTFH01001596v1_decoy,length=1206>
##contig=<ID=chrUn_KN707848v1_decoy,length=1205>
##contig=<ID=chrUn_JTFH01001598v1_decoy,length=1205>
##contig=<ID=chrUn_JTFH01001597v1_decoy,length=1205>
##contig=<ID=chrUn_KN707670v1_decoy,length=1205>
##contig=<ID=chrUn_JTFH01000744v1_decoy,length=1205>
##contig=<ID=chrUn_JTFH01000745v1_decoy,length=1205>
##contig=<ID=chrUn_JTFH01000746v1_decoy,length=1204>
##contig=<ID=chrUn_JTFH01000748v1_decoy,length=1204>
##contig=<ID=chrUn_JTFH01000747v1_decoy,length=1204>
##contig=<ID=chrUn_JTFH01000749v1_decoy,length=1203>
##contig=<ID=chrUn_JTFH01001154v1_decoy,length=1202>
##contig=<ID=chrUn_KI270544v1,length=1202>
##contig=<ID=chrUn_JTFH01001599v1_decoy,length=1202>
##contig=<ID=chrUn_JTFH01000751v1_decoy,length=1201>
##contig=<ID=chrUn_JTFH01000750v1_decoy,length=1201>
##contig=<ID=chrUn_KI270310v1,length=1201>
##contig=<ID=chrUn_JTFH01000752v1_decoy,length=1200>
##contig=<ID=chrUn_JTFH01000753v1_decoy,length=1200>
##contig=<ID=chrUn_JTFH01001600v1_decoy,length=1200>
##contig=<ID=chrUn_KN707931v1_decoy,length=1199>
##contig=<ID=chrUn_JTFH01001155v1_decoy,length=1199>
##contig=<ID=chrUn_JTFH01001601v1_decoy,length=1199>
##contig=<ID=chrUn_JTFH01000754v1_decoy,length=1199>
##contig=<ID=chrUn_JTFH01000755v1_decoy,length=1198>
##contig=<ID=chrUn_JTFH01001604v1_decoy,length=1198>
##contig=<ID=chrUn_JTFH01001603v1_decoy,length=1198>
##contig=<ID=chrUn_JTFH01001602v1_decoy,length=1198>
##contig=<ID=chrUn_JTFH01000756v1_decoy,length=1197>
##contig=<ID=chrUn_JTFH01001156v1_decoy,length=1197>
##contig=<ID=chrUn_JTFH01000757v1_decoy,length=1196>
##contig=<ID=chrUn_JTFH01000758v1_decoy,length=1195>
##contig=<ID=chrUn_JTFH01001605v1_decoy,length=1195>
##contig=<ID=chrUn_JTFH01000760v1_decoy,length=1194>
##contig=<ID=chrUn_JTFH01000759v1_decoy,length=1194>
##contig=<ID=chrUn_JTFH01001606v1_decoy,length=1194>
##contig=<ID=chrUn_JTFH01001157v1_decoy,length=1193>
##contig=<ID=chrUn_JTFH01001607v1_decoy,length=1191>
##contig=<ID=chrUn_JTFH01001158v1_decoy,length=1191>
##contig=<ID=chrUn_JTFH01000761v1_decoy,length=1191>
##contig=<ID=chrUn_JTFH01001608v1_decoy,length=1189>
##contig=<ID=chrUn_JTFH01000762v1_decoy,length=1189>
##contig=<ID=chrUn_JTFH01001609v1_decoy,length=1188>
##contig=<ID=chrUn_JTFH01001159v1_decoy,length=1187>
##contig=<ID=chrUn_JTFH01001160v1_decoy,length=1186>
##contig=<ID=chrUn_KN707697v1_decoy,length=1186>
##contig=<ID=chrUn_JTFH01000763v1_decoy,length=1186>
##contig=<ID=chrUn_JTFH01000764v1_decoy,length=1186>
##contig=<ID=chrUn_KN707669v1_decoy,length=1184>
##contig=<ID=chrUn_JTFH01001162v1_decoy,length=1184>
##contig=<ID=chrUn_JTFH01000765v1_decoy,length=1184>
##contig=<ID=chrUn_JTFH01001161v1_decoy,length=1184>
##contig=<ID=chrUn_JTFH01000767v1_decoy,length=1183>
##contig=<ID=chrUn_JTFH01000766v1_decoy,length=1183>
##contig=<ID=chrUn_JTFH01001163v1_decoy,length=1182>
##contig=<ID=chrUn_JTFH01000768v1_decoy,length=1182>
##contig=<ID=chrUn_JTFH01000770v1_decoy,length=1181>
##contig=<ID=chrUn_JTFH01000772v1_decoy,length=1181>
##contig=<ID=chrUn_JTFH01000771v1_decoy,length=1181>
##contig=<ID=chrUn_JTFH01000769v1_decoy,length=1181>
##contig=<ID=chrUn_JTFH01001611v1_decoy,length=1180>
##contig=<ID=chrUn_JTFH01001610v1_decoy,length=1180>
##contig=<ID=chrUn_JTFH01000773v1_decoy,length=1179>
##contig=<ID=chrUn_KI270412v1,length=1179>
##contig=<ID=chrUn_JTFH01001612v1_decoy,length=1179>
##contig=<ID=chrUn_JTFH01001164v1_decoy,length=1179>
##contig=<ID=chrUn_JTFH01000775v1_decoy,length=1178>
##contig=<ID=chrUn_JTFH01000774v1_decoy,length=1178>
##contig=<ID=chrUn_KN707837v1_decoy,length=1178>
##contig=<ID=chrUn_JTFH01000776v1_decoy,length=1177>
##contig=<ID=chrUn_JTFH01000777v1_decoy,length=1177>
##contig=<ID=chrUn_KN707652v1_decoy,length=1176>
##contig=<ID=chrUn_JTFH01001165v1_decoy,length=1173>
##contig=<ID=chrUn_JTFH01001613v1_decoy,length=1172>
##contig=<ID=chrUn_JTFH01000778v1_decoy,length=1171>
##contig=<ID=chrUn_JTFH01000779v1_decoy,length=1171>
##contig=<ID=chrUn_JTFH01000780v1_decoy,length=1171>
##contig=<ID=chrUn_JTFH01000782v1_decoy,length=1170>
##contig=<ID=chrUn_JTFH01000781v1_decoy,length=1170>
##contig=<ID=chrUn_KN707941v1_decoy,length=1170>
##contig=<ID=chrUn_JTFH01001166v1_decoy,length=1169>
##contig=<ID=chrUn_KN707760v1_decoy,length=1169>
##contig=<ID=chrUn_JTFH01001614v1_decoy,length=1168>
##contig=<ID=chrUn_JTFH01000783v1_decoy,length=1167>
##contig=<ID=chrUn_JTFH01001167v1_decoy,length=1167>
##contig=<ID=chrUn_JTFH01000784v1_decoy,length=1167>
##contig=<ID=chrUn_JTFH01000785v1_decoy,length=1167>
##contig=<ID=chrUn_JTFH01001168v1_decoy,length=1166>
##contig=<ID=chrUn_JTFH01001615v1_decoy,length=1166>
##contig=<ID=chrUn_JTFH01000787v1_decoy,length=1165>
##contig=<ID=chrUn_JTFH01001169v1_decoy,length=1165>
##contig=<ID=chrUn_JTFH01000786v1_decoy,length=1165>
##contig=<ID=chrUn_JTFH01001170v1_decoy,length=1164>
##contig=<ID=chrUn_JTFH01001171v1_decoy,length=1163>
##contig=<ID=chrUn_KN707949v1_decoy,length=1162>
##contig=<ID=chrUn_JTFH01000788v1_decoy,length=1162>
##contig=<ID=chrUn_JTFH01001172v1_decoy,length=1158>
##contig=<ID=chrUn_JTFH01001173v1_decoy,length=1158>
##contig=<ID=chrUn_JTFH01001616v1_decoy,length=1157>
##contig=<ID=chrUn_JTFH01001175v1_decoy,length=1157>
##contig=<ID=chrUn_JTFH01001174v1_decoy,length=1157>
##contig=<ID=chrUn_JTFH01000789v1_decoy,length=1157>
##contig=<ID=chrUn_JTFH01001176v1_decoy,length=1157>
##contig=<ID=chrUn_JTFH01000790v1_decoy,length=1156>
##contig=<ID=chrUn_JTFH01001618v1_decoy,length=1156>
##contig=<ID=chrUn_JTFH01001617v1_decoy,length=1156>
##contig=<ID=chrUn_JTFH01000791v1_decoy,length=1156>
##contig=<ID=chrUn_JTFH01001177v1_decoy,length=1155>
##contig=<ID=chrUn_JTFH01001619v1_decoy,length=1155>
##contig=<ID=chrUn_JTFH01001178v1_decoy,length=1154>
##contig=<ID=chrUn_JTFH01000793v1_decoy,length=1154>
##contig=<ID=chrUn_JTFH01000792v1_decoy,length=1154>
##contig=<ID=chrUn_JTFH01001620v1_decoy,length=1154>
##contig=<ID=chrUn_JTFH01001621v1_decoy,length=1154>
##contig=<ID=chrUn_JTFH01000795v1_decoy,length=1151>
##contig=<ID=chrUn_JTFH01000794v1_decoy,length=1151>
##contig=<ID=chrUn_JTFH01000796v1_decoy,length=1150>
##contig=<ID=chrUn_JTFH01000797v1_decoy,length=1150>
##contig=<ID=chrUn_JTFH01001622v1_decoy,length=1149>
##contig=<ID=chrUn_JTFH01001179v1_decoy,length=1149>
##contig=<ID=chrUn_JTFH01001180v1_decoy,length=1148>
##contig=<ID=chrUn_JTFH01001181v1_decoy,length=1148>
##contig=<ID=chrUn_JTFH01000799v1_decoy,length=1147>
##contig=<ID=chrUn_JTFH01000798v1_decoy,length=1147>
##contig=<ID=chrUn_JTFH01001182v1_decoy,length=1146>
##contig=<ID=chrUn_JTFH01000800v1_decoy,length=1146>
##contig=<ID=chrUn_JTFH01001183v1_decoy,length=1144>
##contig=<ID=chrUn_JTFH01000801v1_decoy,length=1144>
##contig=<ID=chrUn_JTFH01000802v1_decoy,length=1144>
##contig=<ID=chrUn_KI270395v1,length=1143>
##contig=<ID=chrUn_JTFH01000803v1_decoy,length=1143>
##contig=<ID=chrUn_KN707891v1_decoy,length=1143>
##contig=<ID=chrUn_JTFH01001624v1_decoy,length=1143>
##contig=<ID=chrUn_JTFH01001623v1_decoy,length=1143>
##contig=<ID=chrUn_JTFH01000804v1_decoy,length=1142>
##contig=<ID=chrUn_JTFH01000806v1_decoy,length=1141>
##contig=<ID=chrUn_JTFH01000805v1_decoy,length=1141>
##contig=<ID=chrUn_JTFH01001184v1_decoy,length=1140>
##contig=<ID=chrUn_JTFH01000807v1_decoy,length=1140>
##contig=<ID=chrUn_JTFH01001625v1_decoy,length=1140>
##contig=<ID=chrUn_JTFH01000808v1_decoy,length=1138>
##contig=<ID=chrUn_JTFH01001626v1_decoy,length=1137>
##contig=<ID=chrUn_KN707687v1_decoy,length=1136>
##contig=<ID=chrUn_KI270376v1,length=1136>
##contig=<ID=chrUn_JTFH01001185v1_decoy,length=1136>
##contig=<ID=chrUn_JTFH01001629v1_decoy,length=1135>
##contig=<ID=chrUn_JTFH01001627v1_decoy,length=1135>
##contig=<ID=chrUn_JTFH01001628v1_decoy,length=1135>
##contig=<ID=chrUn_JTFH01001186v1_decoy,length=1134>
##contig=<ID=chrUn_JTFH01000809v1_decoy,length=1134>
##contig=<ID=chrUn_JTFH01000810v1_decoy,length=1134>
##contig=<ID=chrUn_JTFH01001187v1_decoy,length=1133>
##contig=<ID=chrUn_JTFH01000811v1_decoy,length=1132>
##contig=<ID=chrUn_JTFH01000813v1_decoy,length=1131>
##contig=<ID=chrUn_JTFH01000812v1_decoy,length=1131>
##contig=<ID=chrUn_KN707930v1_decoy,length=1131>
##contig=<ID=chrUn_JTFH01000814v1_decoy,length=1130>
##contig=<ID=chrUn_JTFH01001188v1_decoy,length=1129>
##contig=<ID=chrUn_JTFH01000815v1_decoy,length=1127>
##contig=<ID=chrUn_JTFH01001189v1_decoy,length=1127>
##contig=<ID=chrUn_JTFH01001190v1_decoy,length=1127>
##contig=<ID=chrUn_JTFH01001631v1_decoy,length=1127>
##contig=<ID=chrUn_JTFH01001630v1_decoy,length=1127>
##contig=<ID=chrUn_JTFH01000816v1_decoy,length=1126>
##contig=<ID=chrUn_KN707976v1_decoy,length=1126>
##contig=<ID=chrUn_JTFH01001632v1_decoy,length=1126>
##contig=<ID=chrUn_JTFH01000817v1_decoy,length=1124>
##contig=<ID=chrUn_KN707723v1_decoy,length=1124>
##contig=<ID=chrUn_JTFH01001633v1_decoy,length=1123>
##contig=<ID=chrUn_JTFH01001635v1_decoy,length=1123>
##contig=<ID=chrUn_JTFH01001634v1_decoy,length=1123>
##contig=<ID=chrUn_JTFH01000818v1_decoy,length=1122>
##contig=<ID=chrUn_JTFH01001637v1_decoy,length=1122>
##contig=<ID=chrUn_JTFH01000819v1_decoy,length=1122>
##contig=<ID=chrUn_JTFH01001636v1_decoy,length=1122>
##contig=<ID=chrUn_JTFH01001639v1_decoy,length=1121>
##contig=<ID=chrUn_JTFH01001638v1_decoy,length=1121>
##contig=<ID=chrUn_JTFH01000820v1_decoy,length=1121>
##contig=<ID=chrUn_KI270337v1,length=1121>
##contig=<ID=chrUn_KN707811v1_decoy,length=1120>
##contig=<ID=chrUn_JTFH01001642v1_decoy,length=1119>
##contig=<ID=chrUn_JTFH01000823v1_decoy,length=1119>
##contig=<ID=chrUn_JTFH01001640v1_decoy,length=1119>
##contig=<ID=chrUn_JTFH01001641v1_decoy,length=1119>
##contig=<ID=chrUn_JTFH01000824v1_decoy,length=1119>
##contig=<ID=chrUn_JTFH01000821v1_decoy,length=1119>
##contig=<ID=chrUn_JTFH01000822v1_decoy,length=1119>
##contig=<ID=chrUn_KN707951v1_decoy,length=1118>
##contig=<ID=chrUn_JTFH01001191v1_decoy,length=1118>
##contig=<ID=chrUn_JTFH01001643v1_decoy,length=1118>
##contig=<ID=chrUn_JTFH01000825v1_decoy,length=1118>
##contig=<ID=chrUn_KN707987v1_decoy,length=1117>
##contig=<ID=chrUn_JTFH01000827v1_decoy,length=1116>
##contig=<ID=chrUn_JTFH01000826v1_decoy,length=1116>
##contig=<ID=chrUn_JTFH01000828v1_decoy,length=1115>
##contig=<ID=chrUn_JTFH01000829v1_decoy,length=1115>
##contig=<ID=chrUn_JTFH01001644v1_decoy,length=1115>
##contig=<ID=chrUn_JTFH01000830v1_decoy,length=1115>
##contig=<ID=chrUn_JTFH01000831v1_decoy,length=1114>
##contig=<ID=chrUn_JTFH01000833v1_decoy,length=1113>
##contig=<ID=chrUn_JTFH01000832v1_decoy,length=1113>
##contig=<ID=chrUn_JTFH01000835v1_decoy,length=1110>
##contig=<ID=chrUn_JTFH01001192v1_decoy,length=1110>
##contig=<ID=chrUn_JTFH01000834v1_decoy,length=1110>
##contig=<ID=chrUn_JTFH01000836v1_decoy,length=1109>
##contig=<ID=chrUn_JTFH01000837v1_decoy,length=1108>
##contig=<ID=chrUn_JTFH01000840v1_decoy,length=1107>
##contig=<ID=chrUn_JTFH01000841v1_decoy,length=1107>
##contig=<ID=chrUn_JTFH01000838v1_decoy,length=1107>
##contig=<ID=chrUn_JTFH01000839v1_decoy,length=1107>
##contig=<ID=chrUn_JTFH01000842v1_decoy,length=1106>
##contig=<ID=chrUn_JTFH01001645v1_decoy,length=1106>
##contig=<ID=chrUn_JTFH01001646v1_decoy,length=1106>
##contig=<ID=chrUn_KN707721v1_decoy,length=1105>
##contig=<ID=chrUn_JTFH01001194v1_decoy,length=1104>
##contig=<ID=chrUn_JTFH01001647v1_decoy,length=1104>
##contig=<ID=chrUn_JTFH01001193v1_decoy,length=1104>
##contig=<ID=chrUn_KN707611v1_decoy,length=1103>
##contig=<ID=chrUn_JTFH01000844v1_decoy,length=1103>
##contig=<ID=chrUn_JTFH01000843v1_decoy,length=1103>
##contig=<ID=chrUn_JTFH01000845v1_decoy,length=1103>
##contig=<ID=chrUn_JTFH01001648v1_decoy,length=1102>
##contig=<ID=chrUn_JTFH01001649v1_decoy,length=1101>
##contig=<ID=chrUn_JTFH01001195v1_decoy,length=1101>
##contig=<ID=chrUn_KN707978v1_decoy,length=1101>
##contig=<ID=chrUn_KN707977v1_decoy,length=1101>
##contig=<ID=chrUn_JTFH01000846v1_decoy,length=1100>
##contig=<ID=chrUn_JTFH01000847v1_decoy,length=1099>
##contig=<ID=chrUn_KN707790v1_decoy,length=1098>
##contig=<ID=chrUn_JTFH01001196v1_decoy,length=1098>
##contig=<ID=chrUn_KN707910v1_decoy,length=1098>
##contig=<ID=chrUn_JTFH01000848v1_decoy,length=1098>
##contig=<ID=chrUn_JTFH01001651v1_decoy,length=1098>
##contig=<ID=chrUn_JTFH01001650v1_decoy,length=1098>
##contig=<ID=chrUn_JTFH01000849v1_decoy,length=1097>
##contig=<ID=chrUn_JTFH01000851v1_decoy,length=1096>
##contig=<ID=chrUn_JTFH01001652v1_decoy,length=1096>
##contig=<ID=chrUn_JTFH01001653v1_decoy,length=1096>
##contig=<ID=chrUn_JTFH01000850v1_decoy,length=1096>
##contig=<ID=chrUn_JTFH01001197v1_decoy,length=1096>
##contig=<ID=chrUn_JTFH01001654v1_decoy,length=1095>
##contig=<ID=chrUn_KN707950v1_decoy,length=1095>
##contig=<ID=chrUn_JTFH01000852v1_decoy,length=1094>
##contig=<ID=chrUn_JTFH01001198v1_decoy,length=1094>
##contig=<ID=chrUn_JTFH01001655v1_decoy,length=1093>
##contig=<ID=chrUn_JTFH01000853v1_decoy,length=1093>
##contig=<ID=chrUn_KN707735v1_decoy,length=1091>
##contig=<ID=chrUn_JTFH01001199v1_decoy,length=1091>
##contig=<ID=chrUn_JTFH01000854v1_decoy,length=1090>
##contig=<ID=chrUn_JTFH01001656v1_decoy,length=1090>
##contig=<ID=chrUn_JTFH01001657v1_decoy,length=1089>
##contig=<ID=chrUn_JTFH01001200v1_decoy,length=1089>
##contig=<ID=chrUn_KN707890v1_decoy,length=1088>
##contig=<ID=chrUn_JTFH01000855v1_decoy,length=1088>
##contig=<ID=chrUn_JTFH01001659v1_decoy,length=1087>
##contig=<ID=chrUn_KN707928v1_decoy,length=1087>
##contig=<ID=chrUn_JTFH01000856v1_decoy,length=1087>
##contig=<ID=chrUn_JTFH01001658v1_decoy,length=1087>
##contig=<ID=chrUn_JTFH01001201v1_decoy,length=1086>
##contig=<ID=chrUn_JTFH01000857v1_decoy,length=1086>
##contig=<ID=chrUn_JTFH01001661v1_decoy,length=1085>
##contig=<ID=chrUn_JTFH01000858v1_decoy,length=1085>
##contig=<ID=chrUn_JTFH01001662v1_decoy,length=1085>
##contig=<ID=chrUn_JTFH01001660v1_decoy,length=1085>
##contig=<ID=chrUn_JTFH01001202v1_decoy,length=1085>
##contig=<ID=chrUn_JTFH01000859v1_decoy,length=1084>
##contig=<ID=chrUn_KN707932v1_decoy,length=1084>
##contig=<ID=chrUn_JTFH01000862v1_decoy,length=1084>
##contig=<ID=chrUn_JTFH01001203v1_decoy,length=1084>
##contig=<ID=chrUn_JTFH01000860v1_decoy,length=1084>
##contig=<ID=chrUn_JTFH01000861v1_decoy,length=1084>
##contig=<ID=chrUn_JTFH01001205v1_decoy,length=1083>
##contig=<ID=chrUn_JTFH01001204v1_decoy,length=1083>
##contig=<ID=chrUn_JTFH01001663v1_decoy,length=1083>
##contig=<ID=chrUn_JTFH01000864v1_decoy,length=1083>
##contig=<ID=chrUn_JTFH01000863v1_decoy,length=1083>
##contig=<ID=chrUn_JTFH01000866v1_decoy,length=1082>
##contig=<ID=chrUn_JTFH01000865v1_decoy,length=1082>
##contig=<ID=chrUn_JTFH01000868v1_decoy,length=1081>
##contig=<ID=chrUn_JTFH01000867v1_decoy,length=1081>
##contig=<ID=chrUn_JTFH01001664v1_decoy,length=1080>
##contig=<ID=chrUn_JTFH01001665v1_decoy,length=1080>
##contig=<ID=chrUn_JTFH01001667v1_decoy,length=1079>
##contig=<ID=chrUn_JTFH01001206v1_decoy,length=1079>
##contig=<ID=chrUn_JTFH01001668v1_decoy,length=1079>
##contig=<ID=chrUn_KN707927v1_decoy,length=1079>
##contig=<ID=chrUn_JTFH01001666v1_decoy,length=1079>
##contig=<ID=chrUn_JTFH01000869v1_decoy,length=1079>
##contig=<ID=chrUn_JTFH01001207v1_decoy,length=1076>
##contig=<ID=chrUn_JTFH01000870v1_decoy,length=1076>
##contig=<ID=chrUn_JTFH01001669v1_decoy,length=1075>
##contig=<ID=chrUn_JTFH01001670v1_decoy,length=1074>
##contig=<ID=chrUn_JTFH01000871v1_decoy,length=1074>
##contig=<ID=chrUn_JTFH01000872v1_decoy,length=1073>
##contig=<ID=chrUn_JTFH01000873v1_decoy,length=1073>
##contig=<ID=chrUn_JTFH01001671v1_decoy,length=1073>
##contig=<ID=chrUn_KN707945v1_decoy,length=1072>
##contig=<ID=chrUn_JTFH01000874v1_decoy,length=1071>
##contig=<ID=chrUn_KN707934v1_decoy,length=1070>
##contig=<ID=chrUn_KN707645v1_decoy,length=1070>
##contig=<ID=chrUn_JTFH01001672v1_decoy,length=1070>
##contig=<ID=chrUn_JTFH01001208v1_decoy,length=1069>
##contig=<ID=chrUn_JTFH01000875v1_decoy,length=1069>
##contig=<ID=chrUn_JTFH01001673v1_decoy,length=1068>
##contig=<ID=chrUn_JTFH01001209v1_decoy,length=1068>
##contig=<ID=chrUn_KN707657v1_decoy,length=1068>
##contig=<ID=chrUn_JTFH01000877v1_decoy,length=1067>
##contig=<ID=chrUn_JTFH01000878v1_decoy,length=1067>
##contig=<ID=chrUn_JTFH01000876v1_decoy,length=1067>
##contig=<ID=chrUn_JTFH01001674v1_decoy,length=1067>
##contig=<ID=chrUn_JTFH01001212v1_decoy,length=1067>
##contig=<ID=chrUn_JTFH01001211v1_decoy,length=1067>
##contig=<ID=chrUn_JTFH01001210v1_decoy,length=1067>
##contig=<ID=chrUn_JTFH01001676v1_decoy,length=1066>
##contig=<ID=chrUn_JTFH01001677v1_decoy,length=1066>
##contig=<ID=chrUn_JTFH01000879v1_decoy,length=1066>
##contig=<ID=chrUn_JTFH01001675v1_decoy,length=1066>
##contig=<ID=chrUn_JTFH01000884v1_decoy,length=1065>
##contig=<ID=chrUn_JTFH01000880v1_decoy,length=1065>
##contig=<ID=chrUn_JTFH01000882v1_decoy,length=1065>
##contig=<ID=chrUn_JTFH01000881v1_decoy,length=1065>
##contig=<ID=chrUn_JTFH01000883v1_decoy,length=1065>
##contig=<ID=chrUn_KN707663v1_decoy,length=1065>
##contig=<ID=chrUn_JTFH01000886v1_decoy,length=1064>
##contig=<ID=chrUn_KN707834v1_decoy,length=1064>
##contig=<ID=chrUn_JTFH01000887v1_decoy,length=1064>
##contig=<ID=chrUn_JTFH01000885v1_decoy,length=1064>
##contig=<ID=chrUn_JTFH01001678v1_decoy,length=1063>
##contig=<ID=chrUn_JTFH01001679v1_decoy,length=1063>
##contig=<ID=chrUn_JTFH01001213v1_decoy,length=1063>
##contig=<ID=chrUn_JTFH01000888v1_decoy,length=1063>
##contig=<ID=chrUn_JTFH01001680v1_decoy,length=1063>
##contig=<ID=chrUn_JTFH01000890v1_decoy,length=1062>
##contig=<ID=chrUn_JTFH01000891v1_decoy,length=1062>
##contig=<ID=chrUn_JTFH01001214v1_decoy,length=1062>
##contig=<ID=chrUn_JTFH01001681v1_decoy,length=1062>
##contig=<ID=chrUn_JTFH01000889v1_decoy,length=1062>
##contig=<ID=chrUn_JTFH01000892v1_decoy,length=1061>
##contig=<ID=chrUn_JTFH01000893v1_decoy,length=1060>
##contig=<ID=chrUn_JTFH01001215v1_decoy,length=1059>
##contig=<ID=chrUn_JTFH01001216v1_decoy,length=1058>
##contig=<ID=chrUn_JTFH01001682v1_decoy,length=1058>
##contig=<ID=chrUn_JTFH01001217v1_decoy,length=1058>
##contig=<ID=chrUn_JTFH01000894v1_decoy,length=1057>
##contig=<ID=chrUn_JTFH01000895v1_decoy,length=1057>
##contig=<ID=chrUn_JTFH01000896v1_decoy,length=1056>
##contig=<ID=chrUn_JTFH01001683v1_decoy,length=1056>
##contig=<ID=chrUn_JTFH01000897v1_decoy,length=1055>
##contig=<ID=chrUn_JTFH01001218v1_decoy,length=1055>
##contig=<ID=chrUn_JTFH01000899v1_decoy,length=1055>
##contig=<ID=chrUn_JTFH01000898v1_decoy,length=1055>
##contig=<ID=chrUn_JTFH01000900v1_decoy,length=1055>
##contig=<ID=chrUn_JTFH01001220v1_decoy,length=1054>
##contig=<ID=chrUn_JTFH01001219v1_decoy,length=1054>
##contig=<ID=chrUn_JTFH01000901v1_decoy,length=1054>
##contig=<ID=chrUn_JTFH01001221v1_decoy,length=1053>
##contig=<ID=chrUn_JTFH01001222v1_decoy,length=1053>
##contig=<ID=chrUn_JTFH01001223v1_decoy,length=1052>
##contig=<ID=chrUn_JTFH01001684v1_decoy,length=1052>
##contig=<ID=chrUn_JTFH01001224v1_decoy,length=1051>
##contig=<ID=chrUn_JTFH01001685v1_decoy,length=1051>
##contig=<ID=chrUn_JTFH01001686v1_decoy,length=1051>
##contig=<ID=chrUn_JTFH01000902v1_decoy,length=1051>
##contig=<ID=chrUn_JTFH01001687v1_decoy,length=1050>
##contig=<ID=chrUn_JTFH01000903v1_decoy,length=1050>
##contig=<ID=chrUn_JTFH01000904v1_decoy,length=1050>
##contig=<ID=chrUn_JTFH01001225v1_decoy,length=1049>
##contig=<ID=chrUn_JTFH01000905v1_decoy,length=1049>
##contig=<ID=chrUn_JTFH01001688v1_decoy,length=1048>
##contig=<ID=chrUn_KI270335v1,length=1048>
##contig=<ID=chrUn_KI270378v1,length=1048>
##contig=<ID=chrUn_JTFH01000906v1_decoy,length=1048>
##contig=<ID=chrUn_JTFH01001226v1_decoy,length=1047>
##contig=<ID=chrUn_JTFH01000907v1_decoy,length=1047>
##contig=<ID=chrUn_JTFH01000908v1_decoy,length=1046>
##contig=<ID=chrUn_JTFH01001689v1_decoy,length=1046>
##contig=<ID=chrUn_JTFH01001690v1_decoy,length=1046>
##contig=<ID=chrUn_JTFH01000910v1_decoy,length=1046>
##contig=<ID=chrUn_JTFH01000909v1_decoy,length=1046>
##contig=<ID=chrUn_KI270379v1,length=1045>
##contig=<ID=chrUn_JTFH01000913v1_decoy,length=1045>
##contig=<ID=chrUn_JTFH01000912v1_decoy,length=1045>
##contig=<ID=chrUn_JTFH01000911v1_decoy,length=1045>
##contig=<ID=chrUn_JTFH01001691v1_decoy,length=1045>
##contig=<ID=chrUn_JTFH01000914v1_decoy,length=1044>
##contig=<ID=chrUn_JTFH01001227v1_decoy,length=1044>
##contig=<ID=chrUn_JTFH01001228v1_decoy,length=1043>
##contig=<ID=chrUn_JTFH01001229v1_decoy,length=1043>
##contig=<ID=chrUn_JTFH01001692v1_decoy,length=1043>
##contig=<ID=chrUn_JTFH01001230v1_decoy,length=1042>
##contig=<ID=chrUn_JTFH01000915v1_decoy,length=1042>
##contig=<ID=chrUn_JTFH01001231v1_decoy,length=1042>
##contig=<ID=chrUn_JTFH01000916v1_decoy,length=1041>
##contig=<ID=chrUn_JTFH01001232v1_decoy,length=1041>
##contig=<ID=chrUn_KI270329v1,length=1040>
##contig=<ID=chrUn_JTFH01001233v1_decoy,length=1040>
##contig=<ID=chrUn_KN707612v1_decoy,length=1039>
##contig=<ID=chrUn_JTFH01000918v1_decoy,length=1039>
##contig=<ID=chrUn_JTFH01001234v1_decoy,length=1039>
##contig=<ID=chrUn_JTFH01000917v1_decoy,length=1039>
##contig=<ID=chrUn_JTFH01000919v1_decoy,length=1038>
##contig=<ID=chrUn_KN707839v1_decoy,length=1038>
##contig=<ID=chrUn_JTFH01001693v1_decoy,length=1038>
##contig=<ID=chrUn_JTFH01001235v1_decoy,length=1038>
##contig=<ID=chrUn_JTFH01001236v1_decoy,length=1037>
##contig=<ID=chrUn_JTFH01001237v1_decoy,length=1037>
##contig=<ID=chrUn_JTFH01000920v1_decoy,length=1036>
##contig=<ID=chrUn_JTFH01000921v1_decoy,length=1036>
##contig=<ID=chrUn_JTFH01001694v1_decoy,length=1036>
##contig=<ID=chrUn_JTFH01001697v1_decoy,length=1035>
##contig=<ID=chrUn_JTFH01000922v1_decoy,length=1035>
##contig=<ID=chrUn_JTFH01001695v1_decoy,length=1035>
##contig=<ID=chrUn_JTFH01000923v1_decoy,length=1035>
##contig=<ID=chrUn_JTFH01001696v1_decoy,length=1035>
##contig=<ID=chrUn_JTFH01001238v1_decoy,length=1035>
##contig=<ID=chrUn_JTFH01001698v1_decoy,length=1033>
##contig=<ID=chrUn_JTFH01000924v1_decoy,length=1033>
##contig=<ID=chrUn_JTFH01001699v1_decoy,length=1032>
##contig=<ID=chrUn_JTFH01000925v1_decoy,length=1032>
##contig=<ID=chrUn_JTFH01000928v1_decoy,length=1031>
##contig=<ID=chrUn_JTFH01000927v1_decoy,length=1031>
##contig=<ID=chrUn_JTFH01000926v1_decoy,length=1031>
##contig=<ID=chrUn_JTFH01001700v1_decoy,length=1031>
##contig=<ID=chrUn_KN707644v1_decoy,length=1030>
##contig=<ID=chrUn_KI270419v1,length=1029>
##contig=<ID=chrUn_JTFH01001239v1_decoy,length=1027>
##contig=<ID=chrUn_JTFH01000930v1_decoy,length=1027>
##contig=<ID=chrUn_JTFH01000929v1_decoy,length=1027>
##contig=<ID=chrUn_JTFH01001703v1_decoy,length=1026>
##contig=<ID=chrUn_JTFH01000931v1_decoy,length=1026>
##contig=<ID=chrUn_JTFH01001702v1_decoy,length=1026>
##contig=<ID=chrUn_KI270336v1,length=1026>
##contig=<ID=chrUn_JTFH01000932v1_decoy,length=1026>
##contig=<ID=chrUn_JTFH01001701v1_decoy,length=1026>
##contig=<ID=chrUn_JTFH01000933v1_decoy,length=1024>
##contig=<ID=chrUn_JTFH01000934v1_decoy,length=1024>
##contig=<ID=chrUn_JTFH01001704v1_decoy,length=1023>
##contig=<ID=chrUn_JTFH01001705v1_decoy,length=1022>
##contig=<ID=chrUn_JTFH01000936v1_decoy,length=1022>
##contig=<ID=chrUn_JTFH01000935v1_decoy,length=1022>
##contig=<ID=chrUn_JTFH01000937v1_decoy,length=1021>
##contig=<ID=chrUn_JTFH01001240v1_decoy,length=1021>
##contig=<ID=chrUn_JTFH01001241v1_decoy,length=1021>
##contig=<ID=chrUn_JTFH01001706v1_decoy,length=1020>
##contig=<ID=chrUn_JTFH01001707v1_decoy,length=1020>
##contig=<ID=chrUn_JTFH01001708v1_decoy,length=1020>
##contig=<ID=chrUn_JTFH01000938v1_decoy,length=1020>
##contig=<ID=chrUn_JTFH01001242v1_decoy,length=1019>
##contig=<ID=chrUn_JTFH01001709v1_decoy,length=1019>
##contig=<ID=chrUn_JTFH01001243v1_decoy,length=1019>
##contig=<ID=chrUn_JTFH01000939v1_decoy,length=1019>
##contig=<ID=chrUn_JTFH01000942v1_decoy,length=1018>
##contig=<ID=chrUn_JTFH01001710v1_decoy,length=1018>
##contig=<ID=chrUn_JTFH01000940v1_decoy,length=1018>
##contig=<ID=chrUn_JTFH01000941v1_decoy,length=1018>
##contig=<ID=chrUn_JTFH01001711v1_decoy,length=1018>
##contig=<ID=chrUn_KN707656v1_decoy,length=1017>
##contig=<ID=chrUn_JTFH01001712v1_decoy,length=1017>
##contig=<ID=chrUn_JTFH01000943v1_decoy,length=1016>
##contig=<ID=chrUn_JTFH01001244v1_decoy,length=1016>
##contig=<ID=chrUn_JTFH01001714v1_decoy,length=1015>
##contig=<ID=chrUn_JTFH01001713v1_decoy,length=1015>
##contig=<ID=chrUn_JTFH01001715v1_decoy,length=1015>
##contig=<ID=chrUn_JTFH01001716v1_decoy,length=1014>
##contig=<ID=chrUn_JTFH01001717v1_decoy,length=1014>
##contig=<ID=chrUn_JTFH01001245v1_decoy,length=1014>
##contig=<ID=chrUn_JTFH01001246v1_decoy,length=1013>
##contig=<ID=chrUn_JTFH01001720v1_decoy,length=1013>
##contig=<ID=chrUn_JTFH01001718v1_decoy,length=1013>
##contig=<ID=chrUn_JTFH01001719v1_decoy,length=1013>
##contig=<ID=chrUn_KN707870v1_decoy,length=1012>
##contig=<ID=chrUn_JTFH01001721v1_decoy,length=1012>
##contig=<ID=chrUn_JTFH01001722v1_decoy,length=1011>
##contig=<ID=chrUn_JTFH01001723v1_decoy,length=1011>
##contig=<ID=chrUn_JTFH01000945v1_decoy,length=1010>
##contig=<ID=chrUn_JTFH01000944v1_decoy,length=1010>
##contig=<ID=chrUn_KN707947v1_decoy,length=1010>
##contig=<ID=chrUn_JTFH01001724v1_decoy,length=1009>
##contig=<ID=chrUn_KN707989v1_decoy,length=1009>
##contig=<ID=chrUn_JTFH01000946v1_decoy,length=1009>
##contig=<ID=chrUn_JTFH01001247v1_decoy,length=1009>
##contig=<ID=chrUn_JTFH01001726v1_decoy,length=1008>
##contig=<ID=chrUn_JTFH01001725v1_decoy,length=1008>
##contig=<ID=chrUn_JTFH01001248v1_decoy,length=1008>
##contig=<ID=chrUn_JTFH01000947v1_decoy,length=1008>
##contig=<ID=chrUn_JTFH01001249v1_decoy,length=1007>
##contig=<ID=chrUn_JTFH01001729v1_decoy,length=1007>
##contig=<ID=chrUn_JTFH01001727v1_decoy,length=1007>
##contig=<ID=chrUn_JTFH01000948v1_decoy,length=1007>
##contig=<ID=chrUn_KN707658v1_decoy,length=1007>
##contig=<ID=chrUn_KN707634v1_decoy,length=1007>
##contig=<ID=chrUn_JTFH01001728v1_decoy,length=1007>
##contig=<ID=chrUn_JTFH01000949v1_decoy,length=1006>
##contig=<ID=chrUn_JTFH01001730v1_decoy,length=1006>
##contig=<ID=chrUn_KN707868v1_decoy,length=1005>
##contig=<ID=chrUn_JTFH01000951v1_decoy,length=1005>
##contig=<ID=chrUn_JTFH01000950v1_decoy,length=1005>
##contig=<ID=chrUn_JTFH01001731v1_decoy,length=1005>
##contig=<ID=chrUn_JTFH01001251v1_decoy,length=1004>
##contig=<ID=chrUn_JTFH01001250v1_decoy,length=1004>
##contig=<ID=chrUn_JTFH01000952v1_decoy,length=1004>
##contig=<ID=chrUn_JTFH01000953v1_decoy,length=1004>
##contig=<ID=chrUn_JTFH01001732v1_decoy,length=1003>
##contig=<ID=chrUn_JTFH01000957v1_decoy,length=1003>
##contig=<ID=chrUn_JTFH01001252v1_decoy,length=1003>
##contig=<ID=chrUn_JTFH01000954v1_decoy,length=1003>
##contig=<ID=chrUn_JTFH01000955v1_decoy,length=1003>
##contig=<ID=chrUn_JTFH01000956v1_decoy,length=1003>
##contig=<ID=chrUn_JTFH01000958v1_decoy,length=1002>
##contig=<ID=chrUn_JTFH01000959v1_decoy,length=1002>
##contig=<ID=chrUn_JTFH01001253v1_decoy,length=1001>
##contig=<ID=chrUn_JTFH01001733v1_decoy,length=1001>
##contig=<ID=chrUn_JTFH01001254v1_decoy,length=1000>
##contig=<ID=chrUn_JTFH01000960v1_decoy,length=1000>
##contig=<ID=chrUn_JTFH01001255v1_decoy,length=1000>
##contig=<ID=chrUn_JTFH01001734v1_decoy,length=1000>
##contig=<ID=chrUn_JTFH01000961v1_decoy,length=1000>
##contig=<ID=chrUn_JTFH01001256v1_decoy,length=1000>
##contig=<ID=chrUn_KI270312v1,length=998>
##contig=<ID=chrUn_KI270539v1,length=993>
##contig=<ID=chrUn_KI270385v1,length=990>
##contig=<ID=chrUn_KI270423v1,length=981>
##contig=<ID=chrUn_KI270392v1,length=971>
##contig=<ID=chrUn_KI270394v1,length=970>
##fileDate=$filedate
##reference=$opts{reference}
##source=gnomit-$version
HEADER

    my @vars = sort { $$infile_hash{$a}{chr} cmp $$infile_hash{$b}{chr} || $$infile_hash{$a}{start} <=> $$infile_hash{$b}{start} } keys %$infile_hash;
    my @samples = sort { $a cmp $b } keys %{ $data_hash{ $vars[0] } };
    print foutname1 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    foreach my $sample (@samples) {
        print foutname1 "\t$sample";
    }
    print foutname1 "\n";
    my $index = 1;
    foreach my $var (@vars) {
        my $chrom  = $$infile_hash{$var}{chr};
        my $id     = $$infile_hash{$var}{id};
        my $alt    = "<INS:MT>";
        my $qual   = $$infile_hash{$var}{qual};
        my $filter = $$infile_hash{$var}{filter};
        my %info   = ();
        my $pos;
        my $end;
        my $ciDelta;

        if ( defined( $$infile_hash{$var}{leftBkpt} ) ) {
            $pos     = $$infile_hash{$var}{leftBkpt} - 1;
            $end     = $$infile_hash{$var}{rightBkpt};
            $ciDelta = $end - $pos + 1;
        }
        else {
            $pos             = $$infile_hash{$var}{start} - 1;
            $end             = $$infile_hash{$var}{end};
            $ciDelta         = $end - $pos + 1;
            $info{IMPRECISE} = undef;
        }
        if ( defined($$infile_hash{$var}{mlen}) ){
            $info{MLEN} = $$infile_hash{$var}{mlen};
            $info{MSTART} = $$infile_hash{$var}{mstart};
            $info{MEND} = $$infile_hash{$var}{mend};
        }
        $info{SAMPLES} = $$infile_hash{$var}{samples};
        $info{CIPOS}  = "0,$ciDelta";
        $info{CIEND}  = "-$ciDelta,0";
        $info{END}    = $end;
        $info{SVTYPE} = "INS";

        my $refline = `$opts{samtools} faidx $opts{reference} $chrom:$pos-$pos`;
        my $ref = ( split( /\n/, $refline ) )[1];
        if ( !defined($ref) ) { $ref = "N"; }
        my $format = "GT:FT:GL0:GQ:PL";

        my $info = "";
        my @sKeys = sort { $a cmp $b } keys %info;
        for ( my $i = 0 ; $i <= $#sKeys ; $i++ ) {
            if ( $i > 0 ) { $info .= ";"; }
            if ( defined( $info{ $sKeys[$i] } ) ) {
                $info .= "$sKeys[$i]=$info{$sKeys[$i]}";
            }
            else {
                $info .= "$sKeys[$i]";
            }
        }
        print foutname1 "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format";
        foreach my $sample (@samples) {
            my @gls = ();
            my @pls = ();
            foreach my $geno ( 0 .. $opts{ploidy} ) {
                push @gls, sprintf( "%.2f", $$data_hash{$var}{$sample}{gl0}{$geno} );
                push @pls, $$data_hash{$var}{$sample}{pl}{$geno};
            }
            my $gl = join( ",", @gls );
            my $pl = join( ",", @pls );
            print foutname1 "\t$$data_hash{$var}{$sample}{gt}:$$data_hash{$var}{$sample}{ft}:$gl:$$data_hash{$var}{$sample}{gq}:$pl";
        }
        print foutname1 "\n";
        $index++;
    }
    close(foutname1);
    print "Exiting report()\n\n" if $opts{verbose};
}

sub calcGl {
    my ( $m, $g, $k, $l, $er, $ea ) = @_;
    
    print "in calcGl():\n"               if $opts{verbose};
    
    if (defined($ea)) {
        print "\t$m\t$g\t$k\t$l\t$er\t$ea\n" if $opts{verbose};
    }
    else {
        print "\t$m\t$g\t$k\t$l\t$er\tN/A\n" if $opts{verbose};
    }

    if ( 1 / $m**$k <= 0 ) { die "problem in calcGL 1, \t$m\t$g\t$k\t$l\t$er\t$ea\n"; }
    my $gl = log10( 1 / ( $m**$k ) );
    foreach my $e ( @{$er} ) {
        if ( ( ( $m - $g ) * $e ) + ( ( 1 - $e ) * $g ) <= 0 ) { die "problem in calcGL 2, \t$m\t$g\t$k\t$l\t$e\n"; }
        $gl += log10( ( ( $m - $g ) * $e ) + ( ( 1 - $e ) * $g ) );
    }
    if (defined($ea)) { 
        foreach my $e ( @{$ea} ) {
            if ( ( $m - $g ) * ( 1 - $e ) + ( $g * $e ) <= 0 ) { die "problem in calcGL 3, \t$m\t$g\t$k\t$l\t$e\n"; }
            $gl += log10( ( $m - $g ) * ( 1 - $e ) + ( $g * $e ) );
        }
    }
    return $gl;
}

sub log10 {
    my $n = shift;
    return log($n) / log(10);
}

sub getInput {
    my ( $infile_hash, $mask_hash, $sample_hash ) = @_;

    print "Entering getInput()\n" if $opts{verbose};

    #input sample information
    open( INFO, $opts{info_filename} ) || die "error opening file $opts{info_filename}, $!\n";
    my @header = ();
    my %fields = ();
    while (<INFO>) {
        chomp;
        if (/^sample/) { @header = split(/\t/); next; }
        my @row = split(/\t/);
        for ( my $i = 0 ; $i <= $#header ; $i++ ) {
            $fields{ $header[$i] } = $row[$i];
        }
        my $sample = $fields{sample};
        $$sample_hash{$sample}{pop}      = $fields{pop};
        $$sample_hash{$sample}{filename} = $fields{filename};
        if ( $fields{median_insert_size} ne "NA" ) {
            $$sample_hash{$sample}{winlen}                    = $fields{median_insert_size} + 3 * $fields{median_absolute_deviation};
            $$sample_hash{$sample}{median_insert_size}        = $fields{median_insert_size};
            $$sample_hash{$sample}{median_absolute_deviation} = $fields{median_absolute_deviation};
        }
        else {
            $$sample_hash{$sample}{winlen}                    = 500;
            $$sample_hash{$sample}{median_insert_size}        = 400;
            $$sample_hash{$sample}{median_absolute_deviation} = 40;
        }
        $$sample_hash{$sample}{meancoverage} = $fields{mean_coverage};
        if ( defined( $fields{read_groups} ) ) {
            my @rgs = split( /,/, $fields{read_groups} );
            %{ $$sample_hash{$sample}{read_groups} } = map { $_, 1 } @rgs;
        }
    }
    close INFO;

    #input mask coordinates
    open( MASK, $opts{mask_filename} ) || die "Could not open $opts{mask_filename} for input, $!\n";
    while (<MASK>) {
        chomp;
        my ( $chr, $start, $end, $id ) = split(/\t/);
        $$mask_hash{$chr}{$start} = $end;
    }
    close MASK;

    #input variant coordinates
    if ( $opts{input_filename} =~ /\.gz$/ ) {
        open( VARS, "zcat $opts{input_filename} |" ) || die "Could not open $opts{input_filename} for input, $!\n";
    }
    else {
        open( VARS, $opts{input_filename} ) || die "Could not open $opts{input_filename} for input, $!\n";
    }
    my $varnum = 1;

    #BED FORMAT
    #while (<VARS>) {
    #    chomp;
    #    my ( $chr, $start, $end ) = split(/\t/);
    #    $chr =~ s/chr//g;
    #    $$infile_hash{$varnum}{chr}   = $chr;
    #    $$infile_hash{$varnum}{start} = $start;
    #    $$infile_hash{$varnum}{end}   = $end;
    #    $varnum++;
    #}

    #VCF FORMAT
    while (<VARS>) {
        next if /^#/;
        chomp;
        my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info ) = split(/\t/);
        if ( defined( $opts{chr} ) ) {
            if ( $chr ne $opts{chr} ) { next; }
        }
        $infile_hash{$varnum}{chr}   = $chr;
        $infile_hash{$varnum}{id}    = $id;
        $infile_hash{$varnum}{start} = $pos + 1;
        my ($end) = $info =~ /END=(\d+)/;
        my ($mstart) = $info =~ /MSTART=(\d+)/;
        my ($mend) = $info =~ /MEND=(\d+)/;
        my ($mlen) = $info =~ /MLEN=(\d+)/;
        my ($samples) = $info =~ /SAMPLES=(\w+?);/;
        if (defined($mlen)) {
            $infile_hash{$varnum}{mlen} = $mlen;
            $infile_hash{$varnum}{mstart} = $mstart;
            $infile_hash{$varnum}{mend} = $mend;
        }
        $infile_hash{$varnum}{samples} = $samples;
        $infile_hash{$varnum}{end}    = $end;
        $infile_hash{$varnum}{filter} = $filter;
        $infile_hash{$varnum}{qual}   = $qual;
        $varnum++;
    }
    close VARS;
    print "Exiting getInput()\n\n" if $opts{verbose};
}

sub getMateInfo {
    my ( $qname, $rnext, $pnext, $readgroup_hash, $filename ) = @_;

    my $command = "";
    if ( $opts{by_chr_dir} ) {
        if ( $opts{ucsc} ) {
            $command = "$opts{samtools} view -T $opts{reference} $filename" . "chr$rnext.*cram chr$rnext:$pnext-$pnext |";
        }
        else {
            $command = "$opts{samtools} view -T $opts{reference} $filename$rnext.*cram $rnext:$pnext-$pnext |";
        }
    }
    else {
        if ( $opts{ucsc} ) {
            $command = "$opts{samtools} view -T $opts{reference} $filename chr$rnext:$pnext-$pnext |";
        }
        else {
            $command = "$opts{samtools} view -T $opts{reference} $filename $rnext:$pnext-$pnext |";
        }
    }
    open( MCI, "$command" ) || die "Could not open $filename, $!\n";

    my $cFlag    = 0;
    my $cPos     = -1;
    my $clipside = "n";
    my $clipsize = -1;
    my $clipseq = "";
    my $seqLen   = 0;
    my $seq      = "";
    my $matchLen = 0;

    while (<MCI>) {
        chomp;
        my ( $m_qname, $m_flag, $m_rname, $m_pos, $m_mapq, $m_cigar, $m_rnext, $m_pnext, $m_tlen, $m_seq, $m_qual, $opt ) = split(/\t/);
        if ( $m_qname ne $qname ) { next; }
        my ($read_group) = $_ =~ /RG:Z:(\S+)/;

        if ( $opts{read_groups} && !defined($read_group) ) { next; }
        elsif ( $opts{read_groups} && !defined( $$readgroup_hash{$read_group} ) ) { next; }

        ( $cFlag, $cPos, $clipside, $clipsize, $clipseq ) = getSoftClipInfo( $m_pos, $m_cigar, $m_qual, $m_seq );
        $seqLen   = 0;
        $matchLen = 0;
        $seq      = $m_seq;
        while ( $m_cigar =~ /(\d+)M/g ) {
            $seqLen   += $1;
            $matchLen += $1;
        }
        while ( $m_cigar =~ /(\d+)N/g ) {
            $seqLen += $1;
        }
        while ( $m_cigar =~ /(\d+)D/g ) {
            $seqLen += $1;
        }
    }
    close MCI;

    return ( $cFlag, $cPos, $clipside, $clipsize, $seqLen, $matchLen, $seq );
}

sub refineData {
    my ( $infile_hash, $mask_hash, $sample_hash, $data_hash ) = @_;
    print "Entering refineData()\n" if $opts{verbose};

    my $numSamp = 0;
    foreach my $var ( keys %{$infile_hash} ) {
        if ( !defined( $$infile_hash{$var}{leftBkpt} ) ) { next; }
        foreach my $sample ( keys %{$sample_hash} ) {
            $numSamp++;

            #last if $numSamp > 50;

            #currently basing coordinates off of ONLY start position, may need to revisit
            my $chr     = $$infile_hash{$var}{chr};
            my $l_start = $$infile_hash{$var}{leftBkpt} - $$sample_hash{$sample}{winlen};
            my $l_end   = $$infile_hash{$var}{leftBkpt};
            my $r_start = $$infile_hash{$var}{rightBkpt};
            my $r_end   = $$infile_hash{$var}{rightBkpt} + $$sample_hash{$sample}{winlen};
            my $command = "";
            my $chrM    = "";

            print "REFINED: $var\t$chr\t$l_start\t$l_end\t$r_start\t$r_end\n" if $opts{verbose};
            if ( $opts{by_chr_dir} ) {
                if ( $opts{ucsc} ) {
                    $command = "$opts{samtools} view -T $opts{reference} $$sample_hash{$sample}{filename}chr$chr.cram chr$chr:$l_start-$r_end |";
                    $chrM    = "M";
                }
                else {
                    $command = "$opts{samtools} view -T $opts{reference} $$sample_hash{$sample}{filename}$chr.*cram $chr:$l_start-$r_end |";
                    $chrM    = "chrM";
                }
            }
            else {
                if ( $opts{ucsc} ) {
                    $command = "$opts{samtools} view -T $opts{reference} $$sample_hash{$sample}{filename} chr$chr:$l_start-$r_end |";
                    $chrM    = "M";
                }
                else {
                    $command = "$opts{samtools} view -T $opts{reference} $$sample_hash{$sample}{filename} $chr:$l_start-$r_end |";
                    $chrM    = "chrM";
                }
            }
            $$data_hash{$var}{$sample}{numRefRP}  = 0;
            $$data_hash{$var}{$sample}{numAltRP}  = 0;
            $$data_hash{$var}{$sample}{numRefSR}  = 0;
            $$data_hash{$var}{$sample}{numAltSR}  = 0;
            $$data_hash{$var}{$sample}{qualRefRP} = ();
            $$data_hash{$var}{$sample}{qualAltRP} = ();
            $$data_hash{$var}{$sample}{qualRefSR} = ();
            $$data_hash{$var}{$sample}{qualAltSR} = ();
            $$data_hash{$var}{$sample}{avgQ}      = 0;

            my %found = ();

            #print "command: $command\n" if $opts{verbose};
            open( SAM, $command ) || die "error in opening file, $!\n";
            while (<SAM>) {
                chomp;
                my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split(/\t/);

                my $mapE = 10**( -1 * $mapq / 10 );
                my ($read_group) = $_ =~ /RG:Z:(\S+)/;

                if ( $opts{read_groups} && !defined($read_group) ) { next; }
                elsif ( $opts{read_groups} && !defined( $$sample_hash{$sample}{read_groups}{$read_group} ) ) { next; }

                if ( $mapq < $opts{min_map_qual} ) { next; }

                my $dir = 0;    #F
                if ( $flag & 16 ) { $dir = 1; }    #R
                my $dnext = 0;                     #mate F
                if ( $flag & 32 ) { $dnext = 1; }  #mate R

                my $seqLen   = 0;
                my $matchLen = 0;
                while ( $cigar =~ /(\d+)M/g ) {
                    $matchLen += $1;
                    $seqLen   += $1;
                }
                while ( $cigar =~ /(\d+)N/g ) {
                    $seqLen += $1;
                }
                while ( $cigar =~ /(\d+)D/g ) {
                    $seqLen += $1;
                }

                my ( $cFlag, $cPos, $clipside, $clipsize, $clipseq ) = getSoftClipInfo( $pos, $cigar, $qual, $seq );
                my ( $m_cFlag, $m_cPos, $m_clipside, $m_clipsize, $m_seqLen, $m_matchLen, $m_seq ) = getMateInfo( $qname, $rname, $pnext, $$sample_hash{$sample}{read_groups}, $$sample_hash{$sample}{filename} );

                if ($opts{breakpoint}) {
                    if ( $cFlag == 1 && ( abs( $cPos - $$infile_hash{$var}{leftBkpt} ) <= $opts{min_clipped_seq} || abs( $cPos - $$infile_hash{$var}{rightBkpt} ) <= $opts{min_clipped_seq} ) ) {
    
                        $$data_hash{$var}{$sample}{numAltSR}++;
                        push @{ $$data_hash{$var}{$sample}{qualAltSR} }, $mapE;
                    }
                    elsif ( $matchLen >= 0.95 * length($seq) && ( ( $pos < $$infile_hash{$var}{leftBkpt} && $pos + $seqLen - 1 > $$infile_hash{$var}{leftBkpt} && abs( $$infile_hash{$var}{leftBkpt} - $pos ) > $opts{min_clipped_seq} && abs( $$infile_hash{$var}{leftBkpt} - ( $pos + $seqLen - 1 ) ) > $opts{min_clipped_seq} ) || ( $pos < $$infile_hash{$var}{rightBkpt} && $pos + $seqLen - 1 > $$infile_hash{$var}{rightBkpt} && abs( $$infile_hash{$var}{rightBkpt} - $pos ) > $opts{min_clipped_seq} && abs( $$infile_hash{$var}{rightBkpt} - ( $pos + $seqLen - 1 ) ) > $opts{min_clipped_seq} ) ) ) {

                        #non-clipped reads with slight overlap with breakpoint may have mismatches instead of clipping, so skip
                        #also require at least 95% of read be 'matched' to reference, as the presence of insertion can cause some wacky mappings
                        $$data_hash{$var}{$sample}{numRefSR}++;
                        push @{ $$data_hash{$var}{$sample}{qualRefSR} }, $mapE;
                    }
                }

                if ( $rnext eq $chrM || checkMaskOverlap( $chr, $pos, $mask_hash ) ) {

                    #mate mapped to mt sequence
                    if ( ( $dir == 0 && $pos >= $l_start && $pos <= $l_end ) || ( $dir == 1 && $pos >= $r_start && $pos <= $r_end ) ) { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                }
                elsif ( $dir == 0 && $dnext == 1 && $rnext eq "=" && !$found{$qname} ) {
                    $found{$qname} = 1;
                    if ( abs($tlen) > $$sample_hash{$sample}{median_insert_size} + 3 * $$sample_hash{$sample}{median_absolute_deviation} ) { next; }    #insert length of potential reference supporting allele out of upper bounds
                    if ( abs($tlen) < $$sample_hash{$sample}{median_insert_size} - 3 * $$sample_hash{$sample}{median_absolute_deviation} ) { next; }    #insert length of potential reference supporting allele out of lower bounds

                    #allow some flexibility for clipped positions around breakpoint due to alignment artifacts (later may want to inplement a realignment step)
                    if    ( $cPos > -1   && abs( $cPos - $l_end ) <= $opts{min_clipped_seq}     && $cFlag == 1 )   { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                    elsif ( $m_cPos > -1 && abs( $m_cPos - $l_end ) <= $opts{min_clipped_seq}   && $m_cFlag == 1 ) { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                    elsif ( $cPos > -1   && abs( $cPos - $r_start ) <= $opts{min_clipped_seq}   && $cFlag == 1 )   { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                    elsif ( $m_cPos > -1 && abs( $m_cPos - $r_start ) <= $opts{min_clipped_seq} && $m_cFlag == 1 ) { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                    elsif ( ( $pos <= $l_end && $pnext + $m_seqLen >= $l_end ) || ( $pos <= $r_start && $pnext + $m_seqLen >= $r_start ) ) {

                        #non-clipped reads with slight overlap with breakpoint may have mismatches instead of clipping
                        if (   abs( $pos - $l_end ) > $opts{min_clipped_seq}
                            && abs( $pos + $seqLen - 1 - $l_end ) > $opts{min_clipped_seq}
                            && abs( $pnext - $l_end ) > $opts{min_clipped_seq}
                            && abs( $pnext + $m_seqLen - 1 - $l_end ) > $opts{min_clipped_seq}
                            && abs( $pos - $r_start ) > $opts{min_clipped_seq}
                            && abs( $pos + $seqLen - 1 - $r_start ) > $opts{min_clipped_seq}
                            && abs( $pnext - $r_start ) > $opts{min_clipped_seq}
                            && abs( $pnext + $m_seqLen - 1 - $r_start ) > $opts{min_clipped_seq}
                            && $matchLen >= 0.95 * length($seq)
                            && $m_matchLen >= 0.95 * length($m_seq) )
                        {
                            $data_hash{$var}{$sample}{numRefRP}++;
                            push @{ $$data_hash{$var}{$sample}{qualRefRP} }, $mapE;
                        }
                    }
                }
                $$data_hash{$var}{$sample}{avgQ} += $mapE;
            }
            close SAM;

            if ( $$data_hash{$var}{$sample}{numRefRP} + $$data_hash{$var}{$sample}{numAltRP} > 0 ) {
                $$data_hash{$var}{$sample}{avgQ} /= ( $$data_hash{$var}{$sample}{numRefRP} + $$data_hash{$var}{$sample}{numAltRP} );
            }
            else {
                $data_hash{$var}{$sample}{avgQ} = 0.00001;
            }
            print "$sample\t$var\t$$data_hash{$var}{$sample}{numRefRP}\t$$data_hash{$var}{$sample}{numRefSR}\t$$data_hash{$var}{$sample}{numAltRP}\t$$data_hash{$var}{$sample}{numAltSR}\t$$data_hash{$var}{$sample}{avgQ}\n" if $opts{verbose};
        }
    }
    print "Exiting refineData()\n\n" if $opts{verbose};
}

sub getData {
    my ( $infile_hash, $mask_hash, $sample_hash, $data_hash ) = @_;
    print "Entering getData()\n" if $opts{verbose};

    my $numSamp = 0;
    foreach my $var ( keys %{$infile_hash} ) {
        foreach my $sample ( keys %{$sample_hash} ) {
            $numSamp++;
            my $chr     = $$infile_hash{$var}{chr};
            my $l_start = $$infile_hash{$var}{start} - $$sample_hash{$sample}{winlen};
            my $l_end   = $$infile_hash{$var}{start};
            my $r_start = $$infile_hash{$var}{end};
            my $r_end   = $$infile_hash{$var}{end} + $$sample_hash{$sample}{winlen};
            my $command = "";
            my $chrM    = "";

            #print "$var\t$chr\t$l_start\t$l_end\t$r_start\t$r_end\n" if $opts{verbose};
            if ( $opts{by_chr_dir} ) {
                if ( $opts{ucsc} ) {
                    $command = "$opts{samtools} view -T $opts{reference} $$sample_hash{$sample}{filename}chr$chr.*cram chr$chr:$l_start-$r_end |";
                    $chrM    = "M";
                }
                else {
                    $command = "$opts{samtools} view -T $opts{reference} $$sample_hash{$sample}{filename}$chr.*cram $chr:$l_start-$r_end |";
                    $chrM    = "chrM";
                }
            }
            else {
                if ( $opts{ucsc} ) {
                    $command = "$opts{samtools} view -T $opts{reference} $$sample_hash{$sample}{filename} chr$chr:$l_start-$r_end |";
                    $chrM    = "M";
                }
                else {
                    $command = "$opts{samtools} view -T $opts{reference} $$sample_hash{$sample}{filename} $chr:$l_start-$r_end |";
                    $chrM    = "chrM";
                }
            }
            $$data_hash{$var}{$sample}{numRefRP}  = 0;
            $$data_hash{$var}{$sample}{numAltRP}  = 0;
            $$data_hash{$var}{$sample}{numRefSR}  = 0;
            $$data_hash{$var}{$sample}{numAltSR}  = 0;
            $$data_hash{$var}{$sample}{qualRefRP} = ();
            $$data_hash{$var}{$sample}{qualAltRP} = ();
            $$data_hash{$var}{$sample}{qualRefSR} = ();
            $$data_hash{$var}{$sample}{qualAltSR} = ();
            $$data_hash{$var}{$sample}{avgQ}      = 0.000001;

            #print "command: $command\n" if $opts{verbose};
            open( SAM, $command ) || die "error in opening file, $!\n";
            while (<SAM>) {
                chomp;
                my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split(/\t/);



                my $mapE = 10**( -1 * $mapq / 10 );
                my ($read_group) = $_ =~ /RG:Z:(\S+)/;

                if ( $opts{read_groups} && !defined($read_group) ) { next; }
                elsif ( $opts{read_groups} && !defined( $$sample_hash{$sample}{read_groups}{$read_group} ) ) { next; }

                if ( $mapq < $opts{min_map_qual} ) { next; }

                my $dir = 0;    #F
                if ( $flag & 16 ) { $dir = 1; }    #R
                my $dnext = 0;                     #mate F
                if ( $flag & 32 ) { $dnext = 1; }  #mate R

                my ( $cFlag, $cPos, $clipside, $clipsize, $clipseq ) = getSoftClipInfo( $pos, $cigar, $qual, $seq );
                if ( $cPos > -1 && $cFlag == 1 && $cPos >= $$infile_hash{$var}{start} - $opts{clipped_flank} && $cPos <= $$infile_hash{$var}{end} + $opts{clipped_flank} ) {

                    #potential breakpoints must be at or between previous breakpoint bounds
                    $$data_hash{$var}{$sample}{clipPos}{$cPos}{$cFlag}++;
                    if ( !defined( $$data_hash{$var}{$sample}{clipSeq}{$cPos} ) || length($clipseq) > length($$data_hash{$var}{$sample}{clipSeq}{$cPos}) ) {
                        $$data_hash{$var}{$sample}{clipSeq}{$cPos} = $clipseq;
                    }
                }
                if ( $rnext eq $chrM || checkMaskOverlap( $chr, $pos, $mask_hash ) ) {

                    #mate mapped to mt sequence
                    if ( ( $dir == 0 && $pos >= $l_start && $pos <= $l_end ) || ( $dir == 1 && $pos >= $r_start && $pos <= $r_end ) ) { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                }
                elsif ( $dir == 0 && $dnext == 1 && $rnext eq "=" ) {
                    if ( abs($tlen) > $$sample_hash{$sample}{median_insert_size} + 3 * $$sample_hash{$sample}{median_absolute_deviation} ) { next; }    #insert length of potential reference supporting allele out of upper bounds
                    if ( abs($tlen) < $$sample_hash{$sample}{median_insert_size} - 3 * $$sample_hash{$sample}{median_absolute_deviation} ) { next; }    #insert length of potential reference supporting allele out of lower bounds
                                                                                                                                                        #"normal" reads which clip at the breakpoint are -not- reference supporting!
                                                                                                                                                        # --------->..............<-------L------ or <------L///////
                    if ( $pos >= $l_start && $pos < $l_end && $pnext < $l_end && $pnext + length($seq) > $l_end ) {

                        #my ( $m_cFlag, $m_cPos, $m_clipside, $m_clipsize ) = getMateInfo( $qname, $rname, $pnext, $readgroup_hash, $$sample_hash{$sample}{filename});
                        #if   ( $m_cPos > -1 ) { $$data_hash{$var}{$sample}{numAltRP}++; }
                        #else                  { $$data_hash{$var}{$sample}{numRefRP}++; }
                    }

                    # --------->..............<-------R------ or <//////R-------
                    elsif ( $pos >= $l_start && $pos < $l_end && $pnext <= $r_end && $pnext + length($seq) > $r_end ) {

                        #my ( $m_cFlag, $m_cPos, $m_clipside, $m_clipsize ) = getMateInfo( $qname, $rname, $pnext, $readgroup_hash, $$sample_hash{$sample}{filename} );
                        #if   ( $m_cPos > -1 ) { $$data_hash{$var}{$sample}{numAltRP}++; }
                        #else                  { $$data_hash{$var}{$sample}{numRefRP}++; }
                    }

                    # --------L//////> or -------L----->..............<----------
                    elsif ( $pos >= $l_start && $pos < $l_end && $pnext > $l_end && $pos + length($seq) > $l_end ) {

                        #if   ( $cPos > -1 ) { $$data_hash{$var}{$sample}{numAltRP}++; }
                        #else                { $$data_hash{$var}{$sample}{numRefRP}++; }
                    }

                    # \\\\\\\\R------> or -------R----->..............<----------
                    elsif ( $pos <= $r_start && $pos + length($seq) > $r_start && $pnext > $r_start ) {

                        #if   ( $cPos > -1 ) { $$data_hash{$var}{$sample}{numAltRP}++; }
                        #else                { $$data_hash{$var}{$sample}{numRefRP}++; }
                    }

                    # ------------->......L......<------------- or ------------->......R......<-------------
                    elsif ( ( $pos >= $l_start && $pos < $l_end && $pnext > $l_end && $pnext <= $r_end ) || ( $pos >= $l_start && $pos < $r_start && $pnext > $r_start && $pnext <= $r_end ) ) {
                        $data_hash{$var}{$sample}{numRefRP}++;
                        push @{ $$data_hash{$var}{$sample}{qualRefRP} }, $mapE;
                    }
                }
                $$data_hash{$var}{$sample}{avgQ} += $mapE;
            }
            close SAM;

            if ( $$data_hash{$var}{$sample}{numRefRP} + $$data_hash{$var}{$sample}{numAltRP} > 0 ) {
                $$data_hash{$var}{$sample}{avgQ} /= ( $$data_hash{$var}{$sample}{numRefRP} + $$data_hash{$var}{$sample}{numAltRP} );
            }
            else {
                $data_hash{$var}{$sample}{avgQ} = 0.99;
            }
            print "\t$sample\t$var\t$$data_hash{$var}{$sample}{numRefRP}\t$$data_hash{$var}{$sample}{numAltRP}\t$$data_hash{$var}{$sample}{avgQ}\n" if $opts{verbose};
        }
    }
    print "Exiting getData()\n\n" if $opts{verbose};
}

sub checkMaskOverlap {
    my ( $chr, $pos, $mask_hash ) = @_;
    my $isMaskOverlap = 0;
    foreach my $maskStart ( keys %{ $$mask_hash{$chr} } ) {
        my $maskEnd = $$mask_hash{$chr}{$maskStart};
        if ( $pos >= $maskStart && $pos <= $maskEnd ) {
            $isMaskOverlap = 1;
            last;
        }
    }
    return $isMaskOverlap;
}

sub getSoftClipInfo {
    my ( $pos, $cigar, $qual, $seq ) = @_;
    my $clipside = "";
    my $clipsize = 0;
    my $clipseq  = "";
    my $cPos     = -1;
    my $avgQual  = -1;
    my $cFlag    = 1;

    if ( $cigar =~ /^(\d+)S.*M.*?(\d+)S$/ ) {
        if ( $1 > $2 ) {
            $cPos     = $pos;
            $clipside = "l";
            $clipsize = $1;

            #consider breakpoint after leftmost soft clipped fragment
        }
        else {
            $cPos     = $pos - 1;
            $clipside = "r";
            $clipsize = $2;
            while ( $cigar =~ /(\d+)M/g ) {    #have to take into account that a CIGAR may contain multiple M's
                $cPos += $1;
            }
            while ( $cigar =~ /(\d+)I/g ) {
                $cPos -= $1;
            }
            while ( $cigar =~ /(\d+)D/g ) {
                $cPos += $1;
            }
        }
    }

    #upstream soft clip only
    elsif ( $cigar =~ /^(\d+)S.*M/ ) {
        $cPos     = $pos;
        $clipside = "l";
        $clipsize = $1;
    }

    #downstream soft clip only
    elsif ( $cigar =~ /M.*?(\d+)S/ ) {
        $cPos     = $pos - 1;
        $clipside = "r";
        $clipsize = $1;
        while ( $cigar =~ /(\d+)M/g ) {    #have to take into account that a CIGAR may contain multiple M's
            $cPos += $1;
        }
        while ( $cigar =~ /(\d+)I/g ) {
            $cPos -= $1;
        }
        while ( $cigar =~ /(\d+)D/g ) {
            $cPos += $1;
        }
    }

    #Check quality of clipped sequence and alignment to reference
    if ( $cPos > -1 ) {
        my $clippedQuals = "";

        if ( $clipside eq "r" ) {
            $clipseq   = substr( $seq,  length($qual) - $clipsize - 1, $clipsize );
            $clippedQuals = substr( $qual, length($qual) - $clipsize - 1, $clipsize );
        }
        else {
            $clipseq   = substr( $seq,  0, $clipsize );
            $clippedQuals = substr( $qual, 0, $clipsize );
        }

        my $avgQualSum = 0;
        my $avgQualNum = 0;
        foreach my $qual ( split( //, $clippedQuals ) ) {
            $avgQualSum += ord($qual) - 33;
            $avgQualNum++;
        }
        $avgQual = $avgQualSum / $avgQualNum;
    }
    if ( $avgQual < 10 ) { $cFlag = 0; }
    if ( $clipsize < $opts{min_clipped_seq} ) { $cFlag = 0; }
    return ( $cFlag, $cPos, $clipside, $clipsize, $clipseq );
}

sub usage {
    my $version = shift;
    printf("\n");
    printf( "%-9s %s\n", "Program:", "gnomit.pl" );
    printf( "%-9s %s\n", "Version:", "$version" );
    printf("\n");
    printf( "%-9s %s\n", "Usage:", "gnomit.pl [options]" );
    printf("\n");
    printf( "%-9s %-35s %s\n", "Options:", "--input_filename=[filename]",     "Input alignment file in BAM format" );
    printf( "%-9s %-35s %s\n", "",         "--info_filename=[filename]",      "Input file wth per-sample information (required)" );
    printf( "%-9s %-35s %s\n", "",         "--output_filename=[filename]",    "Output file (default stdout)" );
    printf( "%-9s %-35s %s\n", "",         "--mask_filename=[filename]",      "Mask file for reference numts in BED format (optional)" );
    printf( "%-9s %-35s %s\n", "",         "--reference=[filename]",          "Reference file");
    printf( "%-9s %-35s %s\n", "",         "--include_mask",                  "Include aberrant reads mapped to mask regions in clustering" );
    printf( "%-9s %-35s %s\n", "",         "--breakpoint",                    "Include soft clipped reads in likelihood calculation" );
    printf( "%-9s %-35s %s\n", "",         "--len_cluster_include=[integer]", "Maximum distance to be included in cluster (default 600)" );
    printf( "%-9s %-35s %s\n", "",         "--len_cluster_link=[integer]",    "Maximum distance to link clusters (default 800)" );
    printf( "%-9s %-35s %s\n", "",         "--min_reads_cluster=[integer]",   "Minimum number of reads to link a cluster (default 1)" );
    printf( "%-9s %-35s %s\n", "",         "--min_evidence=[integer]",        "Minimum evidence to consider an insertion event for genotyping (default 3)" );
    printf( "%-9s %-35s %s\n", "",         "--min_map_qual=[integer]",        "Minimum mapping quality for read consideration (default 10)" );
    printf( "%-9s %-35s %s\n", "",         "--max_read_cov=[integer]",        "Maximum read coverage allowed for breakpoint searching (default 200)" );
    printf( "%-9s %-35s %s\n", "",         "--min_clipped_seq=[integer]",     "Minimum clipped sequence required to consider as putative breakpoint (default 5)" );
    printf( "%-9s %-35s %s\n", "",         "--max_num_clipped=[integer]",     "Maximum number of clipped sequences observed before removing from evidence consideration (default 5)" );
    printf( "%-9s %-35s %s\n", "",         "--read_groups",                   "Stratify analysis to specified read group(s) indicated in info_filename (optional)" );
    printf( "%-9s %-35s %s\n", "",         "--use_priors",                    "Estimate priors using EM framework" );
    printf( "%-9s %-35s %s\n", "",         "--by_chr_dir",                    "If set, expects to find chr specific BAM files info_filename indicated directory" );
    printf( "%-9s %-35s %s\n", "",         "--prefix=[string]",               "Prepend label in report output" );
    printf( "%-9s %-35s %s\n", "",         "--ucsc",                          "Use UCSC genome formatting (e.g. chrM)" );
    printf("\n");
}

sub checkOptions {
    my $optResult = shift;
    my $opts      = shift;
    my $version   = shift;

    if ( !$optResult || $$opts{help} ) {
        usage($version);
        exit;
    }

    if ( !defined( $$opts{input_filename} ) && !defined( $$opts{by_chr_dir} ) ) {
        print "\n***ERROR***\t--input_filename or --by_chr_dir is required\n";
        usage($version);
        exit;
    }
    elsif ( !defined( $$opts{by_chr_dir} ) && !-e $$opts{input_filename} ) {
        print "\n***ERROR***\t--input_filename does not exist\n";
        usage($version);
        exit;
    }
    elsif ( defined( $$opts{by_chr_dir} ) && !-d $$opts{by_chr_dir} ) {
        print "\n***ERROR***\t--by_chr_dir does not exist\n";
        usage($version);
        exit;
    }
    if ( defined( $$opts{mask_filename} ) && !-e ( $$opts{mask_filename} ) ) {
        print "\n***ERROR***\t--mask_filename does not exist\n";
        usage($version);
        exit;
    }
    if ( !$$opts{include_mask} && !defined( $$opts{mask_filename} ) ) {
        print "\n***ERROR***\t--mask_filename is neccessary with --include_mask option\n";
        usage($version);
        exit;
    }
}

sub log2 {
    my $n = shift;
    return log($n) / log(2);
}
