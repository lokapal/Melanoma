#!/usr/bin/perl

    use strict;
    my $ARGC=scalar @ARGV;

    if ($ARGC!=2) { die "Usage: diffexp_select DESeq.full.results.tsv selectionfile\n"; }

    my $infile=$ARGV[0];
    my $listfile=$ARGV[1];

    open (INP, $infile) || die "Can't open input file \"$infile\": $!"; 
    my $basename=substr($infile,0,rindex($infile,"."));
    my $outfile=$basename . ".joined";
    open (OUTP,">$outfile") || die "Can't write to \"$outfile\": $!";

    <INP>; # Omit the first string
    my %expr;
    while (<INP>) {
       chomp;
       chomp;
       $_ =~ s/\"//g;
       my @arr=split('\t');
       my $id=$arr[1];
       my $log2FC=$arr[3];
       my $pvalue=$arr[6];
       my $padj=$arr[7];
       $expr{$id}{log2fc}=$log2FC;
       $expr{$id}{pvalue}=$pvalue;
       $expr{$id}{padj}=$padj;
                  }
    close (INP);

    open GENESF, $listfile or die $!;
    <GENESF>; # omit the first string
    print OUTP "GeneID\tlog2FC1\tpvalue1\tpadj1\tGene\tlog2FC2\tpvalue2\tpadj2\n";
    while (<GENESF>) {
        chomp;
        chomp;
        if (length($_)<5) { next; }
        my @arr=split('\t');
        my $id=$arr[0];
        if (length($id) != 15) { next; }
        if (! defined $expr{$id}) { print "Gene ID $id is missed?!\n"; next; } # exclude empty gene IDs
        my $log2FC=$expr{$id}{log2fc};
        my $pvalue=$expr{$id}{pvalue};
        my $padj=$expr{$id}{padj};
        print OUTP "$_\t$log2FC\t$pvalue\t$padj\n";
                     }

     close (GENESF);
     close (OUTP);
