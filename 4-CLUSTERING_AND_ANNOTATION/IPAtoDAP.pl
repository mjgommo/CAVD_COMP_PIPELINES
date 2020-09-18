#!/usr/bin/perl -w

# For a description, see the subDie subroutine

if (@ARGV < 0) {
	subDie()
}

######################################

open TABLE, "/home/mjgomezr/DATABASES/IPA/IPA_ENSGtoNAME_20181127.txt" or die "Can not open table\n";
while ($line = <TABLE>) {
	chomp ($line);
	@line = split (/\t/,$line);
	$geneName = $line[1];
#	$geneName =~ s/\W/_/g;
	$geneNameTOgeneID{$geneName} = $line[0];
}
close TABLE;

######################################

open FILE, $ARGV[0] or die "Can not open $ARGV[0]\n";

print "Category\tTerm\tCount\tX.\tPValue\tGenes\tList.Total\tPop.Hits\tPop.Total\tFold.Enrichment\tBonferroni\tBenjamini\tFDR\n";

#<FILE>;
#<FILE>;
<FILE>;

#Ingenuity Canonical Pathways
# -log(B-H p-value)
#Ratio
#z-score
#Molecules

while ($line = <FILE>) {
	chomp ($line);
	@line = split (/\t/,$line);
	@molecules = split(",",$line[4]);
	@moleculeID = ();
	foreach $geneName (@molecules) {
		if (exists $geneNameTOgeneID{$geneName}) {
			$geneID = $geneNameTOgeneID{$geneName};
			push(@moleculeID,$geneID); 
		} else {
			$orphanGeneNames{$geneName} = 1;
		}
	}
	$count = @moleculeID;
	if ($count == 0) {
		next;
	}
	$term = $line[0];
	$term =~ s/\W+/_/g;
	$moleculeIDs = join(",",@moleculeID);
	$pval = 10 ** (-1 * $line[1]);
#	print "GOTERM_BP_DIRECT\t$term\t$count\t10\t$pval\t$moleculeIDs\t10\t10\t10\t10\t$pval\t$pval\t$pval\n";
	print "CANONICAL_PATHWAY\t$term\t$count\tNA\t$pval\t$moleculeIDs\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
}

close FILE;

open ORPHAN, ">>orphan_gene_names.txt";
foreach $geneName (keys %orphanGeneNames) {
	print ORPHAN "$geneName\t$ARGV[0]\n";
}
close ORPHAN;

#####################################################################


sub subDie {
die "
#####################################################################
#
# sqm02renameFiles.pl, copyright Manuel J Gomez
#
# You need X arguments.
#
# For example: 
#
# This script 
#
#####################################################################\n\n";
}

#####################################################################
