$/ = "##";

#  This reads in the file that has all the annotated boundaries on the introners, and
#  and compiles them by phase and number of ss nucleotides recruited from the TSD.
#  Command: perl SCRIPT.ROYETAL.PHASERECRUITEDBYNUCS.PL Polarella_glacialis_CCMP1383.IntronersByGroup.BoundariesCalled >
#      
while (<>) {
    ($ic) = /(.+)/;
    while (/>\S+:ph(\d):.+\n.+?([a-z]{0,2})\|[A-Z]*[a-z]+[A-Z]*\|([a-z]{0,2})/g) {
	$tally{$ic}{length $2.$3}[$1]++;
	$tot{$ic}{length $2.$3}++;
    }
}
print "Group\tSSNucsRecruited\tNumberFrames\tPh0\tPh1\tPh2\tEffectivePhases\tNumberIntrons\tPassFilter\n";


#  This reports the data.
for $ic (sort  keys %tally) {
    $numb = 0;
    $tot = 0;
    for $length (keys %{$tally{$ic}}) {
	delete $tally{$ic}{$length} unless $tot{$ic}{$length} >= 50;
    }


    for $length (sort {$a<=>$b} keys %{$tally{$ic}}) {
	$filter =
	    ("N","Y")[int ($tot{$ic}{$length} >= 50)];
	next unless $tot{$ic}{$length} >= 50;
	printf "$ic\t$length\t%d", scalar keys %{$tally{$ic}};
	$effective_phases_inverse = 0;
	for $ph (0..2) {
	    printf "\t%.2f", $tally{$ic}{$length}[$ph]/$tot{$ic}{$length}*100;
	    $effective_phases_inverse += ($tally{$ic}{$length}[$ph]/$tot{$ic}{$length})**2;
	    
	}
	printf "\t%.2f\t$tot{$ic}{$length}\t$filter\n", 1/$effective_phases_inverse;
    }
}
