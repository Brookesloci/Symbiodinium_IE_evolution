$up5 = $down5 = $up3 = $down3 = 8;
open (IN, "SymbC1-SymbF.sharedintrons") ||
    open (IN, "SymbF-SymbC1.sharedintrons");
    
while (<IN>) {
    @F = split;
    push (@pairs, $_);
    $g{$F[0]}++;
    $g{$F[1]}++;
}
delete $g{"-"};


    
for $f ("SymbC1.intronsflanks", "SymbF.intronsflanks") {
    open (IN, $f);
    local $/ = ">";
    while (<IN>) {
	($n,$s) = /(\S+?):.*\n(\S+)/;
	$g{$n} || next;
	next unless $s =~ /([A-Z]{$up5})([a-z]{$down5})[a-z]+([a-z]{$up3})([A-Z]{$down3})/;
	$g{$n} = "$1|$2...$3|$4";
    }
}


open (IN, "SymbC1-SymbF.sharedintrons") ||
    open (IN, "SymbF-SymbC1.sharedintrons");

open (OUT, ">SymbC1-SymbF.sharedintrons.flanks");
while (<IN>) {
    @F = split;
    print OUT "$F[0]:$g{$F[0]}\t$F[1]:$g{$F[1]}\n";
}






open (IN, "SymbC1-SymbF.sharedintrons.flanks");
open (OUT, ">SymbC1-SymbF.ObservedConversions");
while (<IN>) {
    
    next if /\-:\s/;
    
    chomp;
    
    ($ib1,$ib2) = /\|([a-z]{2})/g;
    next if $ib1 eq $ib2;
    
    ($d1,$d2) = /\|([A-Z]{4})/g;
    next unless $d1 eq $d2;


    ($u1,$u2) = /\|([a-z]{4})/g;
    if ($u1 eq lc($d1)) {
	printf OUT "$_\t$ib1-->$ib2\n";
    }
    elsif ($u2 eq lc($d1)) {
	printf OUT "$_\t$ib2-->$ib1\n";
    }

}


    
    
    



