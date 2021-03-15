use English;
$/ = "##";

#  Reads in the training set, aligns each introner to the training set, and marks the
#  boundary of the actual transposable element (as opposed to the intron)
#  Input file should look like:
# ##Polarella_glacialis_CCMP1383.IntronerGroup.+6
# >Pg...evm.model.CCMP1383_scf7180000596960.5.8:ph0:20-104-20:PABCF-PABCF
# ACGGCAACACTTGGAATCACgagc|cttggtagtgtttaggcttacttttttacaatttgttttgttgatgtattgtcttagaataaagtttaatttggcttgggctgctgtaagcattcgtaag|GAGCACATTGCAACGCTGGG
# >Pg...evm.model.CCMP1383_scf7180000541519.3.7:ph2:20-103-20:PABCF-PABCF
# GTGGACTGCGTCCTTGCGTTgtgg|cttggtagtgtttaagctaacttctttactttttgtctggttaattattgtcgtagatgaagttcaatttgtttttgggcttttgtaagcacttgtaag|GTGGGGTCTTTGGTCTCCGT
# >Pg...evm.model.CCMP1383_scf7180000541114.1.45:ph1:20-104-20:P-PABCF
# CGCCTTTGTTTCCAACGTCAgaaacttggtagtgtttagacgtacttttttacattttgtcttgttaattcttgtcttagaaaaaagtttgaacttttttggggctgttgtaagcacttgtaagGAAATGAAGAAAAGCCGACG
    

<>;

#  For each set...
while ($data = <>) {
    chomp ($data);
    $_ = $data;
    s/.+//;
    $n = $&;
    next unless /\S+/;
    # comparisonset will hold the training set, done by hand...
    open (OUT, ">comparisonset");
    $code = 0;
    $trainingset = 0;
    while (/>(\S+)\s*\n(.+)/g) {
	$code[++$code] = $1;
	$backcode{$1} = $code;
	$s = $2;
	if ($s =~ /\|(.+)\|/) {
	    print OUT ">...$code\n$1\n";
	    $trainingset++;
	}
    }
    close OUT;
    next unless $trainingset;

    #  Now we go through and align each intron in the dataset to the comparison set, to give the best guess
    #  of the position of the edge of the TE.
    while ($data =~ />(\S+)\s*\n([A-Za-z\|]+)/g) {
	$code = $backcode{$1};
	$s = $2;
	$s =~ s/|//;
	system "cat comparisonset > temp";
	open (OUT, ">>temp");
	print OUT ">$code\n$s\n";
	close OUT;
	system "clustalw temp -outorder=input -gapext=0.1 -gapopen=30 -output=gde -quiet > /dev/null"; 
	%left = ();
	%right = ();
	{ local $/ = "#";
	  open (IN, "temp.gde");
	  while (<IN>) {
	      ($c) = /\.*(\S+)/;
	      $s = join ("", $POSTMATCH =~ /\S+/g);
	      if (/#$/) {
		  $s =~ /[a-z].+[a-z]/;
		  $left{length $PREMATCH}++;
		  $right{length $PREMATCH.$MATCH}++
	      }
	      else {
		  next unless /\S/;
		  $left = (sort {$left{$a}<=>$left{$b}} keys %left)[-1];
		  $right = (sort {$right{$a}<=>$right{$b}} keys %right)[-1];
		  $te = $right-$left;
		  ($c) = /(\S+)/;
		  $s = $POSTMATCH;
		  $s =~ s/\s//g;
		  $te || next;
		  $left || next;
		  
		  if (($te<32700)&&($left<32700)) {
		      $s =~ /(.{$left})(.{$te})(.*)/;
		      ($up,$insert,$down) = ($1,$2,$3);
		      $up[$c] = $up =~ tr/a-z/a-z/;
		      $down[$c] = $down =~ tr/a-z/a-z/;
		  }


		  }
	      }
	    }
	    

    }
    $data =~ /.+/g;
    print "##$&\n";
    while ($data =~ />(\S+)\s*\n(.+)/g) {
	$c = $backcode{$1};
	$s = $2;
	$s =~ s/[^A-Za-z]+//g;
	$s =~ s/([a-z]{$up[$c]})([a-z]+)([a-z]{$down[$c]})/$1|$2|$3/i;
	print ">$code[$c]\n$s\n";
    }
	      


}
