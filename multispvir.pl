#!/usr/bin/perl
use strict;

my $printspvir = "y"; #If "n" it will only print out the total shared costs among all species, useful if simulating many species

my @ab = (4); #each parameter array must have the same number of items or it will return an error.
my @am = (1);
my @b0 = (0.5);
my @m0 = (0.1);
my $numsp = scalar(@ab);
if ((scalar(@am) != $numsp) or (scalar(@b0) != $numsp) or (scalar(@m0) != $numsp)){
	die "error: parameter arrays of differing length\n";
}
my $accuracy = 0.00005;
my $p = 0;
open (OUTPUT, '>multispvir data mut.txt');
#this creates a tab delineated file with headings
#all the heads should line up in excel unless there are multiple equilibria
#Because it is impossible to tell if there will be multiple equilibria it only prints heads for the first, stable equilibria
for (my $i = 0; $i < $numsp; $i++){
	print "ab$i,   am$i,   b0$i,    m0$i,    ";
	print OUTPUT "ab$i\tam$i\tb0$i\tm0$i\t";
}
print "p, pub_vir, N   ,";
print OUTPUT "p\tpub_vir\tN\t";
for (my $i = 0; $i < $numsp; $i++){
	print "v$i,   I$i,   "; 
	print OUTPUT "v$i\tI$i\tm$i\bI$i\t"; 	
}

print "\n";
print OUTPUT "\n";

for (my $ptemp = 0; $ptemp < 5 ;$ptemp += 0.01){#parameter to be incremented
	for (my $i = 0; $i < $numsp; $i++){ #this loop decreases $am as $p increases
		#$am[$i] = 1 - $ptemp;
	}
	$p=$ptemp;
	for (my $i = 0; $i < $numsp; $i++){
		if ($am[$i] eq "" or $ab[$i] eq "" or $b0[$i] eq "" or $m0[$i] eq ""){
			print "\nerror, parameters not defined for all species"; die;
		}
		printf("%.3f",$ab[$i]);	print ", ";
		printf("%.3f",$am[$i]);	print ", ";	
		printf("%.3f",$b0[$i]);	print ", ";	
		printf("%.3f",$m0[$i]);	print ", ";

		printf OUTPUT ("%.3f",$ab[$i]);	print OUTPUT  "\t";
		printf OUTPUT ("%.3f",$am[$i]);	print OUTPUT  "\t";	
		printf OUTPUT ("%.3f",$b0[$i]);print OUTPUT  "\t";	
		printf OUTPUT ("%.3f",$m0[$i]);print OUTPUT  "\t";
	}
	printf("%.3f",$p);	print ", ";
	printf OUTPUT ("%.3f",$p);	print OUTPUT  "\t";	
		
	my @stability = [];
	my $numunstable = 0;
	my $numstable = 0;
	my @I = [];
	my @v = [];	
	my @sh = [];
	my @shscalc = [];
	
	my $equilibrium = 0;
	
	my $error;		
	my $inbounds = "N";
	my $newerror = 0;
	#beginning search function, goes from 1 to some negative number
	for (my $shared = -.2; $shared < 10; $shared+=$accuracy){
		my @testv;
		my @testI;
		my $outofbounds = "n";
		$shscalc[$equilibrium] = 0;
		for (my $i = 0; $i < $numsp; $i++){
			$testv[$i]= (-2 * $am[$i] *  $b0[$i] + $ab[$i] * ($m0[$i]+$shared))/($ab[$i] * $am[$i]);
			my $m = ($am[$i] * $testv[$i] + $shared + $m0[$i])**2;
			my $b = $ab[$i] * $testv[$i] + $b0[$i];			

			if ($m<0 or $b<=0){
				$testv[$i]= -(($m0[$i]+$shared)/$am[$i]);
			}
			$b = $ab[$i] * $testv[$i] + $b0[$i];
			$m = ($am[$i] * $testv[$i] + $shared + $m0[$i])**2;		
			#print $ab[$i] . "\t" .  $testv[$i] . "\t" .  $b0[$i]. "\t" .  $shared. "\t" . $testv[$i] . "\t" . $m . "\t" .$b ."\n";

			$testI[$i] = 1 - $m / $b;
			if ($testI[$i] < 0){$testI[$i] = 0;}
			$shscalc[$equilibrium] += $testv[$i] * $testI[$i] / $numsp;
			$v[$equilibrium][$i] = $testv[$i];
			$I[$equilibrium][$i] = $testI[$i];	
					
			if ($I[$equilibrium][$i] < 0 or $I[$equilibrium][$i] > 1 or $m < 0 or $b < 0){
				$outofbounds = "y";
			}
		}
		#print $b . "\n";

			
		$shscalc[$equilibrium] = $shscalc[$equilibrium] * $p;		
		for (my $i = 0; $i < $numsp; $i++){
			my $b = $ab[$i] * $testv[$i] + $b0[$i];
			my $m = ($am[$i] * $testv[$i] + $shscalc[$equilibrium] + $m0[$i])**2;
			my $dI = (1-$I[$equilibrium][$i])*$b-$m;
			my $dIdv = (1-$I[$equilibrium][$i])*$ab[$i]-$am[$i];
			if (abs($dI) > 1000000 * $accuracy or abs($dIdv) > 1000000 * $accuracy){
				$outofbounds = "y"; 
			}				
		}
		$newerror =$shscalc[$equilibrium] - $shared;
		$sh[$equilibrium] = $shared;
		if ((($newerror > 0 and $error < 0) or ($newerror < 0 and $error > 0)) and $outofbounds eq "n"){
			if ($newerror < 0 and $error > 0){$stability[$equilibrium] = "S";$numstable++;}						
			if ($newerror > 0 and $error < 0){$stability[$equilibrium] = "U";$numunstable++;}						
			$equilibrium++;
		}
		$error = $newerror;
	}
	print $equilibrium . "\t";
	for (my $count = 0; $count < $equilibrium; $count++){
		printf("%.3f", $shscalc[$count]); print ", ";
		printf OUTPUT ("%.5f", $shscalc[$count]); print OUTPUT  "\t";		
		my $N = 1;
		printf OUTPUT ("%.5f", $N);print OUTPUT "\t";
		printf ("%.3f", $N);print "\t";		
		if ($printspvir eq "y"){
			for (my $i = 0; $i < $numsp; $i++){
				printf("%.3f",$v[$count][$i]); print ", ";
				printf("%.3f",$I[$count][$i]); print ", ";
			
				printf OUTPUT ("%.5f",$v[$count][$i]); print OUTPUT  "\t";			
				printf OUTPUT ("%.5f",$I[$count][$i]); print OUTPUT  "\t";	
				my $b = $ab[$i] * $v[$count][$i] + $b0[$i];
				my $m = ($am[$i] * $v[$count][$i] + $shscalc[$count] + $m0[$i])**2;
				printf OUTPUT ("%.5f",$m); print OUTPUT  "\t";	
				printf OUTPUT ("%.5f",$b); print OUTPUT  "\t";	
				printf("%.3f",$m); print ", ";
			}
			print $stability[$count] . ", ";	#"S" is a stable equilibrium, "U" is an unstable equilibrium.		
			print OUTPUT $stability[$count] . "\t";	
		}
	}
	print "\n";
	print OUTPUT "\n";
}
