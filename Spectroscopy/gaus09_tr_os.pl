#!/usr/local/bin/perl
#;-*- Perl -*-ls
####################################################################
#  This script takes transition energies and oscillator strength from GAUSSIA out;
####################################################################

# Get input

        @args=@ARGV;
#       @args>=4 || die "usage: code.pl name v# ch# units#  \n";
        $name=$args[0];
        $lines=$args[1];
        #$Nst=$args[2];

        print "-------- input parameters ----------- \n";
        print "@args \n";
        print "name: $name \n"; 
        print "# lines: $lines \n";        
        #print "# states: $Nst \n";
        print "------------------------------------- \n";

        $OSout= "Etr_OS_". $name; 
#        $GAUSoutfile= $name . ".log";
         $GAUSoutfile= $name;
        print "working with gaussian out-file: $GAUSoutfile \n";    
#----------
# outfile:
#------------

        #$outfile="intCH_near_at".$pick_at;
        open (OUT0,"> OS_STRENGTH");

     system("grep -A $lines 'Excitation energies and oscillator strengths' $GAUSoutfile > $OSout" );
     system("grep 'Excited State ' $GAUSoutfile > test " );
  
# Find the number of lines in the  file test=#transitions:
        @NL=split(' ',`wc test`);
        $Los= $NL[0];
        print "total # lines in the test-file=# transitions = $Los \n";
        print "total # of transitions = $Los \n";        

        $LL=$Los+1;
# Take Transition Dipole moment:
        @TrD=split(' ',`grep -A $LL "transition electric dipole moments" $GAUSoutfile | tail -n $Los`);
        #print "TrD= @TrD \n";

# Take energies and Oscilator strength:
        @OS=split(' ',`grep 'Excited State ' $GAUSoutfile`);
        #print "OS= @OS \n";

        for ($q = 1; $q < $LL; $q++){
            $q1=10;
            $q2=100; 
            $x=($q-1)*6 + 1;
            $y=($q-1)*6 + 2;
            $z=($q-1)*6 + 3;
            $Dx=$TrD[$x];
            $Dy=$TrD[$y];
            $Dz=$TrD[$z];       
            
            $kE=($q-1)*10 + 4;
            $kf=($q-1)*10 + 8;
            $E=$OS[$kE];
            $F=$OS[$kf];
            @ff= split('=', $F); 
            #print "F:  $F  ff: $ff[1] \n";    
            $fff=$ff[1];


         printf OUT0 "%4d %4d %9.5f  %9.5f    %5.2f %5.2f    %9.5f %9.5f %9.5f \n", $q,$q, $fff,$E, $q1, $q2, $Dx,$Dy,$Dz;
         printf "%4d %4d %9.5f  %9.5f    %5.2f %5.2f    %9.5f %9.5f %9.5f \n", $q,$q, $fff,$E, $q1, $q2, $Dx,$Dy,$Dz;
        }

         print "state state   Osc.      E(eV)      dum   dum        TrDx      TrDy       TrDz \n";

     print "---------------------------------- \n";
     print "  $OSout and OS_STRENGTH are ready \n";
     print OUT0 "state state   Osc.      E(eV)      dum   dum        TrDx      TrDy       TrDz \n";

