#!/usr/bin/perl 
use Math::Trig;
#########################################################################
#	Program: 3pocket.pl					       	#
#									#
#	Mark Glover							#
#	November 2020							#	
#									#	
#	This program scans a series of pdb files to find 		#
#	instances where the 2' and 3' OH of a nucleotide                #
#       H-bond with 2 consecutive backbone NH of a peptide.             #
#                                                                       #
#       The algorithm searches for O2' and O3' atoms in individual      #
#       residues, stores those in two arrays, and also searches for all #
#       backbone NHs and stores those in an array. Next, the O3' atoms  #
#       are scanned for H-bonding with any NHs. When a hit is found,    #
#       the O2' atoms are scanned for atoms in the same residue as the  #
#       O3' atom which also H bond with the NH in the next residue      #
#       from the NH previously identified. When a hit is found,         #
#       lines are printed for a .pml file to load the pdb file,         #
#       align the O2' residue on the 3' terminus of RocR, and show      #
#       lines for the nucleotide, as well as the 2 consecutive amino    #
#       acids in H-bonding distance.                                    #
#########################################################################

if ($#ARGV == -1 )
{
	print "\nUsage perl 3pocket.pl [file with pdb names - one per line] [H-bond cutoff distance] \n\n";
	exit;
}

print "****************************************************************\n\n";
print "				3pocket				\n\n";
print "		    	Mark Glover - November 2020	      	\n\n";
print "                                                                \n\n";
print "****************************************************************\n\n";

#read in file list
open(PDBFILES,"<$ARGV[0]");
@files=<PDBFILES>;

print "\nYour H-bonding distance: $ARGV[1]\n";

#for loop that opens each pdb file sequentially
for ($currentpdb=0; $currentpdb <= $#files; $currentpdb++)
{
	open(CURRENTPDBFILE,"<$files[$currentpdb]");
	my @currentpdb=<CURRENTPDBFILE>;
	close(CURRENTPDBFILE);
 
	#extract residues and put each atom/hetatm line into @residues
	my @residues;
	foreach (@currentpdb)
	{
		my @array=split(" ",$_);

		if ($array[0] eq "ATOM" || $array[0] eq "HETATM")
		{
		 push(@residues,$_);
		}
		@array=();
	}
 
	#variables and array definitions
	my @O2atoms; 	#array for atom lines for all O2' atoms
        my @O3atoms; 	#array for atom lines for all O3' atoms
        my @NHatoms;    #array for atom lines of all mainchain N atoms
	my @O3NHbond;   #array of chain and res IDs for O3'- NH Hbonds
        my @O2NHbond;   #array of chain and res IDs for O2'- NH Hbonds
 
	#go through pdb and generate the O2', O3' mc N arrays
	foreach (@residues)
	{
	 #formatted read for ATOM/HETATM lines
		my $residue = substr $_, 17, 3;
	 	my $atom = substr $_, 12, 4;
		my $chain = substr $_, 21, 1;
		my $resnum = substr $_, 22, 4;
		my $coordx = substr $_, 30, 8;
		my $coordy = substr $_, 38, 8;
	 my $coordz = substr $_, 46, 8;

	 #if found O2' in nuc, push to @O2atoms
	 if (($residue eq 'GMP' || $residue eq 'GDP' || $residue eq 'GTP' || $residue eq 'AMP' || $residue eq 'ADP' || $residue eq 'ATP' || $residue eq 'CMP' || $residue eq 'CDP' || $residue eq 'CTP' ||  $residue eq 'UMP' ||  $residue eq 'UDP' ||  $residue eq 'UTP' ||  $residue eq '  G' || $residue eq 'GUA' || $residue eq '  A' || $residue eq 'ADE' || $residue eq '  C' || $residue eq 'CYT' || $residue eq '  U' || $residue eq 'URA' || $residue eq "  T" || $residue eq "THY") && ($atom eq " O2'" || $atom eq ' O2*'))
		{
		    push(@O2atoms,$_);
		}
	  #if found O3' in nuc, push to @O3atoms
	 if (($residue eq 'GMP' || $residue eq 'GDP' || $residue eq 'GTP' || $residue eq 'AMP' || $residue eq 'ADP' || $residue eq 'ATP' || $residue eq 'CMP' || $residue eq 'CDP' || $residue eq 'CTP' ||  $residue eq 'UMP' ||  $residue eq 'UDP' ||  $residue eq 'UTP' ||  $residue eq '  G' || $residue eq 'GUA' || $residue eq '  A' || $residue eq 'ADE' || $residue eq '  C' || $residue eq 'CYT' || $residue eq '  U' || $residue eq 'URA' || $residue eq "  T" || $residue eq "THY") && ($atom eq " O3'" || $atom eq ' O3*'))

		{
		    push(@O3atoms,$_);
 		}
	  #if found N atom, push to @NHatoms
	 if ($atom eq ' N  ')
		{
		    push(@NHatoms,$_);
 		}
	}

 
#go through O3' array and test for distance to N
#first loop reads through O3'
	foreach (@O3atoms)
	{
	 #formatted read for ATOM/HETATM lines
		my $O3residue = substr $_, 17, 3;
	 	my $O3atom = substr $_, 12, 4;
		my $O3chain = substr $_, 21, 1;
		my $O3resnum = substr $_, 22, 4;
		my $O3coordx = substr $_, 30, 8;
		my $O3coordy = substr $_, 38, 8;
		my $O3coordz = substr $_, 46, 8;
    	 
	 #now read through @NHatoms and test for atoms H-bonding to O3'
         foreach (@NHatoms)
         {
	 #formatted read for ATOM/HETATM lines
		my $NHresidue = substr $_, 17, 3;
	 	my $NHatom = substr $_, 12, 4;
		my $NHchain = substr $_, 21, 1;
		my $NHresnum = substr $_, 22, 4;
		my $NHcoordx = substr $_, 30, 8;
		my $NHcoordy = substr $_, 38, 8;
		my $NHcoordz = substr $_, 46, 8;
    	 
	 #test for H-bonding
	 my $Hbond=sqrt(($NHcoordx-$O3coordx)**2+($NHcoordy-$O3coordy)**2+($NHcoordz-$O3coordz)**2);

	 #add chain and resnums to O3NHbond array if H-bonded
		    if ($Hbond < $ARGV[1])
		    {
			push (@O3NHbond, "$O3chain $O3resnum $NHchain $NHresnum");
		    }
	 }
	}
 

#go through O2' array and test for distance to N
#first loop reads through O2'
	foreach (@O2atoms)
         {
	 #formatted read for ATOM/HETATM lines
		my $O2residue = substr $_, 17, 3;
	 	my $O2atom = substr $_, 12, 4;
		my $O2chain = substr $_, 21, 1;
		my $O2resnum = substr $_, 22, 4;
		my $O2coordx = substr $_, 30, 8;
		my $O2coordy = substr $_, 38, 8;
		my $O2coordz = substr $_, 46, 8;
    	 
	 #now read through @NHatoms and test for atoms H-bonding to O2'
	 foreach (@NHatoms)
         {
	 #formatted read for ATOM/HETATM lines
		my $NHresidue = substr $_, 17, 3;
	 	my $NHatom = substr $_, 12, 4;
		my $NHchain = substr $_, 21, 1;
		my $NHresnum = substr $_, 22, 4;
		my $NHcoordx = substr $_, 30, 8;
		my $NHcoordy = substr $_, 38, 8;
		my $NHcoordz = substr $_, 46, 8;

	 #test for H-bonding
	 my $Hbond=sqrt(($NHcoordx-$O2coordx)**2+($NHcoordy-$O2coordy)**2+($NHcoordz-$O2coordz)**2);

	 #add chain and resnums to O2NHbond array if H-bonded
		    if ($Hbond < $ARGV[1])
		    {
			push (@O2NHbond, "$O2chain $O2resnum $NHchain $NHresnum");
		    }
	 }
	 }
	

 #now find residues that H-bond 2 consecutive NH's
 	$file = $files[$currentpdb];
	chomp $file;

foreach (@O3NHbond)
{
    my @O3array=split(" ",$_);
    foreach (@O2NHbond)
    {
	my @O2array=split(" ",$_);
    
	if (($O3array[0] eq $O2array[0]) && ($O3array[1] == $O2array[1]) && ($O3array[2] eq $O2array[2]) && ($O3array[3] == ($O2array[3]-1)))
	{
			
			print "load $file \n";
			print "align (/$file and c;$O3array[0] and i;$O3array[1]),(/RNA and i;1028)\n";
		        print "show lines, (/$file and ((c;$O3array[0] and i;$O3array[1]) or (c;$O3array[2] and i;$O3array[3]) or (c;$O2array[2] and i;$O2array[3]))\n";
	}
			
     }
}
 
 #clear all arrays/hashes for next PDB file
 @currentpdb=();
 @residues=();
 @O2atoms=();
 @O3atoms=();
 @NHatoms=();
 @O3NHbond=();
 @O2NHbond=();
}

