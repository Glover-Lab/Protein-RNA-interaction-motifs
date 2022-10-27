#!/usr/bin/perl 
use Math::Trig;
#########################################################################
#	Program: Ncap_RNA.pl					       	#
#									#
#	Mark Glover							#
#	April 2021							#	
#									#	
#	This program scans a series of pdb files to find 		#
#	instances where the N-terminus of a helix interacts with        #
#       two consecutive phosphate groups in an RNA chain. The helix     #
#       must be capped with a Ser, the OG H-bonds to the mc N and also  #
#       to one of the phosphates. The N atoms of the the two residues   #
#       C-terminal to the Ser H-bond to the phosphate 3' to the one     #
#       that is recognized by the Ser. This arrangement is found in     #
#       in the RocC-RocR structure and we think it is conserved in the  #
#       FinO and ProQ mechanisms.
#                                                                       #
#       The algorithm searches for O1P, O2P, backbone N, and Ser OG     #
#       and stores those in 4 arrays. Next, it loops through these      #
#       arrays and makes new arrays for OG-OP1, OG-N, OP1-N, and OP2-N  #
#       H bonds where it stores the chain and resids of the bonded      #
#       atoms. Finally, it loops through the bond arrays to find        #
#       instances where all the H bonds are made. takes the 
#########################################################################

if ($#ARGV == -1 )
{
	print "\nUsage perl 3pocket.pl [file with pdb names - one per line] [H-bond cutoff distance] \n\n";
	exit;
}

print "****************************************************************\n\n";
print "				Ncap_RNA				\n\n";
print "		    	Mark Glover - April 2021  	      	\n\n";
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
 #Atom arrays:
	my @OP1; 	#array for atom lines for all OP1 atoms
        my @OP2; 	#array for atom lines for all OP2 atoms
        my @Natoms;     #array for atom lines of all mainchain N atoms
        my @OG;     #array for atom lines for all Ser/Thr OG atoms

 #H bond arrays:
	my @OP1_Nbond;   #array of chain and res IDs for OP1- N Hbonds
        my @OP2_Nbond;   #array of chain and res IDs for OP2- N Hbonds
        my @OP1_OGbond;  #array of chain and res IDs for OP1- OG Hbonds
        my @OG_Nbond;    #array of chain and res IDs for OG- N Hbonds
 
	#go through pdb and generate the OP1, OP2, OG, and mc N arrays
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

	 #if found OP1 in a nucleotide, push to @OP1
	 if (( $residue eq '  G' || $residue eq 'GUA' || $residue eq '  A' || $residue eq 'ADE' || $residue eq '  C' || $residue eq 'CYT' || $residue eq '  U' || $residue eq 'URA' || $residue eq "  T" || $residue eq "THY") && ($atom eq " OP1"))
		{
		    push(@OP1,$_);
		}
	 #if found OP2 in a nucleotide, push to @OP2
	 if (( $residue eq '  G' || $residue eq 'GUA' || $residue eq '  A' || $residue eq 'ADE' || $residue eq '  C' || $residue eq 'CYT' || $residue eq '  U' || $residue eq 'URA' || $residue eq "  T" || $residue eq "THY") && ($atom eq " OP2"))
		{
		    push(@OP2,$_);
		}
	  #if found N atom, push to @NHatoms
	 if ($atom eq ' N  ')
		{
		    push(@Natoms,$_);
 		}
	  #if found OG in ser, push to @Ser_OG
	 if (($residue eq "SER" || $residue eq "THR") && ($atom eq " OG "))

		{
		    push(@OG,$_);
 		}

	}

 
#go through OP1 array and test for distance to N
#first loop reads through OP1
	foreach (@OP1)
	{
	 #formatted read for ATOM/HETATM lines
		my $OP1residue = substr $_, 17, 3;
	 	my $OP1atom = substr $_, 12, 4;
		my $OP1chain = substr $_, 21, 1;
		my $OP1resnum = substr $_, 22, 4;
		my $OP1coordx = substr $_, 30, 8;
		my $OP1coordy = substr $_, 38, 8;
		my $OP1coordz = substr $_, 46, 8;
    	 
	 #now read through @Natoms and test for atoms H-bonding to OP1
         foreach (@Natoms)
         {
	 #formatted read for ATOM/HETATM lines
		my $Nresidue = substr $_, 17, 3;
	 	my $Natom = substr $_, 12, 4;
		my $Nchain = substr $_, 21, 1;
		my $Nresnum = substr $_, 22, 4;
		my $Ncoordx = substr $_, 30, 8;
		my $Ncoordy = substr $_, 38, 8;
		my $Ncoordz = substr $_, 46, 8;
    	 
	 #test for H-bonding
	 my $Hbond=sqrt(($Ncoordx-$OP1coordx)**2+($Ncoordy-$OP1coordy)**2+($Ncoordz-$OP1coordz)**2);

	 #add chain and resnums to OP1Nbond array if H-bonded
		    if ($Hbond < $ARGV[1])
		    {
			push (@OP1_Nbond, "$OP1chain $OP1resnum $Nchain $Nresnum");
		    }
	 }
	}

#go through OP2 array and test for distance to N
#first loop reads through OP2
	foreach (@OP2)
	{
	 #formatted read for ATOM/HETATM lines
		my $OP2residue = substr $_, 17, 3;
	 	my $OP2atom = substr $_, 12, 4;
		my $OP2chain = substr $_, 21, 1;
		my $OP2resnum = substr $_, 22, 4;
		my $OP2coordx = substr $_, 30, 8;
		my $OP2coordy = substr $_, 38, 8;
		my $OP2coordz = substr $_, 46, 8;
    	 
	 #now read through @Natoms and test for atoms H-bonding to OP2
         foreach (@Natoms)
         {
	 #formatted read for ATOM/HETATM lines
		my $Nresidue = substr $_, 17, 3;
	 	my $Natom = substr $_, 12, 4;
		my $Nchain = substr $_, 21, 1;
		my $Nresnum = substr $_, 22, 4;
		my $Ncoordx = substr $_, 30, 8;
		my $Ncoordy = substr $_, 38, 8;
		my $Ncoordz = substr $_, 46, 8;
    	 
	 #test for H-bonding
	 my $Hbond=sqrt(($Ncoordx-$OP2coordx)**2+($Ncoordy-$OP2coordy)**2+($Ncoordz-$OP2coordz)**2);

	 #add chain and resnums to OP2Nbond array if H-bonded
		    if ($Hbond < $ARGV[1])
		    {
			push (@OP2_Nbond, "$OP2chain $OP2resnum $Nchain $Nresnum");
		    }
	 }
	}
 
#go through OP1 array and test for distance to OG
#first loop reads through OP1
	foreach (@OP1)
	{
	 #formatted read for ATOM/HETATM lines
		my $OP1residue = substr $_, 17, 3;
	 	my $OP1atom = substr $_, 12, 4;
		my $OP1chain = substr $_, 21, 1;
		my $OP1resnum = substr $_, 22, 4;
		my $OP1coordx = substr $_, 30, 8;
		my $OP1coordy = substr $_, 38, 8;
		my $OP1coordz = substr $_, 46, 8;
    	 
	 #now read through @OG and test for atoms H-bonding to OP1
         foreach (@OG)
         {
	 #formatted read for ATOM/HETATM lines
		my $OGresidue = substr $_, 17, 3;
	 	my $OGatom = substr $_, 12, 4;
		my $OGchain = substr $_, 21, 1;
		my $OGresnum = substr $_, 22, 4;
		my $OGcoordx = substr $_, 30, 8;
		my $OGcoordy = substr $_, 38, 8;
		my $OGcoordz = substr $_, 46, 8;
    	 
	 #test for H-bonding
	 my $Hbond=sqrt(($OGcoordx-$OP1coordx)**2+($OGcoordy-$OP1coordy)**2+($OGcoordz-$OP1coordz)**2);

	 #add chain and resnums to OP1Nbond array if H-bonded
		    if ($Hbond < $ARGV[1])
		    {
			push (@OP1_OGbond, "$OP1chain $OP1resnum $OGchain $OGresnum");
		    }
	 }
	}

#go through OG array and test for distance to N
#first loop reads through OP1
	foreach (@OG)
	{
	 #formatted read for ATOM/HETATM lines
		my $OGresidue = substr $_, 17, 3;
	 	my $OGatom = substr $_, 12, 4;
		my $OGchain = substr $_, 21, 1;
		my $OGresnum = substr $_, 22, 4;
		my $OGcoordx = substr $_, 30, 8;
		my $OGcoordy = substr $_, 38, 8;
		my $OGcoordz = substr $_, 46, 8;
    	 
	 #now read through @Natoms and test for atoms H-bonding to OP1
         foreach (@Natoms)
         {
	 #formatted read for ATOM/HETATM lines
		my $Nresidue = substr $_, 17, 3;
	 	my $Natom = substr $_, 12, 4;
		my $Nchain = substr $_, 21, 1;
		my $Nresnum = substr $_, 22, 4;
		my $Ncoordx = substr $_, 30, 8;
		my $Ncoordy = substr $_, 38, 8;
		my $Ncoordz = substr $_, 46, 8;
    	 
	 #test for H-bonding
	 my $Hbond=sqrt(($Ncoordx-$OGcoordx)**2+($Ncoordy-$OGcoordy)**2+($Ncoordz-$OGcoordz)**2);

	 #add chain and resnums to OP1Nbond array if H-bonded
		    if ($Hbond < $ARGV[1])
		    {
			push (@OG_Nbond, "$OGchain $OGresnum $Nchain $Nresnum");
		    }
	 }
         }


 #now find residues that H-bond 2 consecutive NH's to OP1 & OP2, plus OG Ncap and OG - OP1
 	$file = $files[$currentpdb];
	chomp $file;

foreach (@OP1_OGbond)
{
  my @OP1_OGarray=split(" ",$_);
  foreach (@OG_Nbond)
  {
    my @OG_Narray=split(" ",$_);
    
    if (($OP1_OGarray[2] eq $OG_Narray[0]) && ($OP1_OGarray[3] == $OG_Narray[1]))
    {
      foreach (@OP2_Nbond)
      {
	my @OP2_Narray=split(" ",$_);
	if (($OP2_Narray[2] eq $OP1_OGarray[2]) && ($OP2_Narray[3] == $OG_Narray[1]+2) && ($OP2_Narray[1] == $OP1_OGarray[1]+1))
	{
	  foreach (@OP1_Nbond)
	  {
	    my @OP1_Narray=split(" ",$_);
	    if (($OP1_Narray[2] eq $OP1_OGarray[2]) && ($OP1_Narray[3] == $OG_Narray[1]+1) && ($OP1_Narray[0] eq $OP2_Narray[0]) && ($OP1_Narray[1] == $OP2_Narray[1]))
	    {
        	print "load $file \n";
	       	print "align (/$file and c;$OP1_Narray[0] and i;$OP1_OGarray[1]),(/RNA and i;1021)\n";
            }
	  }
	}
      }
    }
  }
}
 
 #clear all arrays/hashes for next PDB file
 @currentpdb=();
 @residues=();
 @OP1=();
 @OP2=();
 @Natoms=();
 @OG=();
 @OP1_Nbond=(); 
 @OP2_Nbond=(); 
 @OP1_OGbond=();
 @OG_Nbond=();  
 
}

