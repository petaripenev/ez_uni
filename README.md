# ez_uni.pl - Download rProtein fasta records

Script that accesses Entrez's webservice E-utilities for fetching sequences based on list of txids.
Requires perl and the LWP::Simple module installed.
WARNING: USES SEQUENCES FROM THE EFETCH RESULTS AS HASH IDS TO FILTER OUT REPEAT RESULTS
ALSO: REMOVES ANY RESULTS HAVING MITOCHONDRIAL, PLASTIDIC OR CHLOROPLAST IN THE NAME
It doesn't output clean fasta, since it adds the txids as lines preceding their sequences, for clarity.

Usage:

./ez_uni.pl -db [DATABASE] -prot [PROT NAME] -alt [FILTER PROT NAME] -in [FILE TXIDS] -cut -h [HELP]

-db	Specifies database from NCBI to be searched, I recommend protein. REQUIRED

-prot	Specifies protein name to be searched. Must be formatted as a query if more than one term. REQUIRED
	e.g. \"(term1+OR+term2+NOT+term3)\"
	
-in	File with txids for species, txids should be one per line. REQUIRED

-alt	Used only to filter results, case insensitive, not for actual searching; 
	adding this will decrease the number of your results. Must be guarded by single quotes - \'term1|term2\'
	
-cut	Cut results to one per txid, tries to remove Multispecies, putative, partial entries; takes first entry after sorting alphabetically.

-h	Displays this.
