#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple qw(get);
use HTTP::Request::Common qw(GET);
use LWP::UserAgent 6;
use feature qw(switch);
no warnings qw(experimental);
use Data::Dumper qw(Dumper);

#Help
sub help(){
print "DESCRIPTION:
Script that accesses Entrez's webservice E-utilities for fetching sequences based on list of txids.
Requires perl and the LWP::Simple module installed.
WARNING: USES SEQUENCES FROM THE EFETCH RESULTS AS HASH IDS TO FILTER OUT REPEAT RESULTS
ALSO: REMOVES ANY RESULTS HAVING MITOCHONDRIAL, PLASTIDIC OR CHLOROPLAST IN THE NAME
It doesn't output clean fasta, since it adds the txids as lines preceding the sequences for them, for clarity.
<Add option to stop doing that>
Usage:\n ./ez.pl -db [DATABASE] -prot [PROT NAME] -alt [FILTER PROT NAME] -in [FILE TXIDS] -h [HELP]
-db\tSpecifies database from NCBI to be searched, I reccomend protein. REQUIRED
-prot\tSpecifies protein name to be searched. Can be formatted as a query. REQUIRED
\te.g. \"(term1+OR+term2)\"
-in\tFile with txids for species, txids should be one per line. REQUIRED
-alt\tUsed only to filter results, case insensitive, not for actual searching; 
\tadding this will decrease the number of your results. Must be guarded by single quotes - \'term1|term2\'
-cut\tCut results to one per txid, tries to remove Multispecies, putative, partial entries; rest is random.
-h\tDisplays this.
"
}

#make it so that you'd filter your results by altnam around line 150 or output all




#Options
my ($db, $prot, $file, $altnam);
my $cut = 0;
for (my $i = 0; $i < @ARGV; $i++){
    given ($ARGV[$i]){
	when (/^\-db$/){
	    $i++;
	    $db = $ARGV[$i];
	}
	when (/^\-prot/){
	    $i++;
	    $prot = $ARGV[$i];
	}
	when (/^\-in/){
	    $i++;
	    $file = $ARGV[$i];
	}
	when (/^\-alt/){
	    $i++;
	    $altnam = $ARGV[$i];
	}
	when (/^\-cut/){
	    $i++;
	    $cut = 1;
	}
	when (/^\-h$/){
	    $i++;
	    help;
	    exit;
	}
	default{
	    print "Unrecognized argument ".$_."\n";
	    help;
	    exit;
	}
    }
}
	    
#Argument checks
unless (defined $db && defined $file && defined $prot){
    print STDERR "Missing required option.\n";
    help;
    exit;
}

unless (-e $file){
    print "Input file $file with txids doesn't exist!\n";
    help;
    exit;
}


my @ids;
open(IDS, "<$file");
while(<IDS>){
    push @ids, $_;
}
close IDS;

#chomp @ids;
#print "@ids";
#exit;

#$db = 'protein';
#Loops for each specie txid
my @query;
foreach (@ids){
	my ($base, $url, $web, $key);
    #for whatever reason chomp was misbehaving
    #chomp;
    $_ =~ s/\R//g;
    my $query = "txid".$_."+AND+".$prot;
	#Prints the txid for clarity
    print "%%%$_\n";
	#assemble the esearch URL
    $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

	#post the esearch URL
	my $ua = LWP::UserAgent->new( ssl_opts => { verify_hostname => 1 }, );
    my $req = GET $url;
	my $output = $ua->request($req);
	if ($output->is_success) {
    #print $output->content;
	} else {
	print $output->status_line . "\n";
	}

	#parse WebEnv and QueryKey
    $web = $1 if ($output->content =~ /<WebEnv>(\S+)<\/WebEnv>/);
    $key = $1 if ($output->content =~ /<QueryKey>(\d+)<\/QueryKey>/);
	#print "$web\n";
	
### include this code for ESearch-EFetch
	#assemble the efetch URL
    $url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
    $url .= "&rettype=fasta";
	
	####What the URL looks like
	#print "$url\n";
	#exit;
	####
	
	#post the efetch URL
    my $fetch = GET $url;
	my $data = $ua->request($fetch);
	if ($data->is_success) {
    #print $data->content;
	} else {
	print $data->status_line . "\n";
	}
	
	#Split data into array by entries
    my @seqs = split /[>]+/, $data->content;
	
	#Check for empty results and skip to next txid
	if ( grep( /\?xml/, @seqs ) ) {next;};
	#Clear empty elements in the array
	@seqs = grep /\S/, @seqs;
	my %out;
	#print Dumper \@seqs;
	#next;
	
	#If alternative name has been defined, filter results by it and by search parameter;
	#store them in a hash with keys being the sequences to remove same sequence entries, 
	#otherwise use only the search parameter for filtering. Also removes results from common 
	#organelles, which water down the eukaryotic results.
	foreach (@seqs){
		my @namseq;
		if (defined $altnam){
			if ($_ =~ /$altnam/i && $_ =~ /ribosomal/i && $_ !~ /mitochon/i && $_ !~ /plastid/i && $_ !~ /chloroplast/i){
				@namseq = split (/\n/, $_, 2);
				#Can be used for automatically cutting the names
				#$namseq[0] = "$1_$2" if $namseq[0] =~ /^(\w+).*\[(.*)\]/;
				$out{$namseq[1]}=$namseq[0];
			}
		}else{
			if ($_ =~ /ribosomal/i && $_ !~ /mitochon/i && $_ !~ /plastid/i && $_ !~ /chloroplast/i){
				@namseq = split (/\n/, $_, 2);
				#$namseq[0] = "$1_$2" if $namseq[0] =~ /^(\w+).*\[(.*)\]/;
				$out{$namseq[1]}=$namseq[0];
			}
		}
	}

	#Removing multiple hits for single organism and prioritizing
	#removal of ones that have Multispecies or putative in their names
	if ($cut == 1){
		foreach my $names (keys %out){
			if ($out{$names} =~ /MULTISPECIES/i || $out{$names} =~ /putative/i || $out{$names} =~ /partial/i){
				unless (keys %out < 2){
					delete $out{$names};
				}
			}
		}
		#Sorting and printing only once, ensures the entry with shortest name is outputted
		foreach my $names (sort { $out{$a} cmp $out{$b} } keys %out) {
			print ">$out{$names}\n$names";
			last;
		}
	}else{
		foreach my $names (keys %out) {
			print ">$out{$names}\n$names";
		}
	}
	#print Dumper \%out;
}

#For outputting fasta from the edited files
#sed 's/\%\%\%.*//' uL2.txt | sed 's/>.*\[/>uL2_/' | sed 's/\s/\_/g' | sed 's/\]//' | sed '/^$/d' | sed 's/candidatus_//i'
















