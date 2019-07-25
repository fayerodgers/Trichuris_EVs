#Use the Ensembl perl API to fetch the longest annotated 3' UTR for every mouse gene.

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Transcript;
use Bio::SeqIO;
use Data::Dumper;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);


my $gene_adaptor    = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
my @gene_objects = @{$gene_adaptor->fetch_all()};


while ( my $gene = shift @gene_objects ) {
	my $longest_utr;
	my @transcripts = @{ $gene->get_all_Transcripts()};	#retrieve all transcripts for the gene
	foreach my $transcript (@transcripts){
		next unless ($transcript->is_current() == 1) ;
		my $thr_utr = $transcript->three_prime_utr();
		next unless (defined $thr_utr);
		if (! defined $longest_utr){
			$longest_utr = $thr_utr->seq();
		}
		if (length($thr_utr->seq()) > length($longest_utr)){	#select the longest UTR for each gene 
			$longest_utr = $thr_utr->seq();
		}
	}
	next unless (defined $longest_utr);
	print ">", $gene->stable_id(), "\n",  $longest_utr, "\n";
}
