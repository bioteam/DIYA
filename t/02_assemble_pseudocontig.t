#-*-Perl-*-
# $Id: 02_assemble_pseudocontig.t 290 2008-12-17 19:41:58Z briano $
 
use strict;
use vars qw($NTESTS);

BEGIN: {
    $NTESTS = 9;

    eval { require Test::More;
		     require Test::Exception; };
    if ( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
	 use Test::Exception;

    plan tests => $NTESTS;
}

use_ok('Bio::DB::Fasta');
use_ok('Bio::SeqIO');

my $fasta_file = "t/data/little-buchnera.fa";
my $script = "scripts/diya-assemble_pseudocontig.pl";
my $tmp_dir = "t/tmp";
my $seqid = "test";
my $species = 'Buchnera';
my $strain = 'Buchnera';


SKIP: {

	# Simplest tests, no tiling file or scaffold file

	skip("$fasta_file not found. Skipping next tests.", 7) unless (-e $fasta_file);
	skip("$script not found. Skipping next tests.", 7) unless (-e $script);

	mkdir $tmp_dir unless (-d $tmp_dir);

	system "$script -infile $fasta_file -seqid $seqid -outdir $tmp_dir -species $species -strain $strain"; 

    for my $suffix (qw( fna fna.index gbk fasta )) {
        ok ( -e "$tmp_dir/$seqid.$suffix", "Found $tmp_dir/$seqid.$suffix" );
    }

	my $db = Bio::DB::Fasta->new("$tmp_dir/$seqid.fna");
	isa_ok ( $db, "Bio::DB::Fasta" );
	my $seq = $db->get_Seq_by_id("NC_004061-contig2");
	isa_ok ( $seq, "Bio::PrimarySeq::Fasta" );
	is ( $seq->length, 1200, "Correct sequence length" );

	my $in = Bio::SeqIO->new(-format => 'genbank',
									 -file => "$tmp_dir/$seqid.gbk");

}

END: {
	for my $suffix (qw( fasta fna fna.index gbk )) {
		unlink "$tmp_dir/$seqid.$suffix";
	}
}

1;
