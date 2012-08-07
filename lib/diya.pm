# --------------------------------------------------------------------------
# Copyright 2008
#
# This file is part of DIYA.
#
# DIYA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIYA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the diya package.  If not, see <http://www.gnu.org/licenses/>.
#--------------------------------------------------------------------------

=head1 NAME

DIYA - diya.pm

=head1 VERSION

1.0

=head1 SYNOPSIS

A simple diya script:

  use diya;

  $pipeline = diya->new;
  $pipeline->read_conf;
  $pipeline->run;

The script can be run like this:

  diya-script.pl --conf diya.conf seq1.fa seq2.fa ...

=head1 DESCRIPTION

I<diya> is an open source tool used to build annotation pipelines. A pipeline is
a series of steps linking the various stages of sequence annotation into a 
concise process. The software is designed to use sequences as input. These could 
be complete genomes or the result of shotgun sequencing of a genome library. 
A possible output would be a fully annotated genomic sequence in Genbank format.

You can also use this Genbank file as input and load GFF into a
backend database for viewing with tools like GBrowse.

A pipeline may be executed on a single computer or on a cluster. Currently I<diya> 
only supports the Sun Grid Engine platform if you are using a cluster. 

=head2 In a nutshell

All diya pipelines are made up of B<parser> or B<script> steps that are executed 
in a specific order. The details on the B<parser> and B<script> steps for a given 
pipeline are contained in a single XML configuration, or conf, file.

The C<diya.pm> Perl module is the controller module for a I<diya> annotation pipeline. 
This module reads the configuration file that describes the pipeline, executes each 
step in the pipeline, and launches specific B<parser> modules when required. 
It also keeps track of the input and output files and keeps all these files
in a single output directory.

=head2 B<parser> steps

A B<parser> step is what is doing the analysis in the pipeline. For every I<diya>
B<parser> step there will be a bioinformatics application that will produce
output and a corresponding Perl module that parses the
application output and creates an annotated Genbank file. 
A B<parser> step can act at any time in the pipeline.

=head2 B<script> steps

A B<script> step in I<diya> is simpler than a B<parser> step. Its output is not parsed so it does
not required a corresponding Perl module. A B<script> step may do something
like move a file, format a database, or send an alert. A B<script> step can
act at any time in the pipeline.

=head2 Code functionality

This is a rough description of how a pipeline works in the C<diya.pm> code:

=over 13

=item 1. A I<diya> script is launched from the command-line

=item 2. A pipeline object is created by new()

=item 3. Variables such as the input sequence files and the conf file name are set

=item 4. The I<diya> home directory is set using $DIYAHOME or the current working directory

=item 5. The pipeline object reads the configuration file with read_conf()

=item 6. The run() method is called, it will iterate over all the steps in the pipeline

=item 7. The pipeline object locates or creates the output directory

=item 8. The pipeline object sets the input sequence file

=item 9. The pipeline performs format conversion on the input file, if necessary

=item 10. A command string is constructed using the information in the configuration file

=item 11. The command is executed, creating a program output file

=item 12. The pipeline creates a parser object and the output file is parsed

=item 13. The pipeline proceeds to the next step

=item 14. If the pipeline finishs then run() starts again with the next sequence

=back

A good way to watch what I<diya> is doing is to run it with I<verbose> set
to 1. For example:

  diya.pl --conf diya.conf --verbose 1

=head1 INSTALLATION

The details are in the INSTALL file. I<diya> uses BioPerl, and you will need
to install some other Perl modules from CPAN in addition.

=head2 The $DIYAHOME variable

Consider setting the $DIYAHOME environment variable. By default 
I<diya> uses this directory when it looks for a I<diya> configuration file
and when it creates output directories. If you do not have this set then make
sure to tell I<diya> where your configuration file is using C<-conf> or
C<-use_conf>, see more about this below.

=head2 I<diya> tests

This package comes with a number of test scripts in the C<t/> directory that 
run automatically if you type:

  >perl Makefile.PL
  >make
  >make test

Most of the test scripts run bioinformatics applications, specifically
C<blastall>, C<formatdb>, C<tRNAscan-SE>, 
and C<glimmer3>. The scripts
are written such that they will skip many tests if these applications
are not found in C</usr/local/bin>. If you want to run these tests and
you have these applications installed then you may need to edit the *conf 
files found in C<t/data> to enter the correct paths.

=head1 THE CONFIGURATION FILE

Most of the information about the pipeline is stored in a
configuration file in XML format. The configuration file that comes with the 
package is called C<diya.conf> but you can create your own configuration files 
and call them whatever you want. There are also example *conf files in the 
C<t/data> and C<examples> directory in this package. 

The configuration file contains different sections. These sections can appear 
in any order in the file.

=head2 An example configuration file

 <?xml version="1.0" encoding="UTF-8"?>
 <conf>
   <script>
     <name>download</name>
     <executable>download-genome.pl</executable>
     <command>-id MYID -out OUTPUTFILE</command>
     <home>/Users/bosborne/diya/diya/branches/0.4.0/examples</home>
     <inputfrom></inputfrom>
   </script>
   <parser>
     <name>blastp</name>
     <executable>blastall</executable>
     <command>-i INPUTFILE -d MYDB -p blastp -o OUTPUTFILE</command>
     <home>/usr/local/bin</home>
     <inputformat>fasta</inputformat>
     <inputfrom>download</inputfrom>
   </parser>
   <run>
     <mode>serial</mode>
   </run>
   <order>
     <names>download blastp</names>
   </order>
 </conf>

You might run the pipeline using this configuration file like this:

 diya.pl -conf download-blastp.conf -set MYDB=/opt/gb/at.fa -set MYID=3

=head2 The B<order> section

This section tell I<diya> what B<script> and B<parser> steps to run and in what order.
The B<names> of the steps are separated by spaces. For example:

  <order>
    <names>tRNAscanSE glimmer blastall</names>
  </order>

You will see below that every B<parser> or B<script> section has a line
for its B<name>, like this:

  <name>tRNAscanSE</name>

You use that same B<name> in the B<names> line. This means that
there has to be a corresponding B<parser> or B<script> section for each
B<name> in the B<names> section.

=head2 The B<run> section

This section tells I<diya> whether to run the pipeline on a cluster or not. 
For example:

  <run>
    <mode>serial</mode>
  </run>

You can run I<diya> to run in B<sge> or B<serial> mode.
These are the only 2 possible values.

=head2 The B<parser> sections

This section describes a B<parser> step. An example for 
the application C<tRNAscan-SE>:

  <parser>
    <executable>tRNAscan-SE</executable>
    <home>/usr/local/bin</home>
    <command>-B -o OUTPUTFILE INPUTFILE</command>
    <name>tRNAscanSE</name>
    <inputformat>fasta</inputformat>
    <inputfrom></inputfrom>
  </parser>

=over 6

=item B<executable>

The name of the application. This should be the actual name, not a
synonym. Required.

=item B<home>

The directory where the application is found. Required.

=item B<command>

The command that has to be run, without the application name. Note that
these do not have to be real file names. Instead you can substitute actual 
input and output file names with I<INPUTFILE> and I<OUTPUTFILE>.
See more on this in L<WRITING YOUR OWN conf* FILES>.

=item B<name>

The arbitrary name for the step. It does not have to be the same as the
B<executable> but if the step is a B<parser> then this has to be the same 
as the name of the Perl module that parses the executable output. The only
rule is that no punctuation or spaces are allowed in the B<name>. For example,
a B<name> could be I<tRNAscanSE>, but not I<tRNAscan-SE> (the reason for this is
that spaces and punctuation are not allowed in a Perl module name).
Required.

In addition, you may want to have different steps in a pipeline that
use the same application or script, but in different ways. This way you can assign
a different B<name> to each of these steps.

=item B<inputformat>

The sequence format for the input file. Optional, if there is no B<inputformat> 
set then fasta format is assumed.

If B<inputformat> is set then I<diya> will determine the format of the input 
file for the given step. If this format is different from the B<inputformat> 
of the step then I<diya> will create a new file of the correct format and
make it the new input file for the step.

=item B<inputfrom>

This is optional. Use this if you want the output file from one B<parser>
or B<script> step to be used as the input file for another B<parser> or B<script>
step. For example, if the input file for 'stepA' should be created by 'stepB' 
do this:

  <parser>
    <name>stepA</name>
    <inputfrom>stepB</inputfrom>
    <executable>mixmaster</executable>
    <home>/usr/local/bin</home>
    <command>INPUTFILE</command>
    <inputformat></inputformat>
  </parser>

If you do not specify B<inputfrom> for any step then it is assumed that the 
input file comes from the command-line. For example, if you run I<diya> like this:

   diya-script.pl --conf diya.conf seq1.fa

then the input file will be 'seq1.fa' when there is no B<inputfrom> for a given
step.

=back

=head2 The B<script> sections

A B<script> step simply executes and its output is not parsed. For example, you may 
need to copy a sequence file from some location before running a pipeline. 
Or you may want to send the pipeline output somewhere or do an email alert 
after the pipeline is done, you would write B<script> steps for these purposes. 
An example:

 <script>
   <name>formatdb</name>
   <executable>formatdb.sh</executable>
   <command>INPUTFILE</command>
   <home>/Users/bosborne/diya/branches/0.4.0/examples</home>
   <inputfrom>extractCDS</inputfrom>
 </script>

=over 5

=item B<executable>

The name of the script. This should be the actual name, not a
synonym. Required.

=item B<name>

The B<name> for this step in the pipeline. Required.

=item B<home>

The directory where the script is found. Required.

=item B<command>

The command that has to be run, without the executable name. Optional.

=item B<inputfrom>

This is optional. Use this if you want the output file from a B<parser>
or B<script> step to be used as the input file for a B<script> step.

=back

=head1 INITIALIZATION OPTIONS

You can pass options to your pipeline object when you create it with
the new() method.

=over 4

=item B<verbose>

When you set B<verbose> to 1 the pipeline object will print out useful
diagnostic messages. Set verbose to I<true> like this:

  my $pipeline = diya->new(-verbose => 1)

Setting B<verbose> is optional, the default value is 0 or I<false>

=item B<mode>

The values are I<serial> and I<sge>. Set B<mode> like this:

  my $pipeline = diya->new(-mode => 'sge')

Setting B<mode> is optional, the default value is I<serial>.

=item B<use_conf>

Specify the configuration file for the pipeline. The conf file can have
any sort of name as long as it has the correct format. An example:

  my $pipeline = diya->new( -use_conf => "~/myconf.conf" )

Setting B<use_conf> is optional. If it is not set then I<diya> will
look for a file named C<diya.conf> in your $DIYAHOME directory or in the
current working directory.

=item B<outputdir>

Specify the output directory for the pipeline. An example:

  my $pipeline = diya->new( -outputdir => "~/myfiles" )

Setting B<outputdir> is optional. If it is not set then I<diya> will create
an output directory in your $DIYAHOME directory using a timestamp, 
for example "2012_08_07_15_17_34_diya".

=back

=head1 USING A DIYA SCRIPT

You will run I<diya> using a fairly simple script since most of the details 
are in the configuration file. The C<diya.pl> script that comes with 
I<diya> is an example. 

I<diya> scripts are run like this:

 	% diya.pl [options] [input files]

=head2 Command-line options

=over 4

=item B<--verbose>

Set the verbosity level, 0 or 1.

 diya.pl --verbose 1

=item B<--mode> [serial|sge]

Run the batch in serial mode, or sge mode if SGE is available.

 diya.pl --mode serial

=item B<--conf> 

Use the given conf file. If you use this option then this given
conf file will be used, if there is a conf file specified in the new() 
method it will be ignored.

 diya.pl --conf new-diya.conf

=item B<--outputdir>

Set the output directory.

 diya.pl --outputdir /tmp/mypipeline

=back	

=head2 Using --set to modify a command

You can also modify your commands dynamically from the command-line.
For example, you might want to run C<blastall> and create an output file with
a specific name. Here is an example C<blastall> command from a *conf file:

  <command>-p blastp -d ran.fa -i INPUTFILE -o MYOUTPUTFILE</command>

You could run C<diya.pl> like this:

  diya.pl --set MYOUTPUTFILE=blastp.out

And an output file called C<blastp.out> would be created.

You can add these "wild card" words anywhere you want to in the B<command>
line of the *conf file. The only rule is that you should not use the words 
I<INPUTFILE>, I<OUTPUTFILE>, and I<OUTPUTDIR>. These are already being used by 
I<diya>. One way to make sure your "wild card" is unique is to prefix it with 
'MY'. We also suggest capitalizing these words, for clarity.

=head2 Using --set to set a variable in a Perl module

Suppose you want your Perl module to be able to get some value from
the command-line and use it as a variable, e.g. C<$MYDATABASE>. First add the 
variable name to the C<@EXPORT_OK> array in C<diya.pm>. Then modify the C<use diya;> 
line in your Perl module, for example:

 use diya qw($MYDATABASE);

After these modifications you should be able to do the following:

 diya.pl --set MYDATABASE=ncbi Seq.fa

And the variable C<$MYDATABASE> will have the value I<ncbi> in your Perl module
when C<diya.pl> runs.

When you use --set you are creating global variables that can be used
in your own Perl modules so make sure that your variable names do not collide with 
I<diya> variables. One way to do this is to use variable names that are all 
capitalized, or prefix the name with 'MY'.

=head1 WRITING YOUR OWN conf* FILES

L<THE CONFIGURATION FILE> section discusses the structure of the *conf file
but in order to create your own files you will need to understand some
of the internal details of I<diya>. 

When I<diya> runs it can create the names of input and output files.
This makes it easy for I<diya> to keep track of files since one of its jobs is
to pass the output of one step to the next step as input. I<diya> uses a
timestamp and the B<name> of the step to create file names, for example:

  2008_08_07_10_19_53-create-fasta-db.out

=head2 The meaning of I<INPUTFILE>

The file above was created by the I<create-fasta-db> step, as you can see
from its name. This file could be the input for some other step, and you would 
indicate this by using the B<inputfrom> field. For example:

  <parser>
    <inputfrom>create-fasta-db</inputfrom>
    <executable>blastall</executable>
    <home>/usr/local/bin</home>
    <command>-p blastp -i INPUTFILE -d MYDATABASE -o OUTPUTFILE </command>
    <name>blastpCDS</name>
    <inputformat></inputformat>
  </parser>

The block above says that the I<INPUTFILE> should come from the I<create-fasta-db>
step. When I<diya> runs and the actual command is constructed this part of 
the B<command> line:

    -p blastp -i INPUTFILE

Will be transformed into something like:

    -p blastp -i /tmp/2008_08_07_10_19_53-create-fasta-db.out

I<INPUTFILE> has a second meaning, which is the name of the sequence file passed
to the I<diya> script. Recall that you can run a I<diya> script like this:

  mydiya.pl --conf my.conf NC_123456.fa

If a given step has no B<inputfrom> value then the value of I<INPUTFILE> will be
the name of the sequence file set from the command-line, or "NC_123456.fa" in the
example above.

This does not mean that you have to use I<INPUTFILE> in each step. It means that
when I<INPUTFILE> is present in a B<command> line it will substituted in one
of these 2 ways, depending on whether or not there is an B<inputfrom> value.

=head2 The meaning of I<OUTPUTFILE>

The I<OUTPUTFILE> from one step will frequently be used as the I<INPUTFILE> to
another step. Thus you may need to explicitly create a file with the name
contained in I<OUTPUTFILE> using the B<command> line. An example B<script> block:

 <script>
   <name>create-fasta-db</name>
   <executable>createdb.pl</executable>
   <command>-o OUTPUTFILE</command>
   <home>~/scripts</home>
   <inputfrom></inputfrom>
 </script>

When this step runs a command like this will be run and executed:

  ~/scripts/createdb.pl -o /tmp/2008_08_07_10_19_53-create-fasta-db.out

Here C<~/scripts/createdb.pl> will use the file name provided by I<diya>
and put its output into that file. An alternative is to redirect the 
output of an application into an output file. For example:

 <script>
   <name>create-fasta-db</name>
   <executable>createdb.pl</executable>
   <command> > OUTPUTFILE</command>
   <home>~/scripts</home>
   <inputfrom></inputfrom>
 </script>

=head2 The meaning of I<OUTPUTDIR>

By default I<diya> creates an output directory to contain all the input and output files
created during a pipeline run, something like:

   2008_11_21_22_02_19_diya/

You can get this name in order to use it in your B<script> or B<parser>
steps, something like:

 <script>
   <name>move-file</name>
   <executable>move-file.pl</executable>
   <command>-o OUTPUTDIR</command>
   <home>~/scripts</home>
   <inputfrom></inputfrom>
 </script>

In this example the name of the I<diya> output directory will be passed to 
the script.

=head1 CREATING AND LOADING GFF

The authors use I<diya> to annotate sequences and save these annotations
as GenBank files. They routinely convert these files to GFF, load
the GFF into L<Bio::DB::GFF|Bio::DB::GFF> databases, and visualize the annotations
using GBrowse. The conversion script used is C<diya-genbank2gff3.pl>
in the C<scripts> directory, this script is a modification of the
C<genbank2gff3.pl> script that comes with BioPerl.

=head1 AUTHORS

Andrew Stewart, andrew.stewart@med.navy.mil
Brian Osborne, briano@bioteam.net

=head1 CONTRIBUTORS

Tim Read, timothy.read@med.navy.mil

=cut

package diya;

use strict qw(vars subs);
use vars qw($VERSION @EXPORT_OK);
use base 'Exporter';
use XML::Simple;
use Data::Dumper;
use Cwd qw( cwd );
use Data::Merger qw( merger );
use File::Basename qw( basename fileparse );
use POSIX;
# From BioPerl
use Bio::Tools::GuessSeqFormat;
use Bio::SeqIO;

$VERSION = '1.0';

$Data::Dumper::Purity = 1;

# Add a variable here if you want to set it using --set
@EXPORT_OK = qw($VERSION $ID $SEQFILE $MYSPECIES $MYSTRAIN $PROJECT $CONF $OUTPUT $MYSEQID
					 $MYCLUSTERS $MYCDD $MYLEN );

=head1 METHODS

Public methods are listed first, followed by private methods prefixed with '_'.

=cut

=head2 new

 Name    : new
 Usage   : $diya = diya->new()
 Function: create a diya pipeline object
 Returns : a diya object
 Args    : -verbose   (optional), 1 or 0, set the verbosity level
           -use_conf  (optional), name of the conf file to be used
           -outputdir (optional), the directory where results will reside
           -mode      (optional), 'serial' is default 
 Example : my $pipeline = diya->new(-verbose   => 1,
                                    -use_conf  => "latest.conf",
                                    -outputdir => 'mydir',
                                    -mode      => 'sge' );

=cut

sub new {
	my ($module,@args) = @_;
	my $class = ref($module) || $module;
	my $self = {};
	bless($self,$class);

	# Defaults
	$self->verbose(0);
	$self->mode('serial');

	$self->_initialize(@args);
	$self->_greeting if $self->verbose;
	$self->_diyahome;
	$self->_get_options;

	return $self;
}


=head2 read_conf

 Name    : read_conf
 Usage   : $diya->read_conf("my_conf_file")
 Function: read a diya conf file
 Returns : 1 on success
 Args    : the name of the conf file to be read (optional)
 Example : $pipeline->read_conf()  or  $pipeline->read_conf("latest.conf")

=cut

sub read_conf {
	my ($self, $confpath) = @_;
	my $conf;

	if ( $confpath ) {
		die "Conf file $confpath not found\n" if (! -e $confpath);
		$self->_use_conf($confpath);
	} else {
		$confpath = $self->_use_conf;
	}

	print "Using \'$confpath\'\n" if $self->verbose;

	if ( -e $confpath ) {
		print "Reading \'$confpath\'\n" if $self->verbose;
		$conf = XMLin($confpath, ForceArray => 1 );
	} else {
		die "Could not find conf file $confpath";
	}

	# if someone has already set some values using a method or
	# by using a command-line option then do not overwrite them
	if ( defined $self->_conf ) {
		my $oldconf = $self->_conf;
		merger($conf, $oldconf);
	}

	$self->_conf($conf);
	
	# print Dumper $conf if $self->verbose;

	1;
}

=head2 run

 Name    : run
 Usage   : $diya->run()
 Function: run a diya pipeline
 Returns : 1 on success
 Args    : 
 Example : $pipeline->run

=cut

sub run {
	my $self = shift;
	my $conf = $self->_conf;
	my $filecount;

	die "No configuration found - make sure to call read_conf() before run()" 
	  if ( ! defined $conf );
	
	# The file count is set to the number of input sequence files, if input 
	# files are specified at the command-line. If no input files are specified 
	# it is assumed that the pipeline is just run once and the file count is 1

	if ( @ARGV ) {
		$self->inputfile(@ARGV);
		$filecount = $#ARGV + 1;
	} else {
		$filecount = 1;
	}

	while ( $filecount ) {
		print "File count is $filecount\n" if $self->verbose;

		my $file = $self->_next_inputfile;
		$self->{_lastsgeid} = undef;
		$self->_check_outputdir;

		my $stepcount = 0;

		for my $step ( $self->order ) {
			$stepcount++;
			print "Preparing to run \'$step\', step count of $stepcount\n" 
			  if $self->verbose;

			my $type = $self->_get_type($step);			

			my $modulename = $self->_load_app_module($step) 
			  if ( $type eq 'parser' );

			$self->_check_input_sequence($step,$file) if $file;
			$self->_check_inputfile($step);
			$self->_check_executable($step);
			$self->_make_outputfilename($step);

			my $fullcommand = $self->_make_command($step);
		
			$self->_execute($fullcommand, $stepcount, $modulename, $type, $file);
		}

		$filecount--;
		# Delay 1 second to assure that timestamps on output files differ
		sleep(1);
	}
	1;
}

=head2 new_parser

 Name    : new_parser
 Usage   : $parser = $module->new_parser
 Function: instantiate a new parser object 
 Returns : a new parser object
 Args    : none
 Example : 

=cut

sub new_parser {
	my ($caller,@args) = @_;
	my $class = ref($caller) || $caller;
	my $self = {};

	bless($self,$class);
	return $self;
}

=head2 order

 Name    : order
 Usage   : $diya->order( @array ) or $order = $diya->order
 Function: get or set the order of the steps to be run
 Returns : array of step names
 Args    : To set pass an array of one or more step names the parsers and
           scripts must exist in the conf file in the <parser> and
           <script> sections
 Example : $pipeline->order( qw(tRNAscanSE blastall) )  or 
           $pipeline->order('tRNAscanSE') or
           @my_order = $self->order

=cut

sub order {
	my ($self,@new_order) = @_;
	my $conf = $self->_conf;
	my @current_order;

	if ( ! @new_order ) {
		@current_order = split /\s+/, $conf->{order}->[0]->{names}->[0];
		
		# if <order> is empty then XML::Simple puts a hash there
		if ( $current_order[0] =~ /HASH\(\S+\)/ ) {
			print "No names found in <order>\n" if $self->verbose;
			return;
		}

		print "Current <order> is \'@current_order\'\n\n" if 
		  $self->verbose;
		return @current_order;
	}

	my @conf_scripts = $self->_scripts;
	my @conf_modules = $self->_parsers;
	die "No parsers found - you may need to run read_conf()" if ( ! @conf_modules );
	
	LOOP: foreach my $exe (@new_order) {
		next LOOP if ( grep /^$exe$/, @conf_modules );

		next LOOP if ( grep /^$exe$/, @conf_scripts );

		die "Cannot use new order: the name \'$exe\' can not be found in " .
		  $self->{use_conf};		  
	}

	my $new_order = join ' ', @new_order;
	$conf->{order}->[0]->{names}->[0] = $new_order;
	print "New <order> is \'" . $conf->{order}->[0]->{names}->[0] . "\'\n" 
	  if $self->verbose;

	return @new_order;
}

=head2 write_conf

 Name    : write_conf
 Usage   : $diya->write_conf("my_new_conf_file")
 Function: write a conf file - if no name is supplied then the new file will
           be given a name of format <timestamp>-diya.conf
           (e.g.  2008-06-29-11:35:38-diya.conf )
 Returns : the name of the conf file that was written
 Args    : the name of the conf file that will be written (optional)
 Example : $pipeline->write_conf()  or $pipeline->write_conf("version2.conf")

=cut

sub write_conf {
	my ($self,$confname) = @_;
	my $conf = $self->_conf;

	if ( ! defined $confname ) {
		my $timestamp = strftime("%Y_%m_%d_%H_%M_%S", localtime);
		$confname = $timestamp . "-diya.conf";
		print "No name given for new diya.conf file - " .
		  "will use \'$confname\'\n" if $self->verbose;
	}

	my $header = '<?xml version="1.0" encoding="UTF-8"?>
<!-- $Id: diya.pm 340 2009-04-24 15:03:35Z briano $ -->

';
	open MYOUT,">$confname" or die "Cannot create conf file named $confname";
	print MYOUT $header;
	close MYOUT;

	open MYOUT,">>$confname" or die "Cannot add to conf file named $confname";
	print MYOUT XMLout($conf, noattr => 1, RootName => 'conf');
	close MYOUT;
	$confname;
}

=head2 verbose

 Name    : verbose
 Usage   : $diya->verbose($num) or $verbose_level = $diya->verbose
 Function: get or set the verbose level
 Returns : the verbose level
 Args    : 
 Example : $pipeline->verbose(1)

=cut

sub verbose {
	my ($self,$verbose) = @_;

	if ( ! defined $verbose ) {
		return $self->{verbose};
	} else {
		die "Set verbose() to 1 or 0, not $verbose" if ( $verbose !~ /^[10]$/ );
		$self->{verbose} = $verbose;		
	}

	return $self->{verbose};
}

=head2 project

 Name    : project
 Usage   : $diya->project($num) or $project = $diya->project
 Function: get or set the NCBI project number
 Returns : the NCBI project number
 Args    : 
 Example : $pipeline->project(1355)

=cut

sub project {
	my ($self,$project) = @_;

	if ( ! defined $project ) {
		return $self->{verbose};
	} else {
		die "Set project() to a number, not \'$project\'" if ( $project !~ /^\d+$/ );
		$self->{project} = $project;		
	}

	return $self->{project};
}

=head2 outputdir

 Name    : outputdir
 Usage   : $diya->outputdir()
 Function: get or set the name of the output directory, where all of 
           the files created by the pipeline will be written - the 
           object will try and create the directory if it does not exist

           if no output directory is specified then an output directory will
           be created based on a timestamp, e.g. "2008-06-29-11:35:38-diya"

 Returns : name of the output directory or 0 if not output directory is set
 Args    : 
 Example : $diya->outputdir("pipe-output")

=cut

sub outputdir {
	my ($self,$dir) = @_;

	$self->{outputdir} = $dir if defined $dir;

	$self->{outputdir} .= "/" if ( $self->{outputdir} !~ m|/$| );

	if (! defined $self->{outputdir}) {
		print "No output directory defined - will create one when run() is called\n"
		  if $self->verbose;
		return 0;
	}

	print "Using \'" . $self->{outputdir} . "\' as output directory\n" 
	  if $self->verbose;
	$self->{_outputdir} = $self->{outputdir};
	return $self->{outputdir};
}

=head2 mode

 Name    : mode
 Usage   : 
 Function: get or set the mode corresponding to a pipeline
 Returns : "serial" or "sge"
 Args    : 
 Example : $pipeline->mode("serial")

=cut

sub mode {
	my ($self,$mode) = @_;
	my @modes = qw(sge serial lsf);

	if ( defined $mode ) {
		die "mode() was called with \'$mode\' but the only available modes are: @modes"
		  if ( ! grep /^$mode$/i, @modes );
		$self->{conf}->{run}->{mode} = $mode;
	}

	print "<mode> for the pipeline is \'" . $self->{conf}->{run}->{mode} . "\'\n"
	  if $self->verbose;
	return $self->{conf}->{run}->{mode};
}

=head2 cleanup

 Name    : cleanup
 Usage   : $diya->cleanup
 Function: remove extraneous files created when a pipeline is run
 Returns : 1 on success
 Args    : none
 Example :

=cut

sub cleanup {
	my $self = shift;
	my $dir = $self->outputdir();
	system "rm -fr $dir";
	print "Removed directory \'$dir\'\n" if $self->verbose;
	1;
}

=head2 inputfile

 Name    : inputfile
 Usage   : $diya->inputfile('NC.gbk')
 Function: Get or set the names of the input sequence files
 Returns : 
 Args    : 
 Example : $self->inputfile("234.fa") or $self->inputfile( qw(234.fa AB.fa) )

=cut

sub inputfile {
	my ($self, @files) = @_;
	push @{$self->{inputfiles}}, @files if @files;
	my @arr = @{$self->{inputfiles}};
	
	print "Using \'@arr\' as input files\n" if $self->verbose.
	return @arr;
}

=head2 _next_inputfile

 Name    : _next_inputfile
 Usage   : 
 Function: Get the name of the next input sequence file, remove the 
           last from the queue
 Returns : 
 Args    : 
 Example : 

=cut

sub _next_inputfile {
	my $self = shift;

	if ( ! defined @{$self->{inputfiles}} 
		  || scalar @{$self->{inputfiles}} == 0 ) {
		print "No input files\n" if $self->verbose; 
		return 0;
	}

	$self->{inputfile} = shift @{$self->{inputfiles}};
	print "Next input file is \'" . $self->{inputfile} . "\'\n" if $self->verbose; 
	return $self->{inputfile};
}

=head2 _execute

 Name    : _execute
 Usage   : $self->_execute($command)
 Function: encapsulate the serial and sge execution logic 
 Returns : none
 Args    : command
 Example : 

=cut

sub _execute {
	my ($self, $fullcommand, $count, $modulename, $type, $file) = @_;

	my $output = `$fullcommand`;

	# serial mode execution
	if ( $self->mode ne 'sge' && $self->mode ne 'lsf' ) {
		if ( $? == 0 ) {
			# if $CHILD_ERROR == 0
			print $output, "\n";
			if ( $modulename ) {
				my $parser = $modulename->new_parser;
				print "Created \'$modulename\' object\n" if $self->verbose;
				$parser->parse($self);
			}
			return;
		} else {
			print "Child error: $?\n";
			if ( $output ) {
				die $output, "\n";
			} else {
				exit 1;
			} 
		}
	}
	if ($self->mode eq 'lsf' && $output=~/<(\d+)>/)
	{
		# Extract job id from output
		$output=$1;
	}
	# Handle SGE logic
	if ( $? == 0 ) {
		if ( $output =~ /^(\d+)$/ )	{
			print "sge job id is $1\n" if $self->verbose();
			$self->_lastsgeid($1);
		} else {
			print $output, "\n";
			die "Failed in extracting job id from sge submission output\n";
		}
	} else {
		print $output, "\n";
		die "sge submission fails\n";	
	}
	
	
	if (!defined($modulename))
	{
		return;
	}
	
	#handle the parse logic only if $modulename is defined.
	my $scriptfile = $self->{outputdir} . "/interscript_$count.pl";
	open(OUT, ">$scriptfile") or die "can not open $scriptfile for writing";
	# if the sequence is dumped, the object inheritance is broken when it is
	# revaluated in the generated script below. 
	# to solve this issue, we do not dump the sequence object. Instead
	# the sequence object is rebuilt in the generated script by calling _reconstruct_sequence 
	my $sequence = $self->_sequence();
	$self->{sequence} = undef;
	my $diyastr = Data::Dumper->Dump( [ $self], ["diya"] );
		
	print OUT qq(#!/usr/bin/perl
use strict;
use diya;
use $modulename;
my \$parser={};
bless(\$parser, "$modulename");
my \$diya;
$diyastr;
\$diya->_reconstruct_sequence();
\$parser->parse(\$diya,'$file');
);
	close(OUT);
	$self->_sequence($sequence);

	`chmod uog+x $scriptfile`;

	my $cwd = getcwd();
	my $sgeid = $self->_lastsgeid();
	my $outputdir = $self->outputdir();
	my $intercommand = "qsub -b y -V -cwd -e $outputdir/intercmd_${count}.err -o $outputdir/intercmd_${count}.out -N intercmd -terse -hold_jid $sgeid $scriptfile";
	if ($self->mode eq 'lsf')
	{
		$intercommand = "bsub -e $outputdir/intercmd_${count}.err -o $outputdir/intercmd_${count}.out -J intercmd -w $sgeid $scriptfile";
				
	}
	print "command is \'$intercommand\'\n" if $self->verbose();
	$output = `$intercommand`;
	if ($self->mode eq 'lsf' && $output=~/<(\d+)>/)
	{
		#extract jobid from lsf.
		$output=$1;
	}
	if ( $? == 0 ) {
		if ( $output =~ /^(\d+)$/ ) {
			print "sge job id is $1\n" if $self->verbose();
			$self->_lastsgeid($1);
		} else {
			print $output, "\n";
			die "failed in extracting job id from sge submission output\n";
		}
	} else {
		print $output, "\n";
		die "sge submission fails\n";	
	}	
}

=head2 _reconstruct_sequence

 Usage   : _reconstruct_sequence
 Function: reconstruct the sequence object. Used only when mode is sge. 
           When mode is sge, the _execute() will generate an intermediate
           script performing the $parse->parse($diya). In the intermediate script,
           this method is called. Please see _execute 
 Returns : none
 Args    : none
 Example : $self->_reconstruct_sequence()

=cut

sub _reconstruct_sequence {
	my $self = shift;
	my $in = Bio::SeqIO->new(-format => $self->{__sequenceformat}, 
									 -file   => $self->{__sequencefile});
	my $seq = $in->next_seq;
	$self->_sequence($seq);
}

=head2 _check_executable

 Usage   : _check_executable
 Function: checks to see that the executable exists
 Returns : 1 or die
 Args    : Name of step
 Example : $self->_check_executable($step)

=cut

sub _check_executable {
	my ($self,$step) = @_;
	my $exe = $self->_home($step) . $self->_executable($step);
	die "Executable \'$exe\' not found" if ( ! -e $exe ); 
	1;
}

=head2 _check_input_sequence

 Usage   : _check_input_sequence
 Function: checks to see that the format of the input sequence file is correct,
           if the format is not correct then it creates a sequence file
           of the correct format
 Returns : The name of the file of the correct format
 Args    : Name of step
 Example : $self->_check_input_sequence($step)

=cut

sub _check_input_sequence {
	my ($self,$module,$file) = @_;
	
	my $requiredformat = $self->_inputformat($module);

	my $guesser = Bio::Tools::GuessSeqFormat->new( -file => $file );
	my $inputformat = $guesser->guess;
	print "Format of \'$file\' is \'$inputformat\'\n" if $self->verbose;

	my $in = Bio::SeqIO->new(-format => $inputformat, -file => $file);
	my $seq = $in->next_seq;
	$self->_sequence($seq);

	if (! $requiredformat) {
		print "No <inputformat> found for step \'$module\', no conversion needed\n"
		  if $self->verbose;
		return $file;
	}
	
	print "Step is \'$module\', input file format looks like \'$inputformat\',". 
	  " required format is \'$requiredformat\'\n" if $self->verbose;

	# keep these two values so that _reconstruct_sequence can work.
	$self->{__sequenceformat} = $inputformat;
	$self->{__sequencefile} = $file;

	if ( $requiredformat !~ /$inputformat/i ) {
		
		print "Format of \'$file\' does not match required format, will convert file" .
		  " to \'$requiredformat\'\n" if $self->verbose;

		my $newfile = basename($file);

		if ( $newfile =~ /(\S+)\.\S+$/ ) {
			$newfile = $1 . "." . $requiredformat;
		} else {
			$newfile = $file . "." . $requiredformat;
		}

		$newfile = $self->{outputdir} . $newfile;
		print "New file name is \'$newfile\'\n" if $self->verbose;

		my $out = Bio::SeqIO->new(-format => $requiredformat, -file => ">$newfile" );
		$out->write_seq($seq);

		$self->{inputfile} = $newfile;

		return $newfile;
	}

	$file;
}

=head2 _check_inputfile

 Usage   : _check_inputfile
 Function: records the input file name for the step - this may be
           the input sequence for the entire pipeline but a step may also
           use the output of another step as input
 Returns : the name of the input file for the step
 Args    : name of step
 Example : $self->_check_inputfile($step)

=cut

sub _check_inputfile {
	my ($self,$module) = @_;
	my $conf = $self->_conf;
	my $inputmodule = $self->_inputfrom($module);

	if ( $inputmodule ) {
		if ( $self->{$inputmodule}->{outputfile} ) {
			$self->{$module}->{inputfile} = $self->{$inputmodule}->{outputfile};
			print "Output file from \'$inputmodule\' will be used as input file for \'$module\'\n"
			  if $self->verbose;
			return $self->{$inputmodule}->{outputfile};
		} else {
			die "No output file from \'$inputmodule\' found - " .
			  "make sure that \'$inputmodule\' runs before \'$module\'";
		}
	} else {

		if ( defined $self->{inputfile} ) {
		
			$self->{$module}->{inputfile} = $self->{inputfile};
			print "\'" . $self->{$module}->{inputfile} . "\' will be used as input file for \'$module\'\n"
			  if $self->verbose;

			# Copy input file to output directory if it is not present there
			my $outputdir = $self->outputdir;
			my $inputfile = $self->{$module}->{inputfile};
			my $filename = fileparse($inputfile);

			if ( ! -e "$outputdir$filename" ) {
				print "Copying \'$filename\' to \'$outputdir\'\n" if $self->verbose;
				system "cp $inputfile $outputdir" ;
			}

			$self->{$module}->{inputfile} = "$outputdir$filename";
			print "\'" . $self->{$module}->{inputfile} . "\' will be used as input file for \'$module\'\n"
			  if $self->verbose;
			
			return $self->{$module}->{inputfile};

		}
	}	
}

=head2 _check_outputdir

 Name    : _check_outputdir
 Usage   : $diya->_check_outputdir
 Function: checks that there is an output directory - if no output directory is defined
           then a directory name will be made using a timestamp 
           (e.g.  "2008-06-29-11:35:38-diya") and this directory will be in the diya home
           directory - if an output directory is defined but does not exist then we will 
           attempt to create it
 Returns : name of output directory, on success
 Args    : none
 Example :

=cut

sub _check_outputdir {
	my $self = shift;

	if ( ! defined $self->{outputdir} ) {
		my $timestamp = strftime("%Y_%m_%d_%H_%M_%S", localtime);
		$self->{outputdir} = $self->_diyahome . $timestamp . "_diya";
		$self->{_outputdir} = $self->{outputdir};
		print "No output directory defined, will use \'" . $self->{outputdir} .
		  "\'\n" if $self->verbose;
	}

	if ( $self->mode eq 'sge' ) {
		my $filenum = scalar($self->inputfile); 
		if ( $filenum > 0 ) {
			$self->{outputdir} = $self->{_outputdir} . "_$filenum";
		} else {
			$self->{outputdir} = $self->{_outputdir};
		}
	}

	if (! -d $self->{outputdir}) {
		print "Directory \'" . $self->{outputdir} . "\' does not exist, will create it\n" 
		  if $self->verbose;
		my $dir = $self->{outputdir};
		system("mkdir $dir");
	}
	
	$self->{outputdir} .= "/" if ( $self->{outputdir} !~ m|/$| );

	return $self->{outputdir};
}

=head2 _make_command
 
 Name    : _make_command
 Usage   : $pipeline->_make_command($step)
 Function: to create a complete command using information from the conf file 
           and any command-line options - private method called by run()
 Returns : a command, ready to execute
 Args    : the name of the parser step (e.g. "tRNAscanSE")
 Example : $pipeline->_make_command($parser)

=cut

sub _make_command {
	my ($self,$module) = @_;
	my @substitutions = qw( INPUTFILE OUTPUTFILE );
	
	my $cmd = $self->_home($module) . $self->_executable($module);
	$cmd .=  " " . $self->_command($module) if ( $self->_command($module) ); 
	print "\'$module\' command before substitution is \'$cmd\'\n" if $self->verbose;

	# First substitute any values specific to steps
	for my $word (@substitutions) {
		my $lcword = lc $word;
		my $val = $self->{$module}->{$lcword};
		if ( $val ) {
			print "\'$word\' will be substituted with \'$val\'\n" if $self->verbose;
			# Arguments need to be single-quoted if they contain a space
			if ( $val =~ m|\s| ) {
				$cmd =~ s/$word/'$val'/g;
			} else {
				$cmd =~ s/$word/$val/g;
			}
		}
	}
	print "\'$module\' command after parser-specific substitution is \'$cmd\'\n"
	  if $self->verbose;

	# Second, substitute any values not specific to a specific step, like OUTPUTDIR
	my @global_subs = qw( OUTPUTDIR );
	for my $word ( @global_subs ) {
		my $val = $self->outputdir if ( $word eq 'OUTPUTDIR' ); 
		print "\'$word\' will be substituted with \'$val\'\n" if $self->verbose;
		$cmd =~ s/$word/$val/g;
	}
	print "\'$module\' command after global substitution is \'$cmd\'\n" if $self->verbose;

	# Last, substitute the values of any globals created by using --set NAME=VAL
	for my $word ( keys %{$self->{substitutions}} ) {
		# Arguments need to be single-quoted if they contain a space
		if ( $self->{substitutions}->{$word} =~ m|\s| ) {
			$cmd =~ s/$word/'$self->{substitutions}->{$word}'/g;
		} else {
			$cmd =~ s/$word/$self->{substitutions}->{$word}/g;
		}
		print "\'$word\' will be substituted with \'" . $self->{substitutions}->{$word} . "\'\n"
		  if $self->verbose;
	}
	print "\'$module\' command after variable substitution is \'$cmd\'\n" if $self->verbose;

	my $lastsgeid = $self->_lastsgeid;

	if ( $self->mode eq 'sge' ) {
			my $joberr = $self->outputdir()."/".$self->_executable($module) . ".err";
			my $jobout = $self->outputdir()."/".$self->_executable($module) . ".out";
			my $name   = $self->_executable($module);
			my $holdid = "";
			if ( defined($lastsgeid) )	{
				$holdid = "-hold_jid $lastsgeid";
			}
			$cmd = "qsub -b y -V -cwd -e $joberr -o $jobout -N $name -terse $holdid $cmd";
	} elsif ($self->mode eq 'lsf')
	{
			my $joberr = $self->outputdir()."/".$self->_executable($module) . ".err";
			my $jobout = $self->outputdir()."/".$self->_executable($module) . ".out";
			my $name   = $self->_executable($module);
			my $holdid = "";
			if ( defined($lastsgeid) )	{
				$holdid = "-w $lastsgeid";
			}
			$cmd = "bsub -e $joberr -o $jobout -J $name $holdid $cmd";
	}
	
	print "Command is \'$cmd\'\n" if ( $self->verbose && $self->mode eq 'sge' );
	return $cmd;
}

=head2 _make_outputfilename

 Name    : _make_outputfilename
 Usage   : $diya->_make_outputfilename($parser)
 Function: create an output file name with parser step name and timestamp,
           private method, called by run()
 Returns : output file name
 Args    : step name
 Example :

=cut

sub _make_outputfilename {
	my ($self,$module) = @_;

	my $timestamp = strftime("%Y_%m_%d_%H_%M_%S", localtime);
	my $filename = "$timestamp-$module.out";
	$filename = $self->outputdir . $filename;
	$self->{$module}->{outputfile} = $filename;
	print "Output file name for \'$module\' is \'" . $self->{$module}->{outputfile}  .
	  "\'\n" if $self->verbose;

	return $self->{$module}->{outputfile} ;
}

=head2 _lastsgeid

 Name    : _lastsgeid
 Usage   : $diya->_lastsgeid()
 Function: set or get the last sge job id submitted by current process,
 			  used only internally for job id tracking.
 Returns : last sge job id submitted by current process.
 Args    : 
 Example : $pipeline->_lastsgeid(53)

=cut

sub _lastsgeid {
	my ($self,$id) = @_;

	if ( $id ) {
		return $self->{_lastsgeid} = $id;
	} else {
		 $self->{_lastsgeid};
	}
}

=head2 _outputfile

 Name    : _outputfile
 Usage   : $file = $self->_outputfile($parser)
 Function: get the output file name for a given parser step, 
           private method called by run() or a parser module
 Returns : output file name
 Args    : step name
 Example : $file = $self->_outputfile($step)

=cut

sub _outputfile {
	my ($self,$module) = @_;

	if ( $self->{$module}->{outputfile} ) {
		print "Found output file \'" .  $self->{$module}->{outputfile} .
		  "\' for \'$module\'\n" if $self->verbose;
		return $self->{$module}->{outputfile};
	} else {
		die "No output file name for parser \'$module\' found";
	}
}

=head2 _greeting

 Name    : _greeting
 Usage   : $diya->_greeting
 Function: print a greeting - private method, called by new()
 Returns : nothing
 Args    : 
 Example :

=cut

sub _greeting {
	my $self = shift;
	my $greeting = "
Do It Yourself Annotator (diya) v. $VERSION;
\t--help for usage
\t--license for copyright information

";
	print $greeting;
}


=head2 _help

 Name    : _help
 Usage   : $diya->_help
 Function: print the POD - works only if DIYA is installed
 Returns : nothing
 Args    : 
 Example :

=cut

sub _help {
	my $self = shift;
	eval {
		system "perldoc diya";
	};
	print "Error using _help(), you must install diya first\n" if $@;
}

=head2 _diyahome

 Name    : _diyahome
 Usage   : 
 Function: add the path to the diya home directory to the object -  the path to the 
           diya package comes from the env $DIYAHOME. If this not set then try to 
           use the current working directory. Private method, called by new()
 Returns : the diya home directory
 Args    : none
 Example : 

=cut

sub _diyahome {
	my $self = shift;

	if ( ! defined $ENV{DIYAHOME} or ! $ENV{DIYAHOME} ) {
		$self->{diyahome} = cwd();
	} else {
		$self->{diyahome} = $ENV{DIYAHOME};
	}

	$self->{diyahome} .= "/" if ( $self->{diyahome} !~ m|/$| );
	print "The diya home directory is \'" . $self->{diyahome} . "\'\n"  if $self->verbose;
	return $self->{diyahome};
}

=head2 _conf

 Name    : _conf
 Usage   : 
 Function: get or set a hash representing the conf file - private method, 
           called by read_conf() and write_conf()
 Returns : a hash reference representing the conf file
 Args    : 
 Example : $conf = $diya->_conf

=cut

sub _conf {
	my ($self,$conf) = @_;
	$self->{conf} = $conf if $conf;
	return $self->{conf};
}

=head2 _executable

  Name    : _executable
  Usage   : 
  Function: return the executable name corresponding to a parser
            or script, private method called by _make_command
            and _check_executable
  Returns : executable name
  Args    : parser or script name 
  Example : $exe = $self->_executable($name)

=cut

sub _executable {
 	my ($self, $name) = @_;
	my $conf = $self->_conf;

	my @types = qw(parser script);

	for my $type ( @types ) {
		for my $key ( @{$conf->{$type}} ) {
			if ( $key->{name}->[0] eq $name ) {
				print "Found <executable> \'" . $key->{executable}->[0] . "\' for \'" .
				  $name . "\' in \'" . $self->_use_conf . "\'\n" if $self->verbose;
				return $key->{executable}->[0];
			}
		}
	}
	die "Could not find <executable> for \'$name\'";
}

=head2 _parsers

 Name    : _parsers
 Usage   : 
 Function: return all the parser names in the conf file,
           private method called by order()
 Returns : array of parser names  
 Args    : none
 Example : @parsers = $self->_parsers

=cut

sub _parsers {
	my $self = shift;
	my $conf = $self->_conf;
	my @modules;

	for my $module ( @{$conf->{parser}} ) {
		$module->{name}->[0] =~ s/\s+//g;
		push @modules, $module->{name}->[0];
	}
	print "Parsers in \'" . $self->{use_conf} . "\': \'@modules\'\n" if 
	  $self->verbose;

	return @modules;
}

=head2 _scripts

 Name    : _scripts
 Usage   : 
 Function: return all the script names in the conf file,
           private method called by order()
 Returns : array of script names  
 Args    : none
 Example : @scripts = $self->_scripts

=cut

sub _scripts {
	my $self = shift;
	my $conf = $self->_conf;
	my @scripts;

	for my $script ( @{$conf->{script}} ) {
		$script->{name}->[0] =~ s/\s+//g;
		push @scripts, $script->{name}->[0];
	}
	print "Scripts in \'" . $self->{use_conf} . "\': \'@scripts\'\n" if 
	  $self->verbose;

	return @scripts;
}

=head2 _sequence

 Name    : _sequence
 Usage   : $diya->_sequence($seq) or $seq = $diya->_sequence
 Function: get or set the DNA sequence object, private method used by parser
           modules and by _check_input_sequence() 
 Returns : the sequence object
 Args    : 
 Example : $pipeline->_sequence($seq)

=cut

sub _sequence {
	my ($self,$seq) = @_;

	if ( ! defined $seq ) {
		return $self->{sequence};
	} else {
		$self->{sequence} = $seq;		
	}

	return $self->{sequence};
}

=head2 _home

 Name    : _home
 Usage   : 
 Function: return the home or location corresponding to an executable, 
           private method called by _make_command and _check_executable
 Returns :  
 Args    : parser or script name 
 Example : $path = $self->_home($exe)

=cut

sub _home {
	my ($self, $name) = @_;
	my $conf = $self->_conf;
	my @types = qw(parser script);

	for my $type ( @types ) {
		for my $key ( @{$conf->{$type}} ) {
			if ( $key->{name}->[0] eq $name ) {
				$key->{home}->[0] .= "/" if ( $key->{home}->[0] !~ m|/$| );
				print "Found <home> \'" . $key->{home}->[0] . "\' for \'" .
				  $name . "\'\n" if $self->verbose;
				return $key->{home}->[0];
			}
		}
	}
	die "Could not find <home> for \'$name\'";
}

=head2 _inputfrom

 Name    : _inputfrom
 Usage   : 
 Function: return the inputfrom field corresponding to a parser or
           script
 Returns :  
 Args    : parser or script name
 Example : $path = $self->_inputfrom($module)

=cut

sub _inputfrom {
	my ($self,$name) = @_;
	my $conf = $self->_conf;
	my @types = qw(parser script);

	for my $type ( @types ) {
		for my $key ( @{$conf->{$type}} ) {
			if ( $key->{name}->[0] eq $name ) {
				# if there's an empty hash there, which is what XML::Simple does
				if ( ref $key->{inputfrom}->[0] eq 'HASH' ) {
					print "No <inputfrom> found for \'" .
					  $name . "\'\n" if $self->verbose;
					return 0;
				} else {
					print "Found <inputfrom> \'" . $key->{inputfrom}->[0] . "\' for \'" .
					  $name . "\'\n" if $self->verbose;
					return $key->{inputfrom}->[0];
				}
			}
		}
	}
	die "No <inputfrom> found for \'$name\'";
}

=head2 _command

 Name    : _command
 Usage   : 
 Function: return the command corresponding to a step,
           private method called by _make_command()
 Returns :  
 Args    : step name
 Example : $args = $self->_command($step)

=cut

sub _command {

	my ($self,$name) = @_;
	my $conf = $self->_conf;
	my @types = qw(parser script);
	my $type = "";

	for $type ( @types ) {
		for my $key ( @{$conf->{$type}} ) {
			if ( $key->{name}->[0] eq $name ) {
				# if there is an empty hash, created by XML::Simple from the *conf
				if ( ref $key->{command}->[0] eq 'HASH' ) {
					print "No <command> for $type \'" . $name . "\'\n" if $self->verbose;
					return;
				} elsif (  $key->{command}->[0] ) {
					print "Found <command> \'" . $key->{command}->[0] . "\' for \'" .
					  $name . "\'\n" if $self->verbose;
					return $key->{command}->[0];
				}
			}
		}
	}
}

=head2 _inputformat

 Name    : _inputformat
 Usage   : 
 Function: return the file format required by a parser step,
           private method called by run()
 Returns : format name, or 0 if no value is found
 Args    : parser step name
 Example : $format = $self->_inputformat($step)

=cut

sub _inputformat {
	my ($self,$name) = @_;
	my $conf = $self->_conf;
	
	return 0 if ( $self->_get_type($name) eq 'script' );

	for my $key ( @{$conf->{parser}} ) {
		if ( $key->{name}->[0] eq $name ) {
			# if XML::Simple has created an empty hash
			if ( ref $key->{inputformat}->[0] eq 'HASH' ) {
				print "No <inputformat> for \'" . $name . "\'\n" if $self->verbose;
				return 0;
			} else {
				print "Found <inputformat> \'" . $key->{inputformat}->[0] . 
				  "\' for \'" . $name . "\'\n" if $self->verbose;
				return $key->{inputformat}->[0];	
			} 
		}
	}
	die "Could not find input format for " . $name;
}

=head2 _use_conf

 Name    : _use_conf
 Usage   : $diya->_use_conf("conf_file")
 Function: add the name of the conf file being used to the object - private method, 
           called by read_conf()
 Returns : the name of the conf file being used
 Args    : 
 Example : 

=cut

sub _use_conf {
	my ($self,$use_conf) = @_;
	
	# If a conf file is already specified, e.g. from the command-line
	my $existing_conf = $self->_conf;

	if ( -e $existing_conf ) {
		print "Conf file \'$existing_conf\' will be used\n" 
		  if $self->verbose;
		return;
	}

	if ( $use_conf ) {
		$self->{use_conf} = $use_conf;
		return $self->{use_conf};
	}

	return $self->{use_conf} if $self->{use_conf};

	# if no conf file name is provided look in $DIYAHOME
	my $diyahome = $self->_diyahome;
	if ( -e $diyahome . "diya.conf" ) {
		$self->{use_conf} = $diyahome . "diya.conf";
		print "Will use \'" . $self->{use_conf} . "\'\n" if $self->verbose;
		return $self->{use_conf};
	}
	
	# if no conf file is found in $DIYAHOME look in current directory
	my $cwd = cwd();
	if ( -e $cwd . "/diya.conf" ) {
		$self->{use_conf} = $cwd . "/diya.conf";
		print "Will use \'" . $self->{use_conf} . "\'\n" if $self->verbose;
		return $self->{use_conf};
	}

	die "Could not find diya.conf in $cwd or $diyahome";
}

=head2 _initialize

 Name    : _initialize
 Usage   : $diya->_initialize
 Function: add parameters to the diya object - strips the dash used by named 
           parameters, private method called by new()
 Returns : 
 Args    : 
 Example : 

=cut

sub _initialize {
	my $self = shift;
	# my @params = qw(verbose conf file);

	while ( @_ ) {
		( my $key = shift ) =~ s/^-//;
		$self->{$key} = shift;
	}
 	if (exists($self->{mode}))
 	{
 		$self->mode($self->{mode});
 	}
}

=head2 _get_type

 Name    : _get_type
 Usage   : $type = $diya->_get_type($step)
 Function: return 'script' or 'parser', private method called by run()
 Returns : 'script' or 'parser'
 Args    : 
 Example : 

=cut

sub _get_type {
	my ($self, $name) = @_;
	my $conf = $self->_conf;
	my @types = qw( parser script );

	for my $type ( @types ) {
 		for my $key (  @{$conf->{$type}} ) {			
			if ( $key->{name}->[0] eq $name){
				print "\'$name\' is a \'$type\'\n" if $self->verbose;
				return $type;
			}
 		}
	}

	die "No type found for \'$name\' in " . $self->_use_conf;
}

=head2 _load_app_module

 Name    : _load_app_module
 Usage   : $diya->_load_app_module($module)
 Function: call require on the given module, private method called by run()
 Returns : full name of module, e.g. "diya::tRNAscanSE"
 Args    : name of the module, e.g. 'tRNAscanSE'
 Example : 

=cut

sub _load_app_module {
	my ($self, $module) = @_;

	$module = "diya::" . $module;

	if ($module !~ /^([\w:]+)$/) {
		print "\'$module\' is an illegal Perl package name\n" if $self->verbose;
		$module =~ s/[^\w:]//g;
		print "Will try to use \'$module\' instead\n" if $self->verbose;
	}

	my $modulepath = $module . ".pm";
	$modulepath =~ s|::|/|g;

	# e.g. require "diya/tRNAscanSE.pm"
	eval {
		require "$modulepath";
	};
	
	die "$@" if $@;
	print "Loaded module \'$module\'\n" if $self->verbose;
	
  return $module;
}

sub _license {
	my $self = shift;
print "
#--------------------------------------------------------------------------
# Copyright 2008
#
# This file is part of DIYA.
#
# DIYA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIYA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the diya package.  If not, see <http://www.gnu.org/licenses/>.
#--------------------------------------------------------------------------
";
}

=head2 _get_options

 Name    : _get_options
 Usage   : 
 Function: get command-line options - the --set option is used to
           create globals that can be imported into a 
           parser module or used directly in this module
 Returns : 1 on success 
 Args    : none 
 Example : At the command-line: diya.pl --set REFD=~/mydb.ref

=cut

sub _get_options {
	my $self = shift;
	
	use Getopt::Long;

	GetOptions("mode=s"    => sub { 
					  my ($name,$val) = @_;
					  $self->mode($val); 
				  },

				  "verbose=i" => sub { 
					  my ($name, $val) = @_;
					  $self->verbose($val);
				  },

				  "project=i" => sub { 
					  my ($name, $val) = @_;
					  $self->project($val); 
				  },

				  "outputdir=s" => sub { 
					  my ($name, $val) = @_;
					  $self->outputdir($val); 
				  },

				  "conf=s"    => sub { 
					  my ($name, $val) = @_;
					  $self->_use_conf($val); 
				  },

				  "license"	  => sub { 
					  $self->_license; 
				  },

				  "help|h"	  => sub { 
					  $self->_help; 
				  },

 				  "set=s%"	  => sub {
 					  my ($op, $name, $value) = @_;
 					  if ( defined $name && defined $value ) {
 						  $$name = $value;
 						  $self->{substitutions}->{$name} = $value;
 						  print "--set: the value of the global variable " . '$' . $name . 
 							 " will be \'" . $value . "\'\n" if $self->verbose;
 					  }
 				  }
				 );

	  # "debug"			 => \$DEBUG,
	  # "cleanup|clean" => \$CLEANUP,
	  # "save"			 => \$SAVE,
}

1;

__END__
