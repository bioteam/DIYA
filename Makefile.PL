use ExtUtils::MakeMaker;
# perl Makefile.PL INSTALL_BASE=/usr/local
WriteMakefile(
    NAME         => 'diya',
    VERSION_FROM => 'lib/diya.pm',
    ABSTRACT     => 'Do It Yourself Annotation',
    AUTHOR       => 'Andy Stewart <andrew.stewart@med.navy.mil>',
    PREREQ_PM    => {
        'Bio::Perl'      => '1.005002',
        'XML::Simple'    => '2.0',
        'Data::Merger'   => '0',
        'File::Basename' => '2.0',
        'Getopt::Long'   => '0',
        'FileHandle'     => '0'
    }
);


