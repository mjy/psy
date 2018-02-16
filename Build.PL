use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Psy::Psy',
    license             => 'perl',
    dist_author         => 'Matt <m{j}yoder@{tee}[aye](em)(you).domain4unis>',
    dist_version_from   => 'lib/Psy/Psy.pm',
    requires => {
        'Test::More' => 0,
    },
    add_to_cleanup      => [ 'Psy-*' ],
);

# add the template toolkit files
$builder->add_build_element('tt');

$builder->create_build_script();

