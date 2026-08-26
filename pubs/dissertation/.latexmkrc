use Cwd 'abs_path';
use File::Basename 'dirname';

my $cwd = getcwd; # /work/
my $root = File::Spec->catfile($cwd, 'pubs', 'dissertation'); # /work/pubs/dissertation

# Ensure BibTeX finds elsarticle-num-names.bst
print "BSTINPUTS before: @BSTINPUTS\n";
push @BSTINPUTS, "$root/..";
print "BSTINPUTS after: @BSTINPUTS\n";

# Ensure TeX inputs are found too
print "TEXINPUTS before: @TEXINPUTS\n";
push @TEXINPUTS, "$root/..", "$root/../../mdocean/simulation/modules/OpenFLASH/pubs/JFM";
print "TEXINPUTS after: @TEXINPUTS\n";

# Ensure bib files are found by latexmk
print "BIBINPUTS before: @BIBINPUTS\n";
push @BIBINPUTS, "$root/", "$root/../shared", "$root/../applied-ocean-research-model", "$root/../renewable-energy-mdo", "$root/../../mdocean/simulation/modules/OpenFLASH/pubs/JFM";
print "BIBINPUTS after: @BIBINPUTS\n";

# Ensure bib files are found by BibTeX
$ENV{'BIBINPUTS'} = join(
    ':',
    "$root",
    "$root/../shared",
    "$root/../applied-ocean-research-model",
    "$root/../renewable-energy-mdo",
    "$root/../../mdocean/simulation/modules/OpenFLASH/pubs/JFM",
    $ENV{'BIBINPUTS'} // ''
);

$out_dir = '.';
$aux_dir = 'aux';

# Enable -shell-escape for pdflatex (required for TikZ externalization)
set_tex_cmds( '--shell-escape %O %S' );

add_cus_dep('glo', 'gls', 0, 'run_makeglossaries');
add_cus_dep('glo-abr', 'gls-abr', 0, 'run_makeglossaries');
add_cus_dep('acn', 'acr', 0, 'run_makeglossaries');
add_cus_dep('cld', 'cln', 0, 'run_makeglossaries');
add_cus_dep('mdd', 'mdn', 0, 'run_makeglossaries');
add_cus_dep('ord', 'orn', 0, 'run_makeglossaries');
add_cus_dep('std', 'stn', 0, 'run_makeglossaries');
add_cus_dep('ecd', 'ecn', 0, 'run_makeglossaries');
add_cus_dep('dyd', 'dyn', 0, 'run_makeglossaries');
add_cus_dep('opd', 'opn', 0, 'run_makeglossaries');
add_cus_dep('swd', 'swn', 0, 'run_makeglossaries');
add_cus_dep('wad', 'wan', 0, 'run_makeglossaries');
add_cus_dep('sad', 'san', 0, 'run_makeglossaries');
add_cus_dep('slo', 'sls', 0, 'run_makeglossaries');

# sub run_makeglossaries {
#     return system("makeglossaries -d '$aux_dir' sampleThesis");
# }
$clean_ext .= " acr acn alg glo gls glg glo-abr gls-abr cld cln mdd mdn ord orn std stn ecd ecn dyd dyn opd opn swd swn wad wan sad san slo sls";

sub run_makeglossaries {
    my ($base_name, $path) = fileparse( $_[0] );
    my @args = ( "-q", "-d", $path, $base_name );
    if ($silent) { unshift @args, "-q"; }
    return system "makeglossaries", "-d", $path, $base_name; 
}
$max_repeat = 10;
