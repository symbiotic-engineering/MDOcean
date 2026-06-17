use Cwd 'abs_path';
use File::Basename 'dirname';

my $cwd = getcwd; # /work/
my $root = File::Spec->catfile($cwd, 'pubs', 'renewable-energy-mdo'); # /work/pubs/renewable-energy-mdo

# Ensure BibTeX finds references.bib
push @BIBINPUTS, $root;
print "BIBINPUTS after: @BIBINPUTS\n";

# Ensure BibTeX finds elsarticle-num-names.bst
push @BSTINPUTS, "$root/..";
print "BSTINPUTS after: @BSTINPUTS\n";

# Optional: ensure TeX inputs are found too
push @TEXINPUTS, $root;

$out_dir = '.';
$aux_dir = 'aux';

# Enable -shell-escape for pdflatex (required for bibcop override)
set_tex_cmds( '--shell-escape %O %S' );

$max_repeat = 10;