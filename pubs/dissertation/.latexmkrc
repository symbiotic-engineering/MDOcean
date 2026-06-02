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

# Ensure bib files are found too
print "BIBINPUTS before: @BIBINPUTS\n";
push @BIBINPUTS, "$root/", "$root/../../mdocean/simulation/modules/OpenFLASH/pubs/JFM";
print "BIBINPUTS after: @BIBINPUTS\n";

$out_dir = '.';
$aux_dir = 'aux';

# Enable -shell-escape for pdflatex (required for TikZ externalization)
set_tex_cmds( '--shell-escape %O %S' );

$max_repeat = 10;
