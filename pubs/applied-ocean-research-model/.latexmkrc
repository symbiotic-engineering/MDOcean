use Cwd 'abs_path', 'getcwd';
use File::Basename 'dirname';
use File::Spec;

# method 1: folder of this script (works when latexdiff copies this
# file alongside main.tex into a tmp build folder)
my $script_dir = abs_path(dirname(__FILE__));

# method 2: cwd where latexmk is called from (works in the docker
# container, where cwd is the repo root /work/, two folders above main.tex)
my $cwd_root = File::Spec->catfile(getcwd(), 'pubs', 'applied-ocean-research-model');

# use whichever candidate actually contains main.tex
my $root = -e File::Spec->catfile($script_dir, 'main.tex') ? $script_dir : $cwd_root;

# Ensure BibTeX finds elsarticle-num-names.bst
print "BSTINPUTS before: @BSTINPUTS\n";
push @BSTINPUTS, "$root/..";
print "BSTINPUTS after: @BSTINPUTS\n";

# Ensure TeX inputs are found too
print "TEXINPUTS before: @TEXINPUTS\n";
push @TEXINPUTS, $root;
print "TEXINPUTS after: @TEXINPUTS\n";

$out_dir = '.';
$aux_dir = 'aux';

# Enable -shell-escape for pdflatex (required for bibcop override)
set_tex_cmds( '--shell-escape %O %S' );

$max_repeat = 10;

# Suppress noisy glossary output
$makeglossaries = 'makeglossaries -q %O %S';
