use Cwd 'abs_path';
use File::Basename 'dirname';

# method 1: folder of this script
my $root = dirname(__FILE__);

# method 2: use cwd where script is called from
#my $cwd = getcwd; # /work/
#my $root = File::Spec->catfile($cwd, 'pubs', 'applied-ocean-research-model'); # /work/pubs/applied-ocean-research-model

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
