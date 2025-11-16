use Cwd 'abs_path';
use File::Basename 'dirname';

my $cwd = getcwd; # /work/
my $root = File::Spec->catfile($cwd, 'pubs', 'applied-ocean-research-model'); # /work/pubs/applied-ocean-research-model

# Ensure BibTeX finds references.bib
print "BIBINPUTS before: @BIBINPUTS\n";
push @BIBINPUTS, "$root/..";
print "BIBINPUTS after: @BIBINPUTS\n";

# Ensure BibTeX finds elsarticle-num-names.bst
print "BSTINPUTS before: @BSTINPUTS\n";
push @BSTINPUTS, "$root/..";
print "BSTINPUTS after: @BSTINPUTS\n";

# Optional: ensure TeX inputs are found too
push @TEXINPUTS, $root;

$out_dir = '.';
$aux_dir = 'aux';