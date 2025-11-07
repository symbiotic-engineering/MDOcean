use Cwd;
my $root = getcwd();  # /workdir/pubs/applied-ocean-research-model

# Ensure BibTeX finds references.bib
push @BIBINPUTS, "$root/..";

# Ensure BibTeX finds elsarticle-num-names.bst
push @BSTINPUTS, "$root/..";

# Optional: ensure TeX inputs are found too
push @TEXINPUTS, $root;

$out_dir = 'out';
$aux_dir = 'aux';
