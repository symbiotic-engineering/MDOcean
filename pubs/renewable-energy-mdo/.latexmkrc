use Cwd;
my $root = getcwd();  # /workdir/pubs/renewable-energy-mdo

# Ensure BibTeX finds references.bib
push @BIBINPUTS, $root;

# Ensure BibTeX finds elsarticle-num-names.bst
push @BSTINPUTS, "$root/..";

# Optional: ensure TeX inputs are found too
push @TEXINPUTS, $root;
