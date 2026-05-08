use Cwd;
my $root = getcwd();  # /workdir/pubs/renewable-energy-mdo
print "Current working directory: $root\n";

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