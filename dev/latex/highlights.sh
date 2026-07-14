# uncommented highlights
echo "################################## uncommented highlights ##################################"
grep -r --include="*.tex" -E "^[^%]" pubs/applied-ocean-research-model/ | grep -Fo "\\hl{" | wc -l
grep -rn --include="*.tex" -vE "^\s*%" pubs/applied-ocean-research-model/ | grep -F "\hl{" | sed 'G'

# all highlights
echo "################################## all highlights ##################################"
grep -r --include="*.tex" -Fo "\hl{" pubs/applied-ocean-research-model/ | wc -l
grep -rn --include="*.tex" -F "\hl{" pubs/applied-ocean-research-model/ | sed 'G'