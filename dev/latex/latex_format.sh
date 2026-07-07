 # !/bin/bash

 # Split sentences so there is one line per sentence, while preserving comments and avoiding splitting at common abbreviations. 
 find . -type f -name "*.tex" -exec perl -i -pe '
  next if /^\s*%/;  # Skip full-line comments
  my ($code, $comment) = split(/(?<!\\)%/, $_, 2);  # Split code and comment (handle escaped %)
  $code =~ s{
    (?<!\b(?:e\.g|i\.e|s\.t|etc|vs|Mr|Ms|Dr|Prof|Fig|Eq|Tab|Int|cf|et\ al|U\.S|No|\\left|\\right|[A-Z]))  # avoid common abbreviations and initials
    ([.!?])                           # sentence-ending punctuation
    (\s+)(?=[A-Z\\])                  # followed by space + capital or command
  }{$1\n}gx;
  $_ = defined $comment ? "$code%$comment" : $code;
' {} +

# Replace spaces after newcommands with non-breaking spaces (~) to prevent line breaks after them.
grep -o '\\newcommand{\\[^}]*}' pubs/renewable-energy-mdo/numbers.tex | sed 's/\\newcommand{//;s/}//' | while read cmd; do
  grep -rl --include="*.tex" "$cmd " . | while read file; do
    sed -i "s/${cmd} /${cmd}~/g" "$file"
  done
done
