 find . -type f -name "*.tex" -exec perl -i -pe '
  next if /^\s*%/;  # Skip full-line comments
  my ($code, $comment) = split(/(?<!\\)%/, $_, 2);  # Split code and comment (handle escaped %)
  $code =~ s{
    (?<!\b(?:e\.g|i\.e|etc|vs|Mr|Ms|Dr|Prof|Fig|Eq|Tab|cf|et\ al|U\.S))  # avoid common abbreviations
    ([.!?])                           # sentence-ending punctuation
    (\s+)(?=[A-Z\\])                  # followed by space + capital or command
  }{$1\n}gx;
  $_ = defined $comment ? "$code%$comment" : $code;
' {} +

grep -o '\\newcommand{\\[^}]*}' paper/numbers.tex | sed 's/\\newcommand{//;s/}//' | while read cmd; do
  grep -rl --include="*.tex" "$cmd " . | while read file; do
    sed -i "s/${cmd} /${cmd}~/g" "$file"
  done
done
