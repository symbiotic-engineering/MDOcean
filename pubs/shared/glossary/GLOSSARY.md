# Glossary README

At a high level, how does the glossary work?

- We use the `glossaries-extra` package, which allows a set of glossary entries (including symbols and acronyms) to be defined ahead of time, and then called in-text using the `\gls{}` command.
- The `/gls{}` command looks up the corresponding glossary entry, and inserts the correct symbol/acronym/etc in the text. It also marks that entry as used, which will be relevant for the `\printglossary` command.
- If the `\printglossary` command is used (ie at the beginning of the paper), it prints all glossary entries that have been used in the paper. Due to smart compilation/indexing steps, a glossary can accurately be printed at the beginning before any of the terms have been used, and it will order them alphabetically and include page number references to where each term is used.

Where are the actual glossary entries defined?

- Acronyms are defined in `pubs/shared/glossary/glossary-acronyms.tex` and are shared to all publications
- Symbol definitions that are used in 2+ publications are defined in `pubs/shared/glossary/glossary-symbols-shared.tex` and are shared to all publications
- Symbol definitions unique to a single publication are defined in `pubs/<papername>/glossary/glossary-symbols-<desc>.tex`.

What is the synatx for defining a glossary entry?

- Acronym definition example: `\newabbreviation[type=control]{acr-dp}{DP}{Dynamic Programming}`
- Symbol definition example, standalone (not a parent or a child): `\newsym{vec-hat-F-th}{\vec{\hat{F}}_{\mathrm{th}}}{Complex phasor Thevenin-equivalent force}`
- Symbol definition example, parent (`mu-nu-parent`) with two children (`mu` and `nu`): 
  ```
  \newsym{mu-nu-parent}{\mu,\nu}{Index of constraint in \glsxtrshort{acr-qcqp}}
    \newsym{mu}{\mu}{Index of constraint in \glsxtrshort{acr-qcqp}}[mu-nu-parent]
    \newsym{nu}{\nu}{Index of constraint in \glsxtrshort{acr-qcqp}}[mu-nu-parent]
  ```  

Where is the glossary import/setup/etc defined?

- The `glossaries-extra` package is imported in `pubs/shared/shared-pkg.tex`. This file also imports the setup file (below). Some details of this import are conditional on whether shell escape is enabled (true for local builds, false for arxiv uploads).
- The glossary setup is defined in `pubs/shared/glossary/glossary-setup.tex`. This is where my custom glossary-definition macros live.
- For the AOR journal only, the glossary setup file imports an additional AOR-specific setup file: `pubs/applied-ocean-research-model/glossary/no-appendix-glossary`. This creates a glossary style that does not use page numbers and only includes terms used in the main paper, not terms that only occur in the appendix. It also uses an `alias` entry to make it possible so that when a child glossary term is used, the child symbol appears in-text but the parent symbol is the one that gets indexed so it appears in the printed glossary.

How are glossary files compiled?

- In a local build, the `[automake]` option that `glossaries-extra` is loaded with calls `makeindex` through shell escape. This means it should happen automatically when you compile the document with ie `pdflatex`, `latexmk`, or `calkit latex build`.
- In a local build, if you did not want to use `[automake]`, another option would be to define custom dependencies in the `latekmkrc` and have `latexmk` handle the glossary-making. This was initially tried for the dissertation, but the `automake` option was decided to be simpler.
- In an `ArXiv` build, the glossary files cannot be created, and must be created through a local build and then uploaded to ArXiv. The `calkit latex to-arxiv` command handles this when generating the `arxiv.zip` file. 
