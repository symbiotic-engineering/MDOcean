For Maintainers
===============
These instructions are intended for SEA Lab members who are developing, maintaining, or co-authoring MDOcean.

** Writing: getting github changes into overleaf **

If simulation code (git), figures (dvc), or tex source code (git) have changed 
and you want the most updated version on Overleaf, use the following commands:

.. code-block:: bash

  git checkout main                    # must be main, not another branch
  calkit pull                          # download data - this may take several minutes
  calkit check pipeline -c             # this updates dvc.yaml if needed
  calkit dvc status -d build-AOR-paper # only proceed to overleaf sync step if this step prints "Data and pipelines are up to date."
  calkit overleaf sync --push-only     # this pushes changes to overleaf


** Writing: making changes on overleaf and getting them onto github **

- Before making any changes on Overleaf, do the "pulling changes" steps above.
- Make your changes on Overleaf. Note that any changes you make to `figs/`, `tables/`, and `numeric-results.tex` on overleaf will be lost when you sync.


- When you are done, do the following:

.. code-block:: bash

  git checkout main                    # must be main, not another branch
  calkit pull                          # download data - this may take several minutes
  calkit check pipeline -c             # this updates dvc.yaml if needed
  calkit dvc status -d build-AOR-paper # only proceed to overleaf sync step if this step prints "Data and pipelines are up to date."
  calkit overleaf status               # check that all the changes are as you would expect
  calkit overleaf sync --no-commit     # this will make overleaf changes locally, they are not saved to git yet
  # if the sync produced any merge conflicts, resolve them in the editor before proceeding.
  git checkout -b another-old-or-new-branch # omit the -b if it's an existing branch instead of a new branch.

  # the following steps are only required if you made any changes to references.bib in overleaf
  cp mdocean/pubs/applied-ocean-research-model/references.bib mdocean/pubs/renewable-energy-mdo/references.bib 
  git add mdocean/pubs/renewable-energy-mdo/references.bib 
  git commit -m 'update references.bib from AOR overleaf'

  # the following steps are only required if you changed the order of any figures in overleaf
  # edit mdocean/plots/fig_tab_pub_mapping.m to reflect your new figure order
  git add mdocean/plots/fig_tab_pub_mapping.m
  git commit -m 'update figure order'
  calkit run --log build-AOR-paper # only proceed to next step if this step succeeds (green check mark)
  git add .calkit pubs/applied-ocean-research-model/numeric-results.tex

  # the following steps are always required
  git add 
  git push
  
- Any changes made to `references.bib` on Overleaf need to be handled specially because of this file is used in a `map-paths` stage in Calkit. If you make any changes to this file,  

** CI **

MDOcean uses six self-hosted runners on the lab computer.
To create additional runners, follow these instructions
