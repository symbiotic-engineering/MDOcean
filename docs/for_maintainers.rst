For Maintainers
===============
These instructions are intended for SEA Lab members who are developing, maintaining, or co-authoring MDOcean.

**Writing: getting github changes into overleaf**

If simulation code (git), figures (dvc), or tex source code (git) have changed 
and you want the most updated version on Overleaf, use the following commands:

.. code-block:: bash

  git checkout main                    # must be main, not another branch
  calkit pull                          # download data - this may take several minutes
  calkit check pipeline -c             # this updates dvc.yaml if needed
  calkit dvc status -d build-AOR-paper # only proceed to overleaf sync step if this step prints "Data and pipelines are up to date."
  calkit overleaf sync --push-only     # this pushes changes to overleaf


**Writing: making changes on overleaf and getting them onto github**

- Before making any changes on Overleaf, check if a new version has been released on GitHub. If so, perform the steps above to get those changes onto Overleaf.
- Make your changes on Overleaf. Note that any changes you make to ``figs/``, ``tables/``, and ``numeric-results.tex`` on overleaf will be lost when you sync, so changes to those files should be made in Git only.


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
  code mdocean/plots/fig_tab_pub_mapping.m # opens file in vscode
  # edit mdocean/plots/fig_tab_pub_mapping.m to reflect your new figure order
  git add mdocean/plots/fig_tab_pub_mapping.m
  git commit -m 'update figure order'

  # the following steps are only required if you changed the order of any figures in overleaf.
  # These can be performed locally if you have matlab installed, or skipped and they will be
  # automatically performed on CI. If skipping, you should check the CI for your branch
  # (https://github.com/symbiotic-engineering/MDOcean/actions/workflows/calkit-run.yml)
  #  after pushing and confirm that the 'Run Calkit' step succeeds, and revise if not.
  calkit run make-calkit-stages               # this line is only necessary if you added a new figure, not if you just changed the order
  python mdocean/analysis/update_calkit.py    # this line is only necessary if you added a new figure, not if you just changed the order
  calkit run --log build-AOR-paper            # necessary for figure order changes. Only proceed to next step if this step succeeds (green check mark)
  calkit save dvc.lock dvc.yaml .calkit/ calkit.yaml calkit_stages.yaml pubs/applied-ocean-research-model/numeric-results.tex results/**/end.json results/**/*.tex -m "Run pipeline with updated AOR fig order"
          
  # the following steps are always required
  git add pubs/elsarticle-num-names.bst pubs/applied-ocean-research-model/main.tex pubs/applied-ocean-research-model/sections/ pubs/applied-ocean-research-model/numbers.tex
  git commit -m "Update paper from overleaf"
  git push
  

**CI**

MDOcean uses six self-hosted runners on the lab computer for CI.
To create additional runners, follow these instructions (todo).
