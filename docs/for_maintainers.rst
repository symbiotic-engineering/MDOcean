For Maintainers
===============
These instructions are intended for SEA Lab members who are developing, maintaining, or co-authoring MDOcean.

**One-time setup**

.. code-block:: bash
  #  first install calkit using the appropriate instructions for your system:
  # see https://docs.calkit.org/installation/. You should set up cloud integration.
  git clone --recursive https://github.com/symbiotic-engineering/MDOcean.git

**Writing: getting github changes into overleaf**

If simulation code (git), figures (dvc), or tex source code (git) have changed 
and you want the most updated version on Overleaf, use the following commands:

.. code-block:: bash

  git checkout main                    # must be main, not another branch
  git pull
  git checkout -b overleaf-<desc>      # create new branch called overleaf-<your description here>. Required because main is a protected branch.
  git push --set-upstream origin overleaf-<desc>
  calkit pull -f                       # download data - this may take several minutes
  calkit check pipeline -c             # this updates dvc.yaml if needed
  calkit dvc status -d build-AOR-paper # only proceed to overleaf sync step if this step prints "Data and pipelines are up to date."
  calkit overleaf sync --push-only pubs/applied-ocean-research-model    # this pushes changes to overleaf, and creates/pushes an overleaf sync git commit.


**Writing: making changes on overleaf and getting them onto github**

- Before making any changes on Overleaf, check if a new version has been released on GitHub to the ``main`` branch. If so, perform the steps above to get those changes onto Overleaf.
- Make your changes on Overleaf. Note that any changes you make to ``figs/``, ``tables/``, and ``numeric-results.tex`` on overleaf will be lost when you sync since they are pipeline outputs specified as ``push_paths``, so changes to those files should be made in Git only.


- When you are done, do the following:

.. code-block:: bash

  git checkout main                    # must be main, not another branch
  calkit pull -f                       # download data - this may take several minutes
  calkit check pipeline -c             # this updates dvc.yaml if needed
  calkit dvc status -d build-AOR-paper # only proceed to overleaf sync step if this step prints "Data and pipelines are up to date."
  calkit overleaf status pubs/applied-ocean-research-model               # check that all the changes are as you would expect
  calkit overleaf sync --no-commit pubs/applied-ocean-research-model     # this will make overleaf changes locally, they are not saved to git yet
  # if the sync produced any merge conflicts, resolve them in the editor before proceeding.
  git checkout -b another-old-or-new-branch # omit the -b if it's an existing branch instead of a new branch.
  git push -u origin HEAD              # only necessary if it's a new branch

  # The following steps are only required if you made any changes to references.bib in overleaf,
  # since it is the output of a map-paths calkit stage.
  cp mdocean/pubs/applied-ocean-research-model/references.bib mdocean/pubs/renewable-energy-mdo/references.bib 
  git add mdocean/pubs/renewable-energy-mdo/references.bib 
  git commit -m 'update references.bib from AOR overleaf'

  # The following steps are only required if you changed the order of any figures in overleaf
  code mdocean/plots/fig_tab_pub_mapping.m # opens file in vscode
  # edit mdocean/plots/fig_tab_pub_mapping.m to reflect your new figure order
  git add mdocean/plots/fig_tab_pub_mapping.m
  git commit -m 'update AOR figure order'

  # The following steps are only required if you changed the order of any figures in overleaf.
  # These can be performed locally if you have matlab and docker installed, or skipped and they will be
  # automatically performed on CI. If skipping, you should check the CI for your branch
  # (https://github.com/symbiotic-engineering/MDOcean/actions/workflows/calkit-run.yml)
  #  after pushing and confirm that the 'Run Calkit' step succeeds, and revise if not.
  calkit run make-calkit-stages               # this line is only necessary if you added a new figure, not if you just changed the order. May take a few minutes.
  python mdocean/analysis/update_calkit.py    # this line is only necessary if you added a new figure, not if you just changed the order
  calkit run --log build-AOR-paper            # necessary for figure order changes. May take a few minutes. Only proceed to next step if this step succeeds (green check mark)
  calkit save dvc.lock dvc.yaml .calkit/ calkit.yaml calkit_stages.yaml pubs/applied-ocean-research-model/numeric-results.tex results/**/end.json results/**/*.tex -m "Run pipeline with updated AOR fig order" # pipeline definition and git-tracked pipeline outputs
          
  # the following steps are always required
  git add pubs/elsarticle-num-names.bst pubs/applied-ocean-research-model/main.tex pubs/applied-ocean-research-model/sections/ pubs/applied-ocean-research-model/numbers.tex # git-tracked pipeline inputs
  git commit -m "Update paper from Overleaf sync"
  git push
  

**CI**

MDOcean uses six self-hosted runners on the lab computer for CI.
To create additional runners, follow these instructions (todo).
