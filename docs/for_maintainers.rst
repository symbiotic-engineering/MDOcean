For Maintainers
===============
These instructions are intended for SEA Lab members who are developing, maintaining, or co-authoring MDOcean.

One-time Setup
**************
First install calkit using the appropriate instructions for your system:
see https://docs.calkit.org/installation/. You should set up cloud integration.
Then set up git:

.. code-block:: bash

  git clone https://github.com/symbiotic-engineering/MDOcean.git
  git config merge.dvclock.driver "python dev/dvc_lock_merge.py %O %A %B"

Then, in MATLAB, install the external MATLAB dependencies (run once, from the
repository root):

.. code-block:: matlab

  setup_mip

Writing
*******

.. note::
   If you get ``ERROR: failed to pull data from the cloud - <number> files failed to download``
   or ``ERROR: unexpected error`` on any of the the ``calkit pull`` commands, repeat the command until success. 
   You can safely ignore ``ERROR: failed to pull data from the cloud - Checkout failed for following targets``.

.. note::
   When following these instructions, fill in in ``<paper-stage>`` with either ``build-AOR-paper`` or ``build-RE-paper``, and ``<paper-folder>`` with either ``applied-ocean-research-model`` or ``renewable-energy-mdo`` (no slashes).

.. note::
   For a guided interactive version of these Overleaf sync workflows, run
   ``bash dev/latex/calkit_instructions_interactive.sh`` from the repository root.

**Getting github changes into overleaf**

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


**Making changes on overleaf and getting them onto github**

- Before making any changes on Overleaf, check if a new version has been released on GitHub to the ``main`` branch. If so, perform the steps above to get those changes onto Overleaf.
- Make your changes on Overleaf. Note that any changes you make to ``figs/``, ``tables/``, and ``numeric-results.tex`` on overleaf will be lost when you sync since they are pipeline outputs specified as ``push_paths``, so changes to those files should be made in Git only.
- By default, when you make a change on overleaf, it will apply to both the dissertation version and the journal version. If you want different versions of the text for the journal and dissertation, use the following latex code structure:

.. code-block:: latex

  put content that goes in both dissertation and journal here
  \ifdefined\DISSERTATION
     put content that goes in dissertation here
  \else
     put content that goes in journal here
  \fi
  put more content that goes in both dissertation and journal here


- When you are done, do the following steps:

.. code-block:: bash

  git checkout main                    # must be main, not another branch
  calkit pull -f                       # download data - this may take several minutes
  calkit check pipeline -c             # this updates dvc.yaml if needed
  calkit dvc status -d <paper-stage> # only proceed to overleaf sync step if this step prints "Data and pipelines are up to date."
  calkit overleaf status pubs/<paper-folder>               # check that all the changes are as you would expect
  calkit overleaf sync --no-commit pubs/<paper-folder>     # this will make overleaf changes locally, they are not saved to git yet
  # if the sync produced any merge conflicts, resolve them in the editor before proceeding.
  git checkout -b another-old-or-new-branch # omit the -b if it's an existing branch instead of a new branch. Typically this is overleaf-<desc> from above.
  git push -u origin HEAD              # only necessary if it's a new branch

  # The following steps are only required if you made any changes to references.bib, shared-pkg.tex, 
  # or elsarticle-num-names.bst in overleaf, since they are outputs of a map-paths calkit stage.
  cp pubs/<paper-folder>/<filename> pubs/shared/<filename> # replace <filename> with whichever file edited
  git add pubs/shared/<filename>
  git commit -m 'update shared <filename> from overleaf'

  # The following steps are only required if you changed the order of any figures in overleaf
  code mdocean/plots/fig_tab_pub_mapping.m # opens file in vscode
  # edit mdocean/plots/fig_tab_pub_mapping.m to reflect your new figure order
  git add mdocean/plots/fig_tab_pub_mapping.m
  git commit -m 'update figure order'

  # The following steps are only required if you changed the order of any figures in overleaf.
  # These can be performed locally if you have matlab and docker installed, or skipped and they will be
  # automatically performed on CI. If skipping, you should check the CI for your branch
  # (https://github.com/symbiotic-engineering/MDOcean/actions/workflows/calkit-run.yml)
  #  after pushing and confirm that the 'Run Calkit' step succeeds, and revise if not.
  calkit run make-calkit-stages               # this line is only necessary if you added a new figure, not if you just changed the order. May take a few minutes.
  python mdocean/analysis/update_calkit.py    # this line is only necessary if you added a new figure, not if you just changed the order
  calkit run --log <paper-stage>              # necessary for figure order changes. May take a few minutes. 
  # Only proceed to next step if the previous step  succeeded (green check mark)
  calkit save dvc.lock dvc.yaml .calkit/ calkit.yaml calkit_stages.yaml \
     pubs/<paper-folder>/numeric-results.tex results/**/end.json results/**/*.tex \
     -m "Run pipeline with updated fig order" # pipeline definition and git-tracked pipeline outputs
          
  # the following steps are always required
  git add pubs/<paper-folder>/*.tex pubs/applied-ocean-research-model/sections/ # git-tracked pipeline inputs
  git commit -m "Update paper from Overleaf sync"
  git push



CI
**

MDOcean uses six self-hosted runners on the lab computer for CI.
To create additional runners, follow these instructions (todo).
