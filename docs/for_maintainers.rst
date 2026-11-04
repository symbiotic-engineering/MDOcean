For Maintainers
===============
These instructions are intended for SEA Lab members who are developing, maintaining, or co-authoring MDOcean.

One-time Setup
**************
First install calkit using the appropriate instructions for your system:
see https://docs.calkit.org/installation/. You should set up cloud integration.
Then set up git:

.. code-block:: bash

  git clone --recursive https://github.com/symbiotic-engineering/MDOcean.git # this downloads MDOcean and its submodules
  git config merge.dvclock.driver "python dev/dvc_lock_merge.py %O %A %B"

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

**How to create a self-hosted runner:**

1. Follow the `linked <https://docs.github.com/en/actions/how-tos/manage-runners/self-hosted-runners/add-runners/>`_
   GitHub instructions to add a self-hosted runner to the MDOcean repo.

2. Follow the `linked <https://docs.github.com/en/actions/how-tos/manage-runners/self-hosted-runners/configure-the-application/>`_
   GitHub instructions to configure the self-hosted runner as a service,
   which means it will automatically start the runner every time the computer is restarted
   so you don't have to worry about doing that manually.
   You will need to be careful about sudo and non-sudo users.
   To run the commands, you need to be logged into an account with ``sudo`` privileges,
   but you want to install the service using the ``user`` argument
   corresponding to whichever user you use for development (typically non-sudo).
   On the Cornell lab computer, this meant doing ``su - cornell`` to switch to a user with sudo
   before running the commands in the instructions,
   and using ``becca`` as the ``USERNAME`` argument in instruction number 3 of the "Installing the Service" section of the link.

3. Make sure all the needed MATLAB dependencies are installed.
   These are listed in the ``README.md`` and can be confirmed by running ``tests/dep_check.m``.

4. Make sure MATLAB is on the PATH. On the Cornell lab computer,
   this meant adding the following line to ``.github/workflows/calkit-run.yml``:
   ``echo "/usr/local/MATLAB/R2025a/bin" >> $GITHUB_PATH``.

5. Check if it works: Create a new branch and make an arbitrary change --
   a good example that shouldn't take too long to run would be to change the formatting of the plot in
   ``mdocean/plots/matlab_figs/module_demos/spar_hydro_plot.m``.
   Push the branch and open a PR into ``main``. Go to `GitHub Actions
   <https://github.com/symbiotic-engineering/MDOcean/actions/workflows/calkit-run.yml>`_
   and confirm that the workflow ran and did not encounter any errors.
   After it finishes running, you should see a commit automatically made to your branch by the actions bot.
   The commit should include an update to the Jupyter Notebook showing your plot update, as well as a log and a change to ``dvc.lock``.
   The notebook changes can be easily confirmed visually with
   `ReviewNB <https://app.reviewnb.com/symbiotic-engineering/MDOcean/>`_.

6. Repeat the above steps multiple times to create several runners to allow CI to run on multiple branches at a time.
   For example, on the Cornell lab computer, I set up 6 runners,
   in folders named ``actions-runner``, ``actions-runner-2``, ..., ``actions-runner-6``.


**How to restart the self-hosted runner (can be necessary if all CI is timing out):**

1. log into ssh as ``becca`` user (as normal)

2. ``su - cornell`` to switch to user with sudo

3. ``cd /home/becca/Documents/git/actions-runner`` to change back to req'd directory

4. ``sudo ./svc.sh status`` to show status if desired

5. ``sudo ./svc.sh stop && sudo ./svc.sh start`` to restart.
   Note that you may need to modify the stop script to add the line
   ``pkill -u becca -f matlab`` to get rid of any leftover MATLAB processes under your user.
