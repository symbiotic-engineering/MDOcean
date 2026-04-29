For Maintainers
===============
These instructions are intended for SEA Lab members who are developing, maintaining, or co-authoring MDOcean.

Writing: pulling changes
If simulation code (git), figures (dvc), or tex source code (git) have changed and you want the most updated version on Overleaf, use the following commands:

For main branch (do this by default):
```
git checkout main
git pull
calkit overleaf sync
```
For another branch (only if you have a specific need):
```
git checkout name-of-branch
calkit pull
calkit run --save
calkit overleaf sync
```

Writing: making/pushing changes
- Before making any changes on Overleaf, do the "pulling changes" steps above.
- Make your changes on Overleaf. When you are done, do the following:
  ```
  git checkout overleaf-changes
  calkit pull
  calkit overleaf status # check that all the changes are as you would expect
  calkit overleaf sync
- Any changes made to `references.bib` on Overleaf need to be handled specially because of this file is used in a `map-paths` stage in Calkit. If you make any changes to this file,  

CI
MDOcean uses six self-hosted runners on the lab computer.
To create additional runners, follow these instructions
