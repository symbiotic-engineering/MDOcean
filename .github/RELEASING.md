## Releasing
- Make sure all tests pass. The ready to release workflow makes sure you have updated `__version__.py`, `CHANGELOG.md`, and haven't accidentally downgraded a submodule (a common mistake when merging main into your branch).
- Git tag, GitHub release, and Zenodo release will all happen automatically when the PR is merged to main
- If you want to tag a branch manually, use `git tag -a vX.X -m 'message'` and `git push origin --tags`