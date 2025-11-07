## Releasing
- Make sure all tests pass
- Update version and date-released in CITATION.cff
- Add changes to CHANGELOG.md
- Merge PR
- In the command line, make a tag with `git tag -a vX.X -m 'message'` and `git push origin --tags`
- On GitHub, make a release from the tag. This should automatically trigger a Zenodo update.
