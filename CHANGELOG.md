# Changelog


All notable changes to this project will be documented in this file. The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). This project follows the versions of [SMASH](https://github.com/smash-transport/smash), that means that each new version of SMASH triggers a new version of this project, in order to assure compatibility.
This changelog is in place since version SMASH-analysis-3.0.

The major categories to group changes in this log are:

* `Input / Output` for all, in particular breaking, changes, fixes and additions to the in- and output files.
* `Added` for new features.
* `Changed` for changes in existing functionality.
* `Fixed` for any bug fixes.
* `Removed` for now removed features.

Also possible, but for this project less relevant, is `Deprecated` for soon-to-be removed features.


## Unreleased

## SMASH-analysis-3.0
Date:

### Added
* Energy scan target for 5.02 TeV
* New binary reader version due to added baryon number output
* NNbar-to-5 pions target in the detailed balance analysis

### Fixed
* Correct unstable flag in the cross sections target for NN->NNeta
* Increase equilibration time for detailed balance

### Changed
* ⚠️ The `master` branch has been renamed to `main`
* Use MutableSet from collections.abs instead of collections
