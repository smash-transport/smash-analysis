# Changelog


All notable changes to this project will be documented in this file. The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). This project follows the versions of [SMASH](https://github.com/smash-transport/smash), that means that each new version of SMASH triggers a new version of this project, in order to assure compatibility.
This changelog is in place since version SMASH-analysis-3.0. The naming convention changed from this tag on.

The major categories to group changes in this log are:

* :left_right_arrow: for all, in particular breaking, changes, fixes and additions to the in- and output files.
* :heavy_plus_sign: for new features.
* :recycle: for changes in existing functionality.
* :sos: for any bug fixes.
* :heavy_minus_sign: for now removed features.


Also possible, but for this project less relevant, is `Deprecated` for soon-to-be removed features.


## Unreleased

## SMASH-analysis-3.1
Date: 2024-02-29

* :left_right_arrow: Collision term for boxes uses now bottom up cross sections and pseudoresonances are deactivated
* :heavy_plus_sign: New binary reader version due to added strangeness output

## SMASH-analysis-3.0
Date: 2023-04-28

* :heavy_plus_sign: Energy scan target for 5.02 TeV
* :heavy_plus_sign: New binary reader version due to added baryon number output
* :heavy_plus_sign: NNbar-to-5 pions target in the detailed balance analysis
* :sos: Correct unstable flag in the cross sections target for NN->NNeta
* :sos: Increase equilibration time for detailed balance of N*(1440)
* :recycle: ⚠️ The `master` branch has been renamed to `main`
* :recycle: ⚠️ Development is now performed on the `develop` branch
* :recycle: Use MutableSet from collections.abc instead of collections

[Link to diff from previous version](https://github.com/smash-transport/smash-analysis/compare/SMASH-2.2ana...SMASH-analysis-3.0)

## SMASH-2.2ana
Date: 2022-05-13
