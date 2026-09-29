CI Pipeline Routing
###################

Motivation
==========

Documentation-only changes do not need to compile or test the Spheral source
code. Requiring the full build and test pipeline for these changes consumes
limited CI resources and can make a documentation review-and-update cycle take
days. That delay is a barrier to keeping the documentation current.

The goal of this routing is to recognize a class of documentation-only
changes, validate them with Read the Docs (RTD), and quickly return a
successful GitLab status without running the expensive Spheral build. A branch
that contains any other kind of change continues to run Full CI.

The documentation-only class currently contains:

* All paths under ``docs/``.
* The repository-root ``.readthedocs.yaml``, ``LICENSE``, ``NEWS``, ``NOTICE``,
  ``README.md``, and ``RELEASE_NOTES.md`` files.

This class can be extended when another file does not require source
validation. Adding a path requires updating the classifier and this page.

Deliberate simplifications
==========================

This initial design intentionally does not attempt to:

* Avoid RTD builds when a pull request does not change documentation. RTD is
  comparatively inexpensive, and allowing it to build every pull request
  removes another classification and state-tracking problem.
* Discover and compare against an arbitrary GitHub pull-request target. The
  classifier compares the branch with GitLab's ``CI_DEFAULT_BRANCH``.
* Record the inputs or commit associated with the last successful Full CI or
  RTD build. Classification is stateless and is repeated for each branch
  update.
* Change Hubcast or make GitHub Actions responsible for launching GitLab. The
  existing Hubcast synchronization and status reporting remain in place.
* Select individual source test suites from the paths changed. A branch is
  either docs-only or it receives Full CI.

These choices keep the implementation contained within the GitLab
configuration. They can be revisited if broader path classes, conditional RTD
builds, arbitrary pull-request targets, or more granular source testing become
important.

Overview
========

GitHub is the authoritative repository for Spheral source code and pull
requests. Several services cooperate to validate a change:

* GitHub Actions runs repository-hosted workflows for selected pull-request
  paths and branch pushes.
* Hubcast synchronizes GitHub branch updates to the Spheral project on CZ
  GitLab.
* GitLab runs the Spheral build and test jobs on LC systems. Hubcast reports
  the resulting ``gitlab-ci`` status to the corresponding GitHub commit.
* Read the Docs (RTD) builds the documentation for pull requests and reports
  its status directly to GitHub.

The GitLab configuration gives documentation-only branches a short path
through CI. This avoids waiting for the full Spheral build when a branch only
changes documentation, while still running full CI for any branch containing
a source change.

GitHub Actions
==============

GitHub Actions runs independently of GitLab and RTD. The workflows under
``.github/workflows/`` currently have these responsibilities:

* ``test-tpls.yml`` runs **Test external build** for pull requests that change
  the Spack scripts, developer build scripts, build requirement files, or the
  ``Dockerfile``. GitHub's path filter determines whether this workflow runs.
  It builds the Spheral container images but does not publish them.
* ``docker-image.yml`` runs **Generate build env image** for pushes to
  ``develop`` and ``task/github-actions``. It builds the container images and
  publishes the Spheral image to the GitHub container registry. It is not a
  pull-request workflow, and it currently has no path filter.

These workflows report their results directly through GitHub. The GitLab
docs-only classifier does not enable, disable, or otherwise control them. A
docs-only pull request does not match the paths watched by ``test-tpls.yml``;
however, a later docs-only push to ``develop`` still matches the branch
trigger in ``docker-image.yml``.

GitLab pipeline entry points
============================

GitLab pipelines enter the configuration in three ways:

* An ordinary, non-tag branch push uses the dispatcher and child pipelines
  described below.
* A scheduled pipeline runs the applicable performance, deployment, or
  cleanup jobs directly. It does not use the dispatcher.
* A tag pipeline runs the applicable production jobs directly. It does not
  use the dispatcher.

Branch-push pipelines
=====================

For an ordinary branch push, the first GitLab pipeline is the **CI
dispatcher**. GitLab calls this the *parent pipeline* because its ``run-ci``
trigger job creates a downstream pipeline in the same project. The dispatcher
runs ``classify-ci``, selects one of two configurations, and triggers the
selected pipeline::

  CI dispatcher (parent)
  |
  +-- Full CI child
  |   +-- MPI jobs
  |   +-- sequential jobs
  |   `-- applicable production jobs
  |
  `-- Docs-only CI child
      `-- lightweight success job

``full-ci-child.yml`` reloads the main ``.gitlab-ci.yml`` configuration. In
that child, GitLab sets ``CI_PIPELINE_SOURCE`` to ``parent_pipeline``. The
include rules use this value to activate the normal branch build and test
jobs. The dispatcher jobs themselves only accept a ``push`` source, so they
do not run again inside the child and cannot recursively create more
pipelines.

``docs-only-ci-child.yml`` contains one inexpensive job. Its successful result
means that full Spheral CI was intentionally not required. It does not
validate the documentation; RTD performs that validation separately.

Although GitLab supports a child pipeline that triggers a grandchild, this
workflow uses only one child level.

Branch classification
=====================

The classifier begins with Full CI as the safe default. A feature branch is
classified as docs-only only when every file touched by its commits since it
diverged from ``CI_DEFAULT_BRANCH`` is one of the following:

* A path under ``docs/``.
* The repository-root ``.readthedocs.yaml`` file.

A direct update to the default branch always selects Full CI. Any non-docs
path selects Full CI, including a branch containing both source and
documentation changes. Failure to fetch the default branch, find the merge
base, or inspect the branch history also leaves Full CI selected.

For example, if branch B is created from branch A, and A contains a source
change while B adds only a documentation change, B is not docs-only relative
to the default branch. Full CI is therefore selected for B.

Specialized jobs
================

MPI and sequential jobs are part of the Full CI child. Production job
definitions are also available there, but their individual rules still
determine whether they run. Loading a job definition does not mean that the
job runs for every branch.

Scheduled and tag pipelines remain independent of the dispatcher:

* Performance jobs run in the top-level pipeline for the corresponding
  schedule. An explicit ``test-perf`` commit message can also request them on
  a branch push.
* Deployment and cleanup jobs run in their corresponding scheduled
  pipelines. An explicit ``test-deploy`` commit message can request the
  deployment job on a branch push.
* Tag pipelines make the applicable production jobs available according to
  their existing rules.

Here, **Full CI** means the normal CI applicable to a branch update. It does
not mean that every scheduled, deployment, performance, and tag-only job runs
for every branch.

Status reporting
================

The dispatcher trigger uses ``strategy: depend``, so the parent waits for the
selected child and adopts its result. Hubcast can therefore report the parent
pipeline result to GitHub as the ``gitlab-ci`` check.

GitHub Actions and RTD report their results directly to GitHub. Consequently:

* A docs-only pull request receives the inexpensive GitLab result and the RTD
  documentation result. The existing GitHub pull-request workflow does not
  run because its watched paths are unchanged.
* A source-only or mixed pull request receives the Full CI result and the RTD
  documentation result. It also receives the GitHub **Test external build**
  result when one of that workflow's watched paths changes.

Design assumptions
==================

This routing depends on the following assumptions:

* Hubcast presents synchronized ordinary branch updates to GitLab as
  ``push`` pipelines.
* ``CI_DEFAULT_BRANCH`` is the correct baseline for deciding whether a branch
  is documentation-only.
* The classifier has enough Git history to find the merge base and inspect all
  unique branch commits.
* Parent and child pipelines use the same GitLab project, ref, and commit.
* GitHub Actions continues to apply its own event and path filters
  independently of the GitLab classifier.
* RTD pull-request builds are enabled and report their result directly to
  GitHub.
* New paths that should qualify as documentation-only are added explicitly to
  both the classifier and this documentation.

The file names and job names are project choices. The terms *parent pipeline*
and *child pipeline*, the predefined ``CI_PIPELINE_SOURCE`` variable, and its
``parent_pipeline`` value are GitLab terminology.

Future enhancements
===================

Expand the docs-only path class
-------------------------------

The lightweight path class could include documentation maintained outside
``docs/``. Each addition should be made explicitly in the classifier and
documented in the list under `Motivation`_.

Reuse the last successful Full CI result
----------------------------------------

The current classifier considers every commit unique to a feature branch.
Consequently, once a branch contains a source change, every later update to
that branch selects Full CI, including an update that changes only
documentation.

A future implementation could record the commit associated with the last
successful Full CI run for each branch. If that commit remains an ancestor of
the current branch, the classifier could inspect only the changes made after
it. This would produce the following behavior::

  Source change
  `-- Full CI succeeds and records commit A

  Docs-only changes after A
  `-- Docs-only CI runs

  Source change after A
  `-- Full CI runs again
      +-- success: advance the recorded commit
      `-- failure or cancellation: retain commit A

If no successful commit were recorded, or if a rebase made the recorded
commit no longer an ancestor, classification would fall back to the merge
base with ``CI_DEFAULT_BRANCH``. This feature would require reliable
per-branch state, rules for concurrent pipelines and branch deletion, and a
guarantee that failed or cancelled Full CI runs cannot advance the recorded
commit.

Strengthen RTD validation
-------------------------

RTD currently reports fatal Sphinx build errors to GitHub, but the repository
does not explicitly configure RTD to fail a build for every Sphinx warning. A
future change could enable ``sphinx.fail_on_warning`` in
``.readthedocs.yaml`` after resolving or intentionally suppressing existing
warnings.

Additional documentation checks could include Sphinx's nitpicky mode for
unresolved references and a link-check build for broken links. External-link
checks can fail for temporary network or service problems, so they may be
better suited to a scheduled job or a non-blocking check. Any stricter check
should first be tested against the existing documentation and should produce
a clear failure status on the GitHub pull request.

Refine the GitHub Docker workflows
----------------------------------

``docker-image.yml`` currently builds and publishes container images for every
push to its configured branches, including a push that changes only the
docs-only path class. A future change could add a GitHub path filter so a
docs-only push does not rebuild the images. A mixed push must continue to run
the workflow.

The image-building steps in ``docker-image.yml`` and ``test-tpls.yml`` are
also similar. They could be moved into a reusable workflow or composite action
to prevent the pull-request validation and image-publication implementations
from drifting apart. The refactoring must preserve their different triggers
and permissions: pull-request validation must not publish an image, while the
trusted branch workflow requires permission to publish to the GitHub
container registry.
