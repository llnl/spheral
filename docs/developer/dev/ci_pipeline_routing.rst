CI Pipeline Routing
###################

Motivation
==========

Documentation-only changes do not need to compile or test the Spheral source
code. Requiring the full build and test pipeline for these changes consumes
limited CI resources and can make a documentation review-and-update cycle take
days. That delay is a barrier to keeping the documentation current.

The goal of this routing is to recognize a class of documentation-only
changes and quickly return a successful status without running the expensive
Spheral build. A branch that contains any other kind of change continues to
run Full CI.

The documentation-only class currently contains:

* All paths under ``docs/``.
* The repository-root ``.readthedocs.yaml``, ``LICENSE``, ``NEWS``, ``NOTICE``,
  ``README.md``, and ``RELEASE_NOTES.md`` files.

This class can be extended when another file does not require source
validation. Adding a path requires updating the classifier and this page.

Feature scope
=============

Supported features
------------------

* A feature branch containing only documentation-class changes receives
  Docs-only CI.
* A feature branch containing any source or mixed changes receives Full CI.
* Direct updates to the default branch receive Full CI.
* An inconclusive classification safely selects Full CI.

Features not supported for simplicity
-------------------------------------

* Classifying a branch relative to a merge target other than
  ``CI_DEFAULT_BRANCH``.
* Reusing a successful Full CI result after later documentation-only commits
  on a mixed branch.
* Selecting individual Full CI test suites according to the files changed.
* Automatically treating files outside the documented path class as
  documentation-only.

Branch-push pipelines
=====================

For an ordinary branch push, the first GitLab pipeline is the **CI
dispatcher**. It runs ``classify-ci``, selects one of two purpose-specific
configurations, and triggers the selected pipeline::

  CI dispatcher
  |
  +-- Full CI
  |   +-- MPI jobs
  |   +-- sequential jobs
  |   `-- applicable production jobs
  |
  `-- Docs-only CI
      `-- lightweight success job

The configuration is divided by intent:

* ``ci-common.yml`` defines settings shared by the dispatcher and both
  selected pipelines.
* ``ci-dispatch.yml`` contains ``classify-ci`` and ``run-ci``.
* ``ci-full.yml`` contains the normal branch build and test configuration.
* ``ci-docs-only.yml`` contains one inexpensive success job.

Each purpose-specific configuration sets ``SPHERAL_CI_MODE`` to ``full`` or
``docs-only``. The generated configuration includes ``ci-common.yml`` plus
the selected file. It never includes the root ``.gitlab-ci.yml`` or
``ci-dispatch.yml``, so the selected pipeline cannot recursively dispatch
another pipeline.

The Docs-only CI result means that full Spheral CI was intentionally not
required.

GitLab implements the selected pipeline as one child level beneath the
dispatcher. No further pipeline level is created.

Branch classification
=====================

The classifier begins with Full CI as the safe default. A feature branch is
classified as docs-only only when every file touched by its commits since it
diverged from ``CI_DEFAULT_BRANCH`` is one of the following:

* A path under ``docs/``.
* The repository-root ``.readthedocs.yaml``, ``LICENSE``, ``NEWS``, ``NOTICE``,
  ``README.md``, and ``RELEASE_NOTES.md`` files.

A direct update to the default branch always selects Full CI. Any non-docs
path selects Full CI, including a branch containing both source and
documentation changes. Failure to fetch the default branch, find the merge
base, or inspect the branch history also leaves Full CI selected.

For example, if branch B is created from branch A, and A contains a source
change while B adds only a documentation change, B is not docs-only relative
to the default branch. Full CI is therefore selected for B.

Specialized jobs
================

MPI and sequential jobs are part of the Full CI configuration. Production job
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

Result propagation
==================

The dispatcher trigger uses ``strategy: depend``, so the dispatcher waits for
the selected pipeline and adopts its result.

Design assumptions
==================

This routing depends on the following assumptions:

* ``CI_DEFAULT_BRANCH`` is the correct baseline for deciding whether a branch
  is documentation-only.
* The classifier has enough Git history to find the merge base and inspect all
  unique branch commits.
* The dispatcher and selected pipelines use the same GitLab project, ref, and
  commit.
* New paths that should qualify as documentation-only are added explicitly to
  both the classifier and this documentation.

GitLab implements the dispatcher trigger as a parent/child pipeline
relationship. ``SPHERAL_CI_MODE`` records the routing intent independently of
that execution mechanism.
