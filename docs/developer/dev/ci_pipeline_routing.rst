Full and Docs-Only CI Routing
#############################

Motivation
==========

Documentation-only changes do not need to compile or test the Spheral source
code. Requiring the full build and test pipeline for these changes consumes
limited CI resources and can make a documentation review-and-update cycle take
days. That delay is a barrier to keeping the documentation current.

Spheral therefore provides an explicit branch-naming convention for selecting
a lightweight pipeline. A branch whose name contains ``docs-only`` receives
Docs-Only CI. Every other branch receives Full CI.

Features
========

Supported features
------------------

* A branch whose name contains ``docs-only`` receives Docs-Only CI.
* Every other branch receives Full CI.
* Both configurations produce a single, flat GitLab pipeline.

Features not supported for simplicity
-------------------------------------

* GitLab does not inspect changed paths to choose a pipeline.
* The contents of a ``docs-only`` branch are not validated. The branch name is
  trusted as the developer's declaration that Full CI is unnecessary.
* A Docs-Only CI pipeline cannot switch itself to Full CI. Rename the branch
  without ``docs-only`` and push it when Full CI is required.
* Full CI results are not reused by later pipelines.

Configuration
=============

The root ``.gitlab-ci.yml`` uses ``CI_COMMIT_BRANCH`` and conditional includes
to select the jobs in one flat pipeline. When the branch name contains
``docs-only``, it includes ``.gitlab/docs-only.yml`` and omits the MPI,
sequential, and production job definitions. Otherwise, the existing rules for
those Full CI configurations apply without modification.

Docs-Only CI
============

``.gitlab/docs-only.yml`` defines one inexpensive ``docs-only-ci`` success
job. Its result records that the branch intentionally did not request the
Spheral build and test jobs.

The branch name is the only routing signal. For example,
``mcfadden8/docs-only/update-user-guide`` selects Docs-Only CI. A mixed branch
must not use ``docs-only`` in its name unless the developer has determined that
Full CI is unnecessary.

Full CI configuration
=====================

The existing MPI, sequential, and production configurations provide the normal
Full CI jobs.
