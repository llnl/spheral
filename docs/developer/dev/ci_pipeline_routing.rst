Full and Docs-Only CI Routing
#############################

Motivation
==========

Documentation-only changes do not need to compile or test the Spheral source
code. Requiring the full build and test pipeline for these changes consumes
limited CI resources and can make a documentation review-and-update cycle take
days. That delay is a barrier to keeping the documentation current.

Spheral therefore provides an explicit branch-naming convention for selecting
a lightweight pipeline. A branch whose name begins with ``docs/`` receives
Docs-Only CI. Every other branch receives Full CI.

