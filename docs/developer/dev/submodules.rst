Submodules
##########

Spheral uses submodules to bring other codes/project directly into our source tree, such as BLT, and PYB11Generator. This guide explains some of the common practices one might want to perform with submodules during Spheral development.

.. _update_local_submodules:

Update Local Submodules
=======================

On occassion developers will need to update the submodules of their local Spheral repository. This is common when merging one branch into another. This will be required when a submodule has been bumped to a new verison or different branch. To update your local submodules you will want to run these commands from the Spheral root directory::

  git submodule sync --recursive
  git submodule update --recursive --init


Working Branch of a Submodule
=============================

Developer may want to perform work on one of our submodules to either do development work or fix bugs. When doing this one might want to create a working branch for the given submodule.

 - Enter submodule directoy.::

     cd <submodule_dir>

 - git create a working branch::

     git checkout -b <working_branch_name> 

 - commit and push changes to submodule working branch.::

     git commit -m "<message>"
     git push

 - return to Spheral root directory.::

     cd <Spheral_root_dir>
   
 - Update Spheral to point at submodules new branch.::

     git add <submodule> 
     git commit -m "Updating <submodule> to <working_branch_name>"
     git push

Change Branch of a Submodule
============================

A developer may want to check out a different branch of a submodule to try it out. This might be done to test a bugfix branch or point at a release tag / candidate branch for that submodule. The steps a fairly similar to those above.

 - Enter submodule directoy.::

     cd <submodule_dir>

 - git create a working branch::

     git pull
     git checkout <working_branch_name> 

 - return to Spheral root directory.::

     cd <Spheral_root_dir>
   
 - Update Spheral to point at submodules new branch.::

     git add <submodule> 
     git commit -m "Updating <submodule> to <working_branch_name>"
     git push

Submodule Origin
================

To Change the origin of a submodule you must first edit the ``url`` line of the respective submodule in the ``.gitmodules`` file. 
Next you will need to run the steps from :ref:`Update Local Submodules`. You can then work on the submodule or switch branch normally as detailed in the above sections.
