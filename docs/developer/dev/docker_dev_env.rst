*******************************************
Docker Development Environment
*******************************************

Spheral builds an up-to-date docker container for every merge-request
on ``develop``. Developers can use this container to do development tasks
on local machines.

===========================
Creating a Dev Environment
===========================

We will use ``docker dev create`` with our Spheral docker image and a
local repository. This will allow us to skip setting up a linux system with
external packages, gives us pre-built TPLs and allows us to edit a cloned
repository from our local machines IDE/text editor.bash::

  > rm <path_to_local_repo>/compose-dev.yaml
  > docker dev create --base-image ghcr.io/llnl/spheral --name <name_your_env> --path <path_to_local_repo> -o <path_to_local_repo>

.. note::
   You need to have **Docker Desktop**, **VSCode** and the **VSCode Dev Environment Extension** installed for this to work. You do not need to use VSCode to access the container, but the extension seems to do some of the lifting for us when setting up the volume to our local repo.

**Output** ::

  spheral-recursing_darwin <---- Name of dev environment
  Creating Dev Environment "spheral-recursing_darwin"
  populating volume from /Users/davis291/Projects/spheral
  Creating Dev Environment "spheral-recursing_darwin"
  detecting language
  Detecting main repo language...
  building compose stack
  building compose stack
  starting compose stack
  starting compose stack
   Network spheral-recursing_darwin_default  Creating
   Network spheral-recursing_darwin_default  Created
   Container spheral-recursing_darwin-app-1  Creating
   Container spheral-recursing_darwin-app-1  Created
   Container spheral-recursing_darwin-app-1  Starting
   Container spheral-recursing_darwin-app-1  Started <---- Name of running container to connect to.
  Dev Environment "spheral-recursing_darwin" (5bd37219d27eb68a77ce6fd8fee05a533a52017d8dcc72430867e2471e428e58) is running!%


=============================
Connecting to a Dev Container
=============================

Once the container has been started, you can connect directly through the terminal
with the **Container** name (**NOT** the **Dev Environment** name).::

  > docker exec -it spheral-recursing_darwin-app-1 /bin/bash
  root@671dab5d0b00:/home/spheral/workspace/build_docker-gcc/install#

This drops you into the install location of the ``spheral@develop`` build from
github, this is a fully installed version of the latest ``develop`` Spheral.

.. tip::
  VSCode & Docker Desktop:
    * Open **Docker Desktop** and navigate to the **Dev Environment** tab.
    * Find the container name and select **OPEN IN VSCODE**.


=============================
Development Work
=============================

Your local Spheral repo is mounted from your local filesystem. You can develop directly from your
IDE or text editor of choice. Then you can compile and run from within the container itself.

- The local Spheral repository will be mounted in the container at ``/com.docker.devenvironments.code/``.

- There already exists a full build and install of Spheral at ``develop`` in ``/home/spheral/workspace/build_docker-gcc/install``.

- An updated host config file can be found at ``/home/spheral/wokspace/docker-gcc.cmake``.
