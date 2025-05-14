##########
WSL2 Notes
##########

The Windows Subsystem for Linux (WSL) is a useful method of development on Windows 10 based systems.  If going this route we recommend having at least WSL2 for best results -- the original version of WSL (WSL1) also functioned, but is `significantly` slower for jobs such as compilation.

For the most part using an Ubuntu based WSL environment works just using the Ubuntu notes above.  However, one aspect of WSL2 needs to be adjusted.  The build process requires a fair amount of memory (in particular for a few of the Python binding modules), so we recommend having at least 32GB of swap space available.  On WSL2 this is accomplished by creating a ``.wslconfig`` file in your Windows home directory containing at least the following lines::

    [wsl2]
    swap=32GB

