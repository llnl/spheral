trap 'echo "# $BASH_COMMAND"' DEBUG

run_cmd() {
    if ! eval $1 ; then
        echo "FAILED"
        exit 1
    fi
}

### Create a tar file called $DEV_PKG_NAME.tar.gz containing:
# dev-pkg/
#   *cloned spheral repo.*
#   resources/
#     pip_cache/
#     mirror/
#     build_cache/
#     bootstrap/
#       metadata/
#       bootstrap_cache/

# Also creates a tar file called $DEV_PKG_NAME-spack.tar.gz containing:
# spack
# packages
###############################################################################
###############################################################################

# Where does the spheral pip_cache dir live.
SPHERAL_PIP_CACHE_DIR=${SPHERAL_PIP_CACHE_DIR:-~/.cache/spheral_pip}

# What is the local script path.
SCRIPT_DIR=${SCRIPT_DIR:-'scripts'}

# DEV_PKG_NAME defines the title of the spheral install. This is a combination of
# the system type and the current spheral version string.
DEV_PKG_NAME=${DEV_PKG_NAME:-$SYS_TYPE-spheral-dev-pkg-undefined}

# CI_BUILD_DIR Place to put the build cache tar file
CI_BUILD_DIR=${CI_BUILD_DIR:-$PWD/../}

# TPL_DIR Where to create the TPL directory
TPL_DIR=${TPL_DIR:-$PWD/../spheral-spack-tpls}

###############################################################################
###############################################################################

# Get name of current directory, should be DEV_PKG_NAME
PKG_DIR=${PWD##*/}

# RESOURCE_DIR is a directory created internally to maintain spack & pip
# resources required for building and running Spheral
RESOURCE_DIR=$PWD/resources

# Print for sanity check.
echo $PWD
echo $RESOURCE_DIR
echo $SPHERAL_PIP_CACHE_DIR
echo $TPL_DIR

SPACK_CACHE_DIR=spack-dev-pkg
echo $SPACK_CACHE_DIR

# Clear the stage directory, create resource dir and copy the Spheral repo into
# the current directory
mkdir -p $RESOURCE_DIR

# Copy the SPHERAL_PIP_CACHE_DIR into resource.
mkdir -p $RESOURCE_DIR/pip_cache
cp -a $SPHERAL_PIP_CACHE_DIR/. $RESOURCE_DIR/pip_cache

# tpl-manager --dev-pkg does the following:
# Creates a local Spack repo
# Activates and concretizes the dev_pkg Spheral Spack environment
# Installs the Spheral dependencies for all specs
# --spack-cache-dir causes tpl-manager to create a tar file called
# <tpl_dir>/../<spack-cache-dir>.tar.gz of the --tpl-dir before any installs occur
run_cmd "./$SCRIPT_DIR/devtools/tpl-manager.py --dev-pkg --clean --tpl-dir $TPL_DIR --spack-cache-dir $SPACK_CACHE_DIR"

# Source Spack for the current terminal
source $TPL_DIR/spack/share/spack/setup-env.sh

# Activate our dev spack environment
run_cmd "spack env activate ./scripts/spack/environments/dev_pkg"

# Create a mirror of all tpl specs in our environment
# (should only be our deps for a single spec in the env).
run_cmd "spack mirror create -a -d $RESOURCE_DIR/mirror --exclude-specs 'llnlspheral spheral'"

# Use spack to list all specs in the mirror and push them to the buildcache.
run_cmd "spack buildcache push -u -f $RESOURCE_DIR/mirror $(spack find --format /{hash})"

# Mirror bootstrap packages needed to start a spack instance on an airgapped system.
run_cmd "spack bootstrap mirror --binary-packages $RESOURCE_DIR"

# Tar up everything in the $CI_BUILD_DIR
run_cmd "tar -czf $CI_BUILD_DIR/$DEV_PKG_NAME.tar.gz -C ../ $PKG_DIR"
