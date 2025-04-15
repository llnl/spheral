trap 'echo "# $BASH_COMMAND"' DEBUG

### Create a tar file containing:
# dev-pkg/
#   *cloned spheral repo.*
#   resources/
#     pip_cache/
#     mirror/
#     build_cache/
#     bootstrap/
#       metadata/
#       bootstrap_cache/

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
./$SCRIPT_DIR/devtools/tpl-manager.py --dev-pkg

# Source Spack for the current terminal
source ../spheral-spack-tpls/spack/share/spack/setup-env.sh

# Activate our dev spack environment
spack env activate ./scripts/spack/environments/dev_pkg

# Create a mirror of all tpl specs in our environment
# (should only be our deps for a single spec in the env).
spack mirror create -a -d $RESOURCE_DIR/mirror --exclude-specs "llnlspheral spheral"

# Use spack to list all specs in the mirror and push them to the buildcache.
spack buildcache push -auf $RESOURCE_DIR/mirror $(spack find --format /{hash})

# Mirror bootstrap packages needed to start a spack instance on an airgapped system.
spack bootstrap mirror --binary-packages $RESOURCE_DIR

# Tar up everything in the $CI_BUILD_DIR
tar -czf $CI_BUILD_DIR/$DEV_PKG_NAME.tar.gz -C ../ $PKG_DIR
