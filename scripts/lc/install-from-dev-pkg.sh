set -Eeuo pipefail
trap 'echo "# $BASH_COMMAND"' DEBUG

BUILD_ALLOC=${BUILD_ALLOC:-""}
SCRIPT_DIR=${SCRIPT_DIR:-'scripts'}
SPHERAL_PIP_CACHE_DIR=${SPHERAL_PIP_CACHE_DIR:-~/.cache/spheral_pip}

if [[ -z "${DEV_PKG_SPEC}" ]]; then
  echo "DEV_PKG_SPEC var must be set."
  exit 1
fi

if [[ -z "${INSTALL_DIR}" ]]; then
  echo "INSTALL_DIR var must be set."
  exit 1
fi

mkdir -p $INSTALL_DIR

# Check if a tar file of the spack and spack package repos exists
SPACK_CACHE_DIR=${SPACK_CACHE_DIR:-${DEV_PKG_SPEC}-spack.tar.gz}
TPL_DIR=${TPL_DIR:-spheral-spack-tpls}
if [[ -f "$SPACK_CACHE_DIR" ]]; then
    # Untar into the install directory
    TPL_DIR=$INSTALL_DIR/$TPL_DIR
    if [[ ! -d "${TPL_DIR}" ]]; then
        tar -xzf $SPACK_CACHE_DIR -C $INSTALL_DIR
    fi
fi
# Now ensure the proper directories exist
if [[ ! -d "${TPL_DIR}/spack" ]]; then
    echo "Spack repo must exist in TPL_DIR"
    exit 1
elif [[ ! -d "${TPL_DIR}/packages" ]]; then
    echo "Spack package repo must exist in TPL_DIR"
    exit 1
fi

echo $DEV_PKG_SPEC
echo $INSTALL_DIR
echo $SCRIPT_DIR
echo $BUILD_ALLOC
echo $PWD

cp -a $PWD/resources/pip_cache/. $SPHERAL_PIP_CACHE_DIR

source $TPL_DIR/spack/share/spack/setup-env.sh
spack env activate ./scripts/spack/environments/dev_pkg
spack mirror remove spheral-mirror || true
spack mirror remove spheral-cache || true
spack bootstrap remove spheral-sources || true
spack bootstrap remove spheral-binaries || true
spack bootstrap add --trust spheral-sources $PWD/resources/metadata/sources
spack bootstrap add --trust spheral-binaries $PWD/resources/metadata/binaries
spack mirror add --unsigned spheral-mirror $PWD/resources/mirror
spack mirror add --unsigned spheral-cache $PWD/resources
spack buildcache update-index $PWD/resources/mirror

# With these inputs, tpl-manager will build with --use-buildcache package:never,dependencies:only -u initconfig
# This ensures the TPLs are only built from cache and Spheral isn't built
$BUILD_ALLOC ./$SCRIPT_DIR/devtools/tpl-manager.py --tpl-dir $TPL_DIR --spec $DEV_PKG_SPEC --skip-init --dev-pkg

HOST_CONFIG_FILE=$(ls -t | grep -E "*\.cmake" | head -1)
$BUILD_ALLOC ./$SCRIPT_DIR/devtools/host-config-build.py --host-config $HOST_CONFIG_FILE -i $INSTALL_DIR --build --no-clean -DSPHERAL_PIP_CACHE_DIR=$SPHERAL_PIP_CACHE_DIR -DSPHERAL_NETWORK_CONNECTED=Off
