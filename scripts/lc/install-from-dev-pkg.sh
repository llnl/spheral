set -Eeuo pipefail
trap 'echo "# $BASH_COMMAND"' DEBUG

SPACK_URL=${SPACK_URL:-'https://github.com/spack/spack'}
BUILD_ALLOC=${BUILD_ALLOC}
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

echo $DEV_PKG_SPEC
echo $SPACK_URL
echo $INSTALL_DIR
echo $SCRIPT_DIR
echo $BUILD_ALLOC
echo $PWD

mkdir -p $INSTALL_DIR
cp -a $PWD/resources/pip_cache/. $SPHERAL_PIP_CACHE_DIR

./$SCRIPT_DIR/devtools/tpl-manager.py --spack-url $SPACK_URL --init-only --no-upstream --spack-dir $INSTALL_DIR/spheral-spack-tpls

source $INSTALL_DIR/spheral-spack-tpls/spack/share/spack/setup-env.sh
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
$BUILD_ALLOC ./$SCRIPT_DIR/devtools/tpl-manager.py --no-upstream --spack-dir $INSTALL_DIR/spheral-spack-tpls --spec $DEV_PKG_SPEC --skip-init --dev-pkg

HOST_CONFIG_FILE=$(ls -t | grep -E "*\.cmake" | head -1)
$BUILD_ALLOC ./$SCRIPT_DIR/devtools/host-config-build.py --host-config $HOST_CONFIG_FILE -i $INSTALL_DIR --build --no-clean -DSPHERAL_PIP_CACHE_DIR=$SPHERAL_PIP_CACHE_DIR -DSPHERAL_NETWORK_CONNECTED=Off
