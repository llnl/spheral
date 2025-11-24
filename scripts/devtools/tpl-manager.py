#!/usr/bin/env python3

import argparse, os, sys, re, shutil, glob, time

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from spheralutils import sexe

def get_config_dir(base_dir):
    "Return directory containing the repo.yaml file for a base dir"
    return os.path.join(base_dir, "scripts/spack")

# Spack instance info
default_spack_dir = os.path.join(os.getcwd(), "../spheral-spack-tpls")

# Stay away from ~/.spack!
os.environ["SPACK_USER_CACHE_PATH"] = os.path.join(default_spack_dir, "misc")

default_spack_url = "https://github.com/spack/spack.git"
# Spack version: v1.0.2
# spack_commit = "734c5db2121b01c373eed6538e452f18887e9e44"
# Spack version: v1.1.0
spack_commit = "0c2be44e4ece21eb091ad5de4c97716b7c6d4c87"
# Current repo (either LLNLSpheral or Spheral)
package_name = "spheral"

base_dir = os.getcwd()
package_dirs = {"spheral": base_dir}

default_install_args = dict(stop_at="initconfig", fail_fast=True)
chmod_run = "chmod -Rf g+rwX"

# Find if this repo is LLNLSpheral by checking the submodule list
git_mod_cmd = "git config --file .gitmodules --name-only --get-regexp path$"
git_mod_out = sexe(git_mod_cmd, ret_output=True, echo=False)
if "spheral" in git_mod_out:
    package_name = "llnlspheral"
    package_dirs["spheral"] = os.path.join(base_dir, "spheral")
    package_dirs.update({"llnlspheral": base_dir})
print(f"Managing {package_name} TPLs")

class SpheralTPL:
    def parse_args(self):
        parser = argparse.ArgumentParser()
        group = parser.add_mutually_exclusive_group()
        group.add_argument("--spec", type=str, default=None,
                           help="Install TPLs and create host config file for a given spec.\n"+\
                           "If no spec if given, TPLs for all env specs are installed and "+\
                           "no host config file is created.")
        parser.add_argument("--show-specs", action="store_true",
                            help="Show the specs for the current environment and stop.")
        parser.add_argument("--show-info", action="store_true",
                            help="Show the Spheral Spack info, including dependencies and variants.")
        parser.add_argument("--add-spec", action="store_true",
                            help="Set this flag to add the --spec to the environment.")
        parser.add_argument("--spack-dir", type=str,
                            default=default_spack_dir,
                            help="Directory to install Spack instance and a build directory.")
        parser.add_argument("--update-upstream", action="store_true",
                            help="Install TPLs into the upstream instead of the local build. "+\
                            "Installs all specs in the current environment.")
        parser.add_argument("--spack-url", type=str, default=default_spack_url,
                            help="URL to download spack.")
        parser.add_argument("--clean", action="store_true",
                            help="Set this flag to ensure concretizations/installs are fresh. "+\
                            "If issues arise, try using this flag.")
        parser.add_argument("--init-only", action="store_true",
                            help="Download Spack but do not concretize or install.")
        parser.add_argument("--skip-init", action="store_true",
                            help="Skip downloading Spack repo.")
        parser.add_argument("--ci-run", action="store_true",
                            help="For use by the CI only. Must set a --spec.")
        parser.add_argument("--dry-run", action="store_true",
                            help="Use to do everything but actually install. For testing purposes.")
        parser.add_argument("--id", type=str, default=None,
                            help="ID string to postfix an initconfig file.")
        parser.add_argument("--package-repo", type=str, default=None,
                            help="Specify a location to put the spack-package repo. "+\
                            "Defaults to spheral-spack-tpls/packages.")
        parser.add_argument("--dev-pkg", action="store_true",
                            help="Tells tpl-manager to use the dev_pkg environment. "+\
                            "Assumes TPLs are for buildcache creation if no --spec is provided. "+\
                            "Assumes building from a buildcache if --spec is provided.")

        self.args = parser.parse_args()

        if (not self.args.spec and self.args.ci_run):
            raise Exception("Must specify a --spec if doing --ci-run")

        if (self.args.spec and not self.args.spec.startswith(package_name)):
            raise Exception(f"--spec must start with {package_name}")
        if (self.args.spec):
            print(f"Installing {self.args.spec}")

    def add_spack_paths(self, spack_dir, package_repo):
        "Append spack path to system to use spack python modules"
        spack_path = os.path.join(spack_dir, "lib", "spack")
        sys.path.append(spack_path)
        sys.path.append(os.path.join(spack_path, "spack"))
        sys.path.append(os.path.join(spack_path, "_vendoring"))
        global spack, SpackCommand
        try:
            import spack
            from spack.main import SpackCommand
            spack = spack

            #
            # Workaround: Cause spack to set its internal working directory
            # state before we make any further spack API calls.
            spack.paths.set_working_dir()
        except ImportError as e:
            raise ImportError("Failed to import Spack python module") from e
        repo_yaml = os.path.join(spack_dir, "etc/spack/repos.yaml")
        if (not os.path.exists(repo_yaml)):
            print("Spack instance does not have a builtin package repo set.")
            if package_repo:
                builtin_package = package_repo
                if not os.path.isabs(package_repo):
                    builtin_package = os.path.join(spack_dir, package_repo)
            else:
                builtin_package = os.path.join(spack_dir, "../packages")
            print(f"Creating a builtin package repo at {builtin_package}")
            from spack.main import SpackCommand
            repo_cmd = SpackCommand("repo")
            repo_args = ["set", "--destination", os.path.abspath(builtin_package), "builtin"]
            repo_cmd(*repo_args)
            repo_cmd(*["update"])

    def print_specs(self, specs):
        if (type(specs) != list):
            specs = [specs]
        install_status = spack.spec.Spec.install_status
        print(spack.spec.tree(specs, format=spack.spec.DISPLAY_FORMAT, status_fn=install_status, hashes=True, hashlen=6))

    def clone_spack(self):
        "Clone Spack and add paths to use spack python"
        tpl_root = self.args.spack_dir
        if (not os.path.exists(tpl_root)):
            os.mkdir(tpl_root)
        spack_dir = os.path.join(tpl_root, "spack")
        if (not self.args.skip_init):
            if (not os.path.exists(spack_dir)):
                sexe(f"git -C {tpl_root} clone --depth=2 {self.args.spack_url}")
            # Check commit hash of Spack repo
            cur_hash = sexe(f"git -C {spack_dir} rev-parse HEAD", ret_output=True, echo=False).strip()
            if (cur_hash != spack_commit):
                sexe(f"git -C {spack_dir} fetch --depth=2 origin {spack_commit}")
                sexe(f"git -C {spack_dir} checkout FETCH_HEAD")
        self.add_spack_paths(spack_dir, self.args.package_repo)

    def check_lock_file(self):
        "Check if any files in scripts/spack are newer than the spack.lock file"
        spack_lock = os.path.join(self.env_dir, "spack.lock")
        if (not os.path.exists(spack_lock)):
            return False
        lock_time = os.path.getmtime(spack_lock)
        for package, path in package_dirs.items():
            internal_spack_dir = os.path.join(os.path.abspath(get_config_dir(path)), "**")
            sp_files = glob.glob(internal_spack_dir, recursive=True)
            ctimes = [os.path.getmtime(ff) for ff in sp_files\
                      if os.path.isfile(ff) and "__pycache__" not in ff]
            if (any(x > lock_time for x in ctimes)):
                return True
        return False

    def find_spack_package(self, package_name, req=True):
        """
        Find if the package name/s is already in the environment.
        If any one of them are found, return the name of that package.
        If not, run spack external find and raise an error if none are found.
        """
        cur_packages = spack.config.get("packages")
        if (package_name in cur_packages):
            return True
        ext_cmd = SpackCommand("external")
        if (req):
            find_out = ext_cmd("find", "--not-buildable", package_name)
        else:
            find_out = ext_cmd("find", package_name)
        if (package_name in find_out):
            return True
        return None

    def find_all_spack_packages(self, packages, req=True):
        for i in packages:
            if (not self.find_spack_package(i, req) and req):
                raise Exception(f"System install of {i} not found. "+\
                                "If software is installed, add location to $PATH "+\
                                "environment variable. Otherwise, install package.")

    def modify_env_file(self, env_file, mod_func):
        "Modify the spack.yaml file"
        from spack.util import spack_yaml
        # Load the spack.yaml file
        with open(env_file) as ff:
            try:
                loader = spack_yaml.load(ff)
            except SpackYAMLError as exception:
                print(exception)
        loader = mod_func(loader)
        with open(env_file, 'w') as ff:
            spack_yaml.dump(loader, ff)

    def custom_spack_env(self, env_name):
        "Use/create a custom Spack environment"
        from spack import environment
        if (not self.args.spec):
            raise Exception("Must supply a --spec for a custom environment (IE --spec spheral+mpi%gcc)")
        env_file = os.path.join(self.env_dir, "spack.yaml")
        if (not os.path.exists(env_file)):
            # Create a new environment
            env_cmd = SpackCommand("env")
            env_cmd("create", "--without-view", "-d", self.env_dir)
            def set_concretize(loader):
                new_dict = {"concretizer": {"unify": False}}
                loader["spack"].update(new_dict)
                return loader
            self.modify_env_file(env_file, set_concretize)
        self.spack_env = environment.Environment(self.env_dir)
        environment.activate(self.spack_env)
        # Get all the Spack commands
        repo_cmd = SpackCommand("repo")
        dev_cmd = SpackCommand("develop")
        comp_cmd = SpackCommand("compiler")
        ext_cmd = SpackCommand("external")
        cur_repos = repo_cmd("list") # spack repo list

        # Add the repos and develop paths to the spack environment
        for package, path in package_dirs.items():
            if (package+" " not in cur_repos):
                repo_path = os.path.abspath(os.path.join(get_config_dir(path), f"spack_repo/{package}"))
                repo_cmd("add", f"{repo_path}") # spack repo add <repo_path>
            dev_path = os.path.abspath(path)
            dev_cmd("-p", dev_path, f"{package}@=develop") # spack develop <package>@=develop
        comp_cmd("find") # spack compiler find
        ext_cmd("find") # spack external find
        provider_dict = {}

        # req_packages are packages we refuse to let spack build
        # If they aren't found on the system, an error will be thrown
        # telling the user to install the package or ensure an install is in $PATH
        self.find_all_spack_packages(["python", "perl"])
        # Look for MPI installs
        if (not spack.spec.Spec(self.args.spec).satisfies("~mpi")):
            found_mpi = False
            for i in ["openmpi", "mpich"]:
                if (self.find_spack_package(i)):
                    found_mpi = True
                    provider_dict.update({"mpi": [i]})
                    break
            if (not found_mpi):
                raise Exception(f"System MPI install not found. "+\
                                "If software is installed, add location to $PATH "+\
                                "environment variable. Otherwise, install package.")
        # opt_packages are packages we would like to find on the system
        # but are ok with Spack installing if they are not found
        opt_packages = ["hdf5", "ncurses"]
        for i in opt_packages:
            self.find_spack_package(i, req=False)
        # Hard-coding providers for these packages vastly improves TPL system
        provider_dict.update({"zlib-api": ["zlib"],
                              "blas": ["openblas"],
                              "lapack": ["openblas"]})
        # Always add the spec for a custom environment
        if (provider_dict):
            self.config_env_providers(provider_dict)
        self.args.add_spec = True

    def config_env_providers(self, config_dict):
        env_file = os.path.join(self.env_dir, "spack.yaml")
        def set_providers(loader):
            new_dict = {"all": {"providers": config_dict}}
            loader["spack"]["packages"].update(new_dict)
            return loader
        self.modify_env_file(env_file, set_providers)

    def activate_spack_env(self):
        "Activates a Spack environment or creates and activates one when necessary"
        config_env_dir = os.path.join(get_config_dir(base_dir), "environments")
        # Check if we are on an LC machine and the environment exists
        default_env = os.getenv("SYS_TYPE")
        if (self.args.dev_pkg):
            default_env = "dev_pkg"
        if default_env and os.path.exists(os.path.join(config_env_dir, default_env)):
            # For LC systems
            self.env_dir = os.path.join(config_env_dir, default_env)
            print(f"Activating Spack environment in {self.env_dir}")
            from spack import environment
            self.spack_env = environment.Environment(self.env_dir)
            environment.activate(self.spack_env)
        else:
            # Otherwise, check if environment has been created
            arch_cmd = SpackCommand("arch")
            env_name = arch_cmd().strip()
            self.env_dir = os.path.join(config_env_dir, env_name)
            self.custom_spack_env(env_name)

    def concretize_spec(self, check_spec):
        "Concretize the spec"
        if (self.args.add_spec):
            self.spack_env.add(self.args.spec)
        force_conc = False
        if (self.args.clean):
            force_conc = True
            print("Cleaning and concretizing environment")
        else:
            print("Concretizing environment")
        conc_spec = self.spack_env.concretize(force=force_conc)
        if conc_spec:
            print("Concretized specs")
            for x in conc_spec:
                print(x)
        if (check_spec):
            self.spack_spec = spack.spec.Spec(self.args.spec)
            matches = self.spack_env.matching_spec(self.spack_spec)
            if (not matches):
                raise Exception(f"{self.args.spec} not found in current "+\
                                "environment. Rerun with --add-spec to add it.")
            self.spack_spec = matches
            print(f"Found matching root spec for {self.args.spec}")
            self.print_specs(self.spack_spec)
        else:
            specs = self.spack_env.concrete_roots()
            self.print_specs(specs)
        self.spack_env.write()

    def install_and_config(self):
        "Install TPLs and create host config file for given spec"
        spec = self.args.spec
        # Load the spack package recipe python class
        if (package_name == "llnlspheral"):
            from spack.pkg.llnlspheral.llnlspheral import Llnlspheral
            spack_spheral = Llnlspheral(self.spack_spec)
        else:
            from spack.pkg.spheral.spheral import Spheral
            spack_spheral = Spheral(self.spack_spec)

        # Get host config file name from spack package recipe
        host_config_file = spack_spheral.cache_name
        # If using --id, preserve original host config file
        mod_host_config = False
        if (self.args.id and os.path.exists(host_config_file)):
            # Avoid overwriting existing host config file
            shutil.copyfile(host_config_file, "orig"+host_config_file)
            mod_host_config = True
        if (self.args.dev_pkg):
            self.spack_env.install_specs([self.spack_spec],
                                         package_use_cache=False,
                                         dependencies_cache_only=True,
                                         unsigned=True,
                                         **default_install_args)
        else:
            self.spack_env.install_specs([self.spack_spec], **default_install_args)
        if (self.args.ci_run):
            shutil.copyfile(host_config_file, "gitlab.cmake")
            host_config_file = "gitlab.cmake"
        if (mod_host_config):
            # Apply --id and bring back original host config file
            new_name = host_config_file.replace(".cmake", f"{self.args.id}.cmake")
            os.rename(host_config_file, new_name)
            os.rename("orig"+host_config_file, host_config_file)
            host_config_file = new_name
        print(f"Created {host_config_file}. To configure, run:\n")
        scr_path = os.path.relpath(package_dirs['spheral'], base_dir)
        print(f"{scr_path}/scripts/devtools/host-config-build.py --host-config {host_config_file}\n")

    def __init__(self):
        self.parse_args()
        self.clone_spack()
        if (self.args.init_only):
            return
        self.activate_spack_env()
        if (self.args.show_specs):
            find_cmd = SpackCommand("find")
            find_cmd("-r")
            sys.exit(0)
        if (self.args.show_info):
            info_cmd = SpackCommand("info")
            info_cmd("spheral")
            sys.exit(0)
        # Check if any files in scripts/spack are newer than the spack.lock
        if (not self.args.clean):
            self.args.clean = self.check_lock_file()
        if (self.args.spec):
            # If --spec is given, install TPLs and create host config file
            self.concretize_spec(check_spec=True)
            self.install_and_config()
        else:
            # Concretize the current environment
            self.concretize_spec(check_spec=False)
            if self.args.update_upstream:
                upstream_dir = spack.config.get("upstreams:spheral_shared:install_tree")
                with spack.config.override("config:install_tree", upstream_dir):
                    spack.config.set("config", {"install_tree": {"padded_length": 0}})
                    print("WARNING: Modifying local Spack files, do not commit these changes")
                    with self.spack_env.manifest.use_config():
                        print(spack.config.get("config:install_tree"))
                        print(f"Installing to {upstream_dir}")
                        # Pass None so it installs TPLs for all specs
                        self.spack_env.install_all(install_deps=True, install_package=False, fail_fast=True)
                        # Equivalent of running spack reindex
                        spack.store.STORE.reindex()
                chmod_cmd = chmod_run + f" {upstream_dir}"
                os.system(chmod_cmd)
            else:
                self.spack_env.install_specs(None, **default_install_args)

        # Remove symbolic directory created by Spack
        print("Removing Spack symbolic build directories")
        build_dirs = glob.glob("build-*")
        for i in build_dirs:
            if (os.path.islink(i)):
                os.unlink(i)

if __name__=="__main__":
    spheral_tpl = SpheralTPL()
