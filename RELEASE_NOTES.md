Version v2025.06.0 -- Release date 2025-06-18
==============================================
  * Important Notes:

Notable changes include:

  * New features / API changes:
    * The tpl-manager.py is completely overhauled to include the following:
      * Utilize the Spheral Spack environments.
      * Handle some build cache functionality.
      * Do things Uberenv did like download and install Spack itself.
    * Implicit time integration is now supported.
      * CrankNicolsonIntegrator is the current implicit default
      * The Physics package interface has been augmented to support implicit integration with two important
        methods that must be provide:
        * Physics::dtImplicit to provide a maximum bounding time step
        * Physics::maxResidual should provide a maximum dimensionles residual change to check for convergence
      * Sundials is now an optional (but default) third-party lib in Spheral, and provides a non-linear solver
        we now wrap in a new Solver interface in Spheral (wrapping Sundials KINSOL solver).
    * FSISPH has a new flag (decoupleDamagedMaterial, default True) which can be turned off to more tightly
      couple damaged to undamaged material.
      * planeStrain has been removed as an option in FSISPH as part of unifying deviatoric evolution with other
        hydros.
    * LEOS (Livermore Equation Of State) package now available in Spheral.  Requires access to the LEOS
      package itself, which most folks outside LLNL will not have.
    * gtest suite integration for writing minimal c++ unit tests on the host and device without needing to
      compile large potions of the code.
    * Both ASPH and ASPHClassic now allow the user to override the final H evolution through optional functors added to the classes:
      - HidealFilter
      - RadialFunctor
    * FacetedSurfaceASPHHydro has been removed in favor of providing user filters to the ASPH methods (i.e., the RadialFunctor method).
    * Field resizing operations have been removed from the public interface.
    * Performance analysis tools are improved.
      - Added an "advance" Caliper timer to be used in the future as the default reference timer.
      - Added a deploy CI stage to create a GitLab page with the historical performance benchmarks.

  * Build changes / improvements:
    * Native Spack environments are now being used.
      * Uberenv is no longer used.
      * Adds logic to simplify building on non-LC systems; tries to find existing installed compilers and packages.
      * Adds spack.yaml environment files for current LC systems and a dev_pkg environment, which is used for creating the build cache.
      * Local Spack packages for TPLs are removed or simplified when possible since the builtin Spack packages are no longer replaced.
      * The upstream Spack instance is no longer used when creating the build cache.
      * The package.yaml for Spheral is improved to allow full Spheral installation through Spack.
      * Centralizes things like upstream location, compiler types and versions, and specs in the environments and configs.

  * Bug Fixes / improvements:
    * ATS submodule is updated to fix bug with latest Flux update on LC systems.
    * Update Polytope version.
    * TPL manager removes the symoblic links to the install directory.
    * Consolidated CMake configured files into SpheralConfigs.py.in.
    * Deviatoric stress evolution in lower dimensions (1 and 2D) now consistent with other solid hydros.
    * Changes the `SPHERAL_TEST_INSTALL_PREFIX` to be relative to `CMAKE_INSTALL_PREFIX/tests` directory.
    * Fixed bug where performance tests would incorrectly move a benchmark directory if rerunning failed jobs.
    * Cleaned out some old unused code.
    * Fixed bugs related to DEM timestep, redistribution, and thread safety.

Version v2025.01.0 -- Release date 2025-01-31
==============================================
  * Important Notes:

Notable changes include:

  * New features / API changes:
    * MPI variables are now wrapped as
      ```
      SPHERAL_OP_SUM, SPHERAL_OP_MAX, SPHERAL_OP_MIN
      ```
    * CHAI added as a submodule of Spheral for co-developing features necessary for GPU port.
    * RAJA & Umpire added as first level dependencies.
    * Axom updated to v0.9.0.
    * TPL builds have been split off into a separate Gitlab CI stage to help with timeouts on allocations.
    * Failed ATS runs are automatically retested once in the Gitlab CI.
    * Python execute command is centralized in scripts/spheralutils.py now.
    * Caliper updated v2.11.
    * Adiak added as TPL.
    * Created singleton wrapper for cali::ConfigManger and python wrapped Caliper timer and Adiak routines.
    * New ASPH idealH algorithm implemented, which is much more robust and accurate as H elongations become extreme.
    * New experimental hourglass control algorithm implemented, along with some basic tests/demonstrations.
    * H update algorithms converted to their own independent physics packages, no longer part of the various hydro packages.
    * Physics interface updated slightly:
      * Physics::postStateUpdate now returns a bool indicating if boundary conditions should be enforced again.
      * Physics packages can now have Physics sub-packages, which can be run before or after the main package.  The SpheralController
        now checks for these packages and adds them to the physics package list as needed.
      * Physics packages can indicate if they require Voronoi cell information be available. If so, a new package which computes and
        updates the Voronoi information is automatically added to the package list by the SpheralController (similar to how the
        Reproducing Kernel corrections are handled).
    * Command line options are now consistent. Default values of a string "None" are no longer allowed and any input through the command line of "None" will become the python NoneType None.
    * Cleaned up use of std::any in State objects using a visitor pattern to be rigorous ensuring all state entries are handled properly
      during assignement, equality, and cloning operations. This is intended to help ensure our Physics advance during time integration
      is correct.
    * Performance regression testing is now available. All developers are encouraged to run the performance testing suite for any code changes that might impact performance. See documentation for more details.
    * Added our old ASPH IdealH H update as an option. While it is not as reliable as our current default ASPH, it does not require building the Voronoi and is therefore signifcantly faster.
    * Converted artificial viscosities to Physics packages, and add them as pre-subpackages to Hydro objects.
    * Split artificial viscosities based on the type of pressure they compute (currently Scalar or Tensor), which is slightly more efficient.
      * This required making the hydro packages evaluateDerivatives into templated methods based on the type of Q they are handed.
      * Also introduced a new base class (ArtificialViscosityHandle), which provides a handle class not templated on the type
        of Q pressure for Hydro objects to hold onto.

  * Build changes / improvements:
    * Distributed source directory must always be built now.
    * Git strategies in the Gitlab CI are fixed so a clone only occurs on the first stage for each job, instead of for all stages for each job.
    * New Gitlab CI pipeline cleanup strategy deletes job directories immediately upon successful completion.
    * The FSISPH package is now optional (SPHERAL\_ENABLE\_FSISPH).
    * The GSPH package is now optional (SPHERAL\_ENABLE\_GSPH).
    * The SVPH package is now optional (SPHERAL\_ENABLE\_SVPH).
    * Cleaner Spheral Spack package.
    * ENABLE\_DEV\_BUILD can now export targets properly.
    * Added a GCC flag to prevent building variable tracking symbols when building PYB11 modules.  This is unnecessary, and
      on some platforms trying to build such symbols is very expensive and in some cases fails.
    * Consolidates lcatstest.in and run\_ats.py into a single spheral\_ats.py script.
    * SPHERAL\_TEST\_INSTALL\_PREFIX now includes the tests directory.
    * Removed most configured files and added a SpheralConfigs.py file to use at runtime instead.
    * Python runtime packages are now handled in the Spheral build pipeline with pip.
      * Removed pip package dependencies from spack.
      * Introduced Spheral_Python_Env function to manage Python environments for build and runtime dependencies.
      * spheral-setup-venv now only copies installed Spheral libraries to environments at install time.
      * Added pip cache support to local directory (~/.cache/spheral_pip/), customizable via SPHERAL_PIP_CACHE_DIR.
      * Added ATS as a submodule due to lack of PyPI package.

    * Moved Spheral from BlueOS/NVIDIA systems to support CRAY/AMD.
      * Migrated CI to CRAY/AMD due to pip compatibility issues with BlueOS.
      * Added HIP support for device/offload tests and updated TPLs for HIP-enabled builds.
      * Updated GitLab CI and Developer scripts for flux scheduling system compatibility.

  * Bug Fixes / improvements:
    * Wrappers for MPI calls are simplified and improved.
    * Time step estimate due to velocity divergence in RZ space has been fixed.
    * Fixed tolerances for ANEOS equation of state temperature lookup
    * Clang C++ warnings have eliminated, so the Clang CI tests have been updated to treat warnings as errors.
    * Fix for installing libraries when building individual package with ENABLE\_DEV\_BUILD=On.
    * Bugfix for RZ solid CRKSPH with compatible energy.
    * Parsing of None string now always becomes None python type. Tests have been updated accordingly.
    * IO for checkpoints and visuzalization can now be properly turned off through SpheralController input options.
    * Bugfix for atomicWeight in ANEOS.
    * Fixed porosity model interaction with damage for zero porosity case.

Version v2024.06.1 -- Release date 2024-07-09
==============================================

  * Important Notes:
    * This is a patch release for v2024.06.0.

  * Bug Fixes / improvements:
    * CD pipeline hotfix for installing release builds on LC machines.
    * Fixes an issue with the use of the axom::quest::SignedDistance interface. 

Version v2024.06.0 -- Release date 2024-06-27
==============================================
  * Important Notes:
    * External users of the code will need to supply config files for tpl-manager to find system libraries correctly. Steps to do this are detailed in the external user build guide. 

Notable changes include:

  * New features / API changes:
    * Added MFV hydro from Hopkins 2015 with extension for ALE options.
    * Adding optional user specified smoothing scale method for SPH, FSISPH, and CRKSPH.

  * Build changes / improvements:
    * PYBind11 libraries no longer depend on the structure of the PYB11 source directory.
      * CMake interface for adding PYBind11 target libraries is modified to more closely match how C++ libraries are created.
      * Multiple Spheral Python modules / CMake targets can be specified for a single directory.
      * KernelIntegrator and FieldList directories are divided into 2 modules / targets.
    * tpl-manager.py will no longer use generic x86_64 configs for non LC systems. Users will be required to supply their own configs for pointing spack at external packages.
    * Spack version is increased from 0.19 to 0.22.
    * Spack upstream is updated.
    * Removed the python 3 module load for the Gitlab CI to fix an issue with pkg-config changing.
    * Zlib target and TPL cmake file is removed.
    * PYB11Generator repo is updated.
    * Spack config and package files inside Spheral are updated to accommodate Spack 0.22.
      * Package recipes for py-numpy-stl, py-pillow, py-pipreqs, td, and tk are removed.
      * Versions for python dependencies in the Spheral spack recipe are fixed and updated (in some cases).

  * Bug Fixes / improvements:
    * Corrected an erroneous VERIFY in the P-alpha porosity constructor (with Fields of porosity and sound speed) that forced runs to stop even with correct input parameters
    * Fixed a bug in the standard ASPH hydros (ASPH, SolidASPH, and RZ varieties) that gave incorrect results.  FSI ad CRK models with ASPH smoothing scales were OK, but standard
      SPH using ASPH smoothing scales were simply incorrect for non-unit aspect ratio H's.  Also added ATS tests to help catch such errors going forward.

Version v2024.01.1 -- Release date 2024-02-17
==============================================
  * Important Notes:
    * This is a patch release for v2024.01.0.

Notable changes include:

  * New features/ API changes:
    * Adding an optional second-stage problem start-up hook to the Physics package interface: Physics::initializeProblemStartupDependencies.  The idea is to keep basic sizing
      of arrays and such in the first stage (Physics::initializeProblemStartup), while this new hook is used for updating any initial Physics state (and therefore provides a
      State and StateDerivatives object).
    * DEM
      * new field list to track max particle overlap
      * user can optional turn off fast time stepping
      
  * Build changes / improvements:
    * Improved the target export functionality.

  * Bug Fixes / improvements:
    * Fixed bug with ConstantBoundary in the presence of porosity with the new porosity models introduced in v2024.01.00.
    * Updating header lists for including Spheral modules in external projects.
    * Adding effective viscous pressure back to FSISPH.
    * Initial volumes for damage models were incorrectly not taking into account pore space when computing failure statistics for seeding flaws.  Fixed.
    * DEM
      * fixed bug in solid boundary unique indices that causes particle sticking
      * fixed bug in solid boundary update policies 
      * fixed solid boundary restartability for moving bcs

Version v2024.01.00 -- Release date 2024-01-19
==============================================
  * Important Notes:
    * The PolyClipper, BLT, and PYB11Generator submodules have been modified. Be sure to recursively update the submodules.  

Notable changes include:

  * New features/ API changes:
    * Adding P-alpha porosity model.
    * Updating treatment of various state variables in the presence of porosity.
    * Introduced a new common base class for porosity physics (PorosityModel), which PalphaPorosity and StrainPorosity share.
    * Revamped interaction UpdatePolicies with FieldLists:
      - UpdatePolicies have a new virtual method: clonePerField: True means when registering a FieldList copy the Policy for each Field in the FieldList; False means register the FieldList for update itself with the single instance of the Policy.
      - This change removes most of our redundant Field/FieldList update policies, and allows us to be more granular in applying different policies to single Field values in a FieldList.
    * Adding more Shadow Python interfaces wrapping our C++ classes, in particular PalphaPorosity and StrainPorosity.
    * EquationOfState now requires instances to provide \partial P/\partial \rho and \partial P/\partial \epsilon.  All current equations of state have been updated accordingly.
    * Tillotson and Gruneisen EOSs implementations updated a bit in the revamping.
    * Added more material options to MaterialPropertiesLib.py (mostly from Melosh's 89 book).

  * Build changes / improvements:
    * Spheral now provides First Class CMake support (using the BLT nomenclature). Spheral and its dependencies are now exported to simplify importing the project. To import Spheral into another project using CMake, use:
      ```
      find_package(Spheral_CXX <path_to_spheral_installation>)
      ```
    * CMake variables have a more consistent naming convention. Unused variables are removed.
    * Added ENABLE_DEV_BUILD option to improve build times during code development.
    * Upped our required C++ standard to 17.

  * Bug Fixes / improvements:
    * Fixed melt behavior in Steinberg-Guinan strength model, which was ignoring melt for damaged material.
    * Fixed range of dimensionless melt temperature for Johnson-Cook strength.
    * FSISPH new features and modifications to method. 
      * NOTE constructor inputs have changed.
      * strength implementation modified.
      * new features added including plane strain option and settable minP for interfaces.
      * new, more rigorous, interface and free surface tracking.
    * Fixed initialization of longitudinal sound speed and Youngs modulus for damage models.
    * Corrected some minor bugs/inconsistencies in the Tillotson EOS.
    * lcats updated to work with current TOSS4 machine configurations.
    * Updated various tests to make out automated testing more robust.

Version v2023-06-0 -- Release date 2023-06-20
==============================================

  * Important Notes:
    * This is the first release using Python 3.
    * Restart files are not backwards compatible with prior releases.

Notable changes include:

  * New features/ API changes:
    * We now use Python 3 rather than Python 2. The conversion from Python 2->3 is larger than most prior Python updates, which is part of why we held off for so long on this conversion. Python 2 is no longer supported however, and most other Python packages are only Python 3 now as well. Python comes with a tool (2to3) which can do most of the work converting a Python 2 script to Python 3, so we recommend when updating one of your old Python 2 run scripts or utilities you first let run 2to3 over that script.
    * FileIO has been simplified in how it handles Fields, which has reduced the API of FileIO.  However, this has two compromises:
      * Restart files are not backwards compatible with older releases.
      * Fields in restart files are pretty much binary blobs, and so are not easily inspected by hand.
    * We have removed support for the Python interface to Gnuplot, since the Python Gnuplot package had significant interface changes. Please convert such interactive graphics to Matplotlib, which has very similar support to our older Gnuplot utilities by importing SpheralMatplotlib.

  * Build changes / improvements:
    * PYB11Generator and PolyClipper are now submodules (in extern/)
    * Converted our PYB11 generator build rules to use the newly provided PYB11Generator Cmake rules.
    * toss_4_x86_64_ib system compatibility.
    * Updated spack version to v0.19.1
    * `--debug` and `--no-spec` options added to tpl-manager.py for outputing debug and skipping `spack spec` step.
    * `--init-only` option in tpl-manager.py will only initialize a local spack instance, skipping any TPL configuration.
    * TOSS4 compatibility for LC systems.
    * "risky" builds are installed on LC machines through gitlab CI.

  * Bug Fixes / improvements:
    * Fixed numerous compiler warnings with newer compilers such as G++ 9.4.
    * r-path for additional TPLs can be propogated to Spheral libraries with `SPHERAL_ADDITIONAL_RPATHS`.
    * The DEM package has received a significant updated.
      * Added simple analytic solid boundaries (planes,cylinder,sphere).
      * User given more control over DEM fast time stepping.
      * Simplified initialization and improved robustness.

Version v2023.03.0 -- Release date 2023-03-29
==============================================

This release contains ...

  * Important Notes:
    * Spheral now requires C++ 14
    * Spheral C++ libs are compiled into a single library (Spheral_CXX)

Notable changes include:

  * New features/ API changes:
    * New Discrete Element Model (DEM) physics package with linear-damped spring approach.
    * Adding a CUDA smoke test that can be called from the Spheral python API.
    * NVCC / CUDA 11 gitlab-ci jobs.
    * ATS default filters for non-MPI, debug and CUDA builds are injected into spheral-atstest script.
    * Latest Develop docker containers hosted on ghcr.io/llnl/spheral:latest.
    * External / offline builds are tested through github actions.
    * New polyhedral gravity solver.
    * Improved DEM timestep choice, so that points that do not overlap do not overly constrain the timestep.
    * Point potential gravity solver can now have non-unit metrics.
    * Pair-max damage algorithm now decouples points in different NodeLists or with different fragmentIDs.
    * Added direct support for FacetedVolumes to FileIO interface.
    * ANEOS now allows shallow copies to be made.
    * Polyhedron contain method can now optionally use Axom methods to test for containment.
    * Adding user specified functions for shear modulus and yield strength as a function of damage.

  * Build changes / improvements:
    * The C++ library interface is compiled into a single Spheral_CXX library. 
    * Previous Spheral C++ libraries are still CMAKE targets, but are "ojbect" libraries.
    * ATS bumped to version 7.0.9 for blueos smpi option support.
    * Eigen bumped to 3.4.0 for NVCC compatiblity.
    * C++ flag suppression is gaurded with build time CMake generators to only apply to C++ compilers.
    * Python runtime libraries are now managed through Spack / tpl-manager.
    * Added ENABLE_NAN_EXCEPTIONS (default OFF) Cmake flag to raise an exception when a NAN occurs (Gnu only).
    * Byte-compiling python installed in virtual spheral environment.
    * Invoking spheral no longer byte-compiles Python imported in a spheral script.
    * HDF5 links as a shared library. Static lib usage can be forced with ENABLE_STATIC_TPL.

  * Bug Fixes / improvements:
    * spheral-atstest scripts always point to locally installed ATS instance.
    * gitlab-ci report-results script for analyzing ATS CI runs.
    * Support for offline Spheral builds (provided TPLs are installed).
    * Fixes for restarting without regenerating the original node positions in the Python script.
    * Protected from division by zero in DEM when points coincide.
    * Corrected support for minimum pressure (intact and damaged) with porous materials.
    * Removed term driving damaged material to the reference density in solid hydros.
    * Added verbose flag to EquationOfState::specificThermalEnergyForPressure so users can see how the iterative search proceeds.

Version v2022.06.1 -- Release date 2022-06-24
=============================================

This is a bugfix release, which corrects a path problem that broke our convenient ANEOS
constructors using the provided input for quartz, dunite, and serpentine.


Version v2022.06.0 -- Release date 2022-06-09
=============================================

This release contains new features, bugfixes, and build improvements. Please see the
Spheral Documentation for more information about items in this release.

  * Important Notes: This Spheral release is the first release to support the new Spack based TPL (Third Party Library) management system. The previous CMake based TPL system has been deprecated and removed. The documentation has detailed instructions on how to use this new TPL system.

Notable changes include:

  * New features/ API changes:
    * Spack TPL management system.
    * Support for Gitlab-CI testing with Spheral.
    * Addition of SidreFileIO.
    * Support for Sidre Parallel IO (SPIO).
    * Full Documentation re-organization.
    * GSPH now has it's own pure-virtual GenericHydro class.
    * Adding some persistent state to the hydro objects to remeber criteria for diagnostics.
    * Moved inferace fields of SlideSurface class into the Hydro class.
    * FSISPH handles same-material damaged strength similarly to Spheral's default SPH solver.
    * New ProbabilisticDamageModel, which should be used in place of our prior Grady-Kipp implementations.
    * More Damage application options, and new defaults
    * Artificial viscosity for SPH variants now defaults to LimitedMonaghanGingoldViscosity.

  * Build changes / improvements:
    * Deleting CMake TPL system and all AutoTools Build system artifacts.
    * Removing support for mirrors in spack tpl system.

  * Bug Fixes / improvements:
    * Spheral fixed when running in Debug mode with MPI=Off.
    * Typos fixed in quickstart guide. https://github.com/LLNL/spheral/pull/116
    * Pedantic check for expired pointer to the RestartRegistrar. Ensures we don't call into deleted objects. 
    * Switching GammaLaw and Polytropic EOS to the isentrpic bulk modulus for consistency w/ Solid EOS.
    * Update scalar and tensor damage calc in FSISPH to be more consistent with SolidSPHHydro.
    * CullenDehnen segfault fix.

**Full Changelog**: https://github.com/LLNL/spheral/compare/2022.2.0-pre-spack...v2022.6.0

Version vYYYY.MM.p -- Release date YYYY-MM-DD
==============================================
  * Important Notes:

Notable changes include:

  * New features / API changes:

  * Build changes / improvements:

  * Bug Fixes / improvements:
