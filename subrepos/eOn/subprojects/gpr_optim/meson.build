project('gprd', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++17'])

gprd_deps = []

# Handle Eigen dependencies gracefully
eigen = dependency('Eigen3', method : 'cmake', modules : ['Eigen3::Eigen'], required : false)
if eigen.found()
  eigen_dep = eigen
else
  eigen_proj = subproject(eigen)
  eigen_dep = eigen_proj.get_variable('eigen_dep')
endif

gprd_deps += eigen_dep

if get_option('use_hdf5')
  gprd_deps += dependency('HighFive',
                          method : 'cmake',
                          modules : ['HighFive_HighFive'],
                          required : true
                         )
endif


if get_option('use_mkl')
  message('Be warned, MKL has not been tested extensively!')
  mkldep = dependency('mkl-dynamic-ilp64-iomp', required: true)
  add_project_arguments('-DEIGEN_USE_MKL_ALL', language: 'cpp')
  add_project_link_arguments('-lmkl_intel_lp64', '-lmkl_core', '-lmkl_sequential', language: 'cpp')
endif


libgprd = library('gprdimer',
          'backend/DistributionFunctions.cpp',
          'gpr/auxiliary/Distance.cpp',
          'gpr/auxiliary/ProblemSetUp.cpp',
          'gpr/auxiliary/Gradient.cpp',
          'gpr/covariance_functions/ConstantCF.cpp',
          'gpr/covariance_functions/SexpatCF.cpp',
          'gpr/dimer/Dimer.cpp',
          'gpr/optimization/LBFGS.cpp',
          'gpr/ml/GaussianProcessRegression.cpp',
          'gpr/ml/SCG.cpp',
          'gpr/observation_models/LikGaussian.cpp',
          'gpr/prior/PriorGaussian.cpp',
          'gpr/prior/PriorLogUnif.cpp',
          'gpr/prior/PriorSqrtt.cpp',
          'gpr/prior/PriorT.cpp',
          'gpr/prior/PriorUnif.cpp',
          'gpr/potentials/EAMPotential.cpp',
          'gpr/AtomicDimer.cpp',
          'gpr/AtomicDimerInit.cpp',
          'managers/io/FileManager.cpp',
           dependencies : gprd_deps,
           install : get_option('do_install'))

libgprd_dep = declare_dependency(link_with : libgprd)

if (not meson.is_subproject())
  test_deps = []
  # Handle GTest dependencies gracefully
  gtest = dependency('gtest', required : false)
  if gtest.found()
    test_deps += gtest
  else
    gtest_proj = subproject('gtest')
    test_deps += gtest_proj.get_variable('gtest_dep')
    test_deps += gtest_proj.get_variable('gmock_dep')
  endif

  add_languages('fortran', required: true)
  fc = meson.get_compiler('fortran')

  executable('gprd',
             'main.cpp',
             dependencies : [
               test_deps,
               gprd_deps,
               libgprd_dep
             ],
            install : true)

  # For the tests
  executable('mEAMCUH2',
             'gpr/eam_potential_fortran/main.f90',
             'gpr/eam_potential_fortran/eamroutines.f90',
            )
  # Needed to run the executable
  configure_file(input: 'pot.par',  output: 'pot.par',  copy: true)
  configure_file(input: 'tmp.con',  output: 'tmp.con',  copy: true)
  run_command('bash', '-c', 'cp -r $MESON_SOURCE_ROOT/input $MESON_BUILD_ROOT/input')

endif
