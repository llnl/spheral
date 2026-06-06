include_directories(.)
set(Damage_inst
    DamageModel
    TensorDamageModel
    ProbabilisticDamageModel
    ProbabilisticDamagePolicy
    IvanoviSALEDamageModel
    IvanoviSALEDamagePolicy
    StrainPolicy
    TensorStrainPolicy
    TensorDamagePolicy
    YoungsModulusPolicy
    LongitudinalSoundSpeedPolicy
    weibullFlawDistributionBenzAsphaug
    weibullFlawDistributionOwen
    computeFragmentField
    JohnsonCookDamage
    JohnsonCookFailureStrainPolicy
    JohnsonCookDamagePolicy
    PairMaxDamageNodeCoupling
    ThreePointDamagedNodeCoupling
    DamageGradientNodeCoupling
   )

set(Damage_sources )

instantiate(Damage_inst Damage_sources)

set(Damage_headers
    DamageGradientPolicy.hh
    DamageModel.hh
    EffectiveTensorDamagePolicy.hh
    GradyKippScalarDamage.hh
    JohnsonCookDamage.hh
    JohnsonCookDamageInline.hh
    JohnsonCookDamagePolicy.hh
    JohnsonCookFailureStrainPolicy.hh
    LongitudinalSoundSpeedPolicy.hh
    ScalarDamageModel.hh
    ScalarDamagePolicy.hh
    StrainPolicy.hh
    TensorDamageModel.hh
    TensorDamageModelInline.hh
    ProbabilisticDamageModel.hh
    ProbabilisticDamageModelInline.hh
    IvanoviSALEDamageModel.hh
    IvanoviSALEDamageModelInline.hh
    TensorDamagePolicy.hh
    ProbabilisticDamagePolicy.hh
    IvanoviSALEDamagePolicy.hh
    TensorStrainPolicy.hh
    YoungsModulusPolicy.hh
    computeFragmentField.hh
    oneMinusDamage.hh
    weibullFlawDistributionBenzAsphaug.hh
    weibullFlawDistributionOwen.hh
    DamageGradientNodeCoupling.hh
    PairMaxDamageNodeCoupling.hh
    MaxDamageNodeCoupling.hh
    ThreePointDamagedNodeCoupling.hh
    )

spheral_install_python_files(
  identifyFragments.py
  WeibullDamage.py
  GradyKippTensorDamage.py
  JohnsonCookDamageFactories.py
  ShadowProbabilisticDamageModel.py
  ShadowIvanoviSALEDamageModel.py
  )

spheral_add_obj_library(Damage SPHERAL_OBJ_LIBS)

