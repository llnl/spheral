include_directories(.)
set(Damage_inst
   )

set(Damage_sources )

if(ENABLE_INSTANTIATIONS)
  list(APPEND Damage_sources
       DamageModelInst.cc
       TensorDamageModelInst.cc
       ProbabilisticDamageModelInst.cc
       ProbabilisticDamagePolicyInst.cc
       IvanoviSALEDamageModelInst.cc
       IvanoviSALEDamagePolicyInst.cc
       StrainPolicyInst.cc
       TensorStrainPolicyInst.cc
       TensorDamagePolicyInst.cc
       YoungsModulusPolicyInst.cc
       LongitudinalSoundSpeedPolicyInst.cc
       weibullFlawDistributionBenzAsphaugInst.cc
       weibullFlawDistributionOwenInst.cc
       computeFragmentFieldInst.cc
       JohnsonCookDamageInst.cc
       JohnsonCookFailureStrainPolicyInst.cc
       JohnsonCookDamagePolicyInst.cc
       PairMaxDamageNodeCouplingInst.cc
       ThreePointDamagedNodeCouplingInst.cc
       DamageGradientNodeCouplingInst.cc
       )
endif()

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

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND Damage_headers
     DamageModel.cc
     TensorDamageModel.cc
     ProbabilisticDamageModel.cc
     ProbabilisticDamagePolicy.cc
     IvanoviSALEDamageModel.cc
     IvanoviSALEDamagePolicy.cc
     StrainPolicy.cc
     TensorStrainPolicy.cc
     TensorDamagePolicy.cc
     YoungsModulusPolicy.cc
     LongitudinalSoundSpeedPolicy.cc
     weibullFlawDistributionBenzAsphaug.cc
     weibullFlawDistributionOwen.cc
     computeFragmentField.cc
     JohnsonCookDamage.cc
     JohnsonCookFailureStrainPolicy.cc
     JohnsonCookDamagePolicy.cc
     PairMaxDamageNodeCoupling.cc
     ThreePointDamagedNodeCoupling.cc
     DamageGradientNodeCoupling.cc
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

