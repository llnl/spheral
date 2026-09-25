include_directories(.)
set(Porosity_inst
    PorosityModel
    StrainPorosity
    PalphaPorosity
    PorositySolidMassDensityPolicy
   )

set(Porosity_sources
    )

instantiate(Porosity_inst Porosity_sources)

set(Porosity_headers
    PorosityModel.hh
    PorosityModelInline.hh
    StrainPorosity.hh
    StrainPorosityInline.hh
    PalphaPorosity.hh
    )

spheral_install_python_files(
    ShadowStrainPorosity.py
    ShadowPalphaPorosity.py
    computeHugoniotWithPorosity.py
    PorousEquationOfState.py
    PorousStrengthModel.py
    )

spheral_add_obj_library(Porosity)
