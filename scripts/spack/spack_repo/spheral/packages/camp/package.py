# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.packages.camp.package import Camp as BuiltinCamp


class Camp(BuiltinCamp):
    """CAMP release required by the local RAJA and Umpire recipes."""

    # Backported from spack/spack-packages at
    # db11975340925c44f5e58962c58c7f01bd6bc7e5.
    version(
        "2026.07.2",
        tag="v2026.07.2",
        commit="b4edd077937245c120e25e6080be8caee5926347",
        submodules=False,
    )
