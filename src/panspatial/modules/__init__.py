"""Analysis module registry.

Importing this package registers every analysis layer, so ``describe_registry()`` is the
pipeline's table of contents and ``get_module(name)`` is the single entry point.
"""

from panspatial.modules.base import (
    AnalysisModule,
    DependencyMissing,
    Evidence,
    ModuleContext,
    ModuleError,
    ModuleResult,
    PlatformClass,
    Provenance,
    ReplicationUnit,
    assert_can_validate,
    describe_registry,
    get_module,
    modules_for_phase,
    platform_class,
    register,
    registered_modules,
    require,
)

# Import for side effects: each module registers itself on import.
from panspatial.qc import metrics as _qc                     # noqa: F401
from panspatial.integrate import harmonize as _integrate     # noqa: F401
from panspatial.annotate import celltype as _annotate        # noqa: F401
from panspatial.deconv import consensus as _deconv           # noqa: F401
from panspatial.domains import detect as _domains            # noqa: F401
from panspatial.niches import call as _niches                # noqa: F401
from panspatial.malignant import cnv as _malignant           # noqa: F401
from panspatial.interactions import ccc as _interactions     # noqa: F401
from panspatial.grn import regulons as _grn                  # noqa: F401
from panspatial.comodules import modules as _comodules       # noqa: F401
from panspatial.dynamics import fate as _dynamics            # noqa: F401
from panspatial.alignment import slices as _alignment        # noqa: F401
from panspatial.histology import predict as _histology       # noqa: F401
from panspatial.predict import outcome as _predict           # noqa: F401
from panspatial.translate import targets as _translate       # noqa: F401
from panspatial.nk import module as _nk                       # noqa: F401

__all__ = [
    "AnalysisModule", "DependencyMissing", "Evidence", "ModuleContext", "ModuleError",
    "ModuleResult", "PlatformClass", "Provenance", "ReplicationUnit",
    "assert_can_validate", "describe_registry", "get_module", "modules_for_phase",
    "platform_class", "register", "registered_modules", "require",
]
