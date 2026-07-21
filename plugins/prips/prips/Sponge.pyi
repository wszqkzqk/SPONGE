"""
prips.Sponge
===============
a virtual plugin interface module for communicating with the molecular dynamics simulation package SPONGE

Example - 0
---------------
```python
print("Hello World")
input("press any key to continue...")
```

Example - 1
---------------
```python
from prips import Sponge
Sponge.set_backend("cupy")
print(Sponge.md_info.crd)
input("press any key to continue...")
```

Example - 2
---------------
```python
from prips import Sponge
Sponge.set_backend("jax")

def Calculate_Force():
    force = Sponge.dd.frc.at[:, 2].add(1.0)
    return Sponge.force_result(force)
```
"""

from typing import Any, Literal

_backend: int
"DLPack device type used by the active SPONGE backend"
_device_id: int
"Logical device ordinal used by the active SPONGE process"

FORCE_ENERGY_COMPLETE: int
"Calculate_Force writes its complete contribution to dd.energy when requested"
FORCE_VIRIAL_COMPLETE: int
"Calculate_Force writes its complete contribution to dd.virial when requested"
FORCE_PURE: int
"Calculate_Force has no mutable state that can affect force results"
FORCE_TRANSACTIONAL: int
"Calculate_Force isolates mutable state with transaction lifecycle hooks"

"""
Python force plugins used with the Monte Carlo barostat must define an integer
``SPONGE_FORCE_CAPABILITIES`` in their own script. Exactly one of
``FORCE_PURE`` and ``FORCE_TRANSACTIONAL`` is required. Transactional scripts
must also define callable ``Begin_Force_Transaction()``,
``Commit_Force_Transaction()``, and ``Rollback_Force_Transaction()`` hooks.
"""

def set_backend(backend: Literal["numpy", "jax", "cupy", "pytorch"]) -> None:
    """
    Set the backend for PRIPS.

    Parameters
    -----------
    backend : Literal["numpy", "jax", "cupy", "pytorch"]
        Backend name used to convert DLPack tensors to framework tensors.
        Because JAX arrays are immutable, a JAX ``Calculate_Force`` callback
        must return :func:`force_result` so PRIPS can synchronously write the
        complete updated local buffers back to SPONGE.
    """

class FORCE_RESULT:
    """Complete updated local buffers returned by a functional force hook."""

    force: Any
    energy: Any | None
    virial: Any | None

def force_result(
    force: Any, *, energy: Any | None = None, virial: Any | None = None
) -> FORCE_RESULT:
    """
    Build the result of a functional ``Calculate_Force`` callback.

    ``force`` is the complete updated ``Sponge.dd.frc`` buffer, not an
    isolated contribution. ``energy`` and ``virial``, when supplied, likewise
    replace the complete updated local buffers. This replacement contract is
    exactly equivalent to mutating the exposed buffers in place and preserves
    force contributions already accumulated by the SPONGE core.
    """

class FORCE_EVALUATION_CONTEXT:
    "Read-only context for the currently executing force callback"
    @property
    def commits_sampling_state(self) -> bool:
        "Whether this callback may advance adaptive or history state"

    @property
    def is_exact(self) -> bool:
        "Whether this callback evaluates the current state exactly"

    @property
    def needs_energy(self) -> bool:
        "Whether this callback requires a complete energy contribution"

    @property
    def needs_virial(self) -> bool:
        "Whether this callback requires a complete virial contribution"

force_evaluation = FORCE_EVALUATION_CONTEXT()
" Context for the currently executing Calculate_Force callback "

class MD_INFORMATION:
    "Basic information class for MD"
    class system_information:
        @property
        def steps(self) -> int:
            "The number of the current step"

    @property
    def atom_numbers(self) -> int:
        "The number of atoms"

    @property
    def crd(self) -> SpongeDLPackTensor:
        "The atom coordinates with shape (atom_numbers, 3)"

    @property
    def frc(self) -> SpongeDLPackTensor:
        "The atom forces with shape (atom_numbers, 3)"

md_info = MD_INFORMATION()
" Basic information instance for MD "

class NEIGHBOR_LIST:
    @property
    def index(self) -> SpongeDLPackTensor | None:
        "Neighbor index matrix with shape (atom_numbers, max_neighbor_numbers)"

    @property
    def number(self) -> list[int] | None:
        "Neighbor count for each atom"

    @property
    def max_neighbor_numbers(self) -> int | None:
        "The allocated max number of neighbors for one atom"

neighbor_list: NEIGHBOR_LIST | None
" Neighbor list interface, available after SPONGE finishes initialization "

class DOMAIN_INFORMATION:
    @property
    def atom_numbers(self) -> int | None:
        "Local atom numbers in the current PP rank"

    @property
    def ghost_numbers(self) -> int | None:
        "Ghost atom numbers in the current PP rank"

    @property
    def pp_rank(self) -> int | None:
        "PP rank id of the current domain"

    @property
    def atom_local(self) -> SpongeDLPackTensor | None:
        "Global atom ids for local and ghost atoms"

    @property
    def atom_local_label(self) -> SpongeDLPackTensor | None:
        "Global-to-local membership flags, stored over max_atom_numbers"

    @property
    def atom_local_id(self) -> SpongeDLPackTensor | None:
        "Global-to-local index map, stored over max_atom_numbers"

    @property
    def crd(self) -> SpongeDLPackTensor | None:
        "Local plus ghost coordinates with shape (atom_numbers + ghost_numbers, 3)"

    @property
    def frc(self) -> SpongeDLPackTensor | None:
        "Local plus ghost forces with shape (atom_numbers + ghost_numbers, 3)"

    @property
    def energy(self) -> SpongeDLPackTensor | None:
        "Authoritative owned-atom energy buffer with shape (atom_numbers,)"

    @property
    def virial(self) -> SpongeDLPackTensor | None:
        "Authoritative owned-atom virial with shape (atom_numbers, 6) in a11, a21, a22, a31, a32, a33 order"

dd: DOMAIN_INFORMATION | None
" Current domain interface; it is rebound after DD remeshing, so do not cache its tensors across callbacks "

class CONTROLLER:
    "IO controller"
    def printf(self, *values: any, sep: str = " ", end: str = "\n") -> None:
        """Print the values to the screen and the mdinfo file"""

controller = CONTROLLER()
" IO controller instance"
