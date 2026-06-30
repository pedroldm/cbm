"""Solution validation mirroring the C++ ``Validator``.

A solution is valid iff its permutation contains every column exactly once and
its reported cost equals the independently recomputed 1-block count.
"""

from __future__ import annotations

from .instance import CBMInstance


class ValidationError(Exception):
    """Raised when a produced solution fails validation."""


class SolutionValidator:
    def __init__(self, instance: CBMInstance) -> None:
        self.instance = instance

    def validate(self, perm, cost: int) -> int:
        """Validate ``perm``/``cost`` and return the recomputed block count.

        Raises :class:`ValidationError` on a malformed permutation or a cost
        that disagrees with the recomputed 1-block count.
        """
        self._check_permutation(perm)
        recomputed = self.instance.count_one_blocks(perm)
        if cost != recomputed:
            raise ValidationError(f"reported cost {cost} != recomputed block count {recomputed}")
        return recomputed

    def _check_permutation(self, perm) -> None:
        inst = self.instance
        if len(perm) != inst.cols:
            raise ValidationError(f"permutation has {len(perm)} entries, expected {inst.cols}")
        if set(perm) != set(range(inst.cols)):
            raise ValidationError("permutation does not contain every column exactly once")
