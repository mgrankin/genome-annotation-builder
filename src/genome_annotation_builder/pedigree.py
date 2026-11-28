"""Read and validate the declared PED roles used by the family adapters."""

from dataclasses import dataclass
from pathlib import Path

from .errors import InputValidationError


@dataclass(frozen=True)
class PedigreeMember:
    family_id: str
    individual: str
    father: str | None
    mother: str | None
    sex: int

    @property
    def parents(self) -> tuple[str, ...]:
        return tuple(value for value in (self.father, self.mother) if value is not None)


def parse_pedigree(path: Path) -> dict[str, PedigreeMember]:
    members: dict[str, PedigreeMember] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 6:
                raise InputValidationError(
                    f"PED line {line_number} must contain exactly six fields"
                )
            family_id, individual, father_raw, mother_raw, sex_raw, _phenotype = fields
            if individual in members:
                raise InputValidationError("PED contains a duplicate individual identifier")
            try:
                sex = int(sex_raw)
            except ValueError as exc:
                raise InputValidationError(f"PED line {line_number} has invalid sex") from exc
            if sex not in {0, 1, 2}:
                raise InputValidationError(f"PED line {line_number} has unsupported sex code")
            father = None if father_raw == "0" else father_raw
            mother = None if mother_raw == "0" else mother_raw
            if individual in {father, mother}:
                raise InputValidationError("PED individual cannot be its own parent")
            members[individual] = PedigreeMember(
                family_id=family_id,
                individual=individual,
                father=father,
                mother=mother,
                sex=sex,
            )
    if not members:
        raise InputValidationError("PED contains no individuals")
    for member in members.values():
        for parent in member.parents:
            if parent not in members:
                raise InputValidationError("PED references a parent absent from the pedigree")
            if members[parent].family_id != member.family_id:
                raise InputValidationError("PED parent and child have different family identifiers")

    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(individual: str) -> None:
        if individual in visiting:
            raise InputValidationError("PED contains an ancestry cycle")
        if individual in visited:
            return
        visiting.add(individual)
        for parent in members[individual].parents:
            visit(parent)
        visiting.remove(individual)
        visited.add(individual)

    for individual in members:
        visit(individual)
    return members
