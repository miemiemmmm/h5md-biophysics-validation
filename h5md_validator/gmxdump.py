"""
Parser for ``gmx dump -s run.tpr`` output.

The *source* side of the round-trip: the metadata as GROMACS itself holds it,
read back out of the run input by GROMACS's own dumper. Nothing here reads an
H5MD file, so the two sides of the comparison share no code path.

Used by:
    Validation 2: round-trip test
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# line patterns
# --------------------------------------------------------------------------
_KV = re.compile(
    r"^(?P<indent>\s*)(?P<key>[\w\-\.]+(?:\s*\(3x3\))?)\s*=\s*(?P<value>.*?)\s*$"
)
_MATRIX_ROW = re.compile(r"^\s*(?P<name>[\w\-]+)\[\s*\d+\]=\{(?P<values>[^}]*)\}")
_LIST_LINE = re.compile(r"^\s*(?P<key>[\w\-]+):\s+(?P<values>[-\d\.eE+\s]+)$")
_GRP = re.compile(
    r"^grp\[(?P<kind>[^\]]+)\]\s*nr=(?P<nr>\d+),\s*name=\[(?P<names>[^\]]*)\]"
)

_MOLBLOCK = re.compile(r"^\s*molblock\s*\((?P<index>\d+)\):")
_MOLTYPE_REF = re.compile(r'^\s*moltype\s*=\s*(?P<index>\d+)\s*"(?P<name>.*)"')
_NMOLECULES = re.compile(r"^\s*#molecules\s*=\s*(?P<n>\d+)")
_MOLTYPE = re.compile(r"^\s*moltype\s*\((?P<index>\d+)\):")
_NAME = re.compile(r'^\s*name="(?P<name>.*)"')

_FUNCTYPE = re.compile(
    r"^\s*functype\[(?P<index>\d+)\]=(?P<kind>[A-Z_0-9]+),?\s*(?P<rest>.*)$"
)
_ATOM_PROP = re.compile(r"^\s*atom\[\s*(?P<index>\d+)\]=\{(?P<body>.*)\}")
_ATOM_NAME = re.compile(r'^\s*atom\[\s*(?P<index>\d+)\]=\{name="(?P<name>.*)"\}')
_RESIDUE = re.compile(
    r'^\s*residue\[\s*(?P<index>\d+)\]=\{name="(?P<name>.*?)",\s*nr=(?P<nr>-?\d+)'
)
_EXCL_START = re.compile(
    r"^\s*excls\[(?P<index>\d+)\]\[num=(?P<num>\d+)\]=\{(?P<body>.*)$"
)
_IATOM = re.compile(
    r"^\s*(?P<entry>\d+)\s+type=(?P<type>\d+)\s+\((?P<kind>[A-Z_0-9]+)\)\s*(?P<atoms>[\d\s]*)$"
)
_ILIST_LABEL = re.compile(r"^\s{6}(?P<label>[\w\-\. ]+):\s*$")


class DumpFloat(float):
    """A float that remembers how many significant digits gmx dump printed.

    ``gmx dump`` writes floating-point values as text, and different fields are
    written with different precision (masses with 6 significant digits, LJ
    coefficients with 9).  The round-trip tolerance for a value is bounded by
    whichever is coarser -- the stored HDF5 precision or this printed precision
    -- so it has to be carried along with the value rather than assumed.
    """

    __slots__ = ("sig_digits", "text")

    def __new__(cls, text):
        # Accepts the printed text, or a plain number when reconstructed by
        # copy/pickle -- in which case the printed precision is unknown and the
        # tolerance falls back to the stored precision alone.
        if isinstance(text, str):
            obj = super().__new__(cls, text)
            obj.text = text
            obj.sig_digits = _significant_digits(text)
        else:
            obj = super().__new__(cls, text)
            obj.text = repr(float(text))
            obj.sig_digits = None
        return obj

    def __reduce__(self):
        # copy.deepcopy and pickle must round-trip the printed text, not just
        # the value; dataclasses.asdict deep-copies every field.
        return (DumpFloat, (self.text,))


def _significant_digits(text: str) -> int:
    """Count the significant digits in a printed decimal number."""
    mantissa = text.strip().lower().split("e")[0].lstrip("+-")
    digits = mantissa.replace(".", "")
    stripped = digits.lstrip("0")
    return len(stripped) if stripped else 1


def _coerce(text: str):
    """Turn a dump value into int / float / bool / str."""
    if text in ("true", "false"):
        return text == "true"
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return DumpFloat(text)
    except ValueError:
        return text


@dataclass
class MoleculeType:
    name: str = ""
    #: one dict per atom: type, m, q, resind, atomnumber, ...
    atoms: list = field(default_factory=list)
    atom_names: list = field(default_factory=list)
    #: one (name, nr) per residue
    residues: list = field(default_factory=list)
    #: exclusion list per atom
    excls: list = field(default_factory=list)
    #: GROMACS functype name -> list of (parameter_index, [atom indices])
    ilists: dict = field(default_factory=dict)

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_residues(self) -> int:
        return len(self.residues)


@dataclass
class TprDump:
    inputrec: dict = field(default_factory=dict)
    grpopts: dict = field(default_factory=dict)
    groups: dict = field(default_factory=dict)
    header: dict = field(default_factory=dict)
    system_name: str = ""
    n_atoms: int = 0
    #: [{'moltype': int, 'name': str, 'nmol': int}]
    molblocks: list = field(default_factory=list)
    #: [(functype_name, {parameter: value})] indexed by GROMACS parameter index
    functypes: list = field(default_factory=list)
    moltypes: list = field(default_factory=list)


def run_dump(tpr, gmx="gmx", timeout=1800) -> str:
    """Run ``gmx dump -s tpr`` and return its stdout."""
    proc = subprocess.run(
        [gmx, "-quiet", "-nobackup", "dump", "-s", str(tpr)],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"gmx dump failed for {tpr}:\n{proc.stderr[-2000:]}")
    return proc.stdout


def parse(text: str) -> TprDump:
    """Parse the text of ``gmx dump -s``."""
    dump = TprDump()
    lines = text.splitlines()
    i, n = 0, len(lines)

    section = None  # inputrec | header | topology | None
    matrix_target = None  # name of the 3x3 block currently being read
    current_block = None  # molblock under construction
    mol = None  # MoleculeType under construction
    in_ffparams = False

    while i < n:
        line = lines[i]
        stripped = line.strip()
        i += 1
        if not stripped:
            continue

        # ---------------- top-level section switches ----------------------
        if line.startswith("inputrec:"):
            section, matrix_target = "inputrec", None
            continue
        if line.startswith("grpopts:"):
            section, matrix_target = "grpopts", None
            continue
        if line.startswith("header:"):
            section, matrix_target = "header", None
            continue
        if line.startswith("topology:"):
            section, matrix_target = "topology", None
            continue
        if line.startswith("grp["):
            m = _GRP.match(stripped)
            if m:
                dump.groups[m.group("kind").strip()] = m.group("names").split()
            continue

        # ---------------- inputrec ----------------------------------------
        if section == "inputrec":
            m = _MATRIX_ROW.match(line)
            if m and matrix_target:
                row = [
                    DumpFloat(v.strip())
                    for v in m.group("values").split(",")
                    if v.strip()
                ]
                dump.inputrec.setdefault(matrix_target, []).append(row)
                continue
            if stripped.endswith("(3x3):") or re.match(
                r"^[\w\-]+ \(\dx3\):$", stripped
            ):
                matrix_target = stripped.split(" (")[0]
                dump.inputrec.setdefault(matrix_target, [])
                continue
            m = _KV.match(line)
            if m:
                matrix_target = None
                dump.inputrec[m.group("key").strip()] = _coerce(m.group("value"))
            continue

        if section == "grpopts":
            m = _LIST_LINE.match(line)
            if m:
                dump.grpopts[m.group("key")] = [
                    _coerce(v) for v in m.group("values").split()
                ]
            continue

        if section == "header":
            m = _KV.match(line)
            if m:
                dump.header[m.group("key").strip()] = _coerce(m.group("value"))
            continue

        if section != "topology":
            continue

        # ---------------- topology ----------------------------------------
        if stripped == "ffparams:":
            in_ffparams, mol = True, None
            continue

        m = _MOLTYPE.match(line)
        if m:
            in_ffparams = False
            mol = MoleculeType()
            dump.moltypes.append(mol)
            continue

        m = _MOLBLOCK.match(line)
        if m:
            in_ffparams, mol = False, None
            current_block = {"moltype": None, "name": "", "nmol": 0}
            dump.molblocks.append(current_block)
            continue

        if current_block is not None and mol is None:
            m = _MOLTYPE_REF.match(line)
            if m:
                current_block["moltype"] = int(m.group("index"))
                current_block["name"] = m.group("name")
                continue
            m = _NMOLECULES.match(line)
            if m:
                current_block["nmol"] = int(m.group("n"))
                continue

        if in_ffparams:
            m = _FUNCTYPE.match(line)
            if m:
                params = {}
                for part in m.group("rest").split(","):
                    if "=" in part:
                        k, v = part.split("=", 1)
                        params[k.strip()] = _coerce(v.strip())
                idx = int(m.group("index"))
                while len(dump.functypes) <= idx:
                    dump.functypes.append(None)
                dump.functypes[idx] = (m.group("kind"), params)
                continue
            m = _KV.match(line)
            if m and m.group("key") in ("atnr", "ntypes", "fudgeQQ"):
                dump.header[m.group("key")] = _coerce(m.group("value"))
            continue

        if mol is not None:
            # molecule type name (the first name= inside a moltype block)
            m = _NAME.match(line)
            if m and not mol.name:
                mol.name = m.group("name")
                continue

            m = _ATOM_NAME.match(line)
            if m:
                idx = int(m.group("index"))
                while len(mol.atom_names) <= idx:
                    mol.atom_names.append("")
                mol.atom_names[idx] = m.group("name")
                continue

            m = _ATOM_PROP.match(line)
            if m and "type=" in m.group("body"):
                props = {}
                for part in m.group("body").split(","):
                    if "=" in part:
                        k, v = part.split("=", 1)
                        props[k.strip()] = _coerce(v.strip())
                idx = int(m.group("index"))
                while len(mol.atoms) <= idx:
                    mol.atoms.append({})
                mol.atoms[idx] = props
                continue

            m = _RESIDUE.match(line)
            if m:
                idx = int(m.group("index"))
                while len(mol.residues) <= idx:
                    mol.residues.append(("", 0))
                mol.residues[idx] = (m.group("name"), int(m.group("nr")))
                continue

            m = _EXCL_START.match(line)
            if m:
                body = m.group("body")
                while "}" not in body and i < n:
                    body += lines[i]
                    i += 1
                body = body.split("}")[0]
                idx = int(m.group("index"))
                values = [int(v) for v in body.replace(",", " ").split()]
                while len(mol.excls) <= idx:
                    mol.excls.append([])
                mol.excls[idx] = values
                continue

            m = _ILIST_LABEL.match(line)
            if m and m.group("label") not in ("atoms", "excls"):
                continue

            m = _IATOM.match(line)
            if m:
                kind = m.group("kind")
                atoms = [int(a) for a in m.group("atoms").split()]
                mol.ilists.setdefault(kind, []).append((int(m.group("type")), atoms))
                continue

    # header/topology counts
    dump.system_name = ""
    for idx, line in enumerate(lines):
        if line.startswith("topology:"):
            for probe in lines[idx + 1 : idx + 6]:
                m = _NAME.match(probe)
                if m:
                    dump.system_name = m.group("name")
                    break
            break
    dump.n_atoms = int(dump.header.get("natoms", 0))
    return dump


def from_tpr(tpr, gmx="gmx", cache_dir=None) -> TprDump:
    """Dump and parse a .tpr, optionally caching the dump text on disk."""
    tpr = Path(tpr)
    if cache_dir is not None:
        cache = Path(cache_dir) / f"{tpr.stem}.dump.txt"
        cache.parent.mkdir(parents=True, exist_ok=True)
        if cache.exists():
            return parse(cache.read_text())
        text = run_dump(tpr, gmx=gmx)
        cache.write_text(text)
        return parse(text)
    return parse(run_dump(tpr, gmx=gmx))


# ==========================================================================
# energy file: the source side of the observables comparison
# ==========================================================================
@dataclass
class EdrDump:
    """Per-frame energy terms from ``gmx dump -e``.

    ``terms`` maps a GROMACS energy component name to its value in each frame,
    in frame order; ``steps`` and ``times`` carry the frame index itself, so a
    comparison can align on the step rather than trusting position.
    """

    steps: list = field(default_factory=list)
    times: list = field(default_factory=list)
    units: dict = field(default_factory=dict)
    terms: dict = field(default_factory=dict)

    @property
    def n_frames(self) -> int:
        return len(self.steps)


_EDR_COMPONENT = re.compile(
    r"^\s*\d+\s{2,}(?P<name>\S.*?)\s{2,}\((?P<unit>[^)]*)\)\s*$"
)
_EDR_FRAME = re.compile(r"^\s*time:\s*(?P<time>\S+)\s+step:\s*(?P<step>\d+)\s*$")
_EDR_VALUE = re.compile(
    r"^\s{2,}(?P<name>\S.*?)\s{2,}(?P<value>-?[\d.]+e[+-]\d+|-?[\d.]+)\s*$"
)


def parse_edr(text: str) -> EdrDump:
    """Parse the text of ``gmx dump -e``."""
    dump = EdrDump()
    in_components = False
    for line in text.splitlines():
        if line.startswith("energy components:"):
            in_components = True
            continue

        frame = _EDR_FRAME.match(line)
        if frame:
            in_components = False
            dump.steps.append(int(frame.group("step")))
            dump.times.append(float(frame.group("time")))
            continue

        if in_components:
            m = _EDR_COMPONENT.match(line)
            if m:
                dump.units[m.group("name").strip()] = m.group("unit").strip()
            continue

        if not dump.steps:
            continue
        # Inside a frame. The header line names three columns but only the
        # first, the instantaneous value, is present for every term.
        m = _EDR_VALUE.match(line)
        if m and m.group("name").strip() != "Component":
            name = m.group("name").strip()
            series = dump.terms.setdefault(name, [])
            # Pad a term that first appears after frame 0 so that every series
            # stays indexed by frame.
            while len(series) < len(dump.steps) - 1:
                series.append(None)
            # DumpFloat, not float: gmx dump -e prints 6 significant digits, so
            # the comparison tolerance has to come from the printed text rather
            # than from the stored float width.
            series.append(DumpFloat(m.group("value")))
    return dump


def run_dump_edr(edr, gmx="gmx", timeout=1800) -> str:
    """Run ``gmx dump -e edr`` and return its stdout."""
    proc = subprocess.run(
        [gmx, "-quiet", "-nobackup", "dump", "-e", str(edr)],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"gmx dump -e failed for {edr}:\n{proc.stderr[-2000:]}")
    return proc.stdout


def from_edr(edr, gmx="gmx", cache_dir=None) -> EdrDump:
    """Dump and parse an .edr, optionally caching the dump text on disk."""
    edr = Path(edr)
    if cache_dir is not None:
        cache = Path(cache_dir) / f"{edr.stem}.edr.txt"
        cache.parent.mkdir(parents=True, exist_ok=True)
        if cache.exists():
            return parse_edr(cache.read_text())
        text = run_dump_edr(edr, gmx=gmx)
        cache.write_text(text)
        return parse_edr(text)
    return parse_edr(run_dump_edr(edr, gmx=gmx))
