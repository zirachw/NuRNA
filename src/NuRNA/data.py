import re
import datasets as hf


def load_dataset(family: str | None = None, max_length: int | None = None) -> list[dict]:
    ds = hf.load_dataset("multimolecule/archiveii", split="test")
    records = []
    for row in ds:
        seq = row["sequence"]
        if not validate_sequence(seq):
            continue
        if family and row["family"] != family:
            continue
        if max_length and len(seq) > max_length:
            continue
        records.append({
            "id":        row["id"],
            "family":    row["family"],
            "seq":       seq,
            "reference": parse_dot_bracket(row["secondary_structure"]),
        })
    return records


def parse_dot_bracket(notation: str) -> set[tuple[int, int]]:
    pairs = set()
    stack = []
    for i, ch in enumerate(notation):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack:
                j = stack.pop()
                pairs.add((j, i))
    return pairs


def validate_sequence(seq: str) -> bool:
    return bool(re.fullmatch(r"[AUGC]+", seq))


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(1 for nt in seq if nt in ("G", "C")) / len(seq)
