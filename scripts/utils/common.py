"""Configuration, provenance, TE-annotation and bf16 helpers shared by the
numbered pipeline scripts.

Kept dependency-light on purpose: importing this module must not import
torch or geneformer, so 01_preprocess_data.py runs on a CPU-only machine.
"""
from __future__ import annotations

import hashlib
import json
import logging
import platform
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "config" / "config.yaml"

# Repeat classes that are not transposable elements (strict_te filter,
# TE_Ocean_Acidification 2026-09-06). The config may extend this list.
NON_TE_CLASSES = {
    "Simple_repeat", "Low_complexity", "Satellite", "rRNA", "tRNA",
    "snRNA", "scRNA", "srpRNA", "RNA",
}


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
class Config(dict):
    """dict with attribute access and repo-root-relative path resolution."""

    def __getattr__(self, key: str) -> Any:
        try:
            value = self[key]
        except KeyError as exc:
            raise AttributeError(key) from exc
        return Config(value) if isinstance(value, dict) else value

    def path(self, *keys: str) -> Path:
        """Resolve a nested path entry, e.g. cfg.path("paths", "isp_dir")."""
        node: Any = self
        for key in keys:
            node = node[key]
        if node is None:
            raise ValueError(f"Config path {'.'.join(keys)} is null")
        p = Path(node).expanduser()
        return p if p.is_absolute() else (REPO_ROOT / p).resolve()


def load_config(path: str | Path | None = None) -> Config:
    path = Path(path) if path else DEFAULT_CONFIG
    with open(path) as fh:
        cfg = Config(yaml.safe_load(fh))
    cfg["_config_path"] = str(Path(path).resolve())
    cfg["_config_sha256"] = sha256_file(path)
    return cfg


def setup_logging(cfg: Config, script_name: str) -> logging.Logger:
    log_dir = cfg.path("paths", "logs_dir")
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    logging.basicConfig(
        level=getattr(logging, cfg.runtime.get("log_level", "INFO")),
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_dir / f"{script_name}_{stamp}.log"),
        ],
    )
    return logging.getLogger(script_name)


# ---------------------------------------------------------------------------
# Provenance and resumability
# ---------------------------------------------------------------------------
def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def sha256_dir(path: str | Path) -> str:
    """Order-stable digest over every file under a directory (HF datasets)."""
    h = hashlib.sha256()
    for f in sorted(Path(path).rglob("*")):
        if f.is_file():
            h.update(str(f.relative_to(path)).encode())
            h.update(sha256_file(f).encode())
    return h.hexdigest()


def git_commit(repo: Path = REPO_ROOT) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def provenance(cfg: Config, **extra: Any) -> dict:
    record = {
        "utc": utc_now(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "repo_commit": git_commit(),
        "config_path": cfg.get("_config_path"),
        "config_sha256": cfg.get("_config_sha256"),
    }
    try:  # optional heavy deps — recorded when present
        import torch  # noqa: WPS433

        record["torch"] = torch.__version__
        record["cuda"] = torch.version.cuda
    except ImportError:
        pass
    try:
        import geneformer  # noqa: WPS433

        record["geneformer_path"] = str(Path(geneformer.__file__).parent)
        record["geneformer_commit"] = git_commit(Path(geneformer.__file__).parents[1])
    except ImportError:
        pass
    record.update(extra)
    return record


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n")


def marker_path(directory: Path, unit: str) -> Path:
    return directory / f"{unit}.complete.json"


def is_complete(cfg: Config, directory: Path, unit: str) -> bool:
    return bool(cfg.runtime.get("resume", True)) and marker_path(directory, unit).exists()


def mark_complete(directory: Path, unit: str, payload: dict) -> None:
    """A unit gets its marker only after all of its outputs are written, so an
    interrupted run restarts without repeating finished units."""
    write_json(marker_path(directory, unit), {"completed_utc": utc_now(), **payload})


# ---------------------------------------------------------------------------
# TE annotation
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TEFamily:
    family: str
    superfamily: str
    te_class: str

    @property
    def strict_te(self) -> bool:
        return self.te_class not in NON_TE_CLASSES


_GTF_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def read_te_annotation(path: str | Path) -> dict[str, TEFamily]:
    """Family -> class map from a TEtranscripts-format GTF
    (gene_id = subfamily, family_id = superfamily, class_id = class), as used
    by TEcount in TE_Ocean_Acidification, or from a RepeatMasker .out file
    (repeat in column 10, class/family in column 11)."""
    path = Path(path)
    families: dict[str, TEFamily] = {}
    if path.suffix in {".gtf", ".gz"} or path.name.endswith(".gtf.gz"):
        import gzip

        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                attrs = dict(_GTF_ATTR.findall(fields[8]))
                name = attrs.get("gene_id")
                if name and name not in families:
                    families[name] = TEFamily(
                        family=name,
                        superfamily=attrs.get("family_id", "Unknown"),
                        te_class=attrs.get("class_id", "Unknown"),
                    )
    else:  # RepeatMasker .out
        with open(path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 11 or not parts[0].isdigit():
                    continue
                name, cls_fam = parts[9], parts[10]
                cls, _, fam = cls_fam.partition("/")
                if name not in families:
                    families[name] = TEFamily(name, fam or cls, cls.rstrip("?"))
    if not families:
        raise ValueError(f"No TE families parsed from {path}")
    return families


def normalize_te_class(te_class: str) -> str:
    """'Unknown' stays its own class; uncertain calls ('LINE?') are not upgraded."""
    if te_class in {"", "Unknown", "Unspecified"}:
        return "Unknown"
    return te_class


# ---------------------------------------------------------------------------
# Geneformer dictionaries
# ---------------------------------------------------------------------------
def load_pickle(path: str | Path) -> Any:
    import pickle

    with open(path, "rb") as fh:
        return pickle.load(fh)


def resolve_symbols(symbols: Iterable[str], gene_name_id: dict, token_dict: dict,
                    log: logging.Logger) -> dict[str, str]:
    """HGNC symbol -> Ensembl ID, keeping only IDs present in the token
    dictionary. Unresolved symbols are logged and dropped, never guessed."""
    resolved: dict[str, str] = {}
    for sym in dict.fromkeys(symbols):
        ens = gene_name_id.get(sym)
        if ens is None:
            log.warning("Target %s: no Ensembl ID in gene_name_id dictionary — skipped", sym)
        elif ens not in token_dict:
            log.warning("Target %s (%s): not in Geneformer vocabulary — skipped", sym, ens)
        else:
            resolved[sym] = ens
    return resolved


# ---------------------------------------------------------------------------
# bf16 support (ported from geneformer-lung-tcell/sclc_validation/bf16_bench)
# ---------------------------------------------------------------------------
def install_dtype_cast(dtype: str) -> None:
    """Cast every model loaded through geneformer.perturber_utils.load_model to
    bf16. InSilicoPerturber and EmbExtractor call pu.load_model internally, so
    reassigning the module attribute is visible to them without editing
    upstream source. Must run before the first perturb_data()/extract call."""
    if dtype not in {"fp32", "bf16"}:
        raise ValueError(f"runtime.dtype must be fp32 or bf16, got {dtype!r}")
    import torch
    from geneformer import perturber_utils as pu

    if getattr(pu.load_model, "_te_isp_dtype", None) == dtype:
        return
    original = getattr(pu.load_model, "_te_isp_original", pu.load_model)

    def load_model_cast(*args, **kwargs):
        model = original(*args, **kwargs)
        return model.to(torch.bfloat16) if dtype == "bf16" else model.float()

    load_model_cast._te_isp_dtype = dtype
    load_model_cast._te_isp_original = original
    pu.load_model = load_model_cast


class EvalTieCounter:
    ties = 0
    total = 0


def install_eval_tie_fix() -> None:
    """bf16 TRAINING makes exact logit ties likely; geneformer.evaluation_utils.vote
    then returns the string "tie" and sklearn's confusion_matrix crashes on mixed
    label types. Break ties deterministically (first index, argmax convention)
    and count them so the tie rate can be reported beside any metric."""
    from geneformer import evaluation_utils

    if getattr(evaluation_utils.vote, "_te_isp_patched", False):
        return

    def vote_no_ties(logit_list):
        EvalTieCounter.total += 1
        best = max(logit_list)
        idx = [i for i, x in enumerate(logit_list) if x == best]
        if len(idx) > 1:
            EvalTieCounter.ties += 1
        return idx[0]

    vote_no_ties._te_isp_patched = True
    evaluation_utils.vote = vote_no_ties


def check_bf16_export_patch(log: logging.Logger) -> None:
    """numpy has no bfloat16; Geneformer's embedding export calls .numpy() on
    bf16 tensors unless patches/geneformer-bf16-export.patch is applied."""
    import inspect

    from geneformer import emb_extractor

    if ".float().cpu().numpy()" not in inspect.getsource(emb_extractor.label_cell_embs):
        raise RuntimeError(
            "runtime.dtype=bf16 but the Geneformer checkout lacks the bf16 export fix. "
            "Apply it with: git -C <Geneformer repo> apply "
            f"{REPO_ROOT / 'patches' / 'geneformer-bf16-export.patch'}"
        )
    log.info("bf16 export patch detected in Geneformer checkout")
