import pickle as pkl
from pathlib import Path
from tdsr.types import PathLike
from typing import List, Any


def ensure_dirs(*paths: List[PathLike]):
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


TEST_DIR = Path(__file__).parent
REPO_ROOT = TEST_DIR.parent
DATA_DIR = REPO_ROOT / "data"
PLOT_DIR = TEST_DIR / "plots"
VALUES_DIR = TEST_DIR / "values"

ensure_dirs(VALUES_DIR, PLOT_DIR)


def value_file(key: str) -> PathLike:
    """value file path for key"""
    return (VALUES_DIR / (key + ".pkl")).absolute()


def load_values(key: str) -> Any:
    """unpickle stored values for key"""
    with open(value_file(key), "rb") as f:
        return pkl.load(f)


def dump_values(key: str, values: Any, overwrite: bool = False) -> None:
    """pickle values for key"""
    vf = value_file(key)
    if overwrite or not vf.is_file():
        with open(vf, "wb") as f:
            return pkl.dump(values, f)
