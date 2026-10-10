"""Saved, headless figures; never label preliminary curves as reproduced."""
import os
import tempfile
from pathlib import Path
cache = Path(tempfile.gettempdir()) / "isac-matplotlib"
cache.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

COLORS = ["#222222", "#0072B2", "#D55E00", "#009E73", "#CC79A7"]


def figure():
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .2})
    return plt.subplots(figsize=(7, 4.5), constrained_layout=True)


def save(fig, path):
    fig.savefig(path, dpi=180)
    fig.savefig(Path(path).with_suffix(".svg"))
    plt.close(fig)
