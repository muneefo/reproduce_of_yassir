"""Check scientific dependencies without claiming model validation."""
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from scipy.integrate import quad


def main():
    value, error = quad(lambda x: np.exp(-x), 0, np.inf)
    if abs(value - 1.0) > 1e-10:
        raise RuntimeError(f"Dependency integration check failed: {value}, {error}")
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    fig.canvas.draw()
    plt.close(fig)
    print(f"Dependencies verified: numpy={np.__version__}, scipy={scipy.__version__}, matplotlib={matplotlib.__version__}")
    print("No paper model or reproduction has been validated.")


if __name__ == "__main__":
    main()
