"""Estilo comum das figuras: fonte serifada (combina com abnTeX2), vírgula decimal, PDF vetorial."""
import matplotlib
matplotlib.use("pgf")   # texto composto pelo XeLaTeX: mesmas fontes do documento e camada de texto Unicode correta
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, ScalarFormatter
import os

AZUL, LARANJA, VERDE, AMARELO = "#2a78d6", "#eb6834", "#1baf7a", "#c98500"
CINZA, TINTA = "#6b6a63", "#1f1e1b"
CORES = [AZUL, LARANJA, VERDE, AMARELO]
LS = ["-", "--", "-.", ":"]

plt.rcParams.update({
    "pgf.texsystem": "xelatex", "pgf.rcfonts": False,
    "pgf.preamble": "\n".join([r"\usepackage{fontspec}", r"\usepackage{amsmath}", r"\usepackage{unicode-math}",
                                r"\setmainfont{Latin Modern Roman}", r"\setmathfont{Latin Modern Math}"]),
    "font.family": "serif", "font.serif": ["Latin Modern Roman", "DejaVu Serif"], "mathtext.fontset": "cm",
    "axes.unicode_minus": True, "pdf.fonttype": 3,
    "font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8.5, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.linewidth": 0.7, "lines.linewidth": 1.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e4e3de", "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
    "axes.edgecolor": TINTA, "text.color": TINTA, "axes.labelcolor": TINTA, "xtick.color": TINTA, "ytick.color": TINTA,
})

def _virgula(x, pos=None):
    s = f"{x:g}"
    return s.replace(".", ",").replace("-", "−")

def virgula(ax, eixos="xy"):
    """Vírgula decimal (padrão ABNT) nos eixos lineares."""
    for e in eixos:
        axis = ax.xaxis if e == "x" else ax.yaxis
        scale = ax.get_xscale() if e == "x" else ax.get_yscale()
        if scale == "linear":
            axis.set_major_formatter(FuncFormatter(_virgula))

def rotulo(ax, letra):
    ax.text(0.0, 1.03, r"\textbf{(" + letra + ")}", transform=ax.transAxes, fontsize=10, va="bottom")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras")   # mesma estrutura de pastas do Overleaf

def salvar(fig, nome, **kw):
    """Salva figuras/<nome>.pdf (vetorial; nome pode conter a subpasta, como "6.4/fig08_...").
    Se o pdftoppm (Poppler) estiver instalado, grava também um PNG de 300 dpi para pré-visualização;
    o TCC usa apenas o PDF."""
    import shutil, subprocess
    pdf = os.path.join(OUT, nome + ".pdf")
    os.makedirs(os.path.dirname(pdf), exist_ok=True)
    fig.savefig(pdf, **kw)
    plt.close(fig)
    if shutil.which("pdftoppm"):
        subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile", pdf, os.path.join(OUT, nome)], check=True)
    print("salvo:", nome + ".pdf")
