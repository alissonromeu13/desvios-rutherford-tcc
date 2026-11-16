"""Gera todas as figuras do TCC em figuras/ (mesma estrutura de pastas do Overleaf) e imprime os números
citados no texto. Uso:  python rodar_todas.py        (ou  python rodar_todas.py fig17  para uma só)"""
import glob, os, subprocess, sys, time
aqui = os.path.dirname(os.path.abspath(__file__))
filtro = sys.argv[1] if len(sys.argv) > 1 else ""
scripts = sorted(glob.glob(os.path.join(aqui, "fig*.py"))) + sorted(glob.glob(os.path.join(aqui, "numeros_do_texto", "n*.py")))
falhas = []
for s in scripts:
    if filtro not in os.path.basename(s):
        continue
    t = time.time()
    print(f"\n=== {os.path.relpath(s, aqui)}")
    r = subprocess.run([sys.executable, s], cwd=os.path.dirname(s))
    print(f"--- {time.time() - t:.0f} s")
    if r.returncode:
        falhas.append(os.path.basename(s))
print("\nFalhas:", falhas if falhas else "nenhuma")
