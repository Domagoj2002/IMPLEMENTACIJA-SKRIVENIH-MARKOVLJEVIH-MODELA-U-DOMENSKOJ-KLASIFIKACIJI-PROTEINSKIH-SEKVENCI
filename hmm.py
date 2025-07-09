import subprocess
from graphviz import Digraph
import shutil

def build_hmm(stockholm_file, hmm_file="model.hmm"):
    subprocess.run(["hmmbuild", hmm_file, stockholm_file], check=True)

def parse_hmm(hmm_file):
    with open(hmm_file) as f:
        lines = f.readlines()

    hmm_start = next(i for i, line in enumerate(lines) if line.startswith("HMM"))
    data_lines = lines[hmm_start+1:]
    
    transitions = []
    emissions = {}
    i = 0
    while i < len(data_lines):
        line = data_lines[i].strip()
        if line == "//":
            break
        if line and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 23:  # Match emission + transitions line
                state_index = len(emissions) + 1
                emissions[f"M{state_index}"] = parts[:20]
                i += 1
                transitions.append((f"M{state_index}", parts[20:26]))
            elif len(parts) == 23:  # Insert state emissions
                emissions[f"I{state_index}"] = parts[:20]
        i += 1

    return emissions, transitions

def visualize_hmm(emissions, transitions, output_file="hmm_full_graph"):
    g = Digraph(format='png', engine='dot')
    
    num_states = len([s for s in emissions if s.startswith("M")])
    
    for i in range(1, num_states + 1):
        g.node(f"M{i}", f"M{i}", shape="box", color="black")
        g.node(f"I{i}", f"I{i}", shape="ellipse", style="dashed", color="blue")
        g.node(f"D{i}", f"D{i}", shape="circle", style="dotted", color="red")
        
        if i < num_states:
            g.edge(f"M{i}", f"M{i+1}", label="M→M")
            g.edge(f"M{i}", f"I{i}", label="M→I")
            g.edge(f"M{i}", f"D{i+1}", label="M→D")
            g.edge(f"I{i}", f"M{i+1}", label="I→M")
            g.edge(f"D{i}", f"M{i+1}", label="D→M")

    g.render(output_file, view=False)

def print_emission_matrices(emissions):
    print("\n🔬 Emisijske vjerojatnosti:")
    for state, values in emissions.items():
        vstr = " ".join(f"{v:>5}" for v in values)
        print(f"{state}: {vstr}")

def open_image_safe(path):
    for opener in ["feh", "ristretto", "display", "xdg-open"]:
        if shutil.which(opener):
            subprocess.run([opener, path])
            return
    print("⚠️ Nije pronađen alat za otvaranje slike. Otvori ručno:", path)

def open_image_with_feh(path):
    if shutil.which("feh"):
        subprocess.run(["feh", "--auto-zoom", "--scale-down", "--geometry", "100%x100%", path])
    else:
        print("⚠️ 'feh' nije pronađen. Instaliraj ga s: sudo apt install feh")

# === GLAVNI POZIV ===
stockholm_file = "PF00001.alignment.seed"
hmm_file = "PF00001.hmm"

build_hmm(stockholm_file, hmm_file)
emissions, transitions = parse_hmm(hmm_file)
visualize_hmm(emissions, transitions)
print_emission_matrices(emissions)


# Poziv
open_image_with_feh("hmm_full_graph.png")