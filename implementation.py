import math
from collections import defaultdict, Counter
from graphviz import Digraph

alignment_matrix = []
sequence_order = []
seq_cons_line = None

with open("PF00001.alignment.seed", "r") as f:
    for line in f:
        line = line.rstrip()
        if line.startswith("#=GC seq_cons"):
            seq_cons_line = line.split(None, 2)[-1]
        elif line.startswith("#") or line.strip() == "" or line.startswith("//"):
            continue
        else:
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                sequence_order.append(parts[0])
                alignment_matrix.append(list(parts[1]))

if seq_cons_line:
    sequence_order.append("seq_cons")
    alignment_matrix.append(list(seq_cons_line))

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

num_seqs = len(alignment_matrix) - 1
num_cols = len(alignment_matrix[0])

match_cols = [i for i, c in enumerate(alignment_matrix[-1]) if c.isupper()]

emission_counts = defaultdict(lambda: Counter())
transitions = defaultdict(lambda: Counter())

for seq in alignment_matrix[:-1]:
    prev_state = "S"
    for i in range(num_cols):
        aa = seq[i]
        if i in match_cols:
            if aa == "-":
                curr_state = "D"
            else:
                curr_state = "M"
                emission_counts[f"M{i+1}"][aa] += 1
        else:
            if aa == "-":
                continue
            curr_state = "I"
            emission_counts[f"I{i+1}"][aa] += 1
        transitions[prev_state][curr_state] += 1
        prev_state = curr_state
    transitions[prev_state]["E"] += 1

# Realistična pozadinska distribucija aminokiselina (npr. UniProt frekvencije)
bg_freq = {
    'A': 0.074, 'C': 0.025, 'D': 0.054, 'E': 0.054, 'F': 0.047,
    'G': 0.074, 'H': 0.026, 'I': 0.068, 'K': 0.058, 'L': 0.099,
    'M': 0.025, 'N': 0.045, 'P': 0.039, 'Q': 0.034, 'R': 0.052,
    'S': 0.057, 'T': 0.051, 'V': 0.073, 'W': 0.013, 'Y': 0.032
}

alpha = 1  # pseudobroj (Laplace smoothing)

emissions = {}
log_odds = {}
for state, counts in emission_counts.items():
    total_counts = sum(counts.values())
    if total_counts == 0:
        continue
    K = len(amino_acids)
    emissions[state] = {}
    log_odds[state] = {}
    for aa in amino_acids:
        count = counts[aa] if aa in counts else 0
        p_obs = (count + alpha) / (total_counts + alpha * K)
        emissions[state][aa] = p_obs
        p_bg = bg_freq[aa]
        log_odds[state][aa] = math.log2(p_obs / p_bg)

# Prijelazne vjerojatnosti s pseudobrojevima
all_states = set(transitions.keys()) | {to for d in transitions.values() for to in d.keys()}
transition_probs = {}
for from_state in all_states:
    to_counts = transitions.get(from_state, {})
    total_counts = sum(to_counts.values())
    n_to_states = len(all_states)
    transition_probs[from_state] = {}
    for to_state in all_states:
        count = to_counts.get(to_state, 0)
        if total_counts > 0:
            p = (count + alpha) / (total_counts + alpha * n_to_states)
        else:
            p = 1 / n_to_states
        transition_probs[from_state][to_state] = p

# Export u HMMER-style format
with open("custom_model.hmm", "w") as f:
    f.write("HMMER3/f [Python custom model export]\n")
    f.write("NAME  CustomHMM\n")
    f.write("LENG  {}\n".format(len(match_cols)))
    f.write("ALPH  Amino\n")
    f.write("\nHMM        ")
    f.write(" ".join(amino_acids) + "\n")

    for i in range(1, len(match_cols)+1):
        state = f"M{i}"
        if state in log_odds:
            probs = ["{:.2f}".format(log_odds[state].get(aa, -999.0)) for aa in amino_acids]
        else:
            probs = ["*" for _ in amino_acids]
        f.write("{:5} {}\n".format(i, " ".join(probs)))
        trans = transition_probs.get(state, {})
        f.write("         {}\n".format(" ".join(f"{trans.get(t, 0):.2f}" for t in ["M", "I", "D"])))
    f.write("//\n")

# Vizualizacija HMM
g = Digraph(format='png', engine='dot')
g.attr(rankdir='LR')
g.attr(size='30,10')
g.attr(dpi='600')
g.graph_attr.update(splines="true", nodesep="0.6", ranksep="0.8")

for i in range(1, len(match_cols)+2):
    g.node(f"M{i}", shape='box', style='filled', fillcolor='lightgray')
    g.node(f"I{i}", shape='ellipse', style='filled', fillcolor='lightblue')
    g.node(f"D{i}", shape='circle', style='filled', fillcolor='lightpink')

g.node("S", shape='plaintext')
g.node("E", shape='plaintext')

for from_state, to_dict in transition_probs.items():
    for to_state, prob in to_dict.items():
        g.edge(from_state, to_state, label=f"{prob:.2f}")

g.render("hmm_full_graph", view=False)

# Ispis poravnanja
# for name, row in zip(sequence_order, alignment_matrix):
#     print(f"{name:20} {''.join(row)}")

# print(f"Dimenzije matrice: {len(alignment_matrix)} x {len(alignment_matrix[0])}")
# print(f"Element [0][0]: {alignment_matrix[0][0]}")
# print(f"Element [n-1][m-1]: {alignment_matrix[-1][-1]}")
# print(f"Element [0][1]: {alignment_matrix[0][1]}")
# print(f"Element [n-1][m-2]: {alignment_matrix[-1][-2]}")
# print(f"Element [0][11]: {alignment_matrix[0][11]}")
