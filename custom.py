import math
from collections import defaultdict, Counter
from graphviz import Digraph

# Učitavanje matrice poravnanja
alignment_matrix = []
sequence_order = []
seq_cons_line = None

def apply_laplace_smoothing(transition_counts, all_states, alpha=1):
    smoothed_transitions = defaultdict(lambda: defaultdict(int))

    # 1. Kopiraj postojeće prijelaze
    for from_state, to_dict in transition_counts.items():
        for to_state, count in to_dict.items():
            smoothed_transitions[from_state][to_state] += count

    # 2. Dodaj pseudocounts (Laplace) za moguće prijelaze
    for from_state in all_states:
        for to_state in all_states:
            if from_state == 'S' and to_state.startswith(('M1', 'I0', 'D1')):
                smoothed_transitions[from_state][to_state] += alpha

            elif from_state == 'E':
                continue  # ne dodaj ništa iz E

            elif from_state.startswith('M'):
                idx = int(from_state[1:])
                if to_state in {f'M{idx+1}', f'D{idx+1}', f'I{idx}'}:
                    smoothed_transitions[from_state][to_state] += alpha

            elif from_state.startswith('I'):
                idx = int(from_state[1:])
                if to_state in {f'I{idx}', f'D{idx+1}', f'M{idx+1}'}:
                    smoothed_transitions[from_state][to_state] += alpha

            elif from_state.startswith('D'):
                idx = int(from_state[1:])
                if to_state in {f'D{idx+1}', f'I{idx}', f'M{idx+1}'}:
                    smoothed_transitions[from_state][to_state] += alpha

    # 3. Normalizacija
    transition_probs_smoothed = {
        from_state: {
            to_state: count / sum(to_dict.values())
            for to_state, count in to_dict.items()
        }
        for from_state, to_dict in smoothed_transitions.items()
    }

    return transition_probs_smoothed




def export_hmm_to_r(transition_probs, emission_probs, symbols, filename="hmm_for_r.txt"):
    all_states = sorted(set(transition_probs.keys()) | {s for to_dict in transition_probs.values() for s in to_dict})
    state_indices = {state: idx for idx, state in enumerate(all_states)}

    # Generiraj transition matricu
    num_states = len(all_states)
    trans_matrix = [[0.0] * num_states for _ in range(num_states)]
    for from_state, to_dict in transition_probs.items():
        for to_state, prob in to_dict.items():
            i = state_indices[from_state]
            j = state_indices[to_state]
            trans_matrix[i][j] = prob

    # Generiraj emission matricu
    num_symbols = len(symbols)
    emis_matrix = [[0.0] * num_symbols for _ in range(num_states)]
    for state, emis_dict in emission_probs.items():
        i = state_indices[state]
        for j, sym in enumerate(symbols):
            emis_matrix[i][j] = emis_dict.get(sym, 0.0)

    with open(filename, "w") as f:
        f.write("# Definiraj stanja i simbole\n")
        f.write('states <- c({})\n'.format(", ".join(f'"{s}"' for s in all_states)))
        f.write('symbols <- c({})\n\n'.format(", ".join(f'"{s}"' for s in symbols)))

        f.write("# Prijelazna matrica\n")
        f.write("transProbs <- matrix(c(\n")
        for i, row in enumerate(trans_matrix):
            line = "  " + ", ".join(f"{x:.3f}" for x in row)
            if i != len(trans_matrix) - 1:
                line += ","
            line += "\n"
            f.write(line)
        f.write("), byrow = TRUE, nrow = %d)\n\n" % num_states)

        f.write("# Emisijska matrica\n")
        f.write("emissionProbs <- matrix(c(\n")
        for i, row in enumerate(emis_matrix):
            line = "  " + ", ".join(f"{x:.3f}" for x in row)
            if i != len(emis_matrix) - 1:
                line += ","
            line += "\n"
            f.write(line)
        f.write("), byrow = TRUE, nrow = %d)\n\n" % num_states)

        # Definiranje silent_states automatski:
        silent_states = [s for s in all_states if (s == "S" or s == "E" or s.startswith("D"))]
        f.write("# Silent states (ne emitiraju simbole)\n")
        f.write("silent_states <- c({})\n\n".format(", ".join(f'"{s}"' for s in silent_states)))


    #     f.write("library(HMM)\n")
    #     f.write("hmm <- initHMM(States=states, Symbols=symbols,\n")
        f.write("startProbs = c(" + ", ".join("1" if s == "S" else "0" for s in all_states) + ")\n")
    #     f.write("               transProbs = transProbs,\n")
    #     f.write("               emissionProbs = emissionProbs)\n\n")
    #    # f.write("set.seed(123)\n")
    #     f.write('gen <- simulate_until_end(hmm, start_state = "S", end_states = c("E"))\n')
    #     f.write("cat(\"Stanja:\\n\")\n")
    #     f.write("print(gen$states)\n")
    #     f.write("cat(\"Simboli:\\n\")\n")
    #     f.write("cat(paste(gen$observations[!is.na(gen$observations) & gen$observations != \"\"], collapse = \"\"), \"\\n\")\n")



# SVI simboli (aminokiseline)
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

def correct_insert_indices(temp):
    corrected = temp.copy()
    for i, state in enumerate(temp):
        if state.startswith('M'):
            m_idx = int(state[1:])
            for j in range(i+1, len(temp)):
                if corrected[j].startswith('I'):
                    i_idx = int(corrected[j][1:])
                    if i_idx < m_idx:
                        corrected[j] = f'I{m_idx}'
    return corrected

def visualize_hmm(transition_probs, emission_probs, filename="hmm_profile"):
    dot = Digraph(format="png")
    dot.attr(rankdir='LR')  # Lijevo na desno

    # Prikupi sve prijelaze s njihovim vjerojatnostima da bismo izračunali skaliranje debljine
    all_probs = [prob for to_dict in transition_probs.values() for prob in to_dict.values()]
    max_prob = max(all_probs) if all_probs else 1.0

    # Dodaj čvorove
    all_states = set(transition_probs.keys()) | {s for to_dict in transition_probs.values() for s in to_dict}
    for state in all_states:
        shape = "circle"
        if state == "S":
            shape = "doublecircle"
            #neka se to stanje osjenča u sivo
            dot.node(state, label=state, shape=shape, style="filled", fillcolor="lightgray")
        elif state == "E":
            shape = "doublecircle"
            #neka se to stanje osjenča u sivo
            dot.node(state, label=state, shape=shape, style="filled", fillcolor="lightgray")
        elif state.startswith("D"):
            shape = "circle"
            dot.node(state, label=state, shape=shape, style="filled", fillcolor="lightgray")
        elif state.startswith("I"):
            shape = "diamond"
        elif state.startswith("M"):
            shape = "box"

        label = state
        if state in emission_probs:
            emissions = emission_probs[state]
            top_emissions = sorted(emissions.items(), key=lambda x: -x[1])[:3]  # top 3
            emis_str = "\n".join(f"{aa}:{p:.2f}" for aa, p in top_emissions)
            label = f"{state}\n{emis_str}"

        dot.node(state, label=label, shape=shape)

    # Dodaj prijelaze s debljinom prema vjerojatnosti
    for from_state, to_dict in transition_probs.items():
        for to_state, prob in to_dict.items():
            # skaliraj penwidth između 1.0 i 5.0
            penwidth = 1.0 + (prob / max_prob) * 4.0
            dot.edge(from_state, to_state,
                     label=f"{prob:.2f}",
                     penwidth=str(penwidth))

    # Spremi i prikaži
    dot.render(filename, view=False)

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

# Dodaj zadnji red (seq_cons)
if seq_cons_line:
    sequence_order.append("seq_cons")
    alignment_matrix.append(list(seq_cons_line))

num_seqs = len(alignment_matrix) - 1  # bez seq_cons
num_cols = len(alignment_matrix[0])

seq_cons = alignment_matrix[-1]
#----------------------------------------------------
#smanjena verzija matrice
# Zadnjih 5 redaka poravnanja (bez seq_cons)
# alignment_matrix = alignment_matrix[-6:-1]  # posljednjih 5 sekvenci

# # Smanjimo ih na prvih 20 stupaca
# alignment_matrix = [row[:5] for row in alignment_matrix]

# # Skrati i seq_cons redak na prvih 20 znakova
# seq_cons = list(seq_cons[:5])


# #---------------------------------------------------

# # Ažuriraj osnovne varijable
# num_seqs = len(alignment_matrix)
# num_cols = len(alignment_matrix[0])
#---------------------------------------------------

# Prag za match stanje: broj nepraznih (ne '-' ili '.') > 50%
from collections import defaultdict

alignment_matrix = [
    list(".GAT.T..."),
    list("T...GTG.."),
    list("TG.TG.GAA"),
    list("..A.GC..."),
    list("....GTC..")
]

num_seqs = len(alignment_matrix)
num_cols = len(alignment_matrix[0])

def is_match_column(col_idx):
    threshold = 0.5 * num_seqs
    count = sum(1 for row in alignment_matrix if row[col_idx].isalpha())
    return count >= threshold

# 1. Generiraj temp listu: I* za svaki '' prije/između/nakon *, M* za svaki '*'
indeksi_konz = ['*' if is_match_column(i) else '' for i in range(num_cols)]

temp = []
insert_idx = 0
match_idx = 1

for i, znak in enumerate(indeksi_konz):
    if znak == '*':
        temp.append(f'M{match_idx}')
        match_idx += 1
    else:
        if i > 0 and temp[-1].startswith('I'):
            temp.append(temp[-1])
        else:
            temp.append(f'I{insert_idx}')
            insert_idx += 1

temp = correct_insert_indices(temp)
print("Korigirana temp lista stanja po stupcu:")
print(temp)

# 2. Inicijaliziraj brojače
transition_counts = defaultdict(lambda: defaultdict(int))
emission_counts   = defaultdict(lambda: defaultdict(int))

# 3. Prođi po svakom retku i grade prijelaze/emisije
for r, seq in enumerate(alignment_matrix):
    traj = ['S']
    prev = 'S'
    for c, sym in enumerate(seq):
        state_template = temp[c]
        if sym.isalpha():
            # aminokiselina -> ako temp[c] I* ili M*, prijelaz i emisija
            curr = state_template
            transition_counts[prev][curr] += 1
            emission_counts[curr][sym.upper()] += 1
            prev = curr
            traj.append(curr)

        else:  # sym == '.' ili '-'
            if state_template.startswith('M'):
                # delete stanje s istim indeksom
                idx = state_template[1:]
                curr = f'D{idx}'
                transition_counts[prev][curr] += 1
                prev = curr
                traj.append(curr)
            # ako temp[c] I*, ne radimo ništa
            # (niti prijelaz, niti emisiju)

    # završni prijelaz u E
    transition_counts[prev]['E'] += 1
    traj.append('E')

    print(f"Lista prijelaza reda {r}: {traj}")

# 4. Normaliziraj u vjerojatnosti
transition_probs = {
    frm: {to: cnt / sum(tots.values())
          for to, cnt in tots.items()}
    for frm, tots in transition_counts.items()
}

emission_probs = {
    st: {aa: cnt / sum(ems.values())
         for aa, cnt in ems.items()}
    for st, ems in emission_counts.items()
}

print("\nTransition probabilities:")
for frm, tos in transition_probs.items():
    print(f"{frm}: {tos}")

print("\nEmission probabilities:")
for st, ems in emission_probs.items():
    print(f"{st}: {ems}")

#poziv fje za zaglađivanje
all_states = set(transition_counts.keys()) | {s for to_dict in transition_counts.values() for s in to_dict}
transition_probs = apply_laplace_smoothing(transition_counts, all_states, 0.2)

# Pozovi funkciju za vizualizaciju
#visualize_hmm(transition_probs, emission_probs) #ovo može zaglaviti jer je prevelika slika

export_hmm_to_r(transition_probs, emission_probs, amino_acids)