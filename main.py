import os
import math
from collections import defaultdict, Counter
from graphviz import Digraph

seed_folder = "seed_folder"
output_filename = "hmm_for_r.txt"

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

def apply_laplace_smoothing(transition_counts, all_states, alpha=1):
    smoothed_transitions = defaultdict(lambda: defaultdict(int))
    for from_state, to_dict in transition_counts.items():
        for to_state, count in to_dict.items():
            smoothed_transitions[from_state][to_state] += count
    for from_state in all_states:
        for to_state in all_states:
            if from_state == 'S' and to_state.startswith(('M1', 'I0', 'D1')):
                smoothed_transitions[from_state][to_state] += alpha
            elif from_state == 'E':
                continue
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
    transition_probs_smoothed = {
        from_state: {
            to_state: count / sum(to_dict.values())
            for to_state, count in to_dict.items()
        }
        for from_state, to_dict in smoothed_transitions.items()
    }
    return transition_probs_smoothed

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

def export_hmm_to_r(fh, transition_probs, emission_probs, symbols, domain_name):
    all_states = sorted(set(transition_probs.keys()) | {s for to_dict in transition_probs.values() for s in to_dict})
    state_indices = {state: idx for idx, state in enumerate(all_states)}

    num_states = len(all_states)
    trans_matrix = [[0.0] * num_states for _ in range(num_states)]
    for from_state, to_dict in transition_probs.items():
        for to_state, prob in to_dict.items():
            i = state_indices[from_state]
            j = state_indices[to_state]
            trans_matrix[i][j] = prob

    num_symbols = len(symbols)
    emis_matrix = [[0.0] * num_symbols for _ in range(num_states)]
    for state, emis_dict in emission_probs.items():
        i = state_indices[state]
        for j, sym in enumerate(symbols):
            emis_matrix[i][j] = emis_dict.get(sym, 0.0)

    fh.write(f"# === DOMAIN {domain_name} ===\n")
    fh.write('states_{} <- c({})\n'.format(domain_name, ", ".join(f'"{s}"' for s in all_states)))
    fh.write('symbols_{} <- c({})\n'.format(domain_name, ", ".join(f'"{s}"' for s in symbols)))

    fh.write(f"transProbs_{domain_name} <- matrix(c(\n")
    for i, row in enumerate(trans_matrix):
        line = "  " + ", ".join(f"{x:.3f}" for x in row)
        if i != len(trans_matrix) - 1:
            line += ","
        line += "\n"
        fh.write(line)
    fh.write("), byrow = TRUE, nrow = %d)\n\n" % num_states)

    fh.write(f"emissionProbs_{domain_name} <- matrix(c(\n")
    for i, row in enumerate(emis_matrix):
        line = "  " + ", ".join(f"{x:.3f}" for x in row)
        if i != len(emis_matrix) - 1:
            line += ","
        line += "\n"
        fh.write(line)
    fh.write("), byrow = TRUE, nrow = %d)\n\n" % num_states)

    silent_states = [s for s in all_states if (s == "S" or s == "E" or s.startswith("D"))]
    fh.write("silent_states_{} <- c({})\n\n".format(domain_name, ", ".join(f'"{s}"' for s in silent_states)))
    fh.write("startProbs_{} <- c({})\n\n".format(
        domain_name,
        ", ".join("1" if s == "S" else "0" for s in all_states)
    ))

# Glavna petlja kroz sve datoteke
with open(output_filename, "w") as output_file:
    for filename in os.listdir(seed_folder):
        if not filename.endswith(".alignment.seed"):
            continue

        domain_name = filename.split(".")[0]
        filepath = os.path.join(seed_folder, filename)

        alignment_matrix = []
        sequence_order = []
        seq_cons_line = None

        with open(filepath, "r") as f:
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

        num_seqs = len(alignment_matrix) - 1
        num_cols = len(alignment_matrix[0])
        seq_cons = alignment_matrix[-1]

        def is_match_column(col_idx):
            threshold = 0.5 * num_seqs
            count = sum(1 for row in alignment_matrix[:-1] if row[col_idx].isalpha())
            return count >= threshold

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

        transition_counts = defaultdict(lambda: defaultdict(int))
        emission_counts = defaultdict(lambda: defaultdict(int))

        for seq in alignment_matrix[:-1]:  # bez seq_cons
            traj = ['S']
            prev = 'S'
            for c, sym in enumerate(seq):
                state_template = temp[c]
                if sym.isalpha():
                    curr = state_template
                    transition_counts[prev][curr] += 1
                    emission_counts[curr][sym.upper()] += 1
                    prev = curr
                    traj.append(curr)
                else:
                    if state_template.startswith('M'):
                        idx = state_template[1:]
                        curr = f'D{idx}'
                        transition_counts[prev][curr] += 1
                        prev = curr
                        traj.append(curr)
            transition_counts[prev]['E'] += 1
            traj.append('E')

        all_states = set(transition_counts.keys()) | {s for to_dict in transition_counts.values() for s in to_dict}
        transition_probs = apply_laplace_smoothing(transition_counts, all_states, alpha=0.2)
        emission_probs = {
            st: {aa: cnt / sum(ems.values())
                 for aa, cnt in ems.items()}
            for st, ems in emission_counts.items()
        }

        export_hmm_to_r(output_file, transition_probs, emission_probs, amino_acids, domain_name)
