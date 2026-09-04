#!/usr/bin/env python3

"""
Pure (no Qt, no ChimeraX run()) sequence-motif logic shared by the three
Cell Biology tools (phospho_sites.py, signal_peptide.py, transmembrane_helix.py):
- Kinase phosphorylation-site consensus motifs
- Subcellular targeting-signal motifs
- Kyte & Doolittle hydrophobicity scale for transmembrane-helix detection

Every motif below is a hand-written regex approximation of a published
consensus, not the output of a calibrated predictor (SignalP/TMHMM/NetPhos) -
results are "heuristic motif match", not a prediction confidence. Citations
verified against the literature; anything marked "approximate" loosens the
real consensus (shorter/longer spacer, dropped positional constraint) rather
than reproducing it exactly - see docs/acknowledgements.md for full sources.
"""

import re

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}

# Kyte & Doolittle 1982, J Mol Biol 157(1):105-132.
KYTE_DOOLITTLE = {
    'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
    'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
    'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLU': -3.5,
    'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5,
}

# (kinase, regex, citation, is_approximate) - S/T/Y acceptor residues within
# each match are the reported hits. "Approximate" = looser than the textbook
# consensus (see docs/acknowledgements.md for what's specifically loosened).
PHOSPHO_MOTIFS = [
    ("ATM/ATR", r"[ST]Q", "Kim et al. 1999, J Biol Chem 274:37538-43", False),
    ("PKA", r"[RK]{2}.[ST]", "Kemp & Pearson 1990, Trends Biochem Sci 15:342-6", True),
    ("PKC", r"[RK].{0,2}[ST].[RK]", "Nishikawa et al. 1997, J Biol Chem 272:952-60", True),
    ("CDK", r"[ST]P.[RK]", "Songyang et al. 1994, Curr Biol 4:973-82", False),
    ("CK2", r"[ST].{2}[DE]", "Pinna 1990, Biochim Biophys Acta 1054:267-84", False),
    ("GSK3", r"[ST].{3}[ST]", "Fiol et al. 1987, J Biol Chem; Roach 1990", True),
    ("MAPK", r"P..[ST]P", "Clark-Lewis et al. 1991, J Biol Chem 266:15180-4", True),
    ("AKT", r"R.R..[ST]", "Alessi et al. 1996, FEBS Lett 399:333-8", False),
    ("AMPK", r"[MFLIV]R..[ST]...[MFLIV]", "Gwinn et al. 2008, Mol Cell 30:214-26", True),
]

# (signal, regex, organelle, citation, is_approximate). Regexes ending in "$"
# are C-terminal only (matched against the trailing residues of the sequence).
SIGNAL_MOTIFS = [
    ("NLS (monopartite)", r"[RK]{4,6}", "Nucleus (import)",
     "Dingwall & Laskey 1991, Trends Biochem Sci 16:478-81", True),
    ("NLS (bipartite)", r"[RK]{2}.{10,12}[RK]{3}", "Nucleus (import)",
     "Robbins et al. 1991, Cell 64:615-23", False),
    ("NES (leucine-rich)", r"[LIVFM].{2,3}[LIVFM].{2,3}[LIVFM].[LIVFM]", "Nucleus (export)",
     "Guttler et al. 2010, Nat Struct Mol Biol 17:1367-76", False),
    ("ER retention (KDEL/HDEL)", r"(?:KDEL|HDEL)$", "ER (retention)",
     "Munro & Pelham 1987, Cell 48:899-907", False),
    ("Peroxisomal PTS1", r"[SAC][KRH][LM]$", "Peroxisome (import)",
     "Gould et al. 1989, J Cell Biol 108:1657-64", True),
    ("Peroxisomal PTS2", r"[RK][LVI].{5}[HQ][LA]", "Peroxisome (import)",
     "Swinkels et al. 1991, EMBO J 10:3255-62", True),
]

# Not a regex - N-terminal charge heuristic. Real basis is the amphipathic,
# net-positively-charged mitochondrial presequence (von Heijne 1986, EMBO J
# 5:1335-42; Roise et al. 1986, EMBO J 5:1327-34); the specific thresholds
# below are this project's own heuristic "inspired by" that property, not a
# verbatim published rule.
MTS_WINDOW = 30
MTS_MIN_POSITIVE = 4
MTS_MAX_NEGATIVE = 1
MTS_CITATION = "Heuristic inspired by von Heijne 1986, EMBO J 5:1335-42 (not a published rule)"


def sequence_and_residues(chain):
    """(one-letter sequence, matching list of Residue) for a chimerax.atomic
    Chain, skipping any non-standard-amino-acid residue."""
    seq, residues = [], []
    for res in chain.existing_residues:
        one = THREE_TO_ONE.get(res.name.upper())
        if one is not None:
            seq.append(one)
            residues.append(res)
    return "".join(seq), residues


def residue_range(start_res, end_res, all_residues):
    """All residues (from `all_residues`, in sequence order) between
    start_res and end_res inclusive - both must be members of all_residues."""
    started = False
    for res in all_residues:
        if res is start_res:
            started = True
        if started:
            yield res
        if res is end_res:
            break


def _residue_plddt(residue):
    values = [a.bfactor for a in residue.atoms]
    return sum(values) / len(values) if values else 100.0


def _residue_sasa(residue):
    return sum(getattr(a, 'area', 0.0) or 0.0 for a in residue.atoms)


def scan_phospho_sites(sequence, residues, plddt_threshold=70.0, sasa_threshold=5.0):
    """Return a list of hit dicts: residue, kinase, plddt, sasa, citation,
    approximate. Only S/T/Y residues within a matched motif, filtered to
    disordered (pLDDT < threshold) and surface-exposed (SASA > threshold) -
    the same "phosphorylation needs a flexible, accessible site" filter as
    the original prototype script."""
    hits = {}
    for kinase, pattern, citation, approx in PHOSPHO_MOTIFS:
        for match in re.finditer(pattern, sequence):
            for idx in range(match.start(), match.end()):
                res = residues[idx]
                if res.name.upper() not in ('SER', 'THR', 'TYR'):
                    continue
                plddt = _residue_plddt(res)
                sasa = _residue_sasa(res)
                if plddt < plddt_threshold and sasa > sasa_threshold:
                    entry = hits.setdefault(res, {
                        "residue": res, "plddt": plddt, "sasa": sasa, "kinases": [],
                    })
                    entry["kinases"].append((kinase, citation, approx))
    return sorted(hits.values(), key=lambda h: (h["residue"].chain_id, h["residue"].number))


def scan_signal_motifs(sequence, residues, sasa_threshold=5.0):
    """Return a list of hit dicts: name, organelle, citation, approximate,
    start_residue, end_residue, seq, plddt, sasa. C-terminal-anchored motifs
    (regex ending in "$") only match at the very end of the sequence."""
    hits = []
    for name, pattern, organelle, citation, approx in SIGNAL_MOTIFS:
        for match in re.finditer(pattern, sequence):
            idxs = range(match.start(), match.end())
            plddt_vals = [_residue_plddt(residues[i]) for i in idxs]
            sasa = sum(_residue_sasa(residues[i]) for i in idxs)
            if sasa <= sasa_threshold:
                continue
            hits.append({
                "name": name, "organelle": organelle, "citation": citation, "approximate": approx,
                "start_residue": residues[match.start()], "end_residue": residues[match.end() - 1],
                "seq": match.group(0),
                "plddt": sum(plddt_vals) / len(plddt_vals) if plddt_vals else 0.0,
                "sasa": sasa,
            })
    return hits


def scan_mitochondrial_presequence(sequence):
    """(is_hit, positive_count, negative_count) for the N-terminal MTS_WINDOW
    residues - see MTS_CITATION for what this is and isn't backed by."""
    n_term = sequence[:MTS_WINDOW]
    positive = n_term.count('R') + n_term.count('K')
    negative = n_term.count('D') + n_term.count('E')
    is_hit = positive >= MTS_MIN_POSITIVE and negative <= MTS_MAX_NEGATIVE
    return is_hit, positive, negative


def _mlp_atom_table():
    """Per-atom Fauchere & Pliska 1983 lipophilicity values (Biochem
    Pharmacol 32:2723-8), from ChimeraX's own built-in `mlp` command data -
    same table/method as the rbvi/chimerax-recipes "helixmlp" recipe (SASA-
    weighted per-helix lipophilicity to spot membrane-facing helices).
    Returns None if chimerax.mlp isn't importable (falls back to plain
    Kyte-Doolittle mean below)."""
    try:
        from chimerax.mlp.mlp import Defaults
        return Defaults().fidatadefault
    except Exception:
        return None


def scan_tm_helices(residues, min_length=15, threshold=1.0):
    """Return a list of hit dicts: start_residue, end_residue, length,
    hydrophobicity (mean Kyte & Doolittle 1982 per-residue score), and -
    when SASA has already been measured (atom.area set) and chimerax.mlp is
    available - mlp_score (SASA-weighted Fauchere & Pliska 1983 atomic
    lipophilicity, the citable "helixmlp" recipe method: membrane-facing
    surface should be lipophilic, not just the residue identity). Requires
    DSSP to have already been run on the model (residues carry `ss_type`);
    groups consecutive alpha-helix residues, keeps runs >= min_length whose
    mean Kyte-Doolittle score >= threshold."""
    runs = []
    current = []
    for res in residues:
        if getattr(res, 'ss_type', None) == 1:  # chimerax.atomic.Residue.SS_HELIX
            current.append(res)
        else:
            if len(current) >= min_length:
                runs.append(current)
            current = []
    if len(current) >= min_length:
        runs.append(current)

    mlp_table = _mlp_atom_table()
    results = []
    for helix in runs:
        scores = [KYTE_DOOLITTLE.get(r.name.upper(), 0.0) for r in helix]
        mean_score = sum(scores) / len(scores)
        if mean_score < threshold:
            continue
        hit = {
            "start_residue": helix[0], "end_residue": helix[-1],
            "length": len(helix), "hydrophobicity": mean_score, "mlp_score": None,
        }
        if mlp_table is not None:
            weighted_sum, total_area = 0.0, 0.0
            for res in helix:
                res_table = mlp_table.get(res.name.upper())
                if res_table is None:
                    continue
                for atom in res.atoms:
                    area = getattr(atom, 'area', None)
                    fi = res_table.get(atom.name)
                    if area and fi is not None:
                        weighted_sum += area * fi
                        total_area += area
            if total_area > 0:
                hit["mlp_score"] = weighted_sum / total_area
        results.append(hit)
    return results
