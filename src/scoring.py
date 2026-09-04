# vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

"""
Interface scoring functions for AlphaFold-Multimer predictions.

Pure, ChimeraX-data-driven functions (no Qt/GUI code) so the same implementations
can be reused by both the "Confidence Overview" tab in pae_analysis.py and, later,
a batch/headless command. Only metrics whose formula and reference implementation
have been verified against a published, citable source are included here - see
docs/usage.md and docs/acknowledgements.md for the full citations:

- pDockQ: Bryant, Pozzati & Elofsson, "Improved prediction of protein-protein
  interactions using AlphaFold2", Nat Commun 13:1265 (2022).
  doi:10.1038/s41467-022-28865-w
- LIS / cLIS / iLIS: Kim, Hu, Comjean, Rodiger, Mohr & Perrimon, "Enhanced
  Protein-Protein Interaction Discovery via AlphaFold-Multimer", bioRxiv (2024).
  doi:10.1101/2024.02.19.580970 (github.com/flyark/AFM-LIS)
- Buried interface area and hydrogen bonds use UCSF ChimeraX's own built-in
  chimerax.atomic.buried_area / chimerax.hbonds.find_hbonds (Pettersen et al.,
  Protein Sci. 30:70-82 (2021), doi:10.1002/pro.3943).
- ipSAE (d0chn variant only): Dunbrack, github.com/DunbrackLab/IPSAE (ipsae.py).
  The d0(L) formula is the standard TM-score d0, from Yang & Skolnick, PROTEINS:
  Structure, Function, and Bioinformatics 57:702-710 (2004).

pDockQ2, mpDockQ, ipSAE's adaptive d0dom/d0res variants, and salt bridges are
deliberately not implemented yet - their formulas/constants need dedicated
verification against their original sources before being trusted.
"""

import numpy as np


def _representative_atoms(residues):
    """CB atom (CA for glycine or a residue missing CB) for each residue that has one."""
    residues_out, atoms_out = [], []
    for r in residues:
        atom = r.find_atom("CB") or r.find_atom("CA")
        if atom is not None:
            residues_out.append(r)
            atoms_out.append(atom)
    return residues_out, atoms_out


def _distance_matrix(atoms1, atoms2):
    from chimerax.atomic import Atoms
    xyz1 = Atoms(atoms1).scene_coords
    xyz2 = Atoms(atoms2).scene_coords
    diff = xyz1[:, None, :] - xyz2[None, :, :]
    return np.linalg.norm(diff, axis=-1)


def compute_contact_residues(residues1, residues2, cutoff=8.0):
    """Residues (from both sides) with a CB atom (CA for glycine) within `cutoff`
    Angstrom of a representative atom on the other side - the same "interface
    contact" definition pDockQ uses. Returns (residues_set, n_contacts)."""
    res1, atoms1 = _representative_atoms(residues1)
    res2, atoms2 = _representative_atoms(residues2)
    if not atoms1 or not atoms2:
        return set(), 0

    distances = _distance_matrix(atoms1, atoms2)
    close = distances <= cutoff
    n_contacts = int(close.sum())
    if n_contacts == 0:
        return set(), 0

    rows_with_contact = np.where(close.any(axis=1))[0]
    cols_with_contact = np.where(close.any(axis=0))[0]
    interface_residues = {res1[i] for i in rows_with_contact} | {res2[j] for j in cols_with_contact}
    return interface_residues, n_contacts


def compute_pdockq(residues1, residues2, cutoff=8.0):
    """pDockQ (Bryant, Pozzati & Elofsson, Nat Commun 13:1265, 2022).

    x = mean(interface pLDDT) * log10(n_contacts + 1)
    pDockQ = 0.724 / (1 + exp(-0.052 * (x - 152.611))) + 0.018

    An interface contact is a pair of CB atoms (CA for glycine), one from each
    side, within `cutoff` Angstrom (8 A in the original paper). Interface pLDDT
    is read from the atoms' B-factor column, the standard AlphaFold convention.

    Returns (pdockq, n_contacts).
    """
    interface_residues, n_contacts = compute_contact_residues(residues1, residues2, cutoff)
    if n_contacts == 0:
        return 0.0, 0

    plddt_values = [r.principal_atom.bfactor for r in interface_residues if r.principal_atom is not None]
    if not plddt_values:
        return 0.0, n_contacts

    avg_plddt = float(np.mean(plddt_values))
    x = avg_plddt * np.log10(n_contacts + 1)
    pdockq = 0.724 / (1.0 + np.exp(-0.052 * (x - 152.611))) + 0.018
    return float(pdockq), n_contacts


def _prepare_pae_residues(residues, pae_obj, row_index):
    """(residue, matrix-row-index, representative-atom) for residues that have both."""
    prepared = []
    for r in residues:
        i = row_index.get(r)
        if i is None:
            continue
        atom = r.find_atom("CB") or r.find_atom("CA")
        if atom is None:
            continue
        prepared.append((r, i, atom))
    return prepared


def compute_lis_clis_ilis(residues1, residues2, pae_obj, pae_cutoff=12.0, distance_cutoff=8.0):
    """LIS / cLIS / iLIS (Kim et al., bioRxiv 2024, doi:10.1101/2024.02.19.580970).

    LIS: mean of (pae_cutoff - PAE_ij) / pae_cutoff over residue pairs (i in
         residues1, j in residues2) with PAE_ij <= pae_cutoff (0.0 if none qualify).
    cLIS: the same, but additionally restricted to pairs whose CB (CA for
          glycine) atoms are within `distance_cutoff` Angstrom of each other.
    iLIS: sqrt(LIS * cLIS).

    `pae_obj` is a chimerax.alphafold AlphaFoldPAE instance (`structure.alphafold_pae`).
    Returns (lis, clis, ilis).
    """
    row_index = {ra: i for i, ra in enumerate(pae_obj.row_residues_or_atoms())}
    prepared1 = _prepare_pae_residues(residues1, pae_obj, row_index)
    prepared2 = _prepare_pae_residues(residues2, pae_obj, row_index)
    if not prepared1 or not prepared2:
        return 0.0, 0.0, 0.0

    idx1 = [p[1] for p in prepared1]
    idx2 = [p[1] for p in prepared2]
    pae_sub = pae_obj.pae_matrix[np.ix_(idx1, idx2)]

    atoms1 = [p[2] for p in prepared1]
    atoms2 = [p[2] for p in prepared2]
    dist_sub = _distance_matrix(atoms1, atoms2)

    scaled = np.clip((pae_cutoff - pae_sub) / pae_cutoff, 0.0, None)
    confident = pae_sub <= pae_cutoff
    lis = float(scaled[confident].mean()) if confident.any() else 0.0

    close = dist_sub <= distance_cutoff
    combined = confident & close
    clis = float(scaled[combined].mean()) if combined.any() else 0.0

    ilis = float(np.sqrt(lis * clis)) if lis > 0 and clis > 0 else 0.0
    return lis, clis, ilis


def compute_clis_confident_residues(residues1, residues2, pae_obj, pae_cutoff=12.0, distance_cutoff=8.0):
    """Residues (from both sides) involved in at least one cLIS-qualifying pair:
    PAE <= pae_cutoff AND CB (CA for glycine) distance <= distance_cutoff - the
    same definition compute_lis_clis_ilis() uses for cLIS. Returns a set of residues."""
    row_index = {ra: i for i, ra in enumerate(pae_obj.row_residues_or_atoms())}
    prepared1 = _prepare_pae_residues(residues1, pae_obj, row_index)
    prepared2 = _prepare_pae_residues(residues2, pae_obj, row_index)
    if not prepared1 or not prepared2:
        return set()

    idx1 = [p[1] for p in prepared1]
    idx2 = [p[1] for p in prepared2]
    pae_sub = pae_obj.pae_matrix[np.ix_(idx1, idx2)]

    atoms1 = [p[2] for p in prepared1]
    atoms2 = [p[2] for p in prepared2]
    dist_sub = _distance_matrix(atoms1, atoms2)

    confident = (pae_sub <= pae_cutoff) & (dist_sub <= distance_cutoff)
    rows_with = np.where(confident.any(axis=1))[0]
    cols_with = np.where(confident.any(axis=0))[0]
    res1 = [p[0] for p in prepared1]
    res2 = [p[0] for p in prepared2]
    return {res1[i] for i in rows_with} | {res2[j] for j in cols_with}


def compute_ipsae_d0chn(residues1, residues2, pae_obj, pae_cutoff=10.0):
    """ipSAE_d0chn - the fixed chain-length variant of ipSAE (Dunbrack,
    github.com/DunbrackLab/IPSAE, ipsae.py).

    d0 uses the standard TM-score d0(L) formula (Yang & Skolnick, PROTEINS
    57:702-710, 2004), with L = number of residues in chain1 + chain2 - fixed
    per chain pair, unlike the adaptive d0dom/d0res variants (not implemented
    here). For each residue i on one side, the score is the mean of
    1/(1+(PAE_ij/d0)^2) over all j on the other side with PAE_ij < pae_cutoff;
    the result is the maximum of that row-wise mean over all i, computed in
    both directions and taking the larger of the two (exactly as in ipsae.py).

    Returns 0.0 if no residue pair has PAE < pae_cutoff in either direction.
    """
    row_index = {ra: i for i, ra in enumerate(pae_obj.row_residues_or_atoms())}
    prepared1 = _prepare_pae_residues(residues1, pae_obj, row_index)
    prepared2 = _prepare_pae_residues(residues2, pae_obj, row_index)
    if not prepared1 or not prepared2:
        return 0.0

    idx1 = [p[1] for p in prepared1]
    idx2 = [p[1] for p in prepared2]
    n0 = len(idx1) + len(idx2)
    d0 = max(1.0, 1.24 * (n0 - 15) ** (1.0 / 3.0) - 1.8) if n0 > 27 else 1.0
    pae_matrix = pae_obj.pae_matrix

    def directional(rows, cols):
        sub = pae_matrix[np.ix_(rows, cols)]
        best = 0.0
        for row in sub:
            mask = row < pae_cutoff
            if mask.any():
                best = max(best, float((1.0 / (1.0 + (row[mask] / d0) ** 2.0)).mean()))
        return best

    return max(directional(idx1, idx2), directional(idx2, idx1))


def compute_contact_pae(residues1, residues2, pae_obj, cutoff=8.0):
    """Mean real PAE value (Angstrom) from each contact residue to the
    residues within `cutoff` Angstrom (CB-CB, or CA for glycine - the same
    contact definition compute_contact_residues()/pDockQ use) on the other
    side.

    Unlike pDockQ (contact count + pLDDT) or LIS/cLIS/iLIS/ipSAE (all
    transform/normalize PAE into a 0-1-ish score), this reads the raw PAE_ij
    matrix entries directly - the same underlying numbers ChimeraX's own
    `alphafold contacts` command uses to color its pseudobonds - just
    averaged per residue (one number per node) instead of shown per pair
    (one number per edge), so it can be stored/displayed as a residue
    attribute the same way every other score here is. Uses only the
    residues1-aligned/residues2-scored matrix direction, matching ChimeraX's
    own `alphafold contacts /chain1 to /chain2` (no `flip`) convention -
    not symmetrized, since the PAE matrix itself isn't symmetric.

    Returns a dict {residue: mean_pae}, containing only residues that are
    actually in contact (empty dict if none, or if `pae_obj` is None).
    """
    row_index = {ra: i for i, ra in enumerate(pae_obj.row_residues_or_atoms())}
    prepared1 = _prepare_pae_residues(residues1, pae_obj, row_index)
    prepared2 = _prepare_pae_residues(residues2, pae_obj, row_index)
    if not prepared1 or not prepared2:
        return {}

    idx1 = [p[1] for p in prepared1]
    idx2 = [p[1] for p in prepared2]
    res1 = [p[0] for p in prepared1]
    res2 = [p[0] for p in prepared2]
    atoms1 = [p[2] for p in prepared1]
    atoms2 = [p[2] for p in prepared2]

    dist_sub = _distance_matrix(atoms1, atoms2)
    pae_sub = pae_obj.pae_matrix[np.ix_(idx1, idx2)]
    close = dist_sub <= cutoff

    result = {}
    for i, r in enumerate(res1):
        mask = close[i, :]
        if mask.any():
            result[r] = float(pae_sub[i, mask].mean())
    for j, r in enumerate(res2):
        mask = close[:, j]
        if mask.any():
            result[r] = float(pae_sub[mask, j].mean())
    return result


def compute_buried_area(residues1, residues2, probe_radius=1.4):
    """Buried solvent-accessible surface area between two residue sets, via
    ChimeraX's own chimerax.atomic.buried_area (Pettersen et al., Protein Sci.
    30:70-82, 2021)."""
    from chimerax.atomic import buried_area, Atoms
    atoms1 = Atoms([a for r in residues1 for a in r.atoms])
    atoms2 = Atoms([a for r in residues2 for a in r.atoms])
    if len(atoms1) == 0 or len(atoms2) == 0:
        return 0.0
    ba, _, _, _ = buried_area(atoms1, atoms2, probe_radius)
    return float(ba)


def compute_hbond_count(session, model, chain_id1, chain_id2):
    """Number of hydrogen bonds between two chain groups of one structure, via
    ChimeraX's own chimerax.hbonds.find_hbonds (criteria from Mills & Dean,
    J. Comput.-Aided Mol. Des. 10:607-622, 1996 - the reference ChimeraX's own
    hbonds command cites). chain_id1/chain_id2 may each be a single chain ID
    string (the original, single-pair usage) or an iterable of chain IDs -
    lets a "protein" made of several chains (e.g. a Batch Analysis "first N
    chains vs. last chain" grouping) be treated as one side."""
    from chimerax.hbonds import find_hbonds
    group1 = {chain_id1} if isinstance(chain_id1, str) else set(chain_id1)
    group2 = {chain_id2} if isinstance(chain_id2, str) else set(chain_id2)
    pairs = find_hbonds(session, [model], inter_model=False, intra_model=True, status=False)
    count = 0
    for donor, acceptor in pairs:
        dc, ac = donor.residue.chain_id, acceptor.residue.chain_id
        if (dc in group1 and ac in group2) or (dc in group2 and ac in group1):
            count += 1
    return count
