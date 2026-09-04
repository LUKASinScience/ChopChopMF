# Acknowledgements

## :material-format-quote-close: Citation

If ChopChopMF was useful for your work, please cite the repository - see [`CITATION.cff`](https://github.com/LUKASinScience/ChopChopMF/blob/main/CITATION.cff) for the canonical, machine-readable entry (GitHub renders this as a "Cite this repository" button on the repo page).

## :material-account-group: Contributors

We Thank **Christoph Büschl** and **Mathias Percipalle** (Schur Group, Institute of Science and Technology Austria) for their insightful discussions and contributions toward the refinement of these tools. We are particularly grateful to **Tom Goddard** (University of California, San Francisco) for the release of UCSF ChimeraX version 1.10.dev202502272020, which addressed critical bugs within the Foldseek integration, and for his valuable feedback that helped shape the software's development. Further we want to thank **Pim Huis in ‘t Veld** (Huis Group, Max Perutz Labs), for discussing the idea of a ΔG Filter for the PDBePISA analysis.

We thank **Mariana Pereira Guarda** (de Bono Group, Institute of Science and Technology Austria) for her assistance in selecting the color schemes and for her conceptual ideas that directly influenced the implementation of new features. We also thank **Ana Chegão** for her suggestions regarding the nomenclature of the ChopChop plug-in and **Philipp Bauer** for his expert guidance in designing the intuitive graphical icons for the toolbar. Finally, we acknowledge the use of large language models, specifically ChatGPT-4, Gemini 1.5, and Claude (Anthropic), as computational copilots for scripting, bug identification, and feature development throughout the project.

## :material-flask-outline: Citations for Interface Scores

The "Interface Scores" reported in the **PAE Analysis** tool's Scores tab implement published, peer-reviewed methods. If you use them in your work, please cite the original sources:

- **pDockQ** — Bryant, P., Pozzati, G. & Elofsson, A. Improved prediction of protein-protein interactions using AlphaFold2. *Nat Commun* **13**, 1265 (2022). [doi:10.1038/s41467-022-28865-w](https://doi.org/10.1038/s41467-022-28865-w)
- **LIS / cLIS / iLIS** — Kim, A.-R., Hu, Y., Comjean, A., Rodiger, J., Mohr, S. E. & Perrimon, N. Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer. *bioRxiv* (2024). [doi:10.1101/2024.02.19.580970](https://doi.org/10.1101/2024.02.19.580970) ([github.com/flyark/AFM-LIS](https://github.com/flyark/AFM-LIS))
- **ipSAE (d0chn variant)** — Dunbrack, R. L. Jr. [github.com/DunbrackLab/IPSAE](https://github.com/DunbrackLab/IPSAE). The underlying d0(L) formula is the standard TM-score d0: Yang, Y. & Skolnick, J. Scoring function for automated assessment of protein structure template quality. *Proteins* **57**, 702–710 (2004). Only the fixed chain-length `d0chn` variant is implemented; the adaptive `d0dom`/`d0res` variants are not.
- **Buried interface area** and **hydrogen bond count** use UCSF ChimeraX's own built-in `measure buriedarea` and `hbonds` commands rather than a separate implementation — please also cite ChimeraX itself: Pettersen, E. F. *et al.* UCSF ChimeraX: Structure visualization for researchers, educators, and developers. *Protein Sci.* **30**, 70–82 (2021). [doi:10.1002/pro.3943](https://doi.org/10.1002/pro.3943)

## :material-flask-outline: Citations for Cell Biology Motifs

The **Phospho Sites**, **Signal Peptide**, and **TM Helix** tools flag residues using
hand-written regex approximations of published consensus motifs — a **heuristic
motif match**, not a calibrated prediction (unlike SignalP/TMHMM/NetPhos, which are
trained statistical predictors). Motifs marked *(approximate)* loosen the cited
consensus (a wider spacer, a dropped positional constraint) rather than reproducing
it exactly.

**Phospho Sites** (kinase consensus motifs, S/T/Y acceptor only):

- **ATM/ATR** `[ST]Q` — Kim, S.-T. *et al.* *J Biol Chem* **274**, 37538–37543 (1999).
- **PKA** `[RK]{2}.[ST]` *(approximate)* — Kemp, B. E. & Pearson, R. B. *Trends Biochem Sci* **15**, 342–346 (1990).
- **PKC** `[RK].{0,2}[ST].[RK]` *(approximate)* — Nishikawa, K. *et al.* *J Biol Chem* **272**, 952–960 (1997).
- **CDK** `[ST]P.[RK]` — Songyang, Z. *et al.* *Curr Biol* **4**, 973–982 (1994).
- **CK2** `[ST].{2}[DE]` — Pinna, L. A. *Biochim Biophys Acta* **1054**, 267–284 (1990).
- **GSK3** `[ST].{3}[ST]` *(approximate — omits the "primed", already-phosphorylated downstream residue the real consensus requires)* — Fiol, C. J. *et al.* *J Biol Chem* (1987); Roach, P. J. (1990).
- **MAPK** `P..[ST]P` *(approximate)* — Clark-Lewis, I. *et al.* *J Biol Chem* **266**, 15180–15184 (1991).
- **AKT** `R.R..[ST]` — Alessi, D. R. *et al.* *FEBS Lett* **399**, 333–338 (1996).
- **AMPK** `[MFLIV]R..[ST]...[MFLIV]` *(approximate)* — Gwinn, D. M. *et al.* *Mol Cell* **30**, 214–226 (2008).

**Signal Peptide** (subcellular targeting-signal motifs):

- **NLS, monopartite** `[RK]{4,6}` *(approximate)* — Dingwall, C. & Laskey, R. A. *Trends Biochem Sci* **16**, 478–481 (1991).
- **NLS, bipartite** `[RK]{2}.{10,12}[RK]{3}` — Robbins, J., Dilworth, S. M., Laskey, R. A. & Dingwall, C. *Cell* **64**, 615–623 (1991).
- **NES, leucine-rich** — Güttler, T. *et al.* *Nat Struct Mol Biol* **17**, 1367–1376 (2010).
- **ER retention (KDEL/HDEL)** — Munro, S. & Pelham, H. R. B. *Cell* **48**, 899–907 (1987).
- **Peroxisomal PTS1** `[SAC][KRH][LM]$` *(approximate)* — Gould, S. J. *et al.* *J Cell Biol* **108**, 1657–1664 (1989).
- **Peroxisomal PTS2** *(approximate)* — Swinkels, B. W. *et al.* *EMBO J* **10**, 3255–3262 (1991).
- **Mitochondrial presequence heuristic** (N-terminal positive/negative charge count) — this project's own threshold, inspired by the amphipathic, net-positively-charged presequence property described in von Heijne, G. *EMBO J* **5**, 1335–1342 (1986) and Roise, D. *et al.* *EMBO J* **5**, 1327–1334 (1986); the specific "≥4 positive, ≤1 negative in the first 30 residues" cutoff is not itself a published rule.

**TM Helix**: candidate helices are DSSP-detected alpha-helices scored by mean
Kyte & Doolittle (1982) hydrophobicity. When available, ChimeraX's own built-in `mlp`
lipophilicity table additionally gives a SASA-weighted per-helix score (the
[rbvi/chimerax-recipes "helixmlp"](https://rbvi.github.io/chimerax-recipes/helixmlp/helixmlp.html)
method) — please cite:

- **Kyte & Doolittle hydrophobicity scale** — Kyte, J. & Doolittle, R. F. *J Mol Biol* **157**, 105–132 (1982).
- **Fauchère–Pliska atomic lipophilicity (MLP)** — Fauchère, J.-L. & Pliska, V. *Eur J Med Chem* **18**, 369–375 (1983), as implemented in ChimeraX's `mlp` command.

**License note:** ChopChopMF is licensed under the GNU AGPLv3. None of the code above was copied from the cited sources — each formula was independently reimplemented in `scoring.py` from its published description. For the two GitHub repositories used as a reference ([DunbrackLab/IPSAE](https://github.com/DunbrackLab/IPSAE) and [flyark/AFM-LIS](https://github.com/flyark/AFM-LIS)), both are MIT-licensed, which is compatible with inclusion in an AGPLv3 project either way. The pDockQ formula comes from a Nature Communications article (fully open-access, CC BY 4.0); mathematical formulas and methods described in a paper are in any case not copyrightable subject matter (only the paper's specific text/figures are) — reimplementing them as original code is standard, unrestricted scientific practice, independent of the license of a preprint/article.
