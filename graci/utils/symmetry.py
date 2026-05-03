"""
Point group symmetry correlation tables and utilities.

Provides correlation tables mapping full (non-Abelian) point group irreps to
their decompositions in Abelian computational subgroups, and helper functions
used by find_degen_groups to identify orbitals that are degenerate in the full
point group symmetry.
"""

# ---------------------------------------------------------------------------
# Correlation tables
# ---------------------------------------------------------------------------
# Format:
#   { full_group : { full_irrep : { subgroup : (comp_a,) or (comp_a, comp_b) } } }
#
# A tuple of length 1 means the full-group irrep is 1D in that subgroup (singleton).
# A tuple of length 2 means the full-group irrep is 2D and yields a degenerate pair.
# All labels lowercase.
#
# Currently implemented: Dooh (D∞h) and Coov (C∞v).
# Extend by adding entries for D3h, C3v, D4h, D6h, Oh, Td, etc. as needed.

CORRELATION_TABLES = {

    'dooh': {
        # --- 1D gerade ---
        'A1g': {'d2h': ('ag',),   'c2v': ('a1',)},   # Σg+
        'A2g': {'d2h': ('b1g',),  'c2v': ('a2',)},   # Σg-

        # --- 2D gerade ---
        'E1g': {'d2h': ('b2g', 'b3g'), 'c2v': ('b1', 'b2')},  # Πg
        'E2g': {'d2h': ('ag',  'b1g'), 'c2v': ('a1', 'a2')},  # Δg
        'E3g': {'d2h': ('b2g', 'b3g'), 'c2v': ('b1', 'b2')},  # Φg
        'E4g': {'d2h': ('ag',  'b1g'), 'c2v': ('a1', 'a2')},  # Γg

        # --- 1D ungerade ---
        'A1u': {'d2h': ('b1u',),  'c2v': ('a1',)},   # Σu+
        'A2u': {'d2h': ('au',),   'c2v': ('a2',)},   # Σu-

        # --- 2D ungerade ---
        'E1u': {'d2h': ('b2u', 'b3u'), 'c2v': ('b1', 'b2')},  # Πu
        'E2u': {'d2h': ('au',  'b1u'), 'c2v': ('a1', 'a2')},  # Δu
        'E3u': {'d2h': ('b2u', 'b3u'), 'c2v': ('b1', 'b2')},  # Φu
        'E4u': {'d2h': ('au',  'b1u'), 'c2v': ('a1', 'a2')},  # Γu
    },

    'coov': {
        # --- 1D ---
        'A1': {'c2v': ('a1',)},   # Σ+
        'A2': {'c2v': ('a2',)},   # Σ-

        # --- 2D ---
        'E1': {'c2v': ('b1', 'b2')},  # Π
        'E2': {'c2v': ('a1', 'a2')},  # Δ
        'E3': {'c2v': ('b1', 'b2')},  # Φ
        'E4': {'c2v': ('a1', 'a2')},  # Γ
    },

    # Future entries: 'd3h', 'c3v', 'd4h', 'd6h', 'oh', 'td', ...
}

# ---------------------------------------------------------------------------
# PySCF component label pairs for 2D irreps
# ---------------------------------------------------------------------------
# For Dooh and Coov, PySCF assigns each orbital from a 2D (E-type) irrep a
# component label with an 'x' or 'y' suffix (e.g., 'E1gx', 'E1gy' for Πg).
# These labels appear in scf.orb_irrep (populated via mol.irrep_name).
#
# Using these full labels to identify degenerate partners is more reliable
# than matching via Abelian subgroup irreps + energy, because each component
# label uniquely identifies the full-group irrep (e.g., 'E2gx' is unambiguously
# Δg, while its D2h counterpart 'ag' is shared with Σg+).
#
# Key: full-group 2D irrep name (matching CORRELATION_TABLES keys above).
# Value: (pyscf_label_a, pyscf_label_b) — the two component labels PySCF uses.
#
# NOTE: Coov labels ('E1x'/'E1y' etc.) follow the same x/y convention as Dooh
# but without the g/u suffix. Verify against an actual PySCF Coov calculation
# if unexpected behaviour is seen.

PYSCF_COMPONENTS = {
    'dooh': {
        'E1g': ('E1gx', 'E1gy'),
        'E2g': ('E2gx', 'E2gy'),
        'E3g': ('E3gx', 'E3gy'),
        'E4g': ('E4gx', 'E4gy'),
        'E1u': ('E1ux', 'E1uy'),
        'E2u': ('E2ux', 'E2uy'),
        'E3u': ('E3ux', 'E3uy'),
        'E4u': ('E4ux', 'E4uy'),
    },
    'coov': {
        'E1': ('E1x', 'E1y'),
        'E2': ('E2x', 'E2y'),
        'E3': ('E3x', 'E3y'),
        'E4': ('E4x', 'E4y'),
    },
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_degen_pairs(topgroup, compgroup):
    """Return potentially-degenerate subgroup irrep pairs from CORRELATION_TABLES.

    Queries the table for (topgroup, compgroup) and returns every frozenset of
    two subgroup irrep labels that arises from the same 2D full-group irrep.
    Duplicate pairs (e.g., E1g and E3g both yielding {b2g, b3g}) are
    deduplicated.

    Used as a fallback for groups not in PYSCF_COMPONENTS. Within each returned
    pair, energy matching is required to identify actual degenerate partners.

    Returns
    -------
    list of frozenset of str
        Each frozenset contains two lowercase subgroup irrep labels.
    """
    tbl   = CORRELATION_TABLES.get(topgroup.lower(), {})
    seen  = set()
    pairs = []
    for full_irrep, subgroup_map in tbl.items():
        comp = subgroup_map.get(compgroup.lower())
        if comp is not None and len(comp) == 2:
            key = frozenset(comp)
            if key not in seen:
                seen.add(key)
                pairs.append(key)
    return pairs


def get_pyscf_degen_label_pairs(topgroup):
    """Return PySCF orbital label pairs for 2D irreps in topgroup.

    For Dooh and Coov, PySCF assigns individual component labels (e.g.,
    'E1gx', 'E1gy') to orbitals from 2D full-group irreps.  These appear
    directly in scf.orb_irrep and allow unambiguous identification of
    degenerate partners without requiring energy matching against Abelian
    subgroup irrep labels.

    Returns
    -------
    list of (str, str)
        One tuple per 2D irrep: (label_a, label_b).
        Empty list if topgroup is not in PYSCF_COMPONENTS.
    """
    return list(PYSCF_COMPONENTS.get(topgroup.lower(), {}).values())
