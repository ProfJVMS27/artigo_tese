"""Frozen MCNP6.3 deck builder for Phase-5 validation of the canonical RDD case.

The deck matches the accepted MCNP6.3 campaign: Co-60 source at x=241 cm,
grenade shell and COMP-B cells, FM4=39442, 100 million histories, and explicit
normalization of M6 and M8 fractions.  It creates inputs only; it never runs
MCNP itself.
"""

from __future__ import annotations

import math
from collections.abc import Iterable


N_LAYERS = 7
INNER_RADIUS_CM = 35.0
OUTER_RADIUS_CM = 41.0
TOTAL_THICKNESS_CM = 6.0
MIN_THICKNESS_CM = 0.2
MAX_THICKNESS_CM = 2.0
MATERIALS = {"SiC", "UHMWPE", "Al", "Fe", "B4C-HDPE", "W-Epoxy", "Bi2O3-HDPE"}
MATERIAL_ID_DENSITY = {
    "SiC": (1, 3.2100), "Al": (2, 2.7000), "Fe": (3, 7.8700),
    "Bi2O3-HDPE": (4, 2.5570), "W-Epoxy": (9, 4.6430),
    "B4C-HDPE": (10, 1.4280), "UHMWPE": (11, 0.9600),
}
MATERIAL_CARDS = """M1 6000 0.300200000 $ SiC
       14000 0.699800000
M2 13000 1.000000000 $ Al
M3 26000 1.000000000 $ Fe
M4 83000 0.627900000 $ Bi2O3-HDPE
       8000 0.072100000
       6000 0.256900000
       1001 0.043100000
M5 27060 1.000000000 $ Co60_Source
M6 26000 0.894089409 $ Grenade_Shell
       30000 0.086508651
       13000 0.008400840
       25000 0.003600360
       82000 0.003300330
       24000 0.002500250
       28000 0.000900090
       23000 0.000700070
M7 6000 0.239800000 $ COMP_B
       1001 0.024480000
       7000 0.274950000
       8000 0.460770000
M8 7000 0.781121204 $ Air
       8000 0.209535433
       18000 0.009343363
M9 74000 0.800000000 $ W-Epoxy
       8000 0.040000000
       6000 0.140000000
       1001 0.020000000
M10 5000 0.414500000 $ B4C-HDPE
       6000 0.517900000
       1001 0.067600000
M11 6000 0.856280000 $ UHMWPE
       1001 0.143720000"""


def validate_design_row(row) -> tuple[list[str], list[float]]:
    mats = [str(row[f"material_camada_{index}"]).strip() for index in range(1, N_LAYERS + 1)]
    ths = [float(row[f"espessura_camada_{index}"]) for index in range(1, N_LAYERS + 1)]
    if set(mats) != MATERIALS or len(set(mats)) != N_LAYERS:
        raise ValueError("A validação exige as sete camadas, com um material de cada tipo.")
    if any(not math.isfinite(value) or value < MIN_THICKNESS_CM or value > MAX_THICKNESS_CM for value in ths):
        raise ValueError("Espessura fora dos limites [0,2; 2,0] cm.")
    if not math.isclose(sum(ths), TOTAL_THICKNESS_CM, rel_tol=0.0, abs_tol=1e-6):
        raise ValueError("A soma das espessuras deve ser exatamente 6,0 cm.")
    return mats, ths


def _radii(thicknesses: Iterable[float]) -> list[float]:
    radii = [INNER_RADIUS_CM]
    for thickness in thicknesses:
        radii.append(radii[-1] + float(thickness))
    if not math.isclose(radii[-1], OUTER_RADIUS_CM, rel_tol=0.0, abs_tol=1e-6):
        raise ValueError("Raio externo incompatível com a espessura total de 6,0 cm.")
    return radii


def build_mcnp_input(sim_id: int, materials: list[str], thicknesses: list[float]) -> str:
    if len(materials) != N_LAYERS or len(thicknesses) != N_LAYERS:
        raise ValueError("São necessárias sete camadas.")
    radii = _radii(thicknesses)
    layer_lines = []
    # Surface 4 is the axial PZ boundary, not a radial layer interface.
    # Radial layers must therefore follow 1-2, 2-3, 3-5, 5-6, ..., 8-9.
    surface_pairs = [(1, 2), (2, 3), (3, 5), (5, 6), (6, 7), (7, 8), (8, 9)]
    for index, (material, (left, right)) in enumerate(zip(materials, surface_pairs), start=1):
        mat_id, density = MATERIAL_ID_DENSITY[material]
        layer_lines.append(
            f"{index} {mat_id} -{density:.4f} {left} -{right} 4 -10 11 imp:p=1 $ camada {index}: {material}"
        )
    return "\n".join([
        f"c sim_{sim_id:05d}", "c cells", *layer_lines,
        "8 5 -8.9000 -25 imp:p=1 $ fonte volumetrica Co-60",
        "12 6 -7.6810 -15 16 imp:p=1 $ casco da granada",
        "13 7 -1.6500 -16 imp:p=1 $ explosivo COMP B",
        "20 8 -0.0012 -21 imp:p=1 $ detector",
        "11 0 -30 #1 #2 #3 #4 #5 #6 #7 #8 #12 #13 #20 imp:p=1 $ interior",
        "30 0 30 imp:p=0 $ graveyard", "", "c surfaces",
        f"1 CZ {radii[0]:.3f}", f"2 CZ {radii[1]:.3f}", f"3 CZ {radii[2]:.3f}",
        "4 PZ 90.000", f"5 CZ {radii[3]:.3f}", f"6 CZ {radii[4]:.3f}",
        f"7 CZ {radii[5]:.3f}", f"8 CZ {radii[6]:.3f}", f"9 CZ {radii[7]:.3f}",
        "10 PZ 160.000", "11 PX 0", "15 S 241 10 125 4.5", "16 S 241 10 125 3.5",
        "21 RCC 0 0 120 0 0 10 5", "25 RCC 241 0 124 0 0 2 1", "30 S 125 0 125 150",
        "", "c materials", "MODE P", MATERIAL_CARDS, "c tally and source", "F4:P 20", "FM4 39442",
        "NPS 100000000", "SDEF POS=241 0 125 AXS=0 0 1 EXT=D1 RAD=D2 ERG=D3 PAR=2",
        "SI1 -1 1", "SP1 0 1", "SI2 0 1", "SP2 -21 1", "SI3 L 1.1732 1.3325", "SP3 1 1", "",
    ])


def validate_input_text(text: str) -> None:
    required = ("MODE P", "M6 26000 0.894089409", "M8 7000 0.781121204", "F4:P 20", "FM4 39442", "NPS 100000000", "SDEF POS=241")
    missing = [value for value in required if value not in text]
    if missing:
        raise ValueError(f"Deck MCNP6.3 canônico incompleto: {missing}")
    if "1001 0.143720000\n\nc tally and source" in text:
        raise ValueError("Linha em branco proibida entre os M-cards e o bloco F4/NPS/SDEF.")
    if "4 -6 4 -10" in text:
        raise ValueError("A superfície 4 é axial (PZ) e não pode iniciar uma camada radial.")
