ANIMATE_MAX_LEN = 100

NT_COLORS = {
    "A": "#E8C32A",   # yellow  (hue ~50deg, distinct from FN orange at ~28deg)
    "U": "#2D9B8A",   # teal    (unchanged)
    "G": "#4A6FA5",   # slate-blue (unchanged)
    "C": "#9856A5",   # purple  (hue ~285deg, distinct from FP red at ~8deg)
}

NT_UNKNOWN = "#888888"

BACKBONE_COLOR     = "#AAAAAA"
BACKBONE_COLOR_CMP = "#CCCCCC"
PAIR_COLOR         = "#444444"
ARC_COLOR          = "#555555"

# General chart palette
CHART_COLOR = "#4A6FA5"   # slate-blue primary accent (bars, single-series)
BOX_FILL    = "#B5C9E2"   # tint of slate-blue (box fill)
BOX_MEDIAN  = "#2D9B8A"   # teal (median line)

PAIR_COLORS = {
    "tp": "#2D8C4E",
    "fp": "#C0392B",
    "fn": "#E67E22",
}

# First four anchored to NT_COLORS; remaining six fill out the family palette
FAMILY_PALETTE = [
    "#4A6FA5",  # slate-blue  (G)
    "#2D9B8A",  # teal        (U)
    "#E8C32A",  # yellow      (A)
    "#9856A5",  # purple      (C)
    "#7FBBCF",  # light steel-blue
    "#5FA35F",  # sage green
    "#D4A06B",  # warm sand
    "#B5A642",  # olive
    "#9B6B9B",  # lavender
    "#6BA4AA",  # steel teal
]
