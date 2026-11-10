"""
Generate PowerPoint presentation: Spatial Barcode Layout Optimization
for DNA In-Situ Transcriptomics Arrays
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import copy

# ── Colour palette ──────────────────────────────────────────────────────────
DARK_BG   = RGBColor(0x1A, 0x1A, 0x2E)   # deep navy
MID_BG    = RGBColor(0x16, 0x21, 0x3E)   # mid navy
ACCENT1   = RGBColor(0x0F, 0x3C, 0x96)   # bright blue
ACCENT2   = RGBColor(0xE9, 0x4F, 0x37)   # coral/orange
ACCENT3   = RGBColor(0x39, 0xB5, 0x4A)   # green
ACCENT4   = RGBColor(0xF5, 0xA6, 0x23)   # amber
WHITE     = RGBColor(0xFF, 0xFF, 0xFF)
LIGHT_GRAY= RGBColor(0xCC, 0xD6, 0xE8)
TEXT_DARK = RGBColor(0x1A, 0x1A, 0x2E)

METHOD_COLORS = {
    "Greedy+2-OPT": RGBColor(0x6C, 0x75, 0x7D),
    "Dual":         RGBColor(0x0F, 0x3C, 0x96),
    "Quadtree":     RGBColor(0x39, 0xB5, 0x4A),
    "Zipcodes":     RGBColor(0xE9, 0x4F, 0x37),
}

W = Inches(13.33)   # widescreen width
H = Inches(7.5)     # widescreen height

prs = Presentation()
prs.slide_width  = W
prs.slide_height = H

blank_layout = prs.slide_layouts[6]   # completely blank

# ── Helper utilities ────────────────────────────────────────────────────────

def add_rect(slide, left, top, width, height, fill_color, alpha=None):
    shape = slide.shapes.add_shape(
        1,  # MSO_SHAPE_TYPE.RECTANGLE
        left, top, width, height
    )
    shape.line.fill.background()
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    return shape


def add_textbox(slide, left, top, width, height, text,
                font_size=18, bold=False, color=WHITE,
                align=PP_ALIGN.LEFT, italic=False, wrap=True):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = wrap
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return txBox


def add_para(tf, text, font_size=16, bold=False, color=WHITE,
             align=PP_ALIGN.LEFT, italic=False, space_before=0, bullet=False):
    p = tf.add_paragraph()
    p.alignment = align
    if space_before:
        p.space_before = Pt(space_before)
    run = p.add_run()
    if bullet:
        run.text = "  •  " + text
    else:
        run.text = text
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return p


def slide_background(slide, color=DARK_BG):
    add_rect(slide, 0, 0, W, H, color)


def top_bar(slide, title, subtitle=None, bar_color=ACCENT1):
    """Standard top colour bar with title."""
    bar_h = Inches(1.3)
    add_rect(slide, 0, 0, W, bar_h, bar_color)
    add_textbox(slide, Inches(0.35), Inches(0.15), W - Inches(0.7), Inches(0.75),
                title, font_size=32, bold=True, color=WHITE, align=PP_ALIGN.LEFT)
    if subtitle:
        add_textbox(slide, Inches(0.35), Inches(0.85), W - Inches(0.7), Inches(0.4),
                    subtitle, font_size=16, color=LIGHT_GRAY, align=PP_ALIGN.LEFT)


def bottom_bar(slide, text="Spatial Barcode Layout Optimization"):
    bar_h = Inches(0.35)
    add_rect(slide, 0, H - bar_h, W, bar_h, MID_BG)
    add_textbox(slide, Inches(0.3), H - bar_h + Inches(0.05), Inches(6), bar_h,
                text, font_size=11, color=LIGHT_GRAY, italic=True)


def card(slide, left, top, width, height, title, bullets,
         bg=MID_BG, title_color=ACCENT4, bullet_size=14, title_size=17):
    add_rect(slide, left, top, width, height, bg)
    # title strip
    strip_h = Inches(0.4)
    add_textbox(slide, left + Inches(0.1), top + Inches(0.05),
                width - Inches(0.2), strip_h,
                title, font_size=title_size, bold=True, color=title_color)
    # bullets
    txBox = slide.shapes.add_textbox(
        left + Inches(0.12), top + strip_h + Inches(0.1),
        width - Inches(0.25), height - strip_h - Inches(0.15))
    tf = txBox.text_frame
    tf.word_wrap = True
    first = True
    for b in bullets:
        if first:
            p = tf.paragraphs[0]
            first = False
        else:
            p = tf.add_paragraph()
        p.space_before = Pt(4)
        run = p.add_run()
        run.text = "▸  " + b
        run.font.size = Pt(bullet_size)
        run.font.color.rgb = WHITE


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 1 — Title
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)

# decorative gradient band
add_rect(slide, 0, Inches(2.0), W, Inches(3.2), MID_BG)
add_rect(slide, 0, Inches(2.0), Inches(0.12), Inches(3.2), ACCENT2)

add_textbox(slide, Inches(0.5), Inches(2.15), W - Inches(1.0), Inches(1.1),
            "Spatial Barcode Layout Optimization",
            font_size=42, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
add_textbox(slide, Inches(0.5), Inches(3.3), W - Inches(1.0), Inches(0.7),
            "for DNA In-Situ Transcriptomics Arrays",
            font_size=28, bold=False, color=LIGHT_GRAY, align=PP_ALIGN.CENTER)

add_textbox(slide, Inches(0.5), Inches(4.2), W - Inches(1.0), Inches(0.5),
            "Three Co-Design Strategies: Dual Barcodes  ·  Hybrid Quadtree  ·  Zipcodes",
            font_size=18, color=ACCENT4, bold=True, align=PP_ALIGN.CENTER)

add_textbox(slide, Inches(0.5), Inches(5.2), W - Inches(1.0), Inches(0.4),
            "OncoRNALab  |  Layout Optimization Project",
            font_size=14, color=LIGHT_GRAY, align=PP_ALIGN.CENTER, italic=True)

bottom_bar(slide)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 2 — The Problem
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "The Problem: UV Light Scattering in Photolithographic Synthesis",
        bar_color=ACCENT2)
bottom_bar(slide)

# Left diagram description
lx, ly = Inches(0.3), Inches(1.5)
lw, lh = Inches(4.6), Inches(5.6)
add_rect(slide, lx, ly, lw, lh, MID_BG)

txL = slide.shapes.add_textbox(lx + Inches(0.15), ly + Inches(0.1),
                                lw - Inches(0.3), lh - Inches(0.2))
tfL = txL.text_frame; tfL.word_wrap = True

add_para(tfL, "How Arrays Are Made", 18, bold=True, color=ACCENT4)
add_para(tfL, "Photolithographic arrays synthesize DNA spot-by-spot using UV light to remove "
         "protective groups (deprotection) — allowing the next nucleotide to be added.",
         13, space_before=6)
add_para(tfL, "Each synthesis cycle activates a subset of spots simultaneously, "
         "following a binary schedule (1 = deprotect, 0 = leave protected).",
         13, space_before=6)

add_para(tfL, "The Core Problem", 17, bold=True, color=ACCENT2, space_before=12)
add_para(tfL, "UV light scatters onto neighbouring spots, partially removing "
         "their protective groups → incorrect nucleotide incorporation → insertion errors.",
         13, space_before=6)
add_para(tfL, "The risk scales with feature density: smaller, closer spots suffer more "
         "cross-contamination.", 13, space_before=6)

add_para(tfL, "Synthesis Schedule Representation", 17, bold=True, color=ACCENT3, space_before=12)
add_para(tfL, "Each 36 nt barcode is encoded as a 4×36 = 144-bit one-hot vector "
         "(A=1000, C=0100, G=0010, T=0001).", 13, space_before=6)
add_para(tfL, "Layout cost = sum of Hamming distances between synthesis schedules "
         "of all spatially adjacent pairs.", 13, space_before=6)

# Right: key insight box
rx, ry = Inches(5.2), Inches(1.5)
rw, rh = Inches(7.9), Inches(5.6)
add_rect(slide, rx, ry, rw, rh, MID_BG)

txR = slide.shapes.add_textbox(rx + Inches(0.15), ry + Inches(0.1),
                                rw - Inches(0.3), rh - Inches(0.2))
tfR = txR.text_frame; tfR.word_wrap = True

add_para(tfR, "Key Insight", 20, bold=True, color=ACCENT4)
add_para(tfR, "If two neighbouring spots are BOTH deprotected in the same cycle, "
         "scattered UV hitting spot B while spot A is being activated cannot cause damage "
         "— spot B is already being deprotected anyway.", 14, space_before=6)
add_para(tfR, "→  The MORE similar two neighbours' synthesis schedules are, "
         "the LOWER the synthesis error risk.", 14, bold=True, color=ACCENT3, space_before=8)

add_para(tfR, "Optimization Goal", 18, bold=True, color=ACCENT2, space_before=14)
add_para(tfR, "Minimise: Σ Hamming(schedule_i, schedule_j) for all adjacent pairs (i, j)",
         14, bold=True, color=WHITE, space_before=6)
add_para(tfR, "This problem is MaxSNP-hard (no polynomial-time approximation scheme "
         "exists unless P = NP), motivating heuristic approaches.", 13,
         italic=True, color=LIGHT_GRAY, space_before=6)

add_para(tfR, "Our Approach", 18, bold=True, color=ACCENT4, space_before=14)
add_para(tfR, "Rather than optimizing barcode placement after-the-fact, we explore "
         "co-design strategies that build spatial proximity directly into the barcode "
         "sequences.", 14, space_before=6)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 3 — Project Overview & Array Specs
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Project Overview", "Four layout strategies — all using 36 nt barcodes",
        bar_color=ACCENT1)
bottom_bar(slide)

# array specs card
card(slide, Inches(0.3), Inches(1.45), Inches(3.9), Inches(2.3),
     "Array Specifications",
     ["118 × 182 = 21,476 spots (small array)",
      "235 × 363 = 85,305 spots (large array)",
      "14 µm feature size, 16 µm pitch",
      "All barcodes: 36 nt total length",
      "Full-density & checkerboard (half-density) layouts tested"])

# 4 methods
methods = [
    ("1. Greedy + 2-OPT", ACCENT1,
     ["Post-hoc placement optimisation",
      "GPU-accelerated greedy fill + local swap",
      "Baseline / reference method"]),
    ("2. Dual Barcodes", RGBColor(0x0F, 0x5C, 0xA8),
     ["bar_1 (18 nt) shared per 4×4 block",
      "bar_2 (18 nt) unique within block",
      "Combinatorial clustering"]),
    ("3. Hybrid Quadtree", ACCENT3,
     ["24 nt hierarchical spatial encoding",
      "+ 12 nt high-diversity library",
      "Hamming distance ∝ log(geometric dist)"]),
    ("4. Zipcodes", ACCENT2,
     ["18 nt Gray-coded row index",
      "+ 18 nt Gray-coded col index",
      "Smoothest spatial gradient"]),
]

x0 = Inches(0.3)
for i, (title, col, buls) in enumerate(methods):
    cx = x0 + i * Inches(3.26)
    add_rect(slide, cx, Inches(3.95), Inches(3.1), Inches(0.45), col)
    add_textbox(slide, cx + Inches(0.1), Inches(3.97), Inches(3.0), Inches(0.4),
                title, font_size=15, bold=True, color=WHITE)
    box = slide.shapes.add_textbox(cx + Inches(0.1), Inches(4.45),
                                   Inches(2.9), Inches(2.6))
    tf = box.text_frame; tf.word_wrap = True
    for b in buls:
        p = tf.add_paragraph() if b != buls[0] else tf.paragraphs[0]
        p.space_before = Pt(5)
        run = p.add_run()
        run.text = "▸  " + b
        run.font.size = Pt(13)
        run.font.color.rgb = WHITE

# result teaser
add_rect(slide, Inches(4.2), Inches(1.45), Inches(9.0), Inches(2.3), MID_BG)
txT = slide.shapes.add_textbox(Inches(4.35), Inches(1.5), Inches(8.7), Inches(2.2))
tfT = txT.text_frame; tfT.word_wrap = True
add_para(tfT, "Experimental Snapshot", 17, bold=True, color=ACCENT4)
add_para(tfT, "Arrays were physically synthesised and sequenced. All four methods achieved "
         ">99 % decoding precision.", 13, space_before=5)
add_para(tfT, "Zipcodes showed the highest reads/barcode but exhibited a spatial banding "
         "pattern (lower counts in specific rows/columns).", 13, space_before=5)
add_para(tfT, "Greedy + 2-OPT achieved the lowest synthesis error rate; "
         "it also showed the highest reads/barcode among optimisation-based methods.",
         13, space_before=5)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 4 — Method 1: Greedy + 2-OPT
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Method 1: Greedy + 2-OPT Optimisation",
        "Post-hoc layout optimisation (baseline)", bar_color=ACCENT1)
bottom_bar(slide)

# left
lx, ly = Inches(0.3), Inches(1.5)
lw, lh = Inches(6.0), Inches(5.6)
add_rect(slide, lx, ly, lw, lh, MID_BG)
txL = slide.shapes.add_textbox(lx+Inches(0.15), ly+Inches(0.1), lw-Inches(0.3), lh-Inches(0.2))
tfL = txL.text_frame; tfL.word_wrap = True

add_para(tfL, "Step 1 — Greedy Placement", 17, bold=True, color=ACCENT4)
add_para(tfL, "Fill array column-by-column. For each empty position, compute the Hamming "
         "distance between every remaining barcode's synthesis schedule and those of its "
         "already-placed neighbours. Place the barcode with the minimum total distance.",
         13, space_before=5)
add_para(tfL, "GPU acceleration enables evaluating all ~34,500 candidates at every step.",
         13, italic=True, color=LIGHT_GRAY, space_before=4)

add_para(tfL, "Step 2 — 2-OPT Local Search", 17, bold=True, color=ACCENT4, space_before=10)
add_para(tfL, "Iteratively swap pairs of barcodes if the swap reduces the global layout "
         "cost (total Hamming distance over all adjacent pairs). Continues until no "
         "improving swap exists.", 13, space_before=5)

add_para(tfL, "Barcode Encoding", 17, bold=True, color=ACCENT4, space_before=10)
add_para(tfL, "Each 36 nt barcode → 144-bit one-hot vector (4 bits per base). "
         "Allows uniform-length Hamming distance computation across all barcodes.",
         13, space_before=5)
add_para(tfL, "ACTG → 1000 0100 0001 0010", 12, bold=True,
         color=ACCENT3, space_before=4)

add_para(tfL, "Evaluation Metrics", 17, bold=True, color=ACCENT4, space_before=10)
add_para(tfL, "Mean Hamming distance to 8 neighbours: random baseline → after greedy → "
         "after 2-OPT (cost decreases at each stage).", 13, space_before=5)
add_para(tfL, "Levenshtein distance also computed between neighbouring DNA sequences "
         "to assess real-world sequence similarity.", 13, space_before=5)

# right
rx, ry = Inches(6.5), Inches(1.5)
rw, rh = Inches(6.65), Inches(5.6)
add_rect(slide, rx, ry, rw, rh, MID_BG)
txR = slide.shapes.add_textbox(rx+Inches(0.15), ry+Inches(0.1), rw-Inches(0.3), rh-Inches(0.2))
tfR = txR.text_frame; tfR.word_wrap = True

add_para(tfR, "Synthesis Pattern", 18, bold=True, color=ACCENT4)
add_para(tfR, "A single photolithographic cycle activates a scattered set of small "
         "LOCAL CLUSTERS of spots — groups of features with very similar synthesis "
         "schedules placed next to each other by the greedy algorithm.",
         13, space_before=5)

add_para(tfR, "Advantages", 17, bold=True, color=ACCENT3, space_before=12)
for a in ["Lowest synthesis error rate (fewest cross-activation events)",
          "Works with any barcode set — no structural constraints on sequences",
          "Quantifiable improvement: random → greedy → 2-OPT",
          ">99 % decoding precision in sequencing experiments"]:
    add_para(tfR, a, 13, space_before=3, bullet=True)

add_para(tfR, "Disadvantages", 17, bold=True, color=ACCENT2, space_before=12)
for d in ["Computationally expensive (requires GPU for large arrays)",
          "Post-hoc: barcodes designed first, layout optimised after",
          "Optimisation quality depends heavily on the starting barcode pool",
          "No structural guarantee of spatial similarity"]:
    add_para(tfR, d, 13, space_before=3, bullet=True)

add_para(tfR, "Spatial heatmaps of mean neighbour Hamming distance reveal "
         "local clusters of low-cost regions — confirming the algorithm's effect.",
         12, italic=True, color=LIGHT_GRAY, space_before=10)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 5 — Method 2: Dual Barcodes
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Method 2: Dual Barcodes",
        "Combinatorial block-based clustering (18 nt + 18 nt)", bar_color=ACCENT1)
bottom_bar(slide)

# concept left
lx, ly = Inches(0.3), Inches(1.45)
lw, lh = Inches(6.0), Inches(5.7)
add_rect(slide, lx, ly, lw, lh, MID_BG)
txL = slide.shapes.add_textbox(lx+Inches(0.15), ly+Inches(0.1), lw-Inches(0.3), lh-Inches(0.2))
tfL = txL.text_frame; tfL.word_wrap = True

add_para(tfL, "Core Concept", 18, bold=True, color=ACCENT4)
add_para(tfL, "Partition the array into spatial blocks. Encode position using two "
         "independent 18 nt segments concatenated into one 36 nt barcode.", 13, space_before=5)

add_para(tfL, "bar_1  (18 nt) — Block identifier", 16, bold=True, color=RGBColor(0x6B,0xB8,0xFF), space_before=10)
add_para(tfL, "Shared by ALL spots within the same 4×4 block. "
         "Groups nearby spots under a common synthesis schedule prefix.",
         13, space_before=5)

add_para(tfL, "bar_2  (18 nt) — Position identifier", 16, bold=True, color=ACCENT4, space_before=10)
add_para(tfL, "Unique within each 4×4 block. Distinguishes the 16 individual spots "
         "inside a block from one another.", 13, space_before=5)

add_para(tfL, "Full barcode = bar_1 + bar_2  (36 nt total)", 15,
         bold=True, color=ACCENT3, space_before=8)
add_para(tfL, "Global uniqueness enforced: no two spots share the same bar_1+bar_2 "
         "combination anywhere on the array.", 13, space_before=5)

add_para(tfL, "235 × 363 Array Parameters", 16, bold=True, color=ACCENT4, space_before=10)
for s in ["85,305 total spots",
          "5,369 blocks of 4×4 (edge blocks may have <16 spots)",
          "85,305 barcodes available for both bar_1 and bar_2 pools",
          "Random seed 42 for reproducibility"]:
    add_para(tfL, s, 13, space_before=3, bullet=True)

# block diagram (text art)
rx, ry = Inches(6.5), Inches(1.45)
rw, rh = Inches(6.65), Inches(5.7)
add_rect(slide, rx, ry, rw, rh, MID_BG)
txR = slide.shapes.add_textbox(rx+Inches(0.15), ry+Inches(0.1), rw-Inches(0.3), rh-Inches(0.2))
tfR = txR.text_frame; tfR.word_wrap = True

add_para(tfR, "4 × 4 Block Layout Diagram", 18, bold=True, color=ACCENT4)
add_para(tfR, "All 16 spots share bar_1  ·  Each has a unique bar_2", 13,
         color=LIGHT_GRAY, space_before=4)

block_lines = [
    "┌─────┬─────┬─────┬─────┐",
    "│ b2₁ │ b2₂ │ b2₃ │ b2₄ │  ← bar_1 = X (shared)",
    "├─────┼─────┼─────┼─────┤",
    "│ b2₅ │ b2₆ │ b2₇ │ b2₈ │",
    "├─────┼─────┼─────┼─────┤",
    "│ b2₉ │b2₁₀│b2₁₁│b2₁₂│",
    "├─────┼─────┼─────┼─────┤",
    "│b2₁₃│b2₁₄│b2₁₅│b2₁₆│",
    "└─────┴─────┴─────┴─────┘",
]
for line in block_lines:
    add_para(tfR, line, 12, color=ACCENT3, space_before=1)

add_para(tfR, "Synthesis Pattern", 17, bold=True, color=ACCENT4, space_before=10)
add_para(tfR, "Each photolithographic cycle activates compact 4×4 square blocks — "
         "all spots within a block share bar_1, so they co-activate together in "
         "bar_1 synthesis cycles, minimising cross-block UV scatter.", 13, space_before=5)

add_para(tfR, "Advantages", 15, bold=True, color=ACCENT3, space_before=8)
for a in ["Deterministic spatial structure — no random placement needed",
          "Tolerant to synthesis errors in bar_1 (block-level identification still works)",
          "Efficient: O(rows × cols) time and space"]:
    add_para(tfR, a, 12, space_before=3, bullet=True)

add_para(tfR, "Disadvantages", 15, bold=True, color=ACCENT2, space_before=6)
for d in ["Block boundaries create abrupt transitions between clusters",
          "bar_2 sequences within a block are still random w.r.t. each other",
          "Sequencing quality depends on the 18 nt barcode library design"]:
    add_para(tfR, d, 12, space_before=3, bullet=True)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 6 — Method 3: Hybrid Quadtree
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Method 3: Hybrid Quadtree",
        "Hierarchical spatial encoding with Run-Length Limiting (24 nt + 12 nt)",
        bar_color=ACCENT3)
bottom_bar(slide)

lx, ly = Inches(0.3), Inches(1.45)
lw, lh = Inches(6.3), Inches(5.7)
add_rect(slide, lx, ly, lw, lh, MID_BG)
txL = slide.shapes.add_textbox(lx+Inches(0.15), ly+Inches(0.1), lw-Inches(0.3), lh-Inches(0.2))
tfL = txL.text_frame; tfL.word_wrap = True

add_para(tfL, "Quadtree Spatial Encoding (24 nt)", 18, bold=True, color=ACCENT3)
add_para(tfL, "Recursively subdivide the array into quadrants. At each level, the "
         "quadrant containing a spot (0–3) is mapped to a DNA base (A, C, G, T).",
         13, space_before=5)
add_para(tfL, "3 rounds × 8 cycles = 24 nt quadtree prefix", 13, bold=True,
         color=ACCENT4, space_before=5)
add_para(tfL, "Result: Hamming distance between two barcodes scales logarithmically "
         "with their geometric distance on the array.", 13, space_before=5)

add_para(tfL, "Run-Length Limiting (RLL)", 17, bold=True, color=ACCENT3, space_before=10)
add_para(tfL, "Without constraints, quadtree encoding generates homopolymers "
         "(e.g., AAAA) → synthesis & sequencing problems.", 13, space_before=5)
add_para(tfL, "RLL fix — per synthesis cycle:", 14, bold=True, color=ACCENT4, space_before=5)
for s in ["Build a 4×4 contingency table (quadrant × previous base)",
          "Evaluate all 24 permutations of quadrant→base mapping",
          "Select permutation that minimises homopolymer violations",
          "Forced flips as fallback (<1 % of bases)"]:
    add_para(tfL, s, 12, space_before=3, bullet=True)
add_para(tfL, "max_run = 3  (AAA allowed; AAAA forbidden)", 12, bold=True,
         color=ACCENT3, space_before=5)

add_para(tfL, "Library Extension (12 nt)", 17, bold=True, color=ACCENT3, space_before=10)
add_para(tfL, "12 nt random high-diversity suffix appended to each barcode "
         "to boost sequencing quality and ensure global uniqueness.",
         13, space_before=5)
add_para(tfL, "Source: first 12 bases of the in-house 18 nt barcode library.",
         12, italic=True, color=LIGHT_GRAY, space_before=4)

rx, ry = Inches(6.8), Inches(1.45)
rw, rh = Inches(6.3), Inches(5.7)
add_rect(slide, rx, ry, rw, rh, MID_BG)
txR = slide.shapes.add_textbox(rx+Inches(0.15), ry+Inches(0.1), rw-Inches(0.3), rh-Inches(0.2))
tfR = txR.text_frame; tfR.word_wrap = True

add_para(tfR, "Hierarchical Subdivision Diagram", 18, bold=True, color=ACCENT3)
add_para(tfR, "Level 1: array split into 4 quadrants (A / C / G / T)", 13,
         color=LIGHT_GRAY, space_before=5)
add_para(tfR, "Level 2: each quadrant split again into 4 sub-quadrants", 13,
         color=LIGHT_GRAY, space_before=3)
add_para(tfR, "Level 3 … Level 8: continuing recursion per round", 13,
         color=LIGHT_GRAY, space_before=3)
add_para(tfR,
         "┌──────────┬──────────┐\n"
         "│    A     │    C     │  Level 1\n"
         "│  (AA CC) │  (CA CC) │  Level 2\n"
         "├──────────┼──────────┤\n"
         "│    G     │    T     │\n"
         "│  (GA GC) │  (TA TC) │\n"
         "└──────────┴──────────┘", 11, color=ACCENT3, space_before=8)

add_para(tfR, "Synthesis Pattern", 17, bold=True, color=ACCENT3, space_before=10)
add_para(tfR, "Each cycle activates 16 large quadrant regions — large contiguous "
         "blocks with identical current base, minimising boundary transitions.",
         13, space_before=5)

add_para(tfR, "235 × 363 Config", 15, bold=True, color=ACCENT4, space_before=8)
for s in ["85,305 spots · 24 nt quadtree + 12 nt library = 36 nt",
          "Sequential library assignment (each 12 nt tag used once)",
          "Column-major output ordering"]:
    add_para(tfR, s, 12, space_before=3, bullet=True)

add_para(tfR, "Advantages", 15, bold=True, color=ACCENT3, space_before=8)
for a in ["Mathematically guaranteed spatial similarity gradient",
          "No homopolymers (RLL-controlled)",
          "Combines spatial structure (24 nt) with sequencing quality (12 nt)"]:
    add_para(tfR, a, 12, space_before=3, bullet=True)

add_para(tfR, "Disadvantages", 15, bold=True, color=ACCENT2, space_before=6)
for d in ["Forced flips slightly break the pure hierarchical structure",
          "24 nt alone insufficient for full unique ID — requires library suffix",
          "More complex to implement than block methods"]:
    add_para(tfR, d, 12, space_before=3, bullet=True)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 7 — Method 4: Zipcodes
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Method 4: Zipcodes",
        "Gray-code spatial gradients — the smoothest synthesis pattern (18 nt + 18 nt)",
        bar_color=ACCENT2)
bottom_bar(slide)

lx, ly = Inches(0.3), Inches(1.45)
lw, lh = Inches(6.3), Inches(5.7)
add_rect(slide, lx, ly, lw, lh, MID_BG)
txL = slide.shapes.add_textbox(lx+Inches(0.15), ly+Inches(0.1), lw-Inches(0.3), lh-Inches(0.2))
tfL = txL.text_frame; tfL.word_wrap = True

add_para(tfL, "Gray Code Principle", 18, bold=True, color=ACCENT2)
add_para(tfL, "Gray codes are binary sequences where consecutive integers differ "
         "by exactly ONE bit. This property translates directly to DNA: "
         "spatially adjacent spots differ by exactly ONE base.", 13, space_before=5)

add_para(tfL, "X-code (18 nt) — Column identity", 16, bold=True,
         color=RGBColor(0xFF,0xA0,0x00), space_before=10)
add_para(tfL, "Apply Gray code to column index. Each 2-bit pair → 1 DNA base "
         "(00→A, 01→C, 10→G, 11→T). Adjacent columns differ by exactly 1 base.",
         13, space_before=5)

add_para(tfL, "Y-code (18 nt) — Row identity", 16, bold=True,
         color=RGBColor(0xFF,0xA0,0x00), space_before=10)
add_para(tfL, "Apply Gray code to row index independently. "
         "Adjacent rows differ by exactly 1 base.", 13, space_before=5)

add_para(tfL, "Full barcode = X_code + Y_code  (36 nt, no separator)",
         15, bold=True, color=ACCENT3, space_before=8)

add_para(tfL, "Emergent 2D Spatial Properties", 17, bold=True,
         color=ACCENT2, space_before=10)
for s in ["Horizontal neighbours (same row, ±1 col) → differ by ~2 bases",
          "Vertical neighbours (same col, ±1 row) → differ by ~2 bases",
          "Diagonal neighbours → differ by ~4 bases",
          "Far features (≥5 steps apart on an axis) → differ by ≥8 bases"]:
    add_para(tfL, s, 12, space_before=3, bullet=True)

add_para(tfL, "235 × 363 Config", 16, bold=True, color=ACCENT4, space_before=10)
for s in ["x_bit_length = 36, y_bit_length = 36 (4-bits-per-base scheme)",
          "No spacer sequence between X and Y segments",
          "GC content 40–60 %, max homopolymer = 3",
          "Randomised greedy algorithm with backtracking depth = 5"]:
    add_para(tfL, s, 12, space_before=3, bullet=True)

rx, ry = Inches(6.8), Inches(1.45)
rw, rh = Inches(6.3), Inches(5.7)
add_rect(slide, rx, ry, rw, rh, MID_BG)
txR = slide.shapes.add_textbox(rx+Inches(0.15), ry+Inches(0.1), rw-Inches(0.3), rh-Inches(0.2))
tfR = txR.text_frame; tfR.word_wrap = True

add_para(tfR, "Synthesis Pattern", 18, bold=True, color=ACCENT2)
add_para(tfR, "Because the X-code is determined entirely by column position, "
         "each synthesis cycle activates entire VERTICAL COLUMNS of spots. "
         "This creates the smoothest horizontal gradient but vertical stripes "
         "during synthesis.", 13, space_before=5)

add_para(tfR, "Distance Gradient (per axis)", 16, bold=True, color=ACCENT4, space_before=10)
rows_data = [
    ("Step diff", "Target Hamming"),
    ("1 spot", "~2 bases"),
    ("2 spots", "~4 bases"),
    ("3 spots", "~6 bases"),
    ("4 spots", "~8 bases"),
    ("≥5 spots", "≥8 bases"),
]
for row in rows_data:
    add_para(tfR, f"  {row[0]:<12}  {row[1]}", 12,
             bold=(row[0] == "Step diff"),
             color=ACCENT3 if row[0] == "Step diff" else WHITE,
             space_before=2)

add_para(tfR, "Spatial Error Tolerance", 17, bold=True, color=ACCENT2, space_before=10)
add_para(tfR, "Sequencing errors preferentially misassign to NEARBY positions "
         "(graceful degradation), preserving spatial tissue patterns even under errors.",
         13, space_before=5)

add_para(tfR, "Advantages", 15, bold=True, color=ACCENT3, space_before=8)
for a in ["Smoothest spatial gradient of all methods",
          "Highest reads/barcode in sequencing experiments",
          "Theoretical guarantee: nearby = similar sequences"]:
    add_para(tfR, a, 12, space_before=3, bullet=True)

add_para(tfR, "Disadvantages", 15, bold=True, color=ACCENT2, space_before=6)
for d in ["Vertical stripe synthesis pattern → row/column banding artefact",
          "Observed spatial pattern in read counts: some rows/columns systematically lower",
          "Independent 1D design may not capture all 2D structure"]:
    add_para(tfR, d, 12, space_before=3, bullet=True)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 8 — Synthesis Patterns Comparison
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Synthesis Activation Patterns — One Representative Cycle",
        "How each method distributes deprotection events across the array",
        bar_color=MID_BG)
bottom_bar(slide)

patterns = [
    ("Random Layout", ACCENT1,
     "Spots activated in a cycle are uniformly scattered across the array surface. "
     "Every spot has an equal chance of being a neighbour of an inactive spot → "
     "maximum cross-contamination risk."),
    ("Greedy + 2-OPT", RGBColor(0x6C,0x75,0x7D),
     "Small localised clusters of co-activated spots emerge. The greedy algorithm "
     "groups barcodes with similar synthesis schedules, creating compact 'islands' "
     "of co-activation and reducing inter-cluster scatter."),
    ("Dual Barcodes", RGBColor(0x0F,0x5C,0xA8),
     "Compact 4×4 square blocks are activated uniformly. All 16 spots in a block "
     "share bar_1, so they co-activate in every bar_1 synthesis cycle. Clean, "
     "regular structure with sharp block boundaries."),
    ("Hybrid Quadtree", ACCENT3,
     "16 large quadrant regions are activated. Each cycle separates the array into "
     "4 large zones, then 16, then 64 … The pattern is coarser than Dual Barcodes "
     "but covers very large contiguous regions, minimising total boundary length."),
    ("Zipcodes", ACCENT2,
     "Entire vertical columns are activated together. Because each cycle's base is "
     "determined solely by column position (X-code), whole columns share the same "
     "current base → vertical stripe pattern → explains observed row/column banding."),
]

y0 = Inches(1.45)
bh = Inches(1.06)
gap = Inches(0.04)
for i, (label, col, desc) in enumerate(patterns):
    by = y0 + i * (bh + gap)
    # colour label
    add_rect(slide, Inches(0.3), by, Inches(2.5), bh, col)
    add_textbox(slide, Inches(0.35), by + Inches(0.3), Inches(2.4), Inches(0.5),
                label, font_size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    # description
    add_rect(slide, Inches(2.85), by, Inches(10.15), bh, MID_BG)
    add_textbox(slide, Inches(3.0), by + Inches(0.12),
                Inches(9.9), bh - Inches(0.2),
                desc, font_size=13, color=WHITE, wrap=True)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 9 — Experimental Results
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Experimental Results",
        "Physical synthesis and sequencing of all four layout strategies",
        bar_color=ACCENT1)
bottom_bar(slide)

# Header row
headers = ["Method", "Reads / Barcode", "Decoding Precision", "Synthesis Errors", "Notable Finding"]
col_widths = [Inches(2.0), Inches(2.3), Inches(2.5), Inches(2.5), Inches(3.8)]
header_colors = [ACCENT1, ACCENT1, ACCENT1, ACCENT1, ACCENT1]

y_header = Inches(1.5)
row_h = Inches(0.7)
x_start = Inches(0.3)

# draw header
x = x_start
for w, h_text, hc in zip(col_widths, headers, header_colors):
    add_rect(slide, x, y_header, w - Inches(0.04), row_h, hc)
    add_textbox(slide, x + Inches(0.07), y_header + Inches(0.15),
                w - Inches(0.12), row_h - Inches(0.1),
                h_text, font_size=14, bold=True, color=WHITE)
    x += w

# data rows
rows = [
    ("Greedy + 2-OPT",    "Highest",   "> 99 %",  "Lowest",   "Best synthesis fidelity; small local clusters"),
    ("Dual Barcodes",     "High",      "> 99 %",  "Low",      "Regular 4×4 block pattern; robust decoding"),
    ("Hybrid Quadtree",   "Moderate",  "> 99 %",  "Moderate", "Logarithmic spatial gradient; clean synthesis"),
    ("Zipcodes",          "Highest",   "> 99 %",  "Moderate", "Spatial banding: lower reads in specific rows/cols"),
]
row_colors = [
    RGBColor(0x22,0x2A,0x42),
    RGBColor(0x1A,0x22,0x3A),
    RGBColor(0x22,0x2A,0x42),
    RGBColor(0x1A,0x22,0x3A),
]
for ri, (rdata, rbg) in enumerate(zip(rows, row_colors)):
    y_row = y_header + (ri + 1) * row_h + Inches(0.05)
    x = x_start
    for wi, (cell, cw) in enumerate(zip(rdata, col_widths)):
        add_rect(slide, x, y_row, cw - Inches(0.04), row_h - Inches(0.04), rbg)
        # colour-code reads/errors columns
        txt_col = WHITE
        if wi == 1 and "Highest" in cell:
            txt_col = ACCENT3
        elif wi == 3 and "Lowest" in cell:
            txt_col = ACCENT3
        elif wi == 3 and "Moderate" in cell:
            txt_col = ACCENT4
        elif wi == 4 and "banding" in cell:
            txt_col = ACCENT2
        add_textbox(slide, x + Inches(0.07), y_row + Inches(0.1),
                    cw - Inches(0.12), row_h - Inches(0.1),
                    cell, font_size=13, color=txt_col)
        x += cw

# footnote
add_textbox(slide, Inches(0.3), Inches(5.0), Inches(12.7), Inches(2.1),
            "Key Finding — Zipcodes Spatial Banding: The vertical-stripe synthesis pattern "
            "inherent to Zipcodes (entire columns activated per cycle) leads to systematic "
            "differences in DNA synthesis efficiency between columns. During photolithographic "
            "synthesis, columns with high co-activation frequency may receive slightly "
            "different UV doses from edge effects, resulting in lower read counts for "
            "specific rows and columns in sequencing data. This is a direct consequence "
            "of the column-wise Gray-code structure.",
            font_size=13, color=LIGHT_GRAY, wrap=True)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 10 — Method Comparison
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Method Comparison Summary", bar_color=MID_BG)
bottom_bar(slide)

criteria = [
    "Synthesis error reduction",
    "Spatial similarity guarantee",
    "Sequencing reads/barcode",
    "Computational cost",
    "Implementation complexity",
    "Spatial pattern uniformity",
    "Scalability",
]
scores = {
    # (label, star rating out of 5, color)
    "Greedy+2-OPT": ["★★★★★", "★★★☆☆", "★★★★★", "★★☆☆☆", "★★★☆☆", "★★★★☆", "★★★★☆"],
    "Dual":         ["★★★★☆", "★★★★☆", "★★★★☆", "★★★★★", "★★★★☆", "★★★★☆", "★★★★★"],
    "Quadtree":     ["★★★★☆", "★★★★★", "★★★☆☆", "★★★★☆", "★★★☆☆", "★★★★☆", "★★★★★"],
    "Zipcodes":     ["★★★☆☆", "★★★★★", "★★★★★", "★★★★☆", "★★★★☆", "★★★☆☆", "★★★★★"],
}

methods_order = ["Greedy+2-OPT", "Dual", "Quadtree", "Zipcodes"]
method_labels = ["Greedy + 2-OPT", "Dual Barcodes", "Hybrid Quadtree", "Zipcodes"]
mcols = [RGBColor(0x6C,0x75,0x7D), RGBColor(0x0F,0x5C,0xA8), ACCENT3, ACCENT2]

row_h = Inches(0.53)
col_w = Inches(2.2)
label_w = Inches(3.05)
y0 = Inches(1.55)
x0 = Inches(0.25)

# header
add_rect(slide, x0, y0, label_w, row_h, ACCENT1)
add_textbox(slide, x0+Inches(0.07), y0+Inches(0.1), label_w-Inches(0.1), row_h-Inches(0.1),
            "Criterion", 14, bold=True, color=WHITE)
for mi, (ml, mc) in enumerate(zip(method_labels, mcols)):
    cx = x0 + label_w + mi * col_w
    add_rect(slide, cx, y0, col_w - Inches(0.04), row_h, mc)
    add_textbox(slide, cx+Inches(0.07), y0+Inches(0.1),
                col_w-Inches(0.12), row_h-Inches(0.1),
                ml, 13, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

# data rows
for ri, crit in enumerate(criteria):
    ry = y0 + (ri+1) * row_h + Inches(0.02)
    bg = MID_BG if ri % 2 == 0 else RGBColor(0x22,0x2A,0x42)
    add_rect(slide, x0, ry, label_w, row_h - Inches(0.02), bg)
    add_textbox(slide, x0+Inches(0.07), ry+Inches(0.1),
                label_w-Inches(0.1), row_h-Inches(0.1),
                crit, 13, color=WHITE)
    for mi, mk in enumerate(methods_order):
        cx = x0 + label_w + mi * col_w
        add_rect(slide, cx, ry, col_w - Inches(0.04), row_h - Inches(0.02), bg)
        stars = scores[mk][ri]
        add_textbox(slide, cx+Inches(0.07), ry+Inches(0.1),
                    col_w-Inches(0.12), row_h-Inches(0.1),
                    stars, 13, color=ACCENT4, align=PP_ALIGN.CENTER)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 11 — Future Perspectives
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Future Perspectives & Empirical Validation",
        "From computational design to physical array experiments", bar_color=ACCENT1)
bottom_bar(slide)

steps = [
    ("Physical Synthesis", ACCENT1,
     "Fabricate photolithographic arrays implementing all five layout strategies: "
     "random, greedy + 2-OPT, Dual Barcodes, Quadtree, and Zipcodes. "
     "Use the 235 × 363 design (85,305 spots) for full-scale evaluation."),
    ("Sequencing & QC", ACCENT3,
     "Sequence all arrays under controlled conditions. Measure reads per barcode, "
     "decoding precision, recall, and spatial uniformity. "
     "Directly compare simulation predictions vs. empirical results."),
    ("Synthesis Error Analysis", ACCENT2,
     "Quantify cross-feature contamination and insertion error rates per layout. "
     "Validate that methods with lower Hamming distance between neighbours "
     "genuinely produce fewer synthesis errors."),
    ("Algorithm Refinement", ACCENT4,
     "Use empirical data to refine barcode design parameters, block sizes, "
     "RLL max-run settings, and Gray-code distance targets. "
     "Identify optimal barcode-layout combinations for production-scale arrays."),
]

y0 = Inches(1.5)
sh = Inches(1.35)
gap = Inches(0.05)
for i, (title, col, desc) in enumerate(steps):
    sy = y0 + i * (sh + gap)
    # number circle
    add_rect(slide, Inches(0.3), sy, Inches(0.6), sh, col)
    add_textbox(slide, Inches(0.3), sy + Inches(0.35), Inches(0.6), Inches(0.55),
                str(i+1), font_size=28, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    # title
    add_rect(slide, Inches(0.95), sy, Inches(2.6), sh, col)
    add_textbox(slide, Inches(1.0), sy + Inches(0.4), Inches(2.5), Inches(0.6),
                title, font_size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    # description
    add_rect(slide, Inches(3.6), sy, Inches(9.45), sh, MID_BG)
    add_textbox(slide, Inches(3.75), sy + Inches(0.15),
                Inches(9.2), sh - Inches(0.2),
                desc, font_size=13, color=WHITE, wrap=True)

add_textbox(slide, Inches(0.3), Inches(7.0), Inches(12.5), Inches(0.4),
            "Empirical validation will be critical for confirming whether computational "
            "improvements in layout cost translate to measurable gains in array quality.",
            font_size=13, italic=True, color=LIGHT_GRAY, wrap=True)


# ════════════════════════════════════════════════════════════════════════════
# SLIDE 12 — Conclusions
# ════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(blank_layout)
slide_background(slide)
top_bar(slide, "Conclusions", bar_color=ACCENT1)
bottom_bar(slide)

conclusions = [
    (ACCENT2,  "The Problem",
     "UV light scattering during photolithographic DNA synthesis causes insertion errors "
     "at neighbouring spots. Minimising differences in synthesis schedules between adjacent "
     "features is the key to reducing these errors."),
    (ACCENT1,  "Greedy + 2-OPT",
     "Post-hoc layout optimisation achieves the lowest synthesis error rate and highest "
     "sequencing reads/barcode. It is the most powerful method but requires GPU acceleration "
     "and pre-designed barcode libraries."),
    (RGBColor(0x0F,0x5C,0xA8), "Dual Barcodes",
     "Block-based co-design provides a clean, regular synthesis pattern with 4×4 activation "
     "blocks. Robust decoding and deterministic structure make it an attractive practical option."),
    (ACCENT3,  "Hybrid Quadtree",
     "Hierarchical encoding guarantees a logarithmic spatial similarity gradient. "
     "The RLL algorithm eliminates homopolymers. Combines spatial structure (24 nt) "
     "with high-diversity sequencing quality (12 nt)."),
    (ACCENT2,  "Zipcodes",
     "Gray-code design delivers the smoothest spatial gradient and the highest reads/barcode, "
     "but the column-wise activation pattern introduces systematic spatial banding in "
     "read counts — a direct structural consequence of the encoding scheme."),
]

y0 = Inches(1.45)
bh = Inches(1.05)
gap = Inches(0.03)
for i, (col, title, desc) in enumerate(conclusions):
    by = y0 + i * (bh + gap)
    add_rect(slide, Inches(0.3), by, Inches(0.08), bh, col)
    add_rect(slide, Inches(0.43), by, Inches(12.6), bh, MID_BG)
    add_textbox(slide, Inches(0.55), by + Inches(0.06),
                Inches(2.1), Inches(0.4),
                title, font_size=15, bold=True, color=col)
    add_textbox(slide, Inches(2.7), by + Inches(0.1),
                Inches(10.2), bh - Inches(0.15),
                desc, font_size=13, color=WHITE, wrap=True)


# ════════════════════════════════════════════════════════════════════════════
# Speaker notes — one entry per slide (same order as slides were added)
# ════════════════════════════════════════════════════════════════════════════
speaker_notes = [

    # Slide 1 — Title
    """Welcome, everyone. Today I will walk you through our project on spatial barcode layout
optimization for DNA in-situ transcriptomics arrays. This work tackles a fundamental
challenge in photolithographic array manufacturing: synthesis errors caused by UV light
scattering. We developed and evaluated four distinct strategies — a classical optimization
baseline and three novel co-design methods — all targeting the same 36-nucleotide barcode
format. Let's start with the problem itself.""",

    # Slide 2 — The Problem
    """Photolithographic arrays are built cycle by cycle. In each cycle, a mask selectively
exposes certain spots to UV light, removing a protective group so that the next nucleotide
can be added. The problem is that UV light does not stop exactly at the spot boundary — it
scatters onto adjacent spots, partially deprotecting their surfaces too. This results in
nucleotides being incorrectly incorporated at neighbouring features, predominantly as
insertion errors. As feature sizes shrink and arrays grow denser, this effect becomes
increasingly severe.

The key insight is this: if two neighbouring spots are BOTH being deprotected in the same
cycle, scattered UV is harmless — both spots are already receiving the treatment. The
damage only occurs when one spot is active and its neighbour is not. So our goal is to
design barcode sequences and their spatial arrangement so that neighbours share as many
synthesis cycles as possible. Formally, we minimise the sum of Hamming distances between
the synthesis schedule bit-vectors of all adjacent spot pairs. This problem turns out to
be MaxSNP-hard — there is no efficient exact algorithm — which is why we explore heuristic
and co-design approaches.""",

    # Slide 3 — Project Overview
    """Our project evaluated four methods on two array sizes: a small 118-by-182 array with
21,476 spots, and the full-scale 235-by-363 array with 85,305 spots, at 14-micron feature
size and 16-micron pitch. All barcodes are 36 nucleotides long. We tested both full-density
and checkerboard — half-density — layouts.

The four strategies are: first, Greedy plus 2-OPT, which optimises placement after barcodes
are already designed; second, Dual Barcodes, which splits each 36-mer into two 18-nucleotide
halves with a block structure; third, Hybrid Quadtree, which encodes spatial coordinates
hierarchically; and fourth, Zipcodes, which uses Gray codes to create the smoothest possible
spatial gradient.

Importantly, we were able to physically synthesize and sequence arrays for all four methods.
All achieved greater than 99% decoding precision. Zipcodes gave the highest read counts per
barcode, but showed a spatial banding artefact we will explain shortly. Greedy plus 2-OPT
achieved the lowest synthesis error rate.""",

    # Slide 4 — Greedy + 2-OPT
    """The Greedy plus 2-OPT method is our baseline — it takes a pre-existing library of 36-mer
barcodes and optimises how they are placed on the array.

The greedy step fills the array column by column. At each empty position, it examines all
remaining barcodes and selects the one whose synthesis schedule is most similar to those of
its already-placed neighbours. We measure similarity as the Hamming distance between the
144-bit one-hot encoded vectors of the two barcodes. To do this efficiently for tens of
thousands of candidates at every step, we use GPU acceleration.

After the greedy placement, we apply a 2-OPT local search: iteratively try swapping every
pair of barcodes and keep the swap if it reduces the total layout cost. We continue until no
improving swap exists.

In practice, this creates small local clusters of barcodes with similar synthesis schedules.
Spatial heatmaps of the mean Hamming distance to the 8 nearest neighbours show a clear
reduction in cost compared to a random layout. This method achieved the lowest synthesis
error rate in our experiments, as well as the highest reads per barcode among all tested
strategies.""",

    # Slide 5 — Dual Barcodes
    """The Dual Barcode method takes a different approach: instead of optimising placement after
the fact, it builds spatial structure directly into the barcode design.

The idea is to divide the array into regular 4-by-4 blocks. Every spot in a block shares
the same first 18-nucleotide segment, which we call bar_1. This is the block identifier.
The second 18-nucleotide segment, bar_2, is unique within the block — it tells us which of
the 16 positions inside the block we are reading. The final 36-mer is simply bar_1
concatenated with bar_2.

Because all 16 spots in a block share bar_1, they co-activate together in every synthesis
cycle that involves a base from the bar_1 segment. This means scattered UV within a block
is harmless during those cycles. The block boundaries are the only places where neighbours
have different bar_1 sequences, so that is where cross-contamination risk is concentrated.

For the 235-by-363 array, this gives us 5,369 blocks of 4-by-4, requiring a barcode library
of 85,305 sequences for both the bar_1 and bar_2 pools. Global uniqueness is enforced:
no two spots on the entire array share the same bar_1-plus-bar_2 combination.

The method is deterministic, extremely fast — order N in time and space — and always finds
a valid layout as long as you have at least 4 barcodes.""",

    # Slide 6 — Hybrid Quadtree
    """The Hybrid Quadtree method draws inspiration from spatial indexing data structures.
Instead of grouping spots into flat blocks, it encodes position hierarchically through
recursive subdivision of the array.

In the first round of 8 synthesis cycles, we divide the array into 4 quadrants. Each spot
is assigned the base corresponding to its quadrant — A, C, G, or T. In the second round,
each quadrant is subdivided again into 4 sub-quadrants, and so on. After 3 rounds of 8
cycles, we have a 24-nucleotide prefix that encodes position at multiple scales.
The key mathematical property is that the Hamming distance between two barcodes scales
logarithmically with their geometric distance on the array — nearby spots are more similar.

A serious challenge is homopolymers. Without intervention, the quadtree encoding naturally
produces long runs of the same base. These are problematic both for synthesis and
sequencing. We solve this with a Run-Length Limiting algorithm. At each cycle, we build a
contingency table counting how many spots would get the same base as their previous base
under each of the 24 possible permutations of the quadrant-to-base mapping. We then choose
the permutation that minimises these violations. If a spot still violates the maximum run
length of 3 after the best permutation, we apply a deterministic forced flip. In practice,
fewer than 1% of bases need forced flips.

The final 12 nucleotides are drawn from a high-diversity random library to ensure global
uniqueness and boost sequencing quality. Together, the 24-nucleotide quadtree prefix plus
the 12-nucleotide library suffix give us a 36-mer with both spatial structure and
sequencing robustness.""",

    # Slide 7 — Zipcodes
    """The Zipcode method produces the smoothest possible spatial gradient, using a concept
from coding theory called Gray codes.

A Gray code is a sequence of binary numbers where consecutive values differ by exactly
one bit. If we apply a Gray code to the column index of each spot, consecutive columns
differ by exactly one bit in their code. Translating each 2-bit pair into a DNA base
gives us an 18-nucleotide X-code where adjacent columns differ by exactly one base.
We independently do the same for row indices to get an 18-nucleotide Y-code.
The full barcode is X-code concatenated with Y-code — 36 nucleotides with no separator.

This gives us a beautiful two-dimensional property: horizontal neighbours differ by about
2 bases, vertical neighbours by about 2 bases, diagonal neighbours by about 4 bases,
and features that are 5 or more steps apart differ by at least 8 bases. The gradient is
strictly proportional to Manhattan distance along the array axes.

From a synthesis perspective, because each cycle's base is determined entirely by column
position through the X-code, entire vertical columns of spots are activated together.
This is the smoothest pattern imaginable along the horizontal axis — but it also means
the Y-code synthesis cycles activate entire horizontal rows, creating a striped pattern.

This striped activation is exactly what we see in the sequencing data: certain rows and
columns systematically show lower read counts. This is a direct structural consequence of
the Gray-code design — not a random artefact — and understanding it is key to interpreting
experimental results from Zipcode arrays.""",

    # Slide 8 — Synthesis Patterns
    """Let us look directly at what a single representative photolithographic cycle looks like
for each method — which spots are being deprotected, and which are not.

In a random layout, the activated spots are scattered uniformly across the entire array
surface. Every spot has roughly a 50% chance of being active in any given cycle, and its
neighbours are equally random. This maximises the probability that an active spot sits next
to an inactive one — the worst case for cross-contamination.

After Greedy plus 2-OPT optimisation, we see small localised clusters of co-activated
spots. The algorithm has grouped barcodes with similar synthesis schedules, so features
that co-activate tend to be placed next to each other.

With Dual Barcodes, the pattern is strikingly regular: compact 4-by-4 squares light up
uniformly, separated by dark boundaries. Every spot in a block is either fully active or
fully inactive together — exactly as designed.

The Hybrid Quadtree creates a coarser structure: 16 large rectangular quadrant regions.
Each region is uniformly active or inactive, minimising the total length of the active/
inactive boundary across the array.

Zipcodes produce vertical stripes — entire columns are either on or off. This is the
consequence of the column-wise X-code, and it is the source of the banding pattern
observed in sequencing read counts.""",

    # Slide 9 — Experimental Results
    """Let us now look at the actual sequencing results from physically synthesized arrays.

All four methods achieved greater than 99% decoding precision, meaning that when a read
could be assigned to a barcode, it was overwhelmingly assigned to the correct one. This
validates the barcode design quality across all strategies.

For reads per barcode — the signal strength we get from the array — Greedy plus 2-OPT
and Zipcodes both performed best, while Hybrid Quadtree was more moderate. This reflects
the synthesis error rates: fewer errors mean more reads mapping to the intended barcode.

Synthesis error rates were lowest for Greedy plus 2-OPT, consistent with it having the
lowest layout cost, and highest for Zipcodes, consistent with the striped pattern creating
more boundary events between active and inactive spots at the column edges.

The most striking finding for Zipcodes is the spatial banding pattern. When we plot the
read count per spot on the array, we see clear rows and columns with systematically fewer
reads. This is not noise — it is a direct structural consequence of the column-wise Gray
code. During synthesis, columns that always activate or deactivate together may experience
edge effects in UV dosing, leading to systematic variation in synthesis efficiency along
those axes.""",

    # Slide 10 — Comparison
    """This slide summarises all four methods across seven criteria so you can compare them
at a glance.

Greedy plus 2-OPT wins on synthesis error reduction and reads per barcode, but it is
computationally expensive — it requires GPU acceleration — and it is the only post-hoc
method, meaning barcodes must be designed before layout optimisation can begin.

Dual Barcodes scores highly on spatial similarity guarantee, scalability, and implementation
simplicity. The block structure is intuitive and the algorithm is trivially fast.

Hybrid Quadtree provides the strongest mathematical guarantee of a spatial similarity
gradient — Hamming distance scales logarithmically with distance by construction — but
its implementation is the most complex of the three co-design methods.

Zipcodes deliver the smoothest spatial gradient and the highest read counts, but score
lower on spatial pattern uniformity because of the banding artefact. They are also the
easiest to scale to arbitrary array sizes.

In summary, there is no single winner: the best choice depends on whether you prioritise
synthesis fidelity, sequencing yield, implementation simplicity, or spatial uniformity.""",

    # Slide 11 — Future Perspectives
    """All the results we have shown so far are based on computationally designed and
in-silico evaluated layouts, with the exception of the sequencing experiments just
discussed. The critical next step is a comprehensive empirical validation campaign.

Step one is physical synthesis. We will fabricate arrays implementing all five layout
strategies — including the random baseline — on the full 235-by-363 format, using
identical synthesis conditions for all.

Step two is sequencing and quality control. We will sequence all arrays and measure reads
per barcode, precision, recall, and spatial uniformity maps. Crucially, we will compare
these measurements against our simulation predictions to assess how well the computational
models capture real synthesis behaviour.

Step three is synthesis error analysis. We will quantify cross-feature contamination and
insertion error rates per layout strategy, directly testing whether lower layout cost
translates to fewer synthesis errors.

Step four is algorithm refinement. Armed with empirical data, we can tune design parameters:
block sizes for Dual Barcodes, max-run settings for the Quadtree RLL, and distance targets
for the Zipcode Gray codes. This will allow us to identify the optimal barcode-layout
combination for production-scale spatial transcriptomics arrays.""",

    # Slide 12 — Conclusions
    """To conclude:

UV light scattering during photolithographic synthesis causes insertion errors at
neighbouring spots. The fundamental goal is to ensure that adjacent features share as
many synthesis cycles as possible — quantified as minimising the total Hamming distance
between synthesis schedule vectors of all neighbouring pairs.

Greedy plus 2-OPT is our most powerful method: it achieves the lowest synthesis error
rate and the highest sequencing yield, at the cost of GPU computation and a pre-designed
barcode library.

Dual Barcodes offer a clean, deterministic solution. The 4-by-4 block structure creates
a predictable synthesis pattern, is fast to compute, and is easy to validate.

Hybrid Quadtree provides a mathematically grounded spatial gradient through hierarchical
encoding, with the Run-Length Limiting algorithm preventing homopolymers. It combines
spatial structure with sequencing quality through the hybrid 24-plus-12 architecture.

Zipcodes are theoretically the most elegant design — a smooth, proportional spatial
gradient guaranteed by Gray-code mathematics. In practice, however, the column-wise
activation pattern introduces systematic spatial banding in read counts, which must be
accounted for in downstream analysis.

Together, these four strategies represent a comprehensive exploration of the design space
for spatial barcode layout optimization, and set the stage for empirical validation and
production-scale deployment.

Thank you.""",
]

# Attach notes to slides
for slide_obj, notes_text in zip(prs.slides, speaker_notes):
    notes_slide = slide_obj.notes_slide
    tf = notes_slide.notes_text_frame
    tf.text = notes_text.strip()

# ════════════════════════════════════════════════════════════════════════════
# Save
# ════════════════════════════════════════════════════════════════════════════
outfile = "/user/gent/446/vsc44685/ScratchVO_dir/layout_opt/spatial_barcode_layout_optimization.pptx"
prs.save(outfile)
print(f"Saved: {outfile}")
