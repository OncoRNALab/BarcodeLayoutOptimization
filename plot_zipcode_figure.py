"""
Educational figure: Zipcode barcode strategy
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ── helpers ──────────────────────────────────────────────────────────────────

def gray(n):
    return n ^ (n >> 1)

def to_bits(n, nbits):
    return format(n, f'0{nbits}b')

BASE_MAP = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
BASE_COLOR = {'A': '#4CAF50', 'C': '#2196F3', 'G': '#FF9800', 'T': '#E91E63'}

def bits_to_dna(bits):
    return ''.join(BASE_MAP[bits[i:i+2]] for i in range(0, len(bits), 2))

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

# ── array parameters ──────────────────────────────────────────────────────────
ROWS, COLS = 12, 16
NBITS = 6   # 6 bits → 3-base DNA per axis → 6-base full barcode (for readability)

row_bits = [to_bits(gray(r), NBITS) for r in range(ROWS)]
col_bits = [to_bits(gray(c), NBITS) for c in range(COLS)]
row_dna  = [bits_to_dna(b) for b in row_bits]
col_dna  = [bits_to_dna(b) for b in col_bits]
barcodes = [[col_dna[c] + row_dna[r] for c in range(COLS)] for r in range(ROWS)]

cr, cc = ROWS // 2, COLS // 2
center_bc = barcodes[cr][cc]
ham_grid = np.array([[hamming(barcodes[r][c], center_bc)
                       for c in range(COLS)] for r in range(ROWS)], dtype=float)

# ── color palettes ────────────────────────────────────────────────────────────
DARK   = '#1A1A2E'
MID    = '#16213E'
BLUE   = '#0F3C96'
GREEN  = '#39B54A'
ORANGE = '#F5A623'
CORAL  = '#E94F37'
WHITE  = '#FFFFFF'
LGRAY  = '#CCD6E8'

ham_cmap = LinearSegmentedColormap.from_list(
    'ham', ['#0F3C96', '#39B54A', '#F5A623', '#E94F37'])

# ── figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 13), facecolor=DARK)

gs = gridspec.GridSpec(
    3, 3,
    figure=fig,
    left=0.04, right=0.97,
    top=0.91, bottom=0.04,
    wspace=0.38, hspace=0.48,
    height_ratios=[1, 1.5, 1.4],
    width_ratios=[1.1, 1.3, 1.6],
)

ax_title  = fig.add_subplot(gs[0, :])   # top full-width — title / overview
ax_gray   = fig.add_subplot(gs[1, 0])   # Gray code table
ax_struct = fig.add_subplot(gs[1, 1])   # barcode structure
ax_cycle  = fig.add_subplot(gs[1, 2])   # synthesis cycle
ax_grad   = fig.add_subplot(gs[2, :2])  # Hamming gradient grid
ax_dist   = fig.add_subplot(gs[2, 2])   # distance profile

for ax in fig.axes:
    ax.set_facecolor(MID)
    for spine in ax.spines.values():
        spine.set_edgecolor(LGRAY)
        spine.set_linewidth(0.6)

# ─────────────────────────────────────────────────────────────────────────────
# PANEL 0 — title banner
# ─────────────────────────────────────────────────────────────────────────────
ax_title.set_facecolor(BLUE)
ax_title.axis('off')
ax_title.text(0.5, 0.72, 'Zipcode Barcode Strategy',
              ha='center', va='center', fontsize=26, fontweight='bold',
              color=WHITE, transform=ax_title.transAxes)
ax_title.text(0.5, 0.22,
              'Independent Gray-code encoding of row and column coordinates '
              '→ spatial proximity ≈ sequence similarity',
              ha='center', va='center', fontsize=13, color=LGRAY,
              transform=ax_title.transAxes)

# ─────────────────────────────────────────────────────────────────────────────
# PANEL A — Gray code table (positions 0–7)
# ─────────────────────────────────────────────────────────────────────────────
ax_gray.axis('off')
ax_gray.set_title('A — Gray Code Principle', color=ORANGE, fontsize=12,
                  fontweight='bold', pad=8)

N_SHOW = 8
col_headers = ['Position', 'Binary', 'Gray code', 'DNA']
table_data = []
for i in range(N_SHOW):
    g = gray(i)
    bits = to_bits(g, 4)
    dna = bits_to_dna(bits)
    table_data.append([str(i), to_bits(i, 4), bits, dna])

cell_w = 0.22
cell_h = 0.085
x0, y0 = 0.04, 0.88

# header
for ci, hdr in enumerate(col_headers):
    ax_gray.add_patch(FancyBboxPatch(
        (x0 + ci * cell_w, y0 - cell_h), cell_w - 0.01, cell_h - 0.005,
        boxstyle='round,pad=0.01', facecolor=BLUE, edgecolor='none',
        transform=ax_gray.transAxes))
    ax_gray.text(x0 + ci * cell_w + cell_w / 2, y0 - cell_h / 2,
                 hdr, ha='center', va='center', fontsize=9,
                 fontweight='bold', color=WHITE, transform=ax_gray.transAxes)

# rows
for ri, row in enumerate(table_data):
    y = y0 - (ri + 2) * cell_h
    changed = (ri > 0 and row[2] != table_data[ri - 1][2])
    # highlight changed bit in gray code
    for ci, val in enumerate(row):
        bg = MID if ri % 2 == 0 else '#1e2a40'
        ax_gray.add_patch(FancyBboxPatch(
            (x0 + ci * cell_w, y), cell_w - 0.01, cell_h - 0.005,
            boxstyle='round,pad=0.01', facecolor=bg, edgecolor='none',
            transform=ax_gray.transAxes))
        fc = WHITE
        if ci == 3:   # DNA column — colour per base
            for bi, base in enumerate(val):
                ax_gray.text(x0 + ci * cell_w + 0.015 + bi * 0.045,
                             y + cell_h / 2, base,
                             ha='center', va='center', fontsize=9,
                             fontweight='bold', color=BASE_COLOR[base],
                             transform=ax_gray.transAxes)
            continue
        if ci == 2 and changed:   # highlight the single bit that flipped
            prev_bits = table_data[ri - 1][2]
            cur_bits = val
            for bi, (pb, cb) in enumerate(zip(prev_bits, cur_bits)):
                col = CORAL if pb != cb else WHITE
                ax_gray.text(x0 + ci * cell_w + 0.018 + bi * 0.046,
                             y + cell_h / 2, cb,
                             ha='center', va='center', fontsize=9,
                             fontweight='bold', color=col,
                             transform=ax_gray.transAxes)
            continue
        ax_gray.text(x0 + ci * cell_w + cell_w / 2, y + cell_h / 2,
                     val, ha='center', va='center', fontsize=9,
                     color=fc, transform=ax_gray.transAxes)

ax_gray.text(0.5, 0.03,
             'Red bit = the single bit that changed\n'
             '→ consecutive positions differ by exactly 1 base',
             ha='center', va='bottom', fontsize=8, color=LGRAY,
             transform=ax_gray.transAxes, style='italic')

# ─────────────────────────────────────────────────────────────────────────────
# PANEL B — barcode structure diagram
# ─────────────────────────────────────────────────────────────────────────────
ax_struct.axis('off')
ax_struct.set_title('B — Barcode Construction', color=ORANGE, fontsize=12,
                    fontweight='bold', pad=8)

def draw_segment(ax, x, y, w, h, label, sublabel, color, transform):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle='round,pad=0.02', facecolor=color,
                 edgecolor=WHITE, linewidth=1.2, transform=transform))
    ax.text(x + w / 2, y + h * 0.62, label, ha='center', va='center',
            fontsize=11, fontweight='bold', color=WHITE, transform=transform)
    ax.text(x + w / 2, y + h * 0.22, sublabel, ha='center', va='center',
            fontsize=8.5, color='#ffffffcc', transform=transform)

t = ax_struct.transAxes

# example spot
ex_col, ex_row = 3, 5
x_seq = col_dna[ex_col]
y_seq = row_dna[ex_row]

draw_segment(ax_struct, 0.05, 0.72, 0.40, 0.16,
             f'X-code  (col={ex_col})', f'Gray({ex_col}) → {col_bits[ex_col]} → {x_seq}',
             '#0F5CA8', t)
ax_struct.text(0.50, 0.80, '+', ha='center', va='center',
               fontsize=18, color=WHITE, transform=t)
draw_segment(ax_struct, 0.55, 0.72, 0.40, 0.16,
             f'Y-code  (row={ex_row})', f'Gray({ex_row}) → {row_bits[ex_row]} → {y_seq}',
             CORAL, t)

# arrow down
ax_struct.annotate('', xy=(0.5, 0.60), xytext=(0.5, 0.69),
                   xycoords=t, textcoords=t,
                   arrowprops=dict(arrowstyle='->', color=ORANGE, lw=2))

full_bc = x_seq + y_seq
draw_segment(ax_struct, 0.10, 0.48, 0.80, 0.16,
             'Full barcode (36 nt)', f'{full_bc}  · · ·  (scaled for illustration)',
             GREEN, t)

# distance property boxes
props = [
    ('Horizontal\nneighbours', '~2 bases\ndifferent', BLUE),
    ('Vertical\nneighbours', '~2 bases\ndifferent', BLUE),
    ('Diagonal\nneighbours', '~4 bases\ndifferent', ORANGE),
    ('Far spots\n(≥5 steps)', '≥8 bases\ndifferent', CORAL),
]
bw, bh = 0.20, 0.18
bx0, by0 = 0.03, 0.22
for i, (lbl, val, col) in enumerate(props):
    bx = bx0 + i * (bw + 0.03)
    ax_struct.add_patch(FancyBboxPatch((bx, by0), bw, bh,
                        boxstyle='round,pad=0.02', facecolor=col,
                        edgecolor='none', alpha=0.85, transform=t))
    ax_struct.text(bx + bw / 2, by0 + bh * 0.65, lbl, ha='center', va='center',
                   fontsize=8, fontweight='bold', color=WHITE, transform=t)
    ax_struct.text(bx + bw / 2, by0 + bh * 0.20, val, ha='center', va='center',
                   fontsize=8.5, color='#ffffcc', fontweight='bold', transform=t)

ax_struct.text(0.5, 0.06,
               'Distance increases proportionally with Manhattan distance',
               ha='center', va='bottom', fontsize=8.5, color=LGRAY,
               transform=t, style='italic')

# ─────────────────────────────────────────────────────────────────────────────
# PANEL C — synthesis cycle (vertical stripes)
# ─────────────────────────────────────────────────────────────────────────────
ax_cycle.set_title('C — Single Synthesis Cycle\n(column-wise activation)',
                   color=ORANGE, fontsize=12, fontweight='bold', pad=8)

# cycle: 1st base of x-code determines activation (A → active)
# show which columns are 'active' (base == 'A') vs 'inactive'
active_bases = [col_dna[c][0] for c in range(COLS)]   # first base of x-code
cycle_grid = np.zeros((ROWS, COLS))
for c in range(COLS):
    if active_bases[c] == 'A':
        cycle_grid[:, c] = 1.0

cycle_cmap = LinearSegmentedColormap.from_list('cyc', [MID, BLUE])
ax_cycle.imshow(cycle_grid, cmap=cycle_cmap, aspect='auto',
                vmin=0, vmax=1, interpolation='nearest')

# column labels
for c in range(COLS):
    col_letter = active_bases[c]
    color = WHITE if cycle_grid[0, c] == 1 else LGRAY
    alpha = 1.0 if cycle_grid[0, c] == 1 else 0.5
    ax_cycle.text(c, -0.8, col_letter, ha='center', va='center',
                  fontsize=8, color=color, fontweight='bold', alpha=alpha)

ax_cycle.set_xlim(-0.5, COLS - 0.5)
ax_cycle.set_ylim(ROWS - 0.5, -1.5)
ax_cycle.set_xticks(range(COLS))
ax_cycle.set_xticklabels([str(i) for i in range(COLS)], fontsize=7, color=LGRAY)
ax_cycle.set_yticks([])
ax_cycle.set_xlabel('Column index', fontsize=9, color=LGRAY)

# labels on active columns
for c in range(COLS):
    if cycle_grid[0, c] == 1:
        ax_cycle.axvline(c, color=ORANGE, linewidth=1.2, alpha=0.4)

legend_els = [
    mpatches.Patch(color=BLUE, label='Active (deprotected)'),
    mpatches.Patch(color=MID,  label='Inactive (protected)'),
]
ax_cycle.legend(handles=legend_els, loc='lower right', fontsize=8,
                facecolor=DARK, edgecolor=LGRAY, labelcolor=WHITE)

ax_cycle.text(0.5, 1.08, 'Entire columns activate together\n→ vertical stripe pattern',
              ha='center', va='bottom', fontsize=8.5, color=LGRAY,
              transform=ax_cycle.transAxes, style='italic')

# ─────────────────────────────────────────────────────────────────────────────
# PANEL D — Hamming distance gradient grid
# ─────────────────────────────────────────────────────────────────────────────
ax_grad.set_title(
    f'D — Hamming Distance from Centre Spot ({cr},{cc})  '
    '— colour shows sequence similarity gradient',
    color=ORANGE, fontsize=12, fontweight='bold', pad=8)

im = ax_grad.imshow(ham_grid, cmap=ham_cmap, aspect='equal',
                    interpolation='nearest',
                    vmin=0, vmax=ham_grid.max())

# annotate each cell with distance
for r in range(ROWS):
    for c in range(COLS):
        d = int(ham_grid[r, c])
        fc = WHITE if d > ham_grid.max() * 0.4 else DARK
        ax_grad.text(c, r, str(d), ha='center', va='center',
                     fontsize=9, color=fc, fontweight='bold')

# mark centre
ax_grad.add_patch(plt.Rectangle((cc - 0.5, cr - 0.5), 1, 1,
                                  fill=False, edgecolor=WHITE,
                                  linewidth=3, zorder=5))
ax_grad.text(cc, cr - 0.65, '★', ha='center', va='center',
             fontsize=13, color=WHITE, zorder=6)

# axis labels
ax_grad.set_xticks(range(COLS))
ax_grad.set_xticklabels([str(c) for c in range(COLS)], fontsize=8, color=LGRAY)
ax_grad.set_yticks(range(ROWS))
ax_grad.set_yticklabels([str(r) for r in range(ROWS)], fontsize=8, color=LGRAY)
ax_grad.set_xlabel('Column', fontsize=10, color=LGRAY, labelpad=4)
ax_grad.set_ylabel('Row', fontsize=10, color=LGRAY, labelpad=4)

cbar = fig.colorbar(im, ax=ax_grad, fraction=0.025, pad=0.02)
cbar.set_label('Hamming distance\n(number of differing bases)',
               fontsize=9, color=LGRAY)
cbar.ax.yaxis.set_tick_params(color=LGRAY, labelsize=8)
plt.setp(cbar.ax.yaxis.get_ticklabels(), color=LGRAY)
cbar.outline.set_edgecolor(LGRAY)

# ─────────────────────────────────────────────────────────────────────────────
# PANEL E — distance vs. Manhattan distance profile
# ─────────────────────────────────────────────────────────────────────────────
ax_dist.set_title('E — Distance Profile from Centre',
                  color=ORANGE, fontsize=12, fontweight='bold', pad=8)
ax_dist.set_facecolor(MID)

manhattan = np.array([[abs(r - cr) + abs(c - cc)
                        for c in range(COLS)] for r in range(ROWS)]).ravel()
hamming_flat = ham_grid.ravel()

# scatter
scatter = ax_dist.scatter(manhattan, hamming_flat,
                          c=hamming_flat, cmap=ham_cmap,
                          s=35, alpha=0.7, edgecolors='none',
                          vmin=0, vmax=ham_grid.max())

# trend line (mean Hamming per Manhattan distance)
max_m = int(manhattan.max())
means_x, means_y = [], []
for m in range(max_m + 1):
    mask = manhattan == m
    if mask.sum() > 0:
        means_x.append(m)
        means_y.append(hamming_flat[mask].mean())

ax_dist.plot(means_x, means_y, color=ORANGE, linewidth=2.5,
             marker='o', markersize=5, label='Mean Hamming')

ax_dist.set_xlabel('Manhattan distance (steps)', fontsize=10, color=LGRAY)
ax_dist.set_ylabel('Hamming distance (bases)', fontsize=10, color=LGRAY)
ax_dist.tick_params(colors=LGRAY, labelsize=8)
ax_dist.spines['bottom'].set_color(LGRAY)
ax_dist.spines['left'].set_color(LGRAY)
ax_dist.legend(fontsize=9, facecolor=DARK, edgecolor=LGRAY, labelcolor=WHITE)
ax_dist.text(0.05, 0.92, 'Hamming distance ∝ spatial distance',
             transform=ax_dist.transAxes, fontsize=9, color=LGRAY, style='italic')

# grid
ax_dist.grid(True, color='#2a3a5a', linewidth=0.5, linestyle='--')

# ─────────────────────────────────────────────────────────────────────────────
# Overall figure label
# ─────────────────────────────────────────────────────────────────────────────
fig.text(0.5, 0.975,
         'Figure: Zipcode Spatial Barcode Strategy for Photolithographic Array Synthesis',
         ha='center', va='top', fontsize=15, fontweight='bold',
         color=WHITE)

outfile = '/user/gent/446/vsc44685/ScratchVO_dir/layout_opt/zipcode_strategy.png'
plt.savefig(outfile, dpi=180, bbox_inches='tight', facecolor=DARK)
print(f"Saved: {outfile}")
