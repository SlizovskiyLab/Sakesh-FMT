import matplotlib.pyplot as plt
import cairosvg
from PIL import Image
from io import BytesIO

# Grid config
n_rows = 6
n_cols = 3
row_labels = ["Prep", "PrePostDonor", "Extraction", "Route", "Sequencer", "Study"]
col_labels = ["MDRB", "melanoma", "rCDI"]

# File paths
image_paths = [
    ["C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Prep/pca_mdrb.svg",        "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Prep/pca_melanoma.svg",        "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Prep/pca_rcdi.svg"],
    ["C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/pca_mdrb.svg","C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/pca_melanoma.svg","C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/pca_rcdi.svg"],
    ["C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/pca_mdrb.svg",  "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/pca_melanoma.svg",  "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/pca_rcdi.svg"],
    ["C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Route/pca_mdrb.svg",       "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Route/pca_melanoma.svg",       "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Route/pca_rcdi.svg"],
    ["C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Sequencer/pca_mdrb.svg",   "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Sequencer/pca_melanoma.svg",   "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Sequencer/pca_rcdi.svg"],
    ["C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/pca_mdrb.svg",       "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/pca_melanoma.svg",       "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/pca_rcdi.svg"]
]

# Render SVG to PNG in-memory
def render_svg_to_ax(svg_path, ax):
    try:
        png_bytes = cairosvg.svg2png(url=svg_path, dpi=300)
        image = Image.open(BytesIO(png_bytes))
        ax.imshow(image)
        ax.axis('off')
    except Exception as e:
        ax.set_facecolor("lightgrey")
        ax.text(0.5, 0.5, "Missing", ha='center', va='center', fontsize=8)
        print(f" Failed to render {svg_path}: {e}")

# Create figure and axes
fig, axs = plt.subplots(n_rows, n_cols, figsize=(9, 15))
plt.subplots_adjust(wspace=0.01, hspace=0.01, left=0.15, top=0.92)

# Fill in plots
for row in range(n_rows):
    for col in range(n_cols):
        ax = axs[row, col]
        render_svg_to_ax(image_paths[row][col], ax)

        # Column titles (top row)
        if row == 0:
            ax.set_title(col_labels[col], fontsize=10, fontweight='bold', pad=5)

        # Row labels (horizontal, aligned left of grid)
        if col == 0:
            ax.text(-0.15, 0.5, row_labels[row],
                    va='center', ha='right', fontweight='bold',
                    fontsize=9, transform=ax.transAxes)

# Main figure title
fig.suptitle("PCA of Aitchison Distances for Mobilome Samples",
             fontsize=12, fontweight='bold', y=0.96)

# Save output
fig.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/pca_combined_final.svg",
            format='svg', dpi=300, bbox_inches='tight')

plt.show()
