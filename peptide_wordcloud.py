import pandas as pd
import matplotlib.pyplot as plt
from wordcloud import WordCloud
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

data_df_high_score = pd.read_csv("/Users/chaves/Documents/peptide.csv") # change this path to yours
peptide_counts = data_df_high_score["stripped_sequence"].value_counts() # find the column with the peptide senquences. stripped_sequence in this example
peptide_freq = dict(peptide_counts)
colors = ["#d4e6f1", "#a9cce3", "#7fb3d5", "#5499c7", "#2980b9", "#1f618d", "#154360"]
cmap = LinearSegmentedColormap.from_list("peptide_colors", colors, N=256)


def create_circular_mask(size=800):
    x, y = np.ogrid[:size, :size]
    center = size // 2
    radius = int(size * 0.48)
    mask = (x - center) ** 2 + (y - center) ** 2 > radius**2
    img = np.ones((size, size), dtype=np.uint8)
    img[mask] = 0
    return img


mask = create_circular_mask()

wordcloud = WordCloud(
    width=1000,
    height=1000,
    background_color="white",
    colormap=cmap,
    mask=mask,
    contour_width=1,
    contour_color="#2c3e50",
    min_font_size=6,
    max_font_size=150,
    prefer_horizontal=1,
    relative_scaling=0.7,  # This controls how much bigger frequent words appear
    font_path=None,  # Specify a font path if you want to use a custom font
    random_state=42,
).generate_from_frequencies(peptide_freq)

plt.figure(figsize=(12, 12))
plt.imshow(np.array(wordcloud.to_image()), interpolation="bilinear")
plt.axis("off")
plt.tight_layout(pad=0)

plt.title(
    "Peptide sequence frequency view",
    fontsize=35,
    pad=20,
    fontweight="bold",
    color="#2c3e50",
)

plt.subplots_adjust(bottom=0.1)
cax = plt.axes([0.25, 0.05, 0.5, 0.02])
norm = plt.Normalize(min(peptide_freq.values()), max(peptide_freq.values()))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, cax=cax, orientation="horizontal")
cbar.set_label("Peptide frequency", fontsize=25, color="#2c3e50")

plt.savefig("peptide_wordcloud.png", dpi=300, bbox_inches="tight")
plt.show()
