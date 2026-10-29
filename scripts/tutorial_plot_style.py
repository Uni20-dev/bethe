"""Presentation-only dark companions for Matplotlib tutorial figures.

The caller's light figure is untouched. Work on artists, never on exported
pixels: scalar-mapped images/scatter/meshes keep their data, normalization and
colormaps. Unknown categorical colours are left for the figure author to choose.
"""

from copy import deepcopy


DARK_BACKGROUND = '#1e2129'
DARK_FOREGROUND = '#e5e7eb'
# Brighter counterparts of the tutorial's common categorical colours.
# Hue families and line/marker styles still identify the same physical series.
DARK_COLORS = {
    '#0072b2': '#56b4e9',
    '#d55e00': '#f09956',
    '#009e73': '#62c5a1',
    '#1f77b4': '#69b3e7',
    '#2ca02c': '#74c476',
    '#d62728': '#f08080',
    '#9467bd': '#bc9bdd',
}


def dark_color(color):
    from matplotlib.colors import to_hex, to_rgba

    if isinstance(color, str) and color == 'auto':
        return color
    rgba = to_rgba(color)
    if rgba[3] == 0:
        return rgba
    replacement = DARK_COLORS.get(to_hex(rgba))
    if replacement is not None:
        return to_rgba(replacement, rgba[3])
    r, g, b, alpha = rgba
    if r == g == b and r < .6:
        # Black references and dark-grey grid lines need contrast. White
        # overlays on heatmaps remain white, rather than being inverted.
        return to_rgba(DARK_FOREGROUND if r == 0 else '#9ca3af', alpha)
    return rgba


def dark_figure(figure):
    """Return an independent dark rendering; retain geometry and colour scales."""
    from matplotlib.collections import Collection
    from matplotlib.cm import ScalarMappable
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.text import Text

    result = deepcopy(figure)
    # save_svg has already laid out the light figure. A second constrained
    # layout pass can shift a colourbar by a pixel even though only colours
    # changed. Freeze those positions so the two exports align exactly.
    result.set_layout_engine('none')
    colorbar_axes = {artist.colorbar.ax for artist in result.findobj()
                     if isinstance(artist, ScalarMappable) and artist.colorbar is not None}
    # Include legend samples as well as the plotted artists so keys and curves
    # remain consistent. Scalar-mapped collections are deliberately excluded.
    for artist in result.findobj():
        if isinstance(artist, Text):
            artist.set_color(dark_color(artist.get_color()))
        elif isinstance(artist, Line2D):
            artist.set_color(dark_color(artist.get_color()))
            artist.set_markerfacecolor(dark_color(artist.get_markerfacecolor()))
            artist.set_markeredgecolor(dark_color(artist.get_markeredgecolor()))
        elif isinstance(artist, Collection):
            if artist.get_array() is None:
                artist.set_facecolor([dark_color(c) for c in artist.get_facecolor()])
                artist.set_edgecolor([dark_color(c) for c in artist.get_edgecolor()])
        elif isinstance(artist, Patch):
            # Colourbar extension triangles encode under/over-range values,
            # including black in a greyscale map. They are scientific colours.
            if artist.axes not in colorbar_axes:
                artist.set_facecolor(dark_color(artist.get_facecolor()))
                artist.set_edgecolor(dark_color(artist.get_edgecolor()))
    result.set_facecolor(DARK_BACKGROUND)
    for axes in result.axes:
        axes.set_facecolor(DARK_BACKGROUND)
        # Applies to ticks created lazily during rendering, including log axes.
        axes.tick_params(axis='both', which='both', colors=DARK_FOREGROUND)
        for spine in axes.spines.values():
            spine.set_edgecolor('#9ca3af')
        legend = axes.get_legend()
        if legend is not None:
            legend.get_frame().set_facecolor(DARK_BACKGROUND)
            legend.get_frame().set_edgecolor('#6b7280')
    for legend in result.legends:
        legend.get_frame().set_facecolor(DARK_BACKGROUND)
        legend.get_frame().set_edgecolor('#6b7280')
    return result
