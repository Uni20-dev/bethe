"""Optional plotting checks: presentation must not change scientific colours/data."""

from pathlib import Path
from importlib.util import find_spec
import tempfile
import unittest
import xml.etree.ElementTree as ET

if find_spec('matplotlib') is not None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, to_rgba
    import numpy as np
else:
    plt = None

from tutorial_common import save_svg
from tutorial_plot_style import DARK_BACKGROUND, dark_figure


@unittest.skipIf(plt is None, 'Install docs/tutorials/requirements.txt for rendering checks')
class FigureStyleTests(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_curves_legends_and_original_are_preserved(self):
        fig, ax = plt.subplots(layout='constrained')
        line, = ax.plot([0, 1], [2, 3], 'o-', color='#0072B2', label='Spinon')
        ax.fill_between([0, 1], [1, 1], [2, 3], color='#0072B2', alpha=.2)
        ax.scatter([.5], [2], facecolors='none', edgecolors='#D55E00')
        ax.legend()
        dark = dark_figure(fig)
        plotted = dark.axes[0].lines[0]
        np.testing.assert_array_equal(plotted.get_xydata(), line.get_xydata())
        self.assertEqual(line.get_color(), '#0072B2')
        self.assertEqual(fig.get_facecolor(), to_rgba('white'))
        self.assertNotEqual(to_rgba(plotted.get_color()), to_rgba(line.get_color()))
        self.assertEqual(to_rgba(dark.axes[0].get_legend().get_lines()[0].get_color()),
                         to_rgba(plotted.get_color()))
        self.assertEqual(dark.get_facecolor(), to_rgba(DARK_BACKGROUND))
        self.assertEqual(dark.axes[0].collections[0].get_alpha(), .2)
        self.assertEqual(len(dark.axes[0].collections[1].get_facecolors()), 0)

    def test_heatmap_and_scalar_scatter_keep_identical_colour_mapping(self):
        fig, ax = plt.subplots(figsize=(10, 5.6), layout='constrained')
        values = np.array([[.01, .1], [1, 10]])
        image = ax.imshow(values, cmap='magma', norm=LogNorm(.01, 10))
        cloud = ax.scatter([0, 1], [0, 1], c=[.1, 1], cmap='gray', norm=LogNorm(.01, 10))
        ax.plot([0, 1], [1, 0], color='white')
        fig.colorbar(image, ax=ax, extend='min')
        fig.colorbar(cloud, ax=ax, extend='both')
        dark = dark_figure(fig)
        for original, copied in [(image, dark.axes[0].images[0]),
                                  (cloud, dark.axes[0].collections[0])]:
            np.testing.assert_array_equal(original.get_array(), copied.get_array())
            np.testing.assert_array_equal(original.to_rgba(values), copied.to_rgba(values))
            self.assertEqual(original.get_clim(), copied.get_clim())
        self.assertEqual(to_rgba(dark.axes[0].lines[0].get_color()), to_rgba('white'))
        for original, copied in zip(fig.axes[1:], dark.axes[1:]):
            self.assertEqual([p.get_facecolor() for p in original.patches],
                             [p.get_facecolor() for p in copied.patches])
        # Compare the actual embedded raster after rendering, not only artist settings.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'heatmap.svg'
            save_svg(fig, path)
            ns = {'svg': 'http://www.w3.org/2000/svg'}
            def rasters(file):
                return [x.attrib['{http://www.w3.org/1999/xlink}href']
                        for x in ET.parse(file).findall('.//svg:image', ns)]
            self.assertTrue(rasters(path))
            self.assertEqual(rasters(path), rasters(path.with_name('heatmap-dark.svg')))

    def test_paired_exports_and_custom_figure_opt_out(self):
        with plt.rc_context({'svg.hashsalt': 'bethe-style-test'}), tempfile.TemporaryDirectory() as directory:
            fig, ax = plt.subplots()
            ax.plot([0, 1], [1, 2])
            path = Path(directory) / 'figure.svg'
            save_svg(fig, path)
            first = path.read_text()
            self.assertTrue(path.with_name('figure-dark.svg').is_file())
            save_svg(fig, path)
            self.assertEqual(first, path.read_text())
            self.assertEqual(plt.get_fignums(), [fig.number])
            custom = Path(directory) / 'custom.svg'
            save_svg(fig, custom, dark=False)
            self.assertTrue(custom.is_file())
            self.assertFalse(custom.with_name('custom-dark.svg').exists())


if __name__ == '__main__':
    unittest.main()
