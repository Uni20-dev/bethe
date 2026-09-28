"""Check exported DSF normalization and catch plausible-looking corruptions."""
import math
import unittest
from plot_xxx_structure_factor_tutorial import DATA, N, gaussian_spectrum, load_case, read_case


class XXXStructureFactorTutorial(unittest.TestCase):
    def texts(self):
        return {name: (DATA/f'xxx-dsf-n{N}-{name}.csv').read_text() for name in ('spectrum', 'moments')}

    def test_complete_export(self):
        lines, moments = load_case()
        self.assertEqual(len(lines), 528)
        self.assertEqual(len(moments), 64)
        self.assertAlmostEqual(4*sum(x['weight'] for x in lines)/N, .9770578422424, places=10)

    def test_wrong_channel(self):
        texts = self.texts()
        texts['spectrum'] = texts['spectrum'].replace('# Channel: zz', '# Channel: raising')
        with self.assertRaises(ValueError): read_case(texts)

    def test_wrong_weight_and_incomplete_family(self):
        for broken in ('weight', 'missing'):
            texts = self.texts()
            lines = texts['spectrum'].splitlines()
            first = next(i for i, line in enumerate(lines) if line.startswith('state_id,'))+1
            if broken == 'missing': del lines[first]
            else:
                cells = lines[first].split(','); cells[-1] = str(float(cells[-1])*2)
                lines[first] = ','.join(cells)
            texts['spectrum'] = '\n'.join(lines)
            with self.assertRaises(ValueError): read_case(texts)

    def test_gaussian_area(self):
        lines, _ = load_case()
        omega = [-1+i*.001 for i in range(6001)]
        y = gaussian_spectrum(lines, N//4, omega, .08)
        area = .001*(sum(y)-.5*(y[0]+y[-1]))
        expected = 2*math.pi*sum(x['weight'] for x in lines if x['momentum_index'] == N//4)
        self.assertAlmostEqual(area, expected, places=12)
        with self.assertRaises(ValueError): gaussian_spectrum(lines, 16, omega, 0)


if __name__ == '__main__': unittest.main()
