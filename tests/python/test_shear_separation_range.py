"""Explicit radial controls reach both statistics and their plotted coordinates."""

import sys
from pathlib import Path

# Locate sibling analysis and benchmark helpers for direct script execution.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

import contextlib
import io
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import shear_corr_all_engines as driver


class SeparationRangeTests(unittest.TestCase):
    def test_preset_defaults_and_explicit_controls(self):
        custom = driver.parse_arguments([])
        self.assertEqual((custom.min_sep, custom.max_sep, custom.nbins, custom.sep_units),
                         (.05, 2., 6, 'degree'))
        for preset in ('sofia-fig1', 'paper-8-200-edges'):
            defaults = driver.parse_arguments(['--binning', preset])
            config = driver.RunConfig(min_sep=defaults.min_sep, max_sep=defaults.max_sep,
                                      bins=defaults.nbins, sep_units=defaults.sep_units)
            edges = np.geomspace(*config.native_limits('sphere'), 21)
            measured = np.sqrt(edges[:-1]*edges[1:])[[0,-1]] if preset=='sofia-fig1' else edges[[0,-1]]
            np.testing.assert_allclose(60*np.degrees(2*np.arcsin(measured/2)), [8,200])
            for limits in (['--min-sep','2.5','--max-sep','250','--sep-units','arcmin','--nbins','9'],
                           ['--min-sep=2.5','--max-sep=250','--sep-units=arcmin','--nbins=9']):
                for options in (['--binning',preset]+limits, limits+['--binning',preset]):
                    actual = driver.parse_arguments(options)
                    self.assertEqual((actual.min_sep,actual.max_sep,actual.nbins,actual.sep_units),
                                     (2.5,250.,9,'arcmin'))
            # Explicit values identical to the old custom defaults still override.
            actual = driver.parse_arguments(['--binning',preset,'--min-sep=.05','--max-sep=2',
                                             '--nbins=6','--sep-units=degree'])
            self.assertEqual((actual.min_sep,actual.max_sep,actual.nbins,actual.sep_units),(.05,2.,6,'degree'))
            for units,scale in (('degree',1/60),('radian',np.pi/10800)):
                actual = driver.parse_arguments(['--binning',preset,'--sep-units',units])
                np.testing.assert_allclose([actual.min_sep,actual.max_sep],
                                           np.array([defaults.min_sep,defaults.max_sep])*scale)
                partial = driver.parse_arguments(['--binning',preset,'--sep-units',units,'--min-sep=0.001'])
                self.assertEqual(partial.min_sep,.001)
                self.assertEqual(partial.max_sep,actual.max_sep)
        for options in (['--linear-bins'],['--sep-units=projected']):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                driver.parse_arguments(['--binning','sofia-fig1']+options)
        custom = driver.parse_arguments(['--linear-bins','--min-sep=2.5','--max-sep=250','--sep-units=arcmin'])
        self.assertTrue(custom.linear_bins)

    def test_main_routes_explicit_range_for_each_statistic(self):
        for statistic in ('2pcf','3pcf','both'):
            with patch.object(driver,'inspect_cballs_runtime',return_value={}), \
                 patch.object(driver,'available_engines',return_value=driver.SHEAR_SPHERE_OMP_ENGINES), \
                 patch.object(driver,'run_engine_suite') as run, contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(driver.main(['--synthetic-nbody=32','--geometry=sphere',
                    '--engine=all-omp','--statistics',statistic,'--binning=sofia-fig1',
                    '--min-sep=2.5','--max-sep=250','--sep-units=arcmin',
                    '--no-smooth-pivot','--no-plots']),0)
            config=run.call_args.args[2]
            self.assertEqual((config.min_sep,config.max_sep,config.bins,config.sep_units),(2.5,250.,20,'arcmin'))
            np.testing.assert_allclose(config.native_limits('sphere'),2*np.sin(np.deg2rad([2.5,250])/120))

    def test_plots_use_limits_units_spacing_and_saved_pair_coordinates(self):
        import matplotlib.figure
        for log in (True,False):
            for units,scale in (('arcmin',1.),('degree',1/60),('radian',np.pi/10800)):
                with tempfile.TemporaryDirectory() as directory:
                    config=driver.RunConfig(min_sep=2.5*scale,max_sep=250*scale,sep_units=units,
                        bins=4,multipoles=2,output_dir=Path(directory),use_log_bins=log)
                    edges=(np.geomspace if log else np.linspace)(*config.native_limits('sphere'),5)
                    radius=np.sqrt(edges[:-1]*edges[1:]) if log else (edges[:-1]+edges[1:])/2
                    gamma=np.ones((4,4,4,5),complex)
                    data=dict(radius=radius,geometry='sphere',gamma=gamma,xi_plus=np.ones(4,complex),
                              xi_minus=np.ones(4,complex),pair_valid=np.ones(4,bool),gamma_valid=np.ones((4,4),bool))
                    captures={}
                    def capture(fig,path,**kwargs):captures[Path(path).name]=fig
                    with patch.object(matplotlib.figure.Figure,'savefig',capture):
                        driver.make_plots({driver.SHEAR_SPHERE_TWO_BALLS_ENGINE:data},config)
                    for ax in captures['shear_2pcf_comparison.png'].axes:
                        np.testing.assert_allclose(ax.get_xlim(),[config.min_sep,config.max_sep])
                        self.assertEqual(ax.get_xscale(),'log' if log else 'linear')
                        if ax.get_xlabel():self.assertEqual(ax.get_xlabel(),f'Angular separation [{units}]')
                        np.testing.assert_allclose(ax.lines[0].get_xdata(),60*np.degrees(2*np.arcsin(radius/2))*scale)
                    pairs=np.loadtxt(Path(directory)/'shear_3pcf_radial_pairs.csv',delimiter=',',skiprows=1)
                    np.testing.assert_array_equal(pairs[:,0],np.arange(16))
                    actual=driver._display_radius(radius,config,'sphere')
                    np.testing.assert_allclose(pairs[:,3],actual[pairs[:,1].astype(int)])
                    np.testing.assert_allclose(pairs[:,4],actual[pairs[:,2].astype(int)])
                    self.assertIn('Measured range',captures['shear_3pcf_gamma_n+0.png']._suptitle.get_text())
                    np.testing.assert_array_equal(gamma,np.ones_like(gamma))

if __name__=='__main__':unittest.main()
