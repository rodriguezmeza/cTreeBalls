"""Radial input units, computation edges, and kappa plotting contracts."""

import sys
from pathlib import Path

# Locate sibling analysis and benchmark helpers for direct script execution.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import kappa_corr_all_engines as driver


class KappaRangeTests(unittest.TestCase):
    def config(self, **kw):
        values=dict(engines=('octree-2balls-omp',),output_dir=Path('unused'),bins=4,multipoles=2,
                    theta_min=2.5,theta_max=250,theta_scale='arcmin')
        values.update(kw)
        return driver.RunConfig(**values).normalized()

    def test_defaults_and_aliases_preserve_native_grid(self):
        args=driver.parse_arguments([])
        self.assertEqual((args.theta_min,args.theta_max,args.theta_scale,args.nbins),
                         (0.12250467913471644,3.6279974066581295,'degree',20))
        for names in (('--theta-min','--theta-max','--theta-scale'),('--min-sep','--max-sep','--sep-units')):
            for equal in (False,True):
                flags=[]
                for name,value in zip(names,('2.5','250','arcmin')):
                    flags += [name+'='+value] if equal else [name,value]
                args=driver.parse_arguments(flags)
                self.assertEqual((args.theta_min,args.theta_max,args.theta_scale),(2.5,250.,'arcmin'))
        chord=2*np.sin(np.deg2rad([2.5,250])/120)
        for units,values in (('degree',np.array([2.5,250])/60),('arcmin',[2.5,250]),
                             ('radian',np.deg2rad([2.5,250])/60),('au',chord),('chord',chord)):
            config=self.config(theta_min=values[0],theta_max=values[1],theta_scale=units)
            np.testing.assert_allclose(driver.angular_limits(config),chord,rtol=2e-15)
            np.testing.assert_allclose(driver.display_radius(chord,config),values,rtol=2e-15)

    def test_all_statistics_and_edge_correction_use_requested_edges(self):
        expected=2*np.sin(np.deg2rad([2.5,250])/120)
        for engine in ('octree-2balls-omp','kdtree-2balls-omp','balltree-2balls-omp'):
            for options,edge in ((('only-2pcf',),False),(('only-3pcf',),False),((),False),(('only-3pcf',),True)):
                config=self.config(options=options,edge_corrections=edge)
                params=driver.engine_parameters(config,engine,False)
                np.testing.assert_allclose([params['rminHist'],params['rangeN']],expected)
                self.assertEqual(params['sizeHistN'],4)
                if edge:self.assertIn('edge-corrections',params['options'])
        if hasattr(driver,'dual_node_common_kwargs'):
            params=driver.dual_node_common_kwargs(self.config())
            np.testing.assert_allclose([params['min_sep'],params['max_sep']],expected)

    def test_invalid_angular_limits(self):
        for values in (dict(theta_max=10801),dict(theta_min=0),dict(theta_max=float('nan')),dict(theta_max=2)):
            with self.assertRaises(ValueError):self.config(**values)
        with self.assertRaisesRegex(ValueError,'chord radii'):
            driver.display_radius(np.array([3.]),self.config())

    def test_saved_grids_and_plot_axes_use_actual_units(self):
        import matplotlib.figure
        chord=2*np.sin(np.deg2rad([2.5,250])/120)
        for log in (False,True):
            for units,values in (('degree',np.array([2.5,250])/60),('arcmin',[2.5,250]),
                                 ('radian',np.deg2rad([2.5,250])/60),('au',chord),('chord',chord)):
                with self.subTest(units=units,log=log), tempfile.TemporaryDirectory() as directory:
                    config=self.config(output_dir=Path(directory),theta_min=values[0],theta_max=values[1],
                                       theta_scale=units,use_log_bins=log)
                    edges=(np.geomspace if log else np.linspace)(*chord,5)
                    r=np.sqrt(edges[:-1]*edges[1:]) if log else (edges[:-1]+edges[1:])/2
                    matrix=np.arange(16.).reshape(4,4);original=matrix.copy()
                    name='octree-2balls-omp'
                    result=dict(r=r,xi=np.arange(4.),zeta_m_0=matrix)
                    driver.attach_radial_grid(result,config,name)
                    np.testing.assert_allclose(result['radial_edges'],edges)
                    if units in ('au','chord'):self.assertNotIn('theta_edges_arcmin',result)
                    else:np.testing.assert_allclose(result['theta_edges_arcmin'][[0,-1]],[2.5,250])
                    captures={}
                    def capture(fig,path,**kwargs):captures[Path(path).name]=fig
                    with patch.object(matplotlib.figure.Figure,'savefig',capture):
                        driver.make_plots({name:result,'kdtree-2balls-omp':dict(result,xi=result['xi']+.1)},config)
                    for ax in captures['two_pcf.png'].axes:
                        np.testing.assert_allclose(ax.get_xlim(),values)
                        self.assertEqual(ax.get_xscale(),'log' if log else 'linear')
                        np.testing.assert_allclose(ax.lines[0].get_xdata(),driver.display_radius(r,config))
                    for fname in ('three_pcf_m0_real.png','three_pcf_m0_relative.png'):
                        for ax in captures[fname].axes:
                            if not ax.get_xlabel():continue  # colorbar
                            np.testing.assert_allclose(ax.get_xlim(),values)
                            np.testing.assert_allclose(ax.get_ylim(),values)
                            mesh=ax.collections[0].get_coordinates()
                            np.testing.assert_allclose(mesh[0,:,0],driver.display_radius(edges,config))
                    table=np.loadtxt(Path(directory)/'plots'/f'{name}_radial_pairs.csv',delimiter=',',skiprows=1)
                    np.testing.assert_array_equal(table[:,0],np.arange(16))
                    np.testing.assert_allclose(table[:,3],driver.display_radius(r,config)[table[:,1].astype(int)])
                    np.testing.assert_allclose(table[:,4],driver.display_radius(r,config)[table[:,2].astype(int)])
                    np.testing.assert_array_equal(matrix,original)

    def test_des_presets_only_fill_omitted_values(self):
        if not hasattr(driver,'resolve_binning'):
            self.skipTest('this driver has no DES preset option')
        for name in ('sofia-fig1','paper-8-200-edges'):
            default=driver.parse_arguments(['--binning',name])
            c=self.config(theta_min=default.theta_min,theta_max=default.theta_max,theta_scale=default.theta_scale,bins=default.nbins)
            edges=driver.radial_edges(c)
            selected=np.sqrt(edges[:-1]*edges[1:])[[0,-1]] if name=='sofia-fig1' else edges[[0,-1]]
            np.testing.assert_allclose(60*np.degrees(2*np.arcsin(selected/2)),[8,200])
            for flags in (['--theta-min=2.5','--theta-max=250','--theta-scale=arcmin','--nbins=6'],
                          ['--min-sep','2.5','--max-sep','250','--sep-units','arcmin','--nbins','6']):
                for order in (['--binning',name]+flags,flags+['--binning',name]):
                    a=driver.parse_arguments(order)
                    self.assertEqual((a.theta_min,a.theta_max,a.theta_scale,a.nbins),(2.5,250.,'arcmin',6))
            for units,scale in (('degree',1/60),('radian',np.pi/10800)):
                a=driver.parse_arguments(['--binning',name,'--theta-scale',units,'--theta-min=.001'])
                self.assertEqual(a.theta_min,.001);self.assertAlmostEqual(a.theta_max,default.theta_max*scale)
            a=driver.parse_arguments(['--binning',name,'--theta-min=0.12250467913471644','--theta-scale=degree'])
            self.assertEqual(a.theta_min,0.12250467913471644)

if __name__=='__main__':unittest.main()
