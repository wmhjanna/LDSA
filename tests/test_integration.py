"""Independent analytic checks for missing values, units and input representations."""
import sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from integration import integrate_ldsa,delta_log10
class IntegrationTest(unittest.TestCase):
    def test_area_and_units(self):
        d=np.array([10.,100.,1000.]);x=np.array([[3.,10.,9.]])
        np.testing.assert_allclose(integrate_ldsa(x,d,[.1,.5,.1]),[10*np.pi*.1**2*.5])
    def test_missing_zero_and_bounds(self):
        x=np.array([[np.nan,np.nan],[0.,0.],[np.nan,2.]])
        y=integrate_ldsa(x,[20.,100.],[.5,.5])
        self.assertTrue(np.isnan(y[0]));self.assertEqual(y[1],0.)
        self.assertAlmostEqual(y[2],2*np.pi*.1**2*.5)
    def test_representation(self):
        d=np.array([20.,40.,100.]);x=np.array([[1.,2.,3.]])
        np.testing.assert_allclose(integrate_ldsa(x,d,[.1,.2,.3],representation='dndlog10dp'),integrate_ldsa(x*delta_log10(d),d,[.1,.2,.3]))
    def test_reject_negative(self):
        with self.assertRaises(ValueError):integrate_ldsa([[-1.,2.]],[20.,100.],[.5,.5])
if __name__=='__main__':unittest.main()
