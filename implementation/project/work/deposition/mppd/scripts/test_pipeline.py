"""Synthetic tests for file validation/mapping, not validation of MPPD physics."""
import argparse, json, tempfile, unittest
from pathlib import Path
import df_pipeline as p

class PipelineTests(unittest.TestCase):
    def curve(self,path,rows):
        p.write_csv(path,[dict(zip(p.FIELDS,r)) for r in rows],p.FIELDS)

    def test_exact_and_tolerance(self):
        sizes=[10.,100.];values={10.:[.1,.2,.3,.6],100.:[.2,.1,.1,.4]}
        self.assertEqual(p.lookup(10.,sizes,values)[1],'MPPD_direct')
        self.assertEqual(p.lookup(10.01,sizes,values)[1],'missing')
        self.assertEqual(p.lookup(10.+1e-8,sizes,values)[1],'MPPD_direct')

    def test_log_interpolation_no_extrapolation(self):
        sizes=[10.,100.];values={10.:[.1,.2,.3,.6],100.:[.2,.1,.1,.4]}
        df,method,lo,hi=p.lookup(1000**.5,sizes,values,'interpolate')
        self.assertAlmostEqual(df[-1],.5);self.assertAlmostEqual(sum(df[:3]),df[-1])
        self.assertEqual(method,'log_diameter_linear_interpolation')
        self.assertEqual(p.lookup(5,sizes,values,'interpolate')[1],'missing')
        self.assertEqual(p.lookup(101,sizes,values,'interpolate')[1],'missing')

    def test_bad_rows_rejected(self):
        cases=[[[10,.1,.2,.3,.9]],[[10,10,20,30,60]],[[0,.1,.2,.3,.6]],[[10,float('nan'),.2,.3,.6]],[[10,.1,.2,.3,.6],[10,.2,.2,.3,.7]]]
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'curve.csv'
            for rows in cases:
                self.curve(f,rows)
                with self.assertRaises(ValueError):p.load_curve(f)

    def test_pending_profile_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'profile.json';f.write_text('{"status":"unconfirmed"}')
            with self.assertRaises(ValueError):p.validate_profile(f,Path(d))

    def test_packaged_inputs(self):
        s=p.validate_inputs(p.ROOT)
        self.assertEqual(s['europe_sites'],31)
        self.assertEqual(s['europe_unique_diameters'],3561)
        self.assertEqual(s['europe_unique_diameters_17_500'],2816)

    def test_end_to_end_synthetic(self):
        # Temporary data never enters the delivered results folder.
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'inputs').mkdir();(root/'raw_exports').mkdir()
            channels=[{'scope':'europe31','site_code':'TEST','grid_id':'g1','channel_index':i+1,'diameter_id':f'D{i+1}','diameter_nm':v,'diameter_um':v/1000,'in_17_500_nm':int(17<=v<=500)} for i,v in enumerate([10.,50.,100.])]
            targets=[{k:r[k] for k in ['diameter_id','diameter_nm','diameter_um']} for r in channels]
            p.write_csv(root/'inputs/station_grid_channels.csv',channels,list(channels[0]));p.write_csv(root/'inputs/target_diameters.csv',targets,list(targets[0]))
            (root/'inputs/input_summary.json').write_text(json.dumps({'channel_rows':3,'unique_diameters_with_izana':3}))
            raw=root/'raw_exports/SYNTHETIC_TEST.txt';raw.write_text('synthetic unit test only')
            cfg=json.loads((p.ROOT/'model_profile.json').read_text(encoding='utf8'))
            for k in ['software_name','software_version','parameter_evidence','diameter_interpretation','region_mapping']:cfg[k]='SYNTHETIC TEST'
            cfg.update(status='confirmed',raw_export_files=['raw_exports/SYNTHETIC_TEST.txt']);cfg['parameters']={k:'TEST' for k in cfg['parameters']}
            profile=root/'profile.json';profile.write_text(json.dumps(cfg))
            f=root/'raw_exports/curve.csv';self.curve(f,[[10,.1,.2,.3,.6],[100,.2,.1,.1,.4]])
            args=argparse.Namespace(root=root,profile=profile,curve=f,scope='europe31',window='native',mode='exact',out=root/'out')
            self.assertEqual(p.match(args),2)
            self.assertEqual(len(p.read_csv(root/'out/unmatched_diameters.csv')),1)
            args.mode='interpolate';self.assertEqual(p.match(args),3)
            self.curve(f,[[10,.1,.2,.3,.6],[50,.15,.15,.2,.5],[100,.2,.1,.1,.4]])
            args.mode='exact';self.assertEqual(p.match(args),0)
            self.assertEqual(len(p.read_csv(root/'out/DF_by_station_grid.csv')),3)

if __name__=='__main__':unittest.main(verbosity=2)
