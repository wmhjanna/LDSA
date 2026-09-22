from pathlib import Path
R=Path('@PROJECT_ROOT@');W=Path(__file__).parent
s=(R/'work/network_trends/prepare.py').read_text().replace("W=R/'work/network_trends';A=R/'outputs/network_trends/analysis'","W=R/'work/quality_control/network';A=R/'outputs/quality_control/netcdf/analysis'")
s=s.replace("outputs/network_trends/corrected_LDSA/hourly_ldsa", "outputs/quality_control/netcdf/hourly_ldsa")
s=s[:s.index('old_input=pd.read_csv')]+'''pd.concat(parts,ignore_index=True).to_csv(W/'trend_input_base.csv',index=False)
pd.DataFrame(elig).to_csv(A/'site_eligibility.csv',index=False)
pd.DataFrame(bps).to_csv(A/'breakpoints.csv',index=False)
pd.concat(week,ignore_index=True).to_pickle(W/'weekly.pkl')
(W/'source_manifest.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2))
pd.DataFrame(mapping)[['old','new','type','station_name']].to_csv(A/'new_station_names.csv',index=False)
'''
(W/'network/prepare.py').write_text(s)
s=(R/'work/network_trends/monthly.R').read_text().replace("work/network_trends","work/quality_control/network").replace("outputs/network_trends/analysis","outputs/quality_control/netcdf/analysis")
s=s[:s.index("old<-read.csv")]+"write.csv(extra,file.path(a,'meta_monthly_series.csv'),row.names=FALSE)\n"
(W/'network/monthly.R').write_text(s)
