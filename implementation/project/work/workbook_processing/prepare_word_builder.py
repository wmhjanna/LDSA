from pathlib import Path
import shutil
R=Path('@PROJECT_ROOT@');W=R/'work/table_formatting';W.mkdir(exist_ok=True)
for name in ['SP-nature(2).docx.json']:shutil.copy2(R/'work/word_tables_integration_support'/name,W/name)
t=(R/'work/word_tables_integration_support/build_tables.py').read_text().replace('word_tables_integration_support','table_formatting').replace('integration_support','workbook_processing').replace('LDSA_supplementary_tablesexport_study_method17-500_quality_controlled_data_2026-09-11.docx','LDSA_supplementary_tablesexport_workbook records_PARrestore_PRAcombine_2026-09-13.docx')
t=t.replace("sites=['PRA_2012-2016_SUB','PRA_2017-2019_SUB'] if old=='PRA_SUB' else [oldmaps.get(old,old)]","sites=[oldmaps.get(old,old)]")
t=t.replace("mapped=sorted(mapped", "par_template=list(templates['PRA_SUB'][1]);par_template[0]='PAR_SUB';par_template[1]='FR';templates['PAR_SUB']=(templates['PRA_SUB'][0],par_template);mapped.append('PAR_SUB')\nmapped=sorted(mapped")
t=t.replace('The input records represent observational sites.','The analysis includes observational sites, each represented by one record.').replace('Prague records remain separate.','Paris and Prague configurations are combined after calculation with their respective native diameters and deposition fractions.').replace('Counts refer to independent records, with Prague configurations kept separate.','Counts refer to physical sites, with one combined record each for Paris and Prague.').replace('30 physical sites and 31 independent records','31 physical sites').replace('input records at observational sites','observational sites').replace('Prague is not merged.','Paris and Prague each contribute one combined record.')
(W/'build_tables.py').write_text(t)
print(W/'build_tables.py')
(W/'document_contents.md').write_text('# Supplementary tables\n\nThe six-table document includes one combined record for each physical station. Paris and Prague configurations are integrated on their native grids before combination.\n')
