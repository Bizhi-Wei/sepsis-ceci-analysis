from pathlib import Path
import os
import hashlib,json,duckdb,pandas as pd,numpy as np
P=Path(os.environ["CECI_WORK_DIR"]).resolve(strict=True)
O=P/'eicu_analysis'
O.mkdir(exist_ok=True)
private=O/'restricted_local_only';private.mkdir(exist_ok=True)
root=Path(os.environ["EICU_DIR"]).resolve(strict=True)
a=pd.read_csv(P/'01_data/processed/mimic_sepsis_endothelial_cohort.csv')
mu={};sd={};raw=[]
for name,col,sign,trans in [('lactate','lactate_24h',1,'log2'),('albumin','albumin_24h',-1,'identity'),('platelet','platelet_24h',-1,'log2'),('inr','inr_24h',1,'log2')]:
 v=np.log2(a[col]) if trans=='log2' else a[col]
 mu[name]=float(v.mean());sd[name]=float(v.std(ddof=1));raw.append(sign*(v-mu[name])/sd[name])
r=sum(raw)/4
params={'source_sha256':hashlib.sha256((P/'01_data/processed/mimic_sepsis_endothelial_cohort.csv').read_bytes()).hexdigest(),'source_n':len(a),'component_mean':mu,'component_sample_sd':sd,'composite_mean':float(r.mean()),'composite_sample_sd':float(r.std(ddof=1)),'cutoff_high':0,'external_restandardization':False}
assert np.max(np.abs((r-params['composite_mean'])/params['composite_sample_sd']-a.ceci_z))<1e-12
frozen=O/'frozen_ceci_transform.json'
if frozen.exists():assert json.loads(frozen.read_text())==params
else:frozen.write_text(json.dumps(params,indent=2),encoding='utf-8')
db=duckdb.connect();db.execute("SET threads=6; SET memory_limit='4GB'")
def src(name):
    filename = name + ".csv"
    candidates = (root / filename / filename, root / filename, root / (filename + ".gz"))
    file = next((candidate for candidate in candidates if candidate.is_file()), None)
    if file is None:
        raise FileNotFoundError(f"Missing eICU table: {filename}")
    escaped = file.as_posix().replace("'", "''")
    return f"read_csv_auto('{escaped}',sample_size=20480)"

db.execute(f"CREATE TABLE patient_raw AS SELECT * FROM {src('patient')}")
db.execute("""CREATE TABLE unit_first AS SELECT * EXCLUDE(rn) FROM (SELECT *,row_number() OVER(PARTITION BY patienthealthsystemstayid ORDER BY hospitaladmitoffset DESC,unitvisitnumber,patientunitstayid) rn FROM patient_raw) WHERE rn=1;
CREATE TABLE person_selected AS SELECT * EXCLUDE(rn) FROM (SELECT *,row_number() OVER(PARTITION BY uniquepid ORDER BY sha256('sepsis-external-v1|'||cast(patienthealthsystemstayid AS varchar)),patientunitstayid) rn FROM unit_first WHERE uniquepid IS NOT NULL AND trim(uniquepid)<>'') WHERE rn=1;
""")
db.execute(f"""CREATE TABLE sepsis_flags AS
SELECT patientunitstayid,max(cast(diagnosisoffset BETWEEN -2880 AND 1440 AS INT)) early_code,1 any_code
FROM {src('diagnosis')},unnest(string_split(icd9code,',')) AS tokens(code)
WHERE regexp_matches(upper(replace(trim(code),'.','')),'^(038|A40|A41|R652)') OR upper(replace(trim(code),'.','')) IN ('99591','99592','78552')
GROUP BY patientunitstayid""")
db.execute(f"""CREATE TABLE culture_flags AS SELECT patientunitstayid,1 culture_recorded FROM {src('microLab')} WHERE culturetakenoffset BETWEEN -2880 AND 1440 GROUP BY patientunitstayid""")
db.execute("""CREATE TABLE eligibility AS SELECT p.*,CASE WHEN trim(age)='> 89' THEN 90 ELSE try_cast(age AS DOUBLE) END age_numeric,
coalesce(s.early_code,0) early_sepsis_code,coalesce(s.any_code,0) any_sepsis_code,coalesce(m.culture_recorded,0) culture_recorded,
(lower(coalesce(unitdischargestatus,''))='expired' AND lower(coalesce(hospitaldischargestatus,''))='alive') outcome_contradiction
FROM person_selected p LEFT JOIN sepsis_flags s USING(patientunitstayid) LEFT JOIN culture_flags m USING(patientunitstayid);
CREATE TABLE external_eligible AS SELECT * FROM eligibility WHERE age_numeric>=18 AND unitdischargeoffset>=1440 AND hospitaldischargeoffset>1440 AND early_sepsis_code=1 AND hospitaldischargestatus IN ('Alive','Expired') AND gender IN ('Male','Female') AND NOT outcome_contradiction AND NOT (lower(coalesce(unitdischargestatus,''))='expired' AND unitdischargeoffset<=1440);
""")
flow=[]
queries=[('All ICU unit rows','patient_raw','TRUE'),('First unit per hospital admission','unit_first','TRUE'),('One outcome-blind selected hospital admission per person','person_selected','TRUE'),('Adults','eligibility','age_numeric>=18'),('Adult stays through 24h','eligibility','age_numeric>=18 AND unitdischargeoffset>=1440 AND hospitaldischargeoffset>1440'),('Early coded sepsis among 24h stays','eligibility','age_numeric>=18 AND unitdischargeoffset>=1440 AND hospitaldischargeoffset>1440 AND early_sepsis_code=1'),('Eligible known hospital outcome and sex without contradictory outcome','external_eligible','TRUE')]
for label,table,where in queries:
 n=db.execute(f'SELECT count(*) FROM {table} WHERE {where}').fetchone()[0];flow.append({'Stage':label,'N':n})

pre="age_numeric>=18 AND unitdischargeoffset>=1440 AND hospitaldischargeoffset>1440 AND early_sepsis_code=1"
exclusions=[]
for label,condition in [('Unknown hospital outcome',"hospitaldischargestatus NOT IN ('Alive','Expired') OR hospitaldischargestatus IS NULL"),('Unknown sex after known hospital outcome',"hospitaldischargestatus IN ('Alive','Expired') AND (gender NOT IN ('Male','Female') OR gender IS NULL)"),('ICU death/hospital survival contradiction after known outcome and sex',"hospitaldischargestatus IN ('Alive','Expired') AND gender IN ('Male','Female') AND outcome_contradiction"),('ICU death exactly at 24h after previous checks',"hospitaldischargestatus IN ('Alive','Expired') AND gender IN ('Male','Female') AND NOT outcome_contradiction AND lower(coalesce(unitdischargestatus,''))='expired' AND unitdischargeoffset<=1440")]:
 exclusions.append({'Exclusion':label,'N':db.execute(f'SELECT count(*) FROM eligibility WHERE {pre} AND ({condition})').fetchone()[0]})
pd.DataFrame(exclusions).to_csv(O/'outcome_and_sex_exclusions.csv',index=False)
db.execute(f"""CREATE TABLE marker_events AS SELECT p.patientunitstayid,l.labname,l.labmeasurenamesystem,l.labmeasurenameinterface,l.labresultoffset,l.labresultrevisedoffset,try_cast(l.labresult AS DOUBLE) val,
(l.labresultrevisedoffset<=1440 OR l.labresultrevisedoffset IS NULL) available_by24h
FROM external_eligible p JOIN {src('lab')} l USING(patientunitstayid)
WHERE l.labname IN ('lactate','albumin','platelets x 1000','PT - INR') AND l.labresultoffset BETWEEN -360 AND 1440 AND l.labresultoffset>=p.hospitaladmitoffset AND l.labresultoffset<=p.hospitaldischargeoffset""")
db.execute("""CREATE TABLE score_labs AS SELECT patientunitstayid,
max(val) FILTER(WHERE labname='lactate' AND lower(labmeasurenamesystem)='mmol/l' AND val BETWEEN .1 AND 30 AND available_by24h) lactate,
min(val) FILTER(WHERE labname='albumin' AND lower(labmeasurenamesystem)='g/dl' AND val BETWEEN .5 AND 6.5 AND available_by24h) albumin,
min(val) FILTER(WHERE labname='platelets x 1000' AND lower(labmeasurenamesystem)='k/mcl' AND val BETWEEN 1 AND 1500 AND available_by24h) platelet,
max(val) FILTER(WHERE labname='PT - INR' AND lower(labmeasurenamesystem)='ratio' AND val BETWEEN .5 AND 25 AND available_by24h) inr
FROM marker_events GROUP BY patientunitstayid;
""")
x=db.execute("""SELECT p.patientunitstayid,p.uniquepid,p.patienthealthsystemstayid,p.hospitalid,p.age_numeric AS age,cast(p.gender='Male' AS INTEGER) sex,p.culture_recorded,p.hospitaldischargeoffset,p.unitdischargeoffset,cast(p.hospitaldischargestatus='Expired' AS INTEGER) hospital_mortality,l.lactate,l.albumin,l.platelet,l.inr FROM external_eligible p LEFT JOIN score_labs l USING(patientunitstayid)""").df()
x=x.sort_values('patientunitstayid').reset_index(drop=True)
x['complete_panel']=x[['lactate','albumin','platelet','inr']].notna().all(axis=1)
raw=0
for name,sign in [('lactate',1),('albumin',-1),('platelet',-1),('inr',1)]:
 v=x[name] if name=='albumin' else np.log2(x[name]);raw+=sign*(v-mu[name])/sd[name]/4
x['ceci_frozen']=(raw-params['composite_mean'])/params['composite_sample_sd']
x.to_csv(private/'eicu_eligible_with_frozen_ceci.csv',index=False)
complete=x[x.complete_panel].copy();complete.to_csv(private/'eicu_complete_frozen_ceci.csv',index=False)
flow.append({'Stage':'Complete panel with locked MIMIC transformation','N':len(complete)})
flow.append({'Stage':'Complete panel with recorded early culture sensitivity','N':int(complete.culture_recorded.sum())})
pd.DataFrame(flow).to_csv(O/'cohort_flow.csv',index=False)
quality=db.execute("""SELECT labname,count(*) candidate_records,count(*) FILTER(WHERE NOT available_by24h) revised_after24h,count(*) FILTER(WHERE labresultrevisedoffset IS NULL) revision_offset_missing,count(*) FILTER(WHERE val IS NULL) numeric_missing FROM marker_events GROUP BY labname""").df()
quality.to_csv(O/'laboratory_availability_checks.csv',index=False)
summ=[]
for group,z in [('All eligible',x),('Complete panel',complete),('Incomplete panel',x[~x.complete_panel]),('Complete plus culture record',complete[complete.culture_recorded==1])]:
 summ.append({'Population':group,'N':len(z),'Hospitals':z.hospitalid.nunique(),'Hospital_deaths':int(z.hospital_mortality.sum()),'Hospital_mortality':float(z.hospital_mortality.mean()) if len(z) else None,'Age_median':float(z.age.median()) if len(z) else None,'Male_fraction':float(z.sex.mean()) if len(z) else None,'Recorded_culture_fraction':float(z.culture_recorded.mean()) if len(z) else None})
pd.DataFrame(summ).to_csv(O/'population_summary.csv',index=False)
dis=[]
for name in ['lactate','albumin','platelet','inr','ceci_frozen']:
 vals=complete[name];dis.append({'variable':name,'N':int(vals.notna().sum()),'mean':float(vals.mean()),'sd':float(vals.std()),'p05':float(vals.quantile(.05)),'p25':float(vals.quantile(.25)),'median':float(vals.median()),'p75':float(vals.quantile(.75)),'p95':float(vals.quantile(.95))})
pd.DataFrame(dis).to_csv(O/'marker_and_frozen_score_distributions.csv',index=False)
assert x.patientunitstayid.is_unique and x.uniquepid.is_unique and x.patienthealthsystemstayid.is_unique
assert np.isfinite(complete.ceci_frozen).all()
result={'endpoint':'In-hospital mortality after 24h; not complete 28d vital status','eligible_n':len(x),'complete_n':len(complete),'complete_deaths':int(complete.hospital_mortality.sum()),'complete_hospitals':int(complete.hospitalid.nunique()),'complete_recorded_culture_n':int(complete.culture_recorded.sum()),'high_frozen_ceci_fraction':float((complete.ceci_frozen>=0).mean()),'eligible_to_complete_fraction':len(complete)/len(x),'frozen_transform_sha256':hashlib.sha256(frozen.read_bytes()).hexdigest(),'primary_28d_model_validation_possible':False,'score_association_assessment_possible':bool(len(complete)>=200 and complete.hospital_mortality.sum()>=20 and complete.hospitalid.nunique()>=10),'selection':'one outcome-blind selected hospitalization per person; first ICU unit within that hospitalization','patient_derivatives':'restricted_local_only; not for release'}
(O/'extraction_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2));print(pd.DataFrame(flow).to_string(index=False));print(pd.DataFrame(summ).to_string(index=False))
