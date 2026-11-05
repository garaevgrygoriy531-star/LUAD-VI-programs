#!/usr/bin/env python3
"""Validate the processed LUAD repository files and numerical outputs."""
from __future__ import annotations
import csv, gzip, hashlib, math, re, sys
from collections import Counter
from pathlib import Path
from openpyxl import load_workbook

ROOT=Path(__file__).resolve().parents[1]
checks=[]
def require(ok,message):
    if not ok: raise AssertionError(message)
    checks.append(message)
def rows(path):
    file=ROOT/path
    with file.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t" if file.suffix.lower()==".tsv" else ","))
def sha(path):
    h=hashlib.sha256()
    with (ROOT/path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()
def close(a,b,tol):
    return abs(float(a)-float(b)) <= float(tol)
def boolval(x): return str(x).strip().lower()=="true"

# Core files and directory structure
required=["README.md","LICENSE-CODE","docs/REPRODUCIBILITY.md","docs/DATA_AVAILABILITY.md",
"docs/DATA_PROVENANCE.md","docs/DATA_LICENSE_AND_SOURCES.md","metadata/datasets.csv","metadata/program_definitions.csv",
"metadata/sample_manifest.csv","metadata/column_dictionary.csv","metadata/provenance_manifest.csv",
"metadata/figure_table_map.csv","metadata/panel_manifest.csv","metadata/analysis_dependency_graph.csv",
"data_processed/pathology/gse308103_stage_scores.tsv","data_processed/pathology/gse319666_sample_scores.tsv",
"data_processed/pathology/gse319666_stage_ia_membership.tsv","data_processed/cross_context/gse282617_sample_scores.csv",
"data_processed/supplementary_tables/supplementary_tables_s1_s9.xlsx","tests/expected_results.csv",
"results/stage_architecture/model_selection_bootstrap_recomputed_summary.csv",
"results/stage_architecture/model_selection_bootstrap_recomputed_replicates.csv"]
for rel in required: require((ROOT/rel).is_file(),f"required file exists: {rel}")
if (ROOT/"CITATION.cff").exists():
    citation=(ROOT/"CITATION.cff").read_text(encoding="utf-8-sig")
    require(all(x in citation for x in ["cff-version:","message:","title:","authors:"]),"citation file includes required CFF fields")
else:
    checks.append("citation metadata pending author names and identifiers")
require(not (ROOT/"results_reference").exists(),"obsolete result-reference directory is absent")
for p in ROOT.rglob("*"):
    if p.is_file(): require(p.stat().st_size>0,f"nonempty file: {p.relative_to(ROOT).as_posix()}")

# Current figure files, inputs, and scripts
fig_names={"Figure_1.pdf","Figure_2.pdf","Figure_3.pdf","Figure_4.pdf","Figure_5.pdf",
"Supplementary_Figure_S1.pdf","Supplementary_Figure_S2.pdf","Supplementary_Figure_S3.pdf",
"Supplementary_Figure_S4.pdf","Supplementary_Figure_S5.pdf"}
require({p.name for p in (ROOT/"figures").iterdir() if p.is_file()}==fig_names,"the figures folder contains exactly the ten current PDFs")
canonical=rows("tests/canonical_figure_sha256.csv")
require(len(canonical)==10,"canonical figure manifest contains ten entries")
for entry in canonical:
    name=entry["repository_file"].split("/")[-1]
    require(boolval(entry["byte_identical"]) and entry["repository_sha256"]==entry["source_sha256"],f"canonical source identity is recorded: {name}")
    require(sha(entry["repository_file"])==entry["repository_sha256"],f"canonical figure hash matches manifest: {name}")

provenance=rows("metadata/provenance_manifest.csv")
require(bool(provenance) and set(provenance[0])=={"public_file","source_accession","source_object","transformation_summary","related_script","downstream_use","inferential_unit"},"provenance manifest uses the documented field names")
require("Public analysis or rendering script directly associated with the file; this field does not necessarily denote the original upstream generator." in (ROOT/"docs/DATA_PROVENANCE.md").read_text(encoding="utf-8-sig"),"related-script field is defined in data-provenance documentation")
for r in provenance:
    require((ROOT/r["public_file"]).exists(),f"provenance file exists: {r['public_file']}")
    if r["related_script"].startswith("code/"):
        require((ROOT/r["related_script"]).is_file(),f"related script exists: {r['related_script']}")
for r in rows("metadata/program_definitions.csv"):
    if r["definition_file"]:
        require((ROOT/r["definition_file"]).exists(),f"program definition exists: {r['definition_file']}")
    require(bool(r.get("definition_source")),f"program source is identified: {r['program']}")
for r in rows("metadata/processed_data_catalog.csv"):
    require((ROOT/r["directory"]).is_dir(),f"processed-data directory exists: {r['directory']}")

panel=rows("metadata/panel_manifest.csv")
mapping=rows("metadata/figure_table_map.csv")
panels={x["panel_id"] for x in panel}
mapped={x["panel_or_table"] for x in mapping if x["panel_or_table"] in panels}
require(len(panel)==41 and panels==mapped,"panel manifest contains all 41 current panels")
for x in mapping:
    require((ROOT/x["output_file"]).exists(),f"figure/table output exists: {x['output_file']}")
    if x["generating_script"].startswith("code/"):
        require((ROOT/x["generating_script"]).is_file(),f"generating script exists: {x['generating_script']}")
    else:
        require(x["panel_or_table"]=="Table 1" and "Main Table 1 document" in x["generating_script"],"Main Table 1 transcription source is identified")
    for source in x["processed_inputs"].split(";"):
        source=source.strip()
        require((ROOT/source).exists(),f"mapped processed input exists: {source}")
plot_files={p.relative_to(ROOT).as_posix() for p in (ROOT/"data_processed/plot_ready").iterdir() if p.is_file()}
mapped_plot={s.strip() for x in mapping for s in x["processed_inputs"].split(";") if s.strip().startswith("data_processed/plot_ready/")}
require(plot_files.issubset(mapped_plot),f"every plot-ready file is mapped as an active input ({len(plot_files)} files)")

# Program and cohort constants
defs=rows("metadata/program_definitions.csv")
program_counts={x["program"]:int(x["member_count"]) for x in defs}
require(program_counts=={"V1":115,"V2":37,"V3":182,"V4":140,"ORACLE":23},"published program sizes")
for name,count in {"V1":115,"V2":37,"V3":182,"V4":140}.items():
    genes=[x.strip() for x in (ROOT/"metadata/program_members"/f"{name}.txt").read_text(encoding="utf-8-sig").splitlines() if x.strip()]
    require(len(genes)==count and len(set(genes))==count,f"{name} source gene membership count and uniqueness")
oracle_genes=[x.strip() for x in (ROOT/"metadata/program_members/ORACLE.txt").read_text(encoding="utf-8-sig").splitlines() if x.strip()]
require(len(oracle_genes)==23 and len(set(oracle_genes))==23,"ORACLE has 23 unique published gene symbols")
require(not (ROOT/"data_processed/program_scores/oracle_coefficients.tsv").exists(),"ORACLE coefficient table is not redistributed")
require(not (ROOT/"data_processed/proteomics/pxd046637_measured_protein_gene_ranks.tsv").exists(),"article-derived PXD046637 gene ranks are not redistributed")
manifest=rows("metadata/sample_manifest.csv")
for ds,n,patients in [("GSE308103",75,23),("GSE319666",177,40)]:
    subset=[x for x in manifest if x["dataset"]==ds]
    require(len(subset)==n and len({x["patient_id"] for x in subset})==patients,f"{ds} source specimens and patients")
x=[r for r in manifest if r["dataset"]=="GSE319666"]
primary=[r for r in x if boolval(r["included_in_primary_contrast"])]
paired=[r for r in x if boolval(r["included_in_paired_contrast"])]
ia=[r for r in x if boolval(r["included_in_stage_IA_analysis"])]
require((len(primary),len({r["patient_id"] for r in primary}))==(129,38),"GSE319666 primary contrast membership")
require(len({r["patient_id"] for r in paired})==21,"GSE319666 paired analysis includes 21 patients")
require((len(ia),len({r["patient_id"] for r in ia}))==(92,33),"GSE319666 Stage IA membership")
gse102511=rows("data_processed/cross_context/gse102511_sample_scores.tsv")
source_groups=Counter(r["source_label"].rsplit("-",1)[-1] for r in gse102511)
paired_patients=set()
for patient in {r["patient_id"] for r in gse102511}:
    labels={r["source_label"].rsplit("-",1)[-1] for r in gse102511 if r["patient_id"]==patient}
    if {"AAH","LUAD"}.issubset(labels):paired_patients.add(patient)
require((len(gse102511),len({r["patient_id"] for r in gse102511}),source_groups)==(48,17,{"AAH":17,"LUAD":16,"NL":15}),"GSE102511 source specimens, patients, and source-label groups")
gse102511_meta=[r for r in rows("metadata/datasets.csv") if r["accession"]=="GSE102511"]
require(len(gse102511_meta)==1 and (int(gse102511_meta[0]["analysis_sample_count"]),int(gse102511_meta[0]["analysis_patient_count"]))==(32,16),"GSE102511 primary paired analysis metadata")
require((len(paired_patients),2*len(paired_patients))==(16,32),"GSE102511 primary AAH-LUAD comparison includes 16 pairs")
vi=rows("data_processed/pathology/gse273528_vi_samples.csv")
require((len(vi),len({r["patient_id"] for r in vi}),sum(int(r["VI"]) for r in vi))==(59,58,17),"GSE273528 lesion, patient, and VI-positive counts")
with gzip.open(ROOT/"data_processed/cellular/gse308103_nucleus_annotations.tsv.gz","rt",encoding="utf-8") as f:
    nuclei=sum(1 for _ in f)-1
require(nuclei==798100,"GSE308103 nucleus annotation count")
outcome=rows("data_processed/outcome/stageI_primary_cohort_scores.csv")
surv=[r for r in outcome if r["cohort"] in {"GSE50081","GSE31210","GSE30219"} and boolval(r["primary_eligible"])]
require((len(surv),sum(int(r["os_event"]) for r in surv))==(323,78),"primary stage I survival cohort includes 323 patients and 78 deaths")
cell=rows("data_processed/cellular/scrna_sample_manifest.csv")
require(len({r["patient_id"] for r in cell})==9,"GSE189357 cellular replication includes nine patients")

# Manuscript-facing numerical regression checks
routes={
"figure3a":("results/pathology/v1_adjusted_stage_models.csv","dataset","outcome","adjustment","estimate","bootstrap_ci_low","bootstrap_ci_high","n_specimens","n_patients"),
"v1_refit":("results/pathology/v1_adjusted_stage_models.csv","dataset","outcome","adjustment","estimate",None,None,"n_specimens","n_patients"),
"stage_ia":("results/stage_architecture/stage_ia_contrasts.csv",None,"program",None,"estimate","ci_low","ci_high","n_specimens","n_patients"),
"cross_context":("results/cross_context/pathology_contrasts.csv","dataset","program",None,"estimate","ci_low","ci_high","n_samples","n_patients"),
"os":("results/stage_i_survival/pooled_os_reml_hk.csv",None,"score","model","hr","ci_low","ci_high",None,None),
"vi":("results/vascular_invasion/primary_firth_models.csv",None,None,"model","OR","CI_low","CI_high","n","patients")}
for e in rows("tests/expected_results.csv"):
    kind=e["analysis"]; rel,dscol,pgcol,modcol,estcol,lowcol,highcol,ncol,pcol=routes[kind]
    table=rows(rel); found=[]
    for r in table:
        if dscol and kind not in {"os","vi"} and r.get(dscol)!=e["dataset"]:continue
        if kind=="os" and r.get("cohorts")!=e["dataset"]:continue
        if kind=="vi" and r.get("model")!=e["model"]:continue
        if pgcol and r.get(pgcol)!=e["program"]:continue
        if modcol and r.get(modcol)!=e["model"]:continue
        found.append(r)
    require(len(found)==1,f"one result row matches {kind} {e['dataset']} {e['program']} {e['model']}")
    r=found[0]; tol=float(e["tolerance"])
    require(close(r[estcol],e["estimate"],tol),f"{kind} {e['dataset']} {e['program']} estimate")
    if lowcol and e["ci_low"]:require(close(r[lowcol],e["ci_low"],tol) and close(r[highcol],e["ci_high"],tol),f"{kind} {e['dataset']} {e['program']} interval")
    if ncol and e["n"]:require(int(float(r[ncol]))==int(e["n"]),f"{kind} {e['dataset']} {e['program']} specimen count")
    if pcol and e["patients"]:require(int(float(r[pcol]))==int(e["patients"]),f"{kind} {e['dataset']} {e['program']} patient count")

# Figure 3A bootstrap results match the supplied panel values.
fig3=rows("data_processed/plot_ready/F3a_integrated.csv")
refit=rows("results/pathology/v1_adjusted_stage_models.csv")
require(len(refit)==12 and all(int(r["bootstrap_replicates"])==1000 and int(r["seed"])==20261016 for r in refit),"Figure 3A results include all 12 contrasts with the stated resampling settings")
for r in refit:
    p=[x for x in fig3 if x["dataset"]==r["dataset"] and x["program"]==r["outcome"] and x["series"]==r["adjustment"]]
    require(len(p)==1,f"one Figure 3A panel row matches {r['dataset']} {r['outcome']} {r['adjustment']}")
    require(close(r["estimate"],p[0]["estimate"],1e-8) and close(r["bootstrap_ci_low"],p[0]["ci_low"],1e-8) and close(r["bootstrap_ci_high"],p[0]["ci_high"],1e-8),f"Figure 3A bootstrap result matches panel input: {r['dataset']} {r['outcome']} {r['adjustment']}")

# Current stage model outputs reproduce the stored numerical targets.
target=rows("results/numeric_reference/ARCHITECTURE_PRIMARY.csv")
fitted=rows("results/stage_architecture/primary_estimates.csv")
for t in target:
    z=[r for r in fitted if r["dataset"]==t["dataset"] and r["scenario"]==t["scenario"] and r["program"]==t["program"]]
    require(len(z)==1,f"one stage-model fit for {t['dataset']} {t['scenario']} {t['program']}")
    for col in ["beta_pre","PRE_SPAN","INV_INCREMENT","ALIGNED_PRE_SPAN","ALIGNED_INV_INCREMENT","STEP_MINUS_SPAN_CONTRAST","AICc_hinge","AICc_step","AICc_gradual"]:
        require(close(z[0][col],t[col],1e-6),f"stage-model numeric target {col} {t['dataset']} {t['program']}")

# Recomputed bootstrap outputs preserve the expected model-preference pattern.
boot_summary=rows("results/stage_architecture/model_selection_bootstrap_recomputed_summary.csv")
boot_replicates=rows("results/stage_architecture/model_selection_bootstrap_recomputed_replicates.csv")
require(len(boot_summary)==18,"recomputed bootstrap summary contains all dataset, program, and model combinations")
require(len(boot_replicates)==18000 and len({int(r["replicate"]) for r in boot_replicates})==1000,"recomputed bootstrap replicates contain 1,000 resamples per dataset and program")
for dataset in ["GSE308103","GSE319666"]:
    for program in ["V2","V3","V4"]:
        group=[r for r in boot_summary if r["dataset"]==dataset and r["program"]==program]
        require({r["model"] for r in group}=={"Gradual","Invasive step","Hinge"},"recomputed bootstrap includes all three model classes")
        require(close(sum(float(r["selected"]) for r in group),1,1e-12),"recomputed bootstrap selection fractions sum to one")
        winner=max(group,key=lambda r:float(r["selected"]))
        require(winner["model"]=="Invasive step","recomputed bootstrap preserves the model-preference pattern")

# Supplementary sheet CSVs match their workbook worksheets.
wb=load_workbook(ROOT/"data_processed/supplementary_tables/supplementary_tables_s1_s9.xlsx",read_only=True,data_only=False)
require(len(wb.sheetnames)==9,"supplementary workbook has nine sheets")
def norm(v):
    if v is None:return ""
    if isinstance(v,str):
        v=v.strip()
        if not v:return ""
        try:return ("num",float(v))
        except ValueError:return ("str",v)
    if isinstance(v,(int,float)):return ("num",float(v))
    return ("str",str(v).strip())
for ws in wb.worksheets:
    rel=Path("data_processed/supplementary_tables/csv")/(ws.title+".csv")
    with (ROOT/rel).open(encoding="utf-8-sig",newline="") as f: actual=list(csv.reader(f))
    expected=[[norm(v) for v in row] for row in ws.iter_rows(values_only=True)]
    actual_norm=[[norm(v) for v in row] for row in actual]
    while expected and all(x==("",) for x in expected[-1]):expected.pop()
    while actual_norm and all(x==("",) for x in actual_norm[-1]):actual_norm.pop()
    require(actual_norm==expected,f"exported worksheet matches: {ws.title}")
wb.close()

# Public text and tables contain no local paths or exposed credentials.
private=re.compile(r"(?:[A-Z]:\\\\Users\\\\|[A-Z]:\\\\[^\\\\\\s]+|/home/[^\\s]+|BEGIN [A-Z ]*PRIVATE KEY|api[_-]?key\\s*[:=]|password\\s*[:=]|token\\s*[:=])",re.I)
scan_ext={".md",".txt",".csv",".tsv",".R",".r",".py",".cff",".yml",".yaml",".json"}
for base in ["README.md","docs","metadata","code","environment","results","data_processed"]:
    paths=[ROOT/base] if (ROOT/base).is_file() else (ROOT/base).rglob("*")
    for p in paths:
        if not p.is_file() or p.suffix.lower() not in scan_ext:continue
        text=p.read_text(encoding="utf-8-sig",errors="ignore")
        require(not private.search(text),f"local-path and credential check: {p.relative_to(ROOT).as_posix()}")
# Workbook cells are checked for local paths and exposed credentials as well.
wb=load_workbook(ROOT/"data_processed/supplementary_tables/supplementary_tables_s1_s9.xlsx",read_only=True,data_only=False)
for ws in wb.worksheets:
    for row in ws.iter_rows(values_only=True):
        for value in row:
            if isinstance(value,str):
                require(not private.search(value),f"workbook privacy check: {ws.title}")
wb.close()


# Figure 3C and Supplementary Figure S4D inputs are derived from supplied results.
f3a=rows("data_processed/plot_ready/F3a_integrated.csv")
f3c=rows("data_processed/plot_ready/F3c.csv")
require(len(f3c)==6 and {r["program"] for r in f3c}=={"V2","V3","V4"},"Figure 3C includes six cohort-program attenuation values")
for r in f3c:
    pair=[x for x in f3a if x["dataset"]==r["dataset"] and x["program"]==r["program"]]
    unadjusted=next(x for x in pair if x["series"]=="Unadjusted")
    adjusted=next(x for x in pair if x["series"]=="V1-adjusted")
    value=100*(1-abs(float(adjusted["estimate"]))/abs(float(unadjusted["estimate"])))
    require(close(value,r["attenuation_percent"],0.051),f"Figure 3C attenuation derives from supplied coefficients: {r['dataset']} {r['program']}")
require([r["attenuation_percent"] for r in f3c]==["50.1","52.1","66.6","44.7","73.3","66.2"],"Figure 3C values match the manuscript-reported reductions")

s4d=rows("data_processed/plot_ready/S4d.csv")
s4res=rows("results/cellular/v1_adjusted_decomposition_sensitivity.csv")
require(len(s4d)==16 and len({(r["lineage"],r["program"],r["component"]) for r in s4d})==16,"S4D contains four lineages by two programs by two components")
require(all(int(r["bootstrap_replicates"])==5000 for r in s4d),"S4D confidence intervals use 5,000 bootstrap replicates")
for r in s4d:
    ref=next(x for x in s4res if x["lineage"]==r["lineage"] and x["program"]==r["program"])
    stem="V1_adjusted_within" if r["component"]=="Within-cell-type" else "V1_adjusted_composition"
    require(close(r["estimate"],ref[stem],1e-12),f"S4D estimate matches result table: {r['lineage']} {r['program']} {r['component']}")
    require(close(r["ci_low"],ref[stem+"_CI_low"],1e-12) and close(r["ci_high"],ref[stem+"_CI_high"],1e-12),f"S4D interval endpoints match result table: {r['lineage']} {r['program']} {r['component']}")
require(all(boolval(r["direction_preserved"]) and boolval(r["within_gt_composition"]) for r in s4res),"V1-adjusted total direction and component ordering are preserved in all eight comparisons")
patient_only=rows("results/cellular/v1_adjusted_decomposition_patient_only_robustness.csv")
require(len(patient_only)==1 and int(patient_only[0]["comparison_count"])==8,"patient-only sensitivity summarizes eight comparisons")
require(int(patient_only[0]["total_direction_preserved"])==8 and int(patient_only[0]["within_cell_type_larger_than_composition"])==8,"patient-only sensitivity preserves direction and component ordering in 8/8 comparisons")
require(int(patient_only[0]["minimum_valid_bootstrap_resamples"])==5000,"patient-only sensitivity records the bootstrap count")
s8_text=" ".join(cell for row in rows("data_processed/supplementary_tables/csv/S8.csv") for cell in row.values())
require("F. V1-adjusted decomposition sensitivity" in s8_text and "V1-adjusted" in s8_text,"Supplementary Table S8 includes the V1-adjusted section")

s8grid=list(csv.reader((ROOT/"data_processed/supplementary_tables/csv/S8.csv").open(encoding="utf-8-sig",newline="")))
section=next(i for i,row in enumerate(s8grid) if row and row[0]=="F. V1-adjusted decomposition sensitivity")
header=s8grid[section+1]
frows=s8grid[section+2:section+10]
require(len(frows)==8 and header[:6]==["Lineage","Program","Total effect (95% CI)","Within-cell-type program-score component (95% CI)","Fine-cell-type composition component (95% CI)","Covariance component (95% CI)"],"S8 Section F has the complete decomposition layout")
lineage_names={"epithelial":"Epithelial","fibroblast_stromal":"Fibroblast/stromal","endothelial":"Endothelial","myeloid":"Myeloid"}
for tabrow in frows:
    src=next(r for r in s4res if r["lineage"]==next(k for k,v in lineage_names.items() if v==tabrow[0]) and r["program"]==tabrow[1])
    specs=[("V1_adjusted_total","V1_adjusted_total_CI_low","V1_adjusted_total_CI_high"),
           ("V1_adjusted_within","V1_adjusted_within_CI_low","V1_adjusted_within_CI_high"),
           ("V1_adjusted_composition","V1_adjusted_composition_CI_low","V1_adjusted_composition_CI_high"),
           ("V1_adjusted_covariance","V1_adjusted_covariance_CI_low","V1_adjusted_covariance_CI_high")]
    for cell,(est,low,high) in zip(tabrow[2:6],specs):
        values=[float(x) for x in re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)",cell)]
        expected=[round(float(src[est]),3),round(float(src[low]),3),round(float(src[high]),3)]
        require(len(values)==3 and values==expected,f"S8 Section F values match supplied result: {tabrow[0]} {tabrow[1]} {est}")

coverage=[r for r in rows("data_processed/outcome/program_coverage.csv") if "GSE68465" in str(r.values())]
gse68465=[r for r in rows("metadata/datasets.csv") if r["accession"]=="GSE68465"]
s2_text=(ROOT/"data_processed/supplementary_tables/csv/S2.csv").read_text(encoding="utf-8-sig")
require(bool(coverage) and len(gse68465)==1 and "coverage screening only" in gse68465[0]["role_in_this_study"].lower(),"GSE68465 remains documented for coverage screening")
require("GSE68465" not in s2_text,"GSE68465 is absent from the outcome-coverage table snapshot")
s9_text=(ROOT/"data_processed/supplementary_tables/csv/S9.csv").read_text(encoding="utf-8-sig")
require("Supplementary Table S9. External cellular concordance and spatial support" in s9_text and "Fibroblast/stromal" in s9_text,"S9 title and lineage terminology match the manuscript")
s8_text_lower=s8_text.lower()
require("cell-count weighting" in s8_text_lower and "alternative fine-type scores" in s8_text_lower,"S8 sensitivity labels match the manuscript")
main_table=rows("data_processed/supplementary_tables/main_table_1.csv")
require(len(main_table)==8 and main_table[0]["Evidence layer"]=="Discovery","Main Table 1 matches the accompanying manuscript table structure")

require(any(r["panel_or_table"]=="S4d" for r in mapping) and any(r["panel_or_table"]=="Table S8" and "v1_adjusted_decomposition_sensitivity.csv" in r["processed_inputs"] for r in mapping),"S4D-to-S8 provenance mapping exists")

# Public terminology and filenames follow the manuscript and figure labels.
blocked=["Trans"+"ition-dominant","Pre-"+"LUAD","Dom"+"inance","Within-type/"+"state","Early-"+"lesion","Invasive "+"increment"]
for p in ROOT.rglob("*"):
    if p.is_file() and "dom"+"inance" in p.name.lower():
        require(False,f"obsolete figure filename remains: {p.relative_to(ROOT).as_posix()}")
text_ext={".md",".txt",".csv",".tsv",".R",".r",".py",".cff",".yml",".yaml",".json"}
for p in ROOT.rglob("*"):
    if not p.is_file() or p.suffix.lower() not in text_ext:continue
    txt=p.read_text(encoding="utf-8-sig",errors="ignore")
    for term in blocked:
        require(not re.search(re.escape(term),txt,re.I),f"public terminology check: {p.relative_to(ROOT).as_posix()}")

print(f"PASS: {len(checks)} repository checks")
if not (ROOT/"CITATION.cff").exists():
    print("NOTE: repository citation metadata are pending author names, author identifiers, and an archival DOI.")
