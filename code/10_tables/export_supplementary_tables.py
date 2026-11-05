from pathlib import Path
import csv
from openpyxl import load_workbook
root=Path(__file__).resolve().parents[2]
source=root/"data_processed/supplementary_tables/supplementary_tables_s1_s9.xlsx"
out=root/"data_processed/supplementary_tables/csv"
out.mkdir(parents=True, exist_ok=True)
wb=load_workbook(source,data_only=False,read_only=True)
for ws in wb.worksheets:
    with (out/(ws.title+".csv")).open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.writer(f)
        for row in ws.iter_rows(values_only=True): w.writerow(["" if v is None else v for v in row])
wb.close()
