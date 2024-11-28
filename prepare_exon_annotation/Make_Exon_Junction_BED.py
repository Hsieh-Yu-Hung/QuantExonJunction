import pandas as pd
import argparse
from dataclasses import dataclass

@dataclass
class record:
    chromosome: str
    start: int
    end: int
    transcript_id: str
    gene_name: str
    annotation: str

# 預期 GTF 的格式
EXPECTED_GTF_COLUMNS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# 解析命令列參數
def parse_args():
    parser = argparse.ArgumentParser(description="Make exon junction bed file")
    parser.add_argument("-g", "--gtf", required=True, help="Input Selected GTF file")
    parser.add_argument("-o", "--output", required=True, help="Output exon-junction BED file")
    return parser.parse_args()

# 解析 GTF 的 attribute 欄位
def parse_attribute(attribute:str) -> dict:
    att = attribute.split('; ')
    att_dict = {}
    for item in att:
        key, value = item.split(' ')
        att_dict[key] = value.strip('"')
    return att_dict

# 讀取 GTF 檔案，並解析出 exon 的資訊
def read_gtf_exon(gtf_file:str) -> pd.DataFrame:
    df = pd.read_csv(gtf_file, sep="\t", header=None, names=EXPECTED_GTF_COLUMNS)
    df['transcript_id'] = df['attribute'].apply(lambda x: parse_attribute(x)['transcript_id'])
    df['gene_name'] = df['attribute'].apply(lambda x: parse_attribute(x)['gene_name'])
    df['exon_number'] = df['attribute'].apply(lambda x: int(parse_attribute(x)['exon_number']) if 'exon_number' in parse_attribute(x) else -1)
    return df[df['feature'] == 'exon']

# 製作 exon-junction bed
def make_exon_junction_bed(df:pd.DataFrame, target_transcript_id:str) -> pd.DataFrame:
    df = df[df['transcript_id'] == target_transcript_id]
    df_len = df.shape[0]
    records = []    
    for i in range(df_len - 1):
        rec = record(
            chromosome=df.iloc[i]['chr'],
            start=df.iloc[i]['end'],
            end=df.iloc[i+1]['start'],
            transcript_id=target_transcript_id,
            gene_name=df.iloc[i]['gene_name'],
            annotation=f"exon{df.iloc[i]['exon_number']}-{df.iloc[i+1]['exon_number']}"
        )
        records.append(rec)
    tdf = pd.DataFrame(records)
    return tdf

# 主程式
def main():
    # 解析命令列參數
    args = parse_args()

    # 讀取 GTF 檔案
    df = read_gtf_exon(args.gtf)

    # 讀取 transcript list
    tList = df['transcript_id'].unique()

    # 製作 exon-junction bed
    tdf_list = list(map(lambda t: make_exon_junction_bed(df, t), tList))
    tdf = pd.concat(tdf_list)
    tdf = tdf[tdf['chr'].str.startswith('chr')]
    tdf.to_csv(args.output, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
    