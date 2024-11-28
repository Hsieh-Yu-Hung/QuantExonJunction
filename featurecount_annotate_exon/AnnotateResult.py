import pandas as pd
import argparse
import os

# 中繼檔案名稱
INTERMED_FILE = "annotated_result.tsv"
EXPECTED_COLUMNS = [
    'chr','start','end','transcript_id','gene_name','exon_junction', 'Panel_Target',
    'fc-chr','fc-start','fc-end','fc-transcript_id','fc-count','exon_span'
]

# 解析參數
def parse_args():
    parser = argparse.ArgumentParser(description="Annotate featurecount result table and clean up.")
    parser.add_argument("-b", "--bed", required=True, help="Input featurecount junction-ExonSpan BED file.")
    parser.add_argument("-a", "--annotation", required=True, help="Input Exon junction annotated BED file.")
    parser.add_argument("-o", "--output", required=True, help="Output Annotated featurecount result table, Excel format.")
    parser.add_argument("-k", "--keep", action="store_true", help=f"Keep the intermeiate file as \"{INTERMED_FILE}\" ")
    return parser.parse_args()

# 用 bedtool intersect 做交集
def insetsect(junction_bed:str, junction_anno:str):
    command = (
        f"sed -i -n '/^chr/p' {junction_bed} && "
        f"bedtools intersect -b {junction_bed} -a {junction_anno} -loj > {INTERMED_FILE}"
    )
    os.system(command)

# 讀取中繼檔案做資料清洗
def cleanUpAnnotation(intermed_file:str):
    try:
        df = pd.read_csv(intermed_file, sep='\t', header=None, names=EXPECTED_COLUMNS)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error reading file: {e}\nFile: {intermed_file}")
        exit(1)

    # 將 '.' 替換為 -1，然後將 'exon_span' 轉換為數值型別
    df['exon_span'] = df['exon_span'].replace('.', -1).astype(int)

    # 保留原始 '.' 和等於 1 的資料
    df = df[df['exon_span'] <= 1]

    # 合併 gene_name 和 exon_junction 相同的資料，並將 fc-count 欄位相加
    df_merged = df.groupby(['chr','start','end','transcript_id','gene_name','exon_junction', 'Panel_Target', 'fc-chr'], as_index=False).agg({
        'fc-start':'min',
        'fc-end':'max',
        'fc-transcript_id':'first',
        'fc-count': 'sum',
        'exon_span':'first'
    })

    # 將 fc-count 欄位中的 "." 換成 0
    df_merged['fc-count'] = df_merged['fc-count'].replace('.', 0).astype(int)

    # 將 df_merged 中的 "." 和 -1 替換為 "NA"
    df_merged = df_merged.replace({'.': 'NA', -1: 'NA'})
    
    # 回傳結果
    return df_merged

# 主程式
def main():

    # 解析輸入參數
    args = parse_args()

    # 執行 bedtool intersect
    insetsect(args.bed, args.annotation)

    # 資料清洗
    cleanup_data = cleanUpAnnotation(INTERMED_FILE)

    # 輸出
    cleanup_data.to_excel(args.output,index=False)

    # 清除中繼檔
    if not args.keep:
        os.system(f"rm {INTERMED_FILE}")

if __name__ == "__main__":
    main()
