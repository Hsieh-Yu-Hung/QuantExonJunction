import pandas as pd
import argparse

# 預期 GTF 的格式
EXPECTED_GTF_COLUMNS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# 解析參數
def parse_args():
    parser = argparse.ArgumentParser(description="Annotate Exon span distance.")
    parser.add_argument("-b", "--bed", required=True, help="Input featurecount junction quantification BED file.")
    parser.add_argument("-g", "--gtf", required=True, help="Input genome GTF for annotating exon.")
    return parser.parse_args()

# 解析 GTF 的 attribute 欄位
def parse_attribute(attribute:str) -> dict:
    att = attribute.split('; ')
    att_dict = {}
    for item in att:
        key, value = item.split(' ')
        att_dict[key] = value.strip('"')
    return att_dict

# 讀取 junction count bed 檔
def read_junctions(junction_bed_file:str):
    print(" --> Read and Parsing Junction BED...")
    jcount_df = pd.read_csv(junction_bed_file, sep="\t", header=None)
    column_names = ["chromosome", "start", "end", "transcript_id", "Count"]
    jcount_df.columns = column_names
    return jcount_df

# 讀取 GTF 檔
def read_gtf(gtf_file:str):
    print(" --> Read and Parsing GTF...")
    gtf_df = pd.read_csv(gtf_file, sep='\t', header=None, comment='#', names=EXPECTED_GTF_COLUMNS)
    gtf_df['transcript_id'] = gtf_df['attribute'].apply(lambda x: parse_attribute(x)['transcript_id'])
    gtf_df['gene_name'] = gtf_df['attribute'].apply(lambda x: parse_attribute(x)['gene_name'])
    gtf_df['exon_number'] = gtf_df['attribute'].apply(lambda x: parse_attribute(x)['exon_number'] if 'exon_number' in parse_attribute(x) else "NA")
    return gtf_df

# 搜索 exon 的位置
def search_exon_index(df, position, to_search):
    for _, row in df.iterrows():
        if position <= row[to_search]:
            return row['exon_number']
    return -1

# 計算 Exon-span 距離
def calculate_exon_span(row, gtf_df):
    stid = row['transcript_id']
    start = row['start']
    end = row['end']
    if stid in gtf_df['transcript_id'].values:
        gtf_subdf = gtf_df[gtf_df['transcript_id'] == stid]
        gtf_subdf = gtf_subdf[gtf_subdf['feature'] == 'exon']
        start_exon_number = search_exon_index(gtf_subdf, start, 'start')
        end_exon_number = search_exon_index(gtf_subdf, end, 'end')
        exon_span = abs(int(end_exon_number) - int(start_exon_number)) + 1
        print(f"{stid} exon span distance = {exon_span}.")
        return exon_span
    else:
        return "."

# 主程式
def main():

    # 解析命令列參數
    args = parse_args()

    # 讀取 junction
    jcount_df = read_junctions(args.bed)

    # 讀取 GTF
    gtf_df = read_gtf(args.gtf)

    # 使用 apply 方法計算 'exon_span'
    jcount_df['exon_span'] = jcount_df.apply(calculate_exon_span, axis=1, gtf_df=gtf_df)

    # 輸出
    output_name = args.bed.replace('.bed', '_ExonSpan.bed')
    jcount_df.to_csv(output_name, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()