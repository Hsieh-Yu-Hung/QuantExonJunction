import pandas as pd
import argparse

# GTF 檔案預期的欄位
EXPECTED_COLUMNS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# option parser
def parse_option():
    parser = argparse.ArgumentParser(description="Find the longest transcript for each gene")
    parser.add_argument("-g", "--gtf", required=True, help="Input GTF file")
    parser.add_argument("-o", "--output", required=True, help="Output tsv file")
    parser.add_argument("-l", "--selected_gene", required=True, help="Selected gene list")
    parser.add_argument("-m", "--matrix", action="store_true", help="Output the selected Longest Transcripts-Gene_name matrix, tsv format")
    args = parser.parse_args()
    return args

# 讀取 GTF 檔案
def read_gtf(gtf_file:str) -> pd.DataFrame:
    gtf_df = pd.read_csv(gtf_file, sep='\t', header=None, comment='#')
    try:
        gtf_df.columns = EXPECTED_COLUMNS
    except ValueError:
        print(f"GTF file does not contain the expected columns : {EXPECTED_COLUMNS}")
        exit(1)

    # 解析 attribute 欄位
    gtf_df['Gene_name'] = gtf_df['attribute'].apply(lambda x: parse_attribute(x)['gene_name'])
    gtf_df['transcript_id'] = gtf_df['attribute'].apply(lambda x: parse_attribute(x)['transcript_id'])
    gtf_df['span'] = abs(gtf_df['end'] - gtf_df['start'])
    return gtf_df

# 取得最長的 transcript
def get_longest_transcript(gtf_df:pd.DataFrame) -> pd.DataFrame:
    # 只選出最長的 transcript id
    gtf_df = gtf_df[gtf_df['feature'] != 'gene']
    filtered_df = gtf_df.loc[gtf_df.groupby("Gene_name")["span"].idxmax()]
    filter_transcript_id_list = list(filtered_df['transcript_id'])
    gtf_df = gtf_df[gtf_df['transcript_id'].isin(filter_transcript_id_list)]
    return gtf_df

# 解析 GTF 的 attribute 欄位
def parse_attribute(attribute:str) -> dict:
    att = attribute.split('; ')
    att_dict = {}
    for item in att:
        key, value = item.split(' ')
        att_dict[key] = value.strip('"')
    return att_dict

# 讀取選出的基因列表
def read_selected_gene_list(selected_gene_list_file:str) -> list:
    with open(selected_gene_list_file, 'r') as file:
        selected_gene_list = file.read().splitlines()
    return selected_gene_list

# 驗證 transcript_id 是否為 ENST 開頭
def validate_transcript_id(selected_gtf:pd.DataFrame, transcript_id_list:list) -> bool:
    not_in_gtf = [transcript_id for transcript_id in transcript_id_list if transcript_id not in selected_gtf['transcript_id'].values]
    print(f"{len(transcript_id_list) - len(not_in_gtf)} / {len(transcript_id_list)} transcripts were found in selected GTF.")
    if len(not_in_gtf) > 0:
        print(f"Transcript ID: {not_in_gtf} not found in GTF file.")
    else:
        print("All transcripts were found in selected GTF.")

# 主程式
def main():

    # 讀取參數  
    args = parse_option()

    # 讀取 GTF 檔案
    gtf_df = read_gtf(args.gtf)

    # 取得最長的 transcript
    longest_transcript_df = get_longest_transcript(gtf_df)

    # 選出需要的欄位
    longest_transcript_df = longest_transcript_df[longest_transcript_df['feature'] == 'transcript']

    # 選出指定的基因
    selected_gene_list = read_selected_gene_list(args.selected_gene)
    selected_transcript_list = list(longest_transcript_df[longest_transcript_df['Gene_name'].isin(selected_gene_list)]['transcript_id'])
    
    # 選出指定的 transcript 以及排序
    selected_gtf = gtf_df[gtf_df['transcript_id'].isin(selected_transcript_list)]
    selected_gtf = selected_gtf.sort_values(by=['chr', 'start', 'end'])
    validate_transcript_id(selected_gtf, selected_transcript_list)

    # 輸出 matrix
    if args.matrix:
        export_column = ['chr', 'start', 'end', 'strand', 'Gene_name', 'transcript_id', 'span']
        selected_longest_transcript_df = selected_gtf[selected_gtf['feature'] == 'transcript'][export_column]
        selected_longest_transcript_df.to_csv("Selected_Longest_Transcripts.tsv", sep="\t", index=False, header=True)

    # 輸出結果
    selected_gtf = selected_gtf[EXPECTED_COLUMNS]
    selected_gtf.to_csv(args.output, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()

