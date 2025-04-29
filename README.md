# 定量 exon-splicing events

這個小專案用於定量我們 RNA panel 檢測基因的 exon-splicing 定量, 我們把範圍限縮在 40 個基因, 並且只算這些基因最長的 transcript。使用 Featurecount 算 STAR-align 完之後的 BAM 檔, 他會輸出 sample alignment 中所有splicing event的定量表格, 從中我們在使用bedtool intersect交集找出目標基因的 exon, 去計算他們的 splicing events。

* 2025-04-29 已更新 Target-Exon-Junction BED, 新增最新檔案的被份, 並在 pipeline 中使用, 避免每次更新都要重新建置
* 2025-03-04 加入 Internal Control 基因製作新的 Target-Exon-Junction BED
* 2025-01-08 新增一欄 Kinase Gene 用於近一步篩選
* 如果已經準備好 Target-Exon-Junction BED 且沒有更新, 不用在做一次。
* 目前的 Target_Exon_Junction_BED 為 `prepare_exon_annotation/RNA-Panel_GeneList_ExJ_2025_Feb.bed`

## 安裝軟體

1. 下載最新版本 subread (包含 featurecount)：[連結](https://sourceforge.net/projects/subread/files/subread-2.0.8/subread-2.0.8-Linux-x86_64.tar.gz/download)
2. 解壓縮 subread-2.0.8-Linux-x86_64.tar.gz
3. 執行檔在 subread-2.0.8-Linux-x86_64/bin/ 把他加入環境變數 PATH 就能使用了

## 使用方法

### ＊ Step1: 計算 exon splicing

```bash
bash featureCount.sh <BAM-file> <Genome-GTF> <output-prefix>
```

* 由於是用於 WT exon-junction 的定量作為 control QC, 不算跨越 chromosome 的 splicing
* 不管用哪一種 GTF 都會輸出 BAM 檔中全部的 exon-splicing, 差別是有結果沒有被 Annotate
* 額外輸出 BED 檔格式, 用於後續的交集分析 (output.junction.bed)

### ＊ Step2: 註記 Exon 跨度

```bash
python count_exon_span.py -b <junction.bed> -g <ref_annot.gtf>
```

* <junction.bed> 輸入 Step 的 output BED 檔
* <ref_annot.gtf> 輸入完整 reference genome 的 GTF 檔

為了後續處理我們想看的 splicing 狀況, 在 featurecount 產生 exon-splicong 定量結果的 BED 檔之後, 我們還要幫她找出來的位置註記 exon 跨度, 這 script 參考 hg19 genome GTF 去計算每一個位置的 exon 跨度, 正常來說我們想看得是 exon 跨度為 1 的 splicing event。最後 data clean up 時會一並篩選。

### ＊ Step3: 註記 featurecount 結果

```bash
# --keep 會保留中繼檔, 中繼檔還沒刪掉不符合的 exon-splicing.
python AnnotateResult.py -b <junction_ExonSpan.bed> -a <Target_Exon_Junction_.bed> -o Annotated_featurecount.xlsx --keep
```

### ＊ ！！前準備 Target-Exon-Junction BED！！

1. 複製這個 Repository, 進入到 `prepare_exon_annotation/` 資料夾

   ```
   git clone git@github.com:Hsieh-Yu-Hung/QuantExonJunction.git
   cd prepare_exon_annotation
   ```
2. 準備 hg19 GTF 檔案 (不包含在這個 Repository)
3. 使用  `Find_Longest_Transcript.py` 篩選 GTF

```bash
python Find_Longest_Transcript.py \			# 挑出最長 transcript, 目標基因的 GTF
 -g ref_annot.gtf \ 					# hg19 GTF
 -l GeneList.txt \					# RNA panel 檢測目標基因
 -o output_select.gtf \					# 輸出篩選後的 GTF
 -m							# 輸出最長的 transcript_id - gene_name 對照表
```

3. 使用  `Make_Exon_Junction_BED.py` 製作出定量目標範圍

```bash
python Make_Exon_Junction_BED.py \			  # 製作 Exon-junction BED 檔
 -g test_select.gtf \					  # 請提供篩選後的ＧＴＦ
 -o test_selected_ej.bed				  # 輸出 Exon-junction BED
```

* 將目標基因最長 transcript 的 exons 製作成 BED 檔, 用於後續的結果交集和定量呈現
* 將會列出目標基因所有 exon-junction 的定量結果, **手動標示**出哪一些是我們設計 primer 的地方

4. 手動加入 `ACCUiNPanel` 和 `KinaseGene` 兩個欄位
5. 製作完之後, 上傳至 GCS, Container 在運行時會自動下載, 路徑：
```
gs://accuinbio-core-dev/VEP_hg19_Database/ACCUiN_Panel
```
6. [重要‼️] 更新 `Cloud_JOB_STAR_FUSION_Pipeline` 中的 `main_scripts/step2-STAR_align.sh
```bash
# Annotate featurecount result
python task_scripts/AnnotateResult.py \
 -b /app/$outdir/Exon_Quant.junction_ExonSpan.bed \
 -a RNA-Panel_GeneList_ExJ_2025_Feb.bed \    # 換成新製作的 Target Exon BED
 -o /app/$outdir/Annotated_featurecount.xlsx
```
7. [重要‼️] 重新建置`Cloud_JOB_STAR_FUSION_Pipeline`並上傳到存放區

8. [重要⚠️] 檢查製作好的檔案, 若 PCR 放大目標沒有包含在裡面, 則手動加入❗️
- 因為 GTF 的關係, 有一些基因的 exon junction 和設計有所不同, 以實際上設計的為準❗️

9. 將最新版檔案複製一份為固定名稱, 方便未來更新, 不用再重新建置 STAR Fusion 所有 Images
- 名稱為 `RNA-Panel_GeneList_ExJ_Implemented.bed` 上傳至 GCS.

## 執行範例

* 使用 example/ 底下的檔案
* ref_anno.gtf 需要自行下載, 這是 hg19 的 GTF 檔案

- Step1:

```bash
bash featureCount.sh example/input_files/MGT01-STAR_Aligned.out.bam ref_anno.gtf MGT01_feature_count
# 輸出為 example/output_files/MGT01_feature_count*
```

- Step2:

```bash
python count_exon_span.py -b example/output_files/MGT01_feature_count.junction.bed -g ref_anno.gtf
# 輸出為 example/output_files/MGT01_feature_count.junction_ExonSpan.bed
```

- Step3:

```bash
python AnnotateResult.py \
 -b example/output_files/MGT01_feature_count.junction_ExonSpan.bed \
 -a example/input_files/RNA-Panel_GeneList_ExJ_2024_Nov.bed \
 -o MGT01_Annotated_Feature_count.xlsx --keep
# 輸出為 example/output_files/MGT01_Annotated_Feature_count.xlsx
# 以及 example/output_files/annotated_result.tsv
```

## 額外說明

### **1. featurecount 參數說明**

```bash
featureCounts \
 -T 8 \ 						# 8 CPU
 -a ../ctat_genome_lib_build_dir/ref_annot.gtf \	# GTF
 -o FeatureCount_Result_WholeGenome.txt \		# Output
 -t exon -g transcript_id \				# -t 計算的level -g 計算目標
 -M -J -p \						# -M 允許計算 multimap -J 計算exon-splicing -p Pair-end data
 -C --countReadPairs \ 					# -C 不算跨越chromosome --countReadPairs 算一整個 fragment
 -B STAR-align.bam					# STAR 輸出的 BAM 檔
```

* 輸出範例：

| 檔案                | 說明                            | 輸出         |
| ------------------- | ------------------------------- | ------------ |
| output.txt          | exon 定量結果                   | featurecount |
| output.txt.jcount   | exon-junction 定量結果         | featurecount |
| output.txt.summary  | 統計定量結果                    | featurecount |
| output.junction.bed | exon-junction BED檔 + 定量結果 | 額外         |

```plaintext
# output.txt.jcounts
PrimaryGene	SecondaryGenes	Site1_chr	Site1_location	Site1_strand	Site2_chr	Site2_location	Site2_strand	S1	S2	S3	S4	S5	S6
ENST00000302692.6	NA	chr1	9599806	NA	chr1	9613684	NA	0	0	0	0	2	1
ENST00000253251.8	ENST00000377157.3,ENST00000343090.6,ENST00000377153.1	chr1	10093752	NA	chr1	10132086	NA	0	0	0	0	2	0
ENST00000263934.6	ENST00000377081.1,ENST00000377083.1,ENST00000377086.1,ENST00000377093.4,ENST00000497835.1	chr1	10342591	NA	chr1	10351140	NA	0	0	0	0	1	0
ENST00000263934.6	ENST00000377081.1,ENST00000377083.1,ENST00000377086.1,ENST00000377093.4,ENST00000497835.1	chr1	10351219	NA	chr1	10352105	NA	0	0	0	0	1	0
ENST00000376957.2	ENST00000487300.1,ENST00000459997.1	chr1	11119402	NA	chr1	11119834	NA	0	0	0	0	0	5
ENST00000361445.4	ENST00000495435.1,ENST00000476768.1,ENST00000445982.1,ENST00000420480.1	chr1	11199715	NA	chr1	11204705	NA	0	0	0	0	1	1
ENST00000376369.3	ENST00000491536.1,ENST00000196061.4	chr1	12027148	NA	chr1	12030727	NA	0	0	0	0	1	0
ENST00000376369.3	ENST00000491536.1,ENST00000196061.4,ENST00000481933.1	chr1	12030873	NA	chr1	12032929	NA	0	0	0	0	1	0
ENST00000444836.1	ENST00000235329.5	chr1	12058935	NA	chr1	12059045	NA	0	0	0	0	1	16
ENST00000375799.3	ENST00000375793.2,ENST00000462455.1	chr1	16044487	NA	chr1	16045033	NA	0	0	0	0	1	7
ENST00000375799.3	ENST00000375793.2,ENST00000462455.1	chr1	16045120	NA	chr1	16046229	NA	0	0	0	0	0	2
ENST00000430580.2	ENST00000392963.5	chr1	16918808	NA	chr1	16919936	NA	0	0	0	0	0	1
```

* 欄位說明：

| **欄位名稱** | **說明**                     |
| ------------------ | ---------------------------------- |
| PrimaryGene        | splice-donor gene                  |
| SecondaryGenes     | splice-receptor genes （可能多個） |
| Site1_chr          | dornor 位置 chr                    |
| Site1_location     | dornor 位置 position               |
| Site1_strand       | dornor 基因方向                    |
| Site2_chr          | receptor 位置 chr                  |
| Site2_location     | receptor 位置 position             |
| Site2_strand       | receptor 基因方向                  |
| S1-S6              | 不同 BAM 檔的定量                  |

### **2.  Target-Exon-Junction BED 範例**

* 輸出範例：

```plaintext
chr1	154130197	154142876	ENST00000271850.7	TPM3	exon9-8	Yes
chr1	154142945	154143125	ENST00000271850.7	TPM3	exon8-7	No
chr1	154143187	154144505	ENST00000271850.7	TPM3	exon7-6	No
chr1	154144580	154145384	ENST00000271850.7	TPM3	exon6-5	No
chr1	154145454	154145560	ENST00000271850.7	TPM3	exon5-4	No
chr1	154145677	154148591	ENST00000271850.7	TPM3	exon4-3	Yes
chr1	154148724	154163662	ENST00000271850.7	TPM3	exon3-2	No
chr1	154163787	154164378	ENST00000271850.7	TPM3	exon2-1	No
```

* 欄位說明：

| **欄位名稱** | **說明**                       |
| ------------------ | ------------------------------------ |
| chr                | 染色體編號（BED 格式）               |
| start              | 起始位置（BED 格式）                 |
| end                | 結束位置（BED 格式）                 |
| transcript_id      | 該基因最長的轉錄本 ID                |
| gene_name          | 基因名稱                             |
| exon-junction      | exon-junction 註記                   |
| primer_design      | 手動標註是否為我們設計 primer 的位點 |

### **3.  最終 Annotated featurecount table 範例**

* 輸出範例 (Excel):

```plaintext
chr	start	end	transcript_id	gene_name	exon_junction	Panel_Target	fc-chr	fc-start	fc-end	fc-transcript_id	fc-count	exon_span
chr2	29416788	29419636	ENST00000389048.3	ALK	exon29-28	No	NA	NA	NA	NA	0	NA
chr2	29419726	29420408	ENST00000389048.3	ALK	exon28-27	No	NA	NA	NA	NA	0	NA
chr2	29420542	29430037	ENST00000389048.3	ALK	exon27-26	No	NA	NA	NA	NA	0	NA
chr2	29430138	29432652	ENST00000389048.3	ALK	exon26-25	No	NA	NA	NA	NA	0	NA
chr2	29432744	29436850	ENST00000389048.3	ALK	exon25-24	No	NA	NA	NA	NA	0	NA
chr2	29436947	29443572	ENST00000389048.3	ALK	exon24-23	No	NA	NA	NA	NA	0	NA
chr2	29443701	29445210	ENST00000389048.3	ALK	exon23-22	No	NA	NA	NA	NA	0	NA
chr2	29445274	29445383	ENST00000389048.3	ALK	exon22-21	No	NA	NA	NA	NA	0	NA
chr2	29445473	29446208	ENST00000389048.3	ALK	exon21-20	Yes	NA	NA	NA	NA	0	NA
chr2	29446394	29448327	ENST00000389048.3	ALK	exon20-19	Yes	chr2	29446312	29498275	ENST00000389048.3	21	1
chr2	29448431	29449788	ENST00000389048.3	ALK	exon19-18	Yes	chr2	29446312	29498275	ENST00000389048.3	11	1
chr2	29449940	29450440	ENST00000389048.3	ALK	exon18-17	Yes	chr2	29446312	29498275	ENST00000389048.3	171	1
chr2	29450538	29451750	ENST00000389048.3	ALK	exon17-16	Yes	chr2	29446312	29498275	ENST00000389048.3	1	1
chr2	29451932	29455170	ENST00000389048.3	ALK	exon16-15	No	chr2	29446312	29498275	ENST00000389048.3	1	1
chr2	29455314	29456431	ENST00000389048.3	ALK	exon15-14	No	chr2	29446312	29498275	ENST00000389048.3	1	1
chr2	29456562	29462546	ENST00000389048.3	ALK	exon14-13	No	chr2	29446312	29498275	ENST00000389048.3	1	1
chr2	29462696	29473971	ENST00000389048.3	ALK	exon13-12	No	chr2	29446312	29498275	ENST00000389048.3	1	1
chr2	29474133	29497965	ENST00000389048.3	ALK	exon12-11	No	chr2	29446312	29498275	ENST00000389048.3	1	1
chr2	29498093	29498268	ENST00000389048.3	ALK	exon11-10	Yes	chr2	29446312	29498275	ENST00000389048.3	391	1
chr2	29498362	29519754	ENST00000389048.3	ALK	exon10-9	Yes	chr2	29498362	29519754	ENST00000389048.3	67	1
chr2	29519923	29541170	ENST00000389048.3	ALK	exon9-8	No	NA	NA	NA	NA	0	NA
chr2	29541270	29543617	ENST00000389048.3	ALK	exon8-7	No	NA	NA	NA	NA	0	NA
chr2	29543748	29551216	ENST00000389048.3	ALK	exon7-6	No	NA	NA	NA	NA	0	NA
chr2	29551347	29606598	ENST00000389048.3	ALK	exon6-5	No	NA	NA	NA	NA	0	NA
chr2	29606725	29754781	ENST00000389048.3	ALK	exon5-4	No	NA	NA	NA	NA	0	NA
chr2	29754982	29917716	ENST00000389048.3	ALK	exon4-3	No	NA	NA	NA	NA	0	NA
chr2	29917880	29940444	ENST00000389048.3	ALK	exon3-2	No	NA	NA	NA	NA	0	NA
chr2	29940563	30142859	ENST00000389048.3	ALK	exon2-1	Yes	chr2	29940563	30142859	ENST00000389048.3	20	1
```
