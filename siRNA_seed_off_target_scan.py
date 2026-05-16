"""
siRNA seed区全转录组3'UTR脱靶扫描

主要功能：
1、提取siRNA引导链seed区（第2-8位），并生成反向互补序列
2、扫描人类全转录组3'UTR，寻找精确匹配
3、对命中位置计算杂交自由能，过滤假阳性
4、对高置信命中做功能注释，标记高风险基因

依赖：
   pip install biopython ViennaRNA pandas

数据准备：
    人类3'UTR序列（Ensembl BioMart下载）
    链接：https://www.ensembl.org/biomart/martview/a5362c91587b87582ea6b3653f30104d
    Dataset: Ensembl Genes 115 → Human genes (GRCh38.p14) 
    Attributes: Sequence → 3' UTR sequences
                HEADER INFORMATION → Transcript stable ID, Gene stable ID, Gene name
    保存为：human_3UTR.fa
"""
import RNA
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from pathlib import Path

# ======================================
# Part 1:seed区提取
# ======================================
def get_seed_motif(guide_strand):
    """
    提取siRNA引导链seed区（第2-8位），并生成反向互补序列
    """
    if len(guide_strand) < 8:
        raise ValueError(f"引导链过短（{len(guide_strand)}nt）至少需要8nt")
    
    # 提取seed区
    guide_strand = guide_strand.upper().replace('T', 'U')  # 转换为RNA序列
    seed = guide_strand[1:8]  # 第2-8位
    # 生成反向互补序列
    seed_rc = str(Seq(seed).reverse_complement_rna())
    
    return seed, seed_rc

# ======================================
# Part 2:全转录组3'UTR扫描
# ======================================
def scan_3utr_for_seed(seed_rc, utr_fasta, context_window=10):
    """
    扫描人类全转录组3'UTR，寻找精确匹配
    返回命中列表：[(gene_id, gene_name, utr_seq, match_pos), ...]
    """
    hits = []
    seed_rc_clean = seed_rc.upper().replace('T', 'U')  # 确保seed反向互补序列是RNA格式

    for record in SeqIO.parse(utr_fasta, "fasta"):
        seq = str(record.seq).upper().replace('T', 'U')  # 转换为RNA序列
        
        # 解析Fasta header，假设格式为：>transcript_id gene_id gene_name 
        parts = record.id.split('|')
        transcript_id = parts[0] if len(parts) > 0 else record.id
        gene_id = parts[1] if len(parts) > 1 else "Unknown"
        gene_name = parts[2] if len(parts) > 2 else "Unknown"
        
        # 滑窗搜索精确匹配
        pos = 0
        while True:
            pos = seq.find(seed_rc_clean, pos)
            if pos == -1:
                break
            # 提取命中位置的上下文序列（用于热力学计算）
            # 取命中位置前后各context_window个碱基，确保不越界
            ctx_start = max(0, pos - context_window)
            ctx_end = min(len(seq), pos + len(seed_rc_clean) + context_window) 
            utr_fragment = seq[ctx_start:ctx_end]
            hits.append({
                "transcript_id": transcript_id,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "hit_position": pos + 1,  # 转换为1-based位置
                "utr_seq": utr_fragment,
                "full_utr_seq": seq
            })
            pos += 1  # 继续搜索下一个位置
        
    print(f"精确匹配完成：共命中 {len(hits)} 个位置（来自全转录组3'UTR扫描）)")
    print(f"进入热力学过滤步骤...\n")
    
    return hits

# ======================================
# Part 3:热力学过滤
# ======================================
def calc_seed_binding_energy_vienna(seed_rc, utr_fragment):
    """
    使用ViennaRNA的RNA.cofold()计算seed区与3'UTR片段的杂交自由能
    返回ΔG值（kcal/mol）——越负表示结合越稳定
    参数：
    seed_rc: seed区反向互补序列（7nt）
    utr_fragment: 3'UTR片段序列（上下文）
    """
    # 确保序列格式统一
    seed_rc_clean = seed_rc.upper().replace('T', 'U')
    utr_fragment_clean = utr_fragment.upper().replace('T', 'U')
    # RNA.cofold()需要输入格式为 "seq1&seq2"
    cofold_input = f"{seed_rc_clean}&{utr_fragment_clean}"
    try:
        (structure, energy) = RNA.cofold(cofold_input)
        return energy
    except Exception as e:
        # 
        print(f"警告：计算热力学能量时出错：{e}，跳过次位点")
        return 0.0  # 返回一个默认值，表示无法计算能量，后续会被过滤掉

def calc_seed_binding_energy_rnahybrid(seed_seq, utr_fragment, rnahybrid_path="RNAhybrid"):
    """
    使用RNAhybrid计算seed区与3'UTR片段的杂交自由能
    返回ΔG值（kcal/mol）——越负表示结合越稳定
    参数：
    seed_seq: seed区序列（7nt，非反向互补）
    utr_fragment: 3'UTR片段序列（上下文）
    rnahybrid_path: RNAhybrid可执行文件路径
    """
    import tempfile, os
    
    # 创建临时文件保存输入序列
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_input:
        tmp_input.write(f">seed\n{seed_seq}\n")
        tmp_query  = tmp_input.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_target:
        tmp_target.write(f">utr\n{utr_fragment}\n")
        tmp_target = tmp_target.name  

    mfe = 0.0
    try:
        cmd = [
            rnahybrid_path,
            '-q', tmp_query,
            '-t', tmp_target,
            '-f', '2,8'        # 只考虑seed区（第2-8位）的结合
            '-c'               # 紧凑输出
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split(":")
                if len(parts) >= 5:
                    mfe = float(parts[4])  # MFE值通常在第5列
                    break
    except Exception as e:
        print(f"警告：调用RNAhybrid时出错：{e}，跳过次位点")
    finally:
        # 清理临时文件
        os.remove(tmp_query)
        os.remove(tmp_target)
    
    return mfe

def filter_hits_by_energy(hits, seed_rc, energy_threshold=-8.0, use_rnahybrid=False, seed_seq=""):
    """
    对精确匹配的命中位点进行热力学过滤，保留ΔG值小于energy_threshold的高置信命中
    参数：
    hits: 精确匹配的命中列表
    seed_rc: seed区反向互补序列
    energy_threshold: 能量阈值
    use_rnahybrid: 是否使用RNAhybrid进行计算
    seed_seq: seed区序列（非反向互补）
    """
    print(f"正在对 {len(hits)} 个命中位点进行热力学过滤，能量阈值：{energy_threshold} kcal/mol...")

    results = []
    for i, hit in enumerate(hits):
        if i % 100 == 0:
            print(f"正在处理第 {i+1}/{len(hits)} 个命中位点...")
        
        if use_rnahybrid and seed_seq:
            energy = calc_seed_binding_energy_rnahybrid(seed_seq, hit['utr_seq'])
        else:
            energy = calc_seed_binding_energy_vienna(seed_rc, hit['utr_seq'])

        # 过滤掉能量不满足条件的命中
        if energy < energy_threshold:
            results.append({
                **hit,
                "binding_energy": energy,
                "energy_threshold": energy_threshold
            })
    
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(by="binding_energy")  # 按结合能排序，越负越靠前
    
    print(f"热力学过滤完成：共保留 {len(df)} 个高置信命中位点（ΔG < {energy_threshold} kcal/mol）")
    return df

# ======================================
# Part 4:功能注释与风险评分
# ======================================
# 内置高风险基因列表（可根据项目扩展）
# 来源：心脏安全性（ICH S7A/S7B）、肝毒性（liverTox）、肿瘤抑制基因（COSMIC）、重要转录因子（TFDB）、凋亡通路
HIGH_RISK_GENES_SETS = {
    "cardiotoxicity": {"KCNH2", "SCN5A", "CACNA1C", "KCNE1", "KCNE2"},  # 心脏安全性相关基因
    "hepatotoxicity": {"CYP3A4", "CYP2D6", "UGT1A1", "SULT1A1", "GSTP1", "PCSK9", "ALB", "APOB"},  # 肝毒性相关基因
    "tumor_suppressors": {"TP53", "RB1", "PTEN", "BRCA1", "BRCA2"},  # 肿瘤抑制基因
    "transcription_factors": {"MYC", "NF-kB1", "SP1", "FOXO3", "GATA3"},  # 重要转录因子
    "apoptosis_pathway": {"TP53","BCL2", "BCL2L1", "BAX", "CASP3", "CASP8", "CASP9", "FAS", "APAF1"},  # 凋亡通路相关基因
    "neurotoxicity": {"APP", "PSEN1", "PSEN2", "MAPT", "SNCA", "HTT", "TARDBP", "SOD1"},  # 神经毒性相关基因（阿尔茨海默病、帕金森病相关）
}

def annotate_risk(df, risk_genes_sets=HIGH_RISK_GENES_SETS):
    """
    对高置信命中位点进行功能注释，标记是否命中高风险基因
    标记规则：
        高风险：命中基因属于预设高风险基因集
        中风险：energy < -10 kcal/mol但不属于高风险基因集
        低风险：energy >= -10 kcal/mol且不属于高风险基因集
    参数：
        df: 包含高置信命中位点的DataFrame，来自filter_hits_by_energy()函数的输出
        risk_genes_sets: 高风险基因集合字典，键为风险类型，值为基因名集合
    返回：
        包含风险注释的DataFrame（增加risk_category 和 risk_gene_sets列）
    """
    if df.empty:
        print("没有高置信命中位点需要注释")
        return df
    
    # 建立基因名到风险类型的映射
    gene_to_risk = {}
    for category, genes in risk_genes_sets.items():
        for gene in genes:
            gene_to_risk[gene] = category

    def classify_row(row):
        gene = row.get('gene_name', 'Unknown')
        bind_energy = row.get('binding_energy', 0.0)
        if gene in gene_to_risk:
            return "高风险", gene_to_risk[gene] 
        elif bind_energy < -10.0:
            return "中风险（强结合）", "-"
        else:
            return "低风险", "-"
    
    df[['risk_category', 'risk_gene_sets']] = df.apply(
        lambda row: pd.Series(classify_row(row)), axis=1
    )

    return df

def generate_report(df, guide_strand, seed, seed_rc, energy_threshold):
    """
    生成评估报告
    """
    high_risk_hits = df[df['risk_category'] == "高风险"]
    medium_risk_hits = df[df['risk_category'] == "中风险（强结合）"]

    if len(high_risk_hits) > 0:
        conclusion = "需关注——命中高风险功能基因，建议重新设计siRNA序列"
    elif len(medium_risk_hits) > 5:
        conclusion = "中等风险——多个强结合位点，建议优化siRNA设计"
    else:
        conclusion = "通过"

    report = f"""
================================
siRNA Seed区脱靶评估报告
================================
引导链序列: {guide_strand}
Seed区（2-8位）: {seed}
Seed区反向互补序列: {seed_rc}
能量过滤阈值: {energy_threshold} kcal/mol
--------------------------------

---结果摘要---
总命中位点数: {len(df)}
高风险位点数: {len(high_risk_hits)}
中风险位点数: {len(medium_risk_hits)}
低风险位点数: {len(df) - len(high_risk_hits) - len(medium_risk_hits)}
最强结合能量: {df['binding_energy'].min():.2f if not df.empty else 'N/A'} kcal/mol (df.iloc[0]['gene_name'] if not df.empty else 'N/A')
--------------------------------

---高风险位点详情---
"""
    if high_risk_hits.empty:
        report += "  无\n"
    else:
        for _, row in high_risk_hits.iterrows():
            report += f"  - 基因: {row['gene_name']} (ID: {row['gene_id']}), 位置: {row['hit_position']}, ΔG: {row['binding_energy']:.2f} kcal/mol, 风险类型: {row['risk_gene_sets']}\n"
    
    report += f"""
---结论---
{conclusion}

注：本报告仅覆盖seed区介导的脱靶风险，完整评估还需包括：
1. 非seed区脱靶（全序列匹配、非完全匹配）
2. 免疫激活风险（如TLR7/8识别）
3. 乘客链脱靶风险
"""
    return report

# ======================================
# 主流程函数
# ======================================
if __name__ == "__main__":
    # -- 用户输入 --
    # Inclisiran引导链序列（示例）
    guide_strand = "UUCAAGCCAUAUGAAUUCA"
    utr_fasta = "human_3UTR.fa"  # 全转录组3'UTR序列文件
    energy_threshold = -8.0  # 热力学过滤阈值（kcal/mol）
    OUTPUT_DIR = "siRNA_off_target_report"
    Path(OUTPUT_DIR).mkdir(exist_ok=True)

    # -- 评估流程 --
    print("=== siRNA Seed区脱靶评估开始 ===\n")
    # 1. 提取seed区及反向互补序列
    seed, seed_rc = get_seed_motif(guide_strand)
    print(f"引导链: {guide_strand}")
    print(f"Seed区（2-8位）: {seed}")
    print(f"Seed区反向互补序列: {seed_rc}\n")
    # 2. 全转录组3'UTR扫描
    hits = scan_3utr_for_seed(seed_rc, utr_fasta, context_window=10)
    if not hits:
        print("未发现任何精确匹配位点，脱靶风险较低")
        exit()
    # 3. 热力学过滤
    filtered_df = filter_hits_by_energy(hits, seed_rc, energy_threshold=energy_threshold, 
                                        use_rnahybrid=False, seed_seq=seed)
    if filtered_df.empty:
        print(f"热力学过滤后无高置信命中位点，脱靶风险较低（阈值: {energy_threshold} kcal/mol）")
        print("Seed区脱靶风险评级：低")
        exit()
    # 4. 功能注释与风险评分
    annotated_df = annotate_risk(filtered_df)
    # 5. 生成评估报告
    report = generate_report(annotated_df, guide_strand, seed, seed_rc, energy_threshold)
    print(report)
    # 保存详细结果到CSV文件
    annotated_df.to_csv(f"{OUTPUT_DIR}/siRNA_seed_off_target_hits.csv", index=False)
    print(f"详细命中位点信息已保存到: {OUTPUT_DIR}/siRNA_seed_off_target_hits.csv")

        