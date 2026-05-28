"""
TLR7/TLR8 免疫激活Motif检测

功能：
    1. 检测siRNA乘客链和引导链中的已知TLR7/TLR8激活motif（如UGUGU、GUCCUUCAA等），并标记其位置和来源文献。
    2. 计算GU含量，评估通用免疫激活风险（GU含量>50%为高风险）。
    3. 根据链类型和位置，给出有位置约束的修饰建议
    4. 批量筛查候选序列，输出风险分级报告

关于修饰建议的重要说明：
    第 2 位：
        - AGO2 MID结构域与引导链5'端结合的关键位点
        - 2‘-OMe修饰显著降低AG2装载效率
        文献：Elkayam E et al., 2012, Cell, 150:100-110
    第 10~11 位：
        - AGO2催化切割靶mRNA的活性位点正对引导链的第10~11位
        - 该区域引入2'-OMe直接抑制mRNA切割活性，可能导致脱靶效应增加
        文献： Liu J et al., Science, 2004, 305:1437-1441
    第 2~8 位（seed区）：
        - seed区整体修饰降低RISC对靶mRNA的识别效率，可能显著降低siRNA的基因沉默活性
        - 因此在guide strand的seed区引入修饰需谨慎

乘客链对修饰的耐受性远高于引导链，可在大范围位置引入2'-OMe和2'-F修饰

基于以上限制，本脚本的修饰建议逻辑如下：
    · 乘客链：可大范围修饰，优先处理
    · 引导链：明确标注敏感位点，避开第2位和第10~11位
    · 全链统一修饰的建议不适用于引导链

文献依据：
    TLR7/TLR8 motif 鉴定：
        - Judge et al., 2005, Nat Biotechnol: "UGUGU" motif activates TLR7/8
        - Hornung et al., 2005, Nat Med: "GUCCUUCAA" motif activates TLR7
        - Forsbach et al., 2008, J Immunol: "UUUU"和"UUGU"等motif激活TLR7/8，且GU-rich序列更具刺激性
    GU-rich RNA的TLR识别：
        - Heil F et al., 2004, Science, 303:1526-1529
    2'-OMe修饰降低免疫激活：
        - Sioud M, 2006, Eur J Immunol, 36:1222-1230
    AGO2装载位点限制：
        - Elkayam E et al., 2012, Cell, 150:100-110
    AGO2催化位点限制：
        - Liu J et al., Science, 2004, 305:1437-1441
    化学修饰策略综述：
        - Varley AJ & Desalniers JP, RSC Adv, 2021, 11:2415-2426
    监管依据：
        - FDA Draft Guidance, Docket No. FDA-2024-D-4624 Nov 2024
"""

import pandas as pd
import re

# 已知高风险motif库（可持续扩展）
# 来源：
# 格式：{motif: (来源文献，激活机制)}
TLR_MOTIF_DB = {
    'UGUGU': ('Judge et al., 2005, Nat Biotechnol', 'TLR7/8激活, GU含量以来那个'),
    'GUCCUUCAA': ('Hornung et al., 2005, Nat Med', 'TLR7激活, GU含量非依赖性，激活浆细胞样树突状细胞产生IFN-α'),
    'UUUU': ('Forsbach et al., 2008, J Immunol', 'TLR8激活, AU-rich'),
    'UUGU': ('Forsbach et al., 2008, J Immunol', 'TLR7/8激活, GU-rich'),
    'GUUC': ('Forsbach et al., 2008, J Immunol', 'TLR7/8, GU-rich'),
}

# 引导链修饰敏感位点（1-indexed）
# 这些位置不建议引入2'-OMe修饰
GUIDE_SENSITIVE_POSITIONS = {
    2: 'AGO2 MID结构域结合位点',
    10: 'AGO2催化切割位点',
    11: 'AGO2催化切割位点'
}

# seed区位置（1-indexed）
SEED_REGION_POSITIONS = set(range(2, 8 + 1))  # 2-8位为seed区

def normalize_sequence(seq):
    """将DNA序列转换为RNA序列，并统一大写"""
    return seq.upper().replace('T', 'U')

def calc_gu_content(seq):
    """计算序列的GU含量(TLR7激活的通用风险指标)"""
    seq = normalize_sequence()
    if not seq:
        return 0.0
    gu_count = seq.count('G') + seq.count('U')
    return gu_count / len(seq)

def get_sensitive_pos_in_motif(motif_start:int,
                             motif_len:int,
                             strand_label:str):
    """
    检查motif覆盖的位置中是否包含引导链修饰敏感位点
    """
    if strand_label != 'guide':
        return []  # 乘客链不受这些位点限制

    motif_positions = set(range(motif_start, motif_start + motif_len))
    sensitive_hits = [
        pos for pos in GUIDE_SENSITIVE_POSITIONS
        if pos in motif_positions
    ]
    return sorted(sensitive_hits)

def get_seed_overlap(
    motif_start:int,
    motif_len:int,
    strand_label:str
):
    """
    检查motif是否与引导链seed区（2-8位）重叠
    seed区修饰会降低RISC对靶mRNA的识别效率
    """
    if strand_label != 'guide':
        return False

    motif_positions = set(range(motif_start, motif_start + motif_len))
    return bool(motif_positions & SEED_REGION_POSITIONS)

def build_modification_suggestion(
    strand_label:str,
    motif:str,
    motif_start:int,
    sensitive_positions:list,
    seed_overlap:bool
):
    """
    根据链类型、motif位置和敏感位点，生成有约束的修饰建议
    """
    if strand_label == 'passenger':
        # 乘客链：耐受性高，可直接在motif位置引入修饰
        return(
            f"乘客链第{motif_start}~{motif_start + len(motif) - 1}位含高风险motif"
            f"可在此区域引入2'-OMe修饰以降低TLR激活风险"
            f"乘客链对修饰耐受性高，优先处理"
        )

    # 引导链： 需要根据位置分情况处理
    if sensitive_positions:
        sensitive_desc = "、".join([
            f"第{pos}位（{GUIDE_SENSITIVE_POSITIONS[pos]}）"
            for pos in sensitive_positions
        ])
        return (
            f"引导链第{motif_start}~{motif_start + len(motif) - 1}位含高风险motif"
            f"但该区域覆盖修饰敏感位点：{sensitive_desc}。\n"
            f"      建议：\n"
            f"         · 优先在乘客链对应区域引入2'-OMe修饰\n"
            f"         · 引导链修饰需在活性验证后谨慎引入\n"
            f"         · 若沉默活性可接受，可尝试在敏感位点以外的motif位置引入2'-F修饰\n"
            f"          （2'-F对AGO2活性影响小于2'-OMe）\n"
            f"         ·根本解决方案：重新设计序列"
        )

    if seed_overlap:
        return(
            f"引导链第{motif_start}~{motif_start + len(motif) - 1}位含高风险motif"
            f"与seed区（第2~8位）有重叠。\n"
            f"      建议：\n"
            f"         · 优先在乘客链对应区域引入2'-OMe修饰\n"
            f"         · seed区修饰会降低RISC靶标识别效率，引导链此区域修饰需谨慎\n"
            f"         · 可考虑在motif中非seed区的位置引入2'-F修饰\n"
            f"         ·根本解决方案：重新设计序列"
        )

    # 引导链非敏感位点：可以修饰，但需说明约束
    return(
        f"引导链第{motif_start}~{motif_start + len(motif) - 1}位含高风险motif"
        f"该位置不在已知修饰敏感区域。\n"
        f"      建议：\n"
        f"         · 可在此位置引入引入2'-OMe修饰\n"
        f"         · 同时在乘客链对应区域引入修饰（双重保障）\n"
        f"         · 修饰后需验证沉默活性未显著下降\n"
    )

def build_gu_suggestion(
    strand: str,
    gu_content: float,
    seq_len: int
):
    """
    针对GU含量偏高的修饰建议
    区分乘客链和引导链的不同处理策略
    """
    if strand == "passenger":
        return(
            f"乘客链GU含量偏高（{gu_content:.1%}，>50%警戒线）。\n"
            f"      建议：\n"
            f"      · 乘客链可大范围引入2'-OMe/2'-F交替修饰\n"
            f"      · 乘客链对修饰耐受性高，全链修饰方案可行\n"
            f"      · 若条件允许，优先重新设计序列降低GU含量"
        )
    
    # 引导链GU含量高：不能全链修饰
    return(
        f"引导链GU含量偏高（{gu_content:.1%}，>50%警戒线）。\n"
        f"      注意：引导链不能全链统一引入2'-OMe修饰。\n"
        f"      建议策略（按优先级排序）：\n"
        f"      1. 重新设计序列降低GU含量（根本解决方案）\n"
        f"      2. 若序列不可替换：\n"
        f"          · 在引导链非敏感位点（避开第2位、第10~11位）\n"
        f"            选择性引入2'-OMe修饰\n"
        f"          · 在正义链大范围引入2'-OMe/2'-F修饰\n"
        f"          · 可考虑在引导链采用2'-OMe和2'-F交替修饰模式，\n"
        f"            但须明确跳过第2位和第10~11位\n"
        f"      3. 所有修饰方案均需实验验证沉默活性"
    )

def detect_immune_motifs(
    seq: str, 
    strand_label: str,
    motif_db: dict = TLR_MOTIF_DB,
    gu_threshold: float = 0.5):
    """
    检测单条RNA序列中的免疫激活风险
    
    参数：
    - seq: RNA序列
    - strand_label: "guide"或"passenger"，用于后续修饰建议
    - motif_db: motif数据库
    - gu_threshold: GU含量阈值（可选，默认0.5——文献参考：>50%为高风险）
    
    返回：
        每个风险发现对应一个字典，包含：
        - motif、位置、机制、文献来源
        - 是否覆盖敏感位点（引导链专属）
        - 是否与seed区重叠（引导链专属）
        - 有位置约束的修饰建议
    """
    seq = normalize_sequence(seq)
    detected_motifs = []
    
    # 检测已知高风险motif
    for motif, (source, mechanism) in motif_db.items():
        positions = [m.start() for m in re.finditer(motif, seq)]
        if positions:
            for pos in positions:
                motif_start_1_indexed = pos + 1 #转为 1-indexed
                sensitive_pos = get_sensitive_pos_in_motif(
                    motif_start_1_indexed, len(motif), strand_label
                )
                seed_overlap = get_seed_overlap(
                    motif_start_1_indexed, len(motif), strand_label
                )
                suggestion = build_modification_suggestion(
                    strand_label, motif, motif_start_1_indexed,
                    sensitive_pos, seed_overlap
                )
                detected_motifs.append({
                    "risk_type":        "已知TLR motif",
                    "strand": strand_label,
                    "motif": motif,
                    "position_1indexed": motif_start_1_indexed,
                    "mechanism": mechanism,
                    "literature_source": source,
                    "sensitive_pos_overlap": sensitive_pos,
                    "seed_region_overlap": seed_overlap,
                    "modification_suggestion": suggestion
                })
    
    # GU含量通用风险评估
    gu_content = calc_gu_content(seq)
    if gu_content > gu_threshold:
        detected_motifs.append({
            'risk_type': "GU含量偏高"，
            'strand': strand_label,
            'motif': f"GU-rich (GU含量={gu_content:.1%})",
            'source': '文献综合分析',
            'mechanism': 'TLR7/8 general GU-recongnition',
            'position': '全序列',
            'is_seed_region': False,
            'modification_suggestion': build_gu_suggestion（
                strand_label, gu_content, len(seq)
            ）
        })
    
    return detected_motifs