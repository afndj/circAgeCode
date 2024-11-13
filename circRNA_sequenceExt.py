import pybedtools,sys
from Bio.Seq import Seq

pybedtools.helpers.set_tempdir('/data2/tmpdata/')
def extract_circrna_seq(chrom, coords, strand, ref_genome):
    """
    从参考基因组中提取circRNA序列并拼接

    Args:
        chrom (str): 染色体名称
        coords (str): 用"|"分隔的坐标对字符串
        strand (str): 链方向,'+' 或 '-'
        ref_genome (str): 参考基因组FASTA文件路径

    Returns:
        str: 拼接后的circRNA序列
    """
    circrna_seq = ""
    for coord_pair in coords.split("|"):
        start, end = coord_pair.split(",")
        a = pybedtools.BedTool('%s\t%s\t%s\t.\t.\t%s'%(chrom, start, end, strand), from_string=True)
        a = a.sequence(fi=ref_genome)
        circrna_seq += open(a.seqfn).readlines()[1].strip()
    if strand == '-':
        circrna_seq = str(Seq(circrna_seq).reverse_complement())
    return circrna_seq

def To_circSeq_fa(file, ref_genome):
    with open(file, "r") as f, open(file + '.fa', 'w') as o:
        for line in f:
            fields = line.strip().split("\t")
            chrom, start, end, strand, gene, transcript, region, num_exons, coords = fields
            circrna_seq = extract_circrna_seq(chrom, coords, strand, ref_genome)

            o.write(f">{chrom}:{start}-{end}|{gene}|{transcript}|{region}|{strand}\n")
            o.write(circrna_seq + '\n')

if __name__ == '__main__':
    # 输入文件路径
    input_file = sys.argv[1]
    # 参考基因组FASTA文件路径
    ref_genome = sys.argv[2]
    To_circSeq_fa(input_file, ref_genome)