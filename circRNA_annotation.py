import os,re
import subprocess
import glob

def process_merge_alg_files(directory='.'):
    """
    处理目录中的 _merge.txt 文件，提取指定列并转换为 bed 格式
    
    Args:
        directory (str): 要处理的目录路径，默认为当前目录
    """
    # 切换到指定目录
    original_dir = os.getcwd()
    os.chdir(directory)
    
    try:
        # 获取所有文件
        for file in glob.glob('*_merge.txt'):
            if os.path.isfile(file):
                output_file = f"{file.split('_merge.txt')[0]}_merge_alg.bed"
                
                # 使用 awk 处理文件
                with open(file, 'r') as infile, open(output_file, 'w') as outfile:
                    for line in infile:
                        # 使用冒号或制表符分割
                        fields = line.strip().split('\t')[0]
                        fields = re.split('[:|]', fields)
                        if len(fields) >= 3:
                            outfile.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\n")
    finally:
        # 恢复原始目录
        os.chdir(original_dir)

def convert_genepred_to_bed(input_file, output_file):
    """
    将 genePred 文件转换为 bed 格式
    
    Args:
        input_file (str): 输入的 genePred 文件路径
        output_file (str): 输出的 bed 文件路径
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                # 提取需要的字段并重新排列
                outfile.write(f"{fields[1]}\t{fields[3]}\t{fields[4]}\t{line.strip()}\n")

def perform_bedtools_intersect(directory='.', genepred_bed_file=''):
    """
    对所有 merge_alg.bed 文件执行 bedtools intersect 操作
    
    Args:
        directory (str): 包含 bed 文件的目录路径
        genepred_bed_file (str): genePred bed 文件的路径
    """
    # 切换到指定目录
    original_dir = os.getcwd()
    os.chdir(directory)
    
    try:
        # 获取所有 mergealg.bed 文件
        for bed_file in glob.glob('*_merge_alg.bed'):
            output_file = f"{bed_file.split('merge_alg.bed')[0]}merge_alg.bed.inter"
            
            # 执行 bedtools intersect 命令
            cmd = [
                'bedtools', 'intersect',
                '-a', bed_file,
                '-b', genepred_bed_file,
                '-wao'
            ]
            
            with open(output_file, 'w') as outfile:
                subprocess.run(cmd, stdout=outfile, check=True)
    finally:
        # 恢复原始目录
        os.chdir(original_dir)

import sys


def get_overlapping_exons(circrna_start, circrna_end, exon_starts, exon_ends):
    """
    获取circRNA坐标范围内重叠的外显子索引列表

    Args:
        circrna_start (int): circRNA起始坐标
        circrna_end (int): circRNA终止坐标
        exon_starts (list): 外显子起始坐标列表
        exon_ends (list): 外显子终止坐标列表

    Returns:
        list: circRNA范围内重叠的外显子索引列表
    """
    circrna_range = (circrna_start, circrna_end)
    overlapping_exons = []

    for i, (exon_start, exon_end) in enumerate(zip(exon_starts, exon_ends)):
        exon_range = (exon_start, exon_end)
        if max(circrna_range[0], exon_range[0]) <= min(circrna_range[1], exon_range[1]):
            overlapping_exons.append( (exon_start, exon_end))

    return overlapping_exons

def get_intron_ranges(exon_starts, exon_ends):
    """
    获取基因的内含子坐标范围列表

    Args:
        exon_starts (list): 外显子起始坐标列表
        exon_ends (list): 外显子终止坐标列表

    Returns:
        list: 内含子坐标范围列表,每个元素为一个元组(start, end)
    """
    intron_ranges = []
    for i in range(len(exon_starts) - 1):
        intron_start = exon_ends[i] + 1
        intron_end = exon_starts[i + 1] - 1
        intron_ranges.append((intron_start, intron_end))
    return intron_ranges

def anno_ext(file):
    dic = {}
    output_file = file + '.parser'
    with open(file) as f:
        for i in f.readlines():
            line = i.strip().split()
            if line[3] == '.':
                res = [line[0], line[1], line[2], '+', 'NA', 'NA', 'Intergenic', 0, '%s,%s'%(line[1], line[2])]
                dic.setdefault('%s:%s|%s'%(line[0], line[1], line[2]), []).append(res)
                continue
            chr, start, end, transcript, gene, exon_start, exon_end, strand = line[0], int(line[1]), int(line[2]), line[6], line[17], line[14], line[15], line[8]
            transcript_start, transcript_end, CDS_start, CDS_end = line[4], line[5], line[11], line[12]
            exon_starts = [int(x) for x in exon_start.strip(',').split(',')]
            exon_ends = [int(x) for x in exon_end.strip(',').split(',')]
            if  int(transcript_start) <= int(start) <  int(end) <= int(transcript_end):
                overlap_exons = get_overlapping_exons(start, end, exon_starts, exon_ends)
                if overlap_exons:
                    if overlap_exons[0][0] in range(start-2, start+2) and overlap_exons[-1][1] in range(end-2, end+2):
                        res = [chr, start, end, strand, gene, transcript, 'Exonic', len(overlap_exons),'|'.join(['%s,%s'%(x[0], x[1]) for x in overlap_exons])]
                        dic.setdefault('%s:%s|%s'%(chr, str(start), str(end)), []).append(res)
                else:
                    overlap_intron = get_overlapping_exons(int(start), int(end), [x[0] for x in get_intron_ranges(exon_starts, exon_ends)], [x[1] for x in get_intron_ranges(exon_starts, exon_ends)])
                    res = [chr, start, end, strand, gene, transcript, 'Intronic', 0, '%s,%s'%(start, end)]
                    dic.setdefault('%s:%s|%s'%(chr, str(start), str(end)), []).append(res)
    with open(output_file, 'w') as o: 
        for k in dic.keys():
            circ_anno = sorted(dic[k], key=lambda  x:x[7])
            circ_anno = circ_anno[-1]
            o.write('\t'.join([str(x) for x in circ_anno])+'\n')

def process_all_files(work_directory='.', genepred_file=None):
    """
    执行完整的处理流程
    
    Args:
        work_directory (str): 工作目录路径
        genepred_file (str): genePred 文件路径
    """
    if genepred_file is None:
        raise ValueError("必须提供 genePred 文件路径")
        
    # 1. 处理 merge_alg.txt 文件
    process_merge_alg_files(work_directory)
    
    # 2. 转换 genePred 文件为 bed 格式
    genepred_bed_output = os.path.join(work_directory, os.path.basename(genepred_file) + '.bed')
    convert_genepred_to_bed(genepred_file, genepred_bed_output)
    
    # 3. 执行 bedtools intersect
    perform_bedtools_intersect(work_directory, genepred_bed_output)

    # 4. 提取注释结果
    for bed_file in glob.glob('*merge_alg.bed.inter'):
        anno_ext(bed_file)

# 使用示例
if __name__ == '__main__':
    # 设置文件路径
    work_dir = '.'  # 或指定具体目录
    genepred_file = '/data3/indexs/index2/C_elegans/Caenorhabditis_elegans.WBcel235.111.genePred'
    
    # 执行完整流程
    process_all_files(work_dir, genepred_file)