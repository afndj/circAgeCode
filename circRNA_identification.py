import sys, os
import argparse
import subprocess

def run_cmd(x):
    sys.stderr.write('running %s\n' % x)
    out = subprocess.run(
        x, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if out.returncode != 0:
        sys.stderr.write('cmd running error\n')
        sys.stderr.write('args: %s\n' % out.args)
        sys.stderr.write('log: %s\n' % out.stdout)
        sys.stderr.write('log: %s\n' % out.stderr)
    else:
        sys.stderr.write('done\n')

def align(index_path, output, fq1, fq2=None):
    cmd = 'STAR --chimSegmentMin 10 --runThreadN 30 --genomeDir %s \
        --readFilesIn %s --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
        --chimOutType Junctions SeparateSAMold --outFileNamePrefix %s'%(index_path, ' '.join([fq1, fq2]) if fq2 else fq1, output)
    run_cmd(cmd)

def circExplorer2(Chimeric, ref, fasta, out):
    cmd1 = 'CIRCexplorer2 parse -t STAR -b %s_back_spliced_junction.bed %s'%(out, Chimeric)
    cmd2 = 'CIRCexplorer2 annotate -r %s -g %s -b %s -o %s_circularRNA_known.txt'%(ref, fasta, '%s_back_spliced_junction.bed'%out, out)
    run_cmd(cmd1)
    run_cmd(cmd2)

def circRNA_finder(path, out):
    cmd = 'postProcessStarAlignment.pl --starDir %s --outDir %s'%(path, out)
    run_cmd(cmd)

def main(fq1, fq2, index, out, ref, fasta):
    STARoutDir = out
    os.makedirs('circExplorer2_out', exist_ok=True)
    os.makedirs('circRNA_finder_out', exist_ok=True)
    align(index_path=index, output=f'{STARoutDir}/{out}', fq1=fq1, fq2=fq2)
    circExplorer2(f'{STARoutDir}/{out}Chimeric.out.junction', ref, fasta, f'circExplorer2_out/{out}')
    circRNA_finder(f'{STARoutDir}', 'circRNA_finder_out')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Pipeline for circular RNA analysis using STAR, circExplorer2, and circRNA_finder')
    
    parser.add_argument('-1', '--fq1', required=True,
                        help='Input fastq file 1 (required)')
    parser.add_argument('-2', '--fq2',
                        help='Input fastq file 2 (optional, for paired-end data)')
    parser.add_argument('-i', '--index', required=True,
                        help='Path to STAR index directory')
    parser.add_argument('-o', '--out', required=True,
                        help='Output prefix')
    parser.add_argument('-r', '--ref', required=True,
                        help='Reference gene annotation file')
    parser.add_argument('-f', '--fasta', required=True,
                        help='Reference genome fasta file')

    args = parser.parse_args()

    main(args.fq1, args.fq2, args.index, args.out, args.ref, args.fasta)