import sys, os

def merge_algorithm_res(files, tissue, species, algorithm_name):
    dic = {}
    for sample in files:
        with open(sample) as f:
            for i in f.readlines()[1:]:
                line = i.strip().split('\t')
                circRNA = line[0] + ':' + line[1] + '|' + line[2]
                strand = line[5] if algorithm_name == 'circExplorer2' else 'NA'
                gene = line[14] if algorithm_name == 'circExplorer2' else 'NA'
                circType = line[13] if algorithm_name == 'circExplorer2' else 'NA'
                count = line[12] if algorithm_name == 'circExplorer2' else line[4]
                geneType = 'mRNA' if algorithm_name == 'circExplorer2' else 'NA'
                Conservation = 'Y'
                Agerelated = 'N'
                tissue = tissue
                dic.setdefault(circRNA, {}).setdefault('circinfo', [species, tissue, strand, gene, circType, geneType, Conservation, algorithm_name, Agerelated])
                dic.setdefault(circRNA, {}).setdefault('samples', []).append([os.path.basename(sample), count])
    res = []
    for k in dic.keys():
        output = ''
        output += k + '\t'
        output += '\t'.join(dic[k]['circinfo']) + '\t'
        counts = ''
        for l in dic[k]['samples']:
            counts += l[0] + ',' + l[1] + '|'
        output += counts.strip(',')
        res.append(output)
    return res

def ReadFile(lis):
    dic = {}
    for i in lis:
        line = i.strip().split('\t')
        circRNA, Species, Tissue, Strans, Gene, circType, GeneType, Conservation, Algorithm, Agerelated, Counts = line
        dic.setdefault(circRNA, [circRNA, Species, Tissue, Strans, Gene, circType, GeneType, Conservation, Algorithm, Agerelated, Counts])
    return dic


def format_counts(counts, algorithm_name, suffix):
    formatted_counts = []
    for entry in counts.strip('|').split('|'):
        temp = entry.replace(f'{algorithm_name}/', '').replace(suffix, '').split(',')
        formatted_counts.append(f'{temp[0]},{algorithm_name},{temp[1]}')
    return formatted_counts

def merge_data(circRNA_id, circExplorer2, circRNA_finder):
    circRNA_info = {
        'circRNA': circRNA_id,
        'Species': set(),
        'Tissue': set(),
        'Strans': 'NA',
        'Gene': 'NA',
        'circType': 'NA',
        'GeneType': 'NA',
        'Conservation': 'Y',
        'Algorithm': [],
        'Agerelated': 'N',
        'Counts': []
    }
    
    if circRNA_id in circExplorer2:
        circRNA_info['Species'].add(circExplorer2[circRNA_id][1])
        circRNA_info['Tissue'].add(circExplorer2[circRNA_id][2])
        circRNA_info['Algorithm'].append('circExplorer2')
        circRNA_info['Counts'].extend(format_counts(circExplorer2[circRNA_id][-1], 'circExplorer2', '__circularRNA_known.txt'))
        
    if circRNA_id in circRNA_finder:
        circRNA_info['Species'].add(circRNA_finder[circRNA_id][1])
        circRNA_info['Tissue'].add(circRNA_finder[circRNA_id][2])
        circRNA_info['Algorithm'].append('circRNA_finder')
        circRNA_info['Counts'].extend(format_counts(circRNA_finder[circRNA_id][-1], 'circRNA_finder', '_s_filteredJunctions.bed'))
    
    circRNA_info['Species'] = ','.join(circRNA_info['Species'])
    circRNA_info['Tissue'] = ','.join(circRNA_info['Tissue'])
    circRNA_info['Algorithm'] = ','.join(circRNA_info['Algorithm'])
    circRNA_info['Counts'] = '|'.join(circRNA_info['Counts'])
    
    return circRNA_info

def get_files_from_dir(directory, suffix):
    """
    Returns a list of file paths in a directory that end with the specified suffix.
    """
    return [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(suffix)]


def main():
    circExplorer2 = ReadFile(merge_algorithm_res(get_files_from_dir(sys.argv[3], 'circularRNA_known.txt'), sys.argv[1], sys.argv[2], 'circExplorer2'))
    circRNA_finder = ReadFile(merge_algorithm_res(get_files_from_dir(sys.argv[4], 's_filteredJunctions.bed'), sys.argv[1], sys.argv[2], 'circRNA_finder'))

    all_ids = set(circExplorer2.keys()).union(set(circRNA_finder.keys()))

    for circRNA_id in all_ids:
        merged_data = merge_data(circRNA_id, circExplorer2, circRNA_finder)
        print('\t'.join([
            merged_data['circRNA'], merged_data['Species'], merged_data['Tissue'],
            merged_data['Strans'], merged_data['Gene'], merged_data['circType'],
            merged_data['GeneType'], merged_data['Conservation'], merged_data['Algorithm'],
            merged_data['Agerelated'], merged_data['Counts']
        ]))

if __name__ == "__main__":
    main()
