import json
import os
import pickle
import sys
import pandas
from zipfile import ZipFile

databases = sys.argv[1]
TAX_PICKLE = sys.argv[2]


def successors(g, O):
    """
    Extract ancestors nodes from an starting node
    :param g: starting node name
    :param O: Graph
    :return: a set with node names
    """
    result = {g}
    for o in O.successors(g):
        result.update(successors(o, O))
    return result


if __name__ == '__main__':

    tax = pickle.load(open(TAX_PICKLE, "rb"))
    print('{} taxonomies loaded'.format(len(tax.nodes())))

    taxonomy_groups = {
        'bacteria': {'taxid': '2', 'nodes': [], 'taxa': set()},
        'archaea': {'taxid': '2157', 'nodes': [], 'taxa': set()},
        'viridiplantae': {'taxid': '33090', 'nodes': [], 'taxa': set()},
        'fungi': {'taxid': '4751', 'nodes': [], 'taxa': set()},
        'arthropoda': {'taxid': '6656', 'nodes': [], 'taxa': set()},
        'chordata': {'taxid': '7711', 'nodes': [], 'taxa': set()},
        'metazoa': {'taxid': '33208', 'nodes': [], 'taxa': set()},
        'eukaryota': {'taxid': '2759', 'nodes': [], 'taxa': set()},
        'viruses': {'taxid': '10239', 'nodes': [], 'taxa': set()},
    }
    for k in taxonomy_groups:
        if k != 'metazoa' and k != 'eukaryota':
            taxonomy_groups[k]['nodes'] = [int(i) for i in successors(taxonomy_groups[k]['taxid'], tax)]
            print('{} with {} taxa'.format(k, len(taxonomy_groups[k]['nodes'])))
    k = 'metazoa'
    to_exclude = set(taxonomy_groups['arthropoda']['nodes'] + taxonomy_groups['chordata']['nodes'])
    taxonomy_groups[k]['nodes'] = [int(i) for i in successors(taxonomy_groups[k]['taxid'], tax)]
    taxonomy_groups[k]['nodes'] = list(set(taxonomy_groups[k]['nodes']) - to_exclude)
    print('{} with {} taxa'.format(k, len(taxonomy_groups[k]['nodes'])))
    k = 'eukaryota'
    to_exclude = set(list(to_exclude) + taxonomy_groups['viridiplantae']['nodes'] + taxonomy_groups['fungi']['nodes'] +
                     taxonomy_groups['metazoa']['nodes'])
    taxonomy_groups[k]['nodes'] = [int(i) for i in successors(taxonomy_groups[k]['taxid'], tax)]
    taxonomy_groups[k]['nodes'] = list(set(taxonomy_groups[k]['nodes']) - to_exclude)
    print('{} with {} taxa'.format(k, len(taxonomy_groups[k]['nodes'])))

    databases = databases.split(',')

    for db in databases:
        if not os.path.exists('{}/ncbi_dataset/data'.format(db)):
            os.makedirs('{}/ncbi_dataset/data'.format(db))
        with ZipFile('{}_meta.zip'.format(db), 'r') as zip:
            assemblies = []
            assemblies_tmp = {}
            with zip.open('ncbi_dataset/data/assembly_data_report.jsonl') as fjson, open(
                    '{}/ncbi_dataset/data/assembly_data_report.jsonl'.format(db), 'w') as fjson_out:
                for line in fjson.readlines():
                    d = json.loads(line.decode("utf-8"))
                    v = assemblies_tmp.setdefault(d['taxId'], [])
                    v.append(d)
                for s in assemblies_tmp.keys():
                    rep_genome = []
                    for e in assemblies_tmp[s]:
                        if 'refseqCategory' in e['assemblyInfo']:
                            rep_genome.append(e)
                    if len(rep_genome) == 1:
                        assemblies.append(rep_genome[0]['assemblyInfo']['assemblyAccession'])
                        fjson_out.write('{}\n'.format(json.dumps(rep_genome[0])))
                    else:
                        assemblies.append(assemblies_tmp[s][0]['assemblyInfo']['assemblyAccession'])
                        fjson_out.write('{}\n'.format(json.dumps(assemblies_tmp[s][0])))
                    for k in taxonomy_groups:
                        if int(s) in taxonomy_groups[k]['nodes']:
                            taxonomy_groups[k]['taxa'].add(int(s))
                            break

            print('There are {} assemblies included'.format(len(assemblies)))
            with zip.open('ncbi_dataset/data/dataset_catalog.json') as fjson, open(
                    '{}/ncbi_dataset/data/dataset_catalog.json'.format(db), 'w') as fjson_out:
                d = json.loads(fjson.read().decode("utf-8"))
                catalog = []
                for c in d['assemblies']:
                    if 'accession' in c:
                        if c['accession'] in assemblies:
                            catalog.append(c)
                    else:
                        catalog.append(c)
                d['assemblies'] = catalog
                fjson_out.write(json.dumps(d, indent=2))
            with zip.open('ncbi_dataset/fetch.txt') as fin, open('{}/ncbi_dataset/fetch.txt'.format(db), 'w') as fout:
                for line in fin.readlines():
                    line = line.decode("utf-8")
                    f = os.path.dirname(line.split('\t')[2].replace('data/', ''))
                    if f in assemblies:
                        fout.write(line)

    data = []
    for k in taxonomy_groups:
        data.append([k, len(taxonomy_groups[k]['nodes']), len(taxonomy_groups[k]['taxa'])])
    df = pandas.DataFrame(data, columns=['Tax Group', 'All taxa', 'Included Taxa'])
    df.to_csv('taxa.tsv', index=False, sep='\t')