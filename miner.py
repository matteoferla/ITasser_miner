#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__doc__ = \
    """
A tool to spy on I-TASSER queue.
Download what has been done match it to human protein and see if
* any can be used
* Model for cutting

NB. Written for python 3, not tested under 2.
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "28 Feb 2018"
__license__ = "Cite me!"
__copyright__ = 'GNU'
__version__ = "0"

import argparse, os, csv, time
from pprint import PrettyPrinter
import tarfile
import json
import urllib.request
import urllib.error
from warnings import warn
from Bio import PDB

pprint = PrettyPrinter().pprint

import sys
from warnings import warn

if sys.version_info[0] < 3:
    warn("Sorry man, I told you to use Python 3.")

################## PDB monkeypatch
def get_structure_from_stream(self, id, stream):
    self.header = None
    self.trailer = None
    # Make a StructureBuilder instance (pass id of structure as parameter)
    self.structure_builder.init_structure(id)
    self._parse([x.decode("utf-8") for x in stream.readlines()])
    self.structure_builder.set_header(self.header)
    # Return the Structure instance
    structure = self.structure_builder.get_structure()
    return structure

PDB.PDBParser.get_structure_from_stream = get_structure_from_stream

#####################

"""
* fetch_id(i) - downaloads identifier tarball unless class attribute .can_download is false. Returns self
Getters that return a value or None:
* get_all_seqs(from,to) - given the first and last id, returns a dictionary of ids:seqs
* get_seq_from_id(i) - given an id gets the sequence (from model1.pdb)
* get_seq_from_tar(t) - given a tarball returns a sequence
* get_tar_from_id(i) - given an id returns a tarball
* get_accession_from_seq(s) - given a sequence gets a uniprot accession, requires blast
"""

class Itasser():
    if os.path.isfile('missing.json'):
        invalid = json.load(open('missing.json'))
    else:
        invalid={}
    can_download = True
    verbose = True

    def fetch_id(self, identifier):
        if not self.can_download:
            return None
        url = 'https://zhanglab.ccmb.med.umich.edu/I-TASSER/output/{id}/{id}_results.tar.bz2'.format(id=identifier)
        filename = 'tarballs/{id}.tar.bz2'.format(id=identifier)
        try:
            with urllib.request.urlopen(url) as response:
                with open(filename, 'wb') as w:
                    w.write(response.read())
        except urllib.error.HTTPError as error:
            self.invalid[identifier] = str(error)
            self.store_invalid()
            warn(str(error))
        return self


    def get_seq_from_tar(self, tar):
        for i in tar:
            try:
                stream = tar.extractfile('model1.pdb')
            except Exception as err:
                warn(str(err))
                stream = None
            else:
                if stream: #io.BufferedReader
                    p = PDB.PDBParser().get_structure_from_stream('model', stream)
                    b = PDB.PPBuilder()
                    pp=b.build_peptides(p)[0]
                    tar.close()
                    return pp.get_sequence()
                else:
                    if self.verbose: print('IO Stream is None from tarball.')


    def get_tar_from_id(self,identifier):
        filename = 'tarballs/{id}.tar.bz2'.format(id=identifier)
        if not os.path.isfile(filename):
            self.fetch_id(identifier)
        if identifier in self.invalid:
            if self.verbose: print('\tpreviously excluded.')
            return None
        if not os.path.isfile(filename):
            if self.verbose: print('\tFile absent.')
            return None
        try:
            tar = tarfile.open(filename, "r:bz2")
        except Exception as err:
            warn(str(err))
            os.remove(filename) # not a tar? It was a 403!
            self.invalid[identifier] = str(err)
            self.store_invalid()
            return None
        return tar

    def get_seq_from_id(self, identifier):
        tar= self.get_tar_from_id(identifier)
        if tar:
            return self.get_seq_from_tar(tar)
        else:
            return None

    def iter_ids(self,earliest, latest):
        for i in range(int(earliest[1:]),int(latest[1:])+1):
            yield 'S'+str(i)


    def get_all_seqs(self, latest, earliest):
        data={}
        for identifier in self.iter_ids(earliest, latest):
            if self.verbose: print(identifier)
            seq = self.get_seq_from_id(identifier)
            time.sleep(2)
            if seq:
                data[identifier]= seq
                if self.verbose: print('\thas sequence.')
        return data


    def store_invalid(self):
        json.dump(self.invalid,open('missing.json','w'))

    def get_human_accession_from_seq(self,seq):
        if seq:
            pass
        else:
            return None



######################### MAIN

def part_one():
    it = Itasser()
    it.verbose = True
    it.can_download = True
    os.mkdir('fasta')
    try:
        data = it.get_all_seqs(latest='S449979', earliest='S438773')
        import csv
        with open('data.csv', 'w', newline='') as w:
            sheet = csv.writer(w)
            sheet.writerow(['identifier', 'sequence'])
            sheet.writerows(data.items())
        for k in data:
            with open('fasta/'+k + '.fa', 'w') as w:
                w.write('>{i}\n{s}\n\n'.format(s=data[k], i=k))
    except KeyboardInterrupt:
        it.store_invalid()
        print('User exited')
    else:
        it.store_invalid()
        print('finished')

def part_two():
    data=list(csv.DictReader(open('data.csv')))
    l = [len(s['sequence']) for s in data]
    print(l)
    w = csv.DictWriter(open('data_extra.csv','w',newline=''),fieldnames=('identifier', 'sequence', 'len', 'sheet', 'turn', 'helix'))
    w.writeheader()
    from Bio.SeqUtils import ProtParam
    for gene in data:
        (h,t,s)=ProtParam.ProteinAnalysis(gene['sequence']).secondary_structure_fraction()
        w.writerow({'identifier': gene['identifier'],
                    'sequence': gene['sequence'],
                    'len':len(gene['sequence']),
                    'sheet': s,
                   'turn':t,
                    'helix':h})



def part_blast():
    import csv, os
    #os.system('ncbi-blast-2.8.1+/bin/makeblastdb -in human.fa -parse_seqids -dbtype prot -taxid 9606')
    from Bio.Blast import NCBIXML
    def _get(blast_record):
        for align in blast_record.alignments:
            for hsp in align.hsps:
                d = {'identifier': infile.replace('.fa', ''), 'sequence': 'x', 'match': align.title[0:50], 'match_score': hsp.score, 'match_start': hsp.query_start, 'match_length': hsp.align_length,
                     'match_identity': hsp.identities / hsp.align_length}
                return d
        else:
            return None
    w=csv.DictWriter(open('matched.csv','w',newline=''), fieldnames=('identifier','sequence','match','match_score','match_start','match_length','match_identity'))
    w.writeheader()
    for file in os.listdir('fasta'):
        try:
            if '.fa' not in file:
                continue
            print(file.replace('.fa',''))
            infile = 'fasta/'+file
            blastfile= infile.replace('.fa', '_blast.xml')
            os.system('ncbi-blast-2.8.1+/bin/blastp -query {infile} -db human.fa -outfmt 5 -num_threads 6 > {outfile}'.format(infile=infile, outfile=blastfile))
            blast_record = NCBIXML.read(open(blastfile))
            d = _get(blast_record)
            if d:
                w.writerow(d)
        except Exception as err:
            warn('Caught: '+str(err)+' for file '+file)







if __name__ == "__main__":
    ### ARGPARSER
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('--verbose', action='store_true', help='Runs giving details')
    args = parser.parse_args()
    #part_one()
    part_two()
    #part_blast()
