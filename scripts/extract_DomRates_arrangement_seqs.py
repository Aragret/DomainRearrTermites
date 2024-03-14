#!/usr/bin/python3

import argparse
from collections import defaultdict, namedtuple
from itertools import groupby
import os
import re

__author__ = "Elias Dohmen"
__version__ = "0.2"
__email__ = "e.dohmen@wwu.de"
__institute__ = "IEB Münster"

def extract_rearrangements(domRatesFile: str, nodeids: str, solutions: str, events: str) -> set:
    """
    Extracts all rearrangement events that match the given restrictions and returns these as a list of namedtuples
    :param domRatesFile:
    :param node: node-ID in the phylogenetic tree to consider
    :param solutions: all solution types to consider (comma separated str)
    :param events: all event types to consider (comma separated str)
    :return: set of namedtuples with node, solution, event and arrangement
    """
    DRentry = namedtuple('DRentry', 'node solution event arrangement')
    nodeIDs = nodeids.strip().split(',')
    drset = set()
    solutions_lst = solutions.strip().split(',')
    events_lst = events.strip().split(',')

    with open(domRatesFile, 'r') as drf:
        for raw_line in drf:
            if raw_line[0] not in ('#', '\n'):
                line = raw_line.strip().split('\t')
                if line[0] in nodeIDs and line[1] in solutions_lst and line[2] in events_lst:
                    if line[2] == 'fusion' or line[2] == 'terminal emergence' or line[2] == 'terminal loss':
                        node, solution, event, new_arrangement, old = line
                        drset.add(DRentry(node=node.strip(), solution=solution.strip(), event=event.strip(),
                                          arrangement='-'.join(new_arrangement.strip().split())))
                    elif line[2] == 'fission':
                        node, solution, event, new_arrangement, old = line
                        new_arrgmnts = new_arrangement.strip().split('|')
                        drset.add(DRentry(node=node.strip(), solution=solution.strip(), event=event.strip(),
                                          arrangement='-'.join(new_arrgmnts[0].strip().split())))
                        drset.add(DRentry(node=node.strip(), solution=solution.strip(), event=event.strip(),
                                          arrangement='-'.join(new_arrgmnts[1].strip().split())))
                    elif line[2] == 'single domain emergence':
                        node, solution, event, new_arrangement = line
                        drset.add(DRentry(node=node.strip(), solution=solution.strip(), event=event.strip(),
                                          arrangement='-'.join(new_arrangement.strip().split())))
                    elif line[2] == 'single domain loss':
                        continue

    return drset


def extract_seqids(pfam_scan: str) -> defaultdict:
    """
    Extracts all collapsed domain arrangements and their PfamIDs
    :param pfam_scan: input PfamScan annotation file (string)
    :return: defaultdict(list) with domain arrangement str as keys and list of SequenceIDs as values
    """
    all_arr_dict = defaultdict(list)
    transp_dict = defaultdict(list)

    with open(pfam_scan, 'r') as inf:
        for raw_line in inf:
            if raw_line[0] not in ('#', '\n'):
                line = raw_line.strip().split()
                seq_id = line[0]
                domain_start_pos = int(line[1])
                domain_name = re.match('([^.]+)', line[5]).group(1)
                all_arr_dict[seq_id].append((domain_start_pos, domain_name))

    for seq_id, doms in all_arr_dict.items():
        uncollapsed_arrangement = [x[1] for x in sorted(doms)]
        # collapses consecutive repeats of domains as it is done in DomRates
        arrangement = [k for k, g in groupby(uncollapsed_arrangement)]
        arr_str = '-'.join(arrangement)
        transp_dict[arr_str].append(seq_id)

    return transp_dict


def match_seqs(domRates_set, pf_dir, extens):
    """
    Loops through all PfamScan files with given extension in given directory and matches the domain arrangements from domRates_set to SequenceIDs with same arrangements
    :param domRates_set: set of namedtuples with domRates events (+solutions+arrangements)
    :param pf_dir: directory containing pfamScan files
    :param extens: file name extension of pfamScan files
    :return: defaultdict of defaultdict of list (key1=pfamScanFileName: key2=namedtuple(nodeID,solution,event,arrangement): list of sequenceIDs having this arrangement
    """
    summ_dict = defaultdict(lambda: defaultdict(list))
    if pf_dir == None:
        pf_dir = os.getcwd()
    file_lst = [fi for fi in os.listdir(pf_dir) if fi.endswith(extens)]

    for pf_file in file_lst:
        match_dict = extract_seqids(os.path.join(pf_dir, pf_file))
        for entry in domRates_set:
            summ_dict[pf_file][entry].extend(match_dict[entry.arrangement])

    return summ_dict


def summarise(summary_dict, output):
    """
    Writes a summary to given output file or prints to stdout if no file provided
    :param summary_dict: Dictionary containing output to summarise (see match_seqs() return)
    :param output: output file name
    """
    if output == None:
        print('# filename\tnodeID\tsolutionType\tevenType\tarrangement\tsequenceIDs')
        for pffile, ntdict in summary_dict.items():
            for evs, seqids in ntdict.items():
                if seqids:
                    print(pffile + '\t' + evs.node + '\t' + evs.solution + '\t' + evs.event + '\t' + evs.arrangement +
                          '\t' + ','.join(seqids))
    else:
        with open(output, 'w') as of:
            of.write('# filename\tnodeID\tsolutionType\tevenType\tarrangement\tsequenceIDs')
            for pffile, ntdict in summary_dict.items():
                for evs, seqids in ntdict.items():
                    if seqids:
                        of.write(pffile + '\t' + evs.node + '\t' + evs.solution + '\t' + evs.event + '\t' +
                                 evs.arrangement + '\t' + ','.join(seqids)+'\n')


def main():
    """
    The main function of the program.
    """
    parser = argparse.ArgumentParser(description='Extracts all rearrangement events from a given node of DomRates output and matches sequenceIDs from pfam_scan files to it.')
    parser.add_argument("-f", "--domRatesFile", type=str, required=True, help="The DomRates output file to extract rearrangement events from.")
    parser.add_argument("-n", "--node", type=str, required=True, help="The node-ID(s) for which rearrangements should been extracted (if multiple nodeIDs, provide as comma-separated list).")
    parser.add_argument("-s", "--solutions", type=str, default='exact solution,non-ambiguous solution', help="Comma separated list of solution types to consider.")
    parser.add_argument("-e", "--events", type=str, default='fusion,fission,single domain emergence,single domain loss,terminal emergence,terminal loss', help="Comma separated list of event types to consider")
    parser.add_argument("-d", "--pfam_scan_dir", default=None, type=str, help="The directory containing PfamScan files to match the sequenceIDs to arrangements (if not provided searches in current directory).")
    parser.add_argument("-x", "--extension", type=str, default='.dom', help="The file extension of the PfamScan files (see -d parameter).")
    parser.add_argument("-o", "--output", type=str, default=None, help="The name of the output file (if not provided printed to stdout).")
    args = parser.parse_args()

    drset = extract_rearrangements(args.domRatesFile, args.node, args.solutions, args.events)
    summ_dict = match_seqs(drset, args.pfam_scan_dir, args.extension)
    summarise(summ_dict, args.output)


if __name__ == "__main__":
    main()
