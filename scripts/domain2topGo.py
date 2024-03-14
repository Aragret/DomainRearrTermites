#!/usr/bin/env python3

# edited by Elias Dohmen
# compatibility for DomRates version 1.0.0 added

import argparse
import re
import sys


def read_Pfam(infile, pfam2go, gene2go):
    try:
        with open(infile) as f:
            data = f.readlines()
    except Exception as e:
        sys.exit("failed to open file: %s" % (str(e)))
    for line in data:
        line = line.strip()
        if (len(line) == 0) or (line[0] == "#"):
            continue
        tokens = re.split(" +", line)
        pfam_id = tokens[5].split(".")[0]
        if pfam_id in pfam2go:
            if tokens[0] not in gene2go:
                gene2go[tokens[0]] = set()
            for elem in pfam2go[pfam_id]:
                gene2go[tokens[0]].add(elem)


def read_DAMA(infile, pfam2go, gene2go):
    try:
        with open(infile) as f:
            data = f.readlines()
    except Exception as e:
        sys.exit("failed to open file: %s" % (str(e)))
    for line in data:
        line = line.strip()
        if (len(line) == 0) or (line[0] == "#"):
            continue
        tokens = line.split()
        pfam_id = tokens[4]
        if pfam_id in pfam2go:
            if tokens[3] not in gene2go:
                gene2go[tokens[3]] = set()
            for elem in pfam2go[pfam_id]:
                gene2go[tokens[3]].add(elem)

def read_InterPro(infile, gene2go):
    try:
        with open(infile) as f:
            data = f.readlines()
    except Exception as e:
        sys.exit("failed to open file: %s" % (str(e)))
    for line in data:
        line = line.strip()
        if (len(line) == 0) or (line[0] == "#"):
            continue
        tokens = line.split("\t")
        if len(tokens) >= 14:
            if tokens[0] not in gene2go:
                gene2go[tokens[0]] = set()
            go_terms = tokens[13].split("|")

def read_domRates_statistics(infile, node, pfam2go, gene2go):
    try:
        with open(infile) as f:
            data = f.readlines()
    except Exception as e:
        sys.exit("failed to open file: %s" % (str(e)))
    change_list = []
    for line in data:
        line = line.strip()
        if (len(line) == 0) or (line[0] == "#"):
            continue
        tokens = line.split("\t")
        if (node is None) or (node == tokens[0]):
            name = tokens[0]+tokens[1].strip().replace(" ", "")
            anz=0
            for token in tokens[2:]:
                if token != "0":
                    anz += 1
            domains = tokens[1].split(" ")
            go_terms = set()
            for pfam_id in domains:
                if pfam_id in pfam2go:
                    for elem in pfam2go[pfam_id]:
                        go_terms.add(elem)
            if len(go_terms) != 0:
                gene2go[name] = go_terms
                if (tokens[-1] == "0") and (anz==1) and (tokens[6] == "0"):
                    change_list.append(name)
    return change_list

def read_domRates_v100_statistics(infile, node, pfam2go, gene2go):
    try:
        with open(infile) as f:
            data = f.readlines()
    except Exception as e:
        sys.exit("failed to open file: %s" % (str(e)))
    change_list = []
    for line in data:
        name_list = []
        arr_list = []
        line = line.strip()
        if (len(line) == 0) or (line[0] == "#"):
            continue
        tokens = line.split("\t")
        if (node is None) or (node == tokens[0]):
            arr = None
            if ("loss" in tokens[2]):
                arr = " ".join(re.findall(r'PF\d{5}', tokens[4]))
                arr_list.append(arr)
            elif ("fission" != tokens[2]):
                arr = " ".join(re.findall(r'PF\d{5}', tokens[3]))
                arr_list.append(arr)
            if (arr is not None):
                name = tokens[0]+arr.strip().replace(" ", "")
                name_list.append(name)
            if (tokens[2] == "fission"):
                arrlist = tokens[3].split("|")
                arr1 = " ".join(re.findall(r'PF\d{5}', arrlist[0]))
                arr2 = " ".join(re.findall(r'PF\d{5}', arrlist[1]))
                arr_list.append(arr1)
                arr_list.append(arr2)
                name1 = tokens[0]+arr1.strip().replace(" ", "")
                name2 = tokens[0] + arr2.strip().replace(" ", "")
                name_list.append(name1)
                name_list.append(name2)

            for name, arr in zip(name_list, arr_list):
                domains = arr.split(" ")
                go_terms = set()
                for pfam_id in domains:
                    if pfam_id in pfam2go:
                        for elem in pfam2go[pfam_id]:
                            go_terms.add(elem)
                if len(go_terms) != 0:
                    gene2go[name] = go_terms
                    if tokens[1] in ("exact solution", "non-ambiguous solution"):
                        change_list.append(name)
    return change_list

def read_pfam2go(infile):
    try:
        with open(infile) as f:
            data = f.readlines()
    except Exception as e:
        sys.exit("failed to open file: %s" % (str(e)))
    dom2go = {}
    for line in data:
        if line[0] == "!":
            continue
        line = line.strip()
        tokens = line.split(" ")
        pfam_id = tokens[0].split(":")[1]
        if pfam_id in dom2go:
            dom2go[pfam_id].append(tokens[-1])
        else:
            dom2go[pfam_id] = [tokens[-1]]


    return dom2go

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p','--pfam', help='Pfam domain annotations', nargs='+', metavar='FILE')
    parser.add_argument('-i','--interpro', help='InterPro domain annotations', nargs='+', metavar='FILE')
    parser.add_argument('-d','--dama', help='DAMA domain annotations', nargs='+', metavar='FILE')
    parser.add_argument('-r','--domRates', help='DomRates statistics', metavar='FILE')
    parser.add_argument('-s', '--domRates2', help='DomRates statistics version >= 1.0.0', metavar='FILE')
    parser.add_argument('-g','--pfam2go', help='pfam2go annotation', metavar='FILE')
    parser.add_argument('-n','--node', help='The node to analyse', metavar='FILE')
    parser.add_argument('-o','--out', help='The output file', metavar='FILE', required=True)
    parser.add_argument('-c','--change-out', help='The file with changed arrangements', metavar='FILE')
    args = parser.parse_args()

    if ((args.pfam is not None) or (args.dama is not None)) and (args.pfam2go is None):
        print("Error! You need to provide the pfam2go file!", file=sys.stderr)

    dom2go = {}
    if args.pfam2go is not None:
        dom2go = read_pfam2go(args.pfam2go);

    gene2go = {}
    if args.pfam is not None:
        for elem in args.pfam:
            read_Pfam(elem, dom2go, gene2go)
    if args.dama is not None:
        for elem in args.dama:
            read_DAMA(elem, dom2go, gene2go)

    if args.interpro is not None:
        for elem in args.interpro:
            read_InterPro(elem, gene2go)
    if args.domRates is not None:
        changes = read_domRates_statistics(args.domRates, args.node, dom2go, gene2go)
        with open(args.change_out, "w") as f:
            for elem in changes:
                f.write(elem + "\n")
    if args.domRates2 is not None:
        changes = read_domRates_v100_statistics(args.domRates2, args.node, dom2go, gene2go)
        with open(args.change_out, "w") as f:
            for elem in changes:
                f.write(elem + "\n")

    try:
        with open(args.out, "w") as f:
            for key, val in gene2go.items():
                f.write(key + "\t" + ', '.join(sorted(val)) + "\n")
    except Exception as e:
        sys.exit("failed to open file: %s" % (str(e)))



if __name__ == "__main__":
    main()
