#!/usr/bin/python -u

import os
import re
import sys
import glob
import argparse
import multiprocessing

import cPickle as pickle

#
def Chunk(L, n):
    for i in xrange(0, len(L), n):
        yield L[i:i+n]

#
def runcmd(x):
    print "running: {0}".format(x)

    os.system(x)

#
if __name__ == '__main__':
    #
    parser = argparse.ArgumentParser()

    parser.add_argument('--geno', type=str, required=True)
    parser.add_argument('--meth', type=str, required=True)
    parser.add_argument('--covars', type=str, required=True)
    parser.add_argument('--tp', type=str, required=True, choices=[ "cord", "F7", "15up", "antenatal", "FOM" ])
    parser.add_argument('--fo', type=argparse.FileType('w'), required=True)

    parser.add_argument('--keep-snps', type=str, default=None)
    parser.add_argument('--keep-probes', type=str, default=None)

    parser.add_argument('--pv', type=float, default=1e-5)
    parser.add_argument('--distance', type=float, default=1000000)

    # hidden argument(s) ... used for debugging purposes
    debug = parser.add_argument_group('debugging')

    debug.add_argument('--batch', type=int, default=10000, help=argparse.SUPPRESS)

    debug.add_argument('--skip-annotation', action='store_true', help=argparse.SUPPRESS)
    debug.add_argument('--keep-intermediate-data', action='store_true', help=argparse.SUPPRESS)

    ##

    Args = parser.parse_args()

    try:
        assert Args.geno, "no genotype data specified"

        os.system("rm -rf ./parsed_data")
        os.system("mkdir  ./parsed_data")


        # subset SNP(s)/ALN(s) ...
        cmd = "plink --noweb --bfile {0} --keep {1} --recodeA --freq --out ./parsed_data/ALSPAC ".format(Args.geno, Args.tp)

        if Args.keep_snps:
            cmd += " --extract {0}".format(Args.keep_snps)

        os.system(cmd)

        # batch SNP(s)
        with open("./parsed_data/ALSPAC.raw", "r") as fin:
            print "Batching SNPs ..."

            # Header
            data = [ x.split("_")[0] for x in fin.next()[:-1].split(" ")[6:] ]
            i    = 0
            for x in Chunk(data, Args.batch):
                i += 1

                with open("./parsed_data/ALSPAC." + str(i) + ".tab", "w") as OUT:
                    OUT.write("\t".join( [ "" ] + x ) + "\n")

            # Data
            for record in fin:
                ALN  = record[:-1].split(" ")[0]
                data = record[:-1].split(" ")[6:]

                i    = 0
                for x in Chunk(data, Args.batch):
                    i += 1

                    with open("./parsed_data/ALSPAC." + str(i) + ".tab", "a") as OUT:
                        OUT.write("\t".join( [ ALN ] + x ) + "\n")


        ## Pre-processing ...
        print "Pre-processing ..."

        cmd = "Rscript ./preprocessing.R --meth {0} --covars {1}".format(Args.meth, Args.covars)
        if Args.keep_probes:
            cmd += " --keep-probes {0}".format(Args.keep_probes)

        os.system(cmd)

        ## Matrix EQTL
        cmd  = []
        for p in glob.glob("./parsed_data/ALSPAC.*.tab"):
            cmd.append("Rscript ./MeQTL.R {0} ./parsed_data/probes.tab {1} {2}".format(p, p[:-3] + "meQTL", Args.pv))

        MP   = multiprocessing.Pool()
        Data = MP.map(runcmd, cmd)


        out = glob.glob("./parsed_data/ALSPAC.*.meQTL")
        if len(out) > 0:

            os.system("for x in ./parsed_data/ALSPAC.*.meQTL; do tail -n +2 $x; done > ./parsed_data/MeQTL.tsv")
            os.system("sort -g -k 5 ./parsed_data/MeQTL.tsv -o ./parsed_data/MeQTL.tsv")


            ## Annotate
            if not Args.skip_annotation:
                print "Annotating ..."

                print "reading FREQ info ..."

                frq  = {}
                with open("./parsed_data/ALSPAC.frq", "r") as fin:
                    fin.next()
                    for record in fin:
                        record = re.sub(" +", " ", record[:-1]).strip().split(" ")

                        frq[record[1]] = record # CHR, SNP, A1 (Effect), A2 (Reference), MAF, NCHROBS

                print "reading SNP/probe meta-data ..."

                SNPdata = {} # To add annotations create a dictionary with structure {<rsid>:[SNPchr, SNPpos, SNPgene]}
                CPGdata = {} # To add annotations create a dictionary with structure {<cpgid>:[CPGchr, CPGpos, CPGgene]}

                print "Writing data ..."

                Args.fo.write("\t".join([ "SNP", "SNPchr", "SNPpos", "SNPgene", "CpG", "CPGchr", "CPGpos", "CPGgene", "A1", "A2", "freq", "b", "se", "p", "N", "Trans", "r2" ]) + "\n")

                for record in open("./parsed_data/MeQTL.tsv", "r"):
                    record = record.strip().split("\t")

                    SNP, SNPchr, SNPpos, SNPgene = record[0], ".", ".", "."
                    CpG, CPGchr, CPGpos, CPGgene = record[1], ".", ".", "."

                    if SNPdata.has_key(SNP):
                        SNPchr, SNPpos, SNPgene  = SNPdata[SNP]
                    if CPGdata.has_key(CpG):
                        CPGchr, CPGpos, CPGgene  = CPGdata[CpG]

                    if not frq.has_key(SNP):
                        continue

                    A1   = frq[SNP][2]
                    A2   = frq[SNP][3]
                    freq = frq[SNP][4]

                    b  = float(record[2])
                    t  = float(record[3])
                    p  = record[4]
                    se = b / t
                    N  = int(frq[SNP][5]) / 2
                    r2 = record[6]

                    Trans = "."
                    if not SNPchr == "." and not SNPpos == "." and \
                       not CPGchr == "." and not CPGpos == ".":
                        try:

                            Trans = "N"
                            if SNPchr != CPGchr:
                                Trans = "Y"
                            if abs(int(SNPpos) - int(CPGpos)) > Args.distance:
                                Trans = "Y"
                        except:
                            Trans = "."

                    Args.fo.write("\t".join(map(str, [ SNP, SNPchr, SNPpos, SNPgene, CpG, CPGchr, CPGpos, CPGgene, A1, A2, freq, b, se, p, N, Trans, r2 ])) + "\n")
            else:
                for record in open("./parsed_data/MeQTL.tsv", "r"):
                    Args.fo.write(record)
        else:
            print "No meQTLs passing p-value threshold {0} found".format(Args.pv)
    except:
        raise
    finally:
        if not Args.keep_intermediate_data:
            os.system("rm -rf ./parsed_data")
