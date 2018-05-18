"""
Converts a PDB into a modified PDB format, where the atom type for each atom is indicated
This reduced PDB is what is understood by ATTRACT.
These atoms are parsed from the trans files and topology files derived from the OPLS forcefield
"""
import sys, os
sys.path.append('/home/glenn/Documents/Masterthesis/attract_unchanged/allatom')
import parse_cns_top
import pdbcomplete
from pdbcomplete import run_pdb2pqr, run_whatif, pdbfix, update_patches, pdb_lastresort

has_argparse = False
try:
    import argparse

    has_argparse = True
except ImportError:
    import optparse  # Python 2.6



# this script is a version of aareduce that is not dependend on parsing but can be treated
# as a stand alone version with indivudual functions such that is can be easily imported and used
# Mapping of nucleic-acid codes to DNA/RNA
mapnuc = {
    "A": ["DA", "RA"],
    "A3": ["DA", "RA"],
    "A5": ["DA", "RA"],
    "ADE": ["DA", "RA"],
    "C": ["DC", "RC"],
    "C3": ["DC", "RC"],
    "C5": ["DC", "RC"],
    "CYT": ["DC", "RC"],
    "G": ["DG", "RG"],
    "G3": ["DG", "RG"],
    "G5": ["DG", "RG"],
    "GUA": ["DG", "RG"],
    "T": ["DT", None],
    "T3": ["DT", None],
    "T5": ["DT", None],
    "THY": ["DT", None],
    "U": [None, "RU"],
    "U3": [None, "RU"],
    "U5": [None, "RU"],
    "URA": [None, "RU"],
    "URI": [None, "RU"],
}
mapnucrev = {
    "DA": "A",
    "RA": "A",
    "DC": "C",
    "RC": "C",
    "DG": "G",
    "RG": "G",
    "DT": "T",
    "RU": "U",
}


class PDBres:
    def __init__(self, chain, resid, resname, topology):
        self.chain = chain
        self.resid = resid
        self.resname = resname
        self.coords = {}
        self.chainfirst = False
        self.chainlast = False
        self.nter = False
        self.cter = False
        self.topology = topology


#code_to_type = {}


def parse_transfile(transfile, topname, code_to_type):
    for l in open(transfile):
        ll = l.split()
        type = int(ll[0])
        for code in ll[3:]:
            if code.startswith("#"): break
            assert (code, topname) not in code_to_type, code
            code_to_type[code, topname] = type


#mutations = {}


def read_filelist(filelist):
    ret = []
    for l in open(filelist):
        l = l.strip()
        if not len(l): continue
        assert len(l.split()) == 1, (filelist, l)
        ret.append(l)
    return ret


def read_pdb(pdblines,mutations, add_termini=False, modbase=False, modres=False, args_dna = False, args_rna = False):
    repl = (
        ("H", "HN"),
        ("HT1", "HN"),
        ("OP1", "O1P"),
        ("OP2", "O2P"),
        ("H1", "HN"),
        ("OP3", "O5T"),
        ("HO5'", "H5T"),
        ("HO3'", "H3T"),
    )
    topres, toppatch = parse_cns_top.residues, parse_cns_top.presidues
    pdbres = []
    curr_res = None
    atomlines = []

    if (modbase or modres):
        res0 = {}
        pdblines = list(pdblines)
        for l in pdblines:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                resid = l[22:27]
                if resid not in res0: res0[resid] = set()
                atomcode = l[12:16].strip()
                res0[resid].add(atomcode)
        res_ok = set()
        for r in res0:
            ratoms = res0[r]
            if modres and "CA" in ratoms and "C" in ratoms and "N" in ratoms:
                res_ok.add(r)
            elif modbase:
                cp = 0
                for n in range(5):
                    if ("C%d'" % n) in ratoms: cp += 1
                if cp >= 3:
                    res_ok.add(r)
        for l in pdblines:
            if l.startswith("ATOM"):
                atomlines.append(l)
            elif l.startswith("HETATM"):
                resid = l[22:27]
                if resid in res_ok:
                    atomlines.append(l)
    else:
        atomlines = [l for l in pdblines if l.startswith("ATOM")]

    for l in atomlines:
        atomcode = l[12:16].strip()
        if l[16] not in (" ", "A"): continue  # only keep the first of alternative conformations
        if l[30:38] == " XXXXXXX": continue  # missing atom from --manual mode
        resname = l[17:20].strip()
        if resname == "HIE" or resname == "HIP": resname = "HIS"
        if resname in mutations: resname = mutations[resname]
        if resname in mapnuc:
            if args_dna:
                resname = mapnuc[resname][0]
            elif args_rna:
                resname = mapnuc[resname][1]
            else:
                raise ValueError(
                    "PDB contains a nucleic acid named \"%s\", but it could be either RNA or DNA. Please specify the --dna or --rna option" % resname)

            if resname is None:
                if args_dna: na = "DNA"
                if args_rna: na = "RNA"
                raise ValueError("'%s' can't be %s" % (l[17:20].strip(), na))
        chain = l[21]
        resid = l[22:27]
        x = float(l[30:38])
        y = float(l[38:46])
        z = float(l[46:54])
        newres = False
        nter = False
        chainfirst = False
        if curr_res is None:
            newres = True
            chainfirst = True
            if add_termini: nter = True
        elif chain != curr_res.chain:
            newres = True
            chainfirst = True
            curr_res.chainlast = True
            if add_termini:
                nter = True
                curr_res.cter = True
        elif resid != curr_res.resid or resname != curr_res.resname:
            newres = True
        if newres:
            try:
                if resname is None: raise KeyError
                topr = topres[resname.lower()].copy()
            except KeyError:
                raise KeyError("Residue type %s not known by the topology file" % resname)
            curr_res = PDBres(chain, resid, resname, topr)
            if chainfirst: curr_res.chainfirst = True
            if nter: curr_res.nter = True
            pdbres.append(curr_res)
        curr_res.coords[atomcode] = (x, y, z)
        for pin, pout in repl:
            if atomcode != pin: continue
            curr_res.coords[pout] = (x, y, z)
    if curr_res is not None:
        curr_res.chainlast = True
        if add_termini:
            curr_res.cter = True
    return pdbres


def termini_pdb(pdbres, nter, cter):
    xter = nter, cter
    for n in range(2):
        ter = xter[n]
        for resnr in ter:
            r = [res for res in pdbres if res.resid == resnr]
            if len(r) == 0:
                raise ValueError("Cannot find residue %d" % resnr)
            elif len(r) > 1:
                raise ValueError("Multiple residues %d" % resnr)
            res = r[0]
            if n == 0:
                res.nter = True
            else:
                res.cter = True


def patch_pdb(pdbres, patches):
    topres, toppatch = parse_cns_top.residues, parse_cns_top.presidues
    for res in pdbres:
        if res.resid in patches:
            for p in patches[res.resid]:
                if p is None: continue
                res.topology.patch(toppatch[p])
        elif len(pdbres) > 1 and "ca" in res.topology.atomorder:  # protein
            if res.nter:
                if res.resname == "PRO":
                    res.topology.patch(toppatch["prop"])
                else:
                    res.topology.patch(toppatch["nter"])
            if res.cter:
                res.topology.patch(toppatch["cter2"])
        elif len(pdbres) > 1 and "p" in res.topology.atomorder:  # DNA/RNA
            if res.chainfirst:
                if res.nter:
                    res.topology.patch(toppatch["5pho"])
                else:
                    res.topology.patch(toppatch["5ter"])
            if res.chainlast:
                res.topology.patch(toppatch["3ter"])


def check_pdb(pdbres, heavy=False):
    for res in pdbres:
        top = res.topology
        for a in top.atomorder:
            atom = top.atoms[a]
            if a.lower().startswith("h"):
                if heavy: continue
                if atom.charge == 0: continue
            aa = a.upper()
            if aa.strip() not in res.coords:
                raise ValueError('Missing coordinates for atom "%s" in residue %s %s%s' % (
                aa.strip(), res.resname, res.chain, res.resid))


def write_pdb(pdbres, chain,code_to_type, heavy=False, one_letter_na=False, args_startatom = 1, args_startres = 1):
    pdblines = []
    mapping = []
    atomcounter = args_startatom
    rescounter = args_startres
    for res in pdbres:
        top = res.topology
        for a in top.atomorder:
            atom = top.atoms[a]
            if a.lower().startswith("h"):
                if heavy: continue
                if atom.charge == 0: continue
            aa = a.upper()
            x = " XXXXXXX"
            y = x;
            z = x
            if aa.strip() in res.coords:
                x, y, z = ("%8.3f" % v for v in res.coords[aa.strip()])
            xyz = x + y + z
            type = code_to_type[atom.type.upper(), top.topname]
            a0 = aa
            if len(a0) < 4:
                a0 = " " + a0 + "   "[len(a0):]
            resname = res.resname
            atomname = a0
            if one_letter_na and resname in mapnucrev:
                resname = mapnucrev[resname]
            pdblines.append("ATOM%7d %4s %3s %s%4d    %s %4d %7.3f 0 1.00" % \
                            (atomcounter, atomname, resname, chain, rescounter, xyz, type, atom.charge))
            atomcounter += 1
        mapping.append((res.resid, rescounter))
        rescounter += 1
    return pdblines, mapping


def set_reference(pdbres, pdbreferes):
    if len(pdbres) != len(pdbreferes):
        raise ValueError(
            "PDB and reference do not have the same number of residues, %d vs %s" % (len(pdbres), len(pdbreferes)))
    for n in range(len(pdbres)):
        pdbr, refr = pdbres[n], pdbreferes[n]
        if pdbr.resname != refr.resname:
            rsid = pdbr.resid
            if refr.resid != pdbr.resid: rsid = "%s(%s)" % (pdbr.resid, refr.resid)
            raise ValueError(
                "PDB and reference are different at resid %s: %s vs %s" % (rsid, pdbr.resname, refr.resname))
        pdbr.nter = refr.nter
        pdbr.cter = refr.cter
        pdbr.topology = refr.topology




def main( ):
    if has_argparse:
        parser = argparse.ArgumentParser(description=__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument("pdb", help="PDB file to reduce")
        parser.add_argument("output", help="all-atom reduced output PDB file", nargs="?")
    else:
        parser = optparse.OptionParser()
        parser.add_argument = parser.add_option

    parser.add_argument("--heavy", help="Ignore all hydrogens", action="store_true")
    parser.add_argument("--refe", "--reference",
                        help="Analyze the hydrogens of a reference file to determine histidine/cysteine states")
    parser.add_argument("--autorefe", help="Analyze the hydrogens of the input PDB to determine histidine/cysteine states",
                        action="store_true")
    parser.add_argument("--dna", help="Automatically interpret nucleic acids as DNA", action="store_true")
    parser.add_argument("--rna", help="Automatically interpret nucleic acids as RNA", action="store_true")
    parser.add_argument("--rnalib",
                        help="Use the ATTRACT RNA mononucleotide library to build missing atoms for nucleotides",
                        action="store_true")
    parser.add_argument("--pdb2pqr",
                        help="Use PDB2PQR to complete missing atoms. If no reference has been specified, analyze the hydrogens to determine histidine/cysteine states",
                        action="store_true")
    parser.add_argument("--whatif",
                        help="Use the WHATIF server to complete missing atoms. If no reference has been specified, analyze the hydrogens to determine histidine/cysteine states",
                        action="store_true")
    parser.add_argument("--termini",
                        help="An N-terminus and a C-terminus (5-terminus and 3-terminus for nucleic acids) will be added for each chain",
                        action="store_true")
    parser.add_argument("--nter", "--nterm", dest="nter",
                        help="Add an N-terminus (5-terminus for nucleic acids) for the specified residue number",
                        action="append",
                        type=int, default=[])
    parser.add_argument("--cter", "--cterm", dest="cter",
                        help="Add a C-terminus (3-terminus for nucleic acids) for the specified residue number",
                        action="append",
                        type=int, default=[])
    parser.add_argument("--manual", help="""Enables manual mode. 
    In automatic mode (default), aareduce tries to produce a PDB file that can be used directly by ATTRACT. In case of missing atoms, a number of 
    last-resort fixes are attempted that add pseudo-hydrogens at the position of its connected heavy atom. If there are other missing atoms,
    an exception is raised.
    In manual mode, last-resort fixes are disabled, and missing atoms are simply printed as XXXXXXX in their coordinates. These coordinates
    cannot be read by ATTRACT, they need to be edited manually by the user.
    """, action="store_true")
    parser.add_argument("--trans", "--transfile", dest="transfile",
                        help="Additional trans file that contains additional user-defined atom types (e.g. modified amino acids)",
                        action="append", default=[])
    parser.add_argument("--top", "--topfile", dest="topfile",
                        help="Additional topology file in CNS format that contains additional user-defined atom types (e.g. modified amino acids)",
                        action="append", default=[])
    parser.add_argument("--patch", dest="patches",
                        help="Provide residue number and patch name to apply", nargs=2, action="append", default=[])
    parser.add_argument("--chain", help="Set the chain in the output PDB", default=" ")
    parser.add_argument("--startres", help="Set residue number of first residue", type=int, default=1)
    parser.add_argument("--startatom", help="Set atom number of first atom", type=int, default=1)
    parser.add_argument("--mutate", dest="mutatefiles",
                        help="Provide a 2-column residue mutation file", action="append", default=[])
    parser.add_argument("--modres",
                        help="Interpret HETATM records as ATOM if they have a protein backbone", action="store_true")
    parser.add_argument("--modbase",
                        help="Interpret HETATM records as ATOM if they have at least three sugar atoms",
                        action="store_true")
    parser.add_argument("--batch",
                        help="run aareduce in batch mode. Input and output must be two lists of PDBs", action="store_true")
    parser.add_argument("--dumppatch",
                        help="Dump all applied patches to a file", action="store_true")
    parser.add_argument("--readpatch",
                        help="Read previously applied patches from a file (requires converted input pdb)",
                        action="store_true")
    if has_argparse:
        args = parser.parse_args()
    else:
        args, positional_args = parser.parse_args()
        args.pdb = None
        args.output = None
        if positional_args:
            args.pdb = positional_args[0]
            if len(positional_args) > 1: args.output = positional_args[1]
    output_name_pdb = []
    output_name_mapping = []
    for outfile in output:
        head, tail = os.path.split( outfile)
        output_name_pdb.append(tail)
        output_name_mapping.append(  os.path.splitext(tail)[0]+"-mapping" )
        output_path= head


    aareduce(args.pdb, output_path, output_name_pdb, output_name_mapping, args.chain, args.readpatch, args.patches, args.topfile,
             args.transfile, args.mutatefiles, args.rnalib, args.rna, args.dna,
             args.heavy, args.autorefe, args.refe, args.batch, args.termini,
             args.modbase, args.modres, args.nter, args.cter, args.pdb2pqr, args.whatif,
             args.manual, args.dumppatch, args.startres, args.startatom)

def aareduce( args_pdb ,  output_path, output_name_pdb,output_name_mapping ,args_chain= " ", args_readpatch= False, args_patches=[], args_topfile=[], args_transfile=[], args_mutatefiles=[], args_rnalib = False,args_rna = False, args_dna = False,
              args_heavy = False, args_autorefe = False, args_refe = False,  args_batch = False,args_termini= False,
              args_modbase = False, args_modres = False,args_nter = [], args_cter = [], args_pdb2pqr = False,args_whatif = False, args_manual = False, args_dumppatch = False, args_startres = 1, args_startatom = 1, folder_parameter = os.environ['ATTRACTDIR']+'/../allatom'):
    #currdir = os.path.abspath(os.path.split(__file__)[0])

    topstream = [(folder_parameter + "/topallhdg5.3.pro", "oplsx"),
                 (folder_parameter + "/dna-rna-allatom.top", "dna-rna")
                 ]
    transfiles = [(folder_parameter + "/oplsx.trans", "oplsx"),
                  (folder_parameter + "/dna-rna.trans", "dna-rna")]
    code_to_type={}
    mutations={}
    assert len(args_chain) == 1, args_chain
    
    if args_rna and args_dna:
        raise ValueError("--dna and --rna are mutually incompatible")
    
    if args_rnalib and not args_rna:
        raise ValueError("--rnalib requires option --rna")
    
    if args_heavy and (args_autorefe or args_refe):
        raise ValueError("--(auto)refe and --heavy are mutually incompatible")
    
    if args_autorefe and args_refe:
        raise ValueError("--autorefe and --refe are mutually incompatible")
    if args_autorefe:
        args_refe = args_pdb
    
    if args_batch and (output_name is None and output_path is None):
        raise ValueError("--batch requires a file list as output argument")
    
    if args_readpatch and len(args_patches):
        raise ValueError("--readpatch and explicit patch specification via --patch are mutually incompatible")
    
    for fnr, f in enumerate(args_topfile):
        assert os.path.exists(f), f
        topstream.append((f, "userfile-%d" % (fnr + 1)))
    for f in args_transfile:
        transfiles.append((f, "userfile-%d" % (fnr + 1)))
    
    for f, name in topstream:
        parse_cns_top.parse_stream(open(f), name)


    for f, name in transfiles:
        parse_transfile(f, name, code_to_type)

    for f in args_mutatefiles:
        assert os.path.exists(f), f
        for l in open(f):
            h = l.find("#")
            if h != -1: l = l[:h]
            ll = l.split()
            if len(ll) == 0: continue
            assert len(ll) == 2, l
            mutations[ll[0]] = ll[1]
    rnalib=None
    if args_rnalib:
        rnalib = pdbcomplete.load_rnalib()
    
    if args_batch:
        infiles = read_filelist(args_pdb)
        for f in infiles:
            assert os.path.exists(f), f
        outfiles = read_filelist(output_name)
        for pdb, outfile in zip(infiles, outfiles):
            pdblines, mapping, pdbtop = run(pdb, rnalib,code_to_type,mutations, args_termini, args_modbase, args_modres , args_chain , args_nter, args_cter , args_readpatch, args_patches, args_rnalib , args_heavy ,
                                            args_pdb2pqr , args_refe, args_whatif, args_manual , args_startatom,args_startres    )
            outfile = output_path + outfile
            outf = open(outfile, "w")
            for l in pdblines:
                print >> outf, l
            outf.close()
            if args_dumppatch:
                outfilep = os.path.splitext(outfile)[0] + ".patch"
                outf = open(outfilep, "w")
                for i, res in enumerate(pdbtop):
                    if len(res.topology.patches):
                        for p in res.topology.patches:
                            outf.write(str(i + args_startres) + ' ' + p + '\n')
                outf.close()
    else:
        outfile = os.path.splitext(args_pdb)[0] + "-aa.pdb"
        if output_name_pdb is not None:
            outfile = output_path + output_name_pdb
        pdblines, mapping, pdbtop = run(args_pdb, rnalib,code_to_type,mutations, args_termini, args_modbase, args_modres , args_chain , args_nter, args_cter , args_readpatch, args_patches, args_rnalib , args_heavy ,
                                            args_pdb2pqr , args_refe, args_whatif , args_manual , args_startatom,args_startres  )
        outf = open(outfile, "w")
        for l in pdblines:
            print >> outf, l
        #for v1, v2 in mapping:
         #   print
          #  v1, v2
        outf.close()

        #outf = open(os.path.splitext(outfile)[0] + ".mapping", "w")
        outf = open(output_path + output_name_mapping, "w")
        for v1, v2 in mapping:
            print >> outf,v1, v2
        outf.close()
        if args_dumppatch:
            outfilep = os.path.splitext(outfile)[0] + ".patch"
            outf = open(outfilep, "w")
            for i, res in enumerate(pdbtop):
                if len(res.topology.patches):
                    for p in res.topology.patches:
                        outf.write(str(i + args_startres) + ' ' + p + '\n')
            outf.close()

def run(pdbfile, rnalib,code_to_type,mutations, args_termini= False, args_modbase = False, args_modres = False,
        args_chain = " ", args_nter = [], args_cter = [], args_readpatch = False, args_patches=[],
        args_rnalib = False, args_heavy = False, args_pdb2pqr = False, args_refe= None, args_whatif = False, args_manual = False, args_startatom = 1,args_startres = 1 ):

    
    pdb = read_pdb(open(pdbfile),mutations=mutations, add_termini=args_termini, modbase=args_modbase, modres=args_modres)
    pdblines = write_pdb(pdb, args_chain, code_to_type, args_startatom =args_startatom, args_startres = args_startres)[0]

    termini_pdb(pdb, args_nter, args_cter)
    patches = {}
    if args_readpatch:
        indata = open(os.path.splitext(pdbfile)[0] + '.patch').readlines()
        indata = [line.split() for line in indata]
        args_patches = indata

    for p in args_patches:
        resid = p[0].strip()
        resindices = [ri for ri, r in enumerate(pdb) if r.resid.strip() == resid]
        if len(resindices) == 0:
            raise ValueError("No residues have resid %s" % resid)
        elif len(resindices) > 1:
            raise ValueError("Multiple residues have resid %s" % resid)
        resid2 = pdb[resindices[0]].resid
        if resid2 not in patches: patches[resid2] = []
        pname = p[1].lower()
        if pname == "none": pname = None
        patches[resid2].append(pname)
    patch_pdb(pdb, patches)

    if args_refe:
        refe = read_pdb(open(args_refe),mutations=mutations, add_termini=args_termini)
        patch_pdb(refe, patches)
        if not args_heavy:
            update_patches(refe)
        set_reference(pdb, refe)
    if args_rnalib:
        pdbcomplete.apply_rnalib(pdb, rnalib, args_heavy)
    if args_pdb2pqr:
        pdblines = write_pdb(pdb, args_chain, code_to_type, one_letter_na=True, args_startatom =args_startatom, args_startres = args_startres)[0]
        pqrlines = run_pdb2pqr(pdblines)
        pqr = read_pdb(pqrlines,mutations=mutations)
        pdbcomplete.pdbcomplete(pdb, pqr)
        if not args_heavy and not args_refe:
            update_patches(pdb)
    if args_whatif:
        pdblines = write_pdb(pdb, args_chain, code_to_type, one_letter_na=True, args_startatom =args_startatom, args_startres = args_startres)[0]
        whatiflines = run_whatif(pdblines)
        whatif = read_pdb(whatiflines,mutations=mutations)
        pdbcomplete.pdbcomplete(pdb, whatif)
        if not args_heavy and not args_refe and not args_pdb2pqr:
            update_patches(pdb)

    if args_refe:
        pdbfix(pdb, refe)

    if not args_manual:
        pdb_lastresort(pdb)
        check_pdb(pdb, heavy=args_heavy)
    pdblines, mapping = write_pdb(pdb, args_chain,code_to_type, heavy=args_heavy, args_startatom =args_startatom, args_startres = args_startres)
    return pdblines, mapping, pdb


#aareduce(  args_pdb="/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/receptor.pdb" , output_path="/home/glenn/Documents/Masterthesis/testfolder/benchmark_test/input/", output_name_pdb="receptor-aa.pdb", output_name_mapping="receptor.mapping",
#          args_chain="A", args_readpatch= False, args_patches=[], args_topfile=[], args_transfile=[], args_mutatefiles=[], args_rnalib = False,args_rna = False, args_dna = False,
 #         args_heavy = False, args_autorefe = False, args_refe = False,  args_batch = False ,args_termini= False, args_modbase = False, args_modres = False,
 #         args_nter = [], args_cter = [],args_pdb2pqr = True, args_whatif = False, args_manual = False, args_dumppatch = True, args_startres = 1,args_startatom = 1)
