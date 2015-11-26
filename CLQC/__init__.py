"""
Copyright (c) <2015>,  Sven Giese (sven.giese@tu-berlin.de)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software
   must display the following acknowledgement:
   This product includes software developed by the <organization>.
4. Neither the name of the <organization> nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import pandas as pd
from Bio import PDB
import HTSeq
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator
import brewer2mpl
import matplotlib as mpl
import re

try:
    from Bioinf2 import mpl_overwrite
except:
    print "ERROR! Bioinf2 package could not be found. Rolling back \
    to default plotting settings"


class QC_report():
    """
    Class to store quality control reports on cross-linking MS data.

    Initate the report with the following parameters:

    Parameters:
    ----------------------------------------
    PDB_loc: str,
             file location of the PDB
    links_loc: str,
             file location of the links csv from xiFDR
    fasta_loc: str,
              file location of the fasta file
    """

    def __init__(self, PDB_loc, links_loc, fasta_loc):
        """

        Initate the report with the following parameters:

        Parameters:
        ----------------------------------------
        PDB_loc: str,
                 file location of the PDB
        links_loc: str,
                 file location of the links csv from xiFDR
        fasta_loc: str,
                  file location of the fasta file
        """

        self.PDB_loc = PDB_loc
        self.links_loc = links_loc
        self.fasta_loc = fasta_loc
        self.report = ""
        self.results = None

    def read_links(self, onlyTT=True):
        """
        Reads the link file and stores the information in a dataframe.

        Parameters:
        --------------------------------
        onlyTT: bool,
                  filter Decoy hits automatically.
        """
        #read Xi results
        results = pd.read_csv(self.links_loc)
        #filter decoys and reduce DF
        results = results[results["isTT"] == True]
        return(results)

    def Create_CoverageProfile(self, outdir, winsize=10, norm=True,
                               residue1="K", residue2="STY", annotateK=True,
                               proteinname=""):
        """
        Create a coverageprofile.

        Parameters:
        ----------------------------------------------

        outdir: str,
                 str of residues that are linkable. default: "KSTY"
        winsize: int,
                 size of the windows (to count cross-linked residues)
        norm: bool,
              normalize counts by sum
        residue1: str,
                  to distinguish multiple residue profiles in the protein
                  (e.g."K")
        residue2 str,
                 second part of the residue profile (e.g. "STY")

        annotateK: bool,
                  write all Lysines on the cricular graph (True). False,
                  will annotate every 10th residue
        """
        #read Xi results
        results = self.read_links()

        #transofrm residue strings to lists
        residues1 = list(residue1)
        residues2 = list(residue2)
        #set windowsize
        fasta_seq = self.__get_fasta_obj__().seq
        proteinlength = len(fasta_seq)
        bins = np.arange(0, proteinlength+winsize, winsize)

        #get linked residues
        results["residue1"] = [i[20] for i in results["LinkWindow1"]]
        results["residue2"] = [i[20] for i in results["LinkWindow2"]]

        #all crosslink sites
        all_sites = np.append(results["fromSite"], results["ToSite"])
        x3 = np.bincount(np.digitize(all_sites, bins), minlength=len(bins))
        #plot the profile of KTSY
        # find everey position of the residues in the protein sequence
        starts_r1 = [match.start() for match in
                     re.finditer("|".join(residues1), fasta_seq)]
        starts_r2 = [match.start() for match in
                     re.finditer("|".join(residues2), fasta_seq)]
        #starts = [match.start() for match in re.finditer("K", fasta_seq.seq)]
        x4 = np.bincount(np.digitize(starts_r1, bins), minlength=len(bins))
        x5 = np.bincount(np.digitize(starts_r2, bins), minlength=len(bins))

        #normalize counts
        if norm:
            x3 = 1. * x3 / np.sum(x3)
            x4 = (1. * x4) / np.sum(x4)
            x5 = (1. * x5) / np.sum(x5)
            normed = "(normed)"
        else:
            normed = "(counts)"

        #======================================================================
        #    Actual plot
        #======================================================================
        f, ax = plt.subplots(1, figsize=(14.69, 8.27))
        ax.plot(bins, x3, '-', label="measured links", alpha=0.8)
        ax.plot(bins, x4, '--', label=residue1 + " occs. protein", alpha=0.8)
        ax.plot(bins, x5, '--', label=residue2 + " occs. protein", alpha=0.8)
        ax.legend(loc="upper left")
        ax.set_xlabel("CL pos {}\n(bin size {})".format(proteinname, winsize))
        ax.set_ylabel("frequency (%)".format(normed))
        ax.set_title("cross-link site coverage")
        sns.despine()
        f.savefig("{}coverage_{}.svg".format(outdir, normed),
                  bbox_inches='tight', pad_inches=0.1)
        f.savefig("{}coverage_{}.png".format(outdir, normed),
                  bbox_inches='tight', pad_inches=0.1)
        plt.clf()

        #======================================================================
        # Add a circular plot
        #======================================================================
        #get the graph going
        import networkx as nx
        sites_nx = zip(results["fromSite"], results["ToSite"])
        nnodes = np.arange(0, len(fasta_seq), 1)
        node_labels = {j: "{}\n({})".format(i, j) if (j % 20 == 0) else "" for i, j
                       in zip(list(fasta_seq), nnodes)}

        node_labels = {j: "{}\n({})".format(i, j) if (i in "K") else "" for i, j
                       in zip(list(fasta_seq), nnodes)}
        #create the graph
        G = nx.DiGraph()
        G.add_nodes_from(nnodes)
        pos = nx.circular_layout(G)
        # add the edges and nodes
        nx.draw(G, pos, with_labels=False, linewidths=0.0, node_color="r")
        nx.draw_networkx_labels(G, pos, node_labels, font_color="k")
        # each line adds another round of edges (color)

        nx.draw_networkx_edges(G, pos, edgelist=sites_nx, width=1,
                               edge_color="black", arrows=False)
        plt.savefig("{}all_graph.png".format(outdir),
                    bbox_inches='tight', pad_inches=0.1)
        plt.savefig("{}all_graph.svg".format(outdir),
                    bbox_inches='tight', pad_inches=0.1)


    def Create_ContactMap(self, outfile, maxdistance=25, symmetric=True,
                          proteinname="", title="", annotate=False, chain="A"):
        """
        Generates a contactmap from cross-link data.

        Parameters:
        -----------------------------------
        residues_arr: array,
                      pymol residues
        max_distance: float,
                      distance in Angstrom that is theoretical
                      possible to bridge with the cross-linker
        """
        #open the PDB structure
        p = PDB.PDBParser()
        structure = p.get_structure("HSA", self.PDB_loc)[0][chain]
        residues_arr = self.__get_residues__(structure)

        #create a sequence position, residue (atom) mapping
        res_dic = {i.get_id()[1]: i for i in residues_arr}

        #get FASTA reference sequence
        protein_seq = self.__get_fasta_obj__()

        #place holder for link position in protein1 (x) and protein2 (y)
        xposarr = []
        yposarr = []
        dist = []
        for i, j in enumerate(residues_arr):
            for k, l in enumerate(residues_arr):
                #if distance smaller than constraint, put into contactmap
                if (np.abs(j["CA"] - l["CA"]) <= maxdistance):
                    xposarr.append(j.get_id()[1])
                    yposarr.append(l.get_id()[1])
                    dist.append(np.abs(j["CA"] - l["CA"]))

                else:
                    pass

        #transform tu numpy array
        xposarr = np.array(xposarr)
        yposarr = np.array(yposarr)
        dist = np.array(dist)

        #read Xi results
        results = self.read_links()

        cross_linkdist = self.__get_crosslink_distances__(res_dic,
                                                          results)
        results["distance"] = cross_linkdist
        fps = results[results["distance"] >= maxdistance].copy()
        tps = results[results["distance"] <= maxdistance].copy()
        #======================================================================
        #        Plot the contactmap
        #======================================================================
        n = 5
        steps = np.arange(0, maxdistance+n, maxdistance/n)
        assign_bins = np.digitize(dist, steps)
        blues = brewer2mpl.get_map('Blues', 'Sequential', 9).mpl_colors[::-1]
        bluesmap = brewer2mpl.get_map('Blues', 'Sequential', 9).mpl_colormap

        colors = [blues[i] for i in assign_bins]

        norm = mpl.colors.Normalize(0, 25, 10)
        f, ax = plt.subplots(1, 2, figsize=(11.69, 8.27))
        gs = mpl.gridspec.GridSpec(1, 2, width_ratios=[12, 1])
        ax0 = plt.subplot(gs[0])
        ax0.scatter(xposarr, yposarr, c=colors, alpha=0.8, s=70,
                    edgecolor='none')
        ax0.set_xlim(-10, np.max(tps.fromSite))
        ax0.set_ylim(-10, np.max(tps.ToSite))
        ax0.set_title(title)
        plt.scatter(tps.fromSite, tps.ToSite, s=90, alpha=0.7,
                    linewidth=2, edgecolor='black', color="g")
        plt.scatter(tps.ToSite, tps.fromSite, s=90, alpha=0.7,
                    linewidth=2, edgecolor='black', color="g")
        plt.scatter(fps.fromSite, fps.ToSite, s=90, alpha=0.7,
                    linewidth=2, edgecolor='black', color="r")
        plt.scatter(fps.ToSite, fps.fromSite, s=90, alpha=0.7,
                    linewidth=2, edgecolor='black', color="r")
        if annotate:
            for label, xi, yi in zip(annotate, fps.fromSite.values,
                                     fps.ToSite.values):
                plt.annotate(label, (xi, yi))
        sns.despine()
        ax1 = plt.subplot(gs[1])
        ax0.set(xlabel="Residue Position 1 {}".format(proteinname),
                ylabel="Residue Position 2 {}".format(proteinname))
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=bluesmap,
                                        norm=norm,
                                        orientation='vertical',
                                        ticks=steps[::-1])
        cb1.set_ticklabels(steps[::-1])
        cb1.set_ticks(steps)
        cb1.set_label("distance in AA")
        ax1.yaxis.set_ticks_position('right')
        plt.subplots_adjust(wspace=0.05)
        f.tight_layout()
        f.savefig("{}_contactmap_new_{}A.pdf".format(outfile, maxdistance),
                  bbox_inches='tight', pad_inches=0.1)
        f.savefig("{}_contactmap_new_{}A.png".format(outfile, maxdistance),
                  bbox_inches='tight', pad_inches=0.1)

    def Create_Distogram(self, chain, overlength=25.0, addRandom=True,
                         output="", proteinname="", names=""):
        """
        Plots a link distogram of the cross-link identifications


        A textual report (some stats) and the figures are appended to the QC-
        report.

        Parameters:
        -------------------------------------------
        chain: str,
               Which chain to extract from the PDB, e.g. "A" for HSA
        overlength: float,
                     distance that is considereded not long distance for
                     cross-links
        addRandom: bool,
                   boolean if a random distribution from all links should be
                   added.
        output: str,
                destination of the figure to be plotted (coml. path and name)

        names: list,
               list with short identifiers witht he same length as the #inputs

        """
        self.report += """\r\nQuality report for Distogram:
----------------------------------------------------------------------\r\n"""
        #set default name
        if names == "":
            names = [""]

        #open the PDB structure
        p = PDB.PDBParser()
        structure = p.get_structure("HSA", self.PDB_loc)[0][chain]
        if addRandom:
            residues_arr = self.__get_residues__(structure)
            distances = self.__get_quadratic_distances__(residues_arr)
        else:
            distances = []

        #create a sequence position, residue (atom) mapping
        res_dic = {i.get_id()[1]: i for i in residues_arr}

        #get FASTA reference sequence
        protein_seq = self.__get_fasta_obj__()

        #read Xi results
        results = self.read_links()

        cross_linkdist = self.__get_crosslink_distances__(res_dic,
                                                          results)
        results["distance"] = cross_linkdist
        cross_linkdist = cross_linkdist[np.where(cross_linkdist >= 0)]

        #store new data with distance information
        results.to_csv(output+"Links_distance.csv")
        #plot the figure
        self.__plot_link_histogram__(distances, cross_linkdist,
                                     overlength, title="", outfile=output,
                                     proteinname=proteinname)


    def __get_residues__(self, structure):
        """
        Gets all amino acids residues from a given structure and stores them
        in an array.

        parameters:
        ----------------
        structure: PDB strutore obj,
                   openened PDB structure file object

        Returns:
        ---------------------------------------
        array: np-arr,
               residue objects from Bio.PDB
        """
        residues_arr = []
        for res_i in structure.get_residues():
            if PDB.is_aa(res_i):
                residues_arr.append(res_i)
        return(np.array(residues_arr))

    def __get_fasta_obj__(self):
        """
        Use HTSeq to retrieve the fasta sequence

        Returns:
        ---------------------------
        sequence: str,
                  fsta sequence
        """
        fasta = HTSeq.FastaReader(self.fasta_loc)
        k = 0
        for seq in fasta:
            k += 1
            pass

        #check fasta file
        if k > 1:
            print "Attention: More than 1 sequence in your fasta file! Last \
last sequence is used for the analysis!"
        return(seq)

    def __get_quadratic_distances__(self, residues_arr):
        """
        Computes all pairwise distances from the Carbon atoms in the structure.
        The resulting distribution is used as "random" background
        distribution.

        Parameters:
        -----------------
        residues_arr: arr_like
                  residuesof the structure

        Returns:
        ---------------------------
        dist: np-arr,
              flaot values of distances (once for each residue combination)

        AA_list: np-arr,
                flaot values of distances full quadratic list of all
                combinations
        """
        #create the matrix for storing the distances
        dist = []
        for i, j in enumerate(residues_arr):
            for k, l in enumerate(residues_arr):
                if i >= k:
                    pass
                else:
                    #only use each computed sitance once
                    dist.append(np.abs(j["CA"] - l["CA"]))
        return (dist)

    def __get_crosslink_distances__(self, res_dic, results, bfactors=False):
        """
        Get the pairwise distances from the cross-linked peptides.

        Parameters:
        -----------------------
        res_dic: dictionary,
                 sequence position, residue (atom) mapping

        results: dataframe,
                 XiLink FDR results table. Containg "fromSite" and "ToSite"
                 column

        """
        #compute distances cross-link data
        diff_crosslinks = []
        error_count = 0
        all_counter = 0
        tmp_report = []
        for from_site, to_site in zip(results["fromSite"], results["ToSite"]):
            all_counter += 1
            try:
                #get data
                diff_crosslinks.append(res_dic[from_site]["CA"] -
                                       res_dic[to_site]["CA"])
            except:
                tmp_report.append("error linkapair not found: \
                {} and {}".format(from_site, to_site))
                error_count += 1
                diff_crosslinks.append(-1)

        tmp_report.append("{}/{} ({}%) cross-links could't be mapped to the \
crystal structure\r\n".format(error_count, all_counter,
                          np.round(error_count*1./all_counter, 2)))
        self.report += "\r\n".join(tmp_report)
        return (np.array(diff_crosslinks))

    def __plot_link_histogram__(self, distances, cross_linkdist, overlength,
                                title, outfile, proteinname):
        """
        Plots a Link histogram.

        Paramters:
        -------------------------------------------------------
        distances: np.array,
                   floats distancesof the random distribution
        cross_linkdist: np-arr,
                         floats of measured cross-link distances
        overlength: float,
                    overlength boundary for cross-links considered
                    long-distance
        title: str,
               title of the figure
        outfile: str,
                 outfile location
        proteinname: str,
                     name of the protein
        """

        #weights for true frequency plot
    #    w1 = [1. * np.ones_like(valuesi) / len(distances) for valuesi in distances]
    #    w2 = [1. * np.ones_like(valuesi) / len(cross_linkdist) for valuesi in cross_linkdist]
        noverlength = len(np.where(cross_linkdist >= overlength)[0])
        nlinks = 1. * len(cross_linkdist)
        f, ax = plt.subplots(1, figsize=(11.69, 8.27))
        ax.hist(distances, bins=40, normed=True, alpha=0.9)
        ax.hist(cross_linkdist, bins=40, normed=True, alpha=0.9)
        #ax.hist(distances, weights=w1, bins=40, alpha=0.9)
        #ax.hist(cross_linkdist, weights=w2, bins=40, alpha=0.9)
        ax.axvline(overlength, lw=2, color="k")
        ax.set_xlabel("Angstrom")
        ax.set_ylabel("relative frequency")
        ax.annotate("Number of links: {}\nlong distance (>={}A): {}\nFDR: {:.2f}%".format
                    (int(nlinks), overlength, noverlength,
                     (noverlength/nlinks)*100), xy=(40, 0.06))
        ax.set_title(title + " " + proteinname)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        sns.despine()
        f.tight_layout()
        #estimate fdr with groundtruth data
        maxdist = np.max(cross_linkdist)
        fdr = np.round(noverlength * 1.0 / len(cross_linkdist) * 100, 2)
        #estimate from dsitrubtion right to cl site
        fdr_est = np.round((1. * noverlength / (maxdist - overlength)) *
                            overlength, 2) + noverlength

        self.report += """\r\n    FDR overview:
    ---------------------------------------------------------------------------
    all cross-links: {}
    FDR (hard-cutoff {}A): #{} ({}%)
    true FDR (estimated): #{} ({:.2f}%)\
    """.format(len(cross_linkdist), overlength, noverlength, fdr, fdr_est,
               np.round(fdr_est * 1.0 / len(cross_linkdist) * 100, 2))

        f.savefig("{}_linkhistogram_{}A.pdf".format(outfile, overlength),
                  bbox_inches='tight', pad_inches=0.1)
        f.savefig("{}_linkhistogram_{}A.png".format(outfile, overlength),
                  bbox_inches='tight', pad_inches=0.1)
        return(f)


def test():
    import glob
    PDB_file = "/home/sven/data/Qutsvs.Velos/1AO6_chainA.pdb"
    FASTA_file = "/home/sven/data/Qutsvs.Velos/1AO6.fasta"
    Links = sorted(glob.glob("/home/sven/data/LumosFragmentation2015/*pLink*Links_xiFDR*"))
    PSMs = sorted(glob.glob("/home/sven/data/LumosFragmentation2015/*pLink*PSM_xiFDR*"))

    Links = [i for i in Links if "_distance." not in i]
    PSMs = [i for i in PSMs if "_distance." not in i]

    outfile = "/home/sven/data/LumosFragmentation2015/Distances/"
    outfile = "/home/sven/Dropbox/Fragmentation_Manuscript_privat/Figures_single/"
    Names = [i.split("/")[-1].split("_1.")[0] for i in Links]
    proteinname = "HSA"


    pdb_loc = "/home/sven/data/Qutsvs.Velos/1AO6_chainA.pdb"
    links_loc = Links[0]
    psms_loc = PSMs[0]
    fasta_loc = "/home/sven/data/Qutsvs.Velos/1AO6.fasta"

    QC_rep = QC_report(pdb_loc, links_loc, fasta_loc)
    QC_rep.Create_Distogram("A", overlength=25.0, proteinname="HSA", names="Exp1")
    print QC_rep.report
    QC_rep.Create_ContactMap(outfile="", maxdistance=25, symmetric=True,
                             annotate=None, proteinname="", title="HSA", chain="A")

