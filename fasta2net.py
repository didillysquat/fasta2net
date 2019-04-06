#!#!/usr/bin/env python3
"""A python wrapper for creating splits tree median joining networks.
Inputs will be a fasta file and a name file..
It will be good to add colour functionality and gap functionality to it as well."""

# eventually the fasta and the names file will come from command line arguments but for the time being
# we will hardcode their paths in.
# lets start by running the SP network.

from plumbum import local
import os
import subprocess
import networkx as nx
from collections import defaultdict
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import argparse



class Fasta2Net:
    def __init__(
            self, fasta_file_path ,output_dir=None,  names_file_path=None, color_map_file_path=None,
            spring_pos_iterations=1000, spring_pos_k_value=0.1, labels=True, dpi=1200):
        self.cwd = os.path.dirname(os.path.realpath(__file__))
        if output_dir is None:
            self.output_dir = self.cwd
        else:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            self.output_dir = output_dir
        self.fasta_file_path = fasta_file_path

        self.names_file_path = names_file_path
        self.fasta_dict = self._create_fasta_dict()
        self.names_dict = self._create_names_dict()
        self.color_dict = self._create_color_dict(color_map_file_path)

        # Command line variables
        self.spring_pos_iterations = spring_pos_iterations
        self.spring_pos_k_value = spring_pos_k_value
        self.labels = labels
        self.dpi = dpi

        # attributes for aligning the fasta
        self.aligned_fasta_interleaved_path = self._set_aligned_fasta_interleaved_path()
        self.aligned_fasta_interleaved_as_list = []
        self.aligned_fasta_uncropped = None
        self.aligned_fasta = None
        self.mafft_plum_bum_exe = local["mafft-linsi"]

        # nexus file
        self.nexus_file = None
        self.nexus_file_output_path = os.path.join(self.output_dir, 'splitstree_in.nex')

        # splits tree files
        self.splits_out_path = os.path.join(self.output_dir, 'splitstree_out.nex')
        self.ctrl_path = os.path.join(self.output_dir, 'splitstree_ctrl')
        self.splits_out_file = None
        # we can work with a minimised version of the splits out file that is only the network block
        self.splits_net_block = None

        # Plotting
        self.nex_graph = nx.Graph()
        self.vertice_id_to_seq_id_dict = defaultdict(list)
        self.vertice_list = None
        self.vertice_id_to_abund_dict = {}
        self.vertice_sizes = None
        self.vertice_colors = None
        self.edges_list = []
        self.f, self.ax = plt.subplots(1, 1, figsize=(10,10))
        self.spring_pos = None

    def _set_aligned_fasta_interleaved_path(self):
        fasta_extension = os.path.splitext(self.fasta_file_path)
        if fasta_extension not in ['.fasta', '.fas', '.fa']:
            raise RuntimeError('Unrecognised fasta format')
        return os.path.join(
            self.output_dir,
            self.fasta_file_path.split('/')[-1].replace(fasta_extension, '_aligned{}'.format(fasta_extension)))


    def make_network(self):

        self.aligned_fasta = self._create_aligned_fasta()

        # we have to make the new nexus format by hand as the biopython version was putting out old versions.
        self._make_and_write_splits_tree_nexus()

        self._make_and_write_cntrl_file()

        self.splits_net_block = self._run_splits_trees_and_extract_net_block()

        self._get_vertices()

        self._make_vertice_id_to_abund_dict()

        self._infer_vertice_sizes()

        self._infer_vertice_colors()

        self._make_edge_list()

        self._convert_vert_and_ege_names_to_seq_id_names()

        self._add_vertices_to_graph_object()

        self._add_edges_to_graph_object()

        self._generate_spring_pos()

        self._set_ax_limits()

        self._draw_vertices()

        self._draw_edges()

        self._format_ax()

        self._output_network()

    def _output_network(self):
        plt.savefig(os.path.join(self.output_dir, 'network.png'), dpi=self.dpi)
        plt.savefig(os.path.join(self.output_dir, 'network.svg'))

    def _format_ax(self):
        """The collections that were added from the nx.draw_networkx_nodes and nx.draw_networkx_edges
        called in self._draw_vertices and self._draw_eges can be accessed
        to format the nodes and the edges."""
        # https://stackoverflow.com/questions/22716161/how-can-one-modify-the-outline-color-of-a-node-in-networkx
        # https://matplotlib.org/api/collections_api.html
        self.ax.collections[0].set_edgecolor('grey')
        self.ax.collections[0].set_linewidth(0.5)
        self.ax.collections[1].set_linewidth(0.5)
        self.ax.set_axis_off()

    def _draw_vertices(self):
        nx.draw_networkx_nodes(
            self.nex_graph, pos=self.spring_pos, node_size=self.vertice_sizes, with_labels=self.labels, alpha=1,
            edgecolors='black', ax=self.ax)
    def _draw_edges(self):
        nx.draw_networkx_edges(G=self.nex_graph, pos=self.spring_pos, ax=self.ax)

    def _generate_spring_pos(self):
        self.spring_pos = nx.spring_layout(
            self.nex_graph, k=self.spring_pos_k_value, iterations=self.spring_pos_iterations)

    def _set_ax_limits(self):
        """We will need to set the limits dynamically as we don't know what the positions are going to be.
    I think we will want to have the x and y limits the same so that we end up with a square
    we will therefore just look for the bigest and smallest values, add a buffer and then set the
    axes limits to these
    """
        max_ax_val = 0
        min_ax_val = 9999
        for vert_pos_array in self.spring_pos.values():
            for ind in vert_pos_array:
                if ind > max_ax_val:
                    max_ax_val = ind
                if ind < min_ax_val:
                    min_ax_val = ind
        buffer = 0.2
        self.ax.set_ylim(min_ax_val - buffer, max_ax_val + buffer)
        self.ax.set_xlim(min_ax_val - buffer, max_ax_val + buffer)

    def _add_vertices_to_graph_object(self):
        g.add_nodes_from(self.vertice_list)

    def _add_edges_to_graph_object(self):
        g.add_edges_from(self.edges_list)

    def _convert_vert_and_ege_names_to_seq_id_names(self):
        """Convert the vertices back to the original sequences names so that we can label them on the networks
        do the same for the edges_list"""
        self.vertice_list = ['-'.join(self.vertice_id_to_seq_id_dict[vert]) for vert in self.vertice_list]

        new_edge_list = []
        for tup in self.edges_list:
            new_edge_list.append(
                ('-'.join(self.vertice_id_to_seq_id_dict[tup[0]]), '-'.join(self.vertice_id_to_seq_id_dict[tup[1]])))
        self.edges_list = new_edge_list

    def _make_edge_list(self):
        """The edge list should be a list of tuples
        NB Currently we are not taking into account the length of the edges
        """
        for i in range(len(self.splits_net_block)):
            if self.splits_net_block[i] == 'EDGES':
                # for each line in the edges section
                for j in range(i + 1, len(self.splits_net_block)):
                    if self.splits_net_block[j] == ';':
                        # then we have reached the end of the Edges section
                        break
                    items = self.splits_net_block[j].replace(',', '').split(' ')[1:]
                    self.edges_list.append((items[0], items[1]))

    def _infer_vertice_sizes(self):
        """I haven't worked out what the size units are yet so this may need tweaking accordingly"""
        self.vertice_sizes = [self.vertice_id_to_abund_dict[vert] * 3 for vert in self.vertice_list]

    def _infer_vertice_colors(self):
        self.vertice_colors = [self.color_dict[self.vertice_id_to_seq_id_dict[vert][0]] for vert in self.vertice_list]

    def _make_vertice_id_to_abund_dict(self):
        for vert_id in self.vertice_list:
            count_id_list = self.vertice_id_to_seq_id_dict[vert_id]
            if len(count_id_list) == 1:
                abund = self.names_dict[count_id_list[0]]
                self.vertice_id_to_abund_dict[vert_id] = abund
            elif len(count_id_list) > 1:
                # then we need sum the abundances for each of them
                tot = 0
                for count_id in count_id_list:
                    abund = self.names_dict[count_id]
                    tot += abund
                self.vertice_id_to_abund_dict[vert_id] = tot

    def _get_vertices(self):
        """Get a list of the nodes
        these can be got from the splits_tree_outfile
        we should use the numbered system that the splits tree is using. We can get a list of the different sequences
        that each of the taxa represent form the translate part of the network block.
        This way we can use this to run the sequences against SP to get the names of the sequencs
        we can also make the 'translation dictionary' at the same time"""
        for i in range(len(self.splits_net_block)):
            if self.splits_net_block[i] == 'TRANSLATE':
                # for each line in the translate section
                for j in range(i + 1, len(self.splits_net_block)):
                    if self.splits_net_block[j] == ';':
                        # then we have reached the end of the translation section
                        break
                    items = self.splits_net_block[j].replace('\'', '').replace(',', '').split(' ')
                    self.vertice_id_to_seq_id_dict[items[0]].extend(items[1:])

        self.vertice_list = list(self.vertice_id_to_seq_id_dict.keys())

    def _run_splits_trees_and_extract_net_block(self):
        # Run splitstree
        subprocess.run(['SplitsTree', '-g', '-c', self.ctrl_path])

        # Read in the output file
        # and then we can start making the network finally!
        with open(self.splits_out_path, 'r') as f:
            splits_tree_out_file = [line.rstrip() for line in f]

        for i in range(len(splits_tree_out_file)):
            if splits_tree_out_file[i] == 'BEGIN Network;':
                return splits_tree_out_file[i:]

    def _make_and_write_cntrl_file(self):
        """This creates a control file that can be fed to splits tree on the command line and writes it out
        """
        ctrl_file = []
        ctrl_file.append('BEGIN SplitsTree;')
        ctrl_file.append('EXECUTE FILE={};'.format(self.nexus_file_output_path))
        ctrl_file.append('SAVE FILE={} REPLACE=yes;'.format(self.splits_out_path))
        ctrl_file.append('QUIT;')
        ctrl_file.append('end;')

        # now write out the control file
        with open(self.ctrl_path, 'w') as f:
            for line in ctrl_file:
                f.write('{}\n'.format(line))

    def _make_and_write_splits_tree_nexus(self):
        self.nexus_file = self._splits_tree_nexus_from_fasta()
        # write out the new_nexus
        with open(self.nexus_file_output_path, 'w') as f:
            for line in self.nexus_file:
                f.write('{}\n'.format(line))

    def _create_aligned_fasta(self):
        (self.mafft_plum_bum_exe['--thread', -1, self.fasta_file_path] > self.aligned_fasta_interleaved_path)()
        with open(self.aligned_fasta_interleaved_path, 'r') as f:
            self.aligned_fasta_interleaved_as_list = [line.rstrip() for line in f]
        self.aligned_fasta_uncropped = self._convert_interleaved_to_sequencial_fasta(
            self.aligned_fasta_interleaved_as_list)
        return self._crop_fasta(self.aligned_fasta_uncropped)

    @staticmethod
    def _crop_fasta_df(aligned_fasta_as_pandas_df_to_crop):
        columns_to_drop = []
        for i in list(aligned_fasta_as_pandas_df_to_crop):
            # if there is a gap in the column at the beginning
            if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
                columns_to_drop.append(i)
            else:
                break
        for i in reversed(list(aligned_fasta_as_pandas_df_to_crop)):
            # if there is a gap in the column at the end
            if '-' in list(aligned_fasta_as_pandas_df_to_crop[i]) or '*' in list(aligned_fasta_as_pandas_df_to_crop[i]):
                columns_to_drop.append(i)
            else:
                break

        # get a list that is the columns indices that we want to keep
        col_to_keep = [col_index for col_index in list(aligned_fasta_as_pandas_df_to_crop) if
                       col_index not in columns_to_drop]
        # drop the gap columns
        return aligned_fasta_as_pandas_df_to_crop[col_to_keep]

    def _crop_fasta(self, aligned_fasta):
        # convert each of the sequences in the fasta into a series with the
        # series name as the sequence name from the fasta
        temp_series_list = []
        for i in range(0, len(aligned_fasta), 2):
            temp_series_list.append(pd.Series(list(aligned_fasta[i + 1]), name=aligned_fasta[i][1:]))

        # now create the df from the list of series
        # https://github.com/pandas-dev/pandas/issues/1494
        aligned_fasta_as_df = pd.DataFrame.from_items([(s.name, s) for s in temp_series_list]).T
        # aligned_fasta_as_df = pd.DataFrame(temp_series_list)

        # now do the cropping
        aligned_fasta_as_df_cropped = self._crop_fasta_df(aligned_fasta_as_df)

        # now we need to convert this back to a fasta
        output_fasta = []
        for sequence in aligned_fasta_as_df_cropped.index.values.tolist():
            output_fasta.extend(['>{}'.format(aligned_fasta_as_df_cropped.loc[sequence].name),
                                 ''.join(aligned_fasta_as_df_cropped.loc[sequence].values.tolist())])

        return output_fasta

    @staticmethod
    def _convert_interleaved_to_sequencial_fasta(fasta_in):
        fasta_out = []
        for i in range(len(fasta_in)):
            if fasta_in[i].startswith('>'):
                if fasta_out:
                    # if the fasta is not empty then this is not the first
                    fasta_out.append(temp_seq_str.upper())
                # else then this is the first sequence and there is no need to add the seq.
                temp_seq_str = ''
                fasta_out.append(fasta_in[i])
            else:
                temp_seq_str = temp_seq_str + fasta_in[i]
        # finally we need to add in the last sequence
        fasta_out.append(temp_seq_str)
        return fasta_out

    def _create_color_dict(self, color_map_file_path):
        if color_map_file_path is None:
            # make using the list of colors
            color_list, grey_list = self._get_colour_lists()
            color_len = len(color_list)
            temp_count = 0
            temp_dict = {}
            for seq_key in self.fasta_dict.keys():
                if temp_count < color_len:
                    temp_dict[seq_key] = color_list[temp_count]
                else:
                    temp_dict[seq_key] = grey_list[temp_count%6]
                temp_count += 1
            return temp_dict
        else:
            with open(color_map_file_path, 'r') as f:
                colour_file = [line.rstrip() for line in f]
            temp_dict = {line.split(',')[0]: line.split(',')[1] for line in colour_file}
            # check that the seqs match the seqs in the colour dict
            if set(self.fasta_dict.keys()) != set(temp_dict.keys()):
                raise RuntimeError('Sequence in the fasta file were not found in the color map file.')

    def _create_names_dict(self):
        with open(self.names_file_path, 'r') as f:
            sp_names_file = [line.rstrip() for line in f]
        return {str(line.split('\t')[0]): len(line.split('\t')[1].split(',')) for line in sp_names_file}

    def _create_fasta_dict(self):
        with open(self.fasta_file_path, 'r') as f:
            sp_fasta_file = [line.rstrip() for line in f]
        return {str(sp_fasta_file[i][1:]): sp_fasta_file[i + 1] for i in range(0, len(sp_fasta_file), 2)}

    def _splits_tree_nexus_from_fasta(self):
        new_nexus = []
        new_nexus.append('#NEXUS')
        new_nexus.append('BEGIN taxa;')
        new_nexus.append('\tDIMENSIONS ntax={};'.format(int(len(self.aligned_fasta) / 2)))
        new_nexus.append('TAXLABELS')
        count = 1
        for i in range(0, len(self.aligned_fasta), 2):
            new_nexus.append('[{}]\t{}'.format(count, self.aligned_fasta[i][1:]))
            count += 1
        new_nexus.append(';')
        new_nexus.append('END;')
        new_nexus.append('BEGIN characters;')
        new_nexus.append('\tDIMENSIONS nchar={};'.format(len(self.aligned_fasta[1])))
        new_nexus.append('\tFORMAT')
        new_nexus.append('\t\tdatatype=DNA')
        new_nexus.append('\t\tmissing=?')
        new_nexus.append('\t\tgap=-')
        new_nexus.append('\t\tsymbols="A C G T"')
        new_nexus.append('\t\tlabels=left')
        new_nexus.append('\t\tinterleave=no')
        new_nexus.append('\t;')
        new_nexus.append('MATRIX')
        # now put in the sequences in the style of the mammals.nex example from SplitsTree
        for i in range(0, len(self.aligned_fasta), 2):
            new_nexus.append('{}\t{}'.format(self.aligned_fasta[i][1:], self.aligned_fasta[i + 1].upper()))
        new_nexus.append(';')
        new_nexus.append('END;')
        # finally write in the st_assumption block that will tell SplitsTree to calculate the network
        new_nexus.append('BEGIN st_assumptions;')
        new_nexus.append(
            'CHARTRANSFORM=MedianJoining Epsilon=0 SpringEmbedderIterations=1000 '
            'LabelEdges=false ShowHaplotypes=false SubdivideEdges=false ScaleNodesByTaxa=true;')
        new_nexus.append('end;')
        return new_nexus


    def _get_colour_lists(self):
        colour_palette = self._get_colour_list()
        grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
        return colour_palette, grey_palette

    @staticmethod
    def _get_colour_list():
        colour_list = [
            "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5", "#7A4900",
            "#0000A6",
            "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400",
            "#4FC601",
            "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
            "#B903AA",
            "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
            "#372101",
            "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
            "#0CBD66",
            "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
            "#BEC459",
            "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
            "#FF913F",
            "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7",
            "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
            "#201625",
            "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
            "#CB7E98",
            "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489",
            "#806C66",
            "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5",
            "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58", "#7A7BFF", "#D68E01",
            "#353339",
            "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A",
            "#001325",
            "#02525F", "#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75",
            "#8D8546",
            "#9695C5", "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
            "#00005F",
            "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058",
            "#A45B02",
            "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
            "#F4D749",
            "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE",
            "#C6DC99",
            "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
            "#C6005A",
            "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183", "#AA5199", "#B5D6C3",
            "#A38469",
            "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03",
            "#E4FFFC",
            "#1BE177", "#BCB1E5", "#76912F", "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213",
            "#A76F42",
            "#89412E", "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
            "#BDC9D2",
            "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71",
            "#868E7E",
            "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C",
            "#00B57F",
            "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
        return colour_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='fasta2net is a Python wrapper around SplitsTrees4 '
                    '(https://ab.inf.uni-tuebingen.de/software/splitstree4) and '
                    'networkx (https://networkx.github.io/) to easily plot median joining networks '
                    'from a .fasta or .fasta and .names file input.')
    parser.add_argument("--input_fasta_path",
                        help="The path to the fasta file to be used as input. "
                             "This file can be passed in alone or in conjunction with a .names file.",
                        required=True)
    parser.add_argument("--output_dir", help="Directory in which the network figures should be output", required=False)
    parser.add_argument("--input_names_path",
                        help="The path to a .names file associated with the passed .fasta file.",
                        required=False)
    parser.add_argument("--color_map_path",
                        help="The path to the color map file to be used for coloring the network",
                        required=False)
    parser.add_argument("--spring_pos_iterations", help="The number of iterations used in calculating "
                                                        "vertice (nodes) positions. [1000]", type=int, default=1000)
    parser.add_argument("--spring_pos_k", help='k value used in the splits tree calculation', type=float, default=0.1)
    parser.add_argument("--no_labels", help="If this flag is passed, seq labels will not be plotted. [false]",
                        action="store_true", default=False)
    parser.add_argument("--fig_res", help="The resolution (dpi) to output the .png to. [1200]", type=int, default=1200)
    args = parser.parse_args()
    ftn = Fasta2Net(
        fasta_file_path=args.input_fasta_path,
        output_dir=args.output_dir,
        names_file_path=args.input_names_path,
        color_map_file_path=args.color_map_path,
        spring_pos_iterations=args.spring_pos_iterations,
        spring_pos_k_value=args.spring_pos_k,
        labels=not args.no_labels, dpi=args.fig_res)
    ftn.make_network()