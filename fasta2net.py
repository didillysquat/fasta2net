''' This programe will be aimed at creating splits tree median joining networks form fasta names pair files.
It will be good to add colour functionality and gap functionality to it as well.'''

# eventually the fasta and the names file will come from command line arguments but for the time being
# we will hardcode their paths in.
# lets start by running the SP network.

from plumbum import local
import sys
import os
import subprocess
import networkx as nx
from collections import defaultdict
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

def main():

    output_dir = '/Users/humebc/Google_Drive/projects/barbara_forcioli/sp_networks'
    os.makedirs(output_dir, exist_ok=True)

    in_fasta_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/sp_seqs.fasta'
    with open(in_fasta_path, 'r') as f:
        sp_fasta_file = [line.rstrip() for line in f]
    fasta_dict = {str(sp_fasta_file[i][1:]) : sp_fasta_file[i+1] for i in range(0, len(sp_fasta_file), 2)}

    sp_names_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/sp.names'
    with open(sp_names_path, 'r') as f:
        sp_names_file = [line.rstrip() for line in f]
    names_dict = {str(line.split('\t')[0]) : len(line.split('\t')[1].split(',')) for line in sp_names_file}

    colour_path = '/Users/humebc/Google_Drive/projects/barbara_forcioli/sp_colour.csv'
    with open(colour_path, 'r') as f:
        colour_file = [line.rstrip() for line in f]
    colour_dict = {line.split(',')[0]: line.split(',')[1] for line in colour_file}

    # generate an abundance dict from the fasta and the abund
    seq_to_abund_dict = {}
    for seq_name, seq_seq in fasta_dict.items():
        seq_to_abund_dict[seq_seq] = names_dict[seq_name]



    # now perform the alignment with MAFFT
    mafft = local["mafft-linsi"]

    if in_fasta_path.endswith('.fasta'):
        out_file_path = '{}/{}'.format(output_dir, in_fasta_path.split('/')[-1].replace('.fasta', '_aligned.fasta'))
    elif in_fasta_path.endswith('.fas'):
        out_file_path = '{}/{}'.format(output_dir, in_fasta_path.split('/')[-1].replace('.fas', '_aligned.fasta'))
    elif in_fasta_path.endswith('.fa'):
        out_file_path = '{}/{}'.format(output_dir, in_fasta_path.split('/')[-1].replace('.fa', '_aligned.fasta'))
    else:
        sys.exit('fasta file format extension .{} not recognised'.format(in_fasta_path.split('.')[-1]))

    # now run mafft including the redirect
    (mafft['--thread', -1, in_fasta_path] > out_file_path)()

    with open(out_file_path, 'r') as f:
        aligned_fasta_interleaved = [line.rstrip() for line in f]


    aligned_fasta = convert_interleaved_to_sequencial_fasta(aligned_fasta_interleaved)

    aligned_fasta_cropped = crop_fasta(aligned_fasta)


    # we have to make the new nexus format by hand as the biopython version was putting out old versions.
    new_nexus = splits_tree_nexus_from_fasta(aligned_fasta_cropped)

    # write out the new_nexus
    new_nexus_path = '{}/splitstree_in.nex'.format(output_dir)
    with open(new_nexus_path, 'w') as f:
        for line in new_nexus:
            f.write('{}\n'.format(line))


    # now create the control file that we can use for execution for the no med
    splits_out_path = '{}/splitstree_out.nex'.format(output_dir)
    ctrl_path = '{}/splitstree_ctrl'.format(output_dir)

    make_and_write_cntrl_file(ctrl_path, new_nexus_path, splits_out_path)

    # this creates a control file that can be fed to splits tree on the command line and write it out
    # it then runs splitstrees with the cntrl file before returning the output file
    splits_tree_out_file = run_splits_trees(ctrl_path=ctrl_path, splits_out_path=splits_out_path)


    # networkx graph object
    g = nx.Graph()

    # we can work with a minimised version of the splits out file that is only the network block
    for i in range(len(splits_tree_out_file)):
        if splits_tree_out_file[i] == 'BEGIN Network;':
            network_block = splits_tree_out_file[i:]
            break

    # here we can now work with the network_block rather than the larger splits_tree_out_file
    # get a list of the nodes
    # these can be got from the splits_tree_outfile
    # we should use the numbered system that the splits tree is using. We can get a list of the different sequences
    # that each of the taxa represent form the translate part of the network block.
    # This way we can use this to run the sequences against SP to get the names of the sequencs
    # we can also make the 'translation dictionary' at the same time
    vertice_id_to_seq_id_dict = defaultdict(list)
    for i in range(len(network_block)):
        if network_block[i] == 'TRANSLATE':
            # for each line in the translate section
            for j in range(i + 1, len(network_block)):
                if network_block[j] == ';':
                    # then we have reached the end of the translation section
                    break
                items = network_block[j].replace('\'', '').replace(',', '').split(' ')
                vertice_id_to_seq_id_dict[items[0]].extend(items[1:])

    vertices = list(vertice_id_to_seq_id_dict.keys())

    # here we will use the dictionary that was created in the previous method (passed into this one)
    # to link the count id to the actual sequence, which can then be used to look up the abundance
    # This way we can have a vertice_id_to_rel_abund dict that we can use to put a size to the nodes

    vertice_id_to_abund_dict = {}
    for vert_id in vertices:
        count_id_list = vertice_id_to_seq_id_dict[vert_id]
        if len(count_id_list) == 1:
            abund = names_dict[count_id_list[0]]
            vertice_id_to_abund_dict[vert_id] = abund
        elif len(count_id_list) > 1:
            # then we need sum the abundances for each of them
            tot = 0
            for count_id in count_id_list:
                abund = names_dict[count_id]
                tot += abund
            vertice_id_to_abund_dict[vert_id] = tot

    # the sizes are a little tricky to work out becauase I'm not acutally sure what the units are, possibly pixels
    # lets work where if there was only 1 sequence it would be a size of 1000
    # therefore each vertice will be a node size that is the re_abund * 1000
    # NB the size does appear to be something similar to pixels
    vertices_sizes = [vertice_id_to_abund_dict[vert] * 3 for vert in vertices]

    # need to create a colour list for the nodes as well
    vertices_colours = [colour_dict[vertice_id_to_seq_id_dict[vert][0]] for vert in vertices]


    # the edge list should be a list of tuples
    # currently we are not taking into account the length of the vertices
    edges_list = []
    for i in range(len(network_block)):
        if network_block[i] == 'EDGES':
            # for each line in the edges section
            for j in range(i + 1, len(network_block)):
                if network_block[j] == ';':
                    # then we have reached the end of the Edges section
                    break
                items = network_block[j].replace(',', '').split(' ')[1:]
                edges_list.append((items[0], items[1]))

    # THis is useful, scroll down to where the 'def draw_networkx(......' is
    # https://networkx.github.io/documentation/networkx-1.9/_modules/networkx/drawing/nx_pylab.html#draw_networkx
    #


    g.add_nodes_from(vertices)
    g.add_edges_from(edges_list)
    # we should be able to
    f, ax = plt.subplots(1, 1, figsize=(10,10))
    # we will need to set the limits dynamically as we don't know what the positions are going to be.
    # I think we will want to have the x and y limits the same so that we end up with a square
    # we will therefore just look for the bigest and smallest values, add a buffer and then set the
    # axes limits to thsee
    spring_pos = nx.spring_layout(g, k=0.1, iterations=10000)
    max_ax_val = 0
    min_ax_val = 9999
    for vert_pos_array in spring_pos.values():
        for ind in vert_pos_array:
            if ind > max_ax_val:
                max_ax_val = ind
            if ind < min_ax_val:
                min_ax_val = ind

    ax.set_ylim(min_ax_val - 0.2, max_ax_val + 0.2)
    ax.set_xlim(min_ax_val - 0.2, max_ax_val + 0.2)
    # to set the edge colour of the nodes we need to draw them seperately
    # nx.draw_networkx_nodes(g, pos=spring_pos, node_size=vertices_sizes, with_labels=False, alpha=0.5, edgecolors='black')
    nx.draw_networkx(g, pos=spring_pos, arrows=False, ax=ax, node_color=vertices_colours, alpha=1.0,
                     node_size=vertices_sizes, with_labels=False)

    # https://stackoverflow.com/questions/22716161/how-can-one-modify-the-outline-color-of-a-node-in-networkx
    # https://matplotlib.org/api/collections_api.html
    ax.collections[0].set_edgecolor('grey')
    ax.collections[0].set_linewidth(0.5)
    ax.collections[1].set_linewidth(0.5)
    ax.set_axis_off()
    # ax.set_title(type_id)
    apples = 'asdf'
    plt.savefig('{}/sp_network_s3_k0.1.png'.format(output_dir), dpi=1200)
    plt.savefig('{}/sp_network_s3_k0.1.svg'.format(output_dir))



def make_and_write_cntrl_file(ctrl_path, new_nexus_path, splits_out_path):
    ctrl_file = []
    ctrl_file.append('BEGIN SplitsTree;')
    ctrl_file.append('EXECUTE FILE={};'.format(new_nexus_path))
    ctrl_file.append('SAVE FILE={} REPLACE=yes;'.format(splits_out_path))
    ctrl_file.append('QUIT;')
    ctrl_file.append('end;')
    # now write out the control file
    with open(ctrl_path, 'w') as f:
        for line in ctrl_file:
            f.write('{}\n'.format(line))


def convert_interleaved_to_sequencial_fasta(fasta_in):
    fasta_out = []
    for i in range(len(fasta_in)):
        if fasta_in[i].startswith('>'):
            if fasta_out:
                # if the fasta is not empty then this is not the first
                fasta_out.append(temp_seq_str.upper())
            #else then this is the first sequence and there is no need to add the seq.
            temp_seq_str = ''
            fasta_out.append(fasta_in[i])
        else:
            temp_seq_str = temp_seq_str + fasta_in[i]
    #finally we need to add in the last sequence
    fasta_out.append(temp_seq_str)
    return fasta_out

def crop_fasta(aligned_fasta):
    # convert each of the sequences in the fasta into a series with the series name as the sequence name from the fasta
    temp_series_list = []
    for i in range(0, len(aligned_fasta), 2):
        temp_series_list.append(pd.Series(list(aligned_fasta[i+1]), name=aligned_fasta[i][1:]))

    # now create the df from the list of series
    # https://github.com/pandas-dev/pandas/issues/1494
    aligned_fasta_as_df = pd.DataFrame.from_items([(s.name, s) for s in temp_series_list]).T
    # aligned_fasta_as_df = pd.DataFrame(temp_series_list)

    # now do the cropping
    aligned_fasta_as_df_cropped = crop_fasta_df(aligned_fasta_as_df)

    # now we need to convert this back to a fasta
    output_fasta = []
    for sequence in aligned_fasta_as_df_cropped.index.values.tolist():
        output_fasta.extend(['>{}'.format(aligned_fasta_as_df_cropped.loc[sequence].name), ''.join(aligned_fasta_as_df_cropped.loc[sequence].values.tolist())])

    return output_fasta

def splits_tree_nexus_from_fasta(aligned_fasta):
    new_nexus = []
    new_nexus.append('#NEXUS')
    new_nexus.append('BEGIN taxa;')
    new_nexus.append('\tDIMENSIONS ntax={};'.format(int(len(aligned_fasta) / 2)))
    new_nexus.append('TAXLABELS')
    count = 1
    for i in range(0, len(aligned_fasta), 2):
        new_nexus.append('[{}]\t{}'.format(count, aligned_fasta[i][1:]))
        count += 1
    new_nexus.append(';')
    new_nexus.append('END;')
    new_nexus.append('BEGIN characters;')
    new_nexus.append('\tDIMENSIONS nchar={};'.format(len(aligned_fasta[1])))
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
    for i in range(0, len(aligned_fasta), 2):
        new_nexus.append('{}\t{}'.format(aligned_fasta[i][1:], aligned_fasta[i + 1].upper()))
    new_nexus.append(';')
    new_nexus.append('END;')
    # finally write in the st_assumption block that will tell SplitsTree to calculate the network
    new_nexus.append('BEGIN st_assumptions;')
    new_nexus.append(
        'CHARTRANSFORM=MedianJoining Epsilon=0 SpringEmbedderIterations=1000 LabelEdges=false ShowHaplotypes=false SubdivideEdges=false ScaleNodesByTaxa=true;')
    new_nexus.append('end;')
    return new_nexus

def run_splits_trees(ctrl_path, splits_out_path):


    # now run splitstree
    completedProcess = subprocess.run(
        ['SplitsTree', '-g', '-c', ctrl_path])

    # now we can read in the output file
    # and then we can start making the network finally!
    with open(splits_out_path, 'r') as f:
        splits_tree_out_file = [line.rstrip() for line in f]

    return splits_tree_out_file

def crop_fasta_df(aligned_fasta_as_pandas_df_to_crop):
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
    col_to_keep = [col_index for col_index in list(aligned_fasta_as_pandas_df_to_crop) if col_index not in columns_to_drop]
    # drop the gap columns
    return aligned_fasta_as_pandas_df_to_crop[col_to_keep]

main()