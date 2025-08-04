import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import seaborn as sns
from os import listdir
from scipy.stats import zscore
import networkx as nx
import random
import matplotlib.patches as mpatches

def CreateLinkDf (count_matrix):
    count_matrix_turned = count_matrix.T
    #Calculate the Pearson Correlation
    correlation_matrix = count_matrix_turned.corr().to_numpy().round(decimals=3, out=None)
    #Removing self loops
    np.fill_diagonal(correlation_matrix, 0)
    #Covert the np matrix to padnas df
    df = pd.DataFrame(correlation_matrix, columns = list(count_matrix.index), index= list(count_matrix.index))
    #Define the edge list
    links = df.stack().reset_index()
    #Assign the column name as [source , target , value]
    links.columns = ["source" , "target" , "value"]
    return links

def DensityEvaluation (ths , links , verbose=1):
    density = []
    node_number = []
    edge_number = []
    average_degree = []
    if verbose:
        for th in ths :
            print("starting evaluatiom on {}".format(th))
            links_filtered=links.loc[ (links['value'] > th ) & (links['source'] != links['target']) ]
            print("rows removed")
            g = ig.Graph.TupleList(links_filtered.itertuples(index=False), directed=False, weights=True)
            print("netowrk loaded")
            density.append( g.density(loops=False))
            node_number.append( g.vcount())
            edge_number.append(g.ecount())
            average_degree.append(ig.mean(g.degree()))
            #delete all the nodes in the graph __> free space and avoid other problems
            to_delete_ids = [v.index for v in g.vs ]
            g.delete_vertices(to_delete_ids)
            print("memory cleaned")
            print("ended evaluatiom on {}".format(th))
    else :
        for th in ths :
            links_filtered=links.loc[ (links['value'] > th ) & (links['source'] != links['target']) ]
            g = ig.Graph.TupleList(links_filtered.itertuples(index=False), directed=False, weights=True)
            density.append( g.density(loops=False))
            node_number.append( g.vcount())
            edge_number.append(g.ecount())
            average_degree.append(ig.mean(g.degree()))
            #delete all the nodes in the graph __> free space and avoid other problems
            to_delete_ids = [v.index for v in g.vs ]
            g.delete_vertices(to_delete_ids)
    return density,node_number,edge_number, average_degree


def PlotDensity (ths, density,title) :
    fig = plt.figure(figsize=(7,5))
    color = 'tab:red'
    plt.set_title(title)
    plt.set_xlabel('tresholds')
    plt.set_ylabel('density', color=color)
    plt.plot(ths, density, color=color)
    plt.tick_params(axis='y', labelcolor=color)
    plt.show()


def genomic_position_stackbar (communities : list  , annotation_df : pd.DataFrame) : 
    final_annotation_df = pd.DataFrame()
    for num , community in enumerate(communities) :
        community_anontation = annotation_df.loc[list(communities)[num],:]
        community_anontation["community_numebr"] = "community_" + str(num)
        final_annotation_df = pd.concat([final_annotation_df  , community_anontation] , axis = 0 )

    return final_annotation_df


def transform_annotation_homer (path_file : str , filter_up = 0 , filter_down = 0) :
    annotation_gene = pd.read_csv(path_file,
                         delimiter="\t")
                         
    if filter_up != 0 :
        annotation_gene = annotation_gene.loc[(annotation_gene["Distance to TSS"] <= filter_up) & (annotation_gene["Distance to TSS"] >= filter_down) , :]

    annotation_gene = annotation_gene[["Chr","Start","End" ,'Annotation','Gene Name']]
    annotation_gene[['GenomicRegion', '_']] = annotation_gene.Annotation.str.split(" \\(", expand = True)
    annotation_gene['Start'] =  annotation_gene['Start'] -1
    annotation_gene["Sites"] = annotation_gene['Chr'].astype(str) +"-"+ annotation_gene["Start"].astype(str)+"-"+ annotation_gene["End"].astype(str)
    annotation_gene.drop( ['Chr','Start','End' , 'Annotation' , '_'] , axis = 1 , inplace=True)
    annotation_gene.set_index("Sites" , inplace=True)
    return annotation_gene


# Loading Communities
def loading_comminity_from_folder (path_file : str) :
    onlyfiles = [f for f in listdir(path_file)]
    communities = pd.DataFrame()
    for f in onlyfiles :
        df = pd.read_csv(path_file+f ,delimiter="\t" , names = ['chr','start','end'])
        df["sites"] = df['chr'].astype(str) +"-"+ df["start"].astype(str)+"-"+ df["end"].astype(str)
        
        df.drop(['chr','start','end'], axis  = 1 , inplace = True)
        df['community'] = f.split(".tsv")[0]
        communities = pd.concat([communities , df])
    return communities 



def set_node_community(G , communities) :
    '''Add community to node attributes'''
    for c , v_c in enumerate(communities) :
        for v in v_c :
            #Add 1 to save 0 for external edges
            G.nodes[v]['community']= c + 1

# Create Subgrphs for each community
def communities_subgraphs( G , communities) :
    subgraphs_list = []
    for community in communities :
       subgraphs_list.append(G.subgraph(community))
    return subgraphs_list


def plot_trends (communities_list : list , count_df : pd.DataFrame , dim_x : int , dim_y : int  , custom_ylim = (-2.5,10)) :
    x = 0
    community_number = 0
    fig, axes = plt.subplots(dim_x, dim_y , figsize=(30, 10) , dpi = 300 , constrained_layout=True)
    plt.setp(axes, ylim=custom_ylim)

    while x < dim_x :
        y = 0
        while y < dim_y  :
            #print (x,y,community_number)
            sns.violinplot(count_df.loc[communities_list[community_number],:] , ax = axes[x , y])
            axes[x,y].set_title(  'Community ' + str(community_number), fontstyle='italic')
            axes[x,y].set(xlabel=None)
            axes[x,y].set(ylabel="TMM")
            y = y+1

            community_number = community_number + 1 
            if community_number == len(communities_list) :
                break
        x = x + 1


def plot_trends_zscore (communities_list_toplot : list , data_matrix : pd.DataFrame , custom_ylim = (-2,2)) :
    
    count_df_scored = data_matrix.apply(zscore , axis=1)

    for x in range(0,len(communities_list_toplot)) :
        count_df = count_df_scored.loc[communities_list_toplot[x],:]
        count_df = count_df.stack().reset_index()
        count_df.columns = ['peaks','variable', 'value']
        count_df['community'] = str(x)

        if x == 0 :
            final_df = count_df
        else :
            final_df = pd.concat([final_df , count_df])
        
        plt.figure(figsize=(20, 12))

    # Create the violin plot
    #sns.violinplot(data=final_df,x="variable", y="value", hue="community", scale='width', inner=None)

    # Facet by Cluster
    g = sns.FacetGrid(data=final_df,row='community', sharey=False, height=3, aspect=2, despine=False) 
    g.map_dataframe(sns.violinplot,x="variable", y="value", hue="community", scale='width', fill=False , color = 'gray')
    g.map_dataframe(sns.lineplot,  x = 'variable' , y = 'value' , markers=False , color='black' , lw=2)
    g.figure.subplots_adjust(wspace=0, hspace=0.1)
    g.set_axis_labels("", "Accessibilty Z-score") 
    g.set_titles(row_template="")

    # Customize aesthetics
    for ax in g.axes.flat:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')
        ax.set_xlabel("variable")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_color('black')
        ax.spines['left'].set_color('black')

    g.despine(left=True)
    # plt.savefig('/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/CellReport/pictures/stacked_trends.png', bbox_inches='tight' , dpi=300)
    return (g)



def get_communities_names(G , communities) :
    communities_list_sp = []
    for n , community in enumerate(communities) :
        if len(community) > 100 :
            community_list = [G.vs[node_pos]["name"] for node_pos in community]
            communities_list_sp.append(community_list)
    return(communities_list_sp)


def plot_clusters_resolution (df , steps ) :
    edges = []
    for i in steps:  # Cluster_1 to Cluster_2, Cluster_2 to Cluster_3
        from_col = f"Cluster_{i}"
        to_col = f"Cluster_{round(i+0.1 , 2)}"


        grouped = df.groupby([from_col, to_col]).size().reset_index(name='count')
        
        for _, row in grouped.iterrows():
            source = f"{from_col}:{row[from_col]}"
            target = f"{to_col}:{row[to_col]}"
            edges.append((source, target, row['count']))

    # Build graph
    G = nx.DiGraph()
    for source, target, weight in edges:
        G.add_edge(source, target, weight=weight)

    # Layout positions by level (res_1 → res_3 vertically)
    level_pos = {}
    y_spacing = 10
    for node in G.nodes:
        res, label = node.split(":")
        level = float(res.split("_")[1])
        if level not in level_pos:
            level_pos[level] = []
        level_pos[level].append(node)

    # Assign (x, y) positions
    pos = {}
    for level, nodes in level_pos.items():
        for i, node in enumerate(nodes):
            pos[node] = (i * 4, -level * y_spacing)  # vertical layout


    levels = sorted({float(n.split(":")[0].split("_")[1]) for n in G.nodes})

    # Generate random colors for each level
    def random_color():
        return "#{:06x}".format(random.randint(0, 0xFFFFFF))

    level_colors = {lvl: random_color() for lvl in levels}


    # Node colors per level
    #level_colors = {1: 'lightblue', 2: 'lightgreen', 3: 'salmon'}
    node_colors = [level_colors[float(n.split(":")[0].split("_")[1])] for n in G.nodes]

    # Draw graph
    plt.figure(figsize=(12, 8))
    nx.draw(
        G, pos, with_labels=False, node_size=1000,
        node_color=node_colors, edge_color='gray', arrows=True, font_size=8
    )

    # Create legends
    node_legend = [
        mpatches.Patch(color=clr, label=f"res_{lvl}")
        for lvl, clr in level_colors.items()
    ]

    edge_weights = nx.get_edge_attributes(G, 'weight')
    min_w, max_w = min(edge_weights.values()), max(edge_weights.values())

    # Map edge width to weight
    edge_widths = [1 + 3 * (w - min_w) / (max_w - min_w) if max_w > min_w else 2 for w in edge_weights.values()]
    nx.draw_networkx_edges(G, pos, width=edge_widths)

    # Edge weight legend (example buckets)
    weight_legend = [
        mpatches.Patch(color='gray', label=f"Weight ≈ {w}")
        for w in sorted(set(edge_weights.values()))[:3]  # take up to 3 unique weights
    ]

    plt.legend(handles=node_legend + weight_legend, loc='lower left', bbox_to_anchor=(1, 0.5))
    plt.title("Cluster Transition Tree (Colored by Resolution Level)")
    plt.axis('off')
    plt.tight_layout()
    plt.show()