import networkx as nx
import matplotlib.pyplot as plt
import argparse


def main(filename: str):
    G = nx.Graph()

    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) != 3:
                continue
            node_start, node_end, distance = parts
            G.add_edge(node_start, node_end)
            # , weight=float(distance))

    pos = nx.spring_layout(G, seed=42)
    # edge_labels = nx.get_edge_attributes(G, 'weight')

    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500, font_size=10)
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("txt")
    args = parser.parse_args()
    main(args.txt)