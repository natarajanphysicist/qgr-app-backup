# lqg_simulation/plotting/visualize.py
"""
Provides functions for visualizing spin networks.
"""
import matplotlib.pyplot as plt
import networkx as nx
from lqg_simulation.core.spin_network import SpinNetwork

def plot_spin_network(spin_net: SpinNetwork, title: str = "Spin Network", show_node_names: bool = True, show_link_spins: bool = True, ax=None):
    """
    Plots a 2D representation of a SpinNetwork using NetworkX and Matplotlib.

    Args:
        spin_net: The SpinNetwork object to plot.
        title: The title of the plot.
        show_node_names: Whether to display node names.
        show_link_spins: Whether to display link spin values.
        ax: A Matplotlib Axes object to draw on. If None, a new figure and axes are created.
    """
    if not isinstance(spin_net, SpinNetwork):
        raise TypeError("Input must be a SpinNetwork object.")

    G = nx.Graph()

    # Add nodes to the NetworkX graph
    # We use node.name as the identifier in NetworkX for simplicity if names are unique.
    # Otherwise, node.id should be used. For plotting, names are more readable if available.
    node_map = {node.id: node for node in spin_net.nodes}
    node_labels = {}
    for node_id, node_obj in node_map.items():
        G.add_node(node_obj.name) # Use name for display
        if show_node_names:
            node_labels[node_obj.name] = node_obj.name
        else:
            node_labels[node_obj.name] = ""


    # Add links (edges) to the NetworkX graph
    link_labels = {}
    for link_obj in spin_net.links:
        # Ensure nodes are in G by their names
        node1_name = link_obj.node1.name
        node2_name = link_obj.node2.name

        if not G.has_node(node1_name): G.add_node(node1_name) # Should already be there
        if not G.has_node(node2_name): G.add_node(node2_name) # Should already be there

        G.add_edge(node1_name, node2_name, spin=link_obj.spin_j)
        if show_link_spins:
            # Position link labels at the midpoint of the edge
            link_labels[(node1_name, node2_name)] = f"j={link_obj.spin_j}"

    # Plotting
    if ax is None:
        fig, current_ax = plt.subplots(figsize=(10, 8))
    else:
        current_ax = ax
        fig = current_ax.figure


    # Use a layout algorithm from NetworkX
    # Spring layout is common, but others might be better for specific structures.
    # pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42) # k adjusts spacing, seed for reproducibility
    pos = nx.kamada_kawai_layout(G) # Often good for small to medium graphs
    # pos = nx.shell_layout(G)
    # pos = nx.circular_layout(G)


    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=current_ax, node_color='skyblue', node_size=1500, alpha=0.9)

    # Draw edges
    nx.draw_networkx_edges(G, pos, ax=current_ax, width=1.5, alpha=0.7, edge_color='gray')

    # Draw node labels
    if show_node_names:
        nx.draw_networkx_labels(G, pos, labels=node_labels, ax=current_ax, font_size=10, font_weight='bold')

    # Draw edge (link) labels
    if show_link_spins:
        nx.draw_networkx_edge_labels(G, pos, edge_labels=link_labels, ax=current_ax, font_size=8, font_color='red')

    current_ax.set_title(title, fontsize=15)
    current_ax.axis('off') # Turn off the axis box and ticks

    if ax is None: # If we created the figure, show it.
        plt.tight_layout()
        # Instead of plt.show(), which blocks, we might want to save to a file in some contexts
        # For interactive use, plt.show() is fine. For automated runs, saving is better.
        # For now, let's assume interactive or a context where show() is handled by the caller.
        # If running in a script, plt.show() would typically be called by the script.
        pass # The caller of this function will handle plt.show() or saving.

    return fig, current_ax


if __name__ == '__main__':
    from lqg_simulation.core import SpinNetwork # For example usage
    from sympy import S

    # Create a sample spin network
    sn_example = SpinNetwork()

    n1 = sn_example.add_node(node_name="N1")
    n2 = sn_example.add_node(node_name="N2")
    n3 = sn_example.add_node(node_name="N3")
    n4 = sn_example.add_node(node_name="N4")
    n5 = sn_example.add_node(node_name="N5") # An isolated node

    sn_example.add_link(n1, n2, spin_j=S(1)/2, link_name="L12")
    sn_example.add_link(n2, n3, spin_j=S(1), link_name="L23")
    sn_example.add_link(n3, n1, spin_j=S(3)/2, link_name="L31")
    sn_example.add_link(n3, n4, spin_j=S(1)/2, link_name="L34")
    sn_example.add_link(n2, n4, spin_j=S(1), link_name="L24")

    # Plot the spin network
    fig, ax = plot_spin_network(sn_example, title="Sample Spin Network Visualization")

    # In a script, you would call plt.show() here to display the plot.
    # To save the plot:
    # plt.savefig("sample_spin_network.png")
    # For this __main__ example, let's ensure it can run headlessly if needed.
    # If you run `python visualize.py`, it won't show interactively without plt.show().
    # For testing purposes, we might not want it to block.
    print("Spin network plot generated (not shown interactively in this example run).")
    plt.close(fig) # Close the figure to free memory, especially if generating many plots.

    # Example with an empty spin network
    sn_empty = SpinNetwork()
    fig_empty, _ = plot_spin_network(sn_empty, title="Empty Spin Network")
    print("Empty spin network plot generated.")
    plt.close(fig_empty)

    # Example with a single node
    sn_single_node = SpinNetwork()
    sn_single_node.add_node(node_name="S1")
    fig_single, _ = plot_spin_network(sn_single_node, title="Single Node Spin Network")
    print("Single node spin network plot generated.")
    plt.close(fig_single)

    # Example with two nodes and one link
    sn_two_nodes = SpinNetwork()
    tn1 = sn_two_nodes.add_node(node_name="TN1")
    tn2 = sn_two_nodes.add_node(node_name="TN2")
    sn_two_nodes.add_link(tn1, tn2, spin_j=S(2))
    fig_two, _ = plot_spin_network(sn_two_nodes, title="Two Nodes, One Link")
    print("Two nodes, one link plot generated.")
    plt.close(fig_two)

    print("Completed visualization examples.")
