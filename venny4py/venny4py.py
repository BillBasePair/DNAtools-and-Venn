import matplotlib.pyplot as plt

def venny4py(labels, sets, overlay=True):
    """
    Draws a 4-set Venn-like diagram using matplotlib.
    `labels`: list of 4 strings, one for each set
    `sets`: list of 4 sets of values
    """
    from matplotlib_venn import venn3
    import itertools

    plt.figure(figsize=(9, 9))
    ax = plt.gca()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)

    # Draw four ellipses
    ellipse_props = dict(edgecolor='black', linewidth=1.5, alpha=0.5)
    ellipses = [
        plt.Circle((4, 6), 3.5, color='red', **ellipse_props),
        plt.Circle((6, 6), 3.5, color='green', **ellipse_props),
        plt.Circle((5, 8), 3.5, color='blue', **ellipse_props),
        plt.Circle((5, 4), 3.5, color='yellow', **ellipse_props)
    ]

    for e in ellipses:
        ax.add_patch(e)

    # Optional overlay labels (A–O) for intersection regions
    if overlay:
        grid = {
            'A': (3, 6), 'B': (7, 6), 'C': (5, 9), 'D': (5, 3),
            'E': (4.2, 7.5), 'F': (5.8, 7.5), 'G': (3.6, 4.5),
            'H': (6.4, 4.5), 'I': (4.1, 6.1), 'J': (5.9, 6.1),
            'K': (5, 5), 'L': (4.7, 5.3), 'M': (5.3, 5.3),
            'N': (5.3, 4.7), 'O': (4.7, 4.7)
        }
        for label, (x, y) in grid.items():
            plt.text(x, y, label, fontsize=12, ha='center', va='center', color='black', fontweight='bold')

    for i, label in enumerate(labels):
        plt.text(*[(3, 9), (7, 9), (2.5, 3), (7.5, 3)][i], label, fontsize=13, ha='center', fontweight='bold')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    plt.box(False)
    return plt
