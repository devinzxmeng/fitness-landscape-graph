Here are some of my revision to your plan:
- overview: for printing peak nodes, fitness, and group size (sorted), also printing top-k connection groups with same information set
- for visualization, maybe we can consider plotly, and when mouse cursor hop onto a node, the id will show, then we know what it is. I think this is a better solution than interactive selection
- put the code into a new file so that it doesn't interfere with the current graph_analysis.py. We can always refactor code later.