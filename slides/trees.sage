import sys; sys.path.insert(0, '.')
from feynsage import *
from feynsage.plotting import draw_panels
import matplotlib; matplotlib.use('Agg')
g = diagram("kite")
xs = lambda rem: "$" + "".join("x_{%d}" % (i + 1) for i in rem) + "$"
trees = g.spanning_trees()
fig = draw_panels(g, trees, titles=[xs(r) for r in trees], momenta=False, ncols=4, size=1.6)
fig.savefig('slides/out/kite_trees.svg', bbox_inches='tight')
print(len(trees))
