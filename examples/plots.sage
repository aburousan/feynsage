# Diagrams and plots in the style of the lecture notes.  Run:  sage examples/plots.sage
import sys; sys.path.insert(0, '.')
from feynsage import *
from feynsage.plotting import set_theme, curves, mb_plane, PINK, BLUE
set_theme()
kinE = Kinematics(['p'], ['pp'], {('p','p'): 'pp'}, euclidean=True)
kite = FeynmanGraph([('L','T',0), ('L','B',0), ('T','R',0), ('B','R',0), ('T','B',0)],
                    {'L': mom(p=1), 'R': mom(p=-1)}, kinE)
kite.plot(dots={3: 1}).savefig('examples/kite.svg', bbox_inches='tight')
kinm = Kinematics(['p'], ['pp', 'm2'], {('p','p'): 'pp'}, euclidean=True)
sunset = FeynmanGraph([('A','B','m2'), ('A','B','m2'), ('A','B',0)], {'A': mom(p=1), 'B': mom(p=-1)}, kinm)
sunset.plot().savefig('examples/sunset.svg', bbox_inches='tight')
kin4 = Kinematics(['p1','p2','p3'], ['s','t'],
                  {('p1','p1'): 0, ('p2','p2'): 0, ('p3','p3'): 0, ('p1','p2'): 's/2', ('p2','p3'): 't/2', ('p1','p3'): '-s/2-t/2'},
                  euclidean=False)
box = FeynmanGraph([('a','b',0), ('b','c',0), ('c','d',0), ('d','a',0)],
                   {'a': mom(p1=1), 'b': mom(p2=1), 'c': mom(p3=1), 'd': mom(p1=-1, p2=-1, p3=-1)}, kin4)
box.plot().savefig('examples/box.svg', bbox_inches='tight')
x = var('x')
curves([x/(x^2 + 0.1^2), 1/x], (-3, 3), labels=[r'$x/(x^2+\eta^2)$, $\eta=0.1$', '$1/x$'], ylim=(-6, 6)).savefig('examples/pv.svg')
mb_plane([-1.4, -2.4, -3.4], [0, 1, 2, 3], -0.7, r'$\Gamma(\lambda+z)$', r'$\Gamma(-z)$', close='right').savefig('examples/mb.svg', bbox_inches='tight')
# Sage's own plot() picks up the same fonts and colours after set_theme()
plot(gamma(x), (x, -3.5, 4), ymin=-6, ymax=6, color=PINK, thickness=1.6, detect_poles=True, axes_labels=['$x$', r'$\Gamma(x)$'],
     frame=False).save('examples/gamma.svg')
print("written examples/*.svg")
