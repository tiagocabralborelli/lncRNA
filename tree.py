
from ete3 import Tree, NodeStyle, TreeStyle, NodeStyle, faces, AttrFace, CircleFace
from ete3 import NCBITaxa
ncbi = NCBITaxa()
tree = ncbi.get_topology(['1480154','13451','145481'])
print (tree.get_ascii(attributes=["sci_name", "rank"]))
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
tree.render("teste.png", tree_style=ts)
