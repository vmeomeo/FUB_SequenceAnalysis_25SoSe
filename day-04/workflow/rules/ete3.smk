rule visualize_tree:
    input:
        treefile = f"{config['output_dir_path']}/phylogeny/{{sample}}_{{read_number}}/tree.treefile"
    output:
        image = f"{config['output_dir_path']}/phylogeny/{{sample}}_{{read_number}}/tree.png"
    conda:
        "../env/ete3.yaml"
    log:
        f"{config['output_dir_path']}/phylogeny/ete3_{{sample}}_{{read_number}}.log"
    shell:
#         """
#         python -c "
# import sys
# from ete3 import Tree
# t = Tree('{input.treefile}', format=0)
# t.render('{output.image}', w=800, units='px')
# " > {log} 2>&1
#         """
        """
        mkdir -p $(dirname {output.image})
        xvfb-run -a python -c "
import sys
from ete3 import Tree
t = Tree('{input.treefile}', format=0)
t.render('{output.image}', w=800, units='px')
"
        """ 
