##################################################
###  Genotype-phenotype overlap (GPO):         ###
###  2009, David Ellinghaus                    ###
###  d.ellinghaus(at)ikmb(dot)uni-kiel(dot)de  ###
##################################################


def functionId(nFramesUp):
    """ Create a string naming the function n frames up on the stack.  """
    co = sys._getframe(nFramesUp+1).f_code
    return "%s (%s @ %d)" % (co.co_name, co.co_filename, co.co_firstlineno)
