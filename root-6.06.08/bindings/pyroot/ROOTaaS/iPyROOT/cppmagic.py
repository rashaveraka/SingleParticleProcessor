
from IPython.core.magic import (Magics, magics_class, cell_magic)
from IPython.core.magic_arguments import (argument, magic_arguments, parse_argstring)
import utils


@magics_class
class CppMagics(Magics):
    @cell_magic
    @magic_arguments()
    @argument('-a', '--aclic', action="store_true", help='Compile code with ACLiC.')
    @argument('-d', '--declare', action="store_true", help='Declare functions and/or classes.')
    def cpp(self, line, cell):
        '''Executes the content of the cell as C++ code.'''
        args = parse_argstring(self.cpp, line)
        if args.aclic:
            utils.invokeAclic(cell)
        elif args.declare:
            utils.declareCppCode(cell)
        else:
            utils.processCppCode(cell)

def load_ipython_extension(ipython):
    ipython.register_magics(CppMagics)
