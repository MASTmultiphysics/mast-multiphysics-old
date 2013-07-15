import os, sys

inp = sys.argv[1]
dirval = sys.argv[2]
dirval += '/'

f = open(inp,'r')

strval = ''
for ln in f:
    if (ln.count(':')):
        lns = ln.split(':')
        target_name = lns[0].split('.o')
        strval += (dirval+target_name[0]+'-${METHOD}-cpp.o ')
        strval += (dirval+target_name[0]+'-cpp.dep: ')
        strval += lns[1]
    else:
        strval += ln
f.close();

f = open(inp, 'w')
f.write(strval)
f.close()
