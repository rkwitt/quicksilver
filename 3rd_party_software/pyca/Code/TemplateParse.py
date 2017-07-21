import sys
import re
from optparse import OptionParser

defPatStr = '\s*def\s+(\S*)\s*=\s*\[(.*?)\]'
defPat = re.compile(defPatStr, re.DOTALL)

tplPatStr = '\s*(template.*?;)'
tplPat = re.compile(tplPatStr, re.DOTALL)

replPatStr = '(.*?)\$\{(\S*?)\}(.*)'
replPat = re.compile(replPatStr, re.DOTALL)


def Write(codeStr, outStream=sys.stdout):
    outStream.write(codeStr)

    
def WriteComment(commentStr, outStream=sys.stdout):
    lines = commentStr.splitlines()
    for l in lines:
        outStream.write('// ' + l + '\n')

        
def RecursiveTemplGen(parsed, rest, defDict, outStream=sys.stdout, verbose=False):
    if verbose:
        print 'parsed: ', parsed, 'rest: ', rest
    m = replPat.match(rest)
    if m:
        defName = m.group(2)
        if not defName in defDict:
            raise Exception('Error, definition of ' + defName + ' not found')
        for val in defDict[defName]:
            toRep = "${"+defName+"}"
            RecursiveTemplGen(parsed=parsed + m.group(1) + val,
                              rest=m.group(3).replace(toRep, val),
                              defDict=defDict,
                              outStream=outStream,
                              verbose=verbose)
    else:
        Write('%s%s\n' % (parsed, rest), outStream)

        
def DoParse(text, outStream=sys.stdout, verbose=False, clean=True):
    defList = defPat.findall(text)
    tplList = tplPat.findall(text)

    defDict = {}
    for d in defList:
        defName = d[0]
        l = d[1].split(',')
        defVals = [v.strip() for v in l]
        defDict[defName] = defVals

    if verbose:
        print defDict

    if not clean:
        WriteComment('\nTemplate Expansion Definitions\n\n', outStream)
        for k in defDict:
            defStr = k + ' = [' + ','.join(defDict[k]) + ']\n'
            WriteComment(defStr, outStream)
            WriteComment('\n', outStream)

    for templFmt in tplList:
        if not clean:
            WriteComment('\n', outStream)
            WriteComment(templFmt, outStream)
            WriteComment('\n', outStream)
        RecursiveTemplGen(parsed='', rest=templFmt,
                          defDict=defDict,
                          outStream=outStream,
                          verbose=verbose)
        Write('\n', outStream)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-c", "--clean",
                      action="store_true", dest="clean", default=False,
                      help="output only template definitions, not extra comments")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="print status messages to stdout")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        sys.exit(sys.argv[0] + ': Error, must specify input file')

    infile = args[0]

    # read input text
    try:
        with open(infile, 'r') as f:
            text = f.read()
    except IOError as e:
        sys.exit(sys.argv[0] + ': Error, could not read ' + infile)
    
    try:
        outfile = 'stdout'
        outStream = sys.stdout
        if len(args) > 1:
            outfile = args[1]
            outStream = open(outfile, 'w')
        DoParse(text, outStream=outStream,
                verbose=options.verbose, clean=options.clean)
    except IOError as e:
        sys.exit(sys.argv[0] + ': Error, could not open output file ' + outfile)
    finally:
        if outStream is not sys.stdout:
            outStream.close()
