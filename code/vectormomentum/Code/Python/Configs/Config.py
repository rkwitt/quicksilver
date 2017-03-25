"""Config file spec system using YAML"""

import PyCA.Common as common
import yaml
import os.path


class Param(object):
    """
    Simple parameter class

    Should give a default value which is used for the skeleton file and is the
    default if not required and not provided by the user

    Also, children is either None or a dict containing other Params

    The rare and debug flags just control whether or not this option is shown
    in sample output
    """
    def __init__(self, default=None, required=False, rare=False, debug=False,
                 comment=None):
        self.default = default
        self.required = required
        self.comment = comment
        self.rare = rare  # NOTE: These are currently ignored entirely --JDH
        self.debug = debug

    def decoratedString(self):
        """Output as YAML example"""
        s = str(self.default)
        if self.required or self.comment is not None:
            s += "  # "
        if self.required:
            s += "REQUIRED "
        if self.comment is not None:
            s += self.comment

        return s

    def __str__(self):
        return self.decoratedString()



def SpecToDict(spec):
    """Given a spec, generate a sample dict, as read by LoadYAMLDict()"""
    d = {}
    for k, v in spec.iteritems():
        if k[0] is '_':
            continue  # stuff with _ at the beginning is probably a hook
        elif isinstance(v, Param):
            d[k] = v.default
        elif isinstance(v, dict):
            d[k] = SpecToDict(v)
        else:
            raise Exception("Malformed config spec. " +
                            "Should be dict tree with Param leaves")
    return d


def DictKeysToAttributes(d):
    """Given a dict, convert key/value pairs to attributes"""

    if isinstance(d, dict):
        o = _Object()
        for k, v in d.iteritems():
            if k in o.__dict__:
                raise Exception("Key " + k + " exists as attribute." +
                                " Error in config file spec")
            else:
                o.__dict__[k] = DictKeysToAttributes(v)
    else:  # we're at a leaf
        return d

    return o


def SpecToConfig(spec):
    """Given a spec, generate a sample config object

    This is useful, for instance, in a test script where one might call this to
    get a config object, then customize a few parameters before running the
    test.
    """
    return DictKeysToAttributes(SpecToDict(spec))


def SpecToYAML(spec, level=0):
    """Convert a config spec to YAML"""
    lines = []
    for k, v in spec.iteritems():
        if k[0] is '_':
            continue  # stuff with _ at the beginning is probably a hook
        elif isinstance(v, Param):
            # Only way to do None is missing entry, so comment out if default
            # value is None
            com = "# " if v.default is None else ""

            sv = v.decoratedString()
            lines.append("%s%s%s: %s" % ("    "*level, com, k, sv))
        elif isinstance(v, dict):
            lines.append("%s%s:" % ("    "*level, k))
            lines.append(SpecToYAML(v, level+1))
        else:
            raise Exception("Malformed config spec. " +
                            "Should be dict tree with Param leaves")
    return '\n'.join(lines)


def ConfigToYAML(spec, d, level=0):
    """Convert a config, given spec to YAML"""
    lines = []
    for k, v in spec.iteritems():
        if k[0] is '_':
            continue  # stuff with _ at the beginning is probably a hook
        elif k not in d.__dict__:
            raise Exception("Key %s missing.  Not validated?" % k)
        elif isinstance(v, Param):
            lines.append("%s%s: %s" % ("    "*level, k, d.__dict__[k]))
        elif isinstance(v, dict):
            lines.append("%s%s:" % ("    "*level, k))
            lines.append(ConfigToYAML(spec=spec[k], d=d.__dict__[k],
                                      level=level+1))
        else:
            raise Exception("Malformed config spec. " +
                            "Should be dict tree with Param leaves")
    return '\n'.join(lines)


class IncludeLoader(yaml.Loader):  # pylint: disable-msg=R0901
    """
    Subclassed yaml loader that supports the !include directive
    """
    def __init__(self, stream):
        self._root = os.path.split(stream.name)[0]
        super(IncludeLoader, self).__init__(stream)

    def include(self, node):
        """Include another yaml at this node"""
        filename = os.path.join(self._root, self.construct_scalar(node))
        with open(os.path.expanduser(filename), 'r') as f:
            return yaml.load(f, Loader=IncludeLoader)

IncludeLoader.add_constructor('!include', IncludeLoader.include)


class MissingKeyError(Exception):
    """Exception type for missing config params"""
    def __init__(self, value=None):
        super(MissingKeyError, self).__init__(self, value)
        self.value = value

    def __str__(self):
        return repr(self.value)


def ValidateDict(d, spec, prefix=""):
    """Validate a dict and insert defaults, conforming to spec"""
    for k, v in spec.iteritems():
        if k[0] is '_':
            continue  # stuff with _ at the beginning is probably a hook
        elif isinstance(v, dict):
            # do the recursive thing
            # NOTE: if your spec is a tree then everything below is required
            if k not in d.keys():
                raise MissingKeyError("Required subsection %s is missing" % k)
            else:
                ValidateDict(d[k], spec[k], prefix=prefix + "." + k)
        elif isinstance(v, Param):
            if k not in d.keys():
                if v.required:
                    raise MissingKeyError("Required param %s is missing" % k)
                else:
                    d[k] = v.default
            else:  # they supplied it, we're good
                pass
        else:
            raise Exception("Malformed config spec at key %s.%s" % (prefix, k)
                            + ". Should be dict tree with Param leaves")

    # warn about extra keys (could be misspelled params with default values)
    for k in set(d.keys())-set(spec.keys()):
        # skip stuff starting with underscore
        if k[:1] is not "_":
            print("Warning: config key " + prefix + "." + k +
                  " not covered by spec.  Ignoring...")
        # remove this cruft from the final config
        del d[k]


class _Object(object):
    """Just an empty object that we can set attributes on"""
    def __init__(self):
        pass


def LoadYAMLDict(configFile):
    """Load a yaml file and insert default values"""
    stream = open(configFile, 'r')
    return yaml.load(stream, Loader=IncludeLoader)


class MissingConfigError(Exception):
    """Exception type for missing config file"""
    def __init__(self, value=None):
        super(MissingConfigError, self).__init__(self, value)
        self.value = value

    def __str__(self):
        return repr(self.value)


def RunValidationHooks(c, spec):
    """Recursively search for hooks and run them"""
    for k, v in spec.iteritems():
        if k is "_validation_hook":
            spec["_validation_hook"](c)
        elif isinstance(v, dict):
            RunValidationHooks(c.__dict__[k], v)


def MkConfig(d, spec):
    """
    Given a dict and a spec, validate and return a config object

    If you do not use yaml config files, or are calling functions from a test
    script, this is the preferred way to set up configs.
    """
    ValidateDict(d, spec)

    # convert dict keys to Object attributes (safely)
    c = DictKeysToAttributes(d)

    # run hooks if supplied by spec
    RunValidationHooks(c, spec)

    return c


def Load(spec, argv):
    """Load YAML file, validate it, given spec (inserting defaults) and convert
    to attributes
    """
    if len(argv) < 2:
        print('# Usage: ' + argv[0] + ' <config.yaml>')
        print('# Below is a sample config YAML file')
        print(SpecToYAML(spec))
        print('# PyCA Version: %s (%s)' % (common.PyCAVersion(),
                                           common.PyCABranch()))
        if '_resource' in spec:
            # this won't be printed by default, but let's keep track
            print("_resource: " + spec['_resource'])

        raise MissingConfigError()

    d = LoadYAMLDict(argv[1])

    # TODO: process overrides (given on the cmdline as key=value pairs)

    return MkConfig(d, spec)
