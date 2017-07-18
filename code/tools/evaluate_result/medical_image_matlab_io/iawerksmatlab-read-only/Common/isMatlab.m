function b = isMatlab()
is = getInterpreter();
b = strcmp(is,'matlab');
