static PyObject * hello_wrapper(PyObject * self, PyObject * args)
{
  char * input;
  char * result;
  PyObject * ret;

  // parse arguments
  if (!PyArg_ParseTuple(args, "s", &input)) {
    return NULL;
  }

  // run the actual function
  result = hello(input);

  // build the resulting string into a Python object.
  ret = PyString_FromString(result);
  free(result);

  return ret;
}
