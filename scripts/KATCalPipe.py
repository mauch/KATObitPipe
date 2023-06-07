#! /usr/bin/env python
import ast
import sys
from katim import KATCalibPipe
from argparse import ArgumentParser


def parse_python_assigns(assign_str):
    """
    Parses a string, containing assign statements
    into a dictionary.

    Stolen from https://github.com/ska-sa/katsdpcontim (In util)

    .. code-block:: python

        h5 = katdal.open('123456789.h5')
        kwargs = parse_python_assigns("spw=3; scans=[1,2];"
                                      "targets='bpcal,radec';"
                                      "channels=slice(0,2048)")
        h5.select(**kwargs)

    Parameters
    ----------
    assign_str: str
        Assignment string. Should only contain assignment statements
        assigning python literals or builtin function calls, to variable names.
        Multiple assignment statements should be separated by semi-colons.

    Returns
    -------
    dict
        Dictionary { name: value } containing
        assignment results.
    """

    if not assign_str:
        return {}

    def _eval_value(stmt_value):
        # If the statement value is a call to a builtin, try evaluate it
        if isinstance(stmt_value, ast.Call):
            func_name = stmt_value.func.id

            if func_name not in _BUILTIN_WHITELIST:
                raise ValueError("Function '%s' in '%s' is not builtin. "
                                 "Available builtins: '%s'"
                                 % (func_name, assign_str, list(_BUILTIN_WHITELIST)))

            # Recursively pass arguments through this same function
            if stmt_value.args is not None:
                args = tuple(_eval_value(a) for a in stmt_value.args)
            else:
                args = ()

            # Recursively pass keyword arguments through this same function
            if stmt_value.keywords is not None:
                kwargs = {kw.arg: _eval_value(kw.value) for kw
                          in stmt_value.keywords}
            else:
                kwargs = {}

            return getattr(builtins, func_name)(*args, **kwargs)
        # Try a literal eval
        else:
            return ast.literal_eval(stmt_value)

    # Variable dictionary
    variables = {}

    # Parse the assignment string
    stmts = ast.parse(assign_str, mode='single').body

    for i, stmt in enumerate(stmts):
        if not isinstance(stmt, ast.Assign):
            raise ValueError("Statement %d in '%s' is not a "
                             "variable assignment." % (i, assign_str))

        # Evaluate assignment lhs
        values = _eval_value(stmt.value)

        # "a = b = c" => targets 'a' and 'b' with 'c' as result
        for target in stmt.targets:
            # a = 2
            if isinstance(target, ast.Name):
                variables[target.id] = values

            # Tuple/List unpacking case
            # (a, b) = 2
            elif isinstance(target, (ast.Tuple, ast.List)):
                # Require all tuple/list elements to be variable names,
                # although anything else is probably a syntax error
                if not all(isinstance(e, ast.Name) for e in target.elts):
                    raise ValueError("Tuple unpacking in assignment %d "
                                     "in expression '%s' failed as not all "
                                     "tuple contents are variable names.")

                # Promote for zip and length checking
                if not isinstance(values, (tuple, list)):
                    elements = (values,)
                else:
                    elements = values

                if not len(target.elts) == len(elements):
                    raise ValueError("Unpacking '%s' into a tuple/list in "
                                     "assignment %d of expression '%s' failed. "
                                     "The number of tuple elements did not match "
                                     "the number of values."
                                     % (values, i, assign_str))

                # Unpack
                for variable, value in zip(target.elts, elements):
                    variables[variable.id] = value
            else:
                raise TypeError("'%s' types are not supported"
                                "as assignment targets." % type(target))

    return variables



description = "Calibrate an MVF file."
parser = ArgumentParser(description=description)
parser.add_argument("katdata", nargs='+',
                    help="MVF file our URL to process")
parser.add_argument("--outputdir", default='./',
                    help="Specify the output data directory.")
parser.add_argument("--parms", dest='parmFile',
                    help="Overwrite the default imaging parameters using a parameter file.")
parser.add_argument("--scratchdir",
                    help="Specify the scratch directory.")
parser.add_argument("--targets",
                    help="List of targets to load (You'll need calibrators in this list!!).")
parser.add_argument("--config", dest='configFile',
                    help="Location of .katimrc configuration file.")
parser.add_argument("--timeav", default=1, type=int,
                    help="Number of dumps to average when making uvfits file")
parser.add_argument("--flag", action='store_true', default=False,
                    help="Flag data on the fly during conversion")
parser.add_argument("--reuse", action='store_true', default=False,
                    help="Reuse already loaded data from aipsdisk")
parser.add_argument("--zapraw", action='store_true', default=False,
                    help="zap raw and intermediate uv files")
parser.add_argument("--aipsdisk", default='aipsdisk',
                    help='Name of aipsdisk (in "scratchdir" - or cwd) to use - default is "aipsdisk"')
parser.add_argument("--halfstokes", default=False, action='store_true',
                    help='Only write out HH,VV when saving uv data')
parser.add_argument("--gzip", default=False, action='store_true',
                    help='Gzip the output UV data')
parser.add_argument("--dropants",
                    help='List of antennas to remove from observation.')
parser.add_argument("--blmask", type=float, default=1.e10,
                    help='Baseline length cutoff for the static mask (default apply to all baselines)')
parser.add_argument("--refant", type=str, default=None,
                    help='Reference antenna to use for calibration')
parser.add_argument("--katdal_options", type=parse_python_assigns, default='',
                    help='Options to pass to katdal.open() as kwargs. '
                         'Should only contain python assignment statements to python '
                         'literals, separated by semi-colons. Default: None')
parser.add_argument("--polcal", default=False, action='store_true',
                    help='Switch on polarisation calibration. '
						 'Fix the X & Y gains. This will cause the selected target (via --XTYarg) alone '
 						 'to be used for delay and bandpass then subsequent gains to be computed with avgPol=True')
parser.add_argument("--XYtarg", default=None,
                    help='Name of target to use to fix the X & Y gains (default is 1934-638 or 0408-65)')
parser.add_argument("--delaycal_mvf", default=None, type=str,
                    help='Dataset to use for determining the XY phase calibration. This '
						 'will override the automatic detection of the relevant dataset. NOTE: The dataset MUST have '
						 'been observed using calibrate_delays.py just before the start of the observation. '
						 'This option is only used when --polcal is selected.')
options = parser.parse_args()

kwargs = {}
for k in ['parmFile', 'scratchdir', 'targets', 'configFile', 'timeav', 'flag', 'reuse', 'zapraw', 'aipsdisk', 'halfstokes',
          'gzip', 'dropants', 'blmask', 'refant', 'katdal_options', 'polcal', 'XYtarg', 'delaycal_mvf']:
	if getattr(options,k) != None:
		kwargs[k] =  getattr(options,k)
try:
    KATCalibPipe.MKContPipeline(options.katdata, options.outputdir, **kwargs)
finally:
    pass
