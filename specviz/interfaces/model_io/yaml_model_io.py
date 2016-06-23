#
# Functions in this module support the reading and writing
# of astropy's spectral compound models from/to YAML files.
#

import os
import sys
import re
import dis
import copy
import yaml

from io import StringIO
from ast import literal_eval

from ...third_party.qtpy.QtWidgets import QFileDialog
from ...interfaces.factories import ModelFactory
from ...core.data import GenericSpectrum1DModelLayer

MODEL_FILE_FILTER = "YAML files (*.yaml)"
EXPRESSION_NAME = 'arithmetic behavior'
MODEL_NAME = 'model'


# Helper functions
def _get_model_class_name(function):
    class_string = str(function.__class__)
    return class_string.split('\'>')[0].split(".")[-1]



# ----------  From YAML file to model  ------------------------'

def _ingest_constraints(param_dict):

    bounds = param_dict['constraints']['bounds']
    fixed = param_dict['constraints']['fixed']
    tied = param_dict['constraints']['tied']

    # bounds are tuples stored as strings so the user
    # can read and edit the file using a text editor.
    # They need to be converted back to python tuples.
    for name in bounds:
        bound = literal_eval(bounds[name])
        bounds[name] = (bound[0], bound[1])

    # TODO: re-do this when implementing ties
    # YAML returns different data types depending
    # on the model type. They need to be properly
    # converted.
    for name in fixed:
        if isinstance(fixed[name], str):
            fixed[name] = literal_eval(fixed[name])
            tied[name] = literal_eval(tied[name])

    return bounds, fixed, tied


def _build_single_model(in_map, model_name=None):
    if model_name is None:
        entry_name = list(in_map.keys())[0]
    else:
        entry_name = model_name

    model_name = in_map[entry_name]['class']

    # model names in ModelFactory do not terminate
    # with a redundant '1D' suffix; remove it.
    model_cls = ModelFactory.all_models[model_name[:-2]]

    param_dict = in_map[entry_name]
    name = param_dict['name']

    bounds, fixed, tied = _ingest_constraints(param_dict)

    # the model constructor call can directly handle
    # all parameter constraints, and the name
    model = model_cls(name=name, bounds=bounds, fixed=fixed, tied=tied)

    # parameter values are top level objects in the model
    # instance, unlike other parameter attributes such as
    # bounds and ties. They have to be set explicitly.
    parameters = param_dict['parameters']
    for pname in parameters:
        value = float(param_dict['parameters'][pname])
        setattr(model, pname, value)

    return model


# If no arithmetic behavior expression is know,
# build a compound model by adding together all
# models present in the map.
def _build_additive_model(in_map):
    model = None
    for model_name in in_map[MODEL_NAME]:
        in_model = _build_single_model(in_map[MODEL_NAME], model_name=model_name)
        if model is None:
            model = in_model
        else:
            model += in_model
    return model


# If an arithmetic behavior expression is present,
# use it to build the compound model
def _build_compound_model(in_map):
    model_list = []
    for model_name in in_map[MODEL_NAME]:
        model_list.append(_build_single_model(in_map[MODEL_NAME], model_name=model_name))

    formula = in_map[EXPRESSION_NAME]

    return GenericSpectrum1DModelLayer.from_formula(model_list, formula), formula


def buildModelFromFile(fname):
    """
    Builds a compound model specified in a YAML file.

    This is the main entry point for the 'read from file'
    functionality. The caller is responsible for providing
    the full-path file name.

    Parameters
    ----------
    fname: str
        the ffully qualified file name
    """
    directory = os.path.dirname(fname)

    f = open(fname, "r")
    in_map = yaml.safe_load(f)
    f.close()

    expression = ""

    if MODEL_NAME in in_map:
        # compound model
        if EXPRESSION_NAME in in_map and len(in_map[EXPRESSION_NAME]) > 0:
            model, expression = _build_compound_model(in_map)
        else:
            # add all models together if no formula is present
            model = _build_additive_model(in_map)
    else:
        # single model
        model = _build_single_model(in_map)

    return model, expression, directory



# ----------  From model to YAML file ------------------------'

# Builds a dictionary with model constraints by directly
# referring to the _constraints attribute in a model.
def _build_constraints_dict(model):
    constraints_dict = copy.deepcopy(model._constraints)

    # bounds are stored as strings so
    # they  can be edited by the user.
    for name in constraints_dict['bounds']:
        bound1 = constraints_dict['bounds'][name][0]
        bound2 = constraints_dict['bounds'][name][1]
        constraints_dict['bounds'][name] = "(%s,%s)" % (str(bound1), str(bound2))

    # clean up. This is something that  exists only
    # in single models and is not needed to rebuild
    # the model from its YAML description.
    if 'eqcons' in constraints_dict:
        constraints_dict.pop('eqcons')
        constraints_dict.pop('ineqcons')

    return constraints_dict


# From a single model, builds the dict to be output to YAML file.
def _build_output_dict_single(model):
    model_name = model.name

    param_dict = {}
    for name, value in zip(model.param_names, model.parameters):
        param_dict[name] = str(value)

    constraints = _build_constraints_dict(model)

    model_dict = {
        'name': model_name,
        'class': _get_model_class_name(model),
        'parameters': param_dict,
        'constraints': constraints}

    return model_name, model_dict


# From a compound model, builds the dict to be output to YAML file.
def _build_output_dict_compound(model):
    model_name = model.name

    param_dict = {}
    for parameter_name, value in zip(model.param_names, model.parameters):
        param_dict[parameter_name] = str(value)

    # In a compound model, we don't want the constraints as stored in
    # the _constraints attribute, because these are keyed by the parameter
    # names *in the compound model*, such as 'amplitude_0'. We need keys
    # that relate directly with the underlying Parameter objects. Thus we
    # cannot use method _build_constraints_dict, since it uses a direct
    # copy of the _constraints attribute to drive its machinery.
    constraints_dict = {'bounds':{},'tied':{},'fixed':{}}
    for parameter_name in model.param_names:
        parameter = getattr(model, parameter_name)
        for constraint_name in parameter.constraints:
            constraint = getattr(parameter, constraint_name)
            constraints_dict[constraint_name][parameter_name] = str(constraint)

    model_dict = {
        'name': model_name,
        'class': _get_model_class_name(model),
        'parameters': param_dict,
        'constraints': constraints_dict}

    return model_name, model_dict


# Writes a dict to YAML file.
def _writeToFile(out_model_dict, model_directory, parent):

    fname = QFileDialog.getSaveFileName(parent, 'Save to file', model_directory)[0]

    if len(fname) > 0:
        # enforce correct suffix.
        if not fname.endswith(".yaml"):
            fname += ".yaml"

        f = open(fname, "w")
        yaml.dump(out_model_dict, f,default_flow_style=False)
        f.close()


# Handles the case of a spectral model with a single component. It's
# not strictly a compound model, but saving and retrieving isolated
# components as if they were compound models makes for a simpler
# interface.
def _writeSingleComponentModel(model, model_directory, parent):
    out_model_dict = {}

    model_name, model_dict = _build_output_dict_single(model)
    out_model_dict[model_name] = model_dict

    _writeToFile(out_model_dict, model_directory, parent)


# Handles the case of a compound model
def _writeCompoundModel(compound_model, model_directory, parent, expression):
    out_model_dict = {MODEL_NAME:{}}

    for model in compound_model:
        model_name, model_dict = _build_output_dict_compound(model)
        out_model_dict[MODEL_NAME][model_name] = model_dict

    out_model_dict[EXPRESSION_NAME] = expression

    _writeToFile(out_model_dict, model_directory, parent)


def saveModelToFile(parent, model, model_directory, expression=None):
    """
    Saves spectral model to file.

    This is the main entry point for the 'save to file'
    functionality.

    Parameters
    ----------
    parent : QWidget or None
        optional widget used for screen centering.
    expression: str
        the formula associated with the compound model
    """
    if not hasattr(model, '_format_expression'):
        _writeSingleComponentModel(model, model_directory, parent)
    else:
        _writeCompoundModel(model, model_directory, parent, expression)


#--------------------------------------------------------------------#
#
#  Utility functions that might be used when we implement ties.
#

# Disassembles a tie callable. Ties read from a model
# file are not directly accessible in text form because
# the model file is compiled at import time.
def get_tie_text(tie):
    if tie:
        # dis() only outputs on standard output.....
        keep = sys.stdout
        sys.stdout = StringIO()
        dis.dis(tie)
        assembler_text = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = keep
        result = _parse_assembler_text(assembler_text)
    else:
        result = 'False'
    return result


# This parses the text returned by the disassembler for
# a lambda function that multiplies a constant by a
# variable. That is, we are assuming that ties are coded
# as lambda functions with multiplication by a constant,
# as in the STSDAS' specfit task.
parser = re.compile(r'\(([^)]+)\)') # picks up whatever is enclosed in parenthesis
def _parse_assembler_text(text):
    tokens = parser.findall(text)
    factor = tokens[0]
    lambda_variable_name = tokens[1]
    function_id = tokens[2]
    par_name = tokens[3]

    return "lambda %s: %s *  %s[%s].%s" % \
           (lambda_variable_name,
            factor,
            lambda_variable_name,
            function_id,
            par_name)


