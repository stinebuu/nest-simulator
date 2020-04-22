# -*- coding: utf-8 -*-
#
# hl_api_server.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

import nest

import array
import inspect
import io
import numpy as np
import os
import sys

import flask
from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin

__all__ = [
    'app',
    'do_exec',
    'run_mpi_app',
    'serialize'
]

app = Flask(__name__)
CORS(app)

mpi_comm = None


@app.route('/', methods=['GET'])
@cross_origin()
def nest_index():
    """ Route to fetch metadata of NEST Server.
    """
    data = init_data(request)
    data['response']['version'] = nest.version()
    data['response']['is_mpi'] = mpi_comm is not None
    return jsonify(data)


def do_exec(data, *args, **kwargs):
    with Capturing() as stdout:
        try:
            source = kwargs.get('source', '')
            globals_ = {'__builtins__': None}
            locals_ = {
                'list': list,
                'nest': nest,
                'np': np,
                'print': print,
                'set': set,
                'int': int,
            }
            exec(source, globals_, locals_)
            if 'return' in kwargs:
                if isinstance(kwargs['return'], list):
                    return_data = {}
                    for variable in kwargs['return']:
                        return_data[variable] = locals_.get(variable, None)
                else:
                    return_data = locals_.get(kwargs['return'], None)
                data['response']['data'] = nest.hl_api.serializable(return_data)
            data['response']['status'] = 'ok'
        except Exception as e:
            print(e)
            data['response']['data'] = None
            data['response']['status'] = 'error'
    data['response']['stdout'] = '\n'.join(stdout)


@app.route('/exec', methods=['GET', 'POST'])
@cross_origin()
def route_exec():
    """ Route to execute script in Python.
    """
    data = init_data(request)
    args, kwargs = get_arguments(request, data)
    if mpi_comm is not None:
        print("==> MASTER (exec): sending command bcast")
        mpi_comm.bcast('exec', root=0)
        print("==> MASTER (exec): sending data bcast, data={}".format((data, args, kwargs)))
        mpi_comm.bcast((data, args, kwargs), root=0)
    do_exec(data, args, kwargs)
    worker_responses = [None]
    if mpi_comm is not None:
        print("==> MASTER (exec): waiting for response gather")
        worker_responses = mpi_comm.gather(None, root=0)
    worker_responses[0] = nest.hl_api.serializable(data)
    print("==> MASTER (call): worker_responses={}".format(worker_responses))
    # TODO: combine worker_response data in a meaningful way

    return jsonify(data)


# --------------------------
# RESTful API
# --------------------------

nest_calls = dir(nest)
nest_calls = list(filter(lambda x: not x.startswith('_'), nest_calls))
nest_calls.sort()


@app.route('/api', methods=['GET'])
@cross_origin()
def nest_api():
    """ Route to list call functions in NEST.
    """
    data = init_data(request)
    response = api_client(nest_calls, data)
    return jsonify(response)


@app.route('/api/<call>', methods=['GET', 'POST'])
@cross_origin()
def nest_api_call(call):
    """ Route to call function in NEST.
    """
    data = init_data(request, call)
    args, kwargs = get_arguments(request, data)
    if call in nest_calls:
        call = getattr(nest, call)
        response = api_client(call, data, *args, **kwargs)
    else:
        data['response']['msg'] = 'The request cannot be called in NEST.'
        data['response']['status'] = 'error'
        response = data
    return jsonify(response)


# ----------------------
# Helpers for the server
# ----------------------

class Capturing(list):
    """ Monitor stdout contents i.e. print.
    """
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = io.StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout


def init_data(request, call=None):
    """ Create data variable in dictionary for JSON response.
    """
    url = request.url_rule.rule.split('/')[1]
    data = {
        'request': {
            'url': url,
        },
        'response': {}
    }
    if call:
        data['request']['call'] = call
    return data


def get_arguments(request, data):
    """ Get arguments from the request.
    """
    args, kwargs = [], {}
    if request.is_json:
        json = data['request']['json'] = request.get_json()
        if isinstance(json, list):
            args = json
        elif isinstance(json, dict):
            args = json.get('args', args)
            kwargs = json.get('kwargs', json)
    elif len(request.form) > 0:
        if 'args' in request.form:
            args = data['request']['form'] = request.form.getlist('args')
        else:
            kwargs = data['request']['form'] = request.form.to_dict()
    elif (len(request.args) > 0):
        if 'args' in request.args:
            args = data['request']['args'] = request.args.getlist('args')
        else:
            kwargs = data['request']['args'] = request.args.to_dict()
    return args, kwargs


def get_or_error(func):
    """ Wrapper to get data or error content.
    """
    def func_wrapper(call, data, *args, **kwargs):
        try:
            data = func(call, data, *args, **kwargs)
            if 'data' not in data['response']:
                return data
            response = data['response']['data']
            data['response']['status'] = 'ok'
        except Exception as e:
            data['response']['msg'] = str(e)
            data['response']['status'] = 'error'
        return data
    return func_wrapper


def NodeCollection(kwargs):
    """ Get Node Collection as arguments for NEST functions.
    """
    keys = ['nodes', 'source', 'target', 'pre', 'post']
    for key in keys:
        if key in kwargs:
            kwargs[key] = nest.NodeCollection(kwargs[key])
    return kwargs


def serialize(call, kwargs):
    """ Serialize arguments with keywords for call functions in NEST.
    """
    kwargs = NodeCollection(kwargs)
    if call.startswith('Set'):
        status = {}
        if call == 'SetDefaults':
            status = nest.GetDefaults(kwargs['model'])
        elif call == 'SetKernelStatus':
            status = nest.GetKernelStatus()
        elif call == 'SetStructuralPlasticityStatus':
            status = nest.GetStructuralPlasticityStatus(kwargs['params'])
        elif call == 'SetStatus':
            status = nest.GetStatus(kwargs['nodes'])
        for key, val in kwargs['params'].items():
            if key in status:
                kwargs['params'][key] = type(status[key])(val)
    return kwargs


@get_or_error
def api_client(call, data, *args, **kwargs):
    """ API Client to call function in NEST.
    """
    if callable(call):
        data['request']['call'] = call.__name__
        if str(kwargs.get('return_doc', 'false')) == 'true':
            response = call.__doc__
        elif str(kwargs.get('return_source', 'false')) == 'true':
            response = inspect.getsource(call)
        else:
            if mpi_comm is not None:
                print("==> MASTER (call): sending command bcast")
                mpi_comm.bcast('call', root=0)
                print("==> MASTER (call): sending data bcast, data={}".format((data, args, kwargs)))
                mpi_comm.bcast((data, args, kwargs), root=0)
            response = call(*args, **serialize(call.__name__, kwargs))
            worker_responses = [None]
            if mpi_comm is not None:
                print("==> MASTER (call): waiting for response gather")
                worker_responses = mpi_comm.gather(None, root=0)
            worker_responses[0] = nest.hl_api.serializable(response)
            print("==> MASTER (call): worker_responses={}".format(worker_responses))
            # TODO: combine worker_response in a meaningful way
    else:
        data['request']['call'] = call
        response = call
    data['response']['data'] = nest.hl_api.serializable(response)
    return data


def run_mpi_app(comm):
    global mpi_comm
    mpi_comm = comm
    # NEST segfaults if someone messes with the number of threads.
    # Don't do this!
    app.run(threaded=False)
