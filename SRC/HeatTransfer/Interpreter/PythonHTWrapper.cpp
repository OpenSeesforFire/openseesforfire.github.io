/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

// Written: Minjie

// Description: A python HTwrapper for OpenSees commands
//

#include "PythonHTWrapper.h"
#include "HeatTransferModule.h"
#include <OPS_Globals.h>

static PythonHTWrapper* HTwrapper = 0;

PythonHTWrapper::PythonHTWrapper()
    :currentArgv(0), currentArg(0), numberArgs(0),
     methodsOpenSees(), opensees_docstring(""), currentResult(0)
{
    HTwrapper = this;
}

PythonHTWrapper::~PythonHTWrapper()
{
    HTwrapper = 0;
}

void
PythonHTWrapper::resetCommandLine(int nArgs, int cArg, PyObject* argv)
{
    numberArgs = nArgs;
    currentArg = cArg-1;
    if (currentArg < 0) currentArg = 0;
    currentArgv = argv;
}

void
PythonHTWrapper::resetCommandLine(int cArg)
{
    if (cArg < 0) {
	currentArg += cArg;
    } else {
	currentArg = cArg-1;
    }
    if (currentArg < 0) currentArg = 0;
}

void
PythonHTWrapper::addCommand(const char* name, PyCFunction proc)
{
    PyMethodDef method = {name,proc,METH_VARARGS,opensees_docstring};
    methodsOpenSees.push_back(method);
}

PyMethodDef*
PythonHTWrapper::getMethods()
{
    if (methodsOpenSees.empty()) {
	return 0;
    }

    return &methodsOpenSees[0];
}

void
PythonHTWrapper::setOutputs(int* data, int numArgs)
{
    if (numArgs == 0) return;
    if (numArgs == 1) {
	currentResult = Py_BuildValue("i", data[0]);
	return ;
    }
    currentResult = PyList_New(numArgs);
    for (int i=0; i<numArgs; i++) {
	PyList_SET_ITEM(currentResult, i, Py_BuildValue("i", data[i]));
    }
}

void
PythonHTWrapper::setOutputs(double* data, int numArgs)
{
    if (numArgs == 0) return;
    if (numArgs == 1) {
	currentResult = Py_BuildValue("d", data[0]);
	return ;
    }
    currentResult = PyList_New(numArgs);
    for (int i=0; i<numArgs; i++) {
	PyList_SET_ITEM(currentResult, i, Py_BuildValue("d", data[i]));
    }
}

void
PythonHTWrapper::setOutputs(const char* str)
{
    currentResult = Py_BuildValue("s", str);
}

PyObject*
PythonHTWrapper::getResults()
{
    PyObject* result = currentResult;
    currentResult = 0;

    if (result == 0) {
	Py_INCREF(Py_None);
	result = Py_None;
    }

    return result;
}

//////////////////////////////////////////////
/////// Python HTwrapper functions  ////////////
/////////////////////////////////////////////
static PyObject *Py_ops_addHTMaterial(PyObject *self, PyObject *args)
{
    HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_addHTMaterial() < 0) return NULL;

    return HTwrapper->getResults();
}

static PyObject *Py_ops_addHTEntity(PyObject *self, PyObject *args)
{
    HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_addHTEntity() < 0) return NULL;

    return HTwrapper->getResults();
}

static PyObject* Py_ops_addHTMesh(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHTMesh() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_HTMeshAll(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTMeshAll() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_SetInitialT(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_SetInitialT() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_addHTConstants(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHTConstants() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_HTNodeSet(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTNodeSet() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_HTEleSet(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTEleSet() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_addHTPattern(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHTPattern() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_addFireModel(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addFireModel() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_addHeatFluxBC(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHeatFluxBC() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_addMPTemperatureBC(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addMPTemperatureBC() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_HTRefineMesh(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTRefineMesh() < 0) return NULL;

	return HTwrapper->getResults();
}


static PyObject* Py_ops_HTReset(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTReset() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_HTRecorder(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTRecorder() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_HTAnalysis(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTAnalysis() < 0) return NULL;

	return HTwrapper->getResults();
}

static PyObject* Py_ops_HTAnalyze(PyObject* self, PyObject* args)
{
	HTwrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTAnalyze() < 0) return NULL;

	return HTwrapper->getResults();
}

/////////////////////////////////////////////////
////////////// Add Python commands //////////////
/////////////////////////////////////////////////
void
PythonHTWrapper::addOpenSeesCommands()
{
    addCommand("HTMaterial", &Py_ops_addHTMaterial);
    addCommand("HTEntity", &Py_ops_addHTEntity);
    addCommand("HTMesh", &Py_ops_addHTMesh);
    addCommand("HTMeshAll", &Py_ops_HTMeshAll);
    addCommand("SetInitialT", &Py_ops_SetInitialT);
    addCommand("HTRefineMesh", &Py_ops_HTRefineMesh);
    addCommand("HTConstants", &Py_ops_addHTConstants);
    addCommand("HTEleSet", &Py_ops_HTEleSet);
    addCommand("HTNodeSet", &Py_ops_HTNodeSet);
    addCommand("HTPattern", &Py_ops_addHTPattern);
    addCommand("FireModel", &Py_ops_addFireModel);
    addCommand("HeatFluxBC", &Py_ops_addHeatFluxBC);
    addCommand("timeSeries", &Py_ops_timeSeries);
    addCommand("HTCoupleT", &Py_ops_addMPTemperatureBC);
    addCommand("HTAnalysis", &Py_ops_HTAnalysis);
    addCommand("HTRecorder", &Py_ops_HTRecorder);
    addCommand("HTAnalyze", &Py_ops_HTAnalyze);
    addCommand("HTReset", &Py_ops_HTReset);



    PyMethodDef method = {NULL,NULL,0,NULL};
    methodsOpenSees.push_back(method);
}
