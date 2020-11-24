.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Pipelines
=========

Test Tool
---------
.. automodule:: vre_template_tool

   This is a demonstration pipeline using the testTool.

   Running from the command line
   =============================

   Parameters
   ----------
   config : file
      Location of the config file for the workflow
   in_metadata : file
      Location of the input list of files required by the process
   out_metadata : file
      Location of the output results.json file for returned files

   Returns
   -------
   output : file
      Text file with a single entry

   Example
   -------
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):

   .. code-block:: none
      :linenos:

      cd /home/compss/code/vre_template_tool
      runcompss --lang=python vre_template_tool/VRE_CWL_RUNNER.py --config /vre_template_tool/tests/basic/config.json --in_metadata /vre_template_tool/tests/basic/in_metadata.json --out_metadata /vre_template_tool/tests/basic/results.json
